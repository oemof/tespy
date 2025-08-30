# -*- coding: utf-8

"""Module of class CombustionChamber.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/combustion/base.py

SPDX-License-Identifier: MIT
"""
import itertools

import CoolProp.CoolProp as CP
import numpy as np

from tespy.components.component import Component
from tespy.components.component import component_registry
from tespy.tools import logger
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.fluid_properties import s_mix_pT
from tespy.tools.fluid_properties.helpers import fluid_structure
from tespy.tools.global_vars import FLUID_ALIASES
from tespy.tools.global_vars import combustion_gases
from tespy.tools.helpers import TESPyComponentError
from tespy.tools.helpers import _numeric_deriv
from tespy.tools.helpers import _numeric_deriv_vecvar
from tespy.tools.helpers import fluidalias_in_list


@component_registry
class CombustionChamber(Component):
    r"""
    The class CombustionChamber is parent class of all combustion components.

    **Mandatory Equations**

    - :py:meth:`tespy.components.combustion.base.CombustionChamber.mass_flow_func`
    - :py:meth:`tespy.components.combustion.base.CombustionChamber.combustion_pressure_structure_matrix`
    - :py:meth:`tespy.components.combustion.base.CombustionChamber.stoichiometry`
    - :py:meth:`tespy.components.combustion.base.CombustionChamber.energy_balance_func`

    **Optional Equations**

    - :py:meth:`tespy.components.combustion.base.CombustionChamber.lambda_func`
    - :py:meth:`tespy.components.combustion.base.CombustionChamber.ti_func`

    Available fuels

    - methane, ethane, propane, butane, hydrogen, carbon monoxide, nDodecane

    Inlets/Outlets

    - in1, in2
    - out1

    Image

    .. image:: /api/_images/CombustionChamber.svg
       :alt: flowsheet of the combustion chamber
       :align: center
       :class: only-light

    .. image:: /api/_images/CombustionChamber_darkmode.svg
       :alt: flowsheet of the combustion chamber
       :align: center
       :class: only-dark

    .. note::

        The fuel and the air components can be connected to either of the
        inlets.

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    design_path : str
        Path to the components design case.

    local_offdesign : boolean
        Treat this component in offdesign mode in a design calculation.

    local_design : boolean
        Treat this component in design mode in an offdesign calculation.

    char_warnings : boolean
        Ignore warnings on default characteristics usage for this component.

    printout : boolean
        Include this component in the network's results printout.

    lamb : float, dict
        Actual oxygen to stoichiometric oxygen ratio, :math:`\lambda/1`.

    ti : float, dict
        Thermal input, (:math:`{LHV \cdot \dot{m}_f}`), :math:`ti/\text{W}`.

    Note
    ----
    For more information on the usage of the combustion chamber see the
    examples section on github or look for the combustion chamber tutorials
    at tespy.readthedocs.io.

    Example
    -------
    The combustion chamber calculates energy input due to combustion as well as
    the flue gas composition based on the type of fuel and the amount of
    oxygen supplied. In this example a mixture of methane, hydrogen and
    carbondioxide is used as fuel.

    >>> from tespy.components import Sink, Source, CombustionChamber
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> from tespy.tools.fluid_properties import T_sat_p
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "temperature": "degC"
    ... })
    >>> amb = Source('ambient air')
    >>> sf = Source('fuel')
    >>> fg = Sink('flue gas outlet')
    >>> comb = CombustionChamber('combustion chamber')
    >>> amb_comb = Connection(amb, 'out1', comb, 'in1')
    >>> sf_comb = Connection(sf, 'out1', comb, 'in2')
    >>> comb_fg = Connection(comb, 'out1', fg, 'in1')
    >>> nw.add_conns(sf_comb, amb_comb, comb_fg)

    Specify the thermal input of the combustion chamber. At the given fluid
    compositions this determines the mass flow of the fuel. The outlet
    temperature of the flue gas determines the ratio of oxygen to fuel mass
    flow.

    >>> comb.set_attr(ti=500000)
    >>> amb_comb.set_attr(p=1, T=20, fluid={'Ar': 0.0129, 'N2': 0.7553,
    ... 'CO2': 0.0004, 'O2': 0.2314})
    >>> sf_comb.set_attr(T=25, fluid={'CO2': 0.03, 'H2': 0.01, 'CH4': 0.96})
    >>> comb_fg.set_attr(T=1200)
    >>> nw.solve('design')
    >>> round(comb.lamb.val, 3)
    2.014
    >>> comb.set_attr(lamb=2)
    >>> comb_fg.set_attr(T=None)
    >>> nw.solve('design')
    >>> round(comb_fg.T.val, 1)
    1206.6
    """

    def get_parameters(self):
        return {
            'lamb': dc_cp(
                min_val=1,
                func=self.lambda_func,
                dependents=self.lambda_dependents,
                num_eq_sets=1,
                quantity="ratio"
            ),
            'ti': dc_cp(
                min_val=0,
                deriv=self.ti_deriv,
                func=self.ti_func,
                dependents=self.ti_dependents,
                num_eq_sets=1,
                quantity="heat"
            )
        }

    def _update_num_eq(self):
        inl, outl = self._get_combustion_connections()
        fluid_eqs = set([f for c in inl + outl for f in c.fluid.val])
        self.fluid_eqs_list = list(fluid_eqs)
        self.constraints["stoichiometry_constraints"].num_eq = len(self.fluid_eqs_list)

    def get_mandatory_constraints(self):
        return {
            'mass_flow_constraints': dc_cmc(**{
                'func': self.mass_flow_func,
                'deriv': self.mass_flow_deriv,
                'dependents': self.mass_flow_dependents,
                'constant_deriv': True,
                'num_eq_sets': 1
            }),
            'reactor_pressure_constraints': dc_cmc(**{
                'structure_matrix': self.combustion_pressure_structure_matrix,
                'num_eq_sets': 2
            }),
            'stoichiometry_constraints': dc_cmc(**{
                'func': self.stoichiometry_func,
                'deriv': self.stoichiometry_deriv,
                'dependents': self.stoichiometry_dependents,
                'constant_deriv': False,
                'num_eq_sets': 1
            }),
            'energy_balance_constraints': dc_cmc(**{
                'func': self.energy_balance_func,
                'dependents': self.energy_balance_dependents,
                'constant_deriv': False,
                'num_eq_sets': 1
            })
        }

    @staticmethod
    def inlets():
        return ['in1', 'in2']

    @staticmethod
    def outlets():
        return ['out1']

    def propagate_wrapper_to_target(self, branch):
        if self in branch["components"]:
            return

        outconn = self.outl[0]
        branch["connections"] += [outconn]
        branch["components"] += [self]
        outconn.target.propagate_wrapper_to_target(branch)

    def _preprocess(self, num_nw_vars):
        super()._preprocess(num_nw_vars)
        self.setup_reaction_parameters()

    def _get_combustion_connections(self):
        return (self.inl[:2], [self.outl[0]])

    def setup_reaction_parameters(self):
        r"""Setup parameters for reaction (gas name aliases and LHV)."""
        self.fuel_list = []
        all_fluids = set([f for c in self.inl + self.outl for f in c.fluid.val])
        for f in all_fluids:
            if fluidalias_in_list(f, combustion_gases):
                self.fuel_list += [f]

        self.fuel_list = set(self.fuel_list)

        if len(self.fuel_list) == 0:
            msg = (
                "Your network's fluids do not contain any fuels, that are "
                f"available for the component {self.label} of type "
                f"{self.__class__.__name__}. Available fuels are: "
                f"{', '.join(combustion_gases)}."
            )
            logger.error(msg)
            raise TESPyComponentError(msg)

        else:
            msg = (
                f"The fuels for component {self.label} of type "
                f"{self.__class__.__name__} are: {', '.join(self.fuel_list)}."
            )
            logger.debug(msg)

        for fluid in ["O2", "CO2", "H2O", "N2"]:
            if not fluidalias_in_list(fluid, all_fluids):
                aliases = ", ".join(FLUID_ALIASES.get_fluid(fluid))
                msg = (
                    f"The component {self.label} (class "
                    f"{self.__class__.__name__}) requires that the fluid "
                    f"{fluid} (aliases: {aliases}) is in the network's list of "
                    "fluids."
                )
                logger.error(msg)
                raise TESPyComponentError(msg)
            else:
                setattr(self, fluid.lower(), fluid)

        self.fuels = {}
        for f in self.fuel_list:
            self.fuels[f] = {}
            structure = fluid_structure(f)
            for el in ['C', 'H', 'O']:
                if el in structure:
                    self.fuels[f][el] = structure[el]
                else:
                    self.fuels[f][el] = 0
            self.fuels[f]['LHV'] = self.calc_lhv(f)

    def calc_lhv(self, f):
        r"""
        Calculate the lower heating value of the combustion chamber's fuel.

        - Source for fluids O2, H2O and CO2: :cite:`CODATA1989`
        - Source for all other fluids: :cite:`CRCHandbook2021`

        Parameters
        ----------
        f : str
            Alias of the fuel.

        Returns
        -------
        val : float
            Lower heating value of the combustion chambers fuel.

            .. math::

                LHV = -\frac{\sum_i {\Delta H_f^0}_i -
                \sum_j {\Delta H_f^0}_j }
                {M_{fuel}}\\
                \forall i \in \text{reation products},\\
                \forall j \in \text{reation educts},\\
                \Delta H_f^0: \text{molar formation enthalpy}
        """
        inl, _ = self._get_combustion_connections()
        hf = {}
        hf['hydrogen'] = 0
        hf['methane'] = -74.6
        hf['ethane'] = -84.0
        hf['propane'] = -103.8
        hf['butane'] = -125.7
        hf['nDodecane'] = -289.4
        hf['CO'] = -110.5
        hf[self.o2] = 0
        hf[self.co2] = -393.51
        # water (gaseous)
        hf[self.h2o] = -241.826

        key = set(list(hf.keys())).intersection(FLUID_ALIASES.get_fluid(f))

        val = (
            -(
                self.fuels[f]['H'] / 2 * hf[self.h2o]
                + self.fuels[f]['C'] * hf[self.co2]
                - (
                    (
                        self.fuels[f]['C'] + self.fuels[f]['H'] / 4
                        - self.fuels[f]['O'] / 2
                    ) * hf[self.o2] + hf[list(key)[0]]
                )
            ) / inl[0].fluid.wrapper[f]._molar_mass * 1000
        )

        return val

    def _add_missing_fluids(self, connections):
        inl, outl = self._get_combustion_connections()
        if set(inl + outl) & set(connections):
            return ["H2O", "CO2"]
        else:
            return super()._add_missing_fluids(connections)

    def mass_flow_func(self):
        r"""
        Calculate the residual value for component's mass flow balance.

        Returns
        -------
        residual : list
            Vector with residual value for component's mass flow balance.

            .. math::

                0 = \dot{m}_{in,1} + \dot{m}_{in,2} - \dot{m}_{out,1}
        """
        inl, outl = self._get_combustion_connections()
        return inl[0].m.val_SI + inl[1].m.val_SI - outl[0].m.val_SI

    def mass_flow_deriv(self, increment_filter, k, dependents=None):
        r"""
        Calculate the partial derivatives for all mass flow balance equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the fluid equations.
        """
        inl, outl = self._get_combustion_connections()
        for i in inl:
            if i.m.is_var:
                self.jacobian[k, i.m.J_col] = 1

        if outl[0].m.is_var:
            self.jacobian[k, outl[0].m.J_col] = -1

    def mass_flow_dependents(self):
        inl, outl = self._get_combustion_connections()
        return [[c.m for c in inl + outl]]

    def combustion_pressure_structure_matrix(self, k):
        r"""
        Equations for reactor pressure balance.

        Returns
        -------
        residual : list
            Residual values of equations.

            .. math::

                0 = p_\mathrm{in,1} - p_\mathrm{out,1}\\
                0 = p_\mathrm{in,2} - p_\mathrm{out,1}
        """
        inl, outl = self._get_combustion_connections()

        self._structure_matrix[k, inl[0].p.sm_col] = 1
        self._structure_matrix[k, outl[0].p.sm_col] = -1

        self._structure_matrix[k + 1, inl[1].p.sm_col] = 1
        self._structure_matrix[k + 1, outl[0].p.sm_col] = -1

    def stoichiometry_func(self):
        r"""
        Collect residual values for all fluids in stoichiometry.

        Returns
        -------
        residual : list
            Vector with residual values of equations.
        """
        # calculate equations
        residual = []
        for fluid in self.fluid_eqs_list:
            residual += [self.stoichiometry(fluid)]

        return residual

    def stoichiometry(self, fluid):
        r"""
        Calculate the reaction balance for one fluid.

        - determine molar mass flows of fuel and oxygen
        - calculate mole number of carbon and hydrogen atoms in fuel
        - calculate molar oxygen flow for stoichiometric combustion
        - calculate residual value for the corresponding fluid

        for excess fuel

        - calculate excess carbon and hydrogen in fuels
        - calculate excess fuel shares

        General equations

        .. math::

            res = \sum_i \left(x_{fluid,i} \cdot \dot{m}_{i}\right) -
            \sum_j \left(x_{fluid,j} \cdot \dot{m}_{j}\right)\\
            \forall i \in \text{combustion inlets}\\
            \forall j \in text{flue gas outlet}

            \dot{m}_{fluid,m} = \sum_i \frac{x_{fluid,i} \cdot \dot{m}_{i}}
            {M_{fluid}}\\
            \forall i \in \text{combustion inlets}

            \dot{m}_{O_2,m,stoich}=\frac{\dot{m}_{H_m}}{4} + \dot{m}_{C_m}

            \lambda = \frac{\dot{m}_{O_2,m}}{\dot{m}_{O_2,m,stoich}}

        Excess carbon and hydrogen

        .. math::

           \dot{m}_{H_{exc,m}} = \begin{cases}
           0 & \lambda \geq 1\\
           4 \cdot \left( \dot{m}_{O_2,m,stoich} -
           \dot{m}_{O_2,m}\right) & \lambda < 1
            \end{cases}

           \dot{m}_{C_{exc,m}} = \begin{cases}
           0 & \lambda \geq 1\\
           \dot{m}_{O_2,m,stoich} - \dot{m}_{O_2,m} & \lambda < 1
            \end{cases}

        Equation for fuels

        .. math::

            0 = res - \left(\dot{m}_{f,m} - \dot{m}_{f,exc,m}\right)
            \cdot M_{fuel}\\

            \dot{m}_{f,exc,m} = \begin{cases}
            0 & \lambda \geq 1\\
            \dot{m}_{f,m} - \frac{\dot{m}_{O_2,m}}
            {n_{C,fuel} + 0.25 \cdot n_{H,fuel}}
            \end{cases}

        Equation for oxygen

        .. math::

            0 = res - \begin{cases}
            -\frac{\dot{m}_{O_2,m} \cdot M_{O_2}}{\lambda} &
            \lambda \geq 1\\
            - \dot{m}_{O_2,m} \cdot M_{O_2} & \lambda < 1
            \end{cases}

        Equation for water

        .. math::

            0 = res + \left( \dot{m}_{H_m} - \dot{m}_{H_{exc,m}} \right)
            \cdot 0.5 \cdot M_{H_2O}

        Equation for carbondioxide

        .. math::

            0 = res + \left( \dot{m}_{C_m} - \dot{m}_{C_{exc,m}} \right)
            \cdot M_{CO_2}

        Equation for all other fluids

        .. math::

            0 = res

        Parameters
        ----------
        fluid : str
            Fluid to calculate residual value for.

        Returns
        -------
        residual : float
            Residual value for corresponding fluid.
        """
        # required to work with combustion chamber and engine
        inl, outl = self._get_combustion_connections()

        if fluid in list(self.fuel_list) + [self.co2, self.o2, self.h2o]:
            ###################################################################
            # molar mass flow for fuel and oxygen
            n_fuel = {}
            n_oxy_stoich = {}
            n_h = 0
            n_c = 0
            n_o = 0
            for f in self.fuel_list:
                n_fuel[f] = 0
                for i in inl:
                    n = (
                        i.m.val_SI * i.fluid.val.get(f, 0)
                        / inl[0].fluid.wrapper[f]._molar_mass
                    )
                    n_fuel[f] += n
                    n_h += n * self.fuels[f]['H']
                    n_c += n * self.fuels[f]['C']
                    n_o += n * self.fuels[f]['O']

                # stoichiometric oxygen requirement for each fuel
                n_oxy_stoich[f] = n_fuel[f] * (
                    self.fuels[f]['H'] / 4 + self.fuels[f]['C']
                    - self.fuels[f]['O'] / 2
                )

            ###################################################################
            # calculate stoichiometric oxygen
            n_oxygen_stoich = n_h / 4 + n_c - n_o / 2

            n_oxygen = 0
            for i in inl:
                n_oxygen += (
                    i.m.val_SI
                    * i.fluid.val.get(self.o2, 0)
                    / inl[0].fluid.wrapper[self.o2]._molar_mass
                )

            ###################################################################
            # calculate lambda if not set
            if not self.lamb.is_set:
                self.lamb.val_SI = n_oxygen / n_oxygen_stoich

            ###################################################################
            # calculate excess fuel if lambda is lower than 1
            if self.lamb.val_SI < 1:
                n_h_exc = (n_oxygen_stoich - n_oxygen) * 4
                n_c_exc = (n_oxygen_stoich - n_oxygen)
            else:
                n_h_exc = 0
                n_c_exc = 0

        ###################################################################
        # equation for carbondioxide
        if fluid == self.co2:
            dm = (n_c - n_c_exc) * inl[0].fluid.wrapper[self.co2]._molar_mass

        ###################################################################
        # equation for water
        elif fluid == self.h2o:
            dm = (n_h - n_h_exc) / 2 * inl[0].fluid.wrapper[self.h2o]._molar_mass

        ###################################################################
        # equation for oxygen
        elif fluid == self.o2:
            if self.lamb.val_SI < 1:
                dm = -n_oxygen * inl[0].fluid.wrapper[self.o2]._molar_mass
            else:
                dm = -n_oxygen / self.lamb.val_SI * inl[0].fluid.wrapper[self.o2]._molar_mass

        ###################################################################
        # equation for fuel
        elif fluid in self.fuel_list:
            if self.lamb.val_SI < 1:
                n_fuel_exc = (
                    -(n_oxygen / n_oxygen_stoich - 1) * n_oxy_stoich[fluid]
                    / (self.fuels[fluid]['H'] / 4 + self.fuels[fluid]['C'] - self.fuels[fluid]['O'] / 2)
                )
            else:
                n_fuel_exc = 0
            dm = -(n_fuel[fluid] - n_fuel_exc) * inl[0].fluid.wrapper[fluid]._molar_mass

        ###################################################################
        # equation for other fluids
        else:
            dm = 0

        res = dm
        for i in inl:
            res += i.fluid.val.get(fluid, 0) * i.m.val_SI
        res -= outl[0].fluid.val.get(fluid, 0) * outl[0].m.val_SI
        return res

    def stoichiometry_deriv(self, increment_filter, k, dependents=None):
        r"""
        Calculate partial derivatives of the reaction balance.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of equation in Jacobian matrix.
        """
        # required to work with combustion chamber and engine
        inl, outl = self._get_combustion_connections()
        f = self.stoichiometry
        conns = inl + outl
        for fluid, conn in itertools.product(self.fluid_eqs_list, conns):
            eq_num = self.fluid_eqs_list.index(fluid)
            self._partial_derivative(conn.m, k + eq_num, f, fluid=fluid)
            for fluid_name in conn.fluid.is_var:
                self._partial_derivative_fluid(
                    conn.fluid, k + eq_num, f, fluid_name, fluid=fluid
                )

        # dependency on outlet state is super simple!
        # TODO: make sure inlet and outlet mass flows and fluids are not
        # linear dependent, then it is more difficult!
        for fluid in outl[0].fluid.is_var:
            eq_num = self.fluid_eqs_list.index(fluid)
            self._partial_derivative(outl[0].m, k + eq_num, -outl[0].fluid.val[fluid])
            self.jacobian[k + eq_num, outl[0].fluid.J_col[fluid]] = -outl[0].m.val_SI

    def stoichiometry_dependents(self):
        inl, outl = self._get_combustion_connections()
        mass_flows = [c.m for c in inl + outl]
        return {
            "scalars": [
                mass_flows for f in self.fluid_eqs_list
            ],
            "vectors": [
                {c.fluid: c.fluid.is_var for c in inl}
                | {outl[0].fluid: outl[0].fluid.is_var & {f}}
                for f in self.fluid_eqs_list
            ]
        }

    def energy_balance_func(self):
        r"""
        Calculate the energy balance of the adiabatic combustion chamber.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                \begin{split}
                0 = & \sum_i \dot{m}_{in,i} \cdot
                \left( h_{in,i} - h_{in,i,ref} \right)\\
                & -\dot{m}_{out,2}\cdot\left( h_{out,1}-h_{out,1,ref} \right)\\
                & + LHV_{fuel} \cdot\left(\sum_i\dot{m}_{in,i}\cdot
                x_{fuel,in,i}- \dot{m}_{out,1} \cdot x_{fuel} \right)
                \end{split}\\

                \forall i \in \text{inlets}

        Note
        ----
        The temperature for the reference state is set to 25 °C, thus
        the water may be liquid. In order to make sure, the state is
        referring to the lower heating value, the state of the water in the
        flue gas is fored to gaseous.

        - Reference temperature: 298.15 K.
        - Reference pressure: 1 bar.
        """
        inl, outl = self._get_combustion_connections()
        T_ref = 298.15
        p_ref = 1e5

        res = 0
        for i in inl:
            res += i.m.val_SI * (
                i.h.val_SI
                - h_mix_pT(p_ref, T_ref, i.fluid_data, mixing_rule="forced-gas")
            )

        for o in outl:
            res -= o.m.val_SI * (
                o.h.val_SI
                - h_mix_pT(p_ref, T_ref, o.fluid_data, mixing_rule="forced-gas")
            )

        res += self.calc_ti()
        return res

    def energy_balance_dependents(self):
        inl, outl = self._get_combustion_connections()
        return {
            "scalars": [var for c in inl + outl for var in [c.m, c.h]],
            "vectors": [{
                c.fluid: self.fuel_list & c.fluid.is_var for c in inl + outl
            }]
        }

    def lambda_func(self):
        r"""
        Calculate the residual for specified lambda.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \frac{\dot{m}_{f,m}}{\dot{m}_{O_2,m} \cdot
                \left(n_{C,fuel} + 0.25 \cdot n_{H,fuel}\right)} - \lambda

                \dot{m}_{fluid,m} = \sum_i \frac{x_{fluid,i} \cdot \dot{m}_{i}}
                {M_{fluid}}\\ \forall i \in inlets
        """
        return self.calc_lambda() - self.lamb.val_SI

    def lambda_dependents(self):
        inl, _ = self._get_combustion_connections()
        return {
            "scalars": [c.m for c in inl],
            "vectors": [{
                c.fluid: (self.fuel_list | {self.o2}) & c.fluid.is_var
                for c in inl
            }]
        }

    def calc_lambda(self):
        r"""
        Calculate oxygen to stoichimetric oxygen ration

        Returns
        -------
        lambda : float
            Oxygent to stoichiometric oxygen ratio.

            .. math::

                \lambda = \frac{\dot{m}_{f,m}}{\dot{m}_{O_2,m} \cdot
                \left(n_{C,fuel} + 0.25 \cdot n_{H,fuel}\right)}

                \dot{m}_{fluid,m} = \sum_i \frac{x_{fluid,i} \cdot \dot{m}_{i}}
                {M_{fluid}}\\ \forall i \in inlets
        """
        inl, _ = self._get_combustion_connections()

        n_h = 0
        n_c = 0
        n_o = 0
        for f in self.fuel_list:
            for i in inl:
                n_fuel = (
                    i.m.val_SI * i.fluid.val.get(f, 0)
                    / inl[0].fluid.wrapper[f]._molar_mass
                )
                n_h += n_fuel * self.fuels[f]['H']
                n_c += n_fuel * self.fuels[f]['C']
                n_o += n_fuel * self.fuels[f]['O']

        n_oxygen = 0
        for i in inl:
            n_oxygen += (
                i.m.val_SI * i.fluid.val.get(self.o2, 0)
                / inl[0].fluid.wrapper[self.o2]._molar_mass
            )
        n_oxygen_stoich = n_h / 4 + n_c - n_o / 2
        return n_oxygen / n_oxygen_stoich

    def ti_func(self):
        r"""
        Calculate the residual for specified thermal input.

        Returns
        -------
        residual : float
            Residual value of function.

            .. math::

                0 = ti - \dot{m}_{fuel} \cdot LHV
        """
        return self.ti.val_SI - self.calc_ti()

    def ti_deriv(self, increment_filter, k, dependents=None):
        """
        Calculate partial derivatives of thermal input function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of equation in Jacobian matrix.
        """
        inl, outl = self._get_combustion_connections()
        for i in inl:
            if i.m.is_var:
                deriv = 0
                for f in self.fuel_list:
                    deriv -= i.fluid.val.get(f, 0) * self.fuels[f]['LHV']
                self.jacobian[k, i.m.J_col] = deriv
            for f in (self.fuel_list & i.fluid.is_var):
                self.jacobian[k, i.fluid.J_col[f]] = -i.m.val_SI * self.fuels[f]['LHV']

        o = outl[0]
        if o.m.is_var:
            deriv = 0
            for f in self.fuel_list:
                deriv += o.fluid.val.get(f, 0) * self.fuels[f]['LHV']
            self.jacobian[k, o.m.J_col] = deriv
        for f in (self.fuel_list & o.fluid.is_var):
            self.jacobian[k, o.fluid.J_col[f]] = o.m.val_SI * self.fuels[f]['LHV']

    def ti_dependents(self):
        inl, outl = self._get_combustion_connections()
        return {
            "scalars": [c.m for c in inl + outl],
            "vectors": [{
                c.fluid: self.fuel_list & c.fluid.is_var for c in inl + outl
            }]
        }

    def calc_ti(self):
        r"""
        Calculate the thermal input of the combustion chamber.

        Returns
        -------
        ti : float
            Thermal input.

            .. math::

                ti = LHV \cdot \left[\sum_i \left(\dot{m}_{in,i}
                \cdot x_{fuel,in,i}
                \right) - \dot{m}_{out,1} \cdot x_{fuel,out,1} \right]
                \; \forall i \in [1,2]
        """
        inl, outl = self._get_combustion_connections()
        ti = 0
        for f in self.fuel_list:
            m = 0
            for i in inl:
                m += i.m.val_SI * i.fluid.val.get(f, 0)

            for o in outl:
                m -= o.m.val_SI * o.fluid.val.get(f, 0)

            ti += m * self.fuels[f]['LHV']

        return ti

    def bus_func(self, bus):
        r"""
        Calculate the value of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus.Bus
            TESPy bus object.

        Returns
        -------
        val : float
            Value of energy transfer :math:`\dot{E}`. This value is passed to
            :py:meth:`tespy.components.component.Component.calc_bus_value`
            for value manipulation according to the specified characteristic
            line of the bus.

            .. math::

                \dot{E} = LHV \cdot \dot{m}_{f}
        """
        return self.calc_ti()

    def bus_deriv(self, bus):
        r"""
        Calculate the matrix of partial derivatives of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus.Bus
            TESPy bus object.

        Returns
        -------
        deriv : ndarray
            Matrix of partial derivatives.
        """
        inl, outl = self._get_combustion_connections()
        f = self.calc_bus_value
        for c in inl + outl:
            if c.m.is_var:
                if c.m.J_col not in bus.jacobian:
                    bus.jacobian[c.m.J_col] = 0
                bus.jacobian[c.m.J_col] -= _numeric_deriv(c.m._reference_container, f, bus=bus)

            for fl in (self.fuel_list & c.fluid.is_var):
                if c.fluid.J_col[fl] not in bus.jacobian:
                    bus.jacobian[c.fluid.J_col[fl]] = 0
                bus.jacobian[c.fluid.J_col[fl]] -= _numeric_deriv_vecvar(c.fluid._reference_container, f, fl, bus=bus)

    def convergence_check(self):
        r"""
        Perform a convergence check.

        Note
        ----
        Manipulate enthalpies/pressure at inlet and outlet if not specified
        by user to match physically feasible constraints, keep fluid
        composition within feasible range and then propagates it towards the
        outlet.
        """
        # required to work with combustion chamber and engine
        inl, outl = self._get_combustion_connections()

        m = 0
        for i in inl:
            if i.m.val_SI < 0 and i.m.is_var:
                i.m.set_reference_val_SI(0.01)
            m += i.m.val_SI

        ######################################################################
        # check fluid composition
        outl = outl[0]
        for f in outl.fluid.is_var:
            if f == self.o2:
                if outl.fluid.val[f] > 0.25:
                    outl.fluid._reference_container.val[f] = 0.2
                if outl.fluid.val[f] < 0.001:
                    outl.fluid._reference_container.val[f] = 0.05

            elif f == self.co2:
                if outl.fluid.val[f] > 0.15:
                    outl.fluid._reference_container.val[f] = 0.1

            elif f == self.h2o:
                if outl.fluid.val[f] > 0.15:
                    outl.fluid._reference_container.val[f] = 0.1

            elif f in self.fuel_list:
                if outl.fluid.val[f] > 0:
                    outl.fluid._reference_container.val[f] = 0

            else:
                m_f = 0
                for i in inl:
                    m_f += i.fluid.val.get(f, 0) * i.m.val_SI

                if abs(outl.fluid.val[f] - m_f / m) > 0.03:
                    outl.fluid._reference_container.val[f] = m_f / m

        total_mass_fractions = sum(outl.fluid.val.values())
        for fluid in outl.fluid.is_var:
            outl.fluid._reference_container.val[fluid] /= total_mass_fractions

        if outl.m.val_SI < 0 and outl.m.is_var:
            outl.m.set_reference_val_SI(10)

        if not outl.good_starting_values:
            if outl.h.val_SI < 7.5e5 and outl.h.is_var:
                outl.h.set_reference_val_SI(1e6)

        ######################################################################
        # additional checks for performance improvement
        if not self.lamb.is_set and self.lamb.val_SI < 2:
            # search fuel and air inlet
            for i in inl:
                fuel_found = False
                if not i.good_starting_values:
                    fuel = 0
                    for f in self.fuel_list:
                        fuel += i.fluid.val.get(f, 0)
                    # found the fuel inlet
                    if fuel > 0.75 and i.m.is_var:
                        fuel_found = True
                        fuel_inlet = i

                    # found the air inlet
                    if fuel < 0.75:
                        air_tmp = i.m.val_SI

            if fuel_found:
                fuel_inlet.m.set_reference_val_SI(air_tmp / 25)

    @staticmethod
    def initialise_source(c, key):
        r"""
        Return a starting value for pressure and enthalpy at outlet.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            Connection to perform initialisation on.

        key : str
            Fluid property to retrieve.

        Returns
        -------
        val : float
            Starting value for pressure/enthalpy in SI units.

            .. math::

                val = \begin{cases}
                5 \cdot 10^5 & \text{key = 'p'}\\
                10^6 & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 5e5
        elif key == 'h':
            return 10e5

    @staticmethod
    def initialise_target(c, key):
        r"""
        Return a starting value for pressure and enthalpy at inlet.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            Connection to perform initialisation on.

        key : str
            Fluid property to retrieve.

        Returns
        -------
        val : float
            Starting value for pressure/enthalpy in SI units.

            .. math::

                val = \begin{cases}
                5  \cdot 10^5 & \text{key = 'p'}\\
                5 \cdot 10^5 & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 5e5
        elif key == 'h':
            return 5e5

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        self.ti.val_SI = self.calc_ti()
        self.lamb.val_SI = self.calc_lambda()

    def entropy_balance(self):
        r"""
        Calculate entropy balance of combustion chamber.

        Note
        ----
        The entropy balance makes the following parameter available:

        - :code:`T_mcomb`: Thermodynamic temperature of heat of combustion
        - :code:`S_comb`: Entropy production due to combustion
        - :code:`S_irr`: Entropy production due to irreversibilty

        The methodology for entropy analysis of combustion processes is derived
        from :cite:`Tuschy2001`. Similar to the energy balance of a combustion
        reaction, we need to define the same reference state for the entropy
        balance of the combustion. The temperature for the reference state is
        set to 25 °C and reference pressure is 1 bar. As the water in the flue
        gas may be liquid but the thermodynmic temperature of heat of
        combustion refers to the lower heating value, the water is forced to
        gas at the reference point by considering evaporation.

        - Reference temperature: 298.15 K.
        - Reference pressure: 1 bar.

        .. math::

            T_\mathrm{m,comb}= \frac{\dot{m}_\mathrm{fuel} \cdot LHV}
            {\dot{S}_\mathrm{comb}}\\
            \dot{S}_\mathrm{comb}= \dot{m}_\mathrm{fluegas} \cdot
            \left(s_\mathrm{fluegas}-s_\mathrm{fluegas,ref}\right)
            - \sum_{i=1}^2 \dot{m}_{\mathrm{in,}i} \cdot
            \left( s_{\mathrm{in,}i} - s_{\mathrm{in,ref,}i} \right)\\
            \dot{S}_\mathrm{irr}= 0\\
        """
        T_ref = 298.15
        p_ref = 1e5
        o = self.outl[0]
        self.S_comb = o.m.val_SI * (
            o.s.val_SI - s_mix_pT(p_ref, T_ref, o.fluid_data, "forced-gas")
        )

        for i in self.inl:
            self.S_Qcomb -= i.m.val_SI * (
                i.s.val_SI - s_mix_pT(p_ref, T_ref, i.fluid_data, "forced-gas")
            )

        self.S_irr = 0
        self.T_mcomb = self.calc_ti() / self.S_comb

    def exergy_balance(self, T0):
        self.E_P = self.outl[0].Ex_physical - (
            self.inl[0].Ex_physical + self.inl[1].Ex_physical
        )
        self.E_F = (
            self.inl[0].Ex_chemical + self.inl[1].Ex_chemical
            - self.outl[0].Ex_chemical
        )

        self.E_D = self.E_F - self.E_P
        self.epsilon = self._calc_epsilon()
        self.E_bus = np.nan
