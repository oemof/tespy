# -*- coding: utf-8

"""Module of class CombustionChamber.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/combustion/base.py

SPDX-License-Identifier: MIT
"""

import CoolProp.CoolProp as CP
import numpy as np

from tespy.components.component import Component
from tespy.tools import logger
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.document_models import generate_latex_eq
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.fluid_properties import s_mix_pT
from tespy.tools.global_vars import combustion_gases
from tespy.tools.global_vars import molar_masses
from tespy.tools.helpers import TESPyComponentError
from tespy.tools.helpers import fluid_structure
from tespy.tools.helpers import fluidalias_in_list


class CombustionChamber(Component):
    r"""
    The class CombustionChamber is parent class of all combustion components.

    **Mandatory Equations**

    - :py:meth:`tespy.components.combustion.base.CombustionChamber.mass_flow_func`
    - :py:meth:`tespy.components.combustion.base.CombustionChamber.combustion_pressure_func`
    - :py:meth:`tespy.components.combustion.base.CombustionChamber.stoichiometry`
    - :py:meth:`tespy.components.combustion.base.CombustionChamber.energy_balance_func`

    **Optional Equations**

    - :py:meth:`tespy.components.combustion.base.CombustionChamber.lambda_func`
    - :py:meth:`tespy.components.combustion.base.CombustionChamber.ti_func`

    Available fuels

    - methane, ethane, propane, butane, hydrogen

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
    >>> import shutil
    >>> fluid_list = ['Ar', 'N2', 'H2', 'O2', 'CO2', 'CH4', 'H2O']
    >>> nw = Network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ... iterinfo=False)
    >>> amb = Source('ambient air')
    >>> sf = Source('fuel')
    >>> fg = Sink('flue gas outlet')
    >>> comb = CombustionChamber('combustion chamber')
    >>> comb.component()
    'combustion chamber'
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
    ... 'H2O': 0, 'CH4': 0, 'CO2': 0.0004, 'O2': 0.2314, 'H2': 0})
    >>> sf_comb.set_attr(T=25, fluid={'CO2': 0.03, 'H2': 0.01, 'Ar': 0,
    ... 'N2': 0, 'O2': 0, 'H2O': 0, 'CH4': 0.96})
    >>> comb_fg.set_attr(T=1200)
    >>> nw.solve('design')
    >>> round(comb.lamb.val, 3)
    2.014
    >>> comb.set_attr(lamb=2)
    >>> comb_fg.set_attr(T=np.nan)
    >>> nw.solve('design')
    >>> round(comb_fg.T.val, 1)
    1206.6
    """

    @staticmethod
    def component():
        return 'combustion chamber'

    def get_parameters(self):
        return {
            'lamb': dc_cp(
                min_val=1, deriv=self.lambda_deriv, func=self.lambda_func,
                latex=self.lambda_func_doc, num_eq=1),
            'ti': dc_cp(
                min_val=0, deriv=self.ti_deriv, func=self.ti_func,
                latex=self.ti_func_doc, num_eq=1)
        }

    def get_mandatory_constraints(self):
        self.fluid_eqs = set(
            [
                f for c in self.inl + self.outl
                for f in c.fluid.val
            ]
        )
        self.fluid_eqs_list = list(self.fluid_eqs)
        return {
            'mass_flow_constraints': {
                'func': self.mass_flow_func, 'deriv': self.mass_flow_deriv,
                'constant_deriv': True, 'latex': self.mass_flow_func_doc,
                'num_eq': 1},
            'reactor_pressure_constraints': {
                'func': self.combustion_pressure_func,
                'deriv': self.combustion_pressure_deriv,
                'constant_deriv': True,
                'latex': self.combustion_pressure_func_doc,
                'num_eq': 2},
            'stoichiometry_constraints': {
                'func': self.stoichiometry_func,
                'deriv': self.stoichiometry_deriv,
                'constant_deriv': False,
                'latex': self.stoichiometry_func_doc,
                'num_eq': len(self.fluid_eqs)},
            'energy_balance_constraints': {
                'func': self.energy_balance_func,
                'deriv': self.energy_balance_deriv,
                'constant_deriv': False, 'latex': self.energy_balance_func_doc,
                'num_eq': 1}
        }

    @staticmethod
    def inlets():
        return ['in1', 'in2']

    @staticmethod
    def outlets():
        return ['out1']

    @staticmethod
    def is_branch_source():
        return True

    def start_branch(self):
        outconn = self.outl[0]
        for f in ["H2O", "CO2"]:
            if f not in outconn.fluid.val:
                outconn.fluid.val[f] = 0

        branch = {
            "connections": [outconn],
            "components": [self, outconn.target],
            "subbranches": {}
        }
        outconn.target.propagate_to_target(branch)

        return {outconn.label: branch}

    def propagate_to_target(self, branch):
        return

    def propagate_wrapper_to_target(self, branch):
        if self in branch["components"]:
            return

        outconn = self.outl[0]
        branch["connections"] += [outconn]
        branch["components"] += [self]
        outconn.target.propagate_wrapper_to_target(branch)

    def preprocess(self, num_nw_vars):
        super().preprocess(num_nw_vars)
        self.setup_reaction_parameters()

    def setup_reaction_parameters(self):
        r"""Setup parameters for reaction (gas name aliases and LHV)."""
        self.fuel_list = []
        all_fluids = [f for c in self.inl + self.outl for f in c.fluid.val]
        for f in all_fluids:
            if fluidalias_in_list(f, combustion_gases):
                self.fuel_list += [f]

        self.fuel_list = set(self.fuel_list)

        if len(self.fuel_list) == 0:
            msg = (
                "Your network's fluids do not contain any fuels, that are "
                f"available for the component {self.label} of type "
                f"{self.component()}. Available fuels are: " +
                ", ".join(combustion_gases) + "."
            )
            logger.error(msg)
            raise TESPyComponentError(msg)

        else:
            msg = (
                f"The fuels for component {self.label} of type "
                f"{self.component()} are: " + ", ".join(self.fuel_list) + "."
            )
            logger.debug(msg)

        for fluid in ["O2", "CO2", "H2O", "N2"]:
            if not fluidalias_in_list(fluid, all_fluids):
                aliases = ", ".join(CP.get_aliases(fluid))
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
        hf = {}
        hf['hydrogen'] = 0
        hf['methane'] = -74.6
        hf['ethane'] = -84.0
        hf['propane'] = -103.8
        hf['butane'] = -125.7
        hf['nDodecane'] = -289.4
        hf[self.o2] = 0
        hf[self.co2] = -393.51
        # water (gaseous)
        hf[self.h2o] = -241.826

        key = set(list(hf.keys())).intersection(
                set([a.replace(' ', '')
                     for a in CP.get_aliases(f)]))

        val = (-(self.fuels[f]['H'] / 2 * hf[self.h2o] +
                 self.fuels[f]['C'] * hf[self.co2] -
                 ((self.fuels[f]['C'] + self.fuels[f]['H'] / 4) * hf[self.o2] +
                  hf[list(key)[0]])) /
               self.inl[0].fluid.wrapper[f]._molar_mass * 1000)

        return val

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
        return (self.inl[0].m.val_SI + self.inl[1].m.val_SI -
                self.outl[0].m.val_SI)

    def mass_flow_func_doc(self, label):
        r"""
        Calculate the residual value for component's mass flow balance.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        latex = (
            r'0=\dot{m}_\mathrm{in,1} + \dot{m}_\mathrm{in,2} - '
            r'\dot{m}_\mathrm{out,1}')
        return generate_latex_eq(self, latex, label)

    def mass_flow_deriv(self, k):
        r"""
        Calculate the partial derivatives for all mass flow balance equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the fluid equations.
        """
        for i in self.inl:
            if i.m.is_var:
                self.jacobian[k, i.m.J_col] = 1

        if self.outl[0].m.is_var:
            self.jacobian[k, self.outl[0].m.J_col] = -1

    def combustion_pressure_func(self):
        r"""
        Equations for reactor pressure balance.

        Returns
        -------
        residual : list
            Residual values of equations.

            .. math::

                0 = p_\mathrm{in,3} - p_\mathrm{out,3}\\
                0 = p_\mathrm{in,3} - p_\mathrm{in,4}
        """
        inl = self.inl[::-1][:2][::-1]
        outl = self.outl[::-1][0]

        return [inl[0].p.val_SI - outl.p.val_SI, inl[1].p.val_SI - outl.p.val_SI]

    def combustion_pressure_func_doc(self, label):
        r"""
        Equations for reactor pressure balance.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        inl = self.inl[::-1][:2][::-1]
        outl = self.outl[::-1][0]

        idx_out = outl.source_id[:3] + ',' + outl.source_id[3:]
        idx_in1 = inl[0].target_id[:2] + ',' + inl[0].target_id[2:]
        idx_in2 = inl[1].target_id[:2] + ',' + inl[1].target_id[2:]
        latex = (
            r'\begin{split}' + '\n'
            r'0 = & p_\mathrm{' + idx_in1 + r'} - p_\mathrm{' + idx_out +
            r'}\\' + '\n'
            r'0 = & p_\mathrm{' + idx_in1 + r'} - p_\mathrm{' + idx_in2 +
            r'}\\' + '\n'
            r'\end{split}')
        return generate_latex_eq(self, latex, label)

    def combustion_pressure_deriv(self, k):
        r"""
        Calculate the partial derivatives for combustion pressure equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the fluid equations.
        """
        if self.inl[0].p.is_var:
            self.jacobian[k,  self.inl[0].p.J_col] = 1
        if self.outl[0].p.is_var:
            self.jacobian[k,  self.outl[0].p.J_col] = -1

        if self.inl[1].p.is_var:
            self.jacobian[k + 1,  self.inl[1].p.J_col] = 1
        if self.outl[0].p.is_var:
            self.jacobian[k + 1,  self.outl[0].p.J_col] = -1

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
        inl = self.inl[::-1][:2][::-1]
        outl = self.outl[::-1][0]

        ###################################################################
        # molar mass flow for fuel and oxygen
        n_fuel = {}
        n_oxy_stoich = {}
        n_h = 0
        n_c = 0
        for f in self.fuel_list:
            n_fuel[f] = 0
            for i in inl:
                n = i.m.val_SI * i.fluid.val[f] / self.inl[0].fluid.wrapper[f]._molar_mass
                n_fuel[f] += n
                n_h += n * self.fuels[f]['H']
                n_c += n * self.fuels[f]['C']

            # stoichiometric oxygen requirement for each fuel
            n_oxy_stoich[f] = n_fuel[f] * (
                self.fuels[f]['H'] / 4 + self.fuels[f]['C'])

        n_oxygen = 0
        for i in inl:
            n_oxygen += (
                i.m.val_SI * i.fluid.val[self.o2] / self.inl[0].fluid.wrapper[self.o2]._molar_mass)

        ###################################################################
        # calculate stoichiometric oxygen
        n_oxygen_stoich = n_h / 4 + n_c

        ###################################################################
        # calculate lambda if not set
        if not self.lamb.is_set:
            self.lamb.val = n_oxygen / n_oxygen_stoich

        ###################################################################
        # calculate excess fuel if lambda is lower than 1
        if self.lamb.val < 1:
            n_h_exc = (n_oxygen_stoich - n_oxygen) * 4
            n_c_exc = (n_oxygen_stoich - n_oxygen)
        else:
            n_h_exc = 0
            n_c_exc = 0

        ###################################################################
        # equation for carbondioxide
        if fluid == self.co2:
            dm = (n_c - n_c_exc) * self.inl[0].fluid.wrapper[self.co2]._molar_mass

        ###################################################################
        # equation for water
        elif fluid == self.h2o:
            dm = (n_h - n_h_exc) / 2 * self.inl[0].fluid.wrapper[self.h2o]._molar_mass

        ###################################################################
        # equation for oxygen
        elif fluid == self.o2:
            if self.lamb.val < 1:
                dm = -n_oxygen * self.inl[0].fluid.wrapper[self.o2]._molar_mass
            else:
                dm = -n_oxygen / self.lamb.val * self.inl[0].fluid.wrapper[self.o2]._molar_mass

        ###################################################################
        # equation for fuel
        elif fluid in self.fuel_list:
            if self.lamb.val < 1:
                n_fuel_exc = (
                    -(n_oxygen / n_oxygen_stoich - 1) *
                    n_oxy_stoich[fluid] / (
                        self.fuels[fluid]['H'] / 4 +
                        self.fuels[fluid]['C']))
            else:
                n_fuel_exc = 0
            dm = -(n_fuel[fluid] - n_fuel_exc) * self.inl[0].fluid.wrapper[fluid]._molar_mass

        ###################################################################
        # equation for other fluids
        else:
            dm = 0

        res = dm
        for i in inl:
            res += i.fluid.val[fluid] * i.m.val_SI
        res -= outl.fluid.val[fluid] * outl.m.val_SI
        return res

    def stoichiometry_func_doc(self, label):
        r"""
        Generate stoichiometry LaTeX equations.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        idx_in1 = str(self.inl.index(self.inl[-2]) + 1)
        idx_in2 = str(self.inl.index(self.inl[-1]) + 1)
        idx_out = str(self.outl.index(self.outl[-1]) + 1)

        in1 = (
            r'\dot{m}_\mathrm{in,' + idx_in1 + r'} \cdot '
            r'x_\mathrm{fluid,in,' + idx_in1 + r'} ')
        in2 = (
            r'\dot{m}_\mathrm{in,' + idx_in2 + r'} \cdot '
            r'x_\mathrm{fluid,in,' + idx_in2 + r'}')
        out = (
            r'\dot{m}_\mathrm{out,' + idx_out + r'} \cdot '
            r'x_\mathrm{fluid,out,' + idx_out + r'}')

        equations = ''
        for fluid in self.inl[0].fluid.val:
            if fluid == self.o2:
                latex = (
                    r'0=\Delta\dot{m}_\mathrm{' + fluid + r'}-'
                    r'\dot{m}_\mathrm{' + self.o2 + r',m,stoich} \cdot '
                    r'M_\mathrm{' + self.o2 + r'}'
                )
            elif fluid == self.co2:
                latex = (
                    r'0=\Delta \dot{m}_\mathrm{' + fluid + r'} + '
                    r'\dot{m}_\mathrm{C,m} \cdot M_\mathrm{' +
                    self.co2 + r'} '
                )
            elif fluid == self.h2o:
                latex = (
                    r'0=\Delta \dot{m}_\mathrm{' + fluid + r'} + '
                    r'\frac{\dot{m}_\mathrm{H,m}}{2} \cdot M_\mathrm{' +
                    self.h2o + r'} '
                )
            elif fluid in self.fuel_list:
                latex = (
                    r'0=\Delta\dot{m}_\mathrm{' + fluid + r'}-'
                    r'\dot{m}_\mathrm{' + fluid + r',m} \cdot M_\mathrm{' +
                    fluid + r'}'
                )
            else:
                latex = r'0 = \Delta \dot{m}_\mathrm{' + fluid + '}'

            if fluid == next(iter(self.inl[0].fluid.val)):
                balance = (
                    r'\Delta \dot{m}_\mathrm{fluid} = ' + in1 +
                    '+' + in2 + '-' + out)
                m_fluid_molar = (
                    r'\dot{m}_\mathrm{fluid,m} = \frac{' + in1 + '+' +
                    in2 + r'}{M_\mathrm{fluid}}')
                m_o2_molar_stoich = (
                    r'\dot{m}_\mathrm{' + self.o2 + ',m,stoich}='
                    r'\frac{\dot{m}_\mathrm{H,m}}{4} + \dot{m}_\mathrm{C,m}')
                m_H_molar = r'\dot{m}_\mathrm{H,m}='
                m_C_molar = r'\dot{m}_\mathrm{C,m}='
                for f in self.fuel_list:
                    m_H_molar += (
                        r'\dot{m}_\mathrm{' + f + r',m} \cdot ' +
                        str(self.fuels[f]['H']) + '+')
                    m_C_molar += (
                        r'\dot{m}_\mathrm{' + f + r',m} \cdot ' +
                        str(self.fuels[f]['C']) + '+')
                m_H_molar = m_H_molar[:-1]
                m_C_molar = m_C_molar[:-1]
                latex_general_eq = (
                    r'\begin{split}' + '\n'
                    r'&' + balance + r'\\' + '\n'
                    r'&' + m_fluid_molar + r'\\' + '\n'
                    r'&' + m_H_molar + r'\\' + '\n'
                    r'&' + m_C_molar + r'\\' + '\n'
                    r'&' + m_o2_molar_stoich + r'\\' + '\n'
                    r'\end{split}'
                )
                equations += (
                    generate_latex_eq(
                        self, latex_general_eq, label + '_general_eq') + '\n' +
                    generate_latex_eq(self, latex, label + '_' + fluid) + '\n')
            else:
                equations += (
                    generate_latex_eq(self, latex, label + '_' + fluid) + '\n')
        # remove last newline
        return equations[:-1]

    def stoichiometry_deriv(self, increment_filter, k):
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
        inl = self.inl[::-1][:2][::-1]
        outl = self.outl[::-1][0]
        import itertools
        f = self.stoichiometry
        conns = inl + [outl]
        for fluid, conn in itertools.product(self.fluid_eqs_list, conns):
            eq_num = self.fluid_eqs_list.index(fluid)
            if self.is_variable(conn.m, increment_filter):
                self.jacobian[k + eq_num, conn.m.J_col] = self.numeric_deriv(
                    f, 'm', conn, fluid=fluid
                )
            for fluid_name in conn.fluid.is_var:
                self.jacobian[k + eq_num, conn.fluid.J_col[fluid_name]] = self.numeric_deriv(
                    f, fluid_name, conn, fluid=fluid
                )

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
        T_ref = 298.15
        p_ref = 1e5

        res = 0
        for i in self.inl:
            i.build_fluid_data()
            res += i.m.val_SI * (i.h.val_SI - h_mix_pT(
                p_ref, T_ref, i.fluid_data, mixing_rule="forced-gas"))

        for o in self.outl:
            o.build_fluid_data()
            res -= o.m.val_SI * (o.h.val_SI - h_mix_pT(
                p_ref, T_ref, o.fluid_data, mixing_rule="forced-gas"))

        res += self.calc_ti()
        return res

    def energy_balance_func_doc(self, label):
        r"""
        Calculate the energy balance of the adiabatic combustion chamber.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        latex = (
            r'\begin{split}' + '\n'
            r'0 = & \sum_i \dot{m}_{\mathrm{in,}i} \cdot\left( '
            r'h_{\mathrm{in,}i} - h_{\mathrm{in,}i\mathrm{,ref}} \right) -'
            r'\dot{m}_\mathrm{out,1}\cdot\left( h_\mathrm{out,1}'
            r' - h_\mathrm{out,1,ref}\right)\\' + '\n'
            r'& + LHV_{fuel} \cdot \left(\sum_i \dot{m}_{\mathrm{in,}i} '
            r'\cdot x_{fuel\mathrm{,in,}i} - \dot{m}_\mathrm{out,1} '
            r'\cdot x_{fuel\mathrm{,out,1}} \right)\\' + '\n'
            r'& \forall i \in \text{inlets}\\'
            r'& T_\mathrm{ref}=\unit[298.15]{K}'
            r'\;p_\mathrm{ref}=\unit[10^5]{Pa}\\'
            '\n' + r'\end{split}'
        )
        return generate_latex_eq(self, latex, label)

    def energy_balance_deriv(self, increment_filter, k):
        """
        Calculate partial derivatives of energy balance function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of equation in Jacobian matrix.
        """
        f = self.energy_balance_func
        for c in self.inl + self.outl:
            if self.is_variable(c.m, increment_filter):
                self.jacobian[k, c.m.J_col] = self.numeric_deriv(f, 'm', c)
            if self.is_variable(c.p, increment_filter):
                self.jacobian[k, c.p.J_col] = self.numeric_deriv(f, 'p', c)
            if self.is_variable(c.h, increment_filter):
                if c == self.outl[0]:
                    self.jacobian[k, c.h.J_col] = -c.m.val_SI
                else:
                    self.jacobian[k, c.h.J_col] = c.m.val_SI

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
        # required to work with combustion chamber and engine
        inl = self.inl[::-1][:2][::-1]
        n_h = 0
        n_c = 0
        for f in self.fuel_list:
            n_fuel = 0
            for i in inl:
                n_fuel += i.m.val_SI * i.fluid.val[f] / self.inl[0].fluid.wrapper[f]._molar_mass
                n_h += n_fuel * self.fuels[f]['H']
                n_c += n_fuel * self.fuels[f]['C']

        n_oxygen = 0
        for i in inl:
            n_oxygen += (i.m.val_SI * i.fluid.val[self.o2] /
                         self.inl[0].fluid.wrapper[self.o2]._molar_mass)

        n_oxygen_stoich = n_h / 4 + n_c

        return n_oxygen / n_oxygen_stoich - self.lamb.val

    def lambda_func_doc(self, label):
        r"""
        Calculate the residual for specified lambda.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        latex = (
            r'\begin{split}' + '\n'
            r'0 = \frac{\dot{m}_\mathrm{fuel,m}}{\dot{m}_\mathrm{O_2,m} '
            r'\cdot \left(n_\mathrm{C,fuel} + 0.25 \cdot n_\mathrm{H,fuel}'
            r'\right)} - \lambda \\' + '\n'
            r'\dot{m}_\mathrm{fluid,m} = \frac{x_\mathrm{fluid} \cdot '
            r'\dot{m}}{M_\mathrm{fluid}}\\' + '\n'
            r'\end{split}'
        )
        return generate_latex_eq(self, latex, label)

    def lambda_deriv(self, increment_filter, k):
        """
        Calculate partial derivatives of lambda function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of equation in Jacobian matrix.
        """
        # required to work with combustion chamber and engine
        inl = self.inl[::-1][:2][::-1]
        f = self.lambda_func
        for conn in inl:
            if self.is_variable(conn.m, increment_filter):
                self.jacobian[k, conn.m.J_col] = self.numeric_deriv(f, 'm', conn)
            for fluid in conn.fluid.is_var:
                self.jacobian[k, conn.fluid.J_col[fluid]] = self.numeric_deriv(f, fluid, conn)

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
        return self.ti.val - self.calc_ti()

    def ti_func_doc(self, label):
        r"""
        Calculate the residual for specified thermal input.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        idx = str(self.outl.index(self.outl[-1]) + 1)
        latex = (
            r'\begin{split}' + '\n'
            r'0 = & ti - LHV_\mathrm{fuel} \cdot \left[\sum_i \left('
            r'\dot{m}_{\mathrm{in,}i}\cdot x_{\mathrm{fuel,in,}i}\right)-'
            r' \dot{m}_\mathrm{out,' + idx + r'}\cdot '
            r'x_{\mathrm{fuel,out,}' + idx + r'} \right]\\' + '\n'
            r'& \forall i \in \text{combustion inlets}\\' + '\n'
            r'\end{split}'
        )
        return generate_latex_eq(self, latex, label)

    def ti_deriv(self, increment_filter, k):
        """
        Calculate partial derivatives of thermal input function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of equation in Jacobian matrix.
        """
        for i in self.inl:
            if i.m.is_var:
                deriv = 0
                for f in self.fuel_list:
                    deriv -= i.fluid.val[f] * self.fuels[f]['LHV']
                self.jacobian[k, i.m.J_col] = deriv
            for f in (self.fuel_list & i.fluid.is_var):
                    self.jacobian[k, i.fluid.J_col[f]] = -i.m.val_SI * self.fuels[f]['LHV']

        o = self.outl[0]
        if o.m.is_var:
            deriv = 0
            for f in self.fuel_list:
                deriv += o.fluid.val[f] * self.fuels[f]['LHV']
            self.jacobian[k, o.m.J_col] = deriv
        for f in (self.fuel_list & o.fluid.is_var):
                self.jacobian[k, o.fluid.J_col[f]] = o.m.val_SI * self.fuels[f]['LHV']

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
        ti = 0
        for f in self.fuel_list:
            m = 0
            for i in self.inl:
                m += i.m.val_SI * i.fluid.val[f]

            for o in self.outl:
                m -= o.m.val_SI * o.fluid.val[f]

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

    def bus_func_doc(self, bus):
        r"""
        Return LaTeX string of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus.Bus
            TESPy bus object.

        Returns
        -------
        latex : str
            LaTeX string of bus function.
        """
        idx = str(self.outl.index(self.outl[-1]) + 1)
        return (
            r'LHV_\mathrm{fuel} \cdot \left[\sum_i \left('
            r'\dot{m}_{\mathrm{in,}i}\cdot x_{\mathrm{fuel,in,}i}\right)-'
            r' \dot{m}_\mathrm{out,' + idx + r'}\cdot '
            r'x_{\mathrm{fuel,out,}' + idx + r'} \right]'
        )

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
        deriv = np.zeros((1, 3, self.num_nw_vars))
        f = self.calc_bus_value
        for i in range(3):
            deriv[0, i, 0] = self.numeric_deriv(f, 'm', i, bus=bus)
            deriv[0, i, 3:] = self.numeric_deriv(f, 'fluid', i, bus=bus)

        return deriv

    def initialise_fluids(self):
        """Calculate reaction balance for generic starting values at outlet."""
        N_2 = 0.7655
        O_2 = 0.2345

        n_fuel = 1
        lamb = 3

        fact_fuel = {}
        sum_fuel = 0
        for f in self.fuel_list:
            fact_fuel[f] = 0
            for i in self.inl:
                fact_fuel[f] += i.fluid.val[f] / 2
            sum_fuel += fact_fuel[f]

        for f in self.fuel_list:
            fact_fuel[f] /= sum_fuel

        m_co2 = 0
        m_h2o = 0
        m_fuel = 0
        for f in self.fuel_list:
            m_co2 += (n_fuel * self.fuels[f]['C'] * self.inl[0].fluid.wrapper[self.co2]._molar_mass *
                      fact_fuel[f])
            m_h2o += (n_fuel * self.fuels[f]['H'] / 2 *
                      self.inl[0].fluid.wrapper[self.h2o]._molar_mass * fact_fuel[f])
            m_fuel += n_fuel * self.inl[0].fluid.wrapper[f]._molar_mass * fact_fuel[f]

        n_o2 = (m_co2 / self.inl[0].fluid.wrapper[self.co2]._molar_mass +
                0.5 * m_h2o / self.inl[0].fluid.wrapper[self.h2o]._molar_mass) * lamb

        m_air = n_o2 * self.inl[0].fluid.wrapper[self.o2]._molar_mass / O_2
        m_fg = m_air + m_fuel

        m_o2 = n_o2 * self.inl[0].fluid.wrapper[self.o2]._molar_mass * (1 - 1 / lamb)
        m_n2 = N_2 * m_air

        fg = {
            self.n2: m_n2 / m_fg,
            self.co2: m_co2 / m_fg,
            self.o2: m_o2 / m_fg,
            self.h2o: m_h2o / m_fg
        }

        for o in self.outl:
            for fluid, x in o.fluid.val.items():
                if fluid in o.fluid.is_var and fluid in fg:
                    o.fluid.val[fluid] = fg[fluid]
            o.target.propagate_fluid_to_target(o, o.target)

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
        inl = self.inl[::-1][:2][::-1]
        outl = self.outl[::-1][0]

        m = 0
        for i in inl:
            if not i.good_starting_values:
                if i.m.val_SI < 0 and i.m.is_var:
                    i.m.val_SI = 0.01
                m += i.m.val_SI

        ######################################################################
        # check fluid composition
        if not outl.good_starting_values:
            for f in outl.fluid.is_var:
                if f == self.o2:
                    if outl.fluid.val[f] > 0.25:
                        outl.fluid.val[f] = 0.2
                    if outl.fluid.val[f] < 0.001:
                        outl.fluid.val[f] = 0.05

                elif f == self.co2:
                    if outl.fluid.val[f] > 0.1:
                        outl.fluid.val[f] = 0.075
                    if outl.fluid.val[f] < 0.001:
                        outl.fluid.val[f] = 0.02

                elif f == self.h2o:
                    if outl.fluid.val[f] > 0.1:
                        outl.fluid.val[f] = 0.075
                    if outl.fluid.val[f] < 0.001:
                        outl.fluid.val[f] = 0.02

                elif f in self.fuel_list:
                    if outl.fluid.val[f] > 0:
                        outl.fluid.val[f] = 0

                else:
                    m_f = 0
                    for i in inl:
                        m_f += i.fluid.val[f] * i.m.val_SI

                    if abs(outl.fluid.val[f] - m_f / m) > 0.03:
                        outl.fluid.val[f] = m_f / m

        ######################################################################
        # flue gas propagation
        if not outl.good_starting_values:
            if outl.m.val_SI < 0 and not outl.m.is_set:
                outl.m.val_SI = 10
            outl.target.propagate_fluid_to_target(outl, outl.target)

            if outl.h.val_SI < 7.5e5 and not outl.h.is_set:
                outl.h.val_SI = 1e6

        ######################################################################
        # additional checks for performance improvement
        if self.lamb.val < 2 and not self.lamb.is_set:
            # search fuel and air inlet
            for i in inl:
                fuel_found = False
                if not i.good_starting_values:
                    fuel = 0
                    for f in self.fuel_list:
                        fuel += i.fluid.val[f]
                    # found the fuel inlet
                    if fuel > 0.75 and not i.m.is_set:
                        fuel_found = True
                        fuel_inlet = i

                    # found the air inlet
                    if fuel < 0.75:
                        air_tmp = i.m.val_SI

            if fuel_found:
                fuel_inlet.m.val_SI = air_tmp / 25

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

    def propagate_fluid_to_target(self, inconn, start, entry_point=False):
        r"""
        Fluid propagation to target stops here.

        Parameters
        ----------
        inconn : tespy.connections.connection.Connection
            Connection to initialise.

        start : tespy.components.component.Component
            This component is the fluid propagation starting point.
            The starting component is saved to prevent infinite looping.
        """
        return

    def propagate_fluid_to_source(self, outconn, start, entry_point=False):
        r"""
        Fluid propagation to source stops here.

        Parameters
        ----------
        outconn : tespy.connections.connection.Connection
            Connection to initialise.

        start : tespy.components.component.Component
            This component is the fluid propagation starting point.
            The starting component is saved to prevent infinite looping.
        """
        return

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        self.ti.val = self.calc_ti()

        n_h = 0
        n_c = 0
        for f in self.fuel_list:
            n_fuel = 0
            for i in self.inl:
                n_fuel += i.m.val_SI * i.fluid.val[f] / self.inl[0].fluid.wrapper[f]._molar_mass
                n_h += n_fuel * self.fuels[f]['H']
                n_c += n_fuel * self.fuels[f]['C']

        n_oxygen = 0
        for i in self.inl:
            n_oxygen += (
                i.m.val_SI * i.fluid.val[self.o2]
                / self.inl[0].fluid.wrapper[self.o2]._molar_mass
            )

        n_oxygen_stoich = n_h / 4 + n_c

        self.lamb.val = n_oxygen / n_oxygen_stoich

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
            o.s.val_SI -
            s_mix_pT([0, p_ref, 0, o.fluid.val], T_ref, force_gas=True))
        for c in self.inl:
            self.S_comb -= c.m.val_SI * (
                c.s.val_SI -
                s_mix_pT([0, p_ref, 0, c.fluid.val], T_ref, force_gas=True))

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
        self.epsilon = self.E_P/self.E_F
        self.E_bus = np.nan
