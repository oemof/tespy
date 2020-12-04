# -*- coding: utf-8

"""Module of class CombustionChamber.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/combustion/combustion_chamber.py

SPDX-License-Identifier: MIT
"""

import logging

import CoolProp.CoolProp as CP
import numpy as np

from tespy.components.component import Component
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import DataContainerSimple as dc_simple
from tespy.tools.fluid_properties import h_mix_pQ
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.fluid_properties import h_pT
from tespy.tools.global_vars import err
from tespy.tools.global_vars import molar_masses
from tespy.tools.helpers import TESPyComponentError
from tespy.tools.helpers import fluid_structure
from tespy.tools.helpers import molar_mass_flow


class CombustionChamber(Component):
    r"""
    The class CombustionChamber is parent class of all combustion components.

    Equations

        **mandatory equations**

        - :py:meth:`tespy.components.combustion.combustion_chamber.CombustionChamber.reaction_balance`
        - :py:meth:`tespy.components.component.Component.mass_flow_func`

        .. math::

            0 = p_{in,i} - p_{out} \;
            \forall i \in \mathrm{inlets}

        - :py:meth:`tespy.components.combustion.combustion_chamber.CombustionChamber.energy_balance`

        **optional equations**

        - :py:meth:`tespy.components.combustion.combustion_chamber.CombustionChamber.lambda_func`
        - :py:meth:`tespy.components.combustion.combustion_chamber.CombustionChamber.ti_func`

    Available fuels

        - methane, ethane, propane, butane, hydrogen

    Inlets/Outlets

        - in1, in2
        - out1

    Image

        .. image:: _images/CombustionChamber.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

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

    design_path: str
        Path to the components design case.

    local_offdesign : boolean
        Treat this component in offdesign mode in a design calculation.

    local_design : boolean
        Treat this component in design mode in an offdesign calculation.

    char_warnings: boolean
        Ignore warnings on default characteristics usage for this component.

    printout: boolean
        Include this component in the network's results printout.

    lamb : float, tespy.tools.data_containers.ComponentProperties
        Actual oxygen to stoichiometric oxygen ratio, :math:`\lambda/1`.

    ti : float, tespy.tools.data_containers.ComponentProperties
        Thermal input, (:math:`{LHV \cdot \dot{m}_f}`),
        :math:`ti/\text{W}`.

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
    >>> from tespy.tools.fluid_properties import T_bp_p
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
    2.017
    >>> comb.set_attr(lamb=2)
    >>> comb_fg.set_attr(T=np.nan)
    >>> nw.solve('design')
    >>> round(comb_fg.T.val, 1)
    1208.4
    """

    @staticmethod
    def component():
        return 'combustion chamber'

    @staticmethod
    def attr():
        return {
            'lamb': dc_cp(min_val=1),
            'ti': dc_cp(min_val=0),
            'S': dc_simple()
        }

    @staticmethod
    def inlets():
        return ['in1', 'in2']

    @staticmethod
    def outlets():
        return ['out1']

    def comp_init(self, nw):

        Component.comp_init(self, nw)

        # number of mandatroy equations for
        # reaction balance: num_fl
        # mass flow: 1
        # pressure: 2
        # energy balance: 1
        self.num_eq = self.num_nw_fluids + 4
        for var in [self.lamb, self.ti]:
            if var.is_set is True:
                self.num_eq += 1

        self.jacobian = np.zeros((
            self.num_eq,
            self.num_i + self.num_o + self.num_vars,
            self.num_nw_vars))

        self.residual = np.zeros(self.num_eq)
        self.jacobian[0:1] = self.mass_flow_deriv()
        self.jacobian[1:3] = self.pressure_deriv()

        self.setup_reaction_parameters()

    def setup_reaction_parameters(self):
        r"""Setup parameters for reaction (gas name aliases and LHV)."""
        self.fuel_list = []
        fuels = ['methane', 'ethane', 'propane', 'butane', 'hydrogen']
        for f in fuels:
            self.fuel_list += [x for x in self.nw_fluids if x in [
                a.replace(' ', '') for a in CP.get_aliases(f)]]

        if len(self.fuel_list) == 0:
            msg = ('Your network\'s fluids do not contain any fuels, that are '
                   'available for the component ' + self.label + ' of type ' +
                   self.component() + '. Available fuels are: ' + str(fuels) +
                   '.')
            logging.error(msg)
            raise TESPyComponentError(msg)

        else:
            msg = ('The fuels for component ' + self.label + ' of type ' +
                   self.component() + ' are: ' + str(self.fuel_list) + '.')
            logging.debug(msg)

        self.o2 = [x for x in self.nw_fluids if x in [
            a.replace(' ', '') for a in CP.get_aliases('O2')]][0]
        self.co2 = [x for x in self.nw_fluids if x in [
            a.replace(' ', '') for a in CP.get_aliases('CO2')]][0]
        self.h2o = [x for x in self.nw_fluids if x in [
            a.replace(' ', '') for a in CP.get_aliases('H2O')]][0]
        self.n2 = [x for x in self.nw_fluids if x in [
            a.replace(' ', '') for a in CP.get_aliases('N2')]][0]

        self.fuels = {}
        for f in self.fuel_list:
            self.fuels[f] = {}
            structure = fluid_structure(f)
            for el in ['C', 'H', 'O']:
                if el in structure.keys():
                    self.fuels[f][el] = structure[el]
                else:
                    self.fuels[f][el] = 0
            self.fuels[f]['LHV'] = self.calc_lhv(f)

    def calc_lhv(self, f):
        r"""
        Calculate the lower heating value of the combustion chamber's fuel.

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
        hf['methane'] = -74.85
        hf['ethane'] = -84.68
        hf['propane'] = -103.8
        hf['butane'] = -124.51
        hf[self.o2] = 0
        hf[self.co2] = -393.5
        # water (gaseous)
        hf[self.h2o] = -241.8

        key = set(list(hf.keys())).intersection(
                set([a.replace(' ', '')
                     for a in CP.get_aliases(f)]))

        val = (-(self.fuels[f]['H'] / 2 * hf[self.h2o] +
                 self.fuels[f]['C'] * hf[self.co2] -
                 ((self.fuels[f]['C'] + self.fuels[f]['H'] / 4) * hf[self.o2] +
                  hf[list(key)[0]])) /
               molar_masses[f] * 1000)

        return val

    def equations(self):
        r"""Calculate residual vector with results of equations."""
        k = 0
        ######################################################################
        # eqation for mass flow balance
        self.residual[k] = self.mass_flow_func()
        k += 1

        ######################################################################
        # equations for pressure
        for i in self.inl:
            self.residual[k] = self.outl[0].p.val_SI - i.p.val_SI
            k += 1

        ######################################################################
        # equations for fluids in reaction balance
        for fluid in self.inl[0].fluid.val.keys():
            if (np.absolute(self.residual[k]) > err ** 2 or self.it % 3 == 0 or
                    self.always_all_equations):
                self.residual[k] = self.reaction_balance(fluid)
            k += 1

        ######################################################################
        # equation for energy balance
        if (np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0 or
                self.always_all_equations):
            self.residual[k] = self.energy_balance()
        k += 1

        ######################################################################
        # equation for specified air to stoichiometric air ratio lamb
        if self.lamb.is_set:
            self.residual[k] = self.lambda_func()
            k += 1

        ######################################################################
        # equation for speciified thermal input
        if self.ti.is_set:
            self.residual[k] = self.ti_func()
            k += 1

    def derivatives(self, increment_filter):
        r"""Calculate matrix of partial derivatives for given equations."""
        ######################################################################
        # derivatives for mass flow and pressure are static
        k = 3

        ######################################################################
        # derivatives for reaction balance
        for fluid in self.nw_fluids:
            for i in range(3):
                if not increment_filter[i, 0]:
                    self.jacobian[k, i, 0] = self.rb_numeric_deriv(
                        'm', i, fluid)
                if not all(increment_filter[i, 3:]):
                    self.jacobian[k, i, 3:] = self.rb_numeric_deriv(
                        'fluid', i, fluid)
            k += 1

        ######################################################################
        # derivatives for energy balance equations
        f = self.energy_balance
        for i in range(3):
            if not increment_filter[i, 0]:
                self.jacobian[k, i, 0] = self.numeric_deriv(f, 'm', i)
            if not increment_filter[i, 1]:
                self.jacobian[k, i, 1] = self.numeric_deriv(f, 'p', i)
            if i >= self.num_i:
                self.jacobian[k, i, 2] = -(
                    self.inl + self.outl)[i].m.val_SI
            else:
                self.jacobian[k, i, 2] = (
                    self.inl + self.outl)[i].m.val_SI
        k += 1

        ######################################################################
        # derivatives for specified lamb
        if self.lamb.is_set:
            f = self.lambda_func
            if not increment_filter[0, 0]:
                self.jacobian[k, 0, 0] = self.numeric_deriv(f, 'm', 0)
            if not all(increment_filter[0, 3:]):
                self.jacobian[k, 0, 3:] = self.numeric_deriv(f, 'fluid', 0)
            if not increment_filter[1, 0]:
                self.jacobian[k, 1, 0] = self.numeric_deriv(f, 'm', 1)
            if not all(increment_filter[1, 3:]):
                self.jacobian[k, 1, 3:] = self.numeric_deriv(f, 'fluid', 1)
            k += 1

        ######################################################################
        # derivatives for specified thermal input
        if self.ti.is_set:
            self.jacobian[k, 0, 0] = 0
            self.jacobian[k, 1, 0] = 0
            self.jacobian[k, 2, 0] = 0
            for f in self.fuel_list:
                pos = 3 + self.nw_fluids.index(f)
                lhv = self.fuels[f]['LHV']

                for i in range(2):
                    self.jacobian[k, i, 0] += -(
                        self.inl[i].fluid.val[f] * lhv)
                    self.jacobian[k, i, pos] = -(
                        self.inl[i].m.val_SI * lhv)
                self.jacobian[k, 2, 0] += (
                    self.outl[0].fluid.val[f] * lhv)
                self.jacobian[k, 2, pos] = (
                    self.outl[0].m.val_SI * lhv)
            k += 1

    def pressure_deriv(self):
        r"""
        Calculate the partial derivatives for all pressure equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((2, 3, self.num_nw_vars))
        for k in range(2):
            deriv[k][2][1] = 1
            deriv[k][k][1] = -1
        return deriv

    def reaction_balance(self, fluid):
        r"""
        Calculate the reaction balance for one fluid.

        - determine molar mass flows of fuel and oxygen
        - calculate mole number of carbon and hydrogen atoms in fuel
        - calculate molar oxygen flow for stoichiometric combustion
        - calculate residual value of the fluids balance

        for excess fuel

        - calculate excess carbon and hydrogen in fuels
        - calculate excess fuel shares

        General equations

        .. math::

            i \in [1,2], o \in [1]\\

            res = \sum_i \left(x_{fluid,i} \cdot \dot{m}_{i}\right) -
            \sum_j \left(x_{fluid,j} \cdot \dot{m}_{j}\right) \;
            \forall i, \; \forall j

            \dot{m}_{fluid,m} = \sum_i \frac{x_{fluid,i} \cdot \dot{m}_{i}}
            {M_{fluid}} \; \forall i

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

        Equation for fuel

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
            The fluid to calculate the reation balance for.

        Returns
        -------
        res : float
            Residual value of equation.
        """
        # looks strange, required to work with combustion chamber and engine
        inl = self.inl[::-1][:2][::-1]
        outl = self.outl[::-1][0]

        ######################################################################
        # molar mass flow for fuel and oxygen
        n_fuel = {}
        n_oxy_stoich = {}
        n_h = 0
        n_c = 0
        for f in self.fuel_list:
            n_fuel[f] = 0
            for i in inl:
                n = i.m.val_SI * i.fluid.val[f] / molar_masses[f]
                n_fuel[f] += n
                n_h += n * self.fuels[f]['H']
                n_c += n * self.fuels[f]['C']

            # stoichiometric oxygen requirement for each fuel
            n_oxy_stoich[f] = n_fuel[f] * (self.fuels[f]['H'] / 4 +
                                           self.fuels[f]['C'])

        n_oxygen = 0
        for i in inl:
            n_oxygen += (i.m.val_SI * i.fluid.val[self.o2] /
                         molar_masses[self.o2])

        ######################################################################
        # calculate stoichiometric oxygen
        n_oxygen_stoich = n_h / 4 + n_c

        ######################################################################
        # calculate lambda if not set
        if not self.lamb.is_set:
            self.lamb.val = n_oxygen / n_oxygen_stoich

        ######################################################################
        # calculate excess fuel if lambda is lower than 1
        if self.lamb.val < 1:
            n_h_exc = (n_oxygen_stoich - n_oxygen) * 4
            n_c_exc = (n_oxygen_stoich - n_oxygen)
        else:
            n_h_exc = 0
            n_c_exc = 0

        ######################################################################
        # equation for carbondioxide
        if fluid == self.co2:
            dm = (n_c - n_c_exc) * molar_masses[self.co2]

        ######################################################################
        # equation for water
        elif fluid == self.h2o:
            dm = (n_h - n_h_exc) / 2 * molar_masses[self.h2o]

        ######################################################################
        # equation for oxygen
        elif fluid == self.o2:
            if self.lamb.val < 1:
                dm = -n_oxygen * molar_masses[self.o2]
            else:
                dm = -n_oxygen / self.lamb.val * molar_masses[self.o2]

        ######################################################################
        # equation for fuel
        elif fluid in self.fuel_list:
            if self.lamb.val < 1:
                n_fuel_exc = (-(n_oxygen / n_oxygen_stoich - 1) *
                              n_oxy_stoich[fluid] *
                              (self.fuels[f]['H'] / 4 + self.fuels[f]['C']))
            else:
                n_fuel_exc = 0
            dm = -(n_fuel[fluid] - n_fuel_exc) * molar_masses[fluid]

        ######################################################################
        # equation for other fluids
        else:
            dm = 0

        res = dm
        for i in inl:
            res += i.fluid.val[fluid] * i.m.val_SI
        res -= outl.fluid.val[fluid] * outl.m.val_SI
        return res

    def rb_numeric_deriv(self, dx, pos, fluid):
        r"""
        Calculate derivative of the reaction balance to dx.

        Parameters
        ----------
        dx : str
            Partial derivative.

        pos : int
            Position of connection regarding to inlets and outlet of the
            component, logic: ['in1', 'in2', ..., 'out1', ...] ->
            0, 1, ..., n, n + 1, ..., n + m

        fluid : str
            Fluid to calculate partial derivative of reaction balance for.

        Returns
        -------
        deriv : float, list
            Partial derivative(s) of the function :math:`f` to variable(s)
            :math:`x`.

            .. math::

                \frac{\partial f}{\partial x} = \frac{f(x + d) + f(x - d)}{2 d}
        """
        dm, dp, dh, df = 0, 0, 0, 0
        if dx == 'm':
            dm = 1e-4
        else:
            df = 1e-5

        if dx == 'fluid':
            deriv = []
            for f in self.inl[0].fluid.val.keys():
                val = (self.inl + self.outl)[pos].fluid.val[f]
                exp = 0
                if (self.inl + self.outl)[pos].fluid.val[f] + df <= 1:
                    (self.inl + self.outl)[pos].fluid.val[f] += df
                else:
                    (self.inl + self.outl)[pos].fluid.val[f] = 1
                exp += self.reaction_balance(fluid)
                if (self.inl + self.outl)[pos].fluid.val[f] - 2 * df >= 0:
                    (self.inl + self.outl)[pos].fluid.val[f] -= 2 * df
                else:
                    (self.inl + self.outl)[pos].fluid.val[f] = 0
                exp -= self.reaction_balance(fluid)
                (self.inl + self.outl)[pos].fluid.val[f] = val

                deriv += [exp / (2 * (dm + dp + dh + df))]

        else:
            exp = 0
            (self.inl + self.outl)[pos].m.val_SI += dm
            (self.inl + self.outl)[pos].p.val_SI += dp
            (self.inl + self.outl)[pos].h.val_SI += dh
            exp += self.reaction_balance(fluid)

            (self.inl + self.outl)[pos].m.val_SI -= 2 * dm
            (self.inl + self.outl)[pos].p.val_SI -= 2 * dp
            (self.inl + self.outl)[pos].h.val_SI -= 2 * dh
            exp -= self.reaction_balance(fluid)
            deriv = exp / (2 * (dm + dp + dh + df))

            (self.inl + self.outl)[pos].m.val_SI += dm
            (self.inl + self.outl)[pos].p.val_SI += dp
            (self.inl + self.outl)[pos].h.val_SI += dh

        return deriv

    def energy_balance(self):
        r"""
        Calculate the energy balance of the adiabatic combustion chamber.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                \begin{split}
                res = & \sum_i \dot{m}_{in,i} \cdot
                \left( h_{in,i} - h_{in,i,ref} \right)\\
                & - \sum_j \dot{m}_{out,j} \cdot
                \left( h_{out,j} - h_{out,j,ref} \right)\\
                & + H_{I,f} \cdot \left(\sum_i \dot{m}_{in,i} \cdot x_{f,i} -
                \sum_j \dot{m}_{out,j} \cdot x_{f,j} \right)
                \end{split}\\

                \forall i \in \text{inlets}\; \forall j \in \text{outlets}\\
                x_{f}\text{: mass fraction of fuel}

        Note
        ----
        The temperature for the reference state is set to 20 °C, thus
        the water may be liquid. In order to make sure, the state is
        referring to the lower heating value, the necessary enthalpy
        difference for evaporation is added. The stoichiometric combustion
        chamber uses a different reference, you will find it in the
        :py:meth:`tespy.components.combustion.combustion_chamber_stoich.CombustionChamberStoich.energy_balance`
        documentation.

        - Reference temperature: 293.15 K.
        - Reference pressure: 1 bar.
        """
        T_ref = 293.15
        p_ref = 1e5

        res = 0
        for i in self.inl:
            res += i.m.val_SI * (
                i.h.val_SI - h_mix_pT([0, p_ref, 0, i.fluid.val], T_ref))

        for o in self.outl:
            dh = 0
            n_h2o = o.fluid.val[self.h2o] / molar_masses[self.h2o]
            if n_h2o > 0:
                p = p_ref * n_h2o / molar_mass_flow(o.fluid.val)
                h = h_pT(p, T_ref, self.h2o)
                try:
                    flow = [0, p, 0, {self.h2o: 1}]
                    h_steam = h_mix_pQ(flow, 1)
                except ValueError:
                    flow = [0, 615, 0, {self.h2o: 1}]
                    h_steam = h_mix_pQ(flow, 1)
                if h < h_steam:
                    dh = (h_steam - h) * o.fluid.val[self.h2o]

            res -= o.m.val_SI * (
                o.h.val_SI - h_mix_pT([0, p_ref, 0, o.fluid.val], T_ref) - dh)

        res += self.calc_ti()

        return res

    def lambda_func(self):
        r"""
        Calculate the residual for specified lambda.

        Returns
        -------
        val : float
            Residual value of function.

            .. math::

                \dot{m}_{fluid,m} = \sum_i \frac{x_{fluid,i} \cdot \dot{m}_{i}}
                {M_{fluid}}\\ \forall i \in inlets

                val = \frac{\dot{m}_{f,m}}{\dot{m}_{O_2,m} \cdot
                \left(n_{C,fuel} + 0.25 \cdot n_{H,fuel}\right)} - \lambda
        """
        inl = self.inl[::-1][:2][::-1]

        n_h = 0
        n_c = 0
        for f in self.fuel_list:
            n_fuel = 0
            for i in inl:
                n_fuel += i.m.val_SI * i.fluid.val[f] / molar_masses[f]
                n_h += n_fuel * self.fuels[f]['H']
                n_c += n_fuel * self.fuels[f]['C']

        n_oxygen = 0
        for i in inl:
            n_oxygen += (i.m.val_SI * i.fluid.val[self.o2] /
                         molar_masses[self.o2])

        n_oxygen_stoich = n_h / 4 + n_c

        return n_oxygen / n_oxygen_stoich - self.lamb.val

    def ti_func(self):
        r"""
        Calculate the residual for specified thermal input.

        Returns
        -------
        val : float
            Residual value of function.

            .. math::

                val = ti - \dot{m}_f \cdot LHV
        """
        return self.ti.val - self.calc_ti()

    def calc_ti(self):
        r"""
        Calculate the thermal input of the combustion chamber.

        Returns
        -------
        ti : float
            Thermal input.

            .. math::

                ti = LHV \cdot \left[\sum_i \left(\dot{m}_{in,i} \cdot x_{f,i}
                \right) - \dot{m}_{out,1} \cdot x_{f,1} \right]
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
            m_co2 += (n_fuel * self.fuels[f]['C'] * molar_masses[self.co2] *
                      fact_fuel[f])
            m_h2o += (n_fuel * self.fuels[f]['H'] / 2 *
                      molar_masses[self.h2o] * fact_fuel[f])
            m_fuel += n_fuel * molar_masses[f] * fact_fuel[f]

        n_o2 = (m_co2 / molar_masses[self.co2] +
                0.5 * m_h2o / molar_masses[self.h2o]) * lamb

        m_air = n_o2 * molar_masses[self.o2] / O_2
        m_fg = m_air + m_fuel

        m_o2 = n_o2 * molar_masses[self.o2] * (1 - 1 / lamb)
        m_n2 = N_2 * m_air

        fg = {
            self.n2: m_n2 / m_fg,
            self.co2: m_co2 / m_fg,
            self.o2: m_o2 / m_fg,
            self.h2o: m_h2o / m_fg
        }

        for o in self.outl:
            for fluid, x in o.fluid.val.items():
                if not o.fluid.val_set[fluid] and fluid in fg.keys():
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
        # looks super strange, necessary to work with engine and chamber
        inl = self.inl[::-1][:2][::-1]
        outl = self.outl[::-1][0]

        m = 0
        for i in inl:
            if i.good_starting_values is False:
                if i.m.val_SI < 0 and not i.m.val_set:
                    i.m.val_SI = 0.01
                m += i.m.val_SI

        ######################################################################
        # check fluid composition
        if outl.good_starting_values is False:
            fluids = [f for f in outl.fluid.val.keys()
                      if not outl.fluid.val_set[f]]
            for f in fluids:
                if f not in [self.o2, self.co2, self.h2o] + self.fuel_list:
                    m_f = 0
                    for i in inl:
                        m_f += i.fluid.val[f] * i.m.val_SI

                    if abs(outl.fluid.val[f] - m_f / m) > 0.03:
                        outl.fluid.val[f] = m_f / m

                elif f == self.o2:
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

        ######################################################################
        # flue gas propagation
        if outl.good_starting_values is False:
            if outl.m.val_SI < 0 and not outl.m.val_set:
                outl.m.val_SI = 10
            outl.target.propagate_fluid_to_target(outl, outl.target)

            if outl.h.val_SI < 7.5e5 and not outl.h.val_set:
                outl.h.val_SI = 1e6

        ######################################################################
        # additional checks for performance improvement
        if self.lamb.val < 2 and not self.lamb.is_set:
            # search fuel and air inlet
            for i in inl:
                fuel_found = False
                if i.good_starting_values is False:
                    fuel = 0
                    for f in self.fuel_list:
                        fuel += i.fluid.val[f]
                    # found the fuel inlet
                    if fuel > 0.75 and not i.m.val_set:
                        fuel_found = True
                        fuel_inlet = i

                    # found the air inlet
                    if fuel < 0.75:
                        air_tmp = i.m.val_SI

            if fuel_found is True:
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

    def propagate_fluid_to_target(self, inconn, start):
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

    def propagate_fluid_to_source(self, outconn, start):
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
                n_fuel += i.m.val_SI * i.fluid.val[f] / molar_masses[f]
                n_h += n_fuel * self.fuels[f]['H']
                n_c += n_fuel * self.fuels[f]['C']

        n_oxygen = 0
        for i in self.inl:
            n_oxygen += (i.m.val_SI * i.fluid.val[self.o2] /
                         molar_masses[self.o2])

        n_oxygen_stoich = n_h / 4 + n_c

        self.lamb.val = n_oxygen / n_oxygen_stoich

        self.check_parameter_bounds()
