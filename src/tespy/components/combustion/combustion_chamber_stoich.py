# -*- coding: utf-8

"""Module of class CombustionChamberStoich.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/combustion/combustion_chamber_stoich.py

SPDX-License-Identifier: MIT
"""

import logging

import CoolProp.CoolProp as CP

from tespy.components.combustion.combustion_chamber import CombustionChamber
from tespy.components.component import Component
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import DataContainerSimple as dc_simple
from tespy.tools.document_models import generate_latex_eq
from tespy.tools.fluid_properties import TESPyFluid
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.global_vars import molar_masses
from tespy.tools.helpers import TESPyComponentError
from tespy.tools.helpers import fluid_structure
from tespy.tools.helpers import molar_mass_flow


class CombustionChamberStoich(CombustionChamber):
    r"""
    The class CombustionChamberStoich is a simplified combustion chamber.

    **Mandatory Equations**

    - :py:meth:`tespy.components.combustion.combustion_chamber.CombustionChamber.mass_flow_func`
    - :py:meth:`tespy.components.combustion.combustion_chamber.CombustionChamber.combustion_pressure_func`
    - :py:meth:`tespy.components.combustion.combustion_chamber_stoich.CombustionChamberStoich.stoichiometry`
    - :py:meth:`tespy.components.combustion.combustion_chamber_stoich.CombustionChamberStoich.energy_balance_func`

    **Optional Equations**

    - :py:meth:`tespy.components.combustion.combustion_chamber_stoich.CombustionChamberStoich.lambda_func`
    - :py:meth:`tespy.components.combustion.combustion_chamber_stoich.CombustionChamberStoich.ti_func`

    Available fuels

    - methane, ethane, propane, butane, hydrogen

    Inlets/Outlets

    - in1, in2
    - out1

    Image

    .. image:: _images/CombustionChamber.svg
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

    fuel : dict
        Fuel composition, e.g. :code:`{'CH4': 0.96, 'CO2': 0.04}`.

    fuel_alias : str
        Alias for the fuel.

    air : dict
        Fresh air composition,
        e.g. :code:`{'N2': 0.76, 'O2': 0.23, 'Ar': 0.01}`.

    air_alias : str
        Alias for the fresh air.

    path : str
        Path to existing fluid property table.

    lamb : float, dict
        Air to stoichiometric air ratio, :math:`\lambda/1`.

    ti : float, dict
        Thermal input, (:math:`{LHV \cdot \dot{m}_f}`), :math:`ti/\text{W}`.

    Note
    ----
    This combustion chamber uses fresh air and its fuel as the only
    reactive gas components. Therefore note the following restrictions. You
    are to

    - specify the fluid composition of the fresh air,
    - fully define the fuel's fluid components,
    - provide the aliases of the fresh air and the fuel and
    - make sure, both of the aliases are part of the network fluid vector.

    If you choose 'Air' or 'air' as alias for the fresh air, TESPy will use
    the fluid properties from CoolProp's air. Else, a custom fluid
    'yourairalias' will be created.

    The name of the flue gas will be: 'yourfuelalias_fg'. It is also
    possible to use fluid mixtures for the fuel, e.g.
    :code:`fuel={CH4: 0.9, 'CO2': 0.1}`. If you specify a fluid mixture for
    the fuel, TESPy will automatically create a custom fluid called. For more
    information see the examples section or look for the combustion chamber
    tutorials at tespy.readthedocs.io.

    Example
    -------
    The stoichiometric combustion chamber follows identical physical properties
    as the combustion chamber. The main difference is, that the fuel and the
    air are not stated component wise but are fixed mixtures.
    The main advantage of using the stoichimetric combustion chamber
    comes from a strong improvement in terms of calculation speed.
    This example will show the same calculation as presented in the combustion
    chamber example (see
    :py:meth:`tespy.components.combustion.combustion_chamber.CombustionChamber`
    ).

    >>> from tespy.components import Sink, Source, CombustionChamberStoich
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> from tespy.tools.fluid_properties import T_bp_p
    >>> import numpy as np
    >>> import shutil
    >>> fluid_list = ['myAir', 'myFuel', 'myFuel_fg']
    >>> nw = Network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ... p_range=[1, 10], iterinfo=False)
    >>> amb = Source('ambient air')
    >>> sf = Source('fuel')
    >>> fg = Sink('flue gas outlet')
    >>> comb = CombustionChamberStoich('stoichiometric combustion chamber')
    >>> comb.component()
    'combustion chamber stoichiometric flue gas'
    >>> amb_comb = Connection(amb, 'out1', comb, 'in1')
    >>> sf_comb = Connection(sf, 'out1', comb, 'in2')
    >>> comb_fg = Connection(comb, 'out1', fg, 'in1')
    >>> nw.add_conns(sf_comb, amb_comb, comb_fg)

    Specify the thermal input of the combustion chamber. At the given fluid
    compositions this determines the mass flow of the fuel. The outlet
    temperature of the flue gas determines the ratio of oxygen to fuel mass
    flow. The fluid composition of the fuel and the air are defined, too. The
    results show very small deviation from the actual combustion chamber.

    >>> comb.set_attr(fuel={'CH4': 0.96, 'CO2': 0.03, 'H2': 0.01},
    ... air={'Ar': 0.0129, 'N2': 0.7553, 'H2O': 0, 'CH4': 0, 'CO2': 0.0004,
    ... 'O2': 0.2314}, fuel_alias='myFuel', air_alias='myAir', ti=500000)
    >>> amb_comb.set_attr(T=20, p=1, fluid={'myAir': 1,
    ... 'myFuel': 0,'myFuel_fg': 0})
    >>> sf_comb.set_attr(T=25, fluid={'myAir': 0, 'myFuel': 1,
    ... 'myFuel_fg': 0})
    >>> comb_fg.set_attr(T=1200)
    >>> nw.solve('design')
    >>> round(comb.lamb.val, 3)
    2.015
    >>> comb.set_attr(lamb=2)
    >>> comb_fg.set_attr(T=np.nan)
    >>> nw.solve('design')
    >>> round(comb_fg.T.val, 1)
    1207.1
    >>> shutil.rmtree('./LUT', ignore_errors=True)
    """

    @staticmethod
    def component():
        return 'combustion chamber stoichiometric flue gas'

    def get_variables(self):
        return {
            'fuel': dc_simple(), 'fuel_alias': dc_simple(),
            'air': dc_simple(), 'air_alias': dc_simple(),
            'path': dc_simple(),
            'lamb': dc_cp(
                min_val=1, deriv=self.lambda_deriv, func=self.lambda_func,
                latex=self.lambda_func_doc, num_eq=1),
            'ti': dc_cp(
                min_val=0, deriv=self.ti_deriv, func=self.ti_func,
                latex=self.ti_func_doc, num_eq=1)
        }

    @staticmethod
    def inlets():
        return ['in1', 'in2']

    @staticmethod
    def outlets():
        return ['out1']

    @staticmethod
    def fuels():
        return ['methane', 'ethane', 'propane', 'butane', 'hydrogen']

    def comp_init(self, nw):
        if not self.fuel.is_set or not isinstance(self.fuel.val, dict):
            msg = ('You must specify the fuel composition for stoichimetric '
                   'combustion chamber ' + self.label + '.')
            logging.error(msg)
            raise TESPyComponentError(msg)

        if not self.fuel_alias.is_set:
            msg = ('You must specify a fuel alias for stoichimetric '
                   'combustion chamber ' + self.label + '.')
            logging.error(msg)
            raise TESPyComponentError(msg)

        if not self.air.is_set or not isinstance(self.air.val, dict):
            msg = ('You must specify the air composition for stoichimetric '
                   'combustion chamber ' + self.label + '.')
            logging.error(msg)
            raise TESPyComponentError(msg)

        if not self.air_alias.is_set:
            msg = ('You must specify an air alias for stoichimetric '
                   'combustion chamber ' + self.label + '.')
            logging.error(msg)
            raise TESPyComponentError(msg)

        Component.comp_init(self, nw)

        # adjust the names for required fluids according to naming in the
        # network air
        for f in list(self.air.val.keys()):
            alias = [x for x in self.nw_fluids if x in [
                a.replace(' ', '') for a in CP.get_aliases(f)]]
            if len(alias) > 0:
                self.air.val[alias[0]] = self.air.val.pop(f)

        # fuel
        for f in list(self.fuel.val.keys()):
            alias = [x for x in self.air.val.keys() if x in [
                a.replace(' ', '') for a in CP.get_aliases(f)]]
            if len(alias) > 0:
                self.fuel.val[alias[0]] = self.fuel.val.pop(f)

        # list of all fluids of air and fuel
        fluids = list(self.air.val.keys()) + list(self.fuel.val.keys())

        # oxygen
        alias = [x for x in fluids if x in [
            a.replace(' ', '') for a in CP.get_aliases('O2')]]
        if len(alias) == 0:
            msg = 'Oxygen missing in input fluids.'
            logging.error(msg)
            raise TESPyComponentError(msg)
        else:
            self.o2 = alias[0]

        # carbondioxide
        self.co2 = [x for x in self.nw_fluids if x in [
            a.replace(' ', '') for a in CP.get_aliases('CO2')]]
        if len(self.co2) == 0:
            self.co2 = 'CO2'
        else:
            self.co2 = self.co2[0]

        # water
        self.h2o = [x for x in self.nw_fluids if x in [
            a.replace(' ', '') for a in CP.get_aliases('H2O')]]
        if len(self.h2o) == 0:
            self.h2o = 'H2O'
        else:
            self.h2o = self.h2o[0]

        # calculate lower heating value of specified fuel
        self.lhv = self.calc_lhv()
        msg = ('Combustion chamber fuel (' + self.fuel_alias.val +
               ') LHV is ' + str(self.lhv) + ' for component ' +
               self.label + '.')
        logging.debug(msg)
        # generate fluid properties for stoichiometric flue gas
        self.stoich_flue_gas(nw)

    def calc_lhv(self):
        r"""
        Calculate the lower heating value of the combustion chambers fuel.

        Returns
        -------
        val : float
            Lower heating value of the combustion chambers fuel.

            .. math::

                LHV = \sum_{fuels} \left(-\frac{\sum_i {\Delta H_f^0}_i -
                \sum_j {\Delta H_f^0}_j }
                {M_{fuel}} \cdot x_{fuel} \right)\\
                \forall i \in \text{reation products},\\
                \forall j \in \text{reation educts},\\
                \forall fuel \in \text{fuels},\\
                \Delta H_f^0: \text{molar formation enthalpy},\\
                x_{fuel}: \text{mass fraction of fuel in fuel mixture}
        """
        hf = {}
        hf['hydrogen'] = 0
        hf['methane'] = -74.85
        hf['ethane'] = -84.68
        hf['propane'] = -103.8
        hf['butane'] = -124.51
        hf['O2'] = 0
        hf['CO2'] = -393.5
        # water (gaseous)
        hf['H2O'] = -241.8

        lhv = 0

        for f, x in self.fuel.val.items():
            molar_masses[f] = CP.PropsSI('M', f)
            fl = set(list(hf.keys())).intersection(
                    set([a.replace(' ', '') for a in CP.get_aliases(f)]))
            if len(fl) == 0:
                continue

            if list(fl)[0] in self.fuels():
                structure = fluid_structure(f)

                n = {}
                for el in ['C', 'H', 'O']:
                    if el in structure:
                        n[el] = structure[el]
                    else:
                        n[el] = 0

                lhv += (-(n['H'] / 2 * hf['H2O'] + n['C'] * hf['CO2'] -
                          ((n['C'] + n['H'] / 4) * hf['O2'] +
                           hf[list(fl)[0]])) / molar_masses[f] * 1000) * x

        return lhv

    def stoich_flue_gas(self, nw):
        r"""
        Calculate the fluid composition of the stoichiometric flue gas.

        - uses one mole of fuel as reference quantity and :math:`\lambda=1`
          for stoichiometric flue gas calculation (no oxygen in flue gas)
        - calculate molar quantities of (reactive) fuel components to determine
          water and carbondioxide mass fraction in flue gas
        - calculate required molar quantity for oxygen and required fresh
          air mass
        - calculate residual mass fractions for non reactive components of
          fresh air in the flue gas
        - calculate flue gas fluid composition
        - generate custom fluid porperties

        Reactive components in fuel

        .. math::

            m_{fuel} = \frac{1}{M_{fuel}}\\
            m_{CO_2} = \sum_{i} \frac{x_{i} \cdot m_{fuel} \cdot num_{C,i}
            \cdot M_{CO_{2}}}{M_{i}}\\
            m_{H_{2}O} = \sum_{i} \frac{x_{i} \cdot m_{fuel} \cdot
            num_{H,i} \cdot M_{H_{2}O}}{2 \cdot M_{i}}\\
            \forall i \in \text{fuels in fuel vector},\\
            num = \text{number of atoms in molecule}

        Other components of fuel vector

        .. math::

            m_{fg,j} = x_{j} \cdot m_{fuel}\\
            \forall j \in \text{non fuels in fuel vecotr, e.g. } CO_2,\\
            m_{fg,j} = \text{mass of fluid component j in flue gas}

        Non-reactive components in air

        .. math::

            n_{O_2} = \left( \frac{m_{CO_2}}{M_{CO_2}} +
            \frac{m_{H_{2}O}}
            {0,5 \cdot M_{H_{2}O}} \right) \cdot \lambda,\\
            n_{O_2} = \text{mol of oxygen required}\\
            m_{air} = \frac{n_{O_2} \cdot M_{O_2}}{x_{O_{2}, air}},\\
            m_{air} = \text{required total air mass}\\
            m_{fg,j} = x_{j, air} \cdot m_{air}\\
            m_{fg, O_2} = 0,\\
            m_{fg,j} = \text{mass of fluid component j in flue gas}

        Flue gas composition

        .. math::

            x_{fg,j} = \frac{m_{fg, j}}{m_{air} + m_{fuel}}

        Parameters
        ----------
        nw : tespy.networks.network.Network
            TESPy network to generate stoichiometric flue gas for.
        """
        lamb = 1
        n_fuel = 1
        m_fuel = 1 / molar_mass_flow(self.fuel.val) * n_fuel
        m_fuel_fg = m_fuel
        m_co2 = 0
        m_h2o = 0
        molar_masses[self.h2o] = CP.PropsSI('M', self.h2o)
        molar_masses[self.co2] = CP.PropsSI('M', self.co2)
        molar_masses[self.o2] = CP.PropsSI('M', self.o2)

        self.fg = {}
        self.fg[self.co2] = 0
        self.fg[self.h2o] = 0

        for f, x in self.fuel.val.items():
            fl = set(list(self.fuels())).intersection(
                    set([a.replace(' ', '') for a in CP.get_aliases(f)]))

            if len(fl) == 0:
                if f in self.fg:
                    self.fg[f] += x * m_fuel
                else:
                    self.fg[f] = x * m_fuel
            else:
                n_fluid = x * m_fuel / molar_masses[f]
                m_fuel_fg -= n_fluid * molar_masses[f]
                structure = fluid_structure(f)
                n = {}
                for el in ['C', 'H', 'O']:
                    if el in structure:
                        n[el] = structure[el]
                    else:
                        n[el] = 0

                m_co2 += n_fluid * n['C'] * molar_masses[self.co2]
                m_h2o += n_fluid * n['H'] / 2 * molar_masses[self.h2o]

        self.fg[self.co2] += m_co2
        self.fg[self.h2o] += m_h2o

        n_o2 = (m_co2 / molar_masses[self.co2] +
                0.5 * m_h2o / molar_masses[self.h2o]) * lamb
        m_air = n_o2 * molar_masses[self.o2] / self.air.val[self.o2]

        self.air_min = m_air / m_fuel

        for f, x in self.air.val.items():
            if f != self.o2:
                if f in self.fg:
                    self.fg[f] += m_air * x
                else:
                    self.fg[f] = m_air * x

        m_fg = m_fuel + m_air

        for f in self.fg.keys():
            self.fg[f] /= m_fg

        if not self.path.is_set:
            self.path.val = None

        TESPyFluid(
            self.fuel_alias.val, self.fuel.val, [1000, nw.p_range_SI[1]],
            path=self.path.val)
        TESPyFluid(
            self.fuel_alias.val + '_fg', self.fg, [1000, nw.p_range_SI[1]],
            path=self.path.val)
        msg = (
            'Generated lookup table for ' + self.fuel_alias.val + ' and for '
            'stoichiometric flue gas at component ' + self.label + '.')
        logging.debug(msg)

        if self.air_alias.val not in ['Air', 'air']:
            TESPyFluid(
                self.air_alias.val, self.air.val, [1000, nw.p_range_SI[1]],
                path=self.path.val)
            msg = ('Generated lookup table for ' + self.air_alias.val +
                   ' at stoichiometric combustion chamber ' + self.label + '.')
        else:
            msg = ('Using CoolProp air at stoichiometric combustion chamber ' +
                   self.label + '.')
        logging.debug(msg)

    def stoichiometry(self, fluid):
        r"""
        Calculate the reaction balance for one fluid.

        - determine molar mass flows of fuel and oxygen
        - calculate excess fuel
        - calculate residual value of the fluids balance

        General equations

        .. math::

            res = \sum_i \left(x_{fluid,i} \cdot \dot{m}_{i}\right) -
            \sum_j \left(x_{fluid,j} \cdot \dot{m}_{j}\right)\\
            \forall i \in [1,2], \; \forall j \in [1]

            \dot{m}_{air,min} = \dot{m}_{fuel} \cdot air_{min}

            \lambda = \frac{\dot{m}_{air}}{\dot{m}_{air,min}}

        Equation for fuel

        .. math::

            0 = res - \left(\dot{m}_{f} - \dot{m}_{f,exc}\right)

            \dot{m}_{f,exc} = \begin{cases}
            0 & \lambda \geq 1\\
            \dot{m}_{f} - \frac{\dot{m}_{air}}
            {\lambda \cdot air_{min}} & \lambda < 1
            \end{cases}

        Equation for air

        .. math::

            0 = res - \begin{cases}
            -\dot{m}_{air,min} & \lambda \geq 1\\
            -\dot{m}_{air} & \lambda < 1
            \end{cases}

        Equation for stoichiometric flue gas

        .. math::

            0 = res + \dot{m}_{air,min} + \dot{m}_{f}

        Equation for all other fluids

        .. math::

            0 = res

        Parameters
        ----------
        fluid : str
            The fluid to calculate the reation balance for.

        Returns
        -------
        residual : float
            Residual value of equation.
        """
        air = self.air_alias.val
        fuel = self.fuel_alias.val
        flue_gas = self.fuel_alias.val + '_fg'

        ###################################################################
        # calculate fuel and air mass flow
        m_fuel = 0
        for i in self.inl:
            m_fuel += i.m.val_SI * i.fluid.val[fuel]

        m_air = 0
        for i in self.inl:
            m_air += i.m.val_SI * i.fluid.val[air]

        m_air_min = self.air_min * m_fuel

        ###################################################################
        # calculate lambda if not specified
        if not self.lamb.is_set:
            self.lamb.val = m_air / (self.air_min * m_fuel)

        ###################################################################
        # calculate excess fuel if lambda is smaller than 1
        m_fuel_exc = 0
        if self.lamb.val < 1:
            m_fuel_exc = m_fuel - m_air / (self.lamb.val * self.air_min)

        ###################################################################
        # equation for air
        if fluid == air:
            if self.lamb.val >= 1:
                dm = -m_air_min
            else:
                dm = -m_air

        ###################################################################
        # equation for fuel
        elif fluid == fuel:
            dm = -(m_fuel - m_fuel_exc)

        ###################################################################
        # equation for flue gas
        elif fluid == flue_gas:
            dm = m_air_min + m_fuel

        ###################################################################
        # equation for other components
        else:
            dm = 0

        res = dm
        for i in self.inl:
            res += i.fluid.val[fluid] * i.m.val_SI
        for o in self.outl:
            res -= o.fluid.val[fluid] * o.m.val_SI
        return res

    def stoichiometry_func_doc(self, label):
        r"""
        Calculate the reaction balance for one fluid.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        air = self.air_alias.val
        fuel = self.fuel_alias.val
        flue_gas = self.fuel_alias.val + '_fg'

        equations = ''
        for fluid in self.inl[0].fluid.val.keys():

            in1 = r'\dot{m}_\mathrm{in,1} \cdot x_\mathrm{fluid,in,1}'
            in2 = r'\dot{m}_\mathrm{in,2} \cdot x_\mathrm{fluid,in,2}'
            out = r'\dot{m}_\mathrm{out,1} \cdot x_\mathrm{fluid,out,1}'

            if fluid == air:
                latex = (
                    r'0=\Delta\dot{m}_\mathrm{' + fluid + r'} - '
                    r'\dot{m}_\mathrm{' + fluid + r',stoich}'
                )
            elif fluid == fuel:
                latex = (
                    r'0=\Delta \dot{m}_\mathrm{' + fluid + r'} - '
                    r'\dot{m}_\mathrm{' + fluid + r'}'
                )
            elif fluid == flue_gas:
                latex = (
                    r'0=\Delta \dot{m}_\mathrm{' + fluid.replace('_', ',') +
                    r'} + \dot{m}_\mathrm{' + fuel + r'} +'
                    r'\dot{m}_\mathrm{' + air + r',stoich}'
                )
            else:
                latex = r'0 = \Delta \dot{m}_\mathrm{' + fluid + '}'

            if fluid == next(iter(self.inl[0].fluid.val)):
                balance = (
                    r'\Delta \dot{m}_\mathrm{fluid} = ' + in1 +
                    '+' + in2 + '-' + out)
                m_fluid = r'\dot{m}_\mathrm{fluid} = ' + in1 + '+' + in2
                m_air_stoich = (
                    r'\dot{m}_\mathrm{' + air + ',stoich}='
                    r'\dot{m}_\mathrm{' + fuel + r'} \cdot ' +
                    str(round(self.air_min, 4)))
                latex_general_eq = (
                    r'\begin{split}' + '\n'
                    r'&' + balance + r'\\' + '\n'
                    r'&' + m_fluid + r'\\' + '\n'
                    r'&' + m_air_stoich + r'\\' + '\n'
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
        The temperature for the reference state is set to 100 °C, as the
        custom fluid properties are inacurate at the dew-point of water in
        the flue gas!

        - Reference temperature: 373.15 K.
        - Reference pressure: 1 bar.
        """
        T_ref = 373.15
        p_ref = 1e5

        res = 0
        for i in self.inl:
            res += i.m.val_SI * (
                i.h.val_SI - h_mix_pT([0, p_ref, 0, i.fluid.val], T_ref))
        for o in self.outl:
            res -= o.m.val_SI * (
                o.h.val_SI - h_mix_pT([0, p_ref, 0, o.fluid.val], T_ref))

        return res + self.calc_ti()

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
            r'& T_\mathrm{ref}=\unit[373.15]{K}'
            r'\;p_\mathrm{ref}=\unit[10^5]{Pa}\\'
            '\n' + r'\end{split}'
        )
        return generate_latex_eq(self, latex, label)

    def lambda_func(self):
        r"""
        Calculate the residual for specified lambda.

        Returns
        -------
        residual : float
            Residual value of function.

            .. math::

                0 = \lambda - \frac{\dot{m}_{air}}{\dot{m}_{air,min}}
        """
        air = self.air_alias.val
        fuel = self.fuel_alias.val

        m_air = 0
        m_fuel = 0

        for i in self.inl:
            m_air += (i.m.val_SI * i.fluid.val[air])
            m_fuel += (i.m.val_SI * i.fluid.val[fuel])

        return self.lamb.val - m_air / (m_fuel * self.air_min)

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
            r'0 = \lambda - '
            r'\frac{\dot{m}_\mathrm{air}}{\dot{m}_\mathrm{air,min}}'
        )
        return generate_latex_eq(self, latex, label)

    def ti_func(self):
        r"""
        Calculate the residual for specified thermal input.

        Returns
        -------
        residual : float
            Residual value of function.

            .. math::

                0 = ti - \dot{m}_f \cdot LHV
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
        latex = (
            r'\begin{split}' + '\n'
            r'0 = & ti - LHV_\mathrm{fuel} \cdot \left[\sum_i \left('
            r'\dot{m}_{\mathrm{in,}i}\cdot x_{\mathrm{fuel,in,}i}\right)-'
            r' \dot{m}_\mathrm{out,1}\cdot '
            r'x_{\mathrm{fuel,out,1}} \right]\\' + '\n'
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
        fuel = self.fuel_alias.val
        pos = 3 + self.nw_fluids.index(fuel)

        for i in range(2):
            self.jacobian[k, i, 0] = -self.inl[i].fluid.val[fuel]
            self.jacobian[k, i, pos] = -self.inl[i].m.val_SI
        self.jacobian[k, 2, 0] = self.outl[0].fluid.val[fuel]
        self.jacobian[k, 2, pos] = self.outl[0].m.val_SI
        self.jacobian[k] *= self.lhv
        k += 1

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
        m = 0
        for i in self.inl:
            m += i.m.val_SI * i.fluid.val[self.fuel_alias.val]

        for o in self.outl:
            m -= o.m.val_SI * o.fluid.val[self.fuel_alias.val]

        return m * self.lhv

    def initialise_fluids(self):
        """Calculate reaction balance for generic starting values at outlet."""
        air = self.air_alias.val
        flue_gas = self.fuel_alias.val + '_fg'

        for c in self.outl:
            if not c.fluid.val_set[air]:
                c.fluid.val[air] = 0.8
            if not c.fluid.val_set[flue_gas]:
                c.fluid.val[flue_gas] = 0.2
            c.target.propagate_fluid_to_target(c, c.target)

    def convergence_check(self):
        r"""
        Perform a convergence check.

        Note
        ----
        Manipulate enthalpies/pressure at inlet and outlet if not specified by
        user to match physically feasible constraints, keep fluid composition
        within feasible range and then propagates it towards the outlet.
        """
        air = self.air_alias.val
        flue_gas = self.fuel_alias.val + '_fg'
        fuel = self.fuel_alias.val

        for c in self.outl:
            if not c.fluid.val_set[air]:
                if c.fluid.val[air] > 0.95:
                    c.fluid.val[air] = 0.95
                if c.fluid.val[air] < 0.5:
                    c.fluid.val[air] = 0.5

            if not c.fluid.val_set[flue_gas]:
                if c.fluid.val[flue_gas] > 0.5:
                    c.fluid.val[flue_gas] = 0.5
                if c.fluid.val[flue_gas] < 0.05:
                    c.fluid.val[flue_gas] = 0.05

            if not c.fluid.val_set[fuel]:
                if c.fluid.val[fuel] > 0:
                    c.fluid.val[fuel] = 0

            c.target.propagate_fluid_to_target(c, c.target)

        for i in self.inl:
            if i.m.val_SI < 0 and not i.m.val_set:
                i.m.val_SI = 0.01

        for c in self.outl:
            if c.m.val_SI < 0 and not c.m.val_set:
                c.m.val_SI = 10
            c.target.propagate_fluid_to_target(c, c.target)

        if self.lamb.val < 1 and not self.lamb.is_set:
            self.lamb.val = 2

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        air = self.air_alias.val
        fuel = self.fuel_alias.val

        m_fuel = 0
        for i in self.inl:
            m_fuel += i.m.val_SI * i.fluid.val[fuel]
        m_air = 0
        for i in self.inl:
            m_air += i.m.val_SI * i.fluid.val[air]
        self.lamb.val = (m_air / m_fuel) / self.air_min

        ti = 0
        for i in self.inl:
            ti += i.m.val_SI * i.fluid.val[fuel] * self.lhv
        self.ti.val = ti

    def entropy_balance(self):
        r"""Entropy balance calculation method."""
        return
