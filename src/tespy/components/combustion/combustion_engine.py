# -*- coding: utf-8

"""Module of class CombustionEngine.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/combustion/combustion_engine.py

SPDX-License-Identifier: MIT
"""

import logging

import numpy as np

from tespy.components.combustion.combustion_chamber import CombustionChamber
from tespy.components.component import Component
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import DataContainerSimple as dc_simple
from tespy.tools.document_models import generate_latex_eq
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.fluid_properties import s_mix_ph
from tespy.tools.fluid_properties import s_mix_pT
from tespy.tools.global_vars import molar_masses


class CombustionEngine(CombustionChamber):
    r"""
    An internal combustion engine supplies power and heat cogeneration.

    The combustion engine produces power and heat in cogeneration from fuel
    combustion. The combustion properties are identical to the combustion
    chamber. Thermal input and power output, heat output and heat losses are
    linked with an individual characteristic line for each property.

    **Mandatory Equations**

    - :py:meth:`tespy.components.combustion.combustion_engine.CombustionEngine.fluid_func`
      (for cooling water)
    - :py:meth:`tespy.components.combustion.combustion_engine.CombustionEngine.mass_flow_func`
    - :py:meth:`tespy.components.combustion.combustion_chamber.CombustionChamber.combustion_pressure_func`
    - :py:meth:`tespy.components.combustion.combustion_chamber.CombustionChamber.stoichiometry`
    - :py:meth:`tespy.components.combustion.combustion_engine.CombustionEngine.energy_balance_func`
    - :py:meth:`tespy.components.combustion.combustion_engine.CombustionEngine.tiP_char_func`
    - :py:meth:`tespy.components.combustion.combustion_engine.CombustionEngine.Q1_char_func`
    - :py:meth:`tespy.components.combustion.combustion_engine.CombustionEngine.Q2_char_func`
    - :py:meth:`tespy.components.combustion.combustion_engine.CombustionEngine.Qloss_char_func`

    **Optional Equations**

    - :py:meth:`tespy.components.combustion.combustion_chamber.CombustionChamber.lambda_func`
    - :py:meth:`tespy.components.combustion.combustion_chamber.CombustionChamber.ti_func`
    - :py:meth:`tespy.components.combustion.combustion_engine.CombustionEngine.Q1_func`
    - :py:meth:`tespy.components.combustion.combustion_engine.CombustionEngine.Q2_func`

    - cooling loops:

      - 1 :py:meth:`tespy.components.component.Component.pr_func`
      - 2 :py:meth:`tespy.components.component.Component.pr_func`
      - 1 :py:meth:`tespy.components.component.Component.zeta_func`
      - 2 :py:meth:`tespy.components.component.Component.zeta_func`

    Available fuels

    - methane, ethane, propane, butane, hydrogen

    Inlets/Outlets

    - in1, in2 (cooling water), in3, in4 (air and fuel)
    - out1, out2 (cooling water), out3 (flue gas)

    Image

    .. image:: _images/CombustionEngine.svg
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

    lamb : float, dict
        Air to stoichiometric air ratio, :math:`\lambda/1`.

    ti : float, dict
        Thermal input, (:math:`{LHV \cdot \dot{m}_f}`), :math:`ti/\text{W}`.

    P : float, dict, :code:`"var"`
        Power output, :math:`P/\text{W}`.

    Q1 : float, dict
        Heat output 1, :math:`\dot Q/\text{W}`.

    Q2 : float, dict
        Heat output 2, :math:`\dot Q/\text{W}`.

    Qloss : float, dict, :code:`"var"`
        Heat loss, :math:`\dot Q_{loss}/\text{W}`.

    pr1 : float, dict, :code:`"var"`
        Pressure ratio heat outlet 1, :math:`pr/1`.

    pr2 : float, dict, :code:`"var"`
        Pressure ratio heat outlet 2, :math:`pr/1`.

    zeta1 : float, dict, :code:`"var"`
        Geometry independent friction coefficient heating loop 1,
        :math:`\zeta_1/\frac{1}{\text{m}^4}`.

    zeta2 : float, dict, :code:`"var"`
        Geometry independent friction coefficient heating loop 2,
        :math:`\zeta_2/\frac{1}{\text{m}^4}`.

    tiP_char : tespy.tools.characteristics.CharLine, dict
        Characteristic line linking fuel input to power output.

    Q1_char : tespy.tools.characteristics.CharLine, dict
        Characteristic line linking heat output 1 to power output.

    Q2_char : tespy.tools.characteristics.CharLine, dict
        Characteristic line linking heat output 2 to power output.

    Qloss_char : tespy.tools.characteristics.CharLine, dict
        Characteristic line linking heat loss to power output.

    eta_mech : float
        Value of internal efficiency of the combustion engine. This value is
        required to determine the (virtual) thermodynamic temperature of heat
        inside the combustion engine for the entropy balance calculation.
        Default value is 0.85.

    Note
    ----
    Parameters available through entropy balance are listed in the respective
    method:

    - :py:meth:`tespy.components.combustion.combustion_engine.CombustionEngine.entropy_balance`

    Example
    -------
    The combustion chamber calculates energy input due to combustion as well as
    the flue gas composition based on the type of fuel and the amount of
    oxygen supplied. In this example a mixture of methane, hydrogen and
    carbondioxide is used as fuel. There are two cooling ports, the cooling
    water will flow through them in parallel.

    >>> from tespy.components import (Sink, Source, CombustionEngine, Merge,
    ... Splitter)
    >>> from tespy.connections import Connection, Ref
    >>> from tespy.networks import Network
    >>> import numpy as np
    >>> import shutil
    >>> fluid_list = ['Ar', 'N2', 'O2', 'CO2', 'CH4', 'H2O']
    >>> nw = Network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ... iterinfo=False)
    >>> amb = Source('ambient')
    >>> sf = Source('fuel')
    >>> fg = Sink('flue gas outlet')
    >>> cw_in = Source('cooling water inlet')
    >>> sp = Splitter('cooling water splitter', num_out=2)
    >>> me = Merge('cooling water merge', num_in=2)
    >>> cw_out = Sink('cooling water outlet')
    >>> chp = CombustionEngine(label='internal combustion engine')
    >>> chp.component()
    'combustion engine'
    >>> amb_comb = Connection(amb, 'out1', chp, 'in3')
    >>> sf_comb = Connection(sf, 'out1', chp, 'in4')
    >>> comb_fg = Connection(chp, 'out3', fg, 'in1')
    >>> nw.add_conns(sf_comb, amb_comb, comb_fg)
    >>> cw_sp = Connection(cw_in, 'out1', sp, 'in1')
    >>> sp_chp1 = Connection(sp, 'out1', chp, 'in1')
    >>> sp_chp2 = Connection(sp, 'out2', chp, 'in2')
    >>> chp1_me = Connection(chp, 'out1', me, 'in1')
    >>> chp2_me = Connection(chp, 'out2', me, 'in2')
    >>> me_cw = Connection(me, 'out1', cw_out, 'in1')
    >>> nw.add_conns(cw_sp, sp_chp1, sp_chp2, chp1_me, chp2_me, me_cw)

    The combustion engine produces a power output of 10 MW the oxygen to
    stoichiometric oxygen ratio is set to 1. Only pressure ratio 1 is set as
    we reconnect both cooling water streams. At the merge all pressure values
    will be identical automatically. Reference the mass flow at the splitter
    to be split in half.

    >>> chp.set_attr(pr1=0.99, P=-10e6, lamb=1.0,
    ... design=['pr1'], offdesign=['zeta1'])
    >>> amb_comb.set_attr(p=5, T=30, fluid={'Ar': 0.0129, 'N2': 0.7553,
    ... 'H2O': 0, 'CH4': 0, 'CO2': 0.0004, 'O2': 0.2314})
    >>> sf_comb.set_attr(m0=0.1, T=30, fluid={'CO2': 0, 'Ar': 0, 'N2': 0,
    ... 'O2': 0, 'H2O': 0, 'CH4': 1})
    >>> cw_sp.set_attr(p=3, T=60, m=50, fluid={'CO2': 0, 'Ar': 0, 'N2': 0,
    ... 'O2': 0, 'H2O': 1, 'CH4': 0})
    >>> sp_chp2.set_attr(m=Ref(sp_chp1, 1, 0))
    >>> mode = 'design'
    >>> nw.solve(mode=mode)
    >>> nw.save('tmp')
    >>> round(chp.ti.val, 0)
    25300000.0
    >>> round(chp.Q1.val, 0)
    -4980000.0
    >>> chp.set_attr(Q1=-4e6, P=np.nan)
    >>> mode = 'offdesign'
    >>> nw.solve(mode=mode, init_path='tmp', design_path='tmp')
    >>> round(chp.ti.val, 0)
    17794554.0
    >>> round(chp.P.val / chp.P.design, 3)
    0.617
    >>> chp.set_attr(P=chp.P.design * 0.75, Q1=np.nan)
    >>> mode = 'offdesign'
    >>> nw.solve(mode=mode, init_path='tmp', design_path='tmp')
    >>> round(chp.ti.val, 0)
    20550000.0
    >>> round(chp.P.val / chp.P.design, 3)
    0.75
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    @staticmethod
    def component():
        return 'combustion engine'

    def get_variables(self):
        return {
            'lamb': dc_cp(
                min_val=1, deriv=self.lambda_deriv, func=self.lambda_func,
                latex=self.lambda_func_doc, num_eq=1),
            'ti': dc_cp(
                min_val=0, deriv=self.ti_deriv, func=self.ti_func,
                latex=self.ti_func_doc, num_eq=1),
            'P': dc_cp(val=-1e6, d=1, max_val=-1),
            'Q1': dc_cp(
                max_val=-1, deriv=self.Q1_deriv, func=self.Q1_func,
                num_eq=1, latex=self.Q1_func_doc),
            'Q2': dc_cp(
                max_val=-1, deriv=self.Q2_deriv, func=self.Q2_func,
                num_eq=1, latex=self.Q2_func_doc),
            'Qloss': dc_cp(val=-1e5, d=1, max_val=-1),
            'pr1': dc_cp(
                min_val=1e-4, max_val=1, num_eq=1, deriv=self.pr_deriv,
                latex=self.pr_func_doc,
                func=self.pr_func, func_params={'pr': 'pr1'}),
            'pr2': dc_cp(
                min_val=1e-4, max_val=1, num_eq=1, latex=self.pr_func_doc,
                deriv=self.pr_deriv, func=self.pr_func,
                func_params={'pr': 'pr2', 'inconn': 1, 'outconn': 1}),
            'zeta1': dc_cp(
                min_val=0, max_val=1e15, num_eq=1, latex=self.zeta_func_doc,
                deriv=self.zeta_deriv, func=self.zeta_func,
                func_params={'zeta': 'zeta1'}),
            'zeta2': dc_cp(
                min_val=0, max_val=1e15, num_eq=1, latex=self.zeta_func_doc,
                deriv=self.zeta_deriv, func=self.zeta_func,
                func_params={'zeta': 'zeta2', 'inconn': 1, 'outconn': 1}),
            'tiP_char': dc_cc(), 'Q1_char': dc_cc(), 'Q2_char': dc_cc(),
            'Qloss_char': dc_cc(),
            'eta_mech': dc_simple(val=0.85), 'T_v_inner': dc_simple()}

    def get_mandatory_constraints(self):
        return {
            'mass_flow_constraints': {
                'func': self.mass_flow_func, 'deriv': self.mass_flow_deriv,
                'constant_deriv': True, 'latex': self.mass_flow_func_doc,
                'num_eq': 3},
            'cooling_loop_fuid_constraints': {
                'func': self.fluid_func, 'deriv': self.fluid_deriv,
                'constant_deriv': True, 'latex': self.fluid_func_doc,
                'num_eq': 2 * self.num_nw_fluids},
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
                'num_eq': self.num_nw_fluids},
            'energy_balance_constraints': {
                'func': self.energy_balance_func,
                'deriv': self.energy_balance_deriv,
                'constant_deriv': False, 'latex': self.energy_balance_func_doc,
                'num_eq': 1},
            'power_constraints': {
                'func': self.tiP_char_func,
                'deriv': self.tiP_char_deriv,
                'constant_deriv': False, 'latex': self.tiP_char_func_doc,
                'num_eq': 1, 'char': 'tiP_char'},
            'heat1_constraints': {
                'func': self.Q1_char_func,
                'deriv': self.Q1_char_deriv,
                'constant_deriv': False, 'latex': self.Q1_char_func_doc,
                'num_eq': 1, 'char': 'Q1_char'},
            'heat2_constraints': {
                'func': self.Q2_char_func,
                'deriv': self.Q2_char_deriv,
                'constant_deriv': False, 'latex': self.Q2_char_func_doc,
                'num_eq': 1, 'char': 'Q2_char'},
            'heatloss_constraints': {
                'func': self.Qloss_char_func,
                'deriv': self.Qloss_char_deriv,
                'constant_deriv': False, 'latex': self.Qloss_char_func_doc,
                'num_eq': 1, 'char': 'Qloss_char'},
        }

    @staticmethod
    def inlets():
        return ['in1', 'in2', 'in3', 'in4']

    @staticmethod
    def outlets():
        return ['out1', 'out2', 'out3']

    def comp_init(self, nw):

        if not self.P.is_set:
            self.set_attr(P='var')
            msg = ('The power output of combustion engines must be set! '
                   'We are adding the power output of component ' +
                   self.label + ' as custom variable of the system.')
            logging.info(msg)

        if not self.Qloss.is_set:
            self.set_attr(Qloss='var')
            msg = ('The heat loss of combustion engines must be set! '
                   'We are adding the heat loss of component ' +
                   self.label + ' as custom variable of the system.')
            logging.info(msg)

        Component.comp_init(self, nw)

        self.setup_reaction_parameters()

    def mass_flow_func(self):
        r"""
        Calculate the residual value for component's mass flow balance.

        Returns
        -------
        residual : list
            Vector with residual value for component's mass flow balance.

            .. math::

                0 = \dot{m}_{in,i} - \dot{m}_{out,i}\\
                \forall i \in [1, 2]\\
                0 = \dot{m}_{in,3} + \dot{m}_{in,4} - \dot{m}_{out,3}
        """
        residual = []
        for i in range(2):
            residual += [self.inl[i].m.val_SI - self.outl[i].m.val_SI]
        residual += [self.inl[2].m.val_SI + self.inl[3].m.val_SI -
                     self.outl[2].m.val_SI]
        return residual

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
            r'\begin{split}' + '\n'
            r'0=&\dot{m}_\mathrm{in,1} - \dot{m}_\mathrm{out,1}\\' + '\n'
            r'0=&\dot{m}_\mathrm{in,2} - \dot{m}_\mathrm{out,2}\\' + '\n'
            r'0=&\dot{m}_\mathrm{in,3} + \dot{m}_\mathrm{in,3} - '
            r'\dot{m}_\mathrm{out,3}\\' + '\n'
            r'\end{split}'
        )
        return generate_latex_eq(self, latex, label)

    def mass_flow_deriv(self):
        r"""
        Calculate the partial derivatives for all mass flow balance equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((3, 7 + self.num_vars, self.num_nw_vars))
        for i in range(2):
            deriv[i, i, 0] = 1
        for j in range(2):
            deriv[j, self.num_i + j, 0] = -1
        deriv[2, 2, 0] = 1
        deriv[2, 3, 0] = 1
        deriv[2, 6, 0] = -1
        return deriv

    def fluid_func(self):
        r"""
        Calculate the vector of residual values for cooling loop fluid balance.

        Returns
        -------
        residual : list
            Vector of residual values for component's fluid balance.

            .. math::

                0 = fluid_{i,in_{j}} - fluid_{i,out_{j}}\\
                \forall i \in \mathrm{fluid}, \; \forall j \in [1, 2]
        """
        residual = []
        for i in range(2):
            for fluid, x in self.inl[i].fluid.val.items():
                residual += [x - self.outl[i].fluid.val[fluid]]
        return residual

    def fluid_func_doc(self, label):
        r"""
        Calculate the vector of residual values for cooling loop fluid balance.

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
            r'0=x_{i\mathrm{,in,}j}-x_{i\mathrm{,out,}j}\;'
            r'\forall i \in\text{network fluids,}'
            r'\; \forall j \in [1,2]')
        return generate_latex_eq(self, latex, label)

    def fluid_deriv(self):
        r"""
        Calculate the partial derivatives for cooling loop fluid balance.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros(
            (self.num_nw_fluids * 2, 7 + self.num_vars, self.num_nw_vars))
        for i in range(self.num_nw_fluids):
            deriv[i, 0, i + 3] = 1
            deriv[i, 4, i + 3] = -1
        for j in range(self.num_nw_fluids):
            deriv[i + 1 + j, 1, j + 3] = 1
            deriv[i + 1 + j, 5, j + 3] = -1
        return deriv

    def energy_balance_func(self):
        r"""
        Calculate the energy balance of the combustion engine.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                \begin{split}
                0 = & \sum_i \dot{m}_{in,i} \cdot
                \left( h_{in,i} - h_{in,i,ref} \right)\\
                & - \sum_j \dot{m}_{out,3} \cdot
                \left( h_{out,3} - h_{out,3,ref} \right)\\
                & + LHV_{fuel} \cdot
                \left(\sum_i \left(\dot{m}_{in,i} \cdot x_{fuel,i} \right)-
                \dot{m}_{out,3} \cdot x_{fuel,3} \right)\\
                & - \dot{Q}_1 - \dot{Q}_2 - P - \dot{Q}_{loss}\\
                \end{split}\\
                \forall i \in [3,4]

        Note
        ----
        The temperature for the reference state is set to 25 °C, thus
        the water may be liquid. In order to make sure, the state is
        referring to the lower heating value, the necessary enthalpy
        difference for evaporation is added.

        - Reference temperature: 298.15 K.
        - Reference pressure: 1 bar.
        """
        T_ref = 298.15
        p_ref = 1e5

        res = 0
        for i in self.inl[2:]:
            res += i.m.val_SI * (i.h.val_SI - h_mix_pT(
                [0, p_ref, 0, i.fluid.val], T_ref, force_gas=True))

        for o in self.outl[2:]:
            res -= o.m.val_SI * (o.h.val_SI - h_mix_pT(
                [0, p_ref, 0, o.fluid.val], T_ref, force_gas=True))

        res += self.calc_ti()

        # cooling
        for i in range(2):
            res -= self.inl[i].m.val_SI * (
                self.outl[i].h.val_SI - self.inl[i].h.val_SI)

        # power output and heat loss
        res += self.P.val + self.Qloss.val

        return res

    def energy_balance_func_doc(self, label):
        """
        Calculate the energy balance of the combustion engine.

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
            r'\dot{m}_\mathrm{out,3}\cdot\left( h_\mathrm{out,3}'
            r' - h_\mathrm{out,3,ref}\right)\\' + '\n'
            r'& + LHV_{fuel} \cdot \left(\sum_i \dot{m}_{\mathrm{in,}i} '
            r'\cdot x_{fuel\mathrm{,in,}i} - \dot{m}_\mathrm{out,3} '
            r'\cdot x_{fuel\mathrm{,out,3}} \right)\\' + '\n'
            r'& + \dot{Q}_1 + \dot{Q}_2+P + \dot{Q}_\mathrm{loss}\\' + '\n'
            r'& \forall i \in [3,4]\\'
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
        # mass flow cooling water
        for i in [0, 1]:
            self.jacobian[k, i, 0] = -(
                self.outl[i].h.val_SI - self.inl[i].h.val_SI)

        # mass flow and pressure for combustion reaction
        for i in [2, 3, 6]:
            if not increment_filter[i, 0]:
                self.jacobian[k, i, 0] = self.numeric_deriv(f, 'm', i)
            if not increment_filter[i, 1]:
                self.jacobian[k, i, 1] = self.numeric_deriv(f, 'p', i)

        # enthalpy
        for i in range(4):
            self.jacobian[k, i, 2] = self.inl[i].m.val_SI
        for i in range(3):
            self.jacobian[k, i + 4, 2] = -self.outl[i].m.val_SI

        # fluid composition
        for fl in self.fuel_list:
            pos = 3 + self.nw_fluids.index(fl)
            lhv = self.fuels[fl]['LHV']
            self.jacobian[k, 2, pos] = self.inl[2].m.val_SI * lhv
            self.jacobian[k, 3, pos] = self.inl[3].m.val_SI * lhv
            self.jacobian[k, 6, pos] = -self.outl[2].m.val_SI * lhv

        # power and heat loss
        if self.P.is_var:
            self.jacobian[k, 7 + self.P.var_pos, 0] = (
                self.numeric_deriv(f, 'P', 7))
        if self.Qloss.is_var:
            self.jacobian[k, 7 + self.Qloss.var_pos, 0] = (
                self.numeric_deriv(f, 'Qloss', 7))

    def Q1_func(self):
        r"""
        Calculate residual value of primary heat loop function.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_1 \cdot \left(h_{out,1} +
                h_{in,1} \right) + \dot{Q}_1
        """
        i = self.inl[0]
        o = self.outl[0]
        return i.m.val_SI * (o.h.val_SI - i.h.val_SI) + self.Q1.val

    def Q1_func_doc(self, label):
        r"""
        Calculate residual value of primary heat loop function.

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
            r'0 = \dot{m}_\mathrm{in,1} \cdot \left(h_\mathrm{out,1} +'
            r'h_\mathrm{in,1} \right) + \dot{Q}_1')
        return generate_latex_eq(self, latex, label)

    def Q1_deriv(self, increment_filter, k):
        """
        Calculate partial derivatives of primary heat loop function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of equation in Jacobian matrix.
        """
        self.jacobian[k, 0, 0] = self.outl[0].h.val_SI - self.inl[0].h.val_SI
        self.jacobian[k, 0, 2] = -self.inl[0].m.val_SI
        self.jacobian[k, 4, 2] = self.inl[0].m.val_SI

    def Q2_func(self):
        r"""
        Calculate residual value of secondary heat loop function.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_2 \cdot \left(h_{out,2} - h_{in,2} \right) +
                \dot{Q}_2
        """
        i = self.inl[1]
        o = self.outl[1]
        return i.m.val_SI * (o.h.val_SI - i.h.val_SI) + self.Q2.val

    def Q2_func_doc(self, label):
        r"""
        Calculate residual value of secondary heat loop function.

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
            r'0 = \dot{m}_\mathrm{in,2} \cdot \left(h_\mathrm{out,2} +'
            r'h_\mathrm{in,2} \right) + \dot{Q}_2')
        return generate_latex_eq(self, latex, label)

    def Q2_deriv(self, increment_filter, k):
        """
        Calculate partial derivatives of secondary heat loop function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of equation in Jacobian matrix.
        """
        self.jacobian[k, 1, 0] = self.outl[1].h.val_SI - self.inl[1].h.val_SI
        self.jacobian[k, 1, 2] = -self.inl[1].m.val_SI
        self.jacobian[k, 5, 2] = self.inl[1].m.val_SI

    def tiP_char_func(self):
        r"""
        Calculate the relation of output power and thermal input.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = P \cdot f_{TI}\left(\frac{P}{P_{design}}\right)+ LHV \cdot
                \left[\sum_i \left(\dot{m}_{in,i} \cdot
                x_{f,i}\right) - \dot{m}_{out,3} \cdot x_{f,3} \right]
                \; \forall i \in [3,4]
        """
        if np.isnan(self.P.design):
            expr = 1
        else:
            expr = self.P.val / self.P.design

        return (
            self.calc_ti() + self.tiP_char.char_func.evaluate(expr) *
            self.P.val)

    def tiP_char_func_doc(self, label):
        r"""
        Calculate the relation of output power and thermal input.

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
            r'0=&P\cdot f_\mathrm{TI}\left(\frac{P}{P_\mathrm{design}}'
            r'\right)\\ ' + '\n'
            r'&+ LHV_{fuel} \cdot \left[\sum_i \left('
            r'\dot{m}_{\mathrm{in,}i} \cdot x_{fuel\mathrm{,in,}i}\right)'
            r'-\dot{m}_\mathrm{out,3}\cdot x_{fuel\mathrm{,out,}3}'
            r'\right]\\' + '\n'
            r'&\forall i \in [3,4]\\ ' + '\n'
            r'\end{split}'
        )
        return generate_latex_eq(self, latex, label)

    def tiP_char_deriv(self, increment_filter, k):
        """
        Calculate partial derivatives of power to thermal input characteristic.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of equation in Jacobian matrix.
        """
        f = self.tiP_char_func
        for i in [2, 3, 6]:
            if not increment_filter[i, 0]:
                self.jacobian[k, i, 0] = self.numeric_deriv(f, 'm', i)
            if not all(increment_filter[i, 3:]):
                self.jacobian[k, i, 3:] = self.numeric_deriv(f, 'fluid', i)

        if self.P.is_var:
            self.jacobian[k, 7 + self.P.var_pos, 0] = (
                self.numeric_deriv(f, 'P', 7))

    def Q1_char_func(self):
        r"""
        Calculate the relation of heat output 1 and thermal input.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                \begin{split}
                0 = & \dot{m}_1 \cdot \left(h_{out,1} - h_{in,1} \right) \cdot
                f_{TI}\left(\frac{P}{P_{design}}\right) \\
                & - LHV \cdot \left[\sum_i
                \left(\dot{m}_{in,i} \cdot x_{f,i}\right) -
                \dot{m}_{out,3} \cdot x_{f,3} \right] \cdot
                f_{Q1}\left(\frac{P}{P_{ref}}\right)\\
                \end{split}\\
                \forall i \in [3,4]
        """
        i = self.inl[0]
        o = self.outl[0]

        if np.isnan(self.P.design):
            expr = 1
        else:
            expr = self.P.val / self.P.design

        return (self.calc_ti() * self.Q1_char.char_func.evaluate(expr) -
                self.tiP_char.char_func.evaluate(expr) * i.m.val_SI *
                (o.h.val_SI - i.h.val_SI))

    def Q1_char_func_doc(self, label):
        r"""
        Calculate the relation of heat output 1 and thermal input.

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
            r'0=&LHV_{fuel} \cdot \left[\sum_i \left('
            r'\dot{m}_{\mathrm{in,}i} \cdot x_{fuel\mathrm{,in,}i}\right)'
            r'-\dot{m}_\mathrm{out,3}\cdot x_{fuel\mathrm{,out,}3}'
            r'\right] \cdot f_\mathrm{Q1}\left(\frac{P}{P_\mathrm{design}}'
            r'\right)\\' + '\n'
            r'&-\dot{m}_\mathrm{in,1} \cdot \left( h_\mathrm{out,1} - '
            r'h_\mathrm{in,1}\right) \cdot f_\mathrm{TI}'
            r'\left(\frac{P}{P_\mathrm{design}}'
            r'\right)\\ ' + '\n'
            r'&\forall i \in [3,4]\\ ' + '\n'
            r'\end{split}'
        )
        return generate_latex_eq(self, latex, label)

    def Q1_char_deriv(self, increment_filter, k):
        """
        Calculate partial derivatives of primary heat to thermal input char.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of equation in Jacobian matrix.
        """
        f = self.Q1_char_func
        if not increment_filter[0, 0]:
            self.jacobian[k, 0, 0] = self.numeric_deriv(f, 'm', 0)
        if not increment_filter[0, 2]:
            self.jacobian[k, 0, 2] = self.numeric_deriv(f, 'h', 0)
        if not increment_filter[4, 2]:
            self.jacobian[k, 4, 2] = self.numeric_deriv(f, 'h', 4)
        for i in [2, 3, 6]:
            if not increment_filter[i, 0]:
                self.jacobian[k, i, 0] = self.numeric_deriv(f, 'm', i)
            if not all(increment_filter[i, 3:]):
                self.jacobian[k, i, 3:] = self.numeric_deriv(f, 'fluid', i)

        if self.P.is_var:
            self.jacobian[k, 7 + self.P.var_pos, 0] = (
                self.numeric_deriv(f, 'P', 7))

    def Q2_char_func(self):
        r"""
        Calculate the relation of heat output 2 and thermal input.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                \begin{split}
                0 = & \dot{m}_2 \cdot \left(h_{out,2} - h_{in,2} \right) \cdot
                f_{TI}\left(\frac{P}{P_{design}}\right) \\
                & - LHV \cdot \left[\sum_i
                \left(\dot{m}_{in,i} \cdot x_{f,i}\right) -
                \dot{m}_{out,3} \cdot x_{f,3} \right] \cdot
                f_{Q2}\left(\frac{P}{P_{ref}}\right)\\
                \end{split}\\
                \forall i \in [3,4]
        """
        i = self.inl[1]
        o = self.outl[1]

        if np.isnan(self.P.design):
            expr = 1
        else:
            expr = self.P.val / self.P.design

        return (self.calc_ti() * self.Q2_char.char_func.evaluate(expr) -
                self.tiP_char.char_func.evaluate(expr) * i.m.val_SI *
                (o.h.val_SI - i.h.val_SI))

    def Q2_char_func_doc(self, label):
        r"""
        Calculate the relation of heat output 2 and thermal input.

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
            r'0=&LHV_{fuel} \cdot \left[\sum_i \left('
            r'\dot{m}_{\mathrm{in,}i} \cdot x_{fuel\mathrm{,in,}i}\right)'
            r'-\dot{m}_\mathrm{out,3}\cdot x_{fuel\mathrm{,out,}3}'
            r'\right] \cdot f_\mathrm{Q2}\left(\frac{P}{P_\mathrm{design}}'
            r'\right)\\' + '\n'
            r'&-\dot{m}_\mathrm{in,2} \cdot \left( h_\mathrm{out,2} - '
            r'h_\mathrm{in,2}\right) \cdot f_\mathrm{TI}'
            r'\left(\frac{P}{P_\mathrm{design}}'
            r'\right)\\ ' + '\n'
            r'&\forall i \in [3,4]\\ ' + '\n'
            r'\end{split}'
        )
        return generate_latex_eq(self, latex, label)

    def Q2_char_deriv(self, increment_filter, k):
        """
        Calculate partial derivatives of secondary heat to thermal input char.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of equation in Jacobian matrix.
        """
        f = self.Q2_char_func
        if not increment_filter[1, 0]:
            self.jacobian[k, 1, 0] = self.numeric_deriv(f, 'm', 1)
        if not increment_filter[1, 2]:
            self.jacobian[k, 1, 2] = self.numeric_deriv(f, 'h', 1)
        if not increment_filter[5, 2]:
            self.jacobian[k, 5, 2] = self.numeric_deriv(f, 'h', 5)
        for i in [2, 3, 6]:
            if not increment_filter[i, 0]:
                self.jacobian[k, i, 0] = self.numeric_deriv(f, 'm', i)
            if not all(increment_filter[i, 3:]):
                self.jacobian[k, i, 3:] = self.numeric_deriv(f, 'fluid', i)

        if self.P.is_var:
            self.jacobian[k, 7 + self.P.var_pos, 0] = (
                self.numeric_deriv(f, 'P', 7))

    def Qloss_char_func(self):
        r"""
        Calculate the relation of heat loss and thermal input.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                \begin{split}
                0 = & \dot{Q}_{loss} \cdot
                f_{TI}\left(\frac{P}{P_{design}}\right) \\
                & + LHV \cdot \left[\sum_i
                \left(\dot{m}_{in,i} \cdot x_{f,i}\right) -
                \dot{m}_{out,3} \cdot x_{f,3} \right] \cdot
                f_{QLOSS}\left(\frac{P}{P_{ref}}\right)\\
                \end{split}\\
                \forall i \in [3,4]
        """
        if np.isnan(self.P.design):
            expr = 1
        else:
            expr = self.P.val / self.P.design

        return (self.calc_ti() * self.Qloss_char.char_func.evaluate(expr) +
                self.tiP_char.char_func.evaluate(expr) * self.Qloss.val)

    def Qloss_char_func_doc(self, label):
        r"""
        Calculate the relation of heat loss and thermal input.

        Parameters
        ----------
        label : str
            Label for equation.
        """
        latex = (
            r'\begin{split}' + '\n'
            r'0=&LHV_{fuel} \cdot \left[\sum_i \left('
            r'\dot{m}_{\mathrm{in,}i} \cdot x_{fuel\mathrm{,in,}i}\right)'
            r'-\dot{m}_\mathrm{out,3}\cdot x_{fuel\mathrm{,out,}3}\right]'
            r' \cdot f_\mathrm{QLOSS}\left(\frac{P}{P_\mathrm{design}}'
            r'\right)\\' + '\n'
            r'&+\dot{Q}_\mathrm{loss} \cdot f_\mathrm{TI}'
            r'\left(\frac{P}{P_\mathrm{design}}'
            r'\right)\\ ' + '\n'
            r'&\forall i \in [3,4]\\ ' + '\n'
            r'\end{split}'
        )
        return generate_latex_eq(self, latex, label)

    def Qloss_char_deriv(self, increment_filter, k):
        """
        Calculate partial derivatives of heat loss to thermal input char.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of equation in Jacobian matrix.
        """
        f = self.Qloss_char_func
        for i in [2, 3, 6]:
            if not increment_filter[i, 0]:
                self.jacobian[k, i, 0] = self.numeric_deriv(f, 'm', i)
            if not all(increment_filter[i, 3:]):
                self.jacobian[k, i, 3:] = self.numeric_deriv(f, 'fluid', i)

        if self.P.is_var:
            self.jacobian[k, 7 + self.P.var_pos, 0] = (
                self.numeric_deriv(f, 'P', 7))
        if self.Qloss.is_var:
            self.jacobian[k, 7 + self.Qloss.var_pos, 0] = (
                self.numeric_deriv(f, 'Qloss', 7))

    # ti_func is in class CombustionChamber
    def ti_deriv(self, increment_filter, k):
        """
        Calculate partial derivatives of thermal input equation.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of equation in Jacobian matrix.
        """
        f = self.ti_func
        for i in [2, 3, 6]:
            self.jacobian[k, i, 0] = self.numeric_deriv(f, 'm', i)
            self.jacobian[k, i, 3:] = self.numeric_deriv(f, 'fluid', i)

    def calc_ti(self):
        r"""
        Calculate the thermal input of the combustion engine.

        Returns
        -------
        ti : float
            Thermal input.

            .. math::

                ti = LHV \cdot \left[\sum_i \left(\dot{m}_{in,i} \cdot x_{f,i}
                \right) - \dot{m}_{out,3} \cdot x_{f,3} \right]

                \forall i \in [3,4]
        """
        ti = 0
        for f in self.fuel_list:
            m = 0
            for i in self.inl[2:]:
                m += i.m.val_SI * i.fluid.val[f]

            for o in self.outl[2:]:
                m -= o.m.val_SI * o.fluid.val[f]

            ti += m * self.fuels[f]['LHV']

        return ti

    def calc_P(self):
        r"""
        Calculate the power output of the combustion engine.

        Returns
        -------
        P : float
            Power output.

            .. math::

                P = -\frac{LHV \cdot \dot{m}_{f}}
                {f_{TI}\left(\frac{P}{P_{ref}}\right)}

        """
        if np.isnan(self.P.design):
            expr = 1
        else:
            expr = self.P.val / self.P.design

        return -self.calc_ti() / self.tiP_char.char_func.evaluate(expr)

    def calc_Qloss(self):
        r"""
        Calculate the heat loss of the combustion engine.

        Returns
        -------
        Qloss : float
            Heat loss.

            .. math::

                \dot{Q}_{loss} = -\frac{LHV \cdot \dot{m}_{f} \cdot
                f_{QLOSS}\left(\frac{P}{P_{ref}}\right)}
                {f_{TI}\left(\frac{P}{P_{ref}}\right)}
        """
        if np.isnan(self.P.design):
            expr = 1
        else:
            expr = self.P.val / self.P.design

        return (-self.calc_ti() * self.Qloss_char.char_func.evaluate(expr) /
                self.tiP_char.char_func.evaluate(expr))

    def bus_func(self, bus):
        r"""
        Calculate the value of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus.Bus
            TESPy bus object.

        Returns
        -------
        residual : float
            Value of energy transfer :math:`\dot{E}`. This value is passed to
            :py:meth:`tespy.components.component.Component.calc_bus_value`
            for value manipulation according to the specified characteristic
            line of the bus.

            .. math::

                \dot{E} = \begin{cases}
                LHV \cdot \dot{m}_{f} & \text{key = 'TI'}\\
                P & \text{key = 'P'}\\
                -\dot{m}_1 \cdot \left( h_{1,out} - h_{1,in} \right)
                -\dot{m}_2 \cdot \left( h_{2,out} - h_{2,in} \right) &
                \text{key = 'Q'}\\
                -\dot{m}_1 \cdot \left( h_{1,out} - h_{1,in} \right) &
                \text{key = 'Q1'}\\
                -\dot{m}_2 \cdot \left( h_{2,out} - h_{2,in} \right) &
                \text{key = 'Q2'}\\
                \dot{Q}_{loss} & \text{key = 'Qloss'}
                \end{cases}
        """
        ######################################################################
        # value for bus parameter of thermal input (TI)
        if bus['param'] == 'TI':
            val = self.calc_ti()

        ######################################################################
        # value for bus parameter of power output (P)
        elif bus['param'] == 'P':
            val = self.calc_P()

        ######################################################################
        # value for bus parameter of total heat production (Q)
        elif bus['param'] == 'Q':
            val = 0
            for j in range(2):
                i = self.inl[j]
                o = self.outl[j]
                val -= i.m.val_SI * (o.h.val_SI - i.h.val_SI)

        ######################################################################
        # value for bus parameter of heat production 1 (Q1)
        elif bus['param'] == 'Q1':
            i = self.inl[0]
            o = self.outl[0]
            val = -i.m.val_SI * (o.h.val_SI - i.h.val_SI)

        ######################################################################
        # value for bus parameter of heat production 2 (Q2)
        elif bus['param'] == 'Q2':
            i = self.inl[1]
            o = self.outl[1]
            val = -i.m.val_SI * (o.h.val_SI - i.h.val_SI)

        ######################################################################
        # value for bus parameter of heat loss (Qloss)
        elif bus['param'] == 'Qloss':
            val = self.calc_Qloss()

        ######################################################################
        # missing/invalid bus parameter
        else:
            msg = ('The parameter ' + str(bus['param']) +
                   ' is not a valid parameter for a ' + self.component() + '.')
            logging.error(msg)
            raise ValueError(msg)

        return val

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
        ######################################################################
        # value for bus parameter of thermal input (TI)
        if bus['param'] == 'TI':
            return CombustionChamber.bus_func_doc(self, bus)

        ######################################################################
        # value for bus parameter of power output (P)
        elif bus['param'] == 'P':
            return 'P'

        ######################################################################
        # value for bus parameter of total heat production (Q)
        elif bus['param'] == 'Q':
            return (
                r'-\dot{m}_\mathrm{in,1} \cdot \left( h_\mathrm{out,1} -'
                r'h_\mathrm{in,1} \right) - \dot{m}_\mathrm{in,2} \cdot '
                r'\left( h_\mathrm{out,2} - h_\mathrm{in,2} \right)')

        ######################################################################
        # value for bus parameter of heat production 1 (Q1)
        elif bus['param'] == 'Q1':
            return (
                r'-\dot{m}_\mathrm{in,1} \cdot \left( h_\mathrm{out,1} -'
                r'h_\mathrm{in,1} \right)')

        ######################################################################
        # value for bus parameter of heat production 2 (Q2)
        elif bus['param'] == 'Q2':
            return (
                r'- \dot{m}_\mathrm{in,2} \cdot '
                r'\left( h_\mathrm{out,2} - h_\mathrm{in,2} \right)')

        ######################################################################
        # value for bus parameter of heat loss (Qloss)
        elif bus['param'] == 'Qloss':
            return r'\dot{Q}_\mathrm{loss}'

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
        deriv = np.zeros((1, 7 + self.num_vars, self.num_nw_vars))
        f = self.calc_bus_value
        b = bus.comps.loc[self]

        ######################################################################
        # derivatives for bus parameter of thermal input (TI)
        if b['param'] == 'TI':
            for i in [2, 3, 6]:
                deriv[0, i, 0] = self.numeric_deriv(f, 'm', i, bus=bus)
                deriv[0, i, 3:] = self.numeric_deriv(f, 'fluid', i, bus=bus)

        ######################################################################
        # derivatives for bus parameter of power production (P) or
        # heat loss (Qloss)
        elif b['param'] == 'P' or b['param'] == 'Qloss':
            for i in [2, 3, 6]:
                deriv[0, i, 0] = self.numeric_deriv(f, 'm', i, bus=bus)
                deriv[0, i, 3:] = self.numeric_deriv(f, 'fluid', i, bus=bus)

            # variable power
            if self.P.is_var:
                deriv[0, 7 + self.P.var_pos, 0] = (
                    self.numeric_deriv(f, 'P', 7, bus=bus))

        ######################################################################
        # derivatives for bus parameter of total heat production (Q)
        elif b['param'] == 'Q':
            for i in range(2):
                deriv[0, i, 0] = self.numeric_deriv(f, 'm', i, bus=bus)
                deriv[0, i, 2] = self.numeric_deriv(f, 'h', i, bus=bus)
                deriv[0, i + 4, 2] = self.numeric_deriv(f, 'h', i + 4, bus=bus)

        ######################################################################
        # derivatives for bus parameter of heat production 1 (Q1)
        elif b['param'] == 'Q1':
            deriv[0, 0, 0] = self.numeric_deriv(f, 'm', 0, bus=bus)
            deriv[0, 0, 2] = self.numeric_deriv(f, 'h', 0, bus=bus)
            deriv[0, 4, 2] = self.numeric_deriv(f, 'h', 4, bus=bus)

        ######################################################################
        # derivatives for bus parameter of heat production 2 (Q2)
        elif b['param'] == 'Q2':
            deriv[0, 1, 0] = self.numeric_deriv(f, 'm', 1, bus=bus)
            deriv[0, 1, 2] = self.numeric_deriv(f, 'h', 1, bus=bus)
            deriv[0, 5, 2] = self.numeric_deriv(f, 'h', 5, bus=bus)

        ######################################################################
        # missing/invalid bus parameter
        else:
            msg = ('The parameter ' + str(b['param']) +
                   ' is not a valid parameter for a ' + self.component() + '.')
            logging.error(msg)
            raise ValueError(msg)

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
            m_h2o += (n_fuel * self.fuels[f]['H'] /
                      2 * molar_masses[self.h2o] * fact_fuel[f])
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

        o = self.outl[2]
        for fluid, x in o.fluid.val.items():
            if not o.fluid.val_set[fluid] and fluid in fg:
                o.fluid.val[fluid] = fg[fluid]
        o.target.propagate_fluid_to_target(o, o.target)

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
                5 \cdot 10^5 & \text{key = 'p'}\\
                5 \cdot 10^5 & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 5e5
        elif key == 'h':
            return 5e5

    def propagate_fluid_to_target(self, inconn, start):
        r"""
        Propagate the fluids towards connection's target in recursion.

        Parameters
        ----------
        inconn : tespy.connections.connection.Connection
            Connection to initialise.

        start : tespy.components.component.Component
            This component is the fluid propagation starting point.
            The starting component is saved to prevent infinite looping.
        """
        for outconn in self.outl[:2]:
            for fluid, x in inconn.fluid.val.items():
                if (outconn.fluid.val_set[fluid] is False and
                        outconn.good_starting_values is False):
                    outconn.fluid.val[fluid] = x
            outconn.target.propagate_fluid_to_target(outconn, start)

    def propagate_fluid_to_source(self, outconn, start):
        r"""
        Propagate the fluids towards connection's source in recursion.

        Parameters
        ----------
        outconn : tespy.connections.connection.Connection
            Connection to initialise.

        start : tespy.components.component.Component
            This component is the fluid propagation starting point.
            The starting component is saved to prevent infinite looping.
        """
        for inconn in self.inl[:2]:
            for fluid, x in outconn.fluid.val.items():
                if (inconn.fluid.val_set[fluid] is False and
                        inconn.good_starting_values is False):
                    inconn.fluid.val[fluid] = x

            inconn.source.propagate_fluid_to_source(inconn, start)

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        # Q, pr and zeta
        for i in range(2):
            self.get_attr('Q' + str(i + 1)).val = -self.inl[i].m.val_SI * (
                self.outl[i].h.val_SI - self.inl[i].h.val_SI)
            self.get_attr('pr' + str(i + 1)).val = (
                self.outl[i].p.val_SI / self.inl[i].p.val_SI)
            self.get_attr('zeta' + str(i + 1)).val = (
                (self.inl[i].p.val_SI - self.outl[i].p.val_SI) * np.pi ** 2 / (
                    4 * self.inl[i].m.val_SI ** 2 *
                    (self.inl[i].vol.val_SI + self.outl[i].vol.val_SI)
                ))

        self.P.val = self.calc_P()
        self.Qloss.val = self.calc_Qloss()

        CombustionChamber.calc_parameters(self)

    def check_parameter_bounds(self):
        r"""Check parameter value limits."""
        Component.check_parameter_bounds(self)
        # get bound errors for characteristic lines
        if np.isnan(self.P.design):
            expr = 1
        else:
            expr = self.P.val / self.P.design
        self.tiP_char.char_func.get_domain_errors(expr, self.label)
        self.Qloss_char.char_func.get_domain_errors(expr, self.label)
        self.Q1_char.char_func.get_domain_errors(expr, self.label)
        self.Q2_char.char_func.get_domain_errors(expr, self.label)

    def entropy_balance(self):
        r"""
        Calculate entropy balance of combustion engine.

        For the entropy balance of a combustion engine two additional
        parameters need to be specified:

        - virtual inner temperature :code:`T_v_inner` that is used to determine
          the entropy of heat transferred from the hot side.
        - mechanical efficiency :code:`eta_mech` describing the ratio of power
          output :code:`P` to reversible power of the motor
          :cite:`Zahoransky2019`. It is used to determine the irreversibilty
          inside the motor.

          .. math::

              P_\mathrm{irr,inner}=\left(1 - \frac{1}{\eta_\mathrm{mech}}
              \right) \cdot P

        The default values are:

        - :code:`T_v_inner`: flue gas temperature (result of calculation)
        - :code:`eta_mech`: 0.85

        Note
        ----
        The entropy balance makes the following parameter available:

        - :code:`T_mcomb`: Thermodynamic temperature of heat of combustion
        - :code:`S_comb`: Entropy production due to combustion
        - :code:`T_mQ1`: Thermodynamic temperature of heat at cold side of
          heater 1
        - :code:`S_Q11`: Entropy transport at hot side of heater 1
        - :code:`S_Q12`: Entropy transport at cold side of heater 1
        - :code:`S_Q1irr`: Entropy production due to heat transfer at heater 1
        - :code:`S_irr1`: Entropy production due to pressure losses at heater 1
        - :code:`T_mQ2`: Thermodynamic temperature of heat at cold side of
          heater 2
        - :code:`S_Q21`: Entropy transport at hot side of heater 2
        - :code:`S_Q22`: Entropy transport at cold side of heater 2
        - :code:`S_Q2irr`: Entropy production due to heat transfer at heater 2
        - :code:`S_irr2`: Entropy production due to pressure losses at heater 2
        - :code:`S_irr_i`: Entropy production due to internal irreversibilty
        - :code:`S_Qloss`: Entropy transport with heat loss to ambient
        - :code:`S_Qcomb`: Virtual entropy transport of heat to revert
          combustion gases to reference state
        - :code:`S_irr`: Total entropy production due to irreversibilty

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

            \begin{split}
            T_\mathrm{m,comb}= & \frac{\dot{m}_\mathrm{fuel} \cdot LHV}
            {\dot{S}_\mathrm{comb}}\\
            \dot{S}_\mathrm{comb} =&\dot{S}_\mathrm{Q,comb}-\left(
            \dot{S}_\mathrm{Q,11} + \dot{S}_\mathrm{Q,21} +
            \dot{S}_\mathrm{Q,loss} +\dot{S}_\mathrm{irr,i}\right)\\
            \dot{S}_\mathrm{Q,comb}= & \dot{m}_\mathrm{fluegas} \cdot
            \left(s_\mathrm{fluegas}-s_\mathrm{fluegas,ref}\right)\\
            & - \sum_{i=3}^4 \dot{m}_{\mathrm{in,}i} \cdot
            \left( s_{\mathrm{in,}i} - s_{\mathrm{in,ref,}i} \right)\\
            \dot{S}_\mathrm{Q,11}= & \frac{\dot{Q}_1}{T_\mathrm{v,inner}}\\
            \dot{S}_\mathrm{Q,21}= & \frac{\dot{Q}_2}{T_\mathrm{v,inner}}\\
            \dot{S}_\mathrm{Q,loss}= & \frac{\dot{Q}_\mathrm{loss}}
            {T_\mathrm{v,inner}}\\
            \dot{S}_\mathrm{irr,i}= & \frac{\left(1 -
            \frac{1}{\eta_\mathrm{mech}}\right) \cdot P}{T_\mathrm{v,inner}}\\
            T_\mathrm{Q,12} = &\frac{-\dot{Q}_1}{\dot{m}_1 \cdot \left(
            s_\mathrm{out,1} - s_\mathrm{in,1}\right)}\\
            T_\mathrm{Q,22} = &\frac{-\dot{Q}_2}{\dot{m}_2 \cdot \left(
            s_\mathrm{out,2} - s_\mathrm{in,2}\right)}\\
            \dot{S}_\mathrm{irr} = &\sum \dot{S}_\mathrm{irr}\\
            \end{split}\\
        """
        T_ref = 298.15
        p_ref = 1e5
        o = self.outl[2]
        self.S_Qcomb = o.m.val_SI * (
            o.s.val_SI -
            s_mix_pT([0, p_ref, 0, o.fluid.val], T_ref, force_gas=True))

        for c in self.inl[2:]:
            self.S_Qcomb -= c.m.val_SI * (
                c.s.val_SI -
                s_mix_pT([0, p_ref, 0, c.fluid.val], T_ref, force_gas=True))

        # (virtual) thermodynamic temperature of combustion, use default value
        # if not specified
        if not self.T_v_inner.is_set:
            self.T_v_inner.val = o.T.val_SI

        for i in range(2):
            inl = self.inl[i]
            out = self.outl[i]
            p_star = inl.p.val_SI * (
                self.get_attr('pr' + str(i + 1)).val) ** 0.5
            s_i_star = s_mix_ph(
                [0, p_star, inl.h.val_SI, inl.fluid.val], T0=inl.T.val_SI)
            s_o_star = s_mix_ph(
                [0, p_star, out.h.val_SI, out.fluid.val], T0=out.T.val_SI)

            setattr(self, 'S_Q' + str(i + 1) + '2',
                    inl.m.val_SI * (s_o_star - s_i_star))
            S_Q = self.get_attr('S_Q' + str(i + 1) + '2')
            setattr(self, 'S_irr' + str(i + 1),
                    inl.m.val_SI * (out.s.val_SI - inl.s.val_SI) - S_Q)
            setattr(self, 'T_mQ' + str(i + 1),
                    inl.m.val_SI * (out.h.val_SI - inl.h.val_SI) / S_Q)

        # internal irreversibilty
        self.P_irr_i = (1 / self.eta_mech.val - 1) * self.P.val

        # internal entropy flow and production
        self.S_Q11 = self.Q1.val / self.T_v_inner.val
        self.S_Q21 = self.Q2.val / self.T_v_inner.val
        self.S_Qloss = self.Qloss.val / self.T_v_inner.val
        self.S_irr_i = self.P_irr_i / self.T_v_inner.val

        # entropy production of heaters due to heat transfer
        self.S_Q1irr = self.S_Q12 - self.S_Q11
        self.S_Q2irr = self.S_Q22 - self.S_Q21

        # calculate entropy production of combustion
        self.S_comb = (
            self.S_Qcomb - self.S_Q11 - self.S_Q21 - self.S_Qloss -
            self.S_irr_i)

        # thermodynamic temperature of heat input
        self.T_mcomb = self.calc_ti() / self.S_comb
        # total irreversibilty production
        self.S_irr = (
            self.S_irr_i + self.S_irr2 + self.S_irr1 + self.S_Q1irr +
            self.S_Q2irr)
