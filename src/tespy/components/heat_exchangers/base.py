# -*- coding: utf-8

"""Module of class HeatExchanger.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/heat_exchangers/heat_exchanger.py

SPDX-License-Identifier: MIT
"""
import numpy as np

from tespy.components.component import Component
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentCharacteristics as dc_gcc
from tespy.tools.document_models import generate_latex_eq
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.fluid_properties import s_mix_ph


class HeatExchanger(Component):
    r"""
    Class for counter current heat exchanger.

    The component HeatExchanger is the parent class for the components:

    - :py:class:`tespy.components.heat_exchangers.condenser.Condenser`
    - :py:class:`tespy.components.heat_exchangers.desuperheater.Desuperheater`

    **Mandatory Equations**

    - :py:meth:`tespy.components.component.Component.fluid_func`
    - :py:meth:`tespy.components.component.Component.mass_flow_func`
    - :py:meth:`tespy.components.heat_exchangers.heat_exchanger.HeatExchanger.energy_balance_func`

    **Optional Equations**

    - :py:meth:`tespy.components.heat_exchangers.heat_exchanger.HeatExchanger.energy_balance_hot_func`
    - :py:meth:`tespy.components.heat_exchangers.heat_exchanger.HeatExchanger.kA_func`
    - :py:meth:`tespy.components.heat_exchangers.heat_exchanger.HeatExchanger.kA_char_func`
    - :py:meth:`tespy.components.heat_exchangers.heat_exchanger.HeatExchanger.ttd_u_func`
    - :py:meth:`tespy.components.heat_exchangers.heat_exchanger.HeatExchanger.ttd_l_func`
    - hot side :py:meth:`tespy.components.component.Component.pr_func`
    - cold side :py:meth:`tespy.components.component.Component.pr_func`
    - hot side :py:meth:`tespy.components.component.Component.zeta_func`
    - cold side :py:meth:`tespy.components.component.Component.zeta_func`

    Inlets/Outlets

    - in1, in2 (index 1: hot side, index 2: cold side)
    - out1, out2 (index 1: hot side, index 2: cold side)

    Image

    .. image:: _images/HeatExchanger.svg
       :alt: alternative text
       :align: center

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

    Q : float, dict
        Heat transfer, :math:`Q/\text{W}`.

    pr1 : float, dict, :code:`"var"`
        Outlet to inlet pressure ratio at hot side, :math:`pr/1`.

    pr2 : float, dict, :code:`"var"`
        Outlet to inlet pressure ratio at cold side, :math:`pr/1`.

    zeta1 : float, dict, :code:`"var"`
        Geometry independent friction coefficient at hot side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    zeta2 : float, dict, :code:`"var"`
        Geometry independent friction coefficient at cold side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    ttd_l : float, dict
        Lower terminal temperature difference :math:`ttd_\mathrm{l}/\text{K}`.

    ttd_u : float, dict
        Upper terminal temperature difference :math:`ttd_\mathrm{u}/\text{K}`.

    kA : float, dict
        Area independent heat transfer coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char : dict
        Area independent heat transfer coefficient characteristic.

    kA_char1 : tespy.tools.characteristics.CharLine, dict
        Characteristic line for hot side heat transfer coefficient.

    kA_char2 : tespy.tools.characteristics.CharLine, dict
        Characteristic line for cold side heat transfer coefficient.

    Note
    ----
    The HeatExchanger and subclasses (
    :py:class:`tespy.components.heat_exchangers.condenser.Condenser`,
    :py:class:`tespy.components.heat_exchangers.desuperheater.Desuperheater`)
    are countercurrent heat exchangers. Equations (:code:`kA`, :code:`ttd_u`,
    :code:`ttd_l`) do not work for directcurrent and crosscurrent or
    combinations of different types.

    Example
    -------
    A water cooling is installed to transfer heat from hot exhaust air. The
    heat exchanger is designed for a terminal temperature difference of 5 K.
    From this, it is possible to calculate the heat transfer coefficient and
    predict water and air outlet temperature in offdesign operation.

    >>> from tespy.components import Sink, Source, HeatExchanger
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> from tespy.tools import document_model
    >>> import shutil
    >>> nw = Network(fluids=['water', 'air'], T_unit='C', p_unit='bar',
    ... h_unit='kJ / kg', iterinfo=False)
    >>> exhaust_hot = Source('Exhaust air outlet')
    >>> exhaust_cold = Sink('Exhaust air inlet')
    >>> cw_cold = Source('cooling water inlet')
    >>> cw_hot = Sink('cooling water outlet')
    >>> he = HeatExchanger('waste heat exchanger')
    >>> he.component()
    'heat exchanger'
    >>> ex_he = Connection(exhaust_hot, 'out1', he, 'in1')
    >>> he_ex = Connection(he, 'out1', exhaust_cold, 'in1')
    >>> cw_he = Connection(cw_cold, 'out1', he, 'in2')
    >>> he_cw = Connection(he, 'out2', cw_hot, 'in1')
    >>> nw.add_conns(ex_he, he_ex, cw_he, he_cw)

    The volumetric flow of the air is at 100 l/s. After designing the component
    it is possible to predict the temperature at different flow rates or
    different inlet temperatures of the exhaust air.

    >>> he.set_attr(pr1=0.98, pr2=0.98, ttd_u=5,
    ... design=['pr1', 'pr2', 'ttd_u'], offdesign=['zeta1', 'zeta2', 'kA_char'])
    >>> cw_he.set_attr(fluid={'air': 0, 'water': 1}, T=10, p=3,
    ... offdesign=['m'])
    >>> ex_he.set_attr(fluid={'air': 1, 'water': 0}, v=0.1, T=35)
    >>> he_ex.set_attr(T=17.5, p=1, design=['T'])
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(ex_he.T.val - he_cw.T.val, 0)
    5.0
    >>> ex_he.set_attr(v=0.075)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(he_cw.T.val, 1)
    27.5
    >>> round(he_ex.T.val, 1)
    14.4
    >>> ex_he.set_attr(v=0.1, T=40)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> document_model(nw)
    >>> round(he_cw.T.val, 1)
    33.9
    >>> round(he_ex.T.val, 1)
    18.8
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    @staticmethod
    def component():
        return 'heat exchanger'

    def get_variables(self):
        return {
            'Q': dc_cp(
                max_val=0, func=self.energy_balance_hot_func, num_eq=1,
                deriv=self.energy_balance_hot_deriv,
                latex=self.energy_balance_hot_func_doc),
            'kA': dc_cp(
                min_val=0, num_eq=1, func=self.kA_func, latex=self.kA_func_doc,
                deriv=self.kA_deriv),
            'td_log': dc_cp(min_val=0, is_result=True),
            'ttd_u': dc_cp(
                min_val=0, num_eq=1, func=self.ttd_u_func,
                deriv=self.ttd_u_deriv, latex=self.ttd_u_func_doc),
            'ttd_l': dc_cp(
                min_val=0, num_eq=1, func=self.ttd_l_func,
                deriv=self.ttd_l_deriv, latex=self.ttd_l_func_doc),
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
            'kA_char': dc_gcc(
                elements=['kA_char1', 'kA_char2'],
                num_eq=1, latex=self.kA_char_func_doc, func=self.kA_char_func,
                deriv=self.kA_char_deriv),
            'kA_char1': dc_cc(param='m'),
            'kA_char2': dc_cc(
                param='m', char_params={
                    'type': 'rel', 'inconn': 1, 'outconn': 1})
        }

    def get_mandatory_constraints(self):
        return {
            'mass_flow_constraints': {
                'func': self.mass_flow_func, 'deriv': self.mass_flow_deriv,
                'constant_deriv': True, 'latex': self.mass_flow_func_doc,
                'num_eq': 2},
            'fluid_constraints': {
                'func': self.fluid_func, 'deriv': self.fluid_deriv,
                'constant_deriv': True, 'latex': self.fluid_func_doc,
                'num_eq': self.num_nw_fluids * 2},
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
        return ['out1', 'out2']

    def energy_balance_func(self):
        r"""
        Equation for heat exchanger energy balance.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_{in,1} \cdot \left(h_{out,1} - h_{in,1} \right) +
                \dot{m}_{in,2} \cdot \left(h_{out,2} - h_{in,2} \right)
        """
        return (
            self.inl[0].m.val_SI * (
                self.outl[0].h.val_SI - self.inl[0].h.val_SI) +
            self.inl[1].m.val_SI * (
                self.outl[1].h.val_SI - self.inl[1].h.val_SI))

    def energy_balance_func_doc(self, label):
        r"""
        Equation for heat exchanger energy balance.

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
            r'0 = \dot{m}_\mathrm{in,1} \cdot \left(h_\mathrm{out,1} -'
            r' h_\mathrm{in,1} \right) +\dot{m}_\mathrm{in,2} \cdot '
            r'\left(h_\mathrm{out,2} - h_\mathrm{in,2} \right)')
        return generate_latex_eq(self, latex, label)

    def energy_balance_deriv(self, increment_filter, k):
        r"""
        Partial derivatives of energy balance function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        for i in range(2):
            self.jacobian[k, i, 0] = (
                self.outl[i].h.val_SI - self.inl[i].h.val_SI)
            self.jacobian[k, i, 2] = -self.inl[i].m.val_SI

        self.jacobian[k, 2, 2] = self.inl[0].m.val_SI
        self.jacobian[k, 3, 2] = self.inl[1].m.val_SI

    def energy_balance_hot_func(self):
        r"""
        Equation for hot side heat exchanger energy balance.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 =\dot{m}_{in,1} \cdot \left(h_{out,1}-h_{in,1}\right)-\dot{Q}
        """
        return self.inl[0].m.val_SI * (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI) - self.Q.val

    def energy_balance_hot_func_doc(self, label):
        r"""
        Equation for hot side heat exchanger energy balance.

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
            r'0 =\dot{m}_{in,1} \cdot \left(h_{out,1}-'
            r'h_{in,1}\right)-\dot{Q}')
        return generate_latex_eq(self, latex, label)

    def energy_balance_hot_deriv(self, increment_filter, k):
        r"""
        Partial derivatives for hot side heat exchanger energy balance.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        self.jacobian[k, 0, 0] = self.outl[0].h.val_SI - self.inl[0].h.val_SI
        self.jacobian[k, 0, 2] = -self.inl[0].m.val_SI
        self.jacobian[k, 2, 2] = self.inl[0].m.val_SI

    def kA_func(self):
        r"""
        Calculate heat transfer from heat transfer coefficient.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_{in,1} \cdot \left( h_{out,1} - h_{in,1}\right) +
                kA \cdot \frac{T_{out,1} -
                T_{in,2} - T_{in,1} + T_{out,2}}
                {\ln{\frac{T_{out,1} - T_{in,2}}{T_{in,1} - T_{out,2}}}}
        """
        i1 = self.inl[0]
        i2 = self.inl[1]
        o1 = self.outl[0]
        o2 = self.outl[1]

        # temperature value manipulation for convergence stability
        T_i1 = T_mix_ph(i1.get_flow(), T0=i1.T.val_SI)
        T_i2 = T_mix_ph(i2.get_flow(), T0=i2.T.val_SI)
        T_o1 = T_mix_ph(o1.get_flow(), T0=o1.T.val_SI)
        T_o2 = T_mix_ph(o2.get_flow(), T0=o2.T.val_SI)

        if T_i1 <= T_o2:
            T_i1 = T_o2 + 0.01
        if T_i1 <= T_o2:
            T_o2 = T_i1 - 0.01
        if T_i1 <= T_o2:
            T_o1 = T_i2 + 0.02
        if T_o1 <= T_i2:
            T_i2 = T_o1 - 0.02

        td_log = ((T_o1 - T_i2 - T_i1 + T_o2) /
                  np.log((T_o1 - T_i2) / (T_i1 - T_o2)))

        return i1.m.val_SI * (
            o1.h.val_SI - i1.h.val_SI) + self.kA.val * td_log

    def kA_func_doc(self, label):
        r"""
        Calculate heat transfer from heat transfer coefficient.

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
            r'0 = \dot{m}_\mathrm{in,1} \cdot \left( h_\mathrm{out,1} - '
            r'h_\mathrm{in,1}\right)+ kA \cdot \frac{T_\mathrm{out,1} - '
            r'T_\mathrm{in,2} - T_\mathrm{in,1} + T_\mathrm{out,2}}'
            r'{\ln{\frac{T_\mathrm{out,1} - T_\mathrm{in,2}}'
            r'{T_\mathrm{in,1} - T_\mathrm{out,2}}}}'
        )
        return generate_latex_eq(self, latex, label)

    def kA_deriv(self, increment_filter, k):
        r"""
        Partial derivatives of heat transfer coefficient function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        f = self.kA_func
        self.jacobian[k, 0, 0] = self.outl[0].h.val_SI - self.inl[0].h.val_SI
        for i in range(4):
            if not increment_filter[i, 1]:
                self.jacobian[k, i, 1] = self.numeric_deriv(f, 'p', i)
            if not increment_filter[i, 2]:
                self.jacobian[k, i, 2] = self.numeric_deriv(f, 'h', i)

    def kA_char_func(self):
        r"""
        Calculate heat transfer from heat transfer coefficient characteristic.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_{in,1} \cdot \left( h_{out,1} - h_{in,1}\right) +
                kA_{design} \cdot f_{kA} \cdot \frac{T_{out,1} -
                T_{in,2} - T_{in,1} + T_{out,2}}
                {\ln{\frac{T_{out,1} - T_{in,2}}{T_{in,1} - T_{out,2}}}}

                f_{kA} = \frac{2}{\frac{1}{f_1\left( expr_1\right)} +
                \frac{1}{f_2\left( expr_2\right)}}

        Note
        ----
        For standard functions f\ :subscript:`1` \ and f\ :subscript:`2` \ see
        module :py:mod:`tespy.data`.
        """
        p1 = self.kA_char1.param
        p2 = self.kA_char2.param
        f1 = self.get_char_expr(p1, **self.kA_char1.char_params)
        f2 = self.get_char_expr(p2, **self.kA_char2.char_params)

        i1 = self.inl[0]
        i2 = self.inl[1]
        o1 = self.outl[0]
        o2 = self.outl[1]

        # temperature value manipulation for convergence stability
        T_i1 = T_mix_ph(i1.get_flow(), T0=i1.T.val_SI)
        T_i2 = T_mix_ph(i2.get_flow(), T0=i2.T.val_SI)
        T_o1 = T_mix_ph(o1.get_flow(), T0=o1.T.val_SI)
        T_o2 = T_mix_ph(o2.get_flow(), T0=o2.T.val_SI)

        if T_i1 <= T_o2:
            T_i1 = T_o2 + 0.01
        if T_i1 <= T_o2:
            T_o2 = T_i1 - 0.01
        if T_i1 <= T_o2:
            T_o1 = T_i2 + 0.02
        if T_o1 <= T_i2:
            T_i2 = T_o1 - 0.02

        td_log = ((T_o1 - T_i2 - T_i1 + T_o2) /
                  np.log((T_o1 - T_i2) / (T_i1 - T_o2)))

        fkA1 = self.kA_char1.char_func.evaluate(f1)
        fkA2 = self.kA_char2.char_func.evaluate(f2)
        fkA = 2 / (1 / fkA1 + 1 / fkA2)

        return (
            i1.m.val_SI * (o1.h.val_SI - i1.h.val_SI) +
            self.kA.design * fkA * td_log)

    def kA_char_func_doc(self, label):
        r"""
        Calculate heat transfer from heat transfer coefficient characteristic.

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
            r'0 = & \dot{m}_\mathrm{in,1} \cdot \left( h_\mathrm{out,1} - '
            r'h_\mathrm{in,1}\right)\\' + '\n'
            r'&+kA_\mathrm{design} \cdot '
            r'f_\mathrm{kA} \cdot \frac{T_\mathrm{out,1} - T_\mathrm{in,2}'
            r' - T_\mathrm{in,1} + T_\mathrm{out,2}}{\ln{'
            r'\frac{T_\mathrm{out,1} - T_\mathrm{in,2}}{T_\mathrm{in,1} -'
            r' T_\mathrm{out,2}}}}\\' + '\n'
            r'f_\mathrm{kA}=&\frac{2}{\frac{1}{f\left(X_1\right)}+'
            r'\frac{1}{f\left(X_2\right)}}\\' + '\n'
            r'\end{split}'
        )
        return generate_latex_eq(self, latex, label)

    def kA_char_deriv(self, increment_filter, k):
        r"""
        Partial derivatives of heat transfer coefficient characteristic.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        f = self.kA_char_func
        if not increment_filter[0, 0]:
            self.jacobian[k, 0, 0] = self.numeric_deriv(f, 'm', 0)
        if not increment_filter[1, 0]:
            self.jacobian[k, 1, 0] = self.numeric_deriv(f, 'm', 1)
        for i in range(4):
            if not increment_filter[i, 1]:
                self.jacobian[k, i, 1] = self.numeric_deriv(f, 'p', i)
            if not increment_filter[i, 2]:
                self.jacobian[k, i, 2] = self.numeric_deriv(f, 'h', i)

    def ttd_u_func(self):
        r"""
        Equation for upper terminal temperature difference.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = ttd_{u} - T_{in,1} + T_{out,2}
        """
        T_i1 = T_mix_ph(self.inl[0].get_flow(), T0=self.inl[0].T.val_SI)
        T_o2 = T_mix_ph(self.outl[1].get_flow(), T0=self.outl[1].T.val_SI)
        return self.ttd_u.val - T_i1 + T_o2

    def ttd_u_func_doc(self, label):
        r"""
        Equation for upper terminal temperature difference.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        latex = r'0 = ttd_\mathrm{u} - T_\mathrm{in,1} + T_\mathrm{out,2}'
        return generate_latex_eq(self, latex, label)

    def ttd_u_deriv(self, increment_filter, k):
        """
        Calculate partial derivates of upper terminal temperature function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        f = self.ttd_u_func
        for i in [0, 3]:
            if not increment_filter[i, 1]:
                self.jacobian[k, i, 1] = self.numeric_deriv(f, 'p', i)
            if not increment_filter[i, 2]:
                self.jacobian[k, i, 2] = self.numeric_deriv(f, 'h', i)

    def ttd_l_func(self):
        r"""
        Equation for upper terminal temperature difference.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = ttd_{l} - T_{out,1} + T_{in,2}
        """
        T_i2 = T_mix_ph(self.inl[1].get_flow(), T0=self.inl[1].T.val_SI)
        T_o1 = T_mix_ph(self.outl[0].get_flow(), T0=self.outl[0].T.val_SI)
        return self.ttd_l.val - T_o1 + T_i2

    def ttd_l_func_doc(self, label):
        r"""
        Equation for upper terminal temperature difference.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        latex = r'0 = ttd_\mathrm{l} - T_\mathrm{out,1} + T_\mathrm{in,2}'
        return generate_latex_eq(self, latex, label)

    def ttd_l_deriv(self, increment_filter, k):
        """
        Calculate partial derivates of upper terminal temperature function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        f = self.ttd_l_func
        for i in [1, 2]:
            if not increment_filter[i, 1]:
                self.jacobian[k, i, 1] = self.numeric_deriv(f, 'p', i)
            if not increment_filter[i, 2]:
                self.jacobian[k, i, 2] = self.numeric_deriv(f, 'h', i)

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

                \dot{E} = \dot{m}_{in,1} \cdot \left(
                h_{out,1} - h_{in,1} \right)
        """
        return self.inl[0].m.val_SI * (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI)

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
        return (
            r'\dot{m}_\mathrm{in,1} \cdot \left(h_\mathrm{out,1} - '
            r'h_\mathrm{in,1} \right)')

    def bus_deriv(self, bus):
        r"""
        Calculate partial derivatives of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus.Bus
            TESPy bus object.

        Returns
        -------
        deriv : ndarray
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 4, self.num_nw_vars))
        f = self.calc_bus_value
        deriv[0, 0, 0] = self.numeric_deriv(f, 'm', 0, bus=bus)
        deriv[0, 0, 2] = self.numeric_deriv(f, 'h', 0, bus=bus)
        deriv[0, 2, 2] = self.numeric_deriv(f, 'h', 2, bus=bus)
        return deriv

    def initialise_source(self, c, key):
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
                4 \cdot 10^5 & \text{key = 'p'}\\
                h\left(p, 200 \text{K} \right) & \text{key = 'h' at outlet 1}\\
                h\left(p, 250 \text{K} \right) & \text{key = 'h' at outlet 2}
                \end{cases}
        """
        if key == 'p':
            return 50e5
        elif key == 'h':
            if c.source_id == 'out1':
                T = 200 + 273.15
                return h_mix_pT(c.get_flow(), T)
            else:
                T = 250 + 273.15
                return h_mix_pT(c.get_flow(), T)

    def initialise_target(self, c, key):
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
                4 \cdot 10^5 & \text{key = 'p'}\\
                h\left(p, 300 \text{K} \right) & \text{key = 'h' at inlet 1}\\
                h\left(p, 220 \text{K} \right) & \text{key = 'h' at outlet 2}
                \end{cases}
        """
        if key == 'p':
            return 50e5
        elif key == 'h':
            if c.target_id == 'in1':
                T = 300 + 273.15
                return h_mix_pT(c.get_flow(), T)
            else:
                T = 220 + 273.15
                return h_mix_pT(c.get_flow(), T)

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        # component parameters
        self.Q.val = self.inl[0].m.val_SI * (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI)
        self.ttd_u.val = self.inl[0].T.val_SI - self.outl[1].T.val_SI
        self.ttd_l.val = self.outl[0].T.val_SI - self.inl[1].T.val_SI

        # pr and zeta
        for i in range(2):
            self.get_attr('pr' + str(i + 1)).val = (
                self.outl[i].p.val_SI / self.inl[i].p.val_SI)
            self.get_attr('zeta' + str(i + 1)).val = (
                (self.inl[i].p.val_SI - self.outl[i].p.val_SI) * np.pi ** 2 / (
                    4 * self.inl[i].m.val_SI ** 2 *
                    (self.inl[i].vol.val_SI + self.outl[i].vol.val_SI)
                ))

        # kA and logarithmic temperature difference
        if self.ttd_u.val < 0 or self.ttd_l.val < 0:
            self.td_log.val = np.nan
            self.kA.val = np.nan
        else:
            self.td_log.val = ((self.ttd_l.val - self.ttd_u.val) /
                               np.log(self.ttd_l.val / self.ttd_u.val))
            self.kA.val = -self.Q.val / self.td_log.val

    def entropy_balance(self):
        r"""
        Calculate entropy balance of a heat exchanger.

        The allocation of the entropy streams due to heat exchanged and due to
        irreversibility is performed by solving for T on both sides of the heat
        exchanger:

        .. math::

            h_\mathrm{out} - h_\mathrm{in} = \int_\mathrm{in}^\mathrm{out} v
            \cdot dp - \int_\mathrm{in}^\mathrm{out} T \cdot ds

        As solving :math:`\int_\mathrm{in}^\mathrm{out} v \cdot dp` for non
        isobaric processes would require perfect process knowledge (the path)
        on how specific volume and pressure change throught the component, the
        heat transfer is splitted into three separate virtual processes for
        both sides:

        - in->in*: decrease pressure to
          :math:`p_\mathrm{in*}=p_\mathrm{in}\cdot\sqrt{\frac{p_\mathrm{out}}{p_\mathrm{in}}}`
          without changing enthalpy.
        - in*->out* transfer heat without changing pressure.
          :math:`h_\mathrm{out*}-h_\mathrm{in*}=h_\mathrm{out}-h_\mathrm{in}`
        - out*->out decrease pressure to outlet pressure :math:`p_\mathrm{out}`
          without changing enthalpy.

        Note
        ----
        The entropy balance makes the follwing parameter available:

        .. math::

            \text{S\_Q1}=\dot{m} \cdot \left(s_\mathrm{out*,1}-s_\mathrm{in*,1}
            \right)\\
            \text{S\_Q2}=\dot{m} \cdot \left(s_\mathrm{out*,2}-s_\mathrm{in*,2}
            \right)\\
            \text{S\_Qirr}=\text{S\_Q2} - \text{S\_Q1}\\
            \text{S\_irr1}=\dot{m} \cdot \left(s_\mathrm{out,1}-s_\mathrm{in,1}
            \right) - \text{S\_Q1}\\
            \text{S\_irr2}=\dot{m} \cdot \left(s_\mathrm{out,2}-s_\mathrm{in,2}
            \right) - \text{S\_Q2}\\
            \text{S\_irr}=\sum \dot{S}_\mathrm{irr}\\
            \text{T\_mQ1}=\frac{\dot{Q}}{\text{S\_Q1}}\\
            \text{T\_mQ2}=\frac{\dot{Q}}{\text{S\_Q2}}
        """
        self.S_irr = 0
        for i in range(2):
            inl = self.inl[i]
            out = self.outl[i]
            p_star = inl.p.val_SI * (
                self.get_attr('pr' + str(i + 1)).val) ** 0.5
            s_i_star = s_mix_ph(
                [0, p_star, inl.h.val_SI, inl.fluid.val], T0=inl.T.val_SI)
            s_o_star = s_mix_ph(
                [0, p_star, out.h.val_SI, out.fluid.val], T0=out.T.val_SI)

            setattr(self, 'S_Q' + str(i + 1),
                    inl.m.val_SI * (s_o_star - s_i_star))
            S_Q = self.get_attr('S_Q' + str(i + 1))
            setattr(self, 'S_irr' + str(i + 1),
                    inl.m.val_SI * (out.s.val_SI - inl.s.val_SI) - S_Q)
            setattr(self, 'T_mQ' + str(i + 1),
                    inl.m.val_SI * (out.h.val_SI - inl.h.val_SI) / S_Q)

            self.S_irr += self.get_attr('S_irr' + str(i + 1))

        self.S_irr += self.S_Q1 + self.S_Q2

    def exergy_balance(self, T0):
        r"""
        Calculate exergy balance of a heat exchanger.

        Parameters
        ----------
        T0 : float
            Ambient temperature T0 / K.

        Note
        ----
        .. math::

            \dot{E}_\mathrm{P} =
            \begin{cases}
            \dot{E}_\mathrm{out,2}^\mathrm{T} -
            \dot{E}_\mathrm{in,2}^\mathrm{T}
            & T_\mathrm{in,1}, T_\mathrm{in,2}, T_\mathrm{out,1},
            T_\mathrm{out,2} > T_0\\
            \dot{E}_\mathrm{out,1}^\mathrm{T} -
            \dot{E}_\mathrm{in,1}^\mathrm{T}
            & T_0 \geq  T_\mathrm{in,1}, T_\mathrm{in,2}, T_\mathrm{out,1},
            T_\mathrm{out,2}\\
            \dot{E}_\mathrm{out,1}^\mathrm{T} +
            \dot{E}_\mathrm{out,2}^\mathrm{T}
            & T_\mathrm{in,1}, T_\mathrm{out,2} > T_0 \geq
            T_\mathrm{in,2}, T_\mathrm{out,1}\\
            \dot{E}_\mathrm{out,1}^\mathrm{T}
            & T_\mathrm{in,1} > T_0 \geq
            T_\mathrm{in,2}, T_\mathrm{out,1}, T_\mathrm{out,2}\\
            \text{not defined (nan)}
            & T_\mathrm{in,1}, T_\mathrm{out,1} > T_0 \geq
            T_\mathrm{in,2}, T_\mathrm{out,2}\\
            \dot{E}_\mathrm{out,2}^\mathrm{T}
            & T_\mathrm{in,1}, T_\mathrm{out,1},
            T_\mathrm{out,2} \geq T_0 > T_\mathrm{in,2}\\
            \end{cases}

            \dot{E}_\mathrm{F} =
            \begin{cases}
            \dot{E}_\mathrm{in,1}^\mathrm{PH} -
            \dot{E}_\mathrm{out,1}^\mathrm{PH} +
            \dot{E}_\mathrm{in,2}^\mathrm{M} -
            \dot{E}_\mathrm{out,2}^\mathrm{M}
            & T_\mathrm{in,1}, T_\mathrm{in,2}, T_\mathrm{out,1},
            T_\mathrm{out,2} > T_0\\
            \dot{E}_\mathrm{in,2}^\mathrm{PH} -
            \dot{E}_\mathrm{out,2}^\mathrm{PH} +
            \dot{E}_\mathrm{in,1}^\mathrm{M} -
            \dot{E}_\mathrm{out,1}^\mathrm{M}
            & T_0 \geq T_\mathrm{in,1}, T_\mathrm{in,2}, T_\mathrm{out,1},
            T_\mathrm{out,2}\\
            \dot{E}_\mathrm{in,1}^\mathrm{PH} +
            \dot{E}_\mathrm{in,2}^\mathrm{PH} -
            \dot{E}_\mathrm{out,1}^\mathrm{M} -
            \dot{E}_\mathrm{out,2}^\mathrm{M}
            & T_\mathrm{in,1}, T_\mathrm{out,2} > T_0 \geq
            T_\mathrm{in,2}, T_\mathrm{out,1}\\
            \dot{E}_\mathrm{in,1}^\mathrm{PH} +
            \dot{E}_\mathrm{in,2}^\mathrm{PH} -
            \dot{E}_\mathrm{out,2}^\mathrm{PH} -
            \dot{E}_\mathrm{out,1}^\mathrm{M}
            & T_\mathrm{in,1} > T_0 \geq
            T_\mathrm{in,2}, T_\mathrm{out,1}, T_\mathrm{out,2}\\
            \dot{E}_\mathrm{in,1}^\mathrm{PH} -
            \dot{E}_\mathrm{out,1}^\mathrm{PH} +
            \dot{E}_\mathrm{in,2}^\mathrm{PH} -
            \dot{E}_\mathrm{out,2}^\mathrm{PH}
            & T_\mathrm{in,1}, T_\mathrm{out,1} > T_0 \geq
            T_\mathrm{in,2}, T_\mathrm{out,2}\\
            \dot{E}_\mathrm{in,1}^\mathrm{PH} -
            \dot{E}_\mathrm{out,1}^\mathrm{PH} +
            \dot{E}_\mathrm{in,2}^\mathrm{PH} -
            \dot{E}_\mathrm{out,2}^\mathrm{M}
            & T_\mathrm{in,1}, T_\mathrm{out,1},
            T_\mathrm{out,2} \geq T_0 > T_\mathrm{in,2}\\
            \end{cases}
        """
        if all([c.T.val_SI > T0 for c in self.inl + self.outl]):
            self.E_P = self.outl[1].Ex_therm - self.inl[1].Ex_therm
            self.E_F = self.inl[0].Ex_physical - self.outl[0].Ex_physical + (
                self.inl[1].Ex_mech - self.outl[1].Ex_mech)
        elif all([c.T.val_SI <= T0 for c in self.inl + self.outl]):
            self.E_P = self.outl[0].Ex_therm - self.inl[0].Ex_therm
            self.E_F = self.inl[1].Ex_physical - self.outl[1].Ex_physical + (
                self.inl[0].Ex_mech - self.outl[0].Ex_mech)
        elif (self.inl[0].T.val_SI > T0 and self.outl[1].T.val_SI > T0 and
              self.outl[0].T.val_SI <= T0 and self.inl[1].T.val_SI <= T0):
            self.E_P = self.outl[0].Ex_therm + self.outl[1].Ex_therm
            self.E_F = self.inl[0].Ex_physical + self.inl[1].Ex_physical - (
                self.outl[0].Ex_mech + self.outl[1].Ex_mech)
        elif (self.inl[0].T.val_SI > T0 and self.inl[1].T.val_SI <= T0 and
              self.outl[0].T.val_SI <= T0 and self.outl[1].T.val_SI <= T0):
            self.E_P = self.outl[0].Ex_therm
            self.E_F = self.inl[0].Ex_physical + self.inl[1].Ex_physical - (
                self.outl[1].Ex_physical + self.outl[0].Ex_mech)
        elif (self.inl[0].T.val_SI > T0 and self.outl[0].T.val_SI > T0 and
              self.inl[1].T.val_SI <= T0 and self.outl[1].T.val_SI <= T0):
            self.E_P = np.nan
            self.E_F = self.inl[0].Ex_physical - self.outl[0].Ex_physical + (
                self.inl[1].Ex_physical - self.outl[1].Ex_physical)
        else:
            self.E_P = self.outl[1].Ex_therm
            self.E_F = self.inl[0].Ex_physical - self.outl[0].Ex_physical + (
                self.inl[1].Ex_physical - self.outl[1].Ex_mech)

        self.E_bus = np.nan
        if np.isnan(self.E_P):
            self.E_D = self.E_F
        else:
            self.E_D = self.E_F - self.E_P
        self.epsilon = self.E_P / self.E_F

    def get_plotting_data(self):
        """Generate a dictionary containing FluProDia plotting information.

        Returns
        -------
        data : dict
            A nested dictionary containing the keywords required by the
            :code:`calc_individual_isoline` method of the
            :code:`FluidPropertyDiagram` class. First level keys are the
            connection index ('in1' -> 'out1', therefore :code:`1` etc.).
        """
        return {
            i + 1: {
                'isoline_property': 'p',
                'isoline_value': self.inl[i].p.val,
                'isoline_value_end': self.outl[i].p.val,
                'starting_point_property': 'v',
                'starting_point_value': self.inl[i].vol.val,
                'ending_point_property': 'v',
                'ending_point_value': self.outl[i].vol.val
            } for i in range(2)}
