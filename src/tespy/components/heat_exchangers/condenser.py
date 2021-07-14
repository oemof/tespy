# -*- coding: utf-8

"""Module of class Condenser.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/heat_exchangers/condenser.py

SPDX-License-Identifier: MIT
"""

import numpy as np

from tespy.components.component import Component
from tespy.components.heat_exchangers.heat_exchanger import HeatExchanger
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import DataContainerSimple as dc_simple
from tespy.tools.data_containers import GroupedComponentCharacteristics as dc_gcc
from tespy.tools.document_models import generate_latex_eq
from tespy.tools.fluid_properties import T_bp_p
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.fluid_properties import dh_mix_dpQ
from tespy.tools.fluid_properties import h_mix_pQ


class Condenser(HeatExchanger):
    r"""
    A Condenser cools a fluid until it is in liquid state.

    The condensing fluid is cooled by the cold side fluid. The fluid on the hot
    side of the condenser must be pure. Subcooling is available.

    **Mandatory Equations**

    - :py:meth:`tespy.components.component.Component.fluid_func`
    - :py:meth:`tespy.components.component.Component.mass_flow_func`
    - :py:meth:`tespy.components.heat_exchangers.heat_exchanger.HeatExchanger.energy_balance_func`
    - condensate outlet state, function can be disabled by specifying
      :code:`set_attr(subcooling=True)`
      :py:meth:`tespy.components.heat_exchangers.condenser.Condenser.subcooling_func`

    **Optional Equations**

    - :py:meth:`tespy.components.heat_exchangers.heat_exchanger.HeatExchanger.energy_balance_hot_func`
    - :py:meth:`tespy.components.heat_exchangers.condenser.Condenser.kA_func`
    - :py:meth:`tespy.components.heat_exchangers.condenser.Condenser.kA_char_func`
    - :py:meth:`tespy.components.heat_exchangers.condenser.Condenser.ttd_u_func`
    - :py:meth:`tespy.components.heat_exchangers.heat_exchanger.HeatExchanger.ttd_l_func`
    - hot side :py:meth:`tespy.components.component.Component.pr_func`
    - cold side :py:meth:`tespy.components.component.Component.pr_func`
    - hot side :py:meth:`tespy.components.component.Component.zeta_func`
    - cold side :py:meth:`tespy.components.component.Component.zeta_func`

    Inlets/Outlets

    - in1, in2 (index 1: hot side, index 2: cold side)
    - out1, out2 (index 1: hot side, index 2: cold side)

    Image

    .. image:: _images/Condenser.svg
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
        Upper terminal temperature difference (referring to saturation
        temprature of condensing fluid) :math:`ttd_\mathrm{u}/\text{K}`.

    kA : float, dict
        Area independent heat transfer coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char : tespy.tools.data_containers.DataContainerSimple
        Area independent heat transfer coefficient characteristic.

    kA_char1 : tespy.tools.characteristics.CharLine, dict
        Characteristic line for hot side heat transfer coefficient.

    kA_char2 : tespy.tools.characteristics.CharLine, dict
        Characteristic line for cold side heat transfer coefficient.

    subcooling : boolean
        Enable/disable subcooling, default value: disabled.

    Note
    ----
    The condenser has an additional equation for enthalpy at hot side outlet:
    The fluid leaves the component in saturated liquid state. If subcooling
    is activated, it possible to specify the enthalpy at the outgoing
    connection manually.

    It has different calculation method for given heat transfer coefficient and
    upper terminal temperature dierence: These parameters refer to the
    **condensing** temperature, even if the fluid on the hot side enters the
    component in superheated state.

    Example
    -------
    Air is used to condensate water in a condenser. 1 kg/s waste steam is
    chilled with a terminal temperature difference of 15 K.

    >>> from tespy.components import Sink, Source, Condenser
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> from tespy.tools.fluid_properties import T_bp_p
    >>> import shutil
    >>> nw = Network(fluids=['water', 'air'], T_unit='C', p_unit='bar',
    ... h_unit='kJ / kg', m_range=[0.01, 1000], iterinfo=False)
    >>> amb_in = Sink('ambient air inlet')
    >>> amb_out = Source('air outlet')
    >>> waste_steam = Source('waste steam')
    >>> c = Sink('condensate sink')
    >>> cond = Condenser('condenser')
    >>> cond.component()
    'condenser'
    >>> amb_he = Connection(amb_out, 'out1', cond, 'in2')
    >>> he_amb = Connection(cond, 'out2', amb_in, 'in1')
    >>> ws_he = Connection(waste_steam, 'out1', cond, 'in1')
    >>> he_c = Connection(cond, 'out1', c, 'in1')
    >>> nw.add_conns(amb_he, he_amb, ws_he, he_c)

    The air flow can not be controlled, thus is constant in offdesign
    operation. If the waste steam mass flow or the ambient air temperature
    change, the outlet temperature of the air will change, too.

    >>> cond.set_attr(pr1=0.98, pr2=0.999, ttd_u=15, design=['pr2', 'ttd_u'],
    ... offdesign=['zeta2', 'kA_char'])
    >>> ws_he.set_attr(fluid={'water': 1, 'air': 0}, h=2700, m=1)
    >>> amb_he.set_attr(fluid={'water': 0, 'air': 1}, T=20, offdesign=['v'])
    >>> he_amb.set_attr(p=1, T=40, design=['T'])
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(amb_he.v.val, 2)
    103.17
    >>> round(ws_he.T.val - he_amb.T.val, 1)
    66.9
    >>> round(T_bp_p(ws_he.get_flow()) - 273.15 - he_amb.T.val, 1)
    15.0
    >>> ws_he.set_attr(m=0.7)
    >>> amb_he.set_attr(T=30)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(ws_he.T.val - he_amb.T.val, 1)
    62.5
    >>> round(T_bp_p(ws_he.get_flow()) - 273.15 - he_amb.T.val, 1)
    11.3

    It is possible to activate subcooling. The difference to boiling point
    temperature is specified to 5 K.

    >>> cond.set_attr(subcooling=True)
    >>> he_c.set_attr(Td_bp=-5)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(ws_he.T.val - he_amb.T.val, 1)
    62.5
    >>> round(T_bp_p(ws_he.get_flow()) - 273.15 - he_amb.T.val, 1)
    13.4
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    @staticmethod
    def component():
        return 'condenser'

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
                    'type': 'rel', 'inconn': 1, 'outconn': 1}),
            'subcooling': dc_simple(
                val=False, num_eq=1, latex=self.subcooling_func_doc,
                deriv=self.subcooling_deriv, func=self.subcooling_func)
        }

    def comp_init(self, nw):

        # if subcooling is True, outlet state method must not be calculated
        self.subcooling.is_set = not self.subcooling.val
        Component.comp_init(self, nw)

    def subcooling_func(self):
        r"""
        Equation for hot side outlet state.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0=h_{out,1} -h\left(p_{out,1}, x=0 \right)

        Note
        ----
        This equation is applied in case subcooling is False!
        """
        return self.outl[0].h.val_SI - h_mix_pQ(self.outl[0].get_flow(), 0)

    def subcooling_func_doc(self, label):
        r"""
        Equation for hot side outlet state.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        latex = r'0=h_\mathrm{out,1} -h\left(p_\mathrm{out,1}, x=0 \right)'
        return generate_latex_eq(self, latex, label)

    def subcooling_deriv(self, increment_filter, k):
        """
        Calculate partial derivates of subcooling function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        self.jacobian[k, 2, 1] = -dh_mix_dpQ(self.outl[0].get_flow(), 0)
        self.jacobian[k, 2, 2] = 1

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
                T_{in,2} - T_{sat} \left(p_{in,1}\right) + T_{out,2}}
                {\ln{\frac{T_{out,1} - T_{in,2}}
                {T_{sat} \left(p_{in,1}\right) - T_{out,2}}}}
        """
        i1 = self.inl[0]
        i2 = self.inl[1]
        o1 = self.outl[0]
        o2 = self.outl[1]

        T_i1 = T_bp_p(i1.get_flow())
        T_i2 = T_mix_ph(i2.get_flow(), T0=i2.T.val_SI)
        T_o1 = T_mix_ph(o1.get_flow(), T0=o1.T.val_SI)
        T_o2 = T_mix_ph(o2.get_flow(), T0=o2.T.val_SI)

        if T_i1 <= T_o2 and not i1.T.val_set:
            T_i1 = T_o2 + 0.5
        if T_i1 <= T_o2 and not o2.T.val_set:
            T_o2 = T_i1 - 0.5
        if T_o1 <= T_i2 and not o1.T.val_set:
            T_o1 = T_i2 + 1
        if T_o1 <= T_i2 and not i2.T.val_set:
            T_i2 = T_o1 - 1

        td_log = ((T_o1 - T_i2 - T_i1 + T_o2) /
                  np.log((T_o1 - T_i2) / (T_i1 - T_o2)))

        return i1.m.val_SI * (o1.h.val_SI - i1.h.val_SI) + self.kA.val * td_log

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
            r'T_\mathrm{in,2} -T_\mathrm{sat}\left( p_\mathrm{in,1}\right)'
            r'+ T_\mathrm{out,2}}'
            r'{\ln{\frac{T_\mathrm{out,1} - T_\mathrm{in,2}}'
            r'{T_\mathrm{sat}\left( p_\mathrm{in,1}\right) -'
            r'T_\mathrm{out,2}}}}'
        )
        return generate_latex_eq(self, latex, label)

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
                T_{in,2} - T_{sat} \left(p_{in,1}\right) + T_{out,2}}
                {\ln{\frac{T_{out,1} - T_{in,2}}
                {T_{sat} \left(p_{in,1}\right) - T_{out,2}}}}

                f_{kA} = \frac{2}{\frac{1}{f_1 \left( expr_1\right)} +
                \frac{1}{f_2 \left( expr_2\right)}}

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
        T_i1 = T_bp_p(i1.get_flow())
        T_i2 = T_mix_ph(i2.get_flow(), T0=i2.T.val_SI)
        T_o1 = T_mix_ph(o1.get_flow(), T0=o1.T.val_SI)
        T_o2 = T_mix_ph(o2.get_flow(), T0=o2.T.val_SI)

        if T_i1 <= T_o2 and not i1.T.val_set:
            T_i1 = T_o2 + 0.5
        if T_i1 <= T_o2 and not o2.T.val_set:
            T_o2 = T_i1 - 0.5
        if T_o1 <= T_i2 and not o1.T.val_set:
            T_o1 = T_i2 + 1
        if T_o1 <= T_i2 and not i2.T.val_set:
            T_i2 = T_o1 - 1

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
            r' - T_\mathrm{sat}\left( p_\mathrm{in,1}\right) +'
            r'T_\mathrm{out,2}}{\ln{\frac{T_\mathrm{out,1}-'
            r'T_\mathrm{in,2}}{T_\mathrm{sat}\left( p_\mathrm{in,1}\right)'
            r'- T_\mathrm{out,2}}}}\\' + '\n'
            r'f_\mathrm{kA}=&\frac{2}{\frac{1}{f\left(X_1\right)}+'
            r'\frac{1}{f\left(X_2\right)}}\\' + '\n'
            r'\end{split}'
        )
        return generate_latex_eq(self, latex, label)

    def ttd_u_func(self):
        r"""
        Equation for upper terminal temperature difference.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = ttd_{u} - T_{sat} \left(p_{in,1}\right) + T_{out,2}

        Note
        ----
        The upper terminal temperature difference ttd_u refers to boiling
        temperature at hot side inlet.
        """
        T_i1 = T_bp_p(self.inl[0].get_flow())
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
        latex = (
            r'0=ttd_\mathrm{u}-T_\mathrm{sat}\left(p_\mathrm{in,1}\right)'
            r' + T_\mathrm{out,2}')
        return generate_latex_eq(self, latex, label)

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        # component parameters
        self.Q.val = self.inl[0].m.val_SI * (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI)
        self.ttd_u.val = T_bp_p(self.inl[0].get_flow()) - self.outl[1].T.val_SI
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
