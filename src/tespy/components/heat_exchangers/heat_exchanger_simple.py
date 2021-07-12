# -*- coding: utf-8

"""Module of class HeatExchangerSimple.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/heat_exchangers/heat_exchanger_simple.py

SPDX-License-Identifier: MIT
"""

import logging

import numpy as np

from tespy.components.component import Component
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import DataContainerSimple as dc_simple
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.document_models import generate_latex_eq
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.fluid_properties import s_mix_ph
from tespy.tools.fluid_properties import v_mix_ph
from tespy.tools.fluid_properties import visc_mix_ph
from tespy.tools.helpers import convert_to_SI
from tespy.tools.helpers import darcy_friction_factor as dff


class HeatExchangerSimple(Component):
    r"""
    A basic heat exchanger representing a heat source or heat sink.

    The component HeatExchangerSimple is the parent class for the components:

    - :py:class:`tespy.components.heat_exchangers.solar_collector.SolarCollector`
    - :py:class:`tespy.components.heat_exchangers.parabolic_trough.ParabolicTrough`
    - :py:class:`tespy.components.piping.pipe.Pipe`

    **Mandatory Equations**

    - :py:meth:`tespy.components.component.Component.fluid_func`
    - :py:meth:`tespy.components.component.Component.mass_flow_func`

    **Optional Equations**

    - :py:meth:`tespy.components.component.Component.pr_func`
    - :py:meth:`tespy.components.component.Component.zeta_func`
    - :py:meth:`tespy.components.heat_exchangers.heat_exchanger_simple.HeatExchangerSimple.energy_balance_func`
    - :py:meth:`tespy.components.heat_exchangers.heat_exchanger_simple.HeatExchangerSimple.hydro_group_func`
    - :py:meth:`tespy.components.heat_exchangers.heat_exchanger_simple.HeatExchangerSimple.kA_group_func`
    - :py:meth:`tespy.components.heat_exchangers.heat_exchanger_simple.HeatExchangerSimple.kA_char_group_func`

    Inlets/Outlets

    - in1
    - out1

    Image

    .. image:: _images/Pipe.svg
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

    Q : float, dict, :code:`"var"`
        Heat transfer, :math:`Q/\text{W}`.

    pr : float, dict, :code:`"var"`
        Outlet to inlet pressure ratio, :math:`pr/1`.

    zeta : float, dict, :code:`"var"`
        Geometry independent friction coefficient,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    D : float, dict, :code:`"var"`
        Diameter of the pipes, :math:`D/\text{m}`.

    L : float, dict, :code:`"var"`
        Length of the pipes, :math:`L/\text{m}`.

    ks : float, dict, :code:`"var"`
        Pipe's roughness, :math:`ks/\text{m}` for darcy friction,
        :math:`ks/\text{1}` for hazen-williams equation.

    hydro_group : str, dict
        Parametergroup for pressure drop calculation based on pipes dimensions.
        Choose 'HW' for hazen-williams equation, else darcy friction factor is
        used.

    kA : float, dict, :code:`"var"`
        Area independent heat transfer coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char : tespy.tools.characteristics.CharLine, dict
        Characteristic line for heat transfer coefficient.

    Tamb : float, dict
        Ambient temperature, provide parameter in network's temperature unit.

    kA_group : str, dict
        Parametergroup for heat transfer calculation from ambient temperature
        and area independent heat transfer coefficient kA.

    Example
    -------
    The HeatExchangerSimple can be used as a sink or source of heat. This
    component does not simulate the secondary side of the heat exchanger. It
    is possible to calculate the pressure ratio with the Darcy-Weisbach
    equation or in case of liquid water use the Hazen-Williams equation.
    Also, given ambient temperature and the heat transfer coeffiecient, it is
    possible to predict heat transfer.

    >>> from tespy.components import Sink, Source, HeatExchangerSimple
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> import shutil
    >>> fluids = ['N2']
    >>> nw = Network(fluids=fluids)
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg', iterinfo=False)
    >>> so1 = Source('source 1')
    >>> si1 = Sink('sink 1')
    >>> heat_sink = HeatExchangerSimple('heat sink')
    >>> heat_sink.component()
    'heat exchanger simple'
    >>> heat_sink.set_attr(Tamb=10, pr=0.95, design=['pr'],
    ... offdesign=['zeta', 'kA_char'])
    >>> inc = Connection(so1, 'out1', heat_sink, 'in1')
    >>> outg = Connection(heat_sink, 'out1', si1, 'in1')
    >>> nw.add_conns(inc, outg)

    It is possible to determine the amount of heat transferred when the fluid
    enters the heat sink at a temperature of 200 °C and is cooled down to
    150 °C. Given an ambient temperature of 10 °C this also determines the heat
    transfer coefficient to the ambient. Assuming a characteristic function
    for the heat transfer coefficient we can predict the heat transferred at
    variable flow rates.

    >>> inc.set_attr(fluid={'N2': 1}, m=1, T=200, p=5)
    >>> outg.set_attr(T=150, design=['T'])
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(heat_sink.Q.val, 0)
    -52581.0
    >>> round(heat_sink.kA.val, 0)
    321.0
    >>> inc.set_attr(m=1.25)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(heat_sink.Q.val, 0)
    -56599.0
    >>> round(outg.T.val, 1)
    156.9
    >>> inc.set_attr(m=0.75)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(heat_sink.Q.val, 1)
    -47275.8
    >>> round(outg.T.val, 1)
    140.0
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    @staticmethod
    def component():
        return 'heat exchanger simple'

    def get_variables(self):
        return {
            'Q': dc_cp(
                deriv=self.energy_balance_deriv,
                latex=self.energy_balance_func_doc, num_eq=1,
                func=self.energy_balance_func),
            'pr': dc_cp(
                min_val=1e-4, max_val=1, num_eq=1,
                deriv=self.pr_deriv, latex=self.pr_func_doc,
                func=self.pr_func, func_params={'pr': 'pr'}),
            'zeta': dc_cp(
                min_val=0, max_val=1e15, num_eq=1,
                deriv=self.zeta_deriv, func=self.zeta_func,
                latex=self.zeta_func_doc,
                func_params={'zeta': 'zeta'}),
            'D': dc_cp(min_val=1e-2, max_val=2, d=1e-4),
            'L': dc_cp(min_val=1e-1, d=1e-3),
            'ks': dc_cp(val=1e-4, min_val=1e-7, max_val=1e-3, d=1e-8),
            'kA': dc_cp(min_val=0, d=1),
            'kA_char': dc_cc(param='m'), 'Tamb': dc_cp(),
            'dissipative': dc_simple(val=True),
            'hydro_group': dc_gcp(
                elements=['L', 'ks', 'D'], num_eq=1,
                latex=self.hydro_group_func_doc,
                func=self.hydro_group_func, deriv=self.hydro_group_deriv),
            'kA_group': dc_gcp(
                elements=['kA', 'Tamb'], num_eq=1,
                latex=self.kA_group_func_doc,
                func=self.kA_group_func, deriv=self.kA_group_deriv),
            'kA_char_group': dc_gcp(
                elements=['kA_char', 'Tamb'], num_eq=1,
                latex=self.kA_char_group_func_doc,
                func=self.kA_char_group_func, deriv=self.kA_char_group_deriv)
        }

    @staticmethod
    def inlets():
        return ['in1']

    @staticmethod
    def outlets():
        return ['out1']

    def comp_init(self, nw):
        Component.comp_init(self, nw, num_eq=len(nw.fluids) + 1)

        self.Tamb.val_SI = convert_to_SI('T', self.Tamb.val, nw.T_unit)

    def energy_balance_func(self):
        r"""
        Equation for pressure drop calculation.

        Returns
        -------
        residual : float
            Residual value of equation:

            .. math::

                0 =\dot{m}_{in}\cdot\left( h_{out}-h_{in}\right) -\dot{Q}
        """
        return self.inl[0].m.val_SI * (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI) - self.Q.val

    def energy_balance_func_doc(self, label):
        r"""
        Equation for pressure drop calculation.

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
            r'0 = \dot{m}_\mathrm{in} \cdot \left(h_\mathrm{out} - '
            r'h_\mathrm{in} \right) -\dot{Q}'
        )
        return generate_latex_eq(self, latex, label)

    def energy_balance_deriv(self, increment_filter, k):
        r"""
        Calculate partial derivatives of energy balance.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        self.jacobian[k, 0, 0] = (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI)
        self.jacobian[k, 0, 2] = -self.inl[0].m.val_SI
        self.jacobian[k, 1, 2] = self.inl[0].m.val_SI
        # custom variable Q
        if self.Q.is_var:
            self.jacobian[k, 2 + self.Q.var_pos, 0] = -1

    def hydro_group_func(self):
        r"""
        Equation for pressure drop calculation.

        Returns
        -------
        residual : float
            Residual value of corresponding equation:

            - :py:meth:`tespy.components.heat_exchangers.heat_exchanger_simple.HeatExchangerSimple.darcy_func`
            - :py:meth:`tespy.components.heat_exchangers.heat_exchanger_simple.HeatExchangerSimple.hazen_williams_func`
        """
        # hazen williams equation
        if self.hydro_group.method == 'HW':
            return self.hazen_williams_func()
        # darcy friction factor
        else:
            return self.darcy_func()

    def hydro_group_func_doc(self, label):
        r"""
        Equation for pressure drop calculation.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        # hazen williams equation
        if self.hydro_group.method == 'HW':
            return self.hazen_williams_func_doc(label)
        # darcy friction factor
        else:
            return self.darcy_func_doc(label)

    def hydro_group_deriv(self, increment_filter, k):
        r"""
        Calculate partial derivatives of hydro group (pressure drop).

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        # hazen williams equation
        if self.hydro_group.method == 'HW':
            func = self.hazen_williams_func
        # darcy friction factor
        else:
            func = self.darcy_func
        if not increment_filter[0, 0]:
            self.jacobian[k, 0, 0] = self.numeric_deriv(func, 'm', 0)
        if not increment_filter[0, 1]:
            self.jacobian[k, 0, 1] = self.numeric_deriv(func, 'p', 0)
        if not increment_filter[0, 2]:
            self.jacobian[k, 0, 2] = self.numeric_deriv(func, 'h', 0)
        if not increment_filter[1, 1]:
            self.jacobian[k, 1, 1] = self.numeric_deriv(func, 'p', 1)
        if not increment_filter[1, 2]:
            self.jacobian[k, 1, 2] = self.numeric_deriv(func, 'h', 1)
        # custom variables of hydro group
        for var in self.hydro_group.elements:
            var = self.get_attr(var)
            if var.is_var:
                self.jacobian[k, 2 + var.var_pos, 0] = (
                    self.numeric_deriv(func, self.vars[var], 2))

    def darcy_func(self):
        r"""
        Equation for pressure drop calculation from darcy friction factor.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = p_{in} - p_{out} - \frac{8 \cdot |\dot{m}_{in}| \cdot
                \dot{m}_{in} \cdot \frac{v_{in}+v_{out}}{2} \cdot L \cdot
                \lambda\left(Re, ks, D\right)}{\pi^2 \cdot D^5}\\

                Re = \frac{4 \cdot |\dot{m}_{in}|}{\pi \cdot D \cdot
                \frac{\eta_{in}+\eta_{out}}{2}}\\
                \eta: \text{dynamic viscosity}\\
                v: \text{specific volume}\\
                \lambda: \text{darcy friction factor}
        """
        i, o = self.inl[0].get_flow(), self.outl[0].get_flow()

        if abs(i[0]) < 1e-4:
            return i[1] - o[1]

        visc_i = visc_mix_ph(i, T0=self.inl[0].T.val_SI)
        visc_o = visc_mix_ph(o, T0=self.outl[0].T.val_SI)
        v_i = v_mix_ph(i, T0=self.inl[0].T.val_SI)
        v_o = v_mix_ph(o, T0=self.outl[0].T.val_SI)

        Re = 4 * abs(i[0]) / (np.pi * self.D.val * (visc_i + visc_o) / 2)

        return ((i[1] - o[1]) - 8 * abs(i[0]) * i[0] * (v_i + v_o) / 2 *
                self.L.val * dff(Re, self.ks.val, self.D.val) /
                (np.pi ** 2 * self.D.val ** 5))

    def darcy_func_doc(self, label):
        r"""
        Equation for pressure drop calculation from darcy friction factor.

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
            r'0 = &p_\mathrm{in}-p_\mathrm{out}-'
            r'\frac{8\cdot|\dot{m}_\mathrm{in}| \cdot\dot{m}_\mathrm{in}'
            r'\cdot \frac{v_\mathrm{in}+v_\mathrm{out}}{2} \cdot L \cdot'
            r'\lambda\left(Re, ks, D\right)}{\pi^2 \cdot D^5}\\' + '\n'
            r'Re =&\frac{4 \cdot |\dot{m}_\mathrm{in}|}{\pi \cdot D \cdot'
            r'\frac{\eta_\mathrm{in}+\eta_\mathrm{out}}{2}}\\' + '\n'
            r'\end{split}'
        )
        return generate_latex_eq(self, latex, label)

    def hazen_williams_func(self):
        r"""
        Equation for pressure drop calculation from Hazen-Williams equation.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \left(p_{in} - p_{out} \right) \cdot \left(-1\right)^i -
                \frac{10.67 \cdot |\dot{m}_{in}| ^ {1.852}
                \cdot L}{ks^{1.852} \cdot D^{4.871}} \cdot g \cdot
                \left(\frac{v_{in} + v_{out}}{2}\right)^{0.852}

                i = \begin{cases}
                0 & \dot{m}_{in} \geq 0\\
                1 & \dot{m}_{in} < 0
                \end{cases}

        Note
        ----
        Gravity :math:`g` is set to :math:`9.81 \frac{m}{s^2}`
        """
        i, o = self.inl[0].get_flow(), self.outl[0].get_flow()

        if abs(i[0]) < 1e-4:
            return i[1] - o[1]

        v_i = v_mix_ph(i, T0=self.inl[0].T.val_SI)
        v_o = v_mix_ph(o, T0=self.outl[0].T.val_SI)

        return ((i[1] - o[1]) * np.sign(i[0]) -
                (10.67 * abs(i[0]) ** 1.852 * self.L.val /
                 (self.ks.val ** 1.852 * self.D.val ** 4.871)) *
                (9.81 * ((v_i + v_o) / 2) ** 0.852))

    def hazen_williams_func_doc(self, label):
        r"""
        Equation for pressure drop calculation from Hazen-Williams equation.

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
            r'0 = \left(p_\mathrm{in} - p_\mathrm{out} \right) -'
            r'\frac{10.67 \cdot |\dot{m}_\mathrm{in}| ^ {1.852}'
            r'\cdot L}{ks^{1.852} \cdot D^{4.871}} \cdot g \cdot'
            r'\left(\frac{v_\mathrm{in}+ v_\mathrm{out}}{2}\right)^{0.852}'
        )
        return generate_latex_eq(self, latex, label)

    def kA_group_func(self):
        r"""
        Calculate heat transfer from heat transfer coefficient.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_{in} \cdot \left( h_{out} - h_{in}\right) +
                kA \cdot \Delta T_{log}

                \Delta T_{log} = \begin{cases}
                \frac{T_{in}-T_{out}}{\ln{\frac{T_{in}-T_{amb}}
                {T_{out}-T_{amb}}}} & T_{in} > T_{out} \\
                \frac{T_{out}-T_{in}}{\ln{\frac{T_{out}-T_{amb}}
                {T_{in}-T_{amb}}}} & T_{in} < T_{out}\\
                0 & T_{in} = T_{out}
                \end{cases}

                T_{amb}: \text{ambient temperature}
        """
        i, o = self.inl[0].get_flow(), self.outl[0].get_flow()

        ttd_1 = T_mix_ph(i, T0=self.inl[0].T.val_SI) - self.Tamb.val_SI
        ttd_2 = T_mix_ph(o, T0=self.outl[0].T.val_SI) - self.Tamb.val_SI

        # For numerical stability: If temperature differences have
        # different sign use mean difference to avoid negative logarithm.
        if (ttd_1 / ttd_2) < 0:
            td_log = (ttd_2 + ttd_1) / 2
        elif ttd_1 > ttd_2:
            td_log = (ttd_1 - ttd_2) / np.log(ttd_1 / ttd_2)
        elif ttd_1 < ttd_2:
            td_log = (ttd_2 - ttd_1) / np.log(ttd_2 / ttd_1)
        else:
            td_log = 0

        return i[0] * (o[2] - i[2]) + self.kA.val * td_log

    def kA_group_func_doc(self, label):
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
            r'\begin{split}' + '\n'
            r'0=&\dot{m}_\mathrm{in}\cdot\left(h_\mathrm{out}-'
            r'h_\mathrm{in}\right)+kA \cdot \Delta T_\mathrm{log}\\' + '\n'
            r'\Delta T_\mathrm{log} = &\begin{cases}' + '\n'
            r'\frac{T_\mathrm{in}-T_\mathrm{out}}{\ln{\frac{T_\mathrm{in}-'
            r'T_\mathrm{amb}}{T_\mathrm{out}-T_\mathrm{amb}}}} &'
            r' T_\mathrm{in} > T_\mathrm{out} \\' + '\n'
            r'\frac{T_\mathrm{out}-T_\mathrm{in}}{\ln{\frac{'
            r'T_\mathrm{out}-T_\mathrm{amb}}{T_\mathrm{in}-'
            r'T_\mathrm{amb}}}} & T_\mathrm{in} < T_\mathrm{out}\\' + '\n'
            r'0 & T_\mathrm{in} = T_\mathrm{out}' + '\n'
            r'\end{cases}\\' + '\n'
            r'T_\mathrm{amb} =& \text{ambient temperature}' + '\n'
            r'\end{split}'
        )
        return generate_latex_eq(self, latex, label)

    def kA_group_deriv(self, increment_filter, k):
        r"""
        Calculate partial derivatives of kA group.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        f = self.kA_group_func
        self.jacobian[k, 0, 0] = (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI)
        if not increment_filter[0, 1]:
            self.jacobian[k, 0, 1] = self.numeric_deriv(f, 'p', 0)
        if not increment_filter[0, 2]:
            self.jacobian[k, 0, 2] = self.numeric_deriv(f, 'h', 0)
        if not increment_filter[1, 1]:
            self.jacobian[k, 1, 1] = self.numeric_deriv(f, 'p', 1)
        if not increment_filter[1, 2]:
            self.jacobian[k, 1, 2] = self.numeric_deriv(f, 'h', 1)
        if self.kA.is_var:
            self.jacobian[k, 2 + self.kA.var_pos, 0] = (
                self.numeric_deriv(f, self.vars[self.kA], 2))

    def kA_char_group_func(self):
        r"""
        Calculate heat transfer from heat transfer coefficient characteristic.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_{in} \cdot \left( h_{out} - h_{in}\right) +
                kA_{design} \cdot f_{kA} \cdot \Delta T_{log}

                \Delta T_{log} = \begin{cases}
                \frac{T_{in}-T_{out}}{\ln{\frac{T_{in}-T_{amb}}
                {T_{out}-T_{amb}}}} & T_{in} > T_{out} \\
                \frac{T_{out}-T_{in}}{\ln{\frac{T_{out}-T_{amb}}
                {T_{in}-T_{amb}}}} & T_{in} < T_{out}\\
                0 & T_{in} = T_{out}
                \end{cases}

                f_{kA} = \frac{2}{1 + \frac{1}{f\left( expr\right)}}

                T_{amb}: \text{ambient temperature}

        Note
        ----
        For standard function of f\ :subscript:`kA` \ see module
        :py:mod:`tespy.data`.
        """
        p = self.kA_char.param
        expr = self.get_char_expr(p, **self.kA_char.char_params)
        i, o = self.inl[0].get_flow(), self.outl[0].get_flow()

        # For numerical stability: If temperature differences have
        # different sign use mean difference to avoid negative logarithm.

        ttd_1 = T_mix_ph(i, T0=self.inl[0].T.val_SI) - self.Tamb.val_SI
        ttd_2 = T_mix_ph(o, T0=self.outl[0].T.val_SI) - self.Tamb.val_SI

        if (ttd_1 / ttd_2) < 0:
            td_log = (ttd_2 + ttd_1) / 2
        elif ttd_1 > ttd_2:
            td_log = (ttd_1 - ttd_2) / np.log(ttd_1 / ttd_2)
        elif ttd_1 < ttd_2:
            td_log = (ttd_2 - ttd_1) / np.log(ttd_2 / ttd_1)
        else:
            td_log = 0

        fkA = 2 / (1 + 1 / self.kA_char.char_func.evaluate(expr))

        return i[0] * (o[2] - i[2]) + self.kA.design * fkA * td_log

    def kA_char_group_func_doc(self, label):
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
            r'0=&\dot{m}_\mathrm{in}\cdot\left(h_\mathrm{out}-'
            r'h_\mathrm{in}\right)+kA_\mathrm{design} \cdot f_\mathrm{kA}'
            r' \cdot \Delta T_\mathrm{log}\\' + '\n'
            r'\Delta T_\mathrm{log} = &\begin{cases}' + '\n'
            r'\frac{T_\mathrm{in}-T_\mathrm{out}}{\ln{\frac{T_\mathrm{in}-'
            r'T_\mathrm{amb}}{T_\mathrm{out}-T_\mathrm{amb}}}} &'
            r' T_\mathrm{in} > T_\mathrm{out} \\' + '\n'
            r'\frac{T_\mathrm{out}-T_\mathrm{in}}{\ln{\frac{'
            r'T_\mathrm{out}-T_\mathrm{amb}}{T_\mathrm{in}-'
            r'T_\mathrm{amb}}}} & T_\mathrm{in} < T_\mathrm{out}\\' + '\n'
            r'0 & T_\mathrm{in} = T_\mathrm{out}' + '\n'
            r'\end{cases}\\' + '\n'
            r'f_{kA}=&\frac{2}{1 + \frac{1}{f\left(X\right)}}\\' + '\n'
            r'T_\mathrm{amb} =& \text{ambient temperature}' + '\n'
            r'\end{split}'
        )
        return generate_latex_eq(self, latex, label)

    def kA_char_group_deriv(self, increment_filter, k):
        r"""
        Calculate partial derivatives of kA characteristics.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        f = self.kA_char_group_func
        if not increment_filter[0, 0]:
            self.jacobian[k, 0, 0] = self.numeric_deriv(f, 'm', 0)
        if not increment_filter[0, 1]:
            self.jacobian[k, 0, 1] = self.numeric_deriv(f, 'p', 0)
        if not increment_filter[0, 2]:
            self.jacobian[k, 0, 2] = self.numeric_deriv(f, 'h', 0)
        if not increment_filter[1, 1]:
            self.jacobian[k, 1, 1] = self.numeric_deriv(f, 'p', 1)
        if not increment_filter[1, 2]:
            self.jacobian[k, 1, 2] = self.numeric_deriv(f, 'h', 1)

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

                \dot{E} = \dot{m}_{in} \cdot \left( h_{out} - h_{in} \right)
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
            r'\dot{m}_\mathrm{in} \cdot \left(h_\mathrm{out} - '
            r'h_\mathrm{in} \right)')

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
        deriv = np.zeros((1, 2, self.num_nw_vars))
        f = self.calc_bus_value
        deriv[0, 0, 0] = self.numeric_deriv(f, 'm', 0, bus=bus)
        deriv[0, 0, 2] = self.numeric_deriv(f, 'h', 0, bus=bus)
        deriv[0, 1, 2] = self.numeric_deriv(f, 'h', 1, bus=bus)
        return deriv

    def initialise_source(self, c, key):
        r"""
        Return a starting value for pressure and enthalpy the outlets.

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
                \begin{cases}
                1 \cdot 10^5 \; \frac{\text{J}}{\text{kg}} & \dot{Q} < 0\\
                3 \cdot 10^5 \; \frac{\text{J}}{\text{kg}} & \dot{Q} = 0\\
                5 \cdot 10^5 \; \frac{\text{J}}{\text{kg}} & \dot{Q} > 0
                \end{cases} & \text{key = 'h'}\\
                \; \; \; \; 10^5 \text{Pa} & \text{key = 'p'}
                \end{cases}

        """
        if key == 'p':
            return 1e5
        elif key == 'h':
            if self.Q.val < 0 and self.Q.is_set:
                return 1e5
            elif self.Q.val > 0 and self.Q.is_set:
                return 5e5
            else:
                return 3e5

    def initialise_target(self, c, key):
        r"""
        Return a starting value for pressure and enthalpy the inlets.

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
                1 \cdot 10^5 & \text{key = 'p'}\\
                \begin{cases}
                5 \cdot 10^5 & \dot{Q} < 0\\
                3 \cdot 10^5 & \dot{Q} = 0\\
                1 \cdot 10^5 & \dot{Q} > 0
                \end{cases} & \text{key = 'h'}\\
                \end{cases}
        """
        if key == 'p':
            return 1e5
        elif key == 'h':
            if self.Q.val < 0 and self.Q.is_set:
                return 5e5
            elif self.Q.val > 0 and self.Q.is_set:
                return 1e5
            else:
                return 3e5

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        i = self.inl[0].get_flow()
        o = self.outl[0].get_flow()

        self.Q.val = i[0] * (o[2] - i[2])
        self.pr.val = o[1] / i[1]
        self.zeta.val = ((i[1] - o[1]) * np.pi ** 2 / (
            4 * i[0] ** 2 * (self.inl[0].vol.val_SI + self.outl[0].vol.val_SI)
            ))

        if self.Tamb.is_set:
            ttd_1 = self.inl[0].T.val_SI - self.Tamb.val_SI
            ttd_2 = self.outl[0].T.val_SI - self.Tamb.val_SI

            if ttd_1 > ttd_2:
                td_log = (ttd_1 - ttd_2) / np.log(ttd_1 / ttd_2)
            elif ttd_1 < ttd_2:
                td_log = (ttd_2 - ttd_1) / np.log(ttd_2 / ttd_1)
            else:
                td_log = np.nan

            self.kA.val = abs(i[0] * (o[2] - i[2]) / td_log)
            self.kA.is_result = True
        else:
            self.kA.is_result = False

    def entropy_balance(self):
        r"""
        Calculate entropy balance of a simple heat exchanger.

        The allocation of the entropy streams due to heat exchanged and due to
        irreversibility is performed by solving for T:

        .. math::

            h_\mathrm{out} - h_\mathrm{in} = \int_\mathrm{out}^\mathrm{in}
            v \cdot dp - \int_\mathrm{out}^\mathrm{in} T \cdot ds

        As solving :math:`\int_\mathrm{out}^\mathrm{in} v \cdot dp` for non
        isobaric processes would require perfect process knowledge (the path)
        on how specific volume and pressure change throught the component, the
        heat transfer is splitted into three separate virtual processes:

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

            \text{S\_Q}=\dot{m} \cdot \left(s_\mathrm{out*}-s_\mathrm{in*}
            \right)\\
            \text{S\_irr}=\dot{m} \cdot \left(s_\mathrm{out}-s_\mathrm{in}
            \right) - \text{S\_Q}\\
            \text{T\_mQ}=\frac{\dot{Q}}{\text{S\_Q}}
        """
        i = self.inl[0].get_flow()
        o = self.outl[0].get_flow()

        p1_star = i[1] * (o[1] / i[1]) ** 0.5
        s1_star = s_mix_ph([0, p1_star, i[2], i[3]], T0=self.inl[0].T.val_SI)
        s2_star = s_mix_ph([0, p1_star, o[2], o[3]], T0=self.outl[0].T.val_SI)
        self.S_Q = i[0] * (s2_star - s1_star)
        self.S_irr = i[0] * (
            self.outl[0].s.val_SI - self.inl[0].s.val_SI) - self.S_Q
        self.T_mQ = (o[2] - i[2]) / (s2_star - s1_star)

    def exergy_balance(self, T0):
        r"""
        Calculate exergy balance of a simple heat exchanger.

        The exergy of heat is calculated by allocation of thermal and
        mechanical share of exergy in the physical exergy. Depending on the
        temperature levels at the inlet and outlet of the heat exchanger as
        well as the direction of heat transfer (input or output) fuel and
        product exergy are calculated as follows.

        Parameters
        ----------
        T0 : float
            Ambient temperature T0 / K.

        Note
        ----
        If the fluid transfers heat to the ambient, you can specify
        :code:`mysimpleheatexchanger.set_attr(dissipative=False)` if you do
        NOT want the exergy production nan (only applicable in case
        :math:`\dot{Q}<0`).

        .. math ::

            \dot{E}_\mathrm{P} =
            \begin{cases}
            \begin{cases}
            \begin{cases}
            \text{not defined (nan)} & \text{if dissipative}\\
            \dot{E}_\mathrm{in}^\mathrm{T} - \dot{E}_\mathrm{out}^\mathrm{T} &
            \text{else}\\
            \end{cases}
            & T_\mathrm{in}, T_\mathrm{out} \geq T_0\\
            \dot{E}_\mathrm{out}^\mathrm{T}
            & T_\mathrm{in} \geq T_0 > T_\mathrm{out}\\
            \dot{E}_\mathrm{out}^\mathrm{T} - \dot{E}_\mathrm{in}^\mathrm{T}
            & T_0 \geq T_\mathrm{in}, T_\mathrm{out}\\
            \end{cases} & \dot{Q} < 0\\

            \begin{cases}
            \dot{E}_\mathrm{out}^\mathrm{PH} - \dot{E}_\mathrm{in}^\mathrm{PH}
            & T_\mathrm{in}, T_\mathrm{out} \geq T_0\\
            \dot{E}_\mathrm{in}^\mathrm{T} + \dot{E}_\mathrm{out}^\mathrm{T}
            & T_\mathrm{out} > T_0 \geq T_\mathrm{in}\\
            \dot{E}_\mathrm{in}^\mathrm{T} - \dot{E}_\mathrm{out}^\mathrm{T} +
            \dot{E}_\mathrm{out}^\mathrm{M} - \dot{E}_\mathrm{in}^\mathrm{M} +
            & T_0 \geq T_\mathrm{in}, T_\mathrm{out}\\
            \end{cases} & \dot{Q} > 0\\
            \end{cases}

            \dot{E}_\mathrm{F} =
            \begin{cases}
            \begin{cases}
            \dot{E}_\mathrm{in}^\mathrm{PH} - \dot{E}_\mathrm{out}^\mathrm{PH}
            & T_\mathrm{in}, T_\mathrm{out} \geq T_0\\
            \dot{E}_\mathrm{in}^\mathrm{T} + \dot{E}_\mathrm{in}^\mathrm{M} +
            \dot{E}_\mathrm{out}^\mathrm{T} - \dot{E}_\mathrm{out}^\mathrm{M}
            & T_\mathrm{in} \geq T_0 > T_\mathrm{out}\\
            \dot{E}_\mathrm{out}^\mathrm{T} - \dot{E}_\mathrm{in}^\mathrm{T} +
            \dot{E}_\mathrm{in}^\mathrm{M} - \dot{E}_\mathrm{out}^\mathrm{M} +
            & T_0 \geq T_\mathrm{in}, T_\mathrm{out}\\
            \end{cases} & \dot{Q} < 0\\

            \begin{cases}
            \dot{E}_\mathrm{out}^\mathrm{T} - \dot{E}_\mathrm{in}^\mathrm{T}
            & T_\mathrm{in}, T_\mathrm{out} \geq T_0\\
            \dot{E}_\mathrm{in}^\mathrm{T} + \dot{E}_\mathrm{in}^\mathrm{M} -
            \dot{E}_\mathrm{out}^\mathrm{M}
            & T_\mathrm{out} > T_0 \geq T_\mathrm{in}\\
            \dot{E}_\mathrm{in}^\mathrm{T}-\dot{E}_\mathrm{out}^\mathrm{T}
            & T_0 \geq T_\mathrm{in}, T_\mathrm{out}\\
            \end{cases} & \dot{Q} > 0\\
            \end{cases}

            \dot{E}_\mathrm{bus} =
            \begin{cases}
            \begin{cases}
            \dot{E}_\mathrm{P} & \text{other cases}\\
            \dot{E}_\mathrm{in}^\mathrm{T}
            & T_\mathrm{in} \geq T_0 > T_\mathrm{out}\\
            \end{cases} & \dot{Q} < 0\\
            \dot{E}_\mathrm{F} & \dot{Q} > 0\\
            \end{cases}
        """
        if self.Q.val < 0:
            if self.inl[0].T.val_SI >= T0 and self.outl[0].T.val_SI >= T0:
                if self.dissipative.val:
                    self.E_P = np.nan
                else:
                    self.E_P = self.inl[0].Ex_therm - self.outl[0].Ex_therm
                self.E_F = self.inl[0].Ex_physical - self.outl[0].Ex_physical
                self.E_bus = self.E_P
            elif self.inl[0].T.val_SI >= T0 and self.outl[0].T.val_SI < T0:
                self.E_P = self.outl[0].Ex_therm
                self.E_F = self.inl[0].Ex_therm + self.outl[0].Ex_therm + (
                    self.inl[0].Ex_mech - self.outl[0].Ex_mech)
                self.E_bus = self.inl[0].Ex_therm + self.outl[0].Ex_therm
            elif self.inl[0].T.val_SI <= T0 and self.outl[0].T.val_SI <= T0:
                self.E_P = self.outl[0].Ex_therm - self.inl[0].Ex_therm
                self.E_F = self.outl[0].Ex_therm - self.outl[0].Ex_therm + (
                    self.inl[0].Ex_mech - self.outl[0].Ex_mech)
                self.E_bus = self.E_P
            else:
                msg = ('Exergy balance of simple heat exchangers, where '
                       'outlet temperature is higher than inlet temperature '
                       'with heat extracted is not implmented.')
                logging.warning(msg)
                self.E_P = np.nan
                self.E_F = np.nan
                self.E_bus = np.nan
        elif self.Q.val > 0:
            if self.inl[0].T.val_SI >= T0 and self.outl[0].T.val_SI >= T0:
                self.E_P = self.outl[0].Ex_physical - self.inl[0].Ex_physical
                self.E_F = self.outl[0].Ex_therm - self.inl[0].Ex_therm
                self.E_bus = self.E_F
            elif self.inl[0].T.val_SI <= T0 and self.outl[0].T.val_SI > T0:
                self.E_P = self.outl[0].Ex_therm + self.inl[0].Ex_therm
                self.E_F = self.inl[0].Ex_therm + (
                    self.inl[0].Ex_mech - self.outl[0].Ex_mech)
                self.E_bus = self.inl[0].Ex_therm
            elif self.inl[0].T.val_SI < T0 and self.outl[0].T.val_SI < T0:
                self.E_P = self.inl[0].Ex_therm - self.outl[0].Ex_therm + (
                    self.outl[0].Ex_mech - self.inl[0].Ex_mech)
                self.E_F = self.inl[0].Ex_therm - self.outl[0].Ex_therm
                self.E_bus = self.E_F
            else:
                msg = ('Exergy balance of simple heat exchangers, where '
                       'inlet temperature is higher than outlet temperature '
                       'with heat injected is not implmented.')
                logging.warning(msg)
                self.E_P = np.nan
                self.E_F = np.nan
                self.E_bus = self.E_F
        else:
            # this is basically the exergy balance of a valve
            self.E_P = np.nan
            self.E_F = self.inl[0].Ex_physical - self.outl[0].Ex_physical
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
            1: {
                'isoline_property': 'p',
                'isoline_value': self.inl[0].p.val,
                'isoline_value_end': self.outl[0].p.val,
                'starting_point_property': 's',
                'starting_point_value': self.inl[0].s.val,
                'ending_point_property': 's',
                'ending_point_value': self.outl[0].s.val
            }
        }
