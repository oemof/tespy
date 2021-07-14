# -*- coding: utf-8

"""Module of class SolarCollector.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/heat_exchangers/solar_collector.py

SPDX-License-Identifier: MIT
"""

import numpy as np

from tespy.components.heat_exchangers.heat_exchanger_simple import HeatExchangerSimple
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import DataContainerSimple as dc_simple
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.document_models import generate_latex_eq
from tespy.tools.fluid_properties import T_mix_ph


class SolarCollector(HeatExchangerSimple):
    r"""
    The solar collector calculates heat output from irradiance.

    **Mandatory Equations**

    - :py:meth:`tespy.components.component.Component.fluid_func`
    - :py:meth:`tespy.components.component.Component.mass_flow_func`

    **Optional Equations**

    - :py:meth:`tespy.components.component.Component.pr_func`
    - :py:meth:`tespy.components.component.Component.zeta_func`
    - :py:meth:`tespy.components.heat_exchangers.heat_exchanger_simple.HeatExchangerSimple.energy_balance_func`
    - :py:meth:`tespy.components.heat_exchangers.heat_exchanger_simple.HeatExchangerSimple.hydro_group_func`
    - :py:meth:`tespy.components.heat_exchangers.solar_collector.SolarCollector.energy_group_func`

    Inlets/Outlets

    - in1
    - out1

    Image

    .. image:: _images/SolarCollector.svg
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

    E : float, dict, :code:`"var"`
        irradiance at tilted collector surface area,
        :math:`E/\frac{\text{W}}{\text{m}^2}`.

    eta_opt : float, dict, :code:`"var"`
        optical loss at surface cover,
        :math:`\eta_{opt}`.

    lkf_lin : float, dict, :code:`"var"`
        Linear thermal loss key figure,
        :math:`\alpha_1/\frac{\text{W}}{\text{K} \cdot \text{m}^2}`.

    lkf_quad : float, dict, :code:`"var"`
        Quadratic thermal loss key figure,
        :math:`\alpha_2/\frac{\text{W}}{\text{K}^2 \cdot \text{m}^2}`.

    A : float, dict, :code:`"var"`
        Collector surface area :math:`A/\text{m}^2`.

    Tamb : float, dict
        Ambient temperature, provide parameter in network's temperature unit.

    energy_group : str, dict
        Parametergroup for energy balance of solarthermal collector.

    Example
    -------
    The solar collector is used to calculate heat transferred to the heating
    system from irradiance on a tilted plane. For instance, it is possible to
    calculate the collector surface area required to transfer a specific amount
    of heat at a given irradiance. The collector parameters are the linear and
    the quadratic loss keyfigure as well as the optical effifiency.

    >>> from tespy.components import Sink, Source, SolarCollector
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> import shutil
    >>> fluids = ['H2O']
    >>> nw = Network(fluids=fluids)
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg', iterinfo=False)
    >>> so = Source('source')
    >>> si = Sink('sink')
    >>> sc = SolarCollector('solar collector')
    >>> sc.component()
    'solar collector'
    >>> sc.set_attr(pr=0.95, Q=1e4, design=['pr', 'Q'], offdesign=['zeta'],
    ...     Tamb=25, A='var', eta_opt=0.92, lkf_lin=1, lkf_quad=0.005, E=8e2)
    >>> inc = Connection(so, 'out1', sc, 'in1')
    >>> outg = Connection(sc, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)

    The outlet temperature should be at 90 °C at a constant mass flow, which
    is determined in the design calculation. In offdesign operation (at a
    different irradiance) using the calculated surface area and mass flow, it
    is possible to predict the outlet temperature. It would instead be
    possible to calulate the change in mass flow required to hold the
    specified outlet temperature, too.

    >>> inc.set_attr(fluid={'H2O': 1}, T=40, p=3, offdesign=['m'])
    >>> outg.set_attr(T=90, design=['T'])
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(sc.A.val, 1)
    14.5
    >>> sc.set_attr(A=sc.A.val, E=5e2, Tamb=20)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(sc.Q.val, 1)
    6083.8
    >>> round(outg.T.val, 1)
    70.5
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    @staticmethod
    def component():
        return 'solar collector'

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
            'E': dc_cp(min_val=0), 'A': dc_cp(min_val=0),
            'eta_opt': dc_cp(min_val=0, max_val=1),
            'lkf_lin': dc_cp(min_val=0), 'lkf_quad': dc_cp(min_val=0),
            'Tamb': dc_cp(),
            'Q_loss': dc_cp(max_val=0, val=0),
            'dissipative': dc_simple(val=True),
            'hydro_group': dc_gcp(
                elements=['L', 'ks', 'D'], num_eq=1,
                latex=self.hydro_group_func_doc,
                func=self.hydro_group_func, deriv=self.hydro_group_deriv),
            'energy_group': dc_gcp(
                elements=['E', 'eta_opt', 'lkf_lin', 'lkf_quad', 'A', 'Tamb'],
                num_eq=1, latex=self.energy_group_func_doc,
                func=self.energy_group_func, deriv=self.energy_group_deriv)
        }

    def energy_group_func(self):
        r"""
        Equation for solar collector energy balance.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                \begin{split}
                0 = & \dot{m} \cdot \left( h_{out} - h_{in} \right)\\
                & - A \cdot \left[E \cdot \eta_{opt} - \alpha_1 \cdot
                \left(T_m - T_{amb} \right) - \alpha_2 \cdot
                \left(T_m - T_{amb}\right)^2 \right]\\
                T_m = & \frac{T_{out} + T_{in}}{2}\\
                \end{split}

            Reference: :cite:`Quaschning2013`.
        """
        i = self.inl[0].get_flow()
        o = self.outl[0].get_flow()

        T_m = (T_mix_ph(i, T0=self.inl[0].T.val_SI) +
               T_mix_ph(o, T0=self.outl[0].T.val_SI)) / 2

        return (i[0] * (o[2] - i[2]) -
                self.A.val * (
                    self.E.val * self.eta_opt.val -
                    (T_m - self.Tamb.val_SI) * self.lkf_lin.val -
                    self.lkf_quad.val * (T_m - self.Tamb.val_SI) ** 2))

    def energy_group_func_doc(self, label):
        r"""
        Equation for solar collector energy balance.

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
            r'0 = & \dot{m}_\mathrm{in} \cdot \left( h_\mathrm{out} - '
            r'h_\mathrm{in} \right)\\' + '\n'
            r'& - A \cdot \left[E \cdot \eta_\mathrm{opt} - \alpha_1 \cdot'
            r'\left(T_\mathrm{m} - T_\mathrm{amb} \right) - \alpha_2 \cdot'
            r'\left(T_\mathrm{m} -T_\mathrm{amb}\right)^2 \right]\\' + '\n'
            r'T_\mathrm{m}=&\frac{T_\mathrm{out}+T_\mathrm{in}}{2}\\' +
            '\n'
            r'\end{split}'
        )
        return generate_latex_eq(self, latex, label)

    def energy_group_deriv(self, increment_filter, k):
        r"""
        Calculate partial derivatives of energy group function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        f = self.energy_group_func
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
        # custom variables for the energy-group
        for var in self.energy_group.elements:
            var = self.get_attr(var)
            if var == self.Tamb:
                continue
            if var.is_var:
                self.jacobian[k, 2 + var.var_pos, 0] = (
                    self.numeric_deriv(f, self.vars[var], 2))

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        i = self.inl[0].get_flow()
        o = self.outl[0].get_flow()

        self.Q.val = i[0] * (o[2] - i[2])
        self.pr.val = o[1] / i[1]
        self.zeta.val = ((i[1] - o[1]) * np.pi ** 2 / (
            4 * i[0] ** 2 * (self.inl[0].vol.val_SI + self.outl[0].vol.val_SI)
            ))
        if self.energy_group.is_set:
            self.Q_loss.val = -(self.E.val * self.A.val - self.Q.val)
            self.Q_loss.is_result = True
        else:
            self.Q_loss.is_result = False
