# -*- coding: utf-8

"""Module of class SolarCollector.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/heat_exchangers/solar_collector.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import component_registry
from tespy.components.heat_exchangers.simple import SimpleHeatExchanger
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.document_models import generate_latex_eq


@component_registry
class SolarCollector(SimpleHeatExchanger):
    r"""
    The solar collector calculates heat output from irradiance.

    **Mandatory Equations**

    - :py:meth:`tespy.components.component.Component.fluid_func`
    - :py:meth:`tespy.components.component.Component.mass_flow_func`

    **Optional Equations**

    - :py:meth:`tespy.components.component.Component.pr_func`
    - :py:meth:`tespy.components.component.Component.zeta_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.energy_balance_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.darcy_group_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.hw_group_func`
    - :py:meth:`tespy.components.heat_exchangers.solar_collector.SolarCollector.energy_group_func`

    Inlets/Outlets

    - in1
    - out1

    Image

    .. image:: /api/_images/SolarCollector.svg
       :alt: flowsheet of the solar collector
       :align: center
       :class: only-light

    .. image:: /api/_images/SolarCollector_darkmode.svg
       :alt: flowsheet of the solar collector
       :align: center
       :class: only-dark

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
        Pipe's roughness, :math:`ks/\text{m}`.

    darcy_group : str, dict
        Parametergroup for pressure drop calculation based on pipes dimensions
        using darcy weissbach equation.

    ks_HW : float, dict, :code:`"var"`
        Pipe's roughness, :math:`ks/\text{1}`.

    hw_group : str, dict
        Parametergroup for pressure drop calculation based on pipes dimensions
        using hazen williams equation.

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
    >>> nw = Network()
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
    >>> nw.save('tmp.json')
    >>> round(sc.A.val, 1)
    14.5
    >>> sc.set_attr(A=sc.A.val, E=5e2, Tamb=20)
    >>> nw.solve('offdesign', design_path='tmp.json')
    >>> round(sc.Q.val, 1)
    6083.8
    >>> round(outg.T.val, 1)
    70.5
    >>> shutil.rmtree('./tmp.json', ignore_errors=True)
    """

    @staticmethod
    def component():
        return 'solar collector'

    def get_parameters(self):
        data = super().get_parameters()
        for k in ["kA_group", "kA_char_group", "kA", "kA_char"]:
            del data[k]

        data.update({
            'E': dc_cp(min_val=0), 'A': dc_cp(min_val=0),
            'eta_opt': dc_cp(min_val=0, max_val=1),
            'lkf_lin': dc_cp(min_val=0),
            'lkf_quad': dc_cp(min_val=0),
            'Tamb': dc_cp(),
            'Q_loss': dc_cp(max_val=0, val=0),
            'energy_group': dc_gcp(
                elements=['E', 'eta_opt', 'lkf_lin', 'lkf_quad', 'A', 'Tamb'],
                num_eq=1, latex=self.energy_group_func_doc,
                func=self.energy_group_func, deriv=self.energy_group_deriv
            )
        })
        return data

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
        i = self.inl[0]
        o = self.outl[0]

        T_m = 0.5 * (i.calc_T() + o.calc_T())

        return (
            i.m.val_SI * (o.h.val_SI - i.h.val_SI)
            - self.A.val * (
                self.E.val * self.eta_opt.val
                - (T_m - self.Tamb.val_SI) * self.lkf_lin.val
                - self.lkf_quad.val * (T_m - self.Tamb.val_SI) ** 2
            )
        )

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
        i = self.inl[0]
        o = self.outl[0]
        if self.is_variable(i.m, increment_filter):
            self.jacobian[k, i.m.J_col] = o.h.val_SI - i.h.val_SI
        if self.is_variable(i.p, increment_filter):
            self.jacobian[k, i.p.J_col] = self.numeric_deriv(f, 'p', i)
        if self.is_variable(i.h, increment_filter):
            self.jacobian[k, i.h.J_col] = self.numeric_deriv(f, 'h', i)
        if self.is_variable(o.p, increment_filter):
            self.jacobian[k, o.p.J_col] = self.numeric_deriv(f, 'p', o)
        if self.is_variable(o.h, increment_filter):
            self.jacobian[k, o.h.J_col] = self.numeric_deriv(f, 'h', o)
        # custom variables for the energy-group
        for variable_name in self.energy_group.elements:
            parameter = self.get_attr(variable_name)
            if parameter == self.Tamb:
                continue
            if parameter.is_var:
                self.jacobian[k, parameter.J_col] = (
                    self.numeric_deriv(f, variable_name, None)
                )

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        i = self.inl[0]
        o = self.outl[0]

        self.Q.val = i.m.val_SI * (o.h.val_SI - i.h.val_SI)
        self.pr.val = o.p.val_SI / i.p.val_SI
        self.zeta.val = self.calc_zeta(i, o)

        if self.energy_group.is_set:
            self.Q_loss.val = -(self.E.val * self.A.val - self.Q.val)
            self.Q_loss.is_result = True
        else:
            self.Q_loss.is_result = False
