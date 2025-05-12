# -*- coding: utf-8

"""Module of class ParabolicTrough.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/heat_exchangers/parabolic_trough.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import component_registry
from tespy.components.heat_exchangers.simple import SimpleHeatExchanger
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.document_models import generate_latex_eq


@component_registry
class ParabolicTrough(SimpleHeatExchanger):
    r"""
    The ParabolicTrough calculates heat output from irradiance.

    **Mandatory Equations**

    - :py:meth:`tespy.components.component.Component.fluid_func`
    - :py:meth:`tespy.components.component.Component.mass_flow_func`

    **Optional Equations**

    - :py:meth:`tespy.components.component.Component.pr_func`
    - :py:meth:`tespy.components.component.Component.zeta_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.energy_balance_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.darcy_group_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.hw_group_func`
    - :py:meth:`tespy.components.heat_exchangers.parabolic_trough.ParabolicTrough.energy_group_func`

    Inlets/Outlets

    - in1
    - out1

    Image

    .. image:: /api/_images/ParabolicTrough.svg
       :alt: flowsheet of the parabolic trough
       :align: center
       :class: only-light

    .. image:: /api/_images/ParabolicTrough_darkmode.svg
       :alt: flowsheet of the parabolic trough
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
        Diameter of the absorber tube, :math:`D/\text{m}`.

    L : float, dict, :code:`"var"`
        Length of the absorber tube, :math:`L/\text{m}`.

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
        Direct irradiance to tilted collector,
        :math:`E/\frac{\text{W}}{\text{m}^2}`.

    aoi : float, dict, :code:`"var"`
        Angle of incidience, :math:`aoi/^\circ`.

    doc : float, dict, :code:`"var"`
        Degree of cleanliness (1: full absorption, 0: no absorption),
        :math:`X`.

    eta_opt : float, dict, :code:`"var"`
        (constant) optical losses due to surface reflection,
        :math:`\eta_{opt}`.

    c_1 : float, dict, :code:`"var"`
        Linear thermal loss key figure,
        :math:`c_1/\frac{\text{W}}{\text{K} \cdot \text{m}^2}`.

    c_2 : float, dict, :code:`"var"`
        Quadratic thermal loss key figure,
        :math:`c_2/\frac{\text{W}}{\text{K}^2 \cdot \text{m}^2}`.

    iam_1 : float, dict, :code:`"var"`
        Linear incidence angle modifier,
        :math:`iam_1/\frac{1}{^\circ}`.

    iam_2 : float, dict, :code:`"var"`
        Quadratic incidence angle modifier,
        :math:`iam_2/\left(\frac{1}{^\circ}\right)^2`.

    A : float, dict, :code:`"var"`
        Collector aperture surface area :math:`A/\text{m}^2`.

    Tamb : float, dict
        Ambient temperature, provide parameter in network's temperature unit.

    energy_group : str, dict
        Parametergroup for energy balance of solarthermal collector.

    Example
    -------
    A parabolic trough is installed using S800 as thermo-fluid.
    First, the operation conditions from :cite:`Janotte2014` are reproduced.
    Therefore, the direct normal irradiance :math:`\dot{E}_\mathrm{DNI}` is at
    1000 :math:`\frac{\text{W}}{\text{m}^2}` at an angle of incidence
    :math:`aoi` at 20 °. This means, the direct irradiance to the parabolic
    trough :math:`E` is at
    :math:`\dot{E}_{DNI} \cdot cos\left(20^\circ\right)`.

    >>> from tespy.components import Sink, Source, ParabolicTrough
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> import math
    >>> import shutil
    >>> nw = Network()
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg', iterinfo=False)
    >>> so = Source('source')
    >>> si = Sink('sink')
    >>> pt = ParabolicTrough('parabolic trough collector')
    >>> pt.component()
    'parabolic trough'
    >>> inc = Connection(so, 'out1', pt, 'in1')
    >>> outg = Connection(pt, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)

    The pressure ratio is at a constant level of 1. However, it is possible to
    specify the pressure losses from the absorber tube length, roughness and
    diameter, too. The aperture surface :math:`A` is specified to 1
    :math:`\text{m}^2` for simplicity reasons.

    >>> aoi = 20
    >>> E = 1000 * math.cos(aoi / 180 * math.pi)
    >>> pt.set_attr(
    ... pr=1, aoi=aoi, doc=1,
    ... Tamb=20, A=1, eta_opt=0.816, c_1=0.0622, c_2=0.00023, E=E,
    ... iam_1=-1.59e-3, iam_2=9.77e-5
    ... )
    >>> inc.set_attr(fluid={'INCOMP::S800': 1}, T=220, p=10)
    >>> outg.set_attr(T=260)
    >>> nw.solve('design')
    >>> round(pt.Q.val, 0)
    736.0

    For example, it is possible to calculate the aperture area of the parabolic
    trough given the total heat production, outflow temperature and mass flow.

    >>> pt.set_attr(A='var', Q=5e6, Tamb=25)
    >>> inc.set_attr(T=None)
    >>> outg.set_attr(T=350, m=20)
    >>> nw.solve('design')
    >>> round(inc.T.val, 0)
    229.0
    >>> round(pt.A.val, 0)
    6862.0

    Given this design, it is possible to calculate the outlet temperature as
    well as the heat transfer at different operating points.

    >>> aoi = 30
    >>> E = 800 * math.cos(aoi / 180 * math.pi)
    >>> pt.set_attr(A=pt.A.val, aoi=aoi, Q=None, E=E)
    >>> inc.set_attr(T=150)
    >>> outg.set_attr(T=None)
    >>> nw.solve('design')
    >>> round(outg.T.val, 0)
    244.0
    >>> round(pt.Q.val, 0)
    3602817.0
    """

    @staticmethod
    def component():
        return 'parabolic trough'

    def get_parameters(self):
        data = super().get_parameters()
        for k in ["kA_group", "kA_char_group", "kA", "kA_char"]:
            del data[k]

        data.update({
            'E': dc_cp(min_val=0), 'A': dc_cp(min_val=0),
            'eta_opt': dc_cp(min_val=0, max_val=1),
            'c_1': dc_cp(min_val=0), 'c_2': dc_cp(min_val=0),
            'iam_1': dc_cp(), 'iam_2': dc_cp(),
            'aoi': dc_cp(min_val=-90, max_val=90),
            'doc': dc_cp(min_val=0, max_val=1),
            'Tamb': dc_cp(),
            'Q_loss': dc_cp(max_val=0, val=0),
            'energy_group': dc_gcp(
                elements=[
                    'E', 'eta_opt', 'aoi', 'doc', 'c_1', 'c_2', 'iam_1',
                    'iam_2', 'A', 'Tamb'
                ],
                num_eq=1,
                latex=self.energy_group_func_doc,
                func=self.energy_group_func,
                deriv=self.energy_group_deriv
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
                T_m = & \frac{T_{out} + T_{in}}{2}\\
                iam = & 1 - iam_1 \cdot |aoi| - iam_2 \cdot aoi^2\\
                0 = & \dot{m} \cdot \left( h_{out} - h_{in} \right)\\
                & - A \cdot \left[E \cdot \eta_{opt} \cdot doc^{1.5} \cdot
                iam \right. \\
                & \left. - c_1 \cdot \left(T_m - T_{amb} \right) -
                c_2 \cdot \left(T_m - T_{amb}\right)^2
                \vphantom{ \eta_{opt} \cdot doc^{1.5}} \right]
                \end{split}

            Reference: :cite:`Janotte2014`.
        """
        i = self.inl[0]
        o = self.outl[0]

        T_m = 0.5 * (i.calc_T() + o.calc_T())

        iam = (
            1 - self.iam_1.val * abs(self.aoi.val)
            - self.iam_2.val * self.aoi.val ** 2
        )

        return (
            i.m.val_SI * (o.h.val_SI - i.h.val_SI) - self.A.val * (
                self.E.val * self.eta_opt.val * self.doc.val ** 1.5 * iam
                - (T_m - self.Tamb.val_SI) * self.c_1.val
                - self.c_2.val * (T_m - self.Tamb.val_SI) ** 2
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
            r'& - A \cdot \left[E \cdot \eta_\mathrm{opt} \cdot doc^{1.5}'
            r'\cdot iam \right. \\' + '\n'
            r'&\left. -c_1\cdot\left(T_\mathrm{m}-T_\mathrm{amb}\right) -'
            r'c_2 \cdot \left(T_\mathrm{m} - T_\mathrm{amb}\right)^2'
            r'\vphantom{\eta_\mathrm{opt}\cdot doc^{1.5}}\right]\\' + '\n'
            r'T_\mathrm{m}=&\frac{T_\mathrm{out}+T_\mathrm{in}}{2}\\' +
            '\n'
            r'iam = & 1 - iam_1 \cdot |aoi| - iam_2 \cdot aoi^2\\' + '\n'
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
            self.Q_loss.val = - self.E.val * self.A.val + self.Q.val
            self.Q_loss.is_result = True
        else:
            self.Q_loss.is_result = False
