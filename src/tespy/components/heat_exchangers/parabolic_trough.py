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


@component_registry
class ParabolicTrough(SimpleHeatExchanger):
    r"""
    The ParabolicTrough calculates heat output from irradiance.

    .. image:: /api/_images/components/ParabolicTrough.svg
       :alt: flowsheet of the parabolictrough
       :align: center
       :class: only-light

    .. image:: /api/_images/components/ParabolicTrough_darkmode.svg
       :alt: flowsheet of the parabolictrough
       :align: center
       :class: only-dark

    Ports
    -----

    Fluid inlets: in1

    Fluid outlets: out1

    Power inlets: heat

    Power outlets: heat

    Heat inlets: heat

    Heat outlets: heat

    Mandatory Equations
    -------------------

    - mass flow equality constraint(s): :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
    - fluid composition equality constraint(s): :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`

    When a power or heat connector is attached:

    - energy_connector_balance: :py:meth:`energy_connector_balance_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.energy_connector_balance_func>`

    Parameters
    ----------

    A : float, dict, :code:`"var"`
        Area of the parabolic trough. Quantity: :code:`area`. Can be set as a
        system variable by passing :code:`"var"` as its value.

    aoi : float, dict
        Angle of incidence. Quantity: :code:`angle`.

    c_1 : float, dict
        Thermal loss coefficient 1.

    c_2 : float, dict
        Thermal loss coefficient 2.

    char_warnings : bool
        Ignore warnings on default characteristics usage for this component.

    D : float, dict, :code:`"var"`
        Diameter of channel. Quantity: :code:`length`. Can be set as a system
        variable by passing :code:`"var"` as its value.

    darcy_group : GroupedComponentProperties
        Darcy-Weißbach equation for pressure loss. Elements: :code:`L`,
        :code:`ks`, :code:`D`.
        Equation: :py:meth:`darcy_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.darcy_func>`.

    design : list
        List containing design parameters (stated as String).

    design_path : str
        Path to the components design case.

    dissipative : bool


    doc : float, dict
        Degree of cleanliness. Quantity: :code:`ratio`.

    dp : float, dict
        Inlet to outlet absolute pressure change. Quantity:
        :code:`pressure_difference`.
        Equation: :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`.

    E : float, dict, :code:`"var"`
        Solar irradiation to the parabolic trough. Quantity: :code:`heat`. Can
        be set as a system variable by passing :code:`"var"` as its value.

    energy_group : GroupedComponentProperties
        Energy balance equation of the parabolic trough. Elements: :code:`E`,
        :code:`eta_opt`, :code:`aoi`, :code:`doc`, :code:`c_1`, :code:`c_2`,
        :code:`iam_1`, :code:`iam_2`, :code:`A`, :code:`Tamb`.
        Equation: :py:meth:`energy_group_func <tespy.components.heat_exchangers.parabolic_trough.ParabolicTrough.energy_group_func>`.

    eta_opt : float, dict
        Optical efficiency. Quantity: :code:`efficiency`.

    hw_group : GroupedComponentProperties
        Hazen-Williams equation for pressure loss. Elements: :code:`L`,
        :code:`ks_HW`, :code:`D`.
        Equation: :py:meth:`hazen_williams_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.hazen_williams_func>`.

    iam_1 : float, dict
        Incidence angle modifier 1.

    iam_2 : float, dict
        Incidence angle modifier 2.

    ks : float, dict, :code:`"var"`
        Roughness of wall material. Quantity: :code:`length`. Can be set as a
        system variable by passing :code:`"var"` as its value.

    ks_HW : float, dict, :code:`"var"`
        Hazen-Williams roughness. Can be set as a system variable by passing
        :code:`"var"` as its value.

    L : float, dict, :code:`"var"`
        Length of channel. Quantity: :code:`length`. Can be set as a system
        variable by passing :code:`"var"` as its value.

    label : str
        The label of the component.

    local_design : bool
        Treat this component in design mode in an offdesign calculation.

    local_offdesign : bool
        Treat this component in offdesign mode in a design calculation.

    offdesign : list
        List containing offdesign parameters (stated as String).

    power_connector_location : str


    pr : float, dict
        Outlet to inlet pressure ratio. Quantity: :code:`ratio`.
        Equation: :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`.

    printout : bool
        Include this component in the network's results printout.

    Q : float, dict
        Heat transfer. Quantity: :code:`heat`.
        Equation: :py:meth:`energy_balance_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.energy_balance_func>`.

    Q_loss : float, dict
        Heat dissipation. Quantity: :code:`heat`.

    Tamb : float, dict
        Ambient temperature. Quantity: :code:`temperature`.

    zeta : float, dict
        Non-dimensional friction coefficient for pressure loss calculation.
        Equation: :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`.

    Example
    -------
    A parabolic trough is installed using S800 as thermo-fluid.
    First, the operation conditions from :cite:`Janotte2014` are reproduced.
    Therefore, the direct normal irradiance :math:`\dot{E}_\text{DNI}` is at
    1000 :math:`\frac{\text{W}}{\text{m}^2}` at an angle of incidence
    :math:`aoi` at 20 °. This means, the direct irradiance to the parabolic
    trough :math:`E` is at
    :math:`\dot{E}_{DNI} \cdot cos\left(20^\circ\right)`.

    >>> from tespy.components import Sink, Source, ParabolicTrough
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> import math
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "pressure_difference": "bar",
    ...     "temperature": "degC", "enthalpy": "kJ/kg"
    ... })
    >>> so = Source('source')
    >>> si = Sink('sink')
    >>> pt = ParabolicTrough('parabolic trough collector')
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
    ...     pr=1, aoi=aoi, doc=1,
    ...     Tamb=20, A=1, eta_opt=0.816, c_1=0.0622, c_2=0.00023, E=E,
    ...     iam_1=-1.59e-3, iam_2=9.77e-5
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

    def get_parameters(self):
        data = super().get_parameters()
        for k in ["kA_group", "kA_char_group", "kA", "kA_char"]:
            del data[k]

        data.update({
            'E': dc_cp(
                min_val=0, quantity="heat", _allows_var=True,
                description="solar irradiation to the parabolic trough"
            ),
            'A': dc_cp(
                min_val=0, quantity="area", _allows_var=True,
                description="area of the parabolic trough"
            ),
            'eta_opt': dc_cp(
                min_val=0, max_val=1, quantity="efficiency",
                description="optical efficiency"
            ),
            'c_1': dc_cp(min_val=0, description="thermal loss coefficient 1"),
            'c_2': dc_cp(min_val=0, description="thermal loss coefficient 2"),
            'iam_1': dc_cp(description="incidence angle modifier 1"),
            'iam_2': dc_cp(description="incidence angle modifier 2"),
            'aoi': dc_cp(
                min_val=-90, max_val=90, quantity="angle",
                description="angle of incidence"
            ),
            'doc': dc_cp(
                min_val=0, max_val=1, quantity="ratio",
                description="degree of cleanliness"
            ),
            'Q_loss': dc_cp(
                max_val=0, _val=0, quantity="heat",
                description="heat dissipation",
                calc=self._calc_Q_loss, calc_deps=['Q']
            ),
            'energy_group': dc_gcp(
                elements=[
                    'E', 'eta_opt', 'aoi', 'doc', 'c_1', 'c_2', 'iam_1',
                    'iam_2', 'A', 'Tamb'
                ],
                num_eq_sets=1,
                func=self.energy_group_func,
                dependents=self.energy_group_dependents,
                description="energy balance equation of the parabolic trough"
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
            1 - self.iam_1.val_SI * abs(self.aoi.val_SI)
            - self.iam_2.val_SI * self.aoi.val_SI ** 2
        )

        return (
            i.m.val_SI * (o.h.val_SI - i.h.val_SI) - self.A.val_SI * (
                self.E.val_SI * self.eta_opt.val_SI * self.doc.val_SI ** 1.5 * iam
                - self.c_1.val_SI * (T_m - self.Tamb.val_SI)
                - self.c_2.val_SI * (T_m - self.Tamb.val_SI) ** 2
            )
        )

    def energy_group_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].p,
            self.outl[0].h,
        ] + [self.E, self.A]

    def convergence_check(self):
        pass

    def _calc_Q_loss(self):
        return -self.E.val_SI * self.A.val_SI + self.Q.val_SI

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        super().calc_parameters()
        self.Q_loss.is_result = self.energy_group.is_set
