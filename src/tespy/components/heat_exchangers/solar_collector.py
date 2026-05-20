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


@component_registry
class SolarCollector(SimpleHeatExchanger):
    r"""
    The solar collector calculates heat output from irradiance.

    .. image:: /api/_images/components/SolarCollector.svg
       :alt: flowsheet of the solarcollector
       :align: center
       :class: only-light

    .. image:: /api/_images/components/SolarCollector_darkmode.svg
       :alt: flowsheet of the solarcollector
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
        Area of the solar collector. Quantity: :code:`area`. Can be set as a
        system variable by passing :code:`"var"` as its value.

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


    dp : float, dict
        Inlet to outlet absolute pressure change. Quantity:
        :code:`pressure_difference`.
        Equation: :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`.

    E : float, dict, :code:`"var"`
        Solar irradiation to the solar collector. Quantity: :code:`heat`. Can be
        set as a system variable by passing :code:`"var"` as its value.

    energy_group : GroupedComponentProperties
        Energy balance equation of the solar collector. Elements: :code:`E`,
        :code:`eta_opt`, :code:`lkf_lin`, :code:`lkf_quad`, :code:`A`,
        :code:`Tamb`.
        Equation: :py:meth:`energy_group_func <tespy.components.heat_exchangers.solar_collector.SolarCollector.energy_group_func>`.

    eta_opt : float, dict
        Optical efficiency. Quantity: :code:`efficiency`.

    hw_group : GroupedComponentProperties
        Hazen-Williams equation for pressure loss. Elements: :code:`L`,
        :code:`ks_HW`, :code:`D`.
        Equation: :py:meth:`hazen_williams_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.hazen_williams_func>`.

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

    lkf_lin : float, dict
        Linear heat loss factor.

    lkf_quad : float, dict
        Quadratic heat loss factor.

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
        Ambient air temperature. Quantity: :code:`temperature`.

    zeta : float, dict
        Non-dimensional friction coefficient for pressure loss calculation.
        Equation: :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`.

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
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "pressure_difference": "bar",
    ...     "temperature": "degC", "enthalpy": "kJ/kg"
    ... })
    >>> so = Source('source')
    >>> si = Sink('sink')
    >>> sc = SolarCollector('solar collector')
    >>> sc.set_attr(pr=0.95, Q=1e4, design=['pr', 'Q'], offdesign=['zeta'],
    ...     Tamb=25, A='var', eta_opt=0.92, lkf_lin=1, lkf_quad=0.005, E=8e2)
    >>> inc = Connection(so, 'out1', sc, 'in1')
    >>> outg = Connection(sc, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)

    The outlet temperature should be at 90 °C at a constant mass flow, which
    is determined in the design calculation. In offdesign operation (at a
    different irradiance) using the calculated surface area and mass flow, it
    is possible to predict the outlet temperature. It would instead be
    possible to calculate the change in mass flow required to hold the
    specified outlet temperature, too.

    >>> inc.set_attr(fluid={'H2O': 1}, T=40, p=3, offdesign=['m'])
    >>> outg.set_attr(T=90, design=['T'])
    >>> nw.solve('design')
    >>> design_state = nw.save(as_dict=True)
    >>> round(sc.A.val, 1)
    14.5
    >>> sc.set_attr(A=sc.A.val, E=5e2, Tamb=20)
    >>> nw.solve('offdesign', design_path=design_state)
    >>> round(sc.Q.val, 1)
    6083.8
    >>> round(outg.T.val, 1)
    70.5
    """

    def get_parameters(self):
        data = super().get_parameters()
        for k in ["kA_group", "kA_char_group", "kA", "kA_char"]:
            del data[k]

        data.update({
            'E': dc_cp(
                min_val=0, quantity="heat", _allows_var=True,
                description="solar irradiation to the solar collector"
            ),
            'A': dc_cp(
                min_val=0, quantity="area", _allows_var=True,
                description="area of the solar collector"
            ),
            'eta_opt': dc_cp(
                min_val=0, max_val=1, quantity="efficiency",
                description="optical efficiency"
            ),
            'lkf_lin': dc_cp(
                min_val=0,
                description="linear heat loss factor"
            ),
            'lkf_quad': dc_cp(
                min_val=0,
                description="quadratic heat loss factor"
            ),
            'Tamb': dc_cp(
                quantity="temperature",
                description="ambient air temperature"
            ),
            'Q_loss': dc_cp(
                max_val=0, _val=0, quantity="heat",
                description="heat dissipation",
                calc=self._calc_Q_loss, calc_deps=['Q']
            ),
            'energy_group': dc_gcp(
                elements=['E', 'eta_opt', 'lkf_lin', 'lkf_quad', 'A', 'Tamb'],
                num_eq_sets=1,
                func=self.energy_group_func,
                dependents=self.energy_group_dependents,
                description="energy balance equation of the solar collector"
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
            - self.A.val_SI * (
                self.E.val_SI * self.eta_opt.val_SI
                - self.lkf_lin.val_SI * (T_m - self.Tamb.val_SI)
                - self.lkf_quad.val_SI * (T_m - self.Tamb.val_SI) ** 2
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
        return -(self.E.val_SI * self.A.val_SI - self.Q.val_SI)

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        super().calc_parameters()
        self.Q_loss.is_result = self.energy_group.is_set
