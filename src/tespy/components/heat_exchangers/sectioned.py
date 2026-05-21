# -*- coding: utf-8

"""Module of class SectionedHeatExchanger.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/heat_exchangers/sectioned.py

SPDX-License-Identifier: MIT
"""
import math

import numpy as np
from scipy.optimize import brentq

from tespy.components.component import component_registry
from tespy.components.heat_exchangers.base import HeatExchanger
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentCharacteristics as dc_gcc
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.fluid_properties import h_mix_pQ
from tespy.tools.fluid_properties import phase_mix_ph
from tespy.tools.fluid_properties import single_fluid


@component_registry
class SectionedHeatExchanger(HeatExchanger):
    r"""
    Class for counter flow heat exchanger with UA sections.

    The heat exchanger is internally discretized into 51 sections of equal heat
    transfer. The number of section can be adjusted by the user. It is based on
    the model implemented by :cite:`Quoilin2020`.

    .. image:: /api/_images/components/HeatExchanger.svg
       :alt: flowsheet of the sectionedheatexchanger
       :align: center
       :class: only-light

    .. image:: /api/_images/components/HeatExchanger_darkmode.svg
       :alt: flowsheet of the sectionedheatexchanger
       :align: center
       :class: only-dark

    Ports
    -----

    Fluid inlets: in1, in2

    Fluid outlets: out1, out2

    Mandatory Equations
    -------------------

    - mass flow equality constraint(s): :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
    - fluid composition equality constraint(s): :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
    - hot side to cold side heat transfer equation: :py:meth:`energy_balance_func <tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_func>`

    Parameters
    ----------

    alpha_ratio : float, dict
        Secondary to refrigerant side convective heat transfer coefficient
        ratio. Quantity: :code:`ratio`.

    area_ratio : float, dict
        Secondary to refrigerant side heat transfer area ratio. Quantity:
        :code:`ratio`.

    char_warnings : bool
        Ignore warnings on default characteristics usage for this component.

    design : list
        List containing design parameters (stated as String).

    design_path : str
        Path to the components design case.

    dp1 : float, dict
        Hot side inlet to outlet absolute pressure change. Quantity:
        :code:`pressure_difference`.
        Equation: :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`.

    dp2 : float, dict
        Cold side inlet to outlet absolute pressure change. Quantity:
        :code:`pressure_difference`.
        Equation: :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`.

    eff_cold : float, dict
        Heat exchanger effectiveness for cold side. Quantity:
        :code:`efficiency`.
        Equation: :py:meth:`eff_cold_func <tespy.components.heat_exchangers.base.HeatExchanger.eff_cold_func>`.

    eff_hot : float, dict
        Heat exchanger effectiveness for hot side. Quantity: :code:`efficiency`.
        Equation: :py:meth:`eff_hot_func <tespy.components.heat_exchangers.base.HeatExchanger.eff_hot_func>`.

    eff_max : float, dict
        Maximum heat exchanger effectiveness. Quantity: :code:`efficiency`.
        Equation: :py:meth:`eff_max_func <tespy.components.heat_exchangers.base.HeatExchanger.eff_max_func>`.

    kA : float, dict
        Heat transfer coefficient considering terminal temperature differences.
        Quantity: :code:`heat_transfer_coefficient`.
        Equation: :py:meth:`kA_func <tespy.components.heat_exchangers.base.HeatExchanger.kA_func>`.

    kA_char : GroupedComponentCharacteristics
        Equation for heat transfer based on kA and modification factor.
        Elements: :code:`kA_char1`, :code:`kA_char2`.
        Equation: :py:meth:`kA_char_func <tespy.components.heat_exchangers.base.HeatExchanger.kA_char_func>`.

    kA_char1 : tespy.tools.characteristics.CharLine, dict
        Hot side kA modification lookup table for offdesign.

    kA_char2 : tespy.tools.characteristics.CharLine, dict
        Cold side kA modification lookup table for offdesign.

    label : str
        The label of the component.

    local_design : bool
        Treat this component in design mode in an offdesign calculation.

    local_offdesign : bool
        Treat this component in offdesign mode in a design calculation.

    num_sections : int
        Number of sections of the heat exchanger.

    offdesign : list
        List containing offdesign parameters (stated as String).

    pr1 : float, dict
        Hot side outlet to inlet pressure ratio. Quantity: :code:`ratio`.
        Equation: :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`.

    pr2 : float, dict
        Cold side outlet to inlet pressure ratio. Quantity: :code:`ratio`.
        Equation: :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`.

    printout : bool
        Include this component in the network's results printout.

    Q : float, dict
        Heat transfer from hot side. Quantity: :code:`heat`.
        Equation: :py:meth:`energy_balance_hot_func <tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_hot_func>`.

    re_exp_r : float, dict
        Reynolds exponent for UA modification based on refrigerant side mass
        flow.

    re_exp_sf : float, dict
        Reynolds exponent for UA modification based on secondary fluid side mass
        flow.

    refrigerant_index : int
        Side on which the refrigerant is flowing (0: hot, 1:cold).

    td_log : float, dict
        Logarithmic temperature difference. Quantity:
        :code:`temperature_difference`.

    td_pinch : float, dict
        Equation for minimum pinch. Quantity: :code:`temperature_difference`.
        Equation: :py:meth:`td_pinch_func <tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger.td_pinch_func>`.

    ttd_l : float, dict
        Terminal temperature difference at hot side outlet to cold side inlet.
        Quantity: :code:`temperature_difference`.
        Equation: :py:meth:`ttd_l_func <tespy.components.heat_exchangers.base.HeatExchanger.ttd_l_func>`.

    ttd_min : float, dict
        Minimum terminal temperature difference. Quantity:
        :code:`temperature_difference`.
        Equation: :py:meth:`ttd_min_func <tespy.components.heat_exchangers.base.HeatExchanger.ttd_min_func>`.

    ttd_u : float, dict
        Terminal temperature difference at hot side inlet to cold side outlet.
        Quantity: :code:`temperature_difference`.
        Equation: :py:meth:`ttd_u_func <tespy.components.heat_exchangers.base.HeatExchanger.ttd_u_func>`.

    UA : float, dict
        Sum of UA values of all sections of heat exchanger. Quantity:
        :code:`heat_transfer_coefficient`.
        Equation: :py:meth:`UA_func <tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger.UA_func>`.

    UA_cecchinato : GroupedComponentProperties
        Equation for UA modification in offdesign. Elements: :code:`re_exp_r`,
        :code:`re_exp_sf`, :code:`alpha_ratio`, :code:`area_ratio`.
        Equation: :py:meth:`UA_cecchinato_func <tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger.UA_cecchinato_func>`.

    UA_char : GroupedComponentCharacteristics
        Equation for sectioned UA modification based on characteristic lines.
        Elements: :code:`kA_char1`, :code:`kA_char2`.
        Equation: :py:meth:`UA_char_func <tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger.UA_char_func>`.

    zeta1 : float, dict
        Hot side non-dimensional friction coefficient for pressure loss
        calculation.
        Equation: :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`.

    zeta2 : float, dict
        Cold side non-dimensional friction coefficient for pressure loss
        calculation.
        Equation: :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`.

    Notes
    -----

    .. note::

        The equations only apply to counter-current heat exchangers.

    Example
    -------
    Water vapor should be cooled down, condensed and then further subcooled.
    For this, air is heated up from 15 °C to 25 °C.

    >>> from tespy.components import Source, Sink, SectionedHeatExchanger
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> import numpy as np
    >>> nw = Network()
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "pressure_difference": "bar",
    ...     "temperature": "degC"
    ... })
    >>> nw.iterinfo = False
    >>> so1 = Source("vapor source")
    >>> so2 = Source("air source")
    >>> cd = SectionedHeatExchanger("condenser")
    >>> si1 = Sink("water sink")
    >>> si2 = Sink("air sink")
    >>> c1 = Connection(so1, "out1", cd, "in1", label="1")
    >>> c2 = Connection(cd, "out1", si1, "in1", label="2")
    >>> c11 = Connection(so2, "out1", cd, "in2", label="11")
    >>> c12 = Connection(cd, "out2", si2, "in1", label="12")
    >>> nw.add_conns(c1, c2, c11, c12)

    To generate good guess values, first we run the simulation with fixed
    pressure on the water side. The water enters at superheated vapor state
    with 15 °C superheating and leaves it with 10 °C subcooling.

    >>> c1.set_attr(fluid={"Water": 1}, td_dew=15, m=1)
    >>> c2.set_attr(td_bubble=15, p=1)
    >>> c11.set_attr(fluid={"Air": 1}, T=15)
    >>> c12.set_attr(T=25, p=1)
    >>> cd.set_attr(dp1=0.0, dp2=0.0)
    >>> nw.solve("design")

    Now we can remove the pressure specifications on the air side and impose
    the minimum pinch instead, which will determine the actual water
    condensation pressure.

    >>> c2.set_attr(p=None)
    >>> cd.set_attr(td_pinch=5)
    >>> nw.solve("design")
    >>> round(c1.p.val, 3)
    0.056
    >>> round(c1.T.val, 1)
    50.0

    We can also see the temperature differences in all sections of the heat
    exchanger. Since the water vapor is cooled, condensed and then subcooled,
    while the air does not change phase, three sections will form:

    >>> Q_sections, T_steps_hot, T_steps_cold, Q_per_section, td_log_per_section = cd.calc_sections()
    >>> delta_T_between_sections = T_steps_hot - T_steps_cold
    >>> delta_T_list = [round(float(dT), 2) for dT in delta_T_between_sections]
    >>> delta_T_list[:6]
    [5.0, 16.8, 19.75, 19.6, 19.4, 19.2]

    We can see that the lowest delta T is the first one. This is the delta T
    between the hot side outlet and the cold side inlet, which can also be seen
    if we have a look at the network's results.

    >>> ();nw.print_results();() # doctest: +ELLIPSIS
    (...)

    If we change the subcooling degree at the water outlet, the condensation
    pressure and pinch will move.

    >>> c2.set_attr(td_bubble=5)
    >>> nw.solve("design")
    >>> round(c1.p.val, 3)
    0.042
    >>> Q_sections, T_steps_hot, T_steps_cold, Q_per_section, td_log_per_section = cd.calc_sections()
    >>> delta_T_between_sections = T_steps_hot - T_steps_cold
    >>> delta_T_list = [round(float(dT), 2) for dT in delta_T_between_sections]
    >>> delta_T_list[:6]
    [9.88, 14.8, 14.68, 14.48, 14.28, 14.08]

    Finally, in contrast to the baseclass :code:`HeatExchanger` `kA` value, the
    `UA` value takes into account the heat transfer per section and calculates
    the heat transfer coefficient as the sum of all sections, while the `kA`
    value only takes into account the inlet and outlet temperatures and the
    total heat transfer.

    >>> round(cd.kA.val)
    173307
    >>> round(cd.UA.val)
    273456

    It is also possible to apply a part-load modification to UA following the
    implementation of :cite:`cecchinato2010`. For this you have to specify
    :code:`UA_cecchinato` as offdesign parameter and along with it, values for

    - refrigerant side Reynolds exponent
    - secondary fluid side Reynolds exponent
    - secondary fluid to refrigerant area ratio
    - secondary fluid to refrigerant alpha (heat transfer coefficient) ratio
    - the refrigerant index (which side of the heat exchanger is passed by the
      refrigerant)

    >>> design_state = nw.save(as_dict=True)
    >>> cd.set_attr(
    ...     area_ratio=20,        # typical for a finned heat exchanger
    ...     alpha_ratio=1e-2,     # alpha for water side is higher
    ...     re_exp_r=0.8,
    ...     re_exp_sf=0.55,
    ...     refrigerant_index=0,  # water is refrigerant in this case
    ...     design=["td_pinch"],
    ...     offdesign=["UA_cecchinato"]
    ... )
    >>> nw.solve("offdesign", design_path=design_state)

    Without modifying any parameter, pinch and UA should be identical to
    design conditions.

    >>> round(cd.td_pinch.val, 2)
    5.0
    >>> round(cd.UA.val)
    273456

    With change in operating conditions, e.g. reduction of heat transfer
    we'd typically observe lower pinch, if the heat transfer reduces faster
    than the UA value does.

    >>> c1.set_attr(m=0.8)
    >>> nw.solve("offdesign", design_path=design_state)
    >>> round(cd.Q.val_SI / cd.Q.design, 2)
    0.8
    >>> round(cd.UA.val_SI / cd.UA.design, 2)
    0.88
    >>> round(cd.td_pinch.val, 2)
    4.3

    **Second Example**

    A transcritical gas cooler designed to cool CO2 from 160°C to approximately
    50°C while water is heated from 10°C to 60°C. The heat exchanger uses
    characteristic lines (`kA_char1` and `kA_char2`) to scale the heat transfer
    coefficient in offdesign operation as mass flow varies.

    This two-stage approach improves convergence:

    - **Stage 1**: Design with fixed pressures to establish initial guess
      values
    - **Stage 2**: Offdesign with characteristic line scaling for part-load
      analysis

    >>> from tespy.components import Source, Sink, SectionedHeatExchanger
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> from tespy.tools.characteristics import CharLine, load_default_char

    Set up the network with appropriate units:

    >>> nw = Network()
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "pressure_difference": "bar",
    ...     "temperature": "degC", "mass_flow": "kg/s"
    ... })
    >>> nw.iterinfo = False

    Create network components: two sources (CO2 and water inlets), the heat
    exchanger, and two sinks (outlets):

    >>> so_co2 = Source("CO2 source")
    >>> so_water = Source("Water source")
    >>> hx = SectionedHeatExchanger("transcritical gas cooler")
    >>> si_co2 = Sink("CO2 sink")
    >>> si_water = Sink("Water sink")

    Create connections with counter-current arrangement (CO2 on side 1, water
    on side 2):

    >>> c1 = Connection(so_co2, "out1", hx, "in1", label="CO2_in")
    >>> c2 = Connection(hx, "out1", si_co2, "in1", label="CO2_out")
    >>> c11 = Connection(so_water, "out1", hx, "in2", label="water_in")
    >>> c12 = Connection(hx, "out2", si_water, "in1", label="water_out")
    >>> nw.add_conns(c1, c2, c11, c12)

    **Stage 1: Design calculation with fixed pressures**

    First, we solve with fixed pressures on both sides to generate good initial
    guess values. This improves convergence for the complex transcritical
    cycle.

    Set CO2 inlet at 165 bar and 160°C (transcritical supercritical state):

    >>> c1.set_attr(
    ...     fluid={"CO2": 1},
    ...     p=165,
    ...     T=160,
    ...     m=3.5
    ... )

    Set water inlet at 5 bar and 10°C (cold inlet for cooling):

    >>> c11.set_attr(
    ...     fluid={"Water": 1},
    ...     p=5,
    ...     T=10
    ... )

    Specify water outlet temperature target at 60°C:

    >>> c12.set_attr(T=60)

    Configure heat exchanger with fixed pressures and initial pinch point:

    >>> hx.set_attr(
    ...     td_pinch=20,
    ...     pr1=1,
    ...     pr2=1,
    ...     num_sections=10
    ... )

    Solve the design point and save results:

    >>> nw.solve('design')
    >>> design_state = nw.save(as_dict=True)

    After design computation, the CO2 outlet state is:

    >>> round(c2.p.val, 1)
    165.0
    >>> round(c2.T.val, 1)
    30.0

    **Stage 2: Offdesign analysis with kA_char characteristic scaling**

    Now we activate characteristic line-based scaling. Load the default
    characteristic line for heat exchangers:

    >>> kA_char = load_default_char(
    ...     "HeatExchanger", "kA_char1", "DEFAULT", CharLine
    ... )

    Reconfigure heat exchanger to use characteristic lines for UA scaling in
    offdesign operation:

    >>> hx.set_attr(
    ...     kA_char1=kA_char,
    ...     kA_char2=kA_char,
    ...     design=['td_pinch'],
    ...     offdesign=['UA_char']
    ... )

    When offdesign is set to :code:`['UA_char']`, the solver automatically
    scales the UA value based on the characteristic curve during part-load
    operation. Verify offdesign setup by solving at design conditions. The
    design UA value is approximately 23.4 kW/K and pinch is 20.0 K:

    >>> nw.solve('offdesign', design_path=design_state)
    >>> round(hx.UA.val / 1e3, 2)
    23.42
    >>> round(hx.td_pinch.val, 1)
    20.0

    **Characteristic line scaling at part-load conditions**

    With variable mass flow, the UA value scales according to the
    characteristic curve. At 80 % mass flow, heat transfer reduces to roughly
    83 % while UA reduces by 9.5 % following the characteristic scaling:

    >>> c1.set_attr(m=2.8)
    >>> nw.solve('offdesign', design_path=design_state)
    >>> round(hx.Q.val_SI / hx.Q.design, 2)
    0.83
    >>> round(hx.UA.val_SI / hx.UA.design, 2)
    0.91

    The pinch point decreases to 15.3 K when heat transfer reduces faster than
    the characteristic-based UA scaling:

    >>> round(hx.td_pinch.val, 1)
    15.3

    The :code:`kA_char` parameter allows automatic part-load scaling of UA,
    following the same principle as the standard HeatExchanger component
    (:py:class:`tespy.components.heat_exchangers.base.HeatExchanger`). The
    difference to the :code:`UA_char` usage is that :code:`kA_char` uses a
    characteristic line lookup table to define the scaling relationship.
    :code:`UA_cecchinato` requires the specification of Reynolds number
    exponents, area ratio and alpha ratio of the involved fluids.
    """

    def get_parameters(self):
        params = super().get_parameters()
        params.update({
            'num_sections': dc_simple(
                val=50, dtype="int",
                description="number of sections of the heat exchanger"
            ),
            'UA': dc_cp(
                min_val=0, num_eq_sets=1,
                func=self.UA_func,
                dependents=self.UA_dependents,
                quantity="heat_transfer_coefficient",
                description="sum of UA values of all sections of heat exchanger"
            ),
             'UA_char': dc_gcc(
                    elements=['kA_char1', 'kA_char2'],
                    num_eq_sets=1,
                    func=self.UA_char_func,
                    dependents=self.UA_dependents,
                    description="equation for sectioned UA modification based on characteristic lines"
         ),
            'refrigerant_index': dc_simple(
                val=0, dtype="int",
                description="side on which the refrigerant is flowing (0: hot, 1:cold)"
            ),
            're_exp_r': dc_cp(
                description="Reynolds exponent for UA modification based on refrigerant side mass flow"
            ),
            're_exp_sf': dc_cp(
                description="Reynolds exponent for UA modification based on secondary fluid side mass flow"
            ),
            'alpha_ratio': dc_cp(
                quantity="ratio", min_val=0,
                description="secondary to refrigerant side convective heat transfer coefficient ratio"
            ),
            'area_ratio': dc_cp(
                quantity="ratio", min_val=0,
                description="secondary to refrigerant side heat transfer area ratio"
            ),
            'UA_cecchinato': dc_gcp(
                elements=['re_exp_r', 're_exp_sf', 'alpha_ratio', 'area_ratio'],
                num_eq_sets=1,
                func=self.UA_cecchinato_func,
                dependents=self.UA_dependents,
                description="equation for UA modification in offdesign"
            ),
            'td_pinch': dc_cp(
                min_val=0, num_eq_sets=1,
                func=self.td_pinch_func,
                dependents=self.td_pinch_dependents,
                quantity="temperature_difference",
                description="equation for minimum pinch"
            )
        })
        return params

    @staticmethod
    def _get_steps(num_steps=51):
        """Get the steps as fraction of enthalpy change for either side

        Parameters
        ----------
        c1 : tespy.connections.connection.Connection
            Inlet connection.

        c2 : tespy.connections.connection.Connection
            Outlet connection.

        Returns
        -------
        list
            Steps of enthalpy of the specified connections
        """
        return np.linspace(0, 1, num_steps)

    @staticmethod
    def _get_moving_steps(c1, c2):
        """Get the steps for enthalpy at the boundaries of phases during the
        change of enthalpy from one state to another

        Parameters
        ----------
        c1 : tespy.connections.connection.Connection
            Inlet connection.

        c2 : tespy.connections.connection.Connection
            Outlet connection.

        Returns
        -------
        list
            Steps of enthalpy of the specified connections
        """
        if c1.fluid.val != c2.fluid.val:
            msg = (
                "Both connections need to utilize the same fluid data: "
                f"{c1.fluid.val}, {c2.fluid.val}"
            )
            raise ValueError(msg)

        if c1.p.val_SI != c2.p.val_SI:
            msg = (
                "This method assumes equality of pressure for the inlet and "
                "the outlet connection. The pressure values provided are not "
                "equal, the results may be incorrect."
            )
        # change the order of connections to have c1 as the lower enthalpy
        # connection (enthalpy will be rising in the list)
        if c1.h.val_SI > c2.h.val_SI:
            c1, c2 = c2, c1

        h_at_steps = [0, 1]
        fluid = single_fluid(c1.fluid_data)
        # this should be generalized to "supports two-phase" because it does
        # not work with incompressibles
        is_pure_fluid = fluid is not None

        if is_pure_fluid:
            try:
                phase_h_low = phase_mix_ph(c1.p.val_SI, c1.h.val_SI, c1.fluid_data)
                phase_h_high = phase_mix_ph(c2.p.val_SI, c2.h.val_SI, c2.fluid_data)
            except NotImplementedError:
                return h_at_steps

            delta_h = c2.h.val_SI - c1.h.val_SI
            # we can round delta p here because we only need it in case it is
            # not zero, and then it will never be in the magnitude of subdigit
            # Pascal
            delta_p = round(c2.p.val_SI - c1.p.val_SI, 3)
            if phase_h_high == "g" and phase_h_low == "tp":
                if delta_p == 0:
                    h_sat_gas = h_mix_pQ(c1.p.val_SI, 1, c1.fluid_data)
                    x_gas = (h_sat_gas - c1.h.val_SI) / delta_h
                    h_at_steps = [0, x_gas, 1]
                else:
                    h_sat_gas = h_mix_pQ(c1.p.val_SI, 1, c1.fluid_data)
                    if np.isclose(h_sat_gas, c1.h.val_SI):
                        h_at_steps = [0, 1]
                    elif np.isclose(h_mix_pQ(c2.p.val_SI, 1, c2.fluid_data), c2.h.val_SI):
                        h_at_steps = [0, 1]
                    else:
                        x_gas = brentq(
                            identify_step_at_saturation,
                            0, 1,
                            args=(
                                c1.p.val_SI, c1.h.val_SI,
                                delta_p, delta_h,
                                1, c1.fluid_data#
                            )
                        )
                        h_at_steps = [0, x_gas, 1]

            elif phase_h_high == "g" and phase_h_low == "l":
                if delta_p == 0:
                    h_sat_gas = h_mix_pQ(c1.p.val_SI, 1, c1.fluid_data)
                    x_gas = (h_sat_gas - c1.h.val_SI) / delta_h
                    h_sat_liquid = h_mix_pQ(c1.p.val_SI, 0, c1.fluid_data)
                    x_liq = (h_sat_liquid - c1.h.val_SI) / delta_h
                    h_at_steps = [0, x_liq, x_gas, 1]
                else:
                    # c2 is the higher enthalpy, we have to check if it is
                    # at saturated gas
                    h_sat_gas = h_mix_pQ(c2.p.val_SI, 1, c2.fluid_data)
                    if np.isclose(h_sat_gas, c2.h.val_SI):
                        x_gas = 1
                    else:
                        x_gas = brentq(
                            identify_step_at_saturation,
                            0, 1,
                            args=(
                                c1.p.val_SI, c1.h.val_SI,
                                delta_p, delta_h,
                                1, c1.fluid_data
                            )
                        )
                    # c1 is the lower enthalpy, we have to check if it is
                    # at saturated liquid
                    h_sat_liquid = h_mix_pQ(c1.p.val_SI, 0, c1.fluid_data)
                    if np.isclose(h_sat_liquid, c1.h.val_SI):
                        x_liq = 0
                    else:
                        x_liq = brentq(
                            identify_step_at_saturation,
                            0, x_gas,
                            args=(
                                c1.p.val_SI, c1.h.val_SI,
                                delta_p, delta_h,
                                0, c1.fluid_data#
                            )
                        )
                    h_at_steps = [0, x_liq, x_gas, 1]

            elif phase_h_high == "tp" and phase_h_low == "l":
                if delta_p == 0:
                    h_sat_liquid = h_mix_pQ(c1.p.val_SI, 0, c1.fluid_data)
                    x_liq = (h_sat_liquid - c1.h.val_SI) / delta_h
                    h_at_steps = [0, x_liq, 1]
                else:
                    h_sat_liquid = h_mix_pQ(c1.p.val_SI, 0, c1.fluid_data)
                    if np.isclose(h_sat_liquid, c1.h.val_SI):
                        h_at_steps = [0, 1]
                    elif np.isclose(h_mix_pQ(c2.p.val_SI, 0, c2.fluid_data), c2.h.val_SI):
                        h_at_steps = [0, 1]
                    else:
                        x_liq = brentq(
                            identify_step_at_saturation,
                            0, 1,
                            args=(
                                c1.p.val_SI, c1.h.val_SI,
                                delta_p, delta_h,
                                0, c1.fluid_data#
                            )
                        )
                        h_at_steps = [0, x_liq, 1]

        return h_at_steps

    def _assign_steps(self):
        """Assign the sections of the heat exchanger

        Returns
        -------
        list
            List of cumulative sum of heat exchanged defining the heat exchanger
            sections.
        """
        num_steps = self.num_sections.val + 1
        steps = self._get_steps(num_steps)
        steps_hot = self._get_moving_steps(self.inl[0], self.outl[0])
        steps_cold = self._get_moving_steps(self.inl[1], self.outl[1])

        # unique throws out duplicates and sorts at the same time
        steps = np.unique(np.r_[steps, steps_hot, steps_cold])
        return steps

    def calc_UA(self, sections):
        """Calculate the sum of UA for all sections in the heat exchanger

        Returns
        -------
        float
            Sum of UA values of all heat exchanger sections.
        """
        _, _, _, Q_per_section, td_log_per_section = sections
        UA_sections = Q_per_section / td_log_per_section
        return sum(UA_sections)

    @staticmethod
    def _min_td(sections):
        """Return the minimum hot-minus-cold temperature difference."""
        _, T_hot, T_cold, _, _ = sections
        return np.min(T_hot - T_cold)

    def UA_func(self, **kwargs):
        r"""
        Residual method for fixed heat transfer coefficient UA.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = UA - \sum UA_{i}
        """
        sections = self.calc_sections(False)
        min_td = self._min_td(sections)
        if min_td <= 0.0:
            # Invalid pinch: _calc_td_log_per_section clips negative td to
            # 1e-3 K, making UA_calc >> UA_target (large negative first term).
            # Adding min_td injects a temperature-based gradient that is not
            # proportional to the energy balance row, preventing linear
            # dependency in the Jacobian.  The combined residual is never zero
            # while min_td <= 0, so the solver cannot converge to the false
            # fixed point at the thermodynamic limit.
            return self.UA.val_SI - self.calc_UA(sections) + min_td
        return self.UA.val_SI - self.calc_UA(sections)

    def UA_char_func(self):
        r"""
        Calculate offdesign UA from characteristic lines analogous to standard
        heat exchanger kA_char, but for the sectioned heat exchanger.

        Returns
        -------
        float
            Residual value of equation:

            .. math::

                0 = UA_\text{design} * f_\text{UA} - \sum\left(UA_{i}\right)

        """
        p1 = self.kA_char1.param
        p2 = self.kA_char2.param

        f1 = self.get_char_expr(p1, **self.kA_char1.char_params)
        f2 = self.get_char_expr(p2, **self.kA_char2.char_params)

        fUA1 = self.kA_char1.char_func.evaluate(f1)
        fUA2 = self.kA_char2.char_func.evaluate(f2)

        fUA = 2 / (1 / fUA1 + 1 / fUA2)

        sections = self.calc_sections(False)
        min_td = self._min_td(sections)
        if min_td <= 0:
            return self.UA.design * fUA - self.calc_UA(sections) + min_td
        return self.UA.design * fUA - self.calc_UA(sections)

    def UA_cecchinato_func(self):
        r"""
        Method to calculate heat transfer via UA design with modification
        for part load according to :cite:`cecchinato2010`. UA is determined
        over the UA values of the sections of the heat exchanger.

        .. note::

            You need to specify a couple of parameters to use this method. The
            values depend on the context, as they define relations between the
            refrigerant and the secondary fluid. For an evaporator the
            refrigerant is the cold side, for a condenser it is the hot side.
            You can check the linked publication for reference values.

        - alpha_ratio: Ratio of secondary fluid to refrigerant heat transfer
          coefficient
        - area_ratio: Ratio of secondary fluid to refrigerant area
        - re_exp_sf: Reynolds exponent for the secondary fluid mass flow
        - re_exp_r: Reynolds exponent for the refrigerant mass flow

        The modification factor for UA is calculated as follows

        .. math::

            f_\text{UA}=\frac{
                1 + \frac{\alpha_\text{sf}}{\alpha_\text{r}}
                \cdot\frac{A_\text{sf}}{A_\text{r}}
            }{
                \frac{\dot m_\text{sf}}{\dot m_\text{sf,ref}}^{-Re_\text{sf}} +
                \frac{\alpha_\text{sf}}{\alpha_\text{r}}
                \cdot\frac{A_\text{sf}}{A_\text{r}}
                \cdot\frac{\dot m_\text{r}}{\dot m_\text{r,ref}}^{-Re_\text{r}}
            }

        Returns
        -------
        float
            residual value of equation

            .. math::

                0 = UA_\text{ref} \cdot f_\text{UA} - \sum UA_\text{i}
        """
        alpha_ratio = self.alpha_ratio.val_SI
        area_ratio = self.area_ratio.val_SI
        re_exp_sf = self.re_exp_sf.val_SI
        re_exp_r = self.re_exp_r.val_SI
        refrigerant_index = self.refrigerant_index.val

        if refrigerant_index == 0:
            secondary_index = 1
        else:
            secondary_index = 0

        m_r = self.inl[refrigerant_index].m
        m_ratio_r = max(
            m_r.val_SI / self._conn_design(self.inl[refrigerant_index], 'm'),
            1e-6
        )
        m_sf = self.inl[secondary_index].m
        m_ratio_sf = max(
            m_sf.val_SI / self._conn_design(self.inl[secondary_index], 'm'),
            1e-6
        )

        fUA = (
            (1 + alpha_ratio * area_ratio)
            / (
                m_ratio_sf ** -re_exp_sf
                + alpha_ratio * area_ratio * m_ratio_r ** -re_exp_r
            )
        )
        sections = self.calc_sections(False)
        min_td = self._min_td(sections)
        if min_td <= 0:
            return self.UA.design * fUA - self.calc_UA(sections) + min_td
        return self.UA.design * fUA - self.calc_UA(sections)

    def UA_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].p,
            self.outl[0].h,
            self.inl[1].m,
            self.inl[1].p,
            self.inl[1].h,
            self.outl[1].p,
            self.outl[1].h
        ]

    def calc_td_pinch(self, T_steps_hot, T_steps_cold):
        """Calculate the pinch point temperature difference

        Returns
        -------
        float
            Value of the pinch point temperature difference
        """
        return min(T_steps_hot - T_steps_cold)

    def td_pinch_func(self):
        r"""
        Equation for minimal pinch temperature difference of sections.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = td_\text{pinch} - min(td_\text{i})
        """
        steps = self._assign_steps()
        T_hot, T_cold = self._get_T_at_steps(steps)
        return self.td_pinch.val_SI - self.calc_td_pinch(T_hot, T_cold)

    def td_pinch_dependents(self):
        return [
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].p,
            self.outl[0].h,
            self.inl[1].p,
            self.inl[1].h,
            self.outl[1].p,
            self.outl[1].h
        ]

    def calc_parameters(self):
        super().calc_parameters()

        sections = self.calc_sections()
        self.UA.val_SI = self.calc_UA(sections)
        self.td_pinch.val_SI = self.calc_td_pinch(sections[1], sections[2])


def identify_step_at_saturation(x, p_in, h_in, delta_p, delta_h, Q, fluid_data):
    r"""Method to identify the step corresponding to a saturation line assuming
    the change of pressure delta p is linear to the change of enthalpy delta h.

    .. math::

        \Delta h \cdot \left(p_\text{sat} - p_\text{in}\right) =
        \Delta p \cdot \left(h_\left[p_\text{sat}, Q\right] - h_\text{in}
        \right)


    Parameters
    ----------
    x : float
        Step to solve for
    p_in : float
        pressure at inlet
    h_in : float
        enthalpy at inlet
    delta_p : float
        overall pressure difference
    delta_h : float
        overall enthalpy difference
    Q : float
        vapor mass fraction (0 for bubble line, 1 for dew line)
    fluid_data : dict
        fluid_data dictionary

    Returns
    -------
    float
        residual of equation
    """
    p_sat = p_in + delta_p * x
    h_sat = h_mix_pQ(p_sat, Q, fluid_data)
    return delta_h * (p_sat - p_in) - (h_sat - h_in) * delta_p
