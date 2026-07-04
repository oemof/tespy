# -*- coding: utf-8

"""Module of class MovingBoundaryHeatExchanger.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/heat_exchangers/movingboundary.py

SPDX-License-Identifier: MIT
"""
import numpy as np

from tespy.components.component import component_registry
from tespy.components.heat_exchangers.sectioned import SectionedHeatExchanger


@component_registry
class MovingBoundaryHeatExchanger(SectionedHeatExchanger):
    r"""
    Class for counter flow heat exchanger with UA sections.

    The heat exchanger is internally discretized into multiple sections, which
    are defined by phase changes. The component assumes, that a pressure drop
    is linear to the change in enthalpy, meaning the phase boundary
    identification is done iteratively. In principle the implementations
    follows :cite:`bell2015`.

    .. image:: /api/_images/components/HeatExchanger.svg
       :alt: flowsheet of the movingboundaryheatexchanger
       :align: center
       :class: only-light

    .. image:: /api/_images/components/HeatExchanger_darkmode.svg
       :alt: flowsheet of the movingboundaryheatexchanger
       :align: center
       :class: only-dark

    Ports
    -----

    - Fluid inlets: in1, in2
    - Fluid outlets: out1, out2

    Mandatory Equations
    -------------------

    - mass flow equality constraint(s): :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
    - fluid composition equality constraint(s): :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
    - hot side to cold side heat transfer equation: :py:meth:`energy_balance_func <tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_func>`

    Parameters
    ----------

    alpha1_g : float, dict
        Hot-side heat transfer coefficient in superheated zone. Quantity:
        :code:`heat_transfer_coefficient_per_area`.

    alpha1_l : float, dict
        Hot-side heat transfer coefficient in subcooled zone. Quantity:
        :code:`heat_transfer_coefficient_per_area`.

    alpha1_sc : float, dict
        Hot-side heat transfer coefficient in supercritical zone. Quantity:
        :code:`heat_transfer_coefficient_per_area`.

    alpha1_tp : float, dict
        Hot-side heat transfer coefficient in two-phase zone. Quantity:
        :code:`heat_transfer_coefficient_per_area`.

    alpha2_g : float, dict
        Cold-side heat transfer coefficient in superheated zone. Quantity:
        :code:`heat_transfer_coefficient_per_area`.

    alpha2_l : float, dict
        Cold-side heat transfer coefficient in subcooled zone. Quantity:
        :code:`heat_transfer_coefficient_per_area`.

    alpha2_sc : float, dict
        Cold-side heat transfer coefficient in supercritical zone. Quantity:
        :code:`heat_transfer_coefficient_per_area`.

    alpha2_tp : float, dict
        Cold-side heat transfer coefficient in two-phase zone. Quantity:
        :code:`heat_transfer_coefficient_per_area`.

    alpha_ratio : float, dict
        Secondary to refrigerant side convective heat transfer coefficient
        ratio. Quantity: :code:`ratio`.

    area_hot : float, dict
        Hot-side heat exchange area. Quantity: :code:`area`.

    area_ratio : float, dict
        Heat transfer area ratio; previously defined as secondary to refrigerant
        side ratio, will be defined as hot to cold side ratio in a future
        version. Quantity: :code:`ratio`.

    area_zones : GroupedComponentProperties
        Bell (2015) area-based heat exchanger constraint. All elements must be
        set for the group to activate. For phases that do not occur in your
        application set the corresponding alpha to any value - it only needs to
        be set, as it will not be used. Elements: :code:`area_hot`,
        :code:`area_ratio`, :code:`alpha1_l`, :code:`alpha1_tp`,
        :code:`alpha1_g`, :code:`alpha1_sc`, :code:`alpha2_l`,
        :code:`alpha2_tp`, :code:`alpha2_g`, :code:`alpha2_sc`, :code:`R_cond`.
        Equation: :py:meth:`area_zones_func <tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger.area_zones_func>`.

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
        Deprecated, use :code:`UA` instead. Quantity:
        :code:`heat_transfer_coefficient`.

    kA_char : GroupedComponentCharacteristics
        Deprecated, use :code:`UA_char` instead. Elements: :code:`kA_char1`,
        :code:`kA_char2`.

    kA_char1 : tespy.tools.characteristics.CharLine, dict
        Deprecated, use :code:`UA_char1` instead.

    kA_char2 : tespy.tools.characteristics.CharLine, dict
        Deprecated, use :code:`UA_char2` instead.

    label : str
        The label of the component.

    lmtd : float, dict
        Effective logarithmic mean temperature difference :code:`Q/UA`.
        Quantity: :code:`temperature_difference`.

    lmtd_per_section : numpy.ndarray
        Logarithmic mean temperature difference in each section. Quantity:
        :code:`temperature_difference`. Result only - populated by the network
        after each solve.

    local_design : bool
        Treat this component in design mode in an offdesign calculation.

    local_offdesign : bool
        Treat this component in offdesign mode in a design calculation.

    offdesign : list
        List containing offdesign parameters (stated as String).

    phase_cold_per_section : numpy.ndarray
        Phase index per section on cold side (0=liquid, 1=two-phase, 2=gas,
        3=supercritical). Result only - populated by the network after each
        solve.

    phase_hot_per_section : numpy.ndarray
        Phase index per section on hot side (0=liquid, 1=two-phase, 2=gas,
        3=supercritical). Result only - populated by the network after each
        solve.

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

    Q_per_section : numpy.ndarray
        Heat transferred from hot to cold side in each section. Quantity:
        :code:`heat`. Result only - populated by the network after each solve.

    Q_sections : numpy.ndarray
        Cumulative heat transferred from hot to cold side up to each section
        boundary. Quantity: :code:`heat`. Result only - populated by the network
        after each solve.

    R_cond : float, dict
        Wall conduction thermal resistance. Quantity:
        :code:`thermal_resistance`.

    re_exp_cold : float, dict
        Reynolds exponent for UA modification based on cold side mass flow.

    re_exp_hot : float, dict
        Reynolds exponent for UA modification based on hot side mass flow.

    re_exp_r : float, dict
        Deprecated - Reynolds exponent for refrigerant side mass flow; use
        :code:`re_exp_hot` or :code:`re_exp_cold` depending on which side the
        refrigerant flows on.

    re_exp_sf : float, dict
        Deprecated - Reynolds exponent for secondary fluid side mass flow; use
        :code:`re_exp_hot` or :code:`re_exp_cold` depending on which side the
        secondary fluid flows on.

    refrigerant_index : int
        Deprecated - side on which the refrigerant is flowing (0: hot, 1:cold).

    T_cold_sections : numpy.ndarray
        Cold side temperature at each section boundary. Quantity:
        :code:`temperature`. Result only - populated by the network after each
        solve.

    T_hot_sections : numpy.ndarray
        Hot side temperature at each section boundary. Quantity:
        :code:`temperature`. Result only - populated by the network after each
        solve.

    td_log : float, dict
        Deprecated, use :code:`lmtd` instead. Quantity:
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
        Deprecated - equation for UA modification in offdesign using
        refrigerant/secondary-fluid Reynolds exponents; use
        :code:`UA_cecchinato_hc` with :code:`re_exp_hot` and :code:`re_exp_cold`
        instead. Elements: :code:`re_exp_r`, :code:`re_exp_sf`,
        :code:`alpha_ratio`, :code:`area_ratio`.
        Equation: :py:meth:`UA_cecchinato_legacy_func <tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger.UA_cecchinato_legacy_func>`.

    UA_cecchinato_hc : GroupedComponentProperties
        Temporary name - equation for UA modification in offdesign using
        explicit hot/cold Reynolds exponents; in the next major version
        :code:`UA_cecchinato` will adopt the hot/cold convention of
        :code:`UA_cecchinato_hc`. Elements: :code:`re_exp_hot`,
        :code:`re_exp_cold`, :code:`alpha_ratio`, :code:`area_ratio`.
        Equation: :py:meth:`UA_cecchinato_func <tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger.UA_cecchinato_func>`.

    UA_char : GroupedComponentCharacteristics
        Equation for sectioned UA modification based on characteristic lines.
        Elements: :code:`UA_char1`, :code:`UA_char2`.
        Equation: :py:meth:`UA_char_func <tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger.UA_char_func>`.

    UA_char1 : tespy.tools.characteristics.CharLine, dict
        Hot side UA modification lookup table for offdesign.

    UA_char2 : tespy.tools.characteristics.CharLine, dict
        Cold side UA modification lookup table for offdesign.

    zeta1 : float, dict
        Deprecated, use :code:`zeta1_d4` instead.

    zeta1_d4 : float, dict
        Hot side geometry-independent friction coefficient zeta/D^4 for pressure
        loss calculation.
        Equation: :py:meth:`zeta_d4_func <tespy.components.component.Component.zeta_d4_func>`.

    zeta2 : float, dict
        Deprecated, use :code:`zeta2_d4` instead.

    zeta2_d4 : float, dict
        Cold side geometry-independent friction coefficient zeta/D^4 for
        pressure loss calculation.
        Equation: :py:meth:`zeta_d4_func <tespy.components.component.Component.zeta_d4_func>`.

    Notes
    -----

    .. note::

        The equations only apply to counter-current heat exchangers.

    Example
    -------
    Water vapor should be cooled down, condensed and then further subcooled. For
    this air is heated up from 15 °C to 25 °C.

    >>> from tespy.components import Source, Sink, MovingBoundaryHeatExchanger
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
    >>> cd = MovingBoundaryHeatExchanger("condenser")
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

    >>> c1.set_attr(fluid={"Water": 1}, p=1, td_dew=15, m=1)
    >>> c2.set_attr(td_bubble=15)
    >>> c11.set_attr(fluid={"Air": 1}, p=1, T=15)
    >>> c12.set_attr(T=25)
    >>> cd.set_attr(pr1=1, pr2=1)
    >>> nw.solve("design")

    Now we can remove the pressure specifications on the air side and impose
    the minimum pinch instead, which will determine the actual water
    condensation pressure.

    >>> c1.set_attr(p=None)
    >>> cd.set_attr(td_pinch=5)
    >>> nw.solve("design")
    >>> round(c1.p.val, 3)
    0.056
    >>> round(c1.T.val, 1)
    50.0

    After solving, section data is available directly via the component
    attributes :code:`T_hot_sections`, :code:`T_cold_sections`,
    :code:`Q_sections`, :code:`Q_per_section` and :code:`lmtd_per_section`.
    Since the water vapor is cooled, condensed and then subcooled while the
    air does not change phase, three sections will form:

    >>> delta_T_between_sections = cd.T_hot_sections.val_SI - cd.T_cold_sections.val_SI
    >>> delta_T_between_sections.round(2).tolist()
    [5.0, 19.75, 10.11, 25.0]

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
    >>> delta_T_between_sections = cd.T_hot_sections.val_SI - cd.T_cold_sections.val_SI
    >>> delta_T_between_sections.round(2).tolist()
    [9.88, 14.8, 5.0, 19.88]

    Finally, in contrast to the baseclass :code:`HeatExchanger` `kA` value, the
    `UA` value takes into account the heat transfer per section and calculates
    the heat transfer coefficient as the sum of all sections, while the `kA`
    value only takes into account the inlet and outlet temperatures and the
    total heat transfer.

    >>> round(cd.kA.val)
    173307
    >>> round(cd.UA.val)
    273449

    It is also possible to apply a partload modification to UA following the
    implementation of :cite:`cecchinato2010`. For this you have to specify
    :code:`UA_cecchinato` as offdesign parameter and along with it, values for

    - hot side Reynolds exponent (:code:`re_exp_hot`)
    - cold side Reynolds exponent (:code:`re_exp_cold`)
    - hot to cold side area ratio (:code:`area_ratio`)
    - hot to cold side alpha (heat transfer coefficient) ratio (:code:`alpha_ratio`)

    >>> design_state = nw.save(as_dict=True)
    >>> cd.set_attr(
    ...     area_ratio=20,        # typical for a finned heat exchanger
    ...     alpha_ratio=1e-2,     # alpha for water side is higher
    ...     re_exp_hot=0.8,
    ...     re_exp_cold=0.55,
    ...     design=["td_pinch"],
    ...     offdesign=["UA_cecchinato_hc"]
    ... )
    >>> nw.solve("offdesign", design_path=design_state)

    Without modifying any parameter, pinch and UA should be identical to
    design conditions.

    >>> round(cd.td_pinch.val, 2)
    5.0
    >>> round(cd.UA.val)
    273449

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
    """
    def get_parameters(self):
        params = super().get_parameters()
        del params["num_sections"]
        return params

    def _assign_steps(self, steps_hot=None, steps_cold=None):
        """Assign the sections of the heat exchanger

        Parameters
        ----------
        steps_hot : list, optional
            Pre-computed phase-boundary steps for the hot side. Computed from
            :py:meth:`_get_moving_steps` when not provided.
        steps_cold : list, optional
            Pre-computed phase-boundary steps for the cold side. Computed from
            :py:meth:`_get_moving_steps` when not provided.

        Returns
        -------
        list
            List of cumulative sum of heat exchanged defining the heat exchanger
            sections.
        """
        if steps_hot is None:
            steps_hot, _ = self._get_moving_steps(self.inl[0], self.outl[0])
        if steps_cold is None:
            steps_cold, _ = self._get_moving_steps(self.inl[1], self.outl[1])
        return np.unique(np.r_[steps_hot, steps_cold])
