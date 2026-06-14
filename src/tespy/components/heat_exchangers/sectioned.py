# -*- coding: utf-8

"""Module of class SectionedHeatExchanger.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/heat_exchangers/sectioned.py

SPDX-License-Identifier: MIT
"""
import warnings

import numpy as np
from scipy.optimize import brentq

from tespy.components.component import component_registry
from tespy.components.heat_exchangers.base import HeatExchanger
from tespy.tools.data_containers import ComponentArrayProperties as dc_cap
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentCharacteristics as dc_gcc
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
from tespy.tools.fluid_properties import h_mix_pQ
from tespy.tools.fluid_properties import phase_mix_ph
from tespy.tools.fluid_properties import single_fluid

_PHASE_TO_INT = {"l": 0, "tp": 1, "g": 2, "sc": 3}


@component_registry
class SectionedHeatExchanger(HeatExchanger):
    r"""
    Class for counter flow heat exchanger with UA sections.

    The heat exchanger is internally discretized into 51 sections of equal heat
    transfer. The number of section can be adjusted by the user. It is based a
    combination of the moving boundary approach by :cite:`bell2015` and
    discretization in :cite:`Quoilin2020`.

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

    num_sections : int
        Number of sections of the heat exchanger.

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

    After solving, section data is available directly via the component
    attributes :code:`T_hot_sections`, :code:`T_cold_sections`,
    :code:`Q_sections`, :code:`Q_per_section`, :code:`lmtd_per_section`,
    :code:`phase_hot_per_section` and :code:`phase_cold_per_section`.
    The phase attributes hold integer arrays with one entry per section using
    the mapping 0=liquid, 1=two-phase, 2=gas, 3=supercritical.
    Since the water vapor is cooled, condensed and then subcooled while the
    air does not change phase, three sections will form:

    >>> delta_T_between_sections = cd.T_hot_sections.val_SI - cd.T_cold_sections.val_SI
    >>> delta_T_between_sections[:6].round(2).tolist()
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
    >>> delta_T_between_sections = cd.T_hot_sections.val_SI - cd.T_cold_sections.val_SI
    >>> delta_T_between_sections[:6].round(2).tolist()
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

    >>> nw.solve("design")
    >>> design_state = nw.save(as_dict=True)

    After design computation, the CO2 outlet state is:

    >>> round(c2.p.val, 1)
    165.0
    >>> round(c2.T.val, 1)
    30.0

    **Stage 2: Offdesign analysis with UA_char characteristic scaling**

    Now we activate characteristic line-based scaling. Load the default
    characteristic line for heat exchangers:

    >>> UA_char = load_default_char(
    ...     "HeatExchanger", "UA_char1", "DEFAULT", CharLine
    ... )

    Reconfigure heat exchanger to use characteristic lines for UA scaling in
    offdesign operation:

    >>> hx.set_attr(
    ...     UA_char1=UA_char,
    ...     UA_char2=UA_char,
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

    The :code:`UA_char` parameter allows automatic part-load scaling of UA,
    following the same principle as the standard HeatExchanger component
    (:py:class:`tespy.components.heat_exchangers.base.HeatExchanger`).
    :code:`UA_cecchinato` requires the specification of Reynolds number
    exponents, area ratio and alpha ratio of the involved fluids.
    """

    _parameter_aliases = {
        'kA_char': 'UA_char',
        'kA_char1': 'UA_char1',
        'kA_char2': 'UA_char2',
        'zeta1': 'zeta1_d4',
        'zeta2': 'zeta2_d4',
    }

    def set_attr(self, **kwargs):
        for old, msg in [
            (
                'refrigerant_index',
                f"The parameter 'refrigerant_index' of component {self.label!r} is "
                "deprecated. Use 'UA_cecchinato_hc' with 're_exp_hot' and 're_exp_cold' "
                "instead. In the next major version 'UA_cecchinato' will adopt the "
                "hot/cold convention of 'UA_cecchinato_hc'."
            ),
            (
                're_exp_r',
                f"The parameter 're_exp_r' of component {self.label!r} is deprecated. "
                "Use 'UA_cecchinato_hc' with 're_exp_hot' or 're_exp_cold' instead, "
                "depending on which side the refrigerant flows on. In the next major "
                "version 'UA_cecchinato' will adopt the hot/cold convention of "
                "'UA_cecchinato_hc'."
            ),
            (
                're_exp_sf',
                f"The parameter 're_exp_sf' of component {self.label!r} is deprecated. "
                "Use 'UA_cecchinato_hc' with 're_exp_hot' or 're_exp_cold' instead, "
                "depending on which side the secondary fluid flows on. In the next "
                "major version 'UA_cecchinato' will adopt the hot/cold convention of "
                "'UA_cecchinato_hc'."
            ),
        ]:
            if old in kwargs:
                warnings.warn(msg, FutureWarning, stacklevel=2)
        super().set_attr(**kwargs)

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
                description="sum of UA values of all sections of heat exchanger",
                calc=self._calc_UA_from_sections,
                calc_deps=[]
            ),
            'UA_char': dc_gcc(
                elements=['UA_char1', 'UA_char2'],
                num_eq_sets=1,
                func=self.UA_char_func,
                dependents=self.UA_dependents,
                description="equation for sectioned UA modification based on characteristic lines"
            ),
            'refrigerant_index': dc_simple(
                val=0, dtype="int",
                description="deprecated - side on which the refrigerant is flowing (0: hot, 1:cold)"
            ),
            're_exp_hot': dc_cp(
                description="Reynolds exponent for UA modification based on hot side mass flow"
            ),
            're_exp_cold': dc_cp(
                description="Reynolds exponent for UA modification based on cold side mass flow"
            ),
            're_exp_r': dc_cp(
                description="deprecated - Reynolds exponent for refrigerant side mass flow; "
                "use :code:`re_exp_hot` or :code:`re_exp_cold` depending on which side "
                "the refrigerant flows on"
            ),
            're_exp_sf': dc_cp(
                description="deprecated - Reynolds exponent for secondary fluid side mass flow; "
                "use :code:`re_exp_hot` or :code:`re_exp_cold` depending on which side "
                "the secondary fluid flows on"
            ),
            'alpha_ratio': dc_cp(
                quantity="ratio", min_val=0,
                description="secondary to refrigerant side convective heat transfer coefficient ratio"
            ),
            'area_ratio': dc_cp(
                quantity="ratio", min_val=0,
                description="heat transfer area ratio; previously defined as secondary to "
                "refrigerant side ratio, will be defined as hot to cold side ratio in a "
                "future version"
            ),
            'UA_cecchinato': dc_gcp(
                elements=['re_exp_r', 're_exp_sf', 'alpha_ratio', 'area_ratio'],
                num_eq_sets=1,
                func=self.UA_cecchinato_legacy_func,
                dependents=self.UA_dependents,
                description=(
                    "deprecated - equation for UA modification in offdesign using "
                    "refrigerant/secondary-fluid Reynolds exponents; use "
                    ":code:`UA_cecchinato_hc` with :code:`re_exp_hot` and "
                    ":code:`re_exp_cold` instead"
                )
            ),
            'UA_cecchinato_hc': dc_gcp(
                elements=['re_exp_hot', 're_exp_cold', 'alpha_ratio', 'area_ratio'],
                num_eq_sets=1,
                func=self.UA_cecchinato_func,
                dependents=self.UA_dependents,
                description=(
                    "temporary name - equation for UA modification in offdesign using "
                    "explicit hot/cold Reynolds exponents; in the next major version "
                    ":code:`UA_cecchinato` will adopt the hot/cold convention of "
                    ":code:`UA_cecchinato_hc`"
                )
            ),
            'td_pinch': dc_cp(
                min_val=0, num_eq_sets=1,
                func=self.td_pinch_func,
                dependents=self.td_pinch_dependents,
                quantity="temperature_difference",
                description="equation for minimum pinch",
                calc=self._calc_td_pinch,
                calc_deps=[]
            ),
            'alpha1_l': dc_cp(
                min_val=0, quantity="heat_transfer_coefficient_per_area",
                description="hot-side heat transfer coefficient in subcooled zone"
            ),
            'alpha1_tp': dc_cp(
                min_val=0, quantity="heat_transfer_coefficient_per_area",
                description="hot-side heat transfer coefficient in two-phase zone"
            ),
            'alpha1_g': dc_cp(
                min_val=0, quantity="heat_transfer_coefficient_per_area",
                description="hot-side heat transfer coefficient in superheated zone"
            ),
            'alpha1_sc': dc_cp(
                min_val=0, quantity="heat_transfer_coefficient_per_area",
                description="hot-side heat transfer coefficient in supercritical zone"
            ),
            'alpha2_l': dc_cp(
                min_val=0, quantity="heat_transfer_coefficient_per_area",
                description="cold-side heat transfer coefficient in subcooled zone"
            ),
            'alpha2_tp': dc_cp(
                min_val=0, quantity="heat_transfer_coefficient_per_area",
                description="cold-side heat transfer coefficient in two-phase zone"
            ),
            'alpha2_g': dc_cp(
                min_val=0, quantity="heat_transfer_coefficient_per_area",
                description="cold-side heat transfer coefficient in superheated zone"
            ),
            'alpha2_sc': dc_cp(
                min_val=0, quantity="heat_transfer_coefficient_per_area",
                description="cold-side heat transfer coefficient in supercritical zone"
            ),
            'phase_hot_per_section': dc_cap(
                quantity=None,
                description="phase index per section on hot side (0=liquid, 1=two-phase, 2=gas, 3=supercritical)"
            ),
            'phase_cold_per_section': dc_cap(
                quantity=None,
                description="phase index per section on cold side (0=liquid, 1=two-phase, 2=gas, 3=supercritical)"
            ),
            'R_cond': dc_cp(
                min_val=0, quantity="thermal_resistance",
                description="wall conduction thermal resistance"
            ),
            'area_hot': dc_cp(
                min_val=0, quantity="area",
                description="hot-side heat exchange area"
            ),
            'area_zones': dc_gcp(
                elements=[
                    'area_hot', 'area_ratio',
                    'alpha1_l', 'alpha1_tp', 'alpha1_g', 'alpha1_sc',
                    'alpha2_l', 'alpha2_tp', 'alpha2_g', 'alpha2_sc',
                    'R_cond',
                ],
                num_eq_sets=1,
                func=self.area_zones_func,
                dependents=self.area_zones_dependents,
                description=(
                    "Bell (2015) area-based heat exchanger constraint. All "
                    "elements must be set for the group to activate. For "
                    "phases that do not occur in your application set the "
                    "corresponding alpha to any value - it only needs to be "
                    "set, as it will not be used."
                )
            )
        })
        return params

    def _store_sections(self):
        super()._store_sections()
        steps1, zp1 = self._get_moving_steps(self.inl[0], self.outl[0])
        steps2, zp2 = self._get_moving_steps(self.inl[1], self.outl[1])
        steps_all = self._assign_steps(steps1, steps2)
        self.phase_hot_per_section.val_SI = np.array(
            self._section_phases(steps_all, np.array(steps1), zp1)
        )
        self.phase_cold_per_section.val_SI = np.array(
            self._section_phases(steps_all, np.array(steps2), zp2)
        )

    @staticmethod
    def _get_steps(num_steps=51):
        """Return :code:`num_steps` evenly-spaced fractions in [0, 1].

        Parameters
        ----------
        num_steps : int
            Number of points (equals number of sections + 1).

        Returns
        -------
        numpy.ndarray
            Uniformly spaced step fractions from 0 to 1.
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
        tuple
            Steps of enthalpy of the specified connections and a list of phase
            indices (0=L, 1=TP, 2=G, 3=SC) for each zone between consecutive steps.
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
        zone_phases = [2]
        fluid = single_fluid(c1.fluid_data)
        # this should be generalized to "supports two-phase" because it does
        # not work with incompressibles
        is_pure_fluid = fluid is not None

        if is_pure_fluid:
            phase_h_low = phase_mix_ph(c1.p.val_SI, c1.h.val_SI, c1.fluid_data)
            phase_h_high = phase_mix_ph(c2.p.val_SI, c2.h.val_SI, c2.fluid_data)

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
                    zone_phases = [1, 2]
                else:
                    h_sat_gas = h_mix_pQ(c1.p.val_SI, 1, c1.fluid_data)
                    if np.isclose(h_sat_gas, c1.h.val_SI):
                        h_at_steps = [0, 1]
                        zone_phases = [2]
                    elif np.isclose(h_mix_pQ(c2.p.val_SI, 1, c2.fluid_data), c2.h.val_SI):
                        h_at_steps = [0, 1]
                        zone_phases = [1]
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
                        zone_phases = [1, 2]

            elif phase_h_high == "g" and phase_h_low == "l":
                if delta_p == 0:
                    h_sat_gas = h_mix_pQ(c1.p.val_SI, 1, c1.fluid_data)
                    x_gas = (h_sat_gas - c1.h.val_SI) / delta_h
                    h_sat_liquid = h_mix_pQ(c1.p.val_SI, 0, c1.fluid_data)
                    x_liq = (h_sat_liquid - c1.h.val_SI) / delta_h
                    h_at_steps = [0, x_liq, x_gas, 1]
                    zone_phases = [0, 1, 2]
                else:
                    # c2 is the higher enthalpy, we have to check if it is
                    # at saturated gas
                    h_sat_gas = h_mix_pQ(c2.p.val_SI, 1, c2.fluid_data)
                    x_gas_is_one = np.isclose(h_sat_gas, c2.h.val_SI)
                    if x_gas_is_one:
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
                    x_liq_is_zero = np.isclose(h_sat_liquid, c1.h.val_SI)
                    if x_liq_is_zero:
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
                    if x_liq_is_zero and x_gas_is_one:
                        zone_phases = [1]
                    elif x_liq_is_zero:
                        zone_phases = [1, 2]
                    elif x_gas_is_one:
                        zone_phases = [0, 1]
                    else:
                        zone_phases = [0, 1, 2]

            elif phase_h_high == "tp" and phase_h_low == "l":
                if delta_p == 0:
                    h_sat_liquid = h_mix_pQ(c1.p.val_SI, 0, c1.fluid_data)
                    x_liq = (h_sat_liquid - c1.h.val_SI) / delta_h
                    h_at_steps = [0, x_liq, 1]
                    zone_phases = [0, 1]
                else:
                    h_sat_liquid = h_mix_pQ(c1.p.val_SI, 0, c1.fluid_data)
                    if np.isclose(h_sat_liquid, c1.h.val_SI):
                        h_at_steps = [0, 1]
                        zone_phases = [1]
                    elif np.isclose(h_mix_pQ(c2.p.val_SI, 0, c2.fluid_data), c2.h.val_SI):
                        h_at_steps = [0, 1]
                        zone_phases = [0]
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
                        zone_phases = [0, 1]

            elif "sc" in (phase_h_high, phase_h_low):
                # CoolProp < 8 raises ValueError for h_pT(p, T_crit) when
                # p > p_crit, so the sc/l boundary cannot be located reliably.
                # CoolProp 8 fixes this; re-enable the logic below once the
                # minimum required version is bumped:
                #
                #   elif phase_h_high == "sc" and phase_h_low == "l":
                #       wrapper = c1.fluid_data[fluid]["wrapper"]
                #       h_at_tc = wrapper.h_pT(c1.p.val_SI, wrapper._T_crit)
                #       if np.isclose(h_at_tc, c1.h.val_SI):
                #           h_at_steps = [0, 1]
                #           zone_phases = [3]
                #       elif np.isclose(h_at_tc, c2.h.val_SI):
                #           h_at_steps = [0, 1]
                #           zone_phases = [0]
                #       else:
                #           x_tc = (h_at_tc - c1.h.val_SI) / delta_h
                #           h_at_steps = [0, x_tc, 1]
                #           zone_phases = [0, 3]
                zone_phases = [3]

            else:
                zone_phases = [_PHASE_TO_INT.get(phase_h_low, 2)]

        return h_at_steps, zone_phases

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
        num_steps = self.num_sections.val + 1
        if steps_hot is None:
            steps_hot, _ = self._get_moving_steps(self.inl[0], self.outl[0])
        if steps_cold is None:
            steps_cold, _ = self._get_moving_steps(self.inl[1], self.outl[1])
        return np.unique(np.r_[self._get_steps(num_steps), steps_hot, steps_cold])

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
        sections = self._calc_sections_SI(postprocess=False)
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
        heat exchanger UA_char, but for the sectioned heat exchanger.

        Returns
        -------
        float
            Residual value of equation:

            .. math::

                0 = UA_\text{design} * f_\text{UA} - \sum\left(UA_{i}\right)

        """
        p1 = self.UA_char1.param
        p2 = self.UA_char2.param

        f1 = self.get_char_expr(p1, **self.UA_char1.char_params)
        f2 = self.get_char_expr(p2, **self.UA_char2.char_params)

        fUA1 = self.UA_char1.char_func.evaluate(f1)
        fUA2 = self.UA_char2.char_func.evaluate(f2)

        fUA = 2 / (1 / fUA1 + 1 / fUA2)

        sections = self._calc_sections_SI(postprocess=False)
        min_td = self._min_td(sections)
        if min_td <= 0:
            return self.UA.design * fUA - self.calc_UA(sections) + min_td
        return self.UA.design * fUA - self.calc_UA(sections)

    def _UA_cecchinato_residual(self, re_exp_hot, re_exp_cold, hot_index, cold_index):
        alpha_ratio = self.alpha_ratio.val_SI
        area_ratio = self.area_ratio.val_SI
        m_ratio_hot = max(
            self.inl[hot_index].m.val_SI / self._conn_design(self.inl[hot_index], 'm'),
            1e-6
        )
        m_ratio_cold = max(
            self.inl[cold_index].m.val_SI / self._conn_design(self.inl[cold_index], 'm'),
            1e-6
        )
        fUA = (
            (1 + alpha_ratio * area_ratio)
            / (
                m_ratio_cold ** -re_exp_cold
                + alpha_ratio * area_ratio * m_ratio_hot ** -re_exp_hot
            )
        )
        sections = self._calc_sections_SI(postprocess=False)
        min_td = self._min_td(sections)
        if min_td <= 0:
            return self.UA.design * fUA - self.calc_UA(sections) + min_td
        return self.UA.design * fUA - self.calc_UA(sections)

    def UA_cecchinato_func(self):
        r"""
        Method to calculate heat transfer via UA design with modification
        for part load according to :cite:`cecchinato2010`. UA is determined
        over the UA values of the sections of the heat exchanger.

        Requires :code:`re_exp_hot`, :code:`re_exp_cold`, :code:`alpha_ratio`,
        and :code:`area_ratio`. Hot side is inlet index 0, cold side is inlet
        index 1.

        The modification factor for UA is calculated as follows

        .. math::

            f_\text{UA}=\frac{
                1 + \frac{\alpha_\text{hot}}{\alpha_\text{cold}}
                \cdot\frac{A_\text{hot}}{A_\text{cold}}
            }{
                \frac{\dot m_\text{cold}}{\dot m_\text{cold,ref}}^{-Re_\text{cold}} +
                \frac{\alpha_\text{hot}}{\alpha_\text{cold}}
                \cdot\frac{A_\text{hot}}{A_\text{cold}}
                \cdot\frac{\dot m_\text{hot}}{\dot m_\text{hot,ref}}^{-Re_\text{hot}}
            }

        Returns
        -------
        float
            residual value of equation

            .. math::

                0 = UA_\text{ref} \cdot f_\text{UA} - \sum UA_\text{i}
        """
        return self._UA_cecchinato_residual(
            re_exp_hot=self.re_exp_hot.val_SI,
            re_exp_cold=self.re_exp_cold.val_SI,
            hot_index=0,
            cold_index=1,
        )

    def UA_cecchinato_legacy_func(self):
        r"""
        Deprecated - use :code:`UA_cecchinato_hc` with :code:`re_exp_hot` and
        :code:`re_exp_cold` instead. In the next major version
        :code:`UA_cecchinato` will adopt the new hot/cold parameter convention
        and this group will be removed.

        Requires :code:`re_exp_r`, :code:`re_exp_sf`, :code:`alpha_ratio`,
        :code:`area_ratio`, and :code:`refrigerant_index`.

        Returns
        -------
        float
            residual value of equation
        """
        refrigerant_index = self.refrigerant_index.val
        if refrigerant_index == 0:
            hot_index, cold_index = 0, 1
            re_exp_hot = self.re_exp_r.val_SI
            re_exp_cold = self.re_exp_sf.val_SI
        else:
            hot_index, cold_index = 1, 0
            re_exp_hot = self.re_exp_sf.val_SI
            re_exp_cold = self.re_exp_r.val_SI
        return self._UA_cecchinato_residual(re_exp_hot, re_exp_cold, hot_index, cold_index)

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

    def _calc_UA_from_sections(self):
        return float(sum(self.Q_per_section.val_SI / self.lmtd_per_section.val_SI))

    def _calc_td_pinch(self):
        return float(min(self.T_hot_sections.val_SI - self.T_cold_sections.val_SI))

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

    @staticmethod
    def _section_phases(steps_all, steps_ref, zone_phases):
        """Return the phase index (0=L, 1=TP, 2=G, 3=SC) for each section on one
        fluid side.

        Parameters
        ----------
        steps_all : numpy.ndarray
            Combined step array from :py:meth:`_assign_steps`.
        steps_ref : numpy.ndarray
            Phase-boundary steps for this side from :py:meth:`_get_moving_steps`.
        zone_phases : list
            Zone phase list returned by :py:meth:`_get_moving_steps`.

        Returns
        -------
        list
            Phase index (0=L, 1=TP, 2=G, 3=SC) for each section.
        """
        boundaries = steps_ref[1:-1]
        return [
            zone_phases[int(np.searchsorted(boundaries, (steps_all[i] + steps_all[i + 1]) / 2))]
            for i in range(len(steps_all) - 1)
        ]

    def area_zones_func(self, **kwargs):
        r"""
        Residual for the Bell (2015) area-based heat exchanger constraint.

        Zones both fluid sides independently (SC, TP, SH) and requires that
        the sum of the zone areas equals the specified hot-side area
        :math:`A_h`:

        .. math::

            0 = A_h - \sum_j \frac{\dot{Q}_j}{U_j \cdot \text{LMTD}_j}

        where the overall heat transfer coefficient for zone :math:`j` is
        built from the individual convective coefficients following
        :cite:`bell2015`:

        .. math::

            U_j = \frac{1}{\frac{1}{\alpha_{h,j}} + A_h R_k
            + \frac{A_h}{\alpha_{c,j} A_c}}

        with :math:`R_k` the total wall thermal resistance in
        :math:`\text{K}/\text{W}` (:code:`R_cond`) and
        :math:`A_c = A_h \cdot` :code:`area_ratio`.

        Returns
        -------
        float
            Residual value of the equation.
        """
        steps1, zone_phases1 = self._get_moving_steps(self.inl[0], self.outl[0])
        steps2, zone_phases2 = self._get_moving_steps(self.inl[1], self.outl[1])
        steps_all = self._assign_steps(steps1, steps2)
        T_hot, T_cold = self._get_T_at_steps(steps_all)
        lmtd_per_section = self._calc_lmtd_per_section(T_hot, T_cold, postprocess=False)
        Q_per_section = np.diff(self._get_Q_cumsum_steps(steps_all))
        min_td = float(np.min(T_hot - T_cold))

        area_hot = self.area_hot.val_SI
        if min_td <= 0.0:
            # ×10: Newton overshoots past min_td=0 into the feasible side, giving oscillation_damping a sign-change bracket to bisect; ×1 lands at the branch discontinuity.
            return min_td * 20.0

        phases1 = self._section_phases(steps_all, np.array(steps1), zone_phases1)
        phases2 = self._section_phases(steps_all, np.array(steps2), zone_phases2)

        area_cold = area_hot * self.area_ratio.val_SI
        alpha1 = [self.alpha1_l.val_SI, self.alpha1_tp.val_SI, self.alpha1_g.val_SI, self.alpha1_sc.val_SI]
        alpha2 = [self.alpha2_l.val_SI, self.alpha2_tp.val_SI, self.alpha2_g.val_SI, self.alpha2_sc.val_SI]
        A_req = 0.0
        for Q_j, lmtd_j, ph1, ph2 in zip(Q_per_section, lmtd_per_section, phases1, phases2):
            U_j = 1.0 / (
                1.0 / alpha1[ph1]
                + area_hot * self.R_cond.val_SI
                + area_hot / (alpha2[ph2] * area_cold)
            )
            A_req += abs(Q_j) / (U_j * lmtd_j)

        residual = area_hot - A_req
        return residual

    def area_zones_dependents(self):
        return self.UA_dependents()


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
