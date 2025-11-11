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
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.fluid_properties import h_mix_pQ
from tespy.tools.fluid_properties import single_fluid
from tespy.tools.global_vars import ERR
from tespy.tools.logger import logger


@component_registry
class MovingBoundaryHeatExchanger(SectionedHeatExchanger):
    r"""
    Class for counter flow heat exchanger with UA sections.

    The heat exchanger is internally discretized into multiple sections, which
    are defined by phase changes. The component assumes, that no pressure
    losses occur. In principle the implementations follows :cite:`bell2015`.

    **Mandatory Equations**

    - fluid: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`
    - mass flow: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_func`

    **Optional Equations**

    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_hot_func`
    - :py:meth:`tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger.UA_func`
    - :py:meth:`tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger.td_pinch_func`
    - :py:meth:`tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger.UA_cecchinato_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.ttd_u_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.ttd_l_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.ttd_min_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.eff_cold_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.eff_hot_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.eff_max_func`

    For hot and cold side individually:

    - :py:meth:`tespy.components.component.Component.pr_structure_matrix`
    - :py:meth:`tespy.components.component.Component.dp_structure_matrix`
    - :py:meth:`tespy.components.component.Component.zeta_func`

    Inlets/Outlets

    - in1, in2 (index 1: hot side, index 2: cold side)
    - out1, out2 (index 1: hot side, index 2: cold side)

    Image

    .. image:: /api/_images/HeatExchanger.svg
       :alt: flowsheet of the heat exchanger
       :align: center
       :class: only-light

    .. image:: /api/_images/HeatExchanger_darkmode.svg
       :alt: flowsheet of the heat exchanger
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

    Q : float, dict
        Heat transfer, :math:`Q/\text{W}`.

    pr1 : float, dict, :code:`"var"`
        Outlet to inlet pressure ratio at hot side, :math:`pr/1`.

    pr2 : float, dict, :code:`"var"`
        Outlet to inlet pressure ratio at cold side, :math:`pr/1`.

    dp1 : float, dict, :code:`"var"`
        Inlet to outlet pressure delta at hot side, unit is the network's
        pressure unit!.

    dp2 : float, dict, :code:`"var"`
        Inlet to outlet pressure delta at cold side, unit is the network's
        pressure unit!.

    zeta1 : float, dict, :code:`"var"`
        Geometry independent friction coefficient at hot side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    zeta2 : float, dict, :code:`"var"`
        Geometry independent friction coefficient at cold side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    ttd_l : float, dict
        Lower terminal temperature difference :math:`ttd_\text{l}/\text{K}`.

    ttd_u : float, dict
        Upper terminal temperature difference :math:`ttd_\text{u}/\text{K}`.

    ttd_min : float, dict
        Minimum terminal temperature difference :math:`ttd_\text{min}/\text{K}`.

    eff_cold : float, dict
        Cold side heat exchanger effectiveness :math:`eff_\text{cold}/\text{1}`.

    eff_hot : float, dict
        Hot side heat exchanger effectiveness :math:`eff_\text{hot}/\text{1}`.

    eff_max : float, dict
        Max value of hot and cold side heat exchanger effectiveness values
        :math:`eff_\text{max}/\text{1}`.

    UA : float, dict
        Sum of UA in all sections of the heat exchanger.

    td_pinch : float, dict
        Value of the lowest delta T between hot side and cold side at the
        different sections.

    UA_cecchinato : dict
        Group specification for partload UA modification according to
        :cite:`cecchinato2010`, for usage see details in the
        :py:meth:`tespy.components.heat_exchangers.movingboundary.MovingBoundaryHeatExchanger.UA_cecchinato_func`.
        This method can only be used in offdesign simulations!

    alpha_ration: float
        Secondary fluid to refrigerant heat transfer coefficient ratio.

    area_ration: float
        Secondary fluid to refrigerant heat transfer area ratio.

    re_exp_r: float
        Reynolds exponent for refrigerant side.

    re_exp_sf: float
        Reynolds exponent for secondary fluid side.

    refrigerant_index: int
        Connection index for the refrigerant side, 0 if refrigerant is on hot
        side, 1 if refrigerant is on cold side.

    Note
    ----
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
    ...     "pressure": "bar", "temperature": "degC"
    ... })
    >>> nw.set_attr(iterinfo=False)
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

    We can also see the temperature differences in all sections of the heat
    exchanger. Since the water vapor is cooled, condensed and then subcooled,
    while the air does not change phase, three sections will form:

    >>> Q_sections, T_steps_hot, T_steps_cold, Q_per_section, td_log_per_section = cd.calc_sections()
    >>> T_steps_hot, T_steps_cold = cd._get_T_at_steps(Q_sections)
    >>> delta_T_between_sections = T_steps_hot - T_steps_cold
    >>> [round(float(dT), 2) for dT in delta_T_between_sections]
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
    >>> Q_sections, T_steps_hot, T_steps_cold, Q_per_section, td_log_per_section = cd.calc_sections()
    >>> T_steps_hot, T_steps_cold = cd._get_T_at_steps(Q_sections)
    >>> delta_T_between_sections = T_steps_hot - T_steps_cold
    >>> [round(float(dT), 2) for dT in delta_T_between_sections]
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

    - refrigerant side Reynolds exponent
    - secondary fluid side Reynolds exponent
    - secondary fluid to refrigerant area ratio
    - secondary fluid to refrigerant alpha (heat transfer coefficient) ratio
    - the refrigerant index (which side of the heat exchanger is passed by the
      refrigerant)

    >>> import os
    >>> nw.save("design.json")
    >>> cd.set_attr(
    ...     area_ratio=20,        # typical for a finned heat exchanger
    ...     alpha_ratio=1e-2,     # alpha for water side is higher
    ...     re_exp_r=0.8,
    ...     re_exp_sf=0.55,
    ...     refrigerant_index=0,  # water is refrigerant in this case
    ...     design=["td_pinch"],
    ...     offdesign=["UA_cecchinato"]
    ... )
    >>> nw.solve("offdesign", design_path="design.json")

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
    >>> nw.solve("offdesign", design_path="design.json")
    >>> round(cd.Q.val_SI / cd.Q.design, 2)
    0.8
    >>> round(cd.UA.val_SI / cd.UA.design, 2)
    0.88
    >>> round(cd.td_pinch.val, 2)
    4.3
    >>> os.remove("design.json")
    """
    def get_parameters(self):
        params = super().get_parameters()
        del params["num_sections"]
        return params

    @staticmethod
    def _get_h_steps(c1, c2):
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

        h_at_steps = [c1.h.val_SI, c2.h.val_SI]
        fluid = single_fluid(c1.fluid_data)
        # this should be generalized to "supports two-phase" because it does
        # not work with incompressibles
        is_pure_fluid = fluid is not None

        if is_pure_fluid:
            try:
                h_sat_gas = h_mix_pQ(c1.p.val_SI, 1, c1.fluid_data)
                h_sat_liquid = h_mix_pQ(c1.p.val_SI, 0, c1.fluid_data)
            except (ValueError, NotImplementedError):
                return h_at_steps

            if c1.h.val_SI < h_sat_liquid:
                if c2.h.val_SI > h_sat_gas + ERR:
                    h_at_steps = [c1.h.val_SI, h_sat_liquid, h_sat_gas, c2.h.val_SI]
                elif c2.h.val_SI > h_sat_liquid + ERR:
                    h_at_steps = [c1.h.val_SI, h_sat_liquid, c2.h.val_SI]

            elif c1.h.val_SI < h_sat_gas - ERR:
                if c2.h.val_SI > h_sat_gas + ERR:
                    h_at_steps = [c1.h.val_SI, h_sat_gas, c2.h.val_SI]

        return h_at_steps

    def _assign_sections(self):
        """Assign the sections of the heat exchanger

        Returns
        -------
        list
            List of cumulative sum of heat exchanged defining the heat exchanger
            sections.
        """
        h_steps_hot = self._get_h_steps(self.inl[0], self.outl[0])
        Q_sections_hot = self._get_Q_sections(h_steps_hot, self.inl[0].m.val_SI)
        # do not insert last section, that will come from other side
        Q_sections_hot = np.insert(np.cumsum(Q_sections_hot)[:-1], 0, 0.0)

        h_steps_cold = self._get_h_steps(self.inl[1], self.outl[1])
        Q_sections_cold = self._get_Q_sections(h_steps_cold, self.inl[1].m.val_SI)
        Q_sections_cold = np.cumsum(Q_sections_cold)

        return np.sort(np.r_[Q_sections_cold, Q_sections_hot])

    def _get_T_at_steps(self, Q_sections):
        """Calculate the temperature values for the provided sections.

        Parameters
        ----------
        Q_sections : list
            Cumulative heat exchanged from the hot side to the colde side
            defining the sections of the heat exchanger.

        Returns
        -------
        tuple
            Lists of cold side and hot side temperature
        """
        # now put the Q_sections back on the h_steps on both sides
        # Since Q_sections is defined increasing we have to start back from the
        # outlet of the hot side
        h_steps_hot = self.outl[0].h.val_SI + Q_sections / self.inl[0].m.val_SI
        h_steps_cold = self.inl[1].h.val_SI + Q_sections / self.inl[1].m.val_SI
        T_steps_hot = np.array([
            T_mix_ph(self.inl[0].p.val_SI, h, self.inl[0].fluid_data, self.inl[0].mixing_rule)
            for h in h_steps_hot
        ])
        T_steps_cold = np.array([
            T_mix_ph(self.inl[1].p.val_SI, h, self.inl[1].fluid_data, self.inl[1].mixing_rule)
            for h in h_steps_cold
        ])
        return T_steps_hot, T_steps_cold

    def calc_parameters(self):
        super().calc_parameters()

        if round(self.inl[0].p.val_SI) != round(self.outl[0].p.val_SI):
            msg = (
                f"The {self.__class__.__name__} instance {self.label} is "
                "discovering the phase changes based on constant pressure "
                "assumption. The identification of the heat transfer sections "
                "might be wrong in case phase changes are involved in the "
                "heat transfer process."
            )
            logger.warning(msg)
        if round(self.inl[1].p.val_SI) != round(self.outl[1].p.val_SI):
            msg = (
                f"The {self.__class__.__name__} instance {self.label} is "
                "discovering the phase changes based on constant pressure "
                "assumption. The identification of the heat transfer sections "
                "might be wrong in case phase changes are involved in the "
                "heat transfer process."
            )
            logger.warning(msg)
