# -*- coding: utf-8

"""Module of class MovingBoundaryCondenser.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/heat_exchangers/sectioned.py

SPDX-License-Identifier: MIT
"""
import math

import numpy as np

from tespy.components.heat_exchangers.base import HeatExchanger
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
from tespy.tools.fluid_properties import T_mix_ph


class SectionedHeatExchanger(HeatExchanger):
    r"""
    Class for counter flow heat exchanger with UA sections.

    The heat exchanger is internally discretized into 51 sections of equal heat
    transfer. The number of section can be adjusted by the user. It is based on
    the model implemented by :cite:`quoilin2020`.

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
        Lower terminal temperature difference :math:`ttd_\mathrm{l}/\text{K}`.

    ttd_u : float, dict
        Upper terminal temperature difference :math:`ttd_\mathrm{u}/\text{K}`.

    ttd_min : float, dict
        Minumum terminal temperature difference :math:`ttd_\mathrm{min}/\text{K}`.

    eff_cold : float, dict
        Cold side heat exchanger effectiveness :math:`eff_\text{cold}/\text{1}`.

    eff_hot : float, dict
        Hot side heat exchanger effectiveness :math:`eff_\text{hot}/\text{1}`.

    eff_max : float, dict
        Max value of hot and cold side heat exchanger effectiveness values
        :math:`eff_\text{max}/\text{1}`.

    UA : float, dict
        Sum of UA in all sections of the heat exchanger.

    UA_cecchinato : dict
        Group specification for partload UA modification according to
        :cite:`cecchinato2010`, for usage see details in the
        :py:meth:`tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger.UA_cecchinato_func`.
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

    td_pinch : float, dict
        Value of the lowest delta T between hot side and cold side at the
        different sections.

    Note
    ----
    The equations only apply to counter-current heat exchangers.

    Example
    -------
    Water vapor should be cooled down, condensed and then further subcooled.
    For his air is heated up from 15 °C to 25 °C.

    >>> from tespy.components import Source, Sink, SectionedHeatExchanger
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
    >>> T_steps_hot, T_steps_cold = cd._get_T_at_steps(Q_sections)
    >>> delta_T_between_sections = T_steps_hot - T_steps_cold
    >>> delta_T_list = [round(float(dT), 2) for dT in delta_T_between_sections]
    >>> delta_T_list[::10]
    [5.0, 18.0, 16.0, 14.0, 12.0, 25.0]

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
    >>> delta_T_list = [round(float(dT), 2) for dT in delta_T_between_sections]
    >>> delta_T_list[::10]
    [9.8, 12.8, 10.8, 8.8, 6.8, 19.8]

    Finally, in contrast to the baseclass :code:`HeatExchanger` `kA` value, the
    `UA` value takes into account the heat transfer per section and calculates
    the heat transfer coefficient as the sum of all sections, while the `kA`
    value only takes into account the inlet and outlet temperatures and the
    total heat transfer.

    >>> round(cd.kA.val)
    174399
    >>> round(cd.UA.val)
    274461

    It is also possible to apply a part-load modification to UA following the
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
    274461

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
        params.update({
            'num_sections': dc_simple(val=51),
            'UA': dc_cp(
                min_val=0, num_eq_sets=1,
                func=self.UA_func,
                dependents=self.UA_dependents,
                quantity="heat_transfer_coefficient"
            ),
            'refrigerant_index': dc_simple(val=0),
            're_exp_r': dc_cp(),
            're_exp_sf': dc_cp(),
            'alpha_ratio': dc_cp(),
            'area_ratio': dc_cp(),
            'UA_cecchinato': dc_gcp(
                elements=['re_exp_r', 're_exp_sf', 'alpha_ratio', 'area_ratio'],
                num_eq_sets=1,
                func=self.UA_cecchinato_func,
                dependents=self.UA_dependents,
            ),
            'td_pinch': dc_cp(
                min_val=0, num_eq_sets=1,
                func=self.td_pinch_func,
                dependents=self.td_pinch_dependents,
                quantity="temperature_difference"
            )
        })
        return params

    @staticmethod
    def _get_h_steps(c1, c2, num_sections=51):
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

        # change the order of connections to have c1 as the lower enthalpy
        # connection (enthalpy will be rising in the list)
        if c1.h.val_SI > c2.h.val_SI:
            c1, c2 = c2, c1

        h_at_steps = np.linspace(c1.h.val_SI, c2.h.val_SI, num_sections)

        return h_at_steps

    @staticmethod
    def _get_Q_sections(h_at_steps, mass_flow):
        """Calculate the heat exchange of every section given steps of
        enthalpy and mass flow.

        Parameters
        ----------
        h_at_steps : list
            Enthalpy values at sections (inlet, phase change points, outlet)
        mass_flow : float
            Mass flow value

        Returns
        -------
        float
            Heat exchanged between defined steps of enthalpy.
        """
        return np.diff(h_at_steps) * mass_flow

    def _assign_sections(self):
        """Assign the sections of the heat exchanger

        Returns
        -------
        list
            List of cumulative sum of heat exchanged defining the heat exchanger
            sections.
        """
        h_steps_hot = self._get_h_steps(self.inl[0], self.outl[0], self.num_sections.val)
        Q_sections_hot = self._get_Q_sections(h_steps_hot, self.inl[0].m.val_SI)
        return np.insert(np.cumsum(Q_sections_hot), 0, 0.0)

    def _get_T_at_steps(self, Q_sections):
        """Calculate the temperature values for the provided sections.

        Parameters
        ----------
        Q_sections : list
            Cumulative heat exchanged from the hot side to the cold side
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
        p_steps_hot = np.linspace(
            self.outl[0].p.val_SI, self.inl[0].p.val_SI, self.num_sections.val
        )
        p_steps_cold = np.linspace(
            self.inl[1].p.val_SI, self.outl[1].p.val_SI, self.num_sections.val
        )
        T_steps_hot = np.array([
            T_mix_ph(p, h, self.inl[0].fluid_data, self.inl[0].mixing_rule)
            for p, h in zip(p_steps_hot, h_steps_hot)
        ])
        T_steps_cold = np.array([
            T_mix_ph(p, h, self.inl[1].fluid_data, self.inl[1].mixing_rule)
            for p, h in zip(p_steps_cold, h_steps_cold)
        ])
        return T_steps_hot, T_steps_cold

    @staticmethod
    def _calc_td_log_per_section(T_steps_hot, T_steps_cold, postprocess=False):
        """Calculate the logarithmic temperature difference values per section
        of heat exchanged.

        Parameters
        ----------
        T_steps_hot : list
            Temperature hot side at beginning and end of sections.

        T_steps_cold : list
            Temperature cold side at beginning and end of sections.

        Returns
        -------
        list
            Lists of temperature differences per section of heat exchanged.
        """
        if postprocess:
            td_at_steps = T_steps_hot - T_steps_cold
            if (td_at_steps <= 0).any():
                return np.ones(len(td_at_steps) - 1) * np.nan
        # the temperature ranges both come with increasing values
        td_at_steps = np.abs(T_steps_hot - T_steps_cold)

        return np.array([
            (td_at_steps[i + 1] - td_at_steps[i])
            / math.log(td_at_steps[i + 1] / td_at_steps[i])
            # round is required because tiny differences may cause
            # inconsistencies due to rounding errors
            if round(td_at_steps[i + 1], 6) != round(td_at_steps[i], 6)
            else td_at_steps[i + 1]
            for i in range(len(td_at_steps) - 1)
        ])

    def calc_sections(self, postprocess=True):
        """Calculate the sections of the heat exchanger.

        Returns
        -------
        tuple
            Cumulated heat transfer over sections, temperature at steps hot
            side, temperature at steps cold side, heat transfer per section
        """
        Q_sections = self._assign_sections()
        T_steps_hot, T_steps_cold = self._get_T_at_steps(Q_sections)
        Q_per_section = np.diff(Q_sections)
        td_log_per_section = self._calc_td_log_per_section(T_steps_hot, T_steps_cold, postprocess)
        return Q_sections, T_steps_hot, T_steps_cold, Q_per_section, td_log_per_section

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

    def UA_func(self, **kwargs):
        r"""
        Residual method for fixed heat transfer coefficient UA.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = UA - \sum UA_\text{i}
        """
        sections = self.calc_sections(False)
        return self.UA.val_SI - self.calc_UA(sections)

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
        m_ratio_r = m_r.val_SI / m_r.design
        m_sf = self.inl[secondary_index].m
        m_ratio_sf = m_sf.val_SI / m_sf.design

        fUA = (
            (1 + alpha_ratio * area_ratio)
            / (
                m_ratio_sf ** -re_exp_sf
                + alpha_ratio * area_ratio * m_ratio_r ** -re_exp_r
            )
        )
        sections = self.calc_sections(False)
        return self.UA.design * fUA - self.calc_UA(sections)

    def UA_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].h,
            self.inl[1].m,
            self.inl[1].p,
            self.inl[1].h,
            self.outl[1].h
        ]

    def calc_td_pinch(self, sections):
        """Calculate the pinch point temperature difference

        Returns
        -------
        float
            Value of the pinch point temperature difference
        """
        _, T_steps_hot, T_steps_cold, _, _ = sections

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
        sections = self.calc_sections(False)
        return self.td_pinch.val_SI - self.calc_td_pinch(sections)

    def td_pinch_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].h,
            self.inl[1].m,
            self.inl[1].p,
            self.inl[1].h,
            self.outl[1].h
        ]

    def calc_parameters(self):
        super().calc_parameters()

        sections = self.calc_sections()
        self.UA.val_SI = self.calc_UA(sections)
        self.td_pinch.val_SI = self.calc_td_pinch(sections)
