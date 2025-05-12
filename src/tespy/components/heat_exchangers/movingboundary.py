# -*- coding: utf-8

"""Module of class MovingBoundaryCondenser.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/heat_exchangers/movingboundary.py

SPDX-License-Identifier: MIT
"""
import math

import numpy as np

from tespy.components.heat_exchangers.base import HeatExchanger
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.fluid_properties import h_mix_pQ
from tespy.tools.fluid_properties import single_fluid
from tespy.tools.global_vars import ERR
from tespy.tools.logger import logger


class MovingBoundaryHeatExchanger(HeatExchanger):
    r"""
    Class for counter current heat exchanger with UA sections.

    The heat exchanger is internally discretized into multiple sections, which
    are defined by phase changes. The component assumes, that no pressure
    losses occur. In principle the implementations follows :cite:`bell2015`.

    **Mandatory Equations**

    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_func`

    **Optional Equations**

    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_hot_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.kA_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.kA_char_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.ttd_u_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.ttd_l_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.ttd_min_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.eff_cold_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.eff_hot_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.eff_max_func`
    - hot side :py:meth:`tespy.components.component.Component.pr_func`
    - cold side :py:meth:`tespy.components.component.Component.pr_func`
    - hot side :py:meth:`tespy.components.component.Component.zeta_func`
    - cold side :py:meth:`tespy.components.component.Component.zeta_func`
    - hot side :py:meth:`tespy.components.component.Component.dp_func`
    - cold side :py:meth:`tespy.components.component.Component.dp_func`
    - :py:meth:`tespy.components.heat_exchangers.movingboundary.MovingBoundaryHeatExchanger.UA_func`
    - :py:meth:`tespy.components.heat_exchangers.movingboundary.MovingBoundaryHeatExchanger.td_pinch_func`

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

    kA : float, dict
        Area independent heat transfer coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char : dict
        Area independent heat transfer coefficient characteristic.

    kA_char1 : tespy.tools.characteristics.CharLine, dict
        Characteristic line for hot side heat transfer coefficient.

    kA_char2 : tespy.tools.characteristics.CharLine, dict
        Characteristic line for cold side heat transfer coefficient.

    UA : float, dict
        Sum of UA in all sections of the heat exchanger.

    td_pinch : float, dict
        Value of the lowest delta T between hot side and cold side at the
        different sections.

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
    >>> nw = Network(T_unit="C", p_unit="bar")
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

    >>> c1.set_attr(fluid={"Water": 1}, p=1, Td_bp=15, m=1)
    >>> c2.set_attr(p=1, Td_bp=-15)
    >>> c11.set_attr(fluid={"Air": 1}, p=1, T=15)
    >>> c12.set_attr(p=1, T=25)
    >>> nw.solve("design")

    Now we can remove the pressure specifications on the air side and impose
    the minimum pinch instead, which will determine the actual water
    condensation pressure.

    >>> c1.set_attr(p=None)
    >>> c2.set_attr(p=None)
    >>> c12.set_attr(p=None)
    >>> cd.set_attr(td_pinch=5, pr1=1, pr2=1)
    >>> nw.solve("design")
    >>> round(c1.p.val, 3)
    0.056

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

    >>> c2.set_attr(Td_bp=-5)
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
    """

    @staticmethod
    def component():
        return 'moving boundary heat exchanger'

    def get_parameters(self):
        params = super().get_parameters()
        params.update({
            'U_gas_gas': dc_cp(min_val=0),
            'U_gas_twophase': dc_cp(min_val=0),
            'U_gas_liquid': dc_cp(min_val=0),
            'U_liquid_gas': dc_cp(min_val=0),
            'U_liquid_twophase': dc_cp(min_val=0),
            'U_liquid_liquid': dc_cp(min_val=0),
            'U_twophase_gas': dc_cp(min_val=0),
            'U_twophase_twophase': dc_cp(min_val=0),
            'U_twophase_liquid': dc_cp(min_val=0),
            'A': dc_cp(min_val=0),
            'UA': dc_cp(
                min_val=0, num_eq=1, func=self.UA_func, deriv=self.UA_deriv
            ),
            'td_pinch': dc_cp(
                min_val=0, num_eq=1, func=self.td_pinch_func,
                deriv=self.td_pinch_deriv, latex=None
            )
        })
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
        if c1.fluid_data != c2.fluid_data:
            msg = "Both connections need to utilize the same fluid data."
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
        h_steps_hot = self._get_h_steps(self.inl[0], self.outl[0])
        Q_sections_hot = self._get_Q_sections(h_steps_hot, self.inl[0].m.val_SI)
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

    @staticmethod
    def _calc_td_log_per_section(T_steps_hot, T_steps_cold):
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

    def calc_sections(self):
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
        td_log_per_section = self._calc_td_log_per_section(T_steps_hot, T_steps_cold)
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
        Calculate heat transfer from heat transfer coefficients for
        desuperheating and condensation as well as total heat exchange area.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = UA - \sum UA_\text{i}
        """
        sections = self.calc_sections()
        return self.UA.val - self.calc_UA(sections)

    def UA_deriv(self, increment_filter, k):
        r"""
        Partial derivatives of heat transfer coefficient function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        f = self.UA_func
        for c in self.inl + self.outl:
            if self.is_variable(c.m):
                self.jacobian[k, c.m.J_col] = self.numeric_deriv(f, "m", c)
            if self.is_variable(c.p):
                self.jacobian[k, c.p.J_col] = self.numeric_deriv(f, 'p', c)
            if self.is_variable(c.h):
                self.jacobian[k, c.h.J_col] = self.numeric_deriv(f, 'h', c)

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
        Equation for pinch point temperature difference of a condenser.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = td_\text{pinch} - min(td_\text{sections})
        """
        sections = self.calc_sections()
        return self.td_pinch.val - self.calc_td_pinch(sections)

    def td_pinch_deriv(self, increment_filter, k):
        """
        Calculate partial derivates of upper terminal temperature function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        f = self.td_pinch_func
        for c in self.inl + self.outl:
            if self.is_variable(c.m, increment_filter):
                self.jacobian[k, c.m.J_col] = self.numeric_deriv(f, 'm', c)
            if self.is_variable(c.p, increment_filter):
                self.jacobian[k, c.p.J_col] = self.numeric_deriv(f, 'p', c)
            if self.is_variable(c.h, increment_filter):
                self.jacobian[k, c.h.J_col] = self.numeric_deriv(f, 'h', c)

    def calc_parameters(self):
        super().calc_parameters()

        sections = self.calc_sections()
        self.UA.val = self.calc_UA(sections)
        self.td_pinch.val = self.calc_td_pinch(sections)

        if round(self.inl[0].p.val_SI) != round(self.outl[0].p.val_SI):
            msg = (
                f"The {self.__class__.__name__} instance {self.label} is "
                "defined for constant pressure. The identification of the "
                "heat transfer sections might be wrong in case phase changes "
                "are involved in the heat transfer process."
            )
            logger.warning(msg)
        if round(self.inl[1].p.val_SI) != round(self.outl[1].p.val_SI):
            msg = (
                f"The {self.__class__.__name__} instance {self.label} is "
                "defined for constant pressure. The identification of the "
                "heat transfer sections might be wrong in case phase changes "
                "are involved in the heat transfer process."
            )
            logger.warning(msg)
