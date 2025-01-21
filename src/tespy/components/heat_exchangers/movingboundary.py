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
from tespy.tools.fluid_properties import h_mix_pQ
from tespy.tools.fluid_properties import single_fluid
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.global_vars import ERR


class MovingBoundaryHeatExchanger(HeatExchanger):

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
            # 'U_sections_group': dc_gcp(
            #     elements=['U_desup', 'U_cond', 'U_subcool', 'A'],
            #     func=self.U_sections_func, deriv=self.U_sections_deriv, latex=None,
            #     num_eq=1
            # ),
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
        """Get the steps for enthalpy for a change of state from one connection
        to another

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
        # this should be generalized to "supports two-phase"
        is_pure_fluid = fluid is not None

        if is_pure_fluid:
            try:
                h_sat_gas = h_mix_pQ(c1.p.val_SI, 1, c1.fluid_data)
                h_sat_liquid = h_mix_pQ(c1.p.val_SI, 0, c1.fluid_data)
            except (ValueError, NotImplementedError):
                return h_at_steps

            if c1.h.val_SI < h_sat_liquid:
                if c2.h.val_SI > h_sat_gas:
                    h_at_steps = [c1.h.val_SI, h_sat_liquid, h_sat_gas, c2.h.val_SI]
                elif c2.h.val_SI > h_sat_liquid:
                    h_at_steps = [c1.h.val_SI, h_sat_liquid, c2.h.val_SI]

            elif c1.h.val_SI < h_sat_gas:
                if c2.h.val_SI > h_sat_gas:
                    h_at_steps = [c1.h.val_SI, h_sat_gas, c2.h.val_SI]

        return h_at_steps

    @staticmethod
    def _get_Q_sections(h_at_steps, mass_flow):
        return [
            (h_at_steps[i + 1] - h_at_steps[i]) * mass_flow
            for i in range(len(h_at_steps) - 1)
        ]

    def _assign_sections(self):
        h_steps_hot = self._get_h_steps(self.inl[0], self.outl[0])
        Q_sections_hot = self._get_Q_sections(h_steps_hot, self.inl[0].m.val_SI)
        Q_sections_hot = np.cumsum(Q_sections_hot).round(6).tolist()

        h_steps_cold = self._get_h_steps(self.inl[1], self.outl[1])
        Q_sections_cold = self._get_Q_sections(h_steps_cold, self.inl[1].m.val_SI)
        Q_sections_cold = np.cumsum(Q_sections_cold).round(6).tolist()

        all_sections = [Q for Q in Q_sections_hot + Q_sections_cold + [0.0]]
        return sorted(list(set(all_sections)))

    def _get_T_at_steps(self, Q_sections):
        # now put the Q_sections back on the h_steps on both sides
        h_steps_hot = [self.inl[0].h.val_SI - Q / self.inl[0].m.val_SI for Q in Q_sections]
        h_steps_cold = [self.inl[1].h.val_SI + Q / self.inl[1].m.val_SI for Q in Q_sections]

        T_steps_hot = [
            T_mix_ph(self.inl[0].p.val_SI, h, self.inl[0].fluid_data, self.inl[0].mixing_rule)
            for h in h_steps_hot
        ]
        T_steps_cold = [
            T_mix_ph(self.inl[1].p.val_SI, h, self.inl[1].fluid_data, self.inl[1].mixing_rule)
            for h in h_steps_cold
        ]
        return T_steps_hot, T_steps_cold

    @staticmethod
    def _calc_UA_in_sections(T_steps_hot, T_steps_cold, Q_sections):
        """Calculate logarithmic temperature difference and heat exchange in
        heat exchanger sections.

        Parameters
        ----------
        h_at_steps_1 : list
            Enthalpy values at boundaries of sections.

        Returns
        -------
        tuple
            Lists of logarithmic temperature difference and heat exchange in
            the heat exchanger sections starting from hot side inlet.
        """
        # counter flow version
        td_at_steps = [
            T_hot - T_cold
            for T_hot, T_cold in zip(T_steps_hot, T_steps_cold[::-1])
        ]

        td_at_steps = [abs(td) for td in td_at_steps]
        td_log_in_sections = [
            (td_at_steps[i + 1] - td_at_steps[i])
            / math.log(td_at_steps[i + 1] / td_at_steps[i])
            if td_at_steps[i + 1] != td_at_steps[i] else td_at_steps[i + 1]
            for i in range(len(Q_sections))
        ]
        UA_in_sections = [
            abs(Q) / td_log
            for Q, td_log in zip(Q_sections, td_log_in_sections)
        ]
        return UA_in_sections

    def calc_UA(self):
        Q_sections = self._assign_sections()
        T_steps_hot, T_steps_cold = self._get_T_at_steps(Q_sections)
        Q_per_section = np.diff(Q_sections)
        UA_sections = self._calc_UA_in_sections(T_steps_hot, T_steps_cold, Q_per_section)
        return sum(UA_sections)

    def UA_func(self, **kwargs):
        r"""
        Calculate heat transfer from heat transfer coefficients for
        desuperheating and condensation as well as total heat exchange area.

        Returns
        -------
        residual : float
            Residual value of equation.
        """
        UA_in_sections = self.calc_UA_in_sections()
        return self.UA.val - sum(UA_in_sections)

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

    def U_sections_func(self, **kwargs):
        r"""
        Calculate heat transfer from heat transfer coefficients for
        desuperheating and condensation as well as total heat exchange area.

        Returns
        -------
        residual : float
            Residual value of equation.
        """
        U_in_sections, h_at_steps_1 = self.get_U_sections_and_h_steps(get_U_values=True)
        td_log_in_sections, Q_in_sections = self.calc_td_log_and_Q_in_sections(h_at_steps_1)

        Q_total = sum(Q_in_sections)

        return (
            Q_total
            + self.A.val / Q_total
            * sum([
                    Q * td_log * U
                    for Q, td_log, U
                    in zip(Q_in_sections, td_log_in_sections, U_in_sections)
            ])
        )

    def U_sections_deriv(self, increment_filter, k):
        r"""
        Partial derivatives of heat transfer coefficient function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        f = self.U_sections_func
        for c in self.inl + self.outl:
            if self.is_variable(c.m):
                self.jacobian[k, c.m.J_col] = self.numeric_deriv(f, "m", c)
            if self.is_variable(c.p):
                self.jacobian[k, c.p.J_col] = self.numeric_deriv(f, 'p', c)
            if self.is_variable(c.h):
                self.jacobian[k, c.h.J_col] = self.numeric_deriv(f, 'h', c)

    def calc_td_pinch(self):
        """Calculate the pinch point temperature difference

        Returns
        -------
        float
            Value of the pinch point temperature difference
        """
        o1 = self.outl[0]
        i2 = self.inl[1]

        h_sat = h_mix_pQ(o1.p.val_SI, 1, o1.fluid_data)

        if o1.h.val_SI < h_sat:
            # we have two sections in this case
            Q_cond = o1.m.val_SI * (o1.h.val_SI - h_sat)

            # calculate the intermediate temperatures
            T_cond_i1 =  o1.calc_T_sat()
            h_cond_o2 = i2.h.val_SI + abs(Q_cond) / i2.m.val_SI
            T_cond_o2 = T_mix_ph(i2.p.val_SI, h_cond_o2, i2.fluid_data)

            return T_cond_i1 - T_cond_o2

        else:
            o2 = self.outl[1]
            return o1.calc_T_sat() - o2.calc_T()

    def td_pinch_func(self):
        r"""
        Equation for pinch point temperature difference of a condenser.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = td_\text{pinch} - T_\text{sat,in,1}
                + T_\left(
                    p_\text{in,2},\left[
                        h_\text{in,2}
                        + \frac{\dot Q_\text{cond}}{\dot m_\text{in,2}}
                    \right]
                \right)
        """
        return self.td_pinch.val - self.calc_td_pinch()

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
        for c in [self.outl[0], self.inl[1]]:
            if self.is_variable(c.m, increment_filter):
                self.jacobian[k, c.m.J_col] = self.numeric_deriv(f, 'm', c)
            if self.is_variable(c.p, increment_filter):
                self.jacobian[k, c.p.J_col] = self.numeric_deriv(f, 'p', c)
            if self.is_variable(c.h, increment_filter):
                self.jacobian[k, c.h.J_col] = self.numeric_deriv(f, 'h', c)

    def calc_parameters(self):
        super().calc_parameters()

        UA = self.calc_UA()
        print(UA)
        # U_sections_specified = all([
        #     self.get_attr(f"U_{key}").is_set
        #     for key in ["desup", "cond", "subcool"]
        # ])

        # if U_sections_specified:
        #     U_in_sections, h_at_steps_1 = self.get_U_sections_and_h_steps(get_U_values=True)
        #     td_log_in_sections, Q_in_sections = self.calc_td_log_and_Q_in_sections(h_at_steps_1)
        #     self.A.val = self.Q.val ** 2 / (
        #         sum([
        #             abs(Q) * td_log * U
        #             for Q, td_log, U
        #             in zip(Q_in_sections, td_log_in_sections, U_in_sections)
        #         ])
        #     )
        #     assert abs(abs(self.Q.val) / sum([
        #         ((Q / self.Q.val) * td_log * U)
        #         for Q, td_log, U
        #         in zip(Q_in_sections, td_log_in_sections, U_in_sections)
        #     ]) - self.A.val) < 1e-6
        #     assert round(sum([Q for Q in Q_in_sections]), 3) == round(self.Q.val, 3)

        # self.UA.val = sum(self.calc_UA_in_sections())
        # self.td_pinch.val = self.calc_td_pinch()
