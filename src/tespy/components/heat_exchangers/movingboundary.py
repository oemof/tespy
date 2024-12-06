# -*- coding: utf-8

"""Module of class MovingBoundaryCondenser.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/heat_exchangers/movingboundary.py

SPDX-License-Identifier: MIT
"""
import math

from tespy.components.heat_exchangers.base import HeatExchanger
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.global_vars import ERR
from tespy.tools.fluid_properties import h_mix_pQ


class MovingBoundaryCondenser(HeatExchanger):

    @staticmethod
    def component():
        return 'moving boundary condenser'

    def get_parameters(self):
        params = super().get_parameters()
        params.update({
            'U_desup': dc_cp(min_val=0),
            'U_cond': dc_cp(min_val=0),
            'U_subcool': dc_cp(min_val=0),
            'A': dc_cp(min_val=0),
            'U_sections_group': dc_gcp(
                elements=['U_desup', 'U_cond', 'U_subcool', 'A'],
                func=self.U_sections_func, deriv=self.U_sections_deriv, latex=None,
                num_eq=1
            ),
            'UA': dc_cp(
                min_val=0, num_eq=1, func=self.UA_func, deriv=self.UA_deriv
            ),
            'td_pinch': dc_cp(
                min_val=0, num_eq=1, func=self.td_pinch_func,
                deriv=self.td_pinch_deriv, latex=None
            )
        })
        return params

    def get_U_sections_and_h_steps(self, get_U_values=False):
        """Get the U values of the sections and the boundary hot side enthalpies

        Parameters
        ----------
        get_U_values : boolean
            Also return the U values for the sections of the heat exchanger.

        Returns
        -------
        tuple
            U values in the heat exchange sections and boundary hot side
            enthalpies
        """
        i1, _ = self.inl
        o1, _ = self.outl
        U_in_sections = []

        h_sat_gas = h_mix_pQ(o1.p.val_SI, 1, o1.fluid_data)
        h_sat_liquid = h_mix_pQ(o1.p.val_SI, 0, o1.fluid_data)

        if i1.h.val_SI > h_sat_gas + ERR and o1.h.val_SI < h_sat_liquid - ERR:
            h_at_steps_1 = [i1.h.val_SI, h_sat_gas, h_sat_liquid, o1.h.val_SI]
            if get_U_values:
                U_in_sections = [self.U_desup.val, self.U_cond.val, self.U_subcool.val]

        elif ((i1.h.val_SI > h_sat_gas + ERR) ^ (o1.h.val_SI < h_sat_liquid - ERR)):
            if i1.h.val_SI > h_sat_gas + ERR:
                h_at_steps_1 = [i1.h.val_SI, h_sat_gas, o1.h.val_SI]
                if get_U_values:
                    U_in_sections = [self.U_desup.val, self.U_cond.val]
            else:
                h_at_steps_1 = [i1.h.val_SI, h_sat_liquid, o1.h.val_SI]
                if get_U_values:
                    U_in_sections = [self.U_cond.val, self.U_subcool.val]

        else:
            h_at_steps_1 = [i1.h.val_SI, o1.h.val_SI]

            if get_U_values:
                if i1.h.val_SI > h_sat_gas + ERR:
                    U_in_sections = [self.U_desup.val]
                elif i1.h.val_SI > h_sat_liquid - ERR:
                    U_in_sections = [self.U_cond.val]
                else:
                    U_in_sections = [self.U_subcool.val]

        return U_in_sections, h_at_steps_1

    def calc_td_log_and_Q_in_sections(self, h_at_steps_1):
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
        i1, i2 = self.inl
        steps = len(h_at_steps_1)
        sections = steps - 1

        p_at_steps_1 = [i1.p.val_SI for _ in range(steps)]
        T_at_steps_1 = [
            T_mix_ph(p, h, i1.fluid_data, i1.mixing_rule)
            for p, h in zip(p_at_steps_1, h_at_steps_1)
        ]

        Q_in_sections = [
            i1.m.val_SI * (h_at_steps_1[i + 1] - h_at_steps_1[i])
            for i in range(sections)
        ]

        h_at_steps_2 = [i2.h.val_SI]
        for Q in Q_in_sections[::-1]:
            h_at_steps_2.append(h_at_steps_2[-1] + abs(Q) / i2.m.val_SI)

        p_at_steps_2 = [i2.p.val_SI for _ in range(sections + 1)]
        T_at_steps_2 = [
            T_mix_ph(p, h, i2.fluid_data, i2.mixing_rule)
            for p, h in zip(p_at_steps_2, h_at_steps_2)
        ]

        # counter flow version
        td_at_steps = [
            T1 - T2 for T1, T2 in zip(T_at_steps_1, T_at_steps_2[::-1])
        ]

        td_at_steps = [abs(td) for td in td_at_steps]
        td_log_in_sections = [
            (td_at_steps[i + 1] - td_at_steps[i])
            / math.log(td_at_steps[i + 1] / td_at_steps[i])
            for i in range(sections)
        ]
        return td_log_in_sections, Q_in_sections

    def calc_UA_in_sections(self):
        """Calc UA values for all sections.

        Returns
        -------
        list
            List of UA values starting from hot side inlet.
        """
        _, h_at_steps_1 = self.get_U_sections_and_h_steps()
        td_log_in_sections, Q_in_sections = self.calc_td_log_and_Q_in_sections(h_at_steps_1)
        UA_in_sections = [
            abs(Q) / td_log
            for Q, td_log in zip(Q_in_sections, td_log_in_sections)
        ]

        return UA_in_sections

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

        U_sections_specified = all([
            self.get_attr(f"U_{key}").is_set
            for key in ["desup", "cond", "subcool"]
        ])

        if U_sections_specified:
            U_in_sections, h_at_steps_1 = self.get_U_sections_and_h_steps(get_U_values=True)
            td_log_in_sections, Q_in_sections = self.calc_td_log_and_Q_in_sections(h_at_steps_1)
            self.A.val = self.Q.val ** 2 / (
                sum([
                    abs(Q) * td_log * U
                    for Q, td_log, U
                    in zip(Q_in_sections, td_log_in_sections, U_in_sections)
                ])
            )
            assert abs(abs(self.Q.val) / sum([
                ((Q / self.Q.val) * td_log * U)
                for Q, td_log, U
                in zip(Q_in_sections, td_log_in_sections, U_in_sections)
            ]) - self.A.val) < 1e-6
            assert round(sum([Q for Q in Q_in_sections]), 3) == round(self.Q.val, 3)

        self.UA.val = sum(self.calc_UA_in_sections())
        self.td_pinch.val = self.calc_td_pinch()
