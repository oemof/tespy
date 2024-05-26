# -*- coding: utf-8

"""Module of class MovingBoundaryCondenser.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/heat_exchangers/movingboundary.py

SPDX-License-Identifier: MIT
"""
import numpy as np

from tespy.components.heat_exchangers.base import HeatExchanger
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.fluid_properties import T_mix_ph
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
            'A': dc_cp(min_val=0),
            'UA_group': dc_gcp(
                elements=['U_desup', 'U_cond', 'A'],
                func=self.UA_func, deriv=self.UA_deriv, latex=None,
                num_eq=1
            ),
            'td_pinch': dc_cp(
                min_val=0, num_eq=1, func=self.td_pinch_func,
                deriv=self.td_pinch_deriv, latex=None
            )
        })
        return params

    def UA_func(self, **kwargs):
        r"""
        Calculate heat transfer from heat transfer coefficients for
        desuperheating and condensation as well as total heat exchange area.

        Returns
        -------
        residual : float
            Residual value of equation.
        """
        Q_total = self.inl[0].m.val_SI * (self.outl[0].h.val_SI - self.inl[0].h.val_SI)
        h_sat = h_mix_pQ(self.outl[0].p.val_SI, 1, self.outl[0].fluid_data)

        if self.outl[0].h.val_SI < h_sat:
            # we have two sections in this case
            Q_desup = self.inl[0].m.val_SI * (h_sat - self.inl[0].h.val_SI)
            Q_cond = self.inl[0].m.val_SI * (self.outl[0].h.val_SI - h_sat)

            # calculate the intermediate temperatures
            T_desup_i1 = self.inl[0].calc_T()
            T_desup_o1 =  self.outl[0].calc_T_sat()
            T_cond_i1 = T_desup_o1
            # not considering any pressure loss yet
            T_cond_o1 = T_desup_o1

            T_cond_i2 = self.inl[1].calc_T()
            h_cond_o2 = self.inl[1].h.val_SI + abs(Q_cond) / self.inl[1].m.val_SI
            T_cond_o2 = T_mix_ph(self.inl[1].p.val_SI, h_cond_o2, self.inl[1].fluid_data)
            T_desup_i2 = T_cond_o2
            T_desup_o2 = self.outl[1].calc_T()

            ttd_desup_u = T_desup_i1 - T_desup_o2
            ttd_desup_l = T_desup_o1 - T_desup_i2

            ttd_cond_u = T_cond_i1 - T_cond_o2
            ttd_cond_l = T_cond_o1 - T_cond_i2

            if ttd_desup_u < 0:
                ttd_desup_u = abs(ttd_desup_u)
            if ttd_desup_l < 0:
                ttd_desup_l = abs(ttd_desup_l)

            if ttd_cond_u < 0:
                ttd_cond_u = abs(ttd_cond_u)
            if ttd_cond_l < 0:
                ttd_cond_l = abs(ttd_cond_l)

            td_log_desup = (ttd_desup_l - ttd_desup_u) / np.log(ttd_desup_l / ttd_desup_u)
            td_log_cond = (ttd_cond_l - ttd_cond_u) / np.log(ttd_cond_l / ttd_cond_u)

            residual = (
                Q_total
                + self.A.val * (Q_desup / Q_total) * self.U_desup.val * td_log_desup
                + self.A.val * (Q_cond / Q_total) * self.U_cond.val * td_log_cond
            )

        else:
            # only condensation is happening
            residual = Q_total + self.A.val * self.U_cond.val * self.calculate_td_log()

        return residual

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
                self.jacobian[k, c.h.J_col] = self.numeric_deriv(f, 'h', c, d=1e-5)

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

            return self.td_pinch.val - T_cond_i1 + T_cond_o2

        else:
            o2 = self.outl[1]
            return self.td_pinch.val - o1.calc_T_sat() + o2.calc_T()

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

        # this should be exported to a function, which can be called from here
        # and from the UA_func.
        h_sat = h_mix_pQ(self.outl[0].p.val_SI, 1, self.outl[0].fluid_data)
        self.Q_desup = self.inl[0].m.val_SI * (h_sat - self.inl[0].h.val_SI)
        Q_desup = self.Q_desup
        self.Q_cond = self.inl[0].m.val_SI * (self.outl[0].h.val_SI - h_sat)
        Q_cond = self.Q_cond
        T_desup_i1 = self.inl[0].calc_T()
        self.T_desup_o1 =  self.outl[0].calc_T_sat()
        T_desup_o1 = self.T_desup_o1
        T_cond_i1 = T_desup_o1
        T_cond_o1 = T_desup_o1

        T_cond_i2 = self.inl[1].calc_T()
        h_cond_o2 = self.inl[1].h.val_SI + abs(Q_cond) / self.inl[1].m.val_SI
        T_cond_o2 = T_mix_ph(self.inl[1].p.val_SI, h_cond_o2, self.inl[1].fluid_data)
        T_desup_i2 = T_cond_o2
        self.T_desup_i2 = T_desup_i2
        T_desup_o2 = self.outl[1].calc_T()

        ttd_desup_u = T_desup_i1 - T_desup_o2
        ttd_desup_l = T_desup_o1 - T_desup_i2

        ttd_cond_u = T_cond_i1 - T_cond_o2
        ttd_cond_l = T_cond_o1 - T_cond_i2

        td_log_desup = (ttd_desup_l - ttd_desup_u) / np.log(ttd_desup_l / ttd_desup_u)
        td_log_cond = (ttd_cond_l - ttd_cond_u) / np.log(ttd_cond_l / ttd_cond_u)

        self.td_pinch.val = T_desup_o1 - T_desup_i2

        self.A.val = abs(self.Q.val) / ((Q_desup / self.Q.val) * self.U_desup.val * td_log_desup + (Q_cond / self.Q.val) * self.U_cond.val * td_log_cond)

        # some intermediate tests for consistency
        assert abs(abs(self.Q.val) / ((Q_desup / self.Q.val) * self.U_desup.val * td_log_desup + (Q_cond / self.Q.val) * self.U_cond.val * td_log_cond) - self.A.val) < 1e-6
        assert round(Q_cond + Q_desup, 3) == round(self.Q.val, 3)
