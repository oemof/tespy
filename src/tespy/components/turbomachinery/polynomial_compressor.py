
# -*- coding: utf-8

"""Module of class PolynomialCompressor.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/turbomachinery/polynomial_compressor.py

SPDX-License-Identifier: MIT
"""
from CoolProp.CoolProp import PropsSI as PSI
from tespy.tools.fluid_properties import isentropic


from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.components.component import Component
from tespy.components.turbomachinery.base import Turbomachine


# the polynomial compressor is a fundamentally different component
class PolynomialCompressor(Turbomachine):

    @staticmethod
    def powerinlets():
        return ["power"]

    def _preprocess(self, num_nw_vars):
        return Component._preprocess(self, num_nw_vars)

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        if len(self.power_inl) > 0:
            constraints["energy_connector_balance"] = dc_cmc(**{
                "func": self.energy_connector_balance_func,
                "dependents": self.energy_connector_dependents,
                "num_eq_sets": 1
            })

        return constraints

    def get_parameters(self):
        return {
            "Q_diss": dc_cp(max_val=0),
            "P": dc_cp(min_val=0),
            "eta_vol": dc_cp(min_val=0, max_val=1),
            "Q_diss_rel": dc_cp(min_val=0, max_val=1),
            "rpm": dc_cp(min_val=0),
            "reference_state": dc_simple(),
            "eta_vol_group": dc_gcp(
                elements=["reference_state", "eta_vol", "rpm"],
                func=self.eta_vol_group_func,
                dependents=self.eta_vol_group_dependents,
                num_eq_sets=1
            ),
            "eta_s": dc_cp(min_val=0, max_val=1),
            "eta_s_group": dc_gcp(
                elements=["eta_s", "Q_diss_rel"],
                func=self.eta_s_group_func,
                dependents=self.eta_s_group_dependents,
                num_eq_sets=1
            )
        }

    # this is a bit different that in other cases, because the power cannot
    # directly be deduced from the change in enthalpy
    def energy_connector_balance_func(self):
        return (
            self.inl[0].m.val_SI
            * (self.outl[0].h.val_SI - self.inl[0].h.val_SI)
            / (1 - self.Q_diss_rel.val)
            - self.power_inl[0].E.val_SI
        )

    def energy_connector_dependents(self):
        return [self.power_inl[0].E, self.inl[0].m, self.inl[0].h, self.outl[0].h]

    def eta_s_group_func(self):
        i = self.inl[0]
        o = self.outl[0]
        h_out_s = isentropic(
            i.p.val_SI,
            i.h.val_SI,
            o.p.val_SI,
            i.fluid_data,
            i.mixing_rule,
            T0=None
        )
        return (
            self.eta_s.val
            * (o.h.val_SI - i.h.val_SI) / (1 - self.Q_diss_rel.val)
            - (h_out_s - i.h.val_SI)
        )

    def eta_s_deriv(self, increment_filter, k, dependents=None):
        r"""
        Partial derivatives for isentropic efficiency.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        from tespy.tools.helpers import _get_dependents
        dependents = dependents["scalars"][0]
        i = self.inl[0]
        o = self.outl[0]
        f = self.eta_s_group_func

        if o.h.is_var and not i.h.is_var:
            self._partial_derivative(o.h, k, self.eta_s.val / (1 - self.Q_diss_rel.val), increment_filter)
            # remove o.h from the dependents
            dependents = dependents.difference(_get_dependents([o.h])[0])

        for dependent in dependents:
            self._partial_derivative(dependent, k, f, increment_filter)

    def eta_s_group_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].p,
            self.outl[0].h,
        ]

    def eta_vol_group_func(self):
        displacement_ref = self.reference_state.val["displacement_ref"]
        displacement = (
            self.reference_state.val["displacement"] * self.rpm.val
            / (displacement_ref * 3600)
        )
        # calculation of mass flow rate
        return (
            self.inl[0].m.val_SI
            - self.eta_vol.val * displacement / self.inl[0].calc_vol()
        )

    def eta_vol_group_dependents(self):
        return [self.inl[0].m, self.inl[0].p, self.inl[0].h, self.rpm]

    def calc_parameters(self):
        i = self.inl[0]
        o = self.outl[0]
        if self.Q_diss_rel.is_set:
            h_2 = (o.h.val_SI - i.h.val_SI * self.Q_diss_rel.val) / (1 - self.Q_diss_rel.val)
            self.P.val = i.m.val_SI * (h_2 - i.h.val_SI)
            self.Q_diss.val = i.m.val_SI * (o.h.val_SI - h_2)
            self.eta_s.val = (
                isentropic(
                    i.p.val_SI,
                    i.h.val_SI,
                    o.p.val_SI,
                    i.fluid_data,
                    i.mixing_rule,
                    T0=None
                ) - i.h.val_SI
            ) / (
                (h_2 - i.h.val_SI)
            )

    @staticmethod
    def calc_etas_from_polynome(fluid, reference_state, polynomes):

        T_evap = reference_state["T_evap"]
        T_cond = reference_state["T_cond"]
        P_comp = calc_EN12900(polynomes["power"], T_evap, T_cond) * 1000
        Q_evap = calc_EN12900(polynomes["heat"], T_evap, T_cond) * 1000

        p_evap = PSI("P", "T", T_evap, "Q", 1, fluid)
        p_cond = PSI("P", "T", T_cond, "Q", 0, fluid)
        T_sh = reference_state["T_sh"]
        T_sc = reference_state["T_sc"]

        if T_sh > 0:
            h_evap_out = PSI("H", "P", p_evap, "T", T_evap + T_sh, fluid)
        else:
            h_evap_out = PSI("H", "P", p_evap, "Q", 1, fluid)

        if T_sc > 0:
            h_cond_out = PSI("H", "P", p_cond, "T", T_cond - T_sc, fluid)
        else:
            h_cond_out = PSI("H", "P", p_cond, "Q", 1, fluid)

        s_comp_in = PSI("S", "P", p_evap, "H", h_evap_out, fluid)
        h_comp_s = PSI("H", "P", p_cond, "S", s_comp_in, fluid)

        dot_m = Q_evap / (h_evap_out - h_cond_out)
        eta_s = dot_m * (h_comp_s - h_evap_out) / P_comp

        rpm_ref = reference_state["rpm_ref"]
        displacement_ref = reference_state["displacement_ref"]
        displacement = (
            reference_state["displacement"] * rpm_ref
            / (displacement_ref * 3600)
        )
        rho_comp_in = PSI("D", "P", p_evap, "H", h_evap_out, fluid)
        eta_vol = dot_m / (rho_comp_in * displacement)

        return {
            "eta_s": eta_s,
            "eta_vol": eta_vol
        }


@staticmethod
def calc_EN12900(c, t_evap, t_cond):
    t_evap = t_evap - 273.15
    t_cond = t_cond - 273.15
    return (
        c[0]
        + c[1] * t_evap + c[2] * t_cond
        + c[3] * t_evap**2 + c[4] * t_evap * t_cond + c[5] * t_cond ** 2
        + c[6] * t_evap**3 + c[7] * t_evap ** 2 * t_cond
        + c[8] * t_evap * t_cond ** 2 + c[9] * t_cond ** 3
    )
