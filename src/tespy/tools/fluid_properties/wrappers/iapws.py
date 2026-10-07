# -*- coding: utf-8

"""iapws wrapper for fluid property calls.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/tools/fluid_properties/wrappers/iapws.py

SPDX-License-Identifier: MIT
"""

from .base import FluidPropertyWrapper
from .base import wrapper_registry


@wrapper_registry
class IAPWSWrapper(FluidPropertyWrapper):


    def __init__(self, fluid, back_end=None, **kwargs) -> None:
        """Wrapper for iapws library calls

        Parameters
        ----------
        fluid : str
            Name of the fluid
        back_end : str, optional
            CoolProp back end for the AbstractState object, by default "IF97"
        """
        # avoid unnecessary loading time if not used
        try:
            import iapws
        except ModuleNotFoundError:
            msg = (
                "To use the iapws fluid properties you need to install "
                "iapws."
            )
            raise ModuleNotFoundError(msg)

        if back_end is None:
            back_end = "IF97"
        super().__init__(fluid, back_end)

        if self.back_end == "IF97":
            self.AS = iapws.IAPWS97
        elif self.back_end == "IF95":
            self.AS = iapws.IAPWS95
        else:
            msg = f"The specified back_end {self.back_end} is not available."
            raise NotImplementedError(msg)
        self._set_constants(iapws)

    def _set_constants(self, iapws):
        self._T_min = iapws._iapws.Tt
        self._T_max = 2000
        self._p_min = iapws._iapws.Pt * 1e6
        self._p_max = 100e6
        self._p_crit = iapws._iapws.Pc * 1e6
        self._T_crit = iapws._iapws.Tc
        self._molar_mass = iapws._iapws.M

    def _is_below_T_critical(self, T):
        return T < self._T_crit

    def isentropic(self, p_1, h_1, p_2):
        return self.h_ps(p_2, self.s_ph(p_1, h_1))

    def T_ph(self, p, h):
        return self.AS(h=h / 1e3, P=p / 1e6).T

    def T_ps(self, p, s):
        return self.AS(s=s / 1e3, P=p / 1e6).T

    def h_pQ(self, p, Q):
        return self.AS(P=p / 1e6, x=Q).h * 1e3

    def h_ps(self, p, s):
        return self.AS(P=p / 1e6, s=s / 1e3).h * 1e3

    def h_pT(self, p, T):
        return self.AS(P=p / 1e6, T=T).h * 1e3

    def h_QT(self, Q, T):
        return self.AS(T=T, x=Q).h * 1e3

    def s_QT(self, Q, T):
        return self.AS(T=T, x=Q).s * 1e3

    def T_sat(self, p):
        return self.AS(P=p / 1e6, x=0).T

    def p_sat(self, T):
        if T > self._T_crit:
            T = self._T_crit * 0.99

        return self.AS(T=T, x=0).P * 1e6

    def Q_ph(self, p, h):
        return self.AS(h=h / 1e3, P=p / 1e6).x

    def phase_ph(self, p, h):
        # (h, P) correctly identifies two-phase but gives unreliable .phase for
        # single-phase in some regions (e.g. IF95 superheated steam). Use it
        # only for the two-phase check; re-evaluate via (T, P) for everything else.
        state_hp = self.AS(h=h / 1e3, P=p / 1e6)
        if state_hp.phase in ["Two phases", "Saturated vapor", "Saturated liquid"]:
            return "tp"
        phase = self.AS(T=state_hp.T, P=p / 1e6).phase
        if phase in ["Liquid", "Compressible liquid"]:
            return "l"
        elif phase in ["Vapour", "Gas"]:
            return "g"
        else:
            return "sc"

    def d_ph(self, p, h):
        return self.AS(h=h / 1e3, P=p / 1e6).rho

    def d_pT(self, p, T):
        return self.AS(T=T, P=p / 1e6).rho

    def d_QT(self, Q, T):
        return self.AS(T=T, x=Q).rho

    def viscosity_ph(self, p, h):
        return self.AS(P=p / 1e6, h=h / 1e3).mu

    def viscosity_pT(self, p, T):
        return self.AS(T=T, P=p / 1e6).mu

    def s_ph(self, p, h):
        return self.AS(P=p / 1e6, h=h / 1e3).s * 1e3

    def s_pT(self, p, T):
        return self.AS(P=p / 1e6, T=T).s * 1e3
