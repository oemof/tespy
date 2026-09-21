# -*- coding: utf-8

"""pyromat wrapper for fluid property calls.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/tools/fluid_properties/wrappers/pyromat.py

SPDX-License-Identifier: MIT
"""

from .base import FluidPropertyWrapper
from .base import wrapper_registry


@wrapper_registry
class PyromatWrapper(FluidPropertyWrapper):

    def __init__(self, fluid, back_end=None, **kwargs) -> None:
        """Wrapper for the Pyromat fluid property library

        Parameters
        ----------
        fluid : str
            Name of the fluid
        back_end : str, optional
            CoolProp back end for the AbstractState object, by default None
        """
        # avoid unnecessary loading time if not used
        try:
            import pyromat as pm
            pm.config['unit_energy'] = "J"
            pm.config['unit_pressure'] = "Pa"
            pm.config['unit_molar'] = "mol"
        except ModuleNotFoundError:
            msg = (
                "To use the pyromat fluid properties you need to install "
                "pyromat."
            )
            raise ModuleNotFoundError(msg)

        super().__init__(fluid, back_end)
        self._create_AS(pm)
        self._set_constants()

    def _create_AS(self, pm):
        self.AS = pm.get(f"{self.back_end}.{self.fluid}")

    def _set_constants(self):
        self._p_min, self._p_max = 100, 1000e5
        self._T_min, self._T_max = self.AS.Tlim()

        if self.back_end == "mp":
            self._T_crit, self._p_crit = self.AS.critical()
        else:
            self._T_crit = self._T_max
            self._p_crit = self._p_max

        self._molar_mass = self.AS.mw()

    def isentropic(self, p_1, h_1, p_2):
        return self.h_ps(p_2, self.s_ph(p_1, h_1))

    def T_ph(self, p, h):
        return self.AS.T(p=p, h=h)[0]

    def T_ps(self, p, s):
        return self.AS.T(p=p, s=s)[0]

    def h_pT(self, p, T):
        return self.AS.h(p=p, T=T)[0]

    def h_ps(self, p, s):
        return self.AS.h(p=p, s=s)[0]

    def d_ph(self, p, h):
        return self.AS.d(p=p, h=h)[0]

    def d_pT(self, p, T):
        return self.AS.d(p=p, T=T)[0]

    def s_ph(self, p, h):
        return self.AS.s(p=p, h=h)[0]

    def s_pT(self, p, T):
        return self.AS.s(p=p, T=T)[0]

    def h_QT(self, Q, T):
        if self.back_end == "ig":
            self._not_implemented()
        return self.AS.h(x=Q, T=T)[0]

    def h_pQ(self, p, Q):
        if self.back_end == "ig":
            self._not_implemented()
        return self.AS.h(p=p, x=Q)[0]

    def s_QT(self, Q, T):
        if self.back_end == "ig":
            self._not_implemented()
        return self.AS.s(x=Q, T=T)[0]

    def Q_ph(self, p, h):
        if self.back_end == "ig":
            self._not_implemented()
        return self.AS.x(p=p, h=h)[0]

    def d_QT(self, Q, T):
        if self.back_end == "ig":
            self._not_implemented()
        return self.AS.d(x=Q, T=T)[0]

    def phase_ph(self, p, h):
        if self.back_end == "ig":
            return "g"
        if p >= self._p_crit:
            if self.T_ph(p, h) >= self._T_crit:
                return "sc"
            else:
                return "l"
        h_bubble = self.h_pQ(p, 0)
        h_dew = self.h_pQ(p, 1)
        if h <= h_bubble:
            return "l"
        elif h >= h_dew:
            return "g"
        else:
            return "tp"
