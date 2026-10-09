# -*- coding: utf-8

"""CoolProp wrapper for fluid property calls.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/tools/fluid_properties/wrappers/coolprop.py

SPDX-License-Identifier: MIT
"""
import CoolProp as CP

from .base import FluidPropertyWrapper
from .base import wrapper_registry


class SerializableAbstractState(CP.AbstractState):

    def __init__(self, back_end, fluid_name):
        self.back_end = back_end
        self.fluid_name = fluid_name

    def __reduce__(self):
        return (self.__class__, (self.back_end, self.fluid_name))


@wrapper_registry
class CoolPropWrapper(FluidPropertyWrapper):

    def __init__(self, fluid, back_end=None, **kwargs) -> None:
        """Wrapper for CoolProp.CoolProp.AbstractState instance calls

        Parameters
        ----------
        fluid : str
            Name of the fluid
        back_end : str, optional
            CoolProp back end for the AbstractState object, by default "HEOS"
        """
        super().__init__(fluid, back_end)

        if self.back_end is None:
            self.back_end = "HEOS"

        self._identify_mixture()
        self.AS = SerializableAbstractState(self.back_end, self.fluid)
        self._set_mixture_fractions()
        self._set_constants()
        self._last_ip = None
        self._last_a = None
        self._last_b = None

    def _update(self, input_pair, a, b):
        if input_pair == self._last_ip and a == self._last_a and b == self._last_b:
            return
        self._last_ip = None
        self.AS.update(input_pair, a, b)
        self._last_ip = input_pair
        self._last_a = a
        self._last_b = b

    def _set_mixture_fractions(self):
        """Set the fractions for provided mixture"""
        if self.mixture_type == "mass":
            self.AS.set_mass_fractions(self.fractions)
        elif self.mixture_type == "molar":
            self.AS.set_mole_fractions(self.fractions)
        elif self.mixture_type == "volume":
            self.AS.set_volu_fractions(self.fractions)

    def _set_constants(self):
        """Setup constants for later quick access, e.g. mixture fractions
        minimum/maximum pressure/temperature and critical point properties
        """
        self._T_min = self.AS.trivial_keyed_output(CP.iT_min)
        self._T_max = self.AS.trivial_keyed_output(CP.iT_max)

        if self.back_end == "INCOMP":
            self._p_min = 1e2
            self._p_max = 1e8
            self._p_crit = 1e8
            self._T_crit = None
            self._molar_mass = 1
            if self.mixture_type is not None:
                try:
                    self._T_min = max(
                        self.AS.trivial_keyed_output(CP.iT_freeze),
                        self._T_min
                    )
                except ValueError:
                    pass
        else:
            if self.back_end == "HEOS":
                # see https://github.com/CoolProp/CoolProp/discussions/2443
                self._T_max *= 1.45

            if self.back_end == "REFPROP":
                if self.mixture_type is not None:
                    self._T_min += 5
                self._p_min = 1e1
            else:
                self._p_min = self.AS.trivial_keyed_output(CP.iP_min)
            self._p_max = self.AS.trivial_keyed_output(CP.iP_max)
            self._p_crit = self.AS.trivial_keyed_output(CP.iP_critical)
            self._T_crit = self.AS.trivial_keyed_output(CP.iT_critical)
            self._molar_mass = self.AS.trivial_keyed_output(CP.imolar_mass)

    def _is_below_T_critical(self, T):
        return T < self._T_crit

    def get_T_max(self, p):
        if self.back_end == "INCOMP":
            return self.T_sat(p)
        else:
            return self._T_max

    def isentropic(self, p_1, h_1, p_2):
        return self.h_ps(p_2, self.s_ph(p_1, h_1))

    def T_ph(self, p, h):
        self._update(CP.HmassP_INPUTS, h, p)
        return self.AS.T()

    def T_ps(self, p, s):
        self._update(CP.PSmass_INPUTS, p, s)
        return self.AS.T()

    def h_pQ(self, p, Q):
        self._update(CP.PQ_INPUTS, p, Q)
        return self.AS.hmass()

    def h_ps(self, p, s):
        self._update(CP.PSmass_INPUTS, p, s)
        return self.AS.hmass()

    def h_pT(self, p, T):
        self._update(CP.PT_INPUTS, p, T)
        return self.AS.hmass()

    def h_QT(self, Q, T):
        self._update(CP.QT_INPUTS, Q, T)
        return self.AS.hmass()

    def s_QT(self, Q, T):
        self._update(CP.QT_INPUTS, Q, T)
        return self.AS.smass()

    def T_sat(self, p):
        self._update(CP.PQ_INPUTS, p, 0)
        return self.AS.T()

    def T_dew(self, p):
        self._update(CP.PQ_INPUTS, p, 1)
        return self.AS.T()

    def T_bubble(self, p):
        self._update(CP.PQ_INPUTS, p, 0)
        return self.AS.T()

    def p_sat(self, T):
        self._update(CP.QT_INPUTS, 0.5, T)
        return self.AS.p()

    def p_dew(self, T):
        self._update(CP.QT_INPUTS, 1, T)
        return self.AS.p()

    def p_bubble(self, T):
        self._update(CP.QT_INPUTS, 0, T)
        return self.AS.p()

    def p_sat_TQ(self, T, Q):
        self._update(CP.QT_INPUTS, Q, T)
        return self.AS.p()

    def Q_ph(self, p, h):
        self._update(CP.HmassP_INPUTS, h, p)
        if len(self.fractions) > 1:
            return self.AS.Q()

        phase = self.AS.phase()
        if phase == CP.iphase_twophase:
            return self.AS.Q()
        elif phase == CP.iphase_liquid:
            return 0
        elif phase == CP.iphase_gas:
            return 1
        else:  # all other phases - though this should be unreachable as p is sub-critical
            return -1

    def phase_ph(self, p, h):
        if self.back_end == "INCOMP":
            return "l"

        self._update(CP.HmassP_INPUTS, h, p)
        if self.mixture_type is not None:
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

        phase = self.AS.phase()
        if phase == CP.iphase_twophase:
            return "tp"
        elif phase in (CP.iphase_liquid, CP.iphase_supercritical_liquid):
            return "l"
        elif phase in (CP.iphase_gas, CP.iphase_supercritical_gas):
            return "g"
        else:
            return "sc"

    def d_ph(self, p, h):
        self._update(CP.HmassP_INPUTS, h, p)
        return self.AS.rhomass()

    def d_pT(self, p, T):
        self._update(CP.PT_INPUTS, p, T)
        return self.AS.rhomass()

    def d_QT(self, Q, T):
        self._update(CP.QT_INPUTS, Q, T)
        return self.AS.rhomass()

    def viscosity_ph(self, p, h):
        self._update(CP.HmassP_INPUTS, h, p)
        return self.AS.viscosity()

    def viscosity_pT(self, p, T):
        self._update(CP.PT_INPUTS, p, T)
        return self.AS.viscosity()

    def conductivity_ph(self, p, h):
        self._update(CP.HmassP_INPUTS, h, p)
        return self.AS.conductivity()

    def conductivity_pT(self, p, T):
        self._update(CP.PT_INPUTS, p, T)
        return self.AS.conductivity()

    def s_ph(self, p, h):
        self._update(CP.HmassP_INPUTS, h, p)
        return self.AS.smass()

    def s_pT(self, p, T):
        self._update(CP.PT_INPUTS, p, T)
        return self.AS.smass()
