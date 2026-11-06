# -*- coding: utf-8

"""CoolProp wrapper for fluid property calls.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/tools/fluid_properties/wrappers/coolprop.py

SPDX-License-Identifier: MIT
"""
from itertools import islice

import CoolProp as CP
import numpy as np

from .base import FluidPropertyWrapper
from .base import cache_generation
from .base import wrapper_registry

MEMO_ENABLED = True
# the memo is cleared at the start of every solve (up to 8000 states per
# solve for the benchmark models, 40000 for a 100 section heat exchanger
# model), the cap bounds the memory at about 600 bytes per state, the
# oldest states are evicted first
MEMO_MAX_STATES = 20000
_EAGER_OUTPUTS = ("T", "p", "hmass", "smass", "rhomass")
# every connection holds its own wrapper instance, the memo is shared by all
# instances of the same fluid and back end
_MEMOS = {}


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
        self._memo = _MEMOS.setdefault(
            (self.back_end, self.fluid, self.mixture_type, tuple(self.fractions)),
            {"generation": -1, "data": {}}
        )
        self._eager_phase = True

    def _get(self, input_pair, a, b, output):
        """Output of the state, memoized per solve: equations and their
        derivatives request the same states repeatedly within one iteration,
        and the states of connections solved in an earlier block repeat in
        every later block iteration, the memo answers the repeats without a
        flash."""
        if not MEMO_ENABLED:
            self._update(input_pair, a, b)
            return getattr(self.AS, output)()
        memo = self._memo
        if memo["generation"] != cache_generation.value:
            memo["data"].clear()
            memo["generation"] = cache_generation.value
        key = (input_pair, a, b)
        outputs = memo["data"].get(key)
        if outputs is None:
            if len(memo["data"]) >= MEMO_MAX_STATES:
                # the oldest tenth goes (insertion order), one entry at a
                # time would scan the deleted slots at the front of the dict
                for old in list(islice(memo["data"], MEMO_MAX_STATES // 10)):
                    del memo["data"][old]
            outputs = memo["data"][key] = {}
        elif output in outputs:
            return outputs[output]
        self._update(input_pair, a, b)
        # the getters of a set state are free compared to a flash, so the
        # outputs the equations usually ask for are collected right away
        AS = self.AS
        for name in _EAGER_OUTPUTS:
            outputs[name] = getattr(AS, name)()
        # the phase is read along for the back ends that provide it, the
        # quality only for two phase states, it raises otherwise
        if self._eager_phase:
            try:
                phase = outputs["phase"] = AS.phase()
            except ValueError:
                self._eager_phase = False
            else:
                if phase == CP.iphase_twophase:
                    outputs["Q"] = AS.Q()
        if output not in outputs:
            outputs[output] = getattr(AS, output)()
        return outputs[output]

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
        return self._get(CP.HmassP_INPUTS, h, p, "T")

    def T_ps(self, p, s):
        return self._get(CP.PSmass_INPUTS, p, s, "T")

    def h_pQ(self, p, Q):
        return self._get(CP.PQ_INPUTS, p, Q, "hmass")

    def h_ps(self, p, s):
        return self._get(CP.PSmass_INPUTS, p, s, "hmass")

    def h_pT(self, p, T):
        return self._get(CP.PT_INPUTS, p, T, "hmass")

    def h_QT(self, Q, T):
        return self._get(CP.QT_INPUTS, Q, T, "hmass")

    def s_QT(self, Q, T):
        return self._get(CP.QT_INPUTS, Q, T, "smass")

    def T_sat(self, p):
        return self._get(CP.PQ_INPUTS, p, 0, "T")

    def T_dew(self, p):
        return self._get(CP.PQ_INPUTS, p, 1, "T")

    def T_bubble(self, p):
        return self._get(CP.PQ_INPUTS, p, 0, "T")

    def p_sat(self, T):
        return self._get(CP.QT_INPUTS, 0.5, T, "p")

    def p_dew(self, T):
        return self._get(CP.QT_INPUTS, 1, T, "p")

    def p_bubble(self, T):
        return self._get(CP.QT_INPUTS, 0, T, "p")

    def p_sat_TQ(self, T, Q):
        return self._get(CP.QT_INPUTS, Q, T, "p")

    def Q_ph(self, p, h):
        # the state is established first as before, invalid inputs raise here
        self._get(CP.HmassP_INPUTS, h, p, "T")
        if len(self.fractions) > 1:
            return self._get(CP.HmassP_INPUTS, h, p, "Q")

        if p >= self._p_crit:
            return np.nan

        phase = self._get(CP.HmassP_INPUTS, h, p, "phase")
        if phase == CP.iphase_twophase:
            return self._get(CP.HmassP_INPUTS, h, p, "Q")
        elif phase == CP.iphase_liquid:
            return 0
        elif phase in (CP.iphase_gas, CP.iphase_supercritical_gas):
            return 1
        else:
            return np.nan

    def phase_ph(self, p, h):
        if self.back_end == "INCOMP":
            return "l"

        self._get(CP.HmassP_INPUTS, h, p, "T")
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

        phase = self._get(CP.HmassP_INPUTS, h, p, "phase")
        if phase == CP.iphase_twophase:
            return "tp"
        elif phase in (CP.iphase_liquid, CP.iphase_supercritical_liquid):
            return "l"
        elif phase in (CP.iphase_gas, CP.iphase_supercritical_gas):
            return "g"
        else:
            return "sc"

    def d_ph(self, p, h):
        return self._get(CP.HmassP_INPUTS, h, p, "rhomass")

    def d_pT(self, p, T):
        return self._get(CP.PT_INPUTS, p, T, "rhomass")

    def d_QT(self, Q, T):
        return self._get(CP.QT_INPUTS, Q, T, "rhomass")

    def viscosity_ph(self, p, h):
        return self._get(CP.HmassP_INPUTS, h, p, "viscosity")

    def viscosity_pT(self, p, T):
        return self._get(CP.PT_INPUTS, p, T, "viscosity")

    def viscosity_pQ(self, p, Q):
        return self._get(CP.PQ_INPUTS, p, Q, "viscosity")

    def conductivity_ph(self, p, h):
        return self._get(CP.HmassP_INPUTS, h, p, "conductivity")

    def conductivity_pT(self, p, T):
        return self._get(CP.PT_INPUTS, p, T, "conductivity")

    def s_ph(self, p, h):
        return self._get(CP.HmassP_INPUTS, h, p, "smass")

    def s_pT(self, p, T):
        return self._get(CP.PT_INPUTS, p, T, "smass")
