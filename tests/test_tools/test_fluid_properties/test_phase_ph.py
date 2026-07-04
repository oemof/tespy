# -*- coding: utf-8

"""Unit tests for phase_ph across fluid property wrappers.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_tools/test_fluid_properties/test_phase_ph.py

SPDX-License-Identifier: MIT
"""
import pytest

from tespy.tools.fluid_properties.functions import phase_mix_ph
from tespy.tools.fluid_properties.wrappers import CoolPropWrapper
from tespy.tools.fluid_properties.wrappers import IAPWSWrapper
from tespy.tools.fluid_properties.wrappers import PyromatWrapper

# (p, T, expected_phase) states that every water wrapper must handle identically
WATER_PHASE_CASES = [
    (1e5,   300.0, "l"),   # subcooled liquid  (P < Pc, T < T_sat)
    (1e5,   450.0, "g"),   # superheated steam (P < Pc, T > T_sat)
    (1e6,   700.0, "g"),   # supercritical gas (P < Pc, T > Tc)
    (30e6,  300.0, "l"),   # supercritical liquid (P > Pc, T < Tc)
    (30e6,  700.0, "sc"),  # true supercritical (P > Pc, T > Tc)
]


class TestCoolPropWrapperPhasePh:
    """phase_ph covers all six CoolProp phase regions for water."""

    def setup_method(self):
        self.w = CoolPropWrapper("Water")

    @pytest.mark.parametrize("p,T,expected", WATER_PHASE_CASES)
    def test_phase_from_pT(self, p, T, expected):
        h = self.w.h_pT(p, T)
        assert self.w.phase_ph(p, h) == expected

    def test_two_phase(self):
        h = self.w.h_pQ(1e5, 0.5)
        assert self.w.phase_ph(1e5, h) == "tp"

    def test_saturated_liquid_edge(self):
        h = self.w.h_pQ(1e5, 0.0)
        assert self.w.phase_ph(1e5, h) == "tp"

    def test_saturated_vapor_edge(self):
        h = self.w.h_pQ(1e5, 1.0)
        assert self.w.phase_ph(1e5, h) == "tp"


class TestCoolPropWrapperIncompPhasePh:
    """INCOMP back-end always returns liquid."""

    def setup_method(self):
        self.w = CoolPropWrapper("DowJ", back_end="INCOMP")

    def test_always_liquid(self):
        h = self.w.h_pT(1e5, 300.0)
        assert self.w.phase_ph(1e5, h) == "l"


class TestIAPWSWrapperPhasePh:
    """phase_ph covers the main thermodynamic regions for water via IAPWSWrapper."""

    @pytest.mark.parametrize("back_end", ["IF97", "IF95"])
    @pytest.mark.parametrize("p,T,expected", WATER_PHASE_CASES)
    def test_phase_from_pT(self, back_end, p, T, expected):
        w = IAPWSWrapper("H2O", back_end=back_end)
        h = w.h_pT(p, T)
        assert w.phase_ph(p, h) == expected

    @pytest.mark.parametrize("back_end", ["IF97", "IF95"])
    def test_two_phase(self, back_end):
        w = IAPWSWrapper("H2O", back_end=back_end)
        h = w.h_pQ(1e5, 0.5)
        assert w.phase_ph(1e5, h) == "tp"


class TestPyromatWrapperPhasePh:
    """phase_ph for PyromatWrapper: ig always gas, mp covers all water regions."""

    def test_ig_always_gas(self):
        w = PyromatWrapper("N2", back_end="ig")
        h = w.h_pT(1e5, 400.0)
        assert w.phase_ph(1e5, h) == "g"

    @pytest.mark.parametrize("p,T,expected", WATER_PHASE_CASES)
    def test_mp_phase_from_pT(self, p, T, expected):
        w = PyromatWrapper("H2O", back_end="mp")
        h = w.h_pT(p, T)
        assert w.phase_ph(p, h) == expected

    def test_mp_two_phase(self):
        w = PyromatWrapper("H2O", back_end="mp")
        h = w.h_pQ(1e5, 0.5)
        assert w.phase_ph(1e5, h) == "tp"


# Air composition from the gas turbine tutorial
_AIR_FLUID_DATA = {
    "Ar":  {"wrapper": CoolPropWrapper("Ar"),  "mass_fraction": 0.0129},
    "N2":  {"wrapper": CoolPropWrapper("N2"),  "mass_fraction": 0.7553},
    "CO2": {"wrapper": CoolPropWrapper("CO2"), "mass_fraction": 0.0004},
    "O2":  {"wrapper": CoolPropWrapper("O2"),  "mass_fraction": 0.2314},
}


class TestPhaseMixPhIdealCond:
    """phase_mix_ph returns 'g' for ideal and ideal-cond air mixtures."""

    @pytest.mark.parametrize("mixing_rule", ["ideal", "ideal-cond"])
    def test_gas_at_ambient(self, mixing_rule):
        assert phase_mix_ph(1e5, 3e5, _AIR_FLUID_DATA, mixing_rule) == "g"

    @pytest.mark.parametrize("mixing_rule", ["ideal", "ideal-cond"])
    def test_gas_at_high_temperature(self, mixing_rule):
        assert phase_mix_ph(1e5, 1.4e6, _AIR_FLUID_DATA, mixing_rule) == "g"


class TestPhaseMixPhIncompressible:
    """phase_mix_ph returns 'l' for incompressible water-ethanol mixture."""

    def setup_method(self):
        self.fluid_data = {
            "Water":   {"wrapper": CoolPropWrapper("Water",   back_end="INCOMP"), "mass_fraction": 0.6},
            "Ethanol": {"wrapper": CoolPropWrapper("Ethanol", back_end="INCOMP"), "mass_fraction": 0.4},
        }

    def test_liquid_at_ambient(self):
        assert phase_mix_ph(1e5, 1e5, self.fluid_data, "incompressible") == "l"

    def test_liquid_at_high_pressure(self):
        assert phase_mix_ph(50e6, 1e5, self.fluid_data, "incompressible") == "l"
