# -*- coding: utf-8

"""Tests for the MixingRuleRegistry.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_tools/test_fluid_properties/test_mixing_rule_registry.py

SPDX-License-Identifier: MIT
"""
import pytest

from tespy.tools.fluid_properties import MIXING_RULES
from tespy.tools.fluid_properties import CoolPropWrapper
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.fluid_properties import s_mix_pT
from tespy.tools.fluid_properties import v_mix_pT
from tespy.tools.fluid_properties.mixtures import MixingRuleRegistry

FLUID_DATA = {
    "N2": {"wrapper": CoolPropWrapper("N2"), "mass_fraction": 0.75},
    "O2": {"wrapper": CoolPropWrapper("O2"), "mass_fraction": 0.25},
}


class TestMixingRuleRegistryBuiltins:

    def test_h_pT_registered_for_all_rules(self):
        for rule in ("ideal", "ideal-cond", "incompressible", "forced-gas", "humidair"):
            assert rule in MIXING_RULES._h_pT

    def test_forced_gas_excluded_from_T_ph_inversion(self):
        assert "forced-gas" not in MIXING_RULES._T_ph

    def test_all_standard_rules_support_T_ph_inversion(self):
        for rule in ("ideal", "ideal-cond", "incompressible", "humidair"):
            assert rule in MIXING_RULES._T_ph

    def test_exergy_chemical_only_for_ideal_cond(self):
        assert "ideal-cond" in MIXING_RULES._exergy_chemical
        assert "ideal" not in MIXING_RULES._exergy_chemical


class TestMixingRuleRegistryUnknownRule:

    def test_unknown_rule_h_pT_raises(self):
        with pytest.raises(KeyError, match="not available"):
            h_mix_pT(1e5, 300, FLUID_DATA, mixing_rule="nonexistent")

    def test_unknown_rule_s_pT_raises(self):
        with pytest.raises(KeyError, match="not available"):
            s_mix_pT(1e5, 300, FLUID_DATA, mixing_rule="nonexistent")

    def test_unknown_rule_v_pT_raises(self):
        with pytest.raises(KeyError, match="not available"):
            v_mix_pT(1e5, 300, FLUID_DATA, mixing_rule="nonexistent")

    def test_forced_gas_T_ph_raises(self):
        from tespy.tools.fluid_properties.mixtures import MIXING_RULES as MR
        with pytest.raises(KeyError, match="not available"):
            MR.T_ph("forced-gas")


class TestMixingRuleRegistryCustomRule:

    def setup_method(self):
        self.registry = MixingRuleRegistry()

    def test_register_h_pT(self):
        sentinel = object()
        self.registry.register("my-rule", h_pT=lambda p, T, fd, **kw: sentinel)
        assert self.registry.h_pT("my-rule")(None, None, None) is sentinel

    def test_register_h_pT_also_registers_T_ph_by_default(self):
        func = lambda p, T, fd, **kw: 0.0
        self.registry.register("my-rule", h_pT=func)
        assert self.registry.T_ph("my-rule") is func

    def test_register_h_pT_with_inversion_false(self):
        self.registry.register("no-inv", h_pT=lambda p, T, fd, **kw: 0.0, T_ph_inversion=False)
        assert "no-inv" in self.registry._h_pT
        assert "no-inv" not in self.registry._T_ph

    def test_register_partial(self):
        self.registry.register("partial", h_pT=lambda p, T, fd, **kw: 0.0)
        with pytest.raises(KeyError, match="not available"):
            self.registry.s_pT("partial")

    def test_custom_rule_called_via_h_mix_pT(self):
        def mass_weighted_h(p, T, fluid_data, **kwargs):
            return sum(
                d["wrapper"].h_pT(p, T) * d["mass_fraction"]
                for d in fluid_data.values()
            )

        MIXING_RULES.register("test-mass-weighted", h_pT=mass_weighted_h)
        result = h_mix_pT(1e5, 300, FLUID_DATA, mixing_rule="test-mass-weighted")
        assert isinstance(result, float)
