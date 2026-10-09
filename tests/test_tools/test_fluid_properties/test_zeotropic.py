# -*- coding: utf-8

"""Tests for zeotropic mixture fluid property support.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_tools/test_fluid_properties/test_zeotropic.py

SPDX-License-Identifier: MIT
"""
from CoolProp.CoolProp import get_global_param_string
from pytest import approx
from pytest import mark

from tespy.tools.fluid_properties.functions import T_bubble_p
from tespy.tools.fluid_properties.functions import T_dew_p
from tespy.tools.fluid_properties.functions import p_sat_TQ
from tespy.tools.fluid_properties.wrappers import CoolPropWrapper

# TESPy fluid name format — NOT the raw CoolProp name "ISOBUTAN&IPENTANE".
# CoolPropWrapper._identify_mixture() parses the fractions and |mass suffix
# and calls set_mass_fractions() on the AbstractState internally.
ZEOTROPIC_FLUID = "ISOBUTAN[0.5]&IPENTANE[0.5]|mass"
PURE_FLUID = "R134A"
BACKEND = "REFPROP"

skipif_no_refprop = mark.skipif(
    get_global_param_string("REFPROP_version") == "n/a",
    reason="This test requires REFPROP.",
)


@skipif_no_refprop
class TestZeotropicWrapperProperties:

    def setup_method(self):
        self.wrapper = CoolPropWrapper(ZEOTROPIC_FLUID, BACKEND)
        self.fluid_data = {ZEOTROPIC_FLUID: {"wrapper": self.wrapper, "mass_fraction": 1.0}}
        self.pure_wrapper = CoolPropWrapper(PURE_FLUID, BACKEND)
        self.pure_fluid_data = {PURE_FLUID: {"wrapper": self.pure_wrapper, "mass_fraction": 1.0}}
        self.T = 300.0  # K, well inside the two-phase region
        self.p = 5e5   # Pa

    def test_p_sat_TQ_pure_fluid_all_Q_equal(self):
        """For a pure fluid p_sat_TQ is independent of Q."""
        p0 = p_sat_TQ(self.T, 0, self.pure_fluid_data)
        p05 = p_sat_TQ(self.T, 0.5, self.pure_fluid_data)
        p1 = p_sat_TQ(self.T, 1, self.pure_fluid_data)
        assert p0 == approx(p05)
        assert p0 == approx(p1)

    def test_p_sat_TQ_zeotropic_bubble_dew_differ(self):
        """Bubble and dew pressures differ for a zeotropic mixture."""
        p_bubble = p_sat_TQ(self.T, 0, self.fluid_data)
        p_dew = p_sat_TQ(self.T, 1, self.fluid_data)
        assert p_bubble != approx(p_dew, rel=1e-3)

    def test_p_sat_TQ_zeotropic_intermediate_x_between_bubble_and_dew(self):
        """Intermediate quality gives pressure strictly between bubble and dew."""
        p_bubble = p_sat_TQ(self.T, 0, self.fluid_data)
        p_dew = p_sat_TQ(self.T, 1, self.fluid_data)
        p_mid = p_sat_TQ(self.T, 0.5, self.fluid_data)
        p_lo, p_hi = sorted([p_bubble, p_dew])
        assert p_lo < p_mid < p_hi

    def test_T_dew_bubble_pure_equal(self):
        """Dew and bubble temperatures coincide for a pure fluid."""
        assert T_dew_p(self.p, self.pure_fluid_data) == approx(T_bubble_p(self.p, self.pure_fluid_data))

    def test_T_dew_bubble_zeotropic_glide_present(self):
        """Dew temperature exceeds bubble temperature (temperature glide)."""
        assert T_dew_p(self.p, self.fluid_data) > T_bubble_p(self.p, self.fluid_data)

    def test_phase_ph_liquid(self):
        h_bubble = self.wrapper.h_pQ(self.p, 0)
        assert self.wrapper.phase_ph(self.p, h_bubble - 5000) == "l"

    def test_phase_ph_two_phase(self):
        h_bubble = self.wrapper.h_pQ(self.p, 0)
        h_dew = self.wrapper.h_pQ(self.p, 1)
        assert self.wrapper.phase_ph(self.p, (h_bubble + h_dew) / 2) == "tp"

    def test_phase_ph_gas(self):
        h_dew = self.wrapper.h_pQ(self.p, 1)
        assert self.wrapper.phase_ph(self.p, h_dew + 5000) == "g"
