# -*- coding: utf-8

"""Module for testing fluid properties of gas mixtures.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_tools/test_fluid_properties/test_humid_mixtures.py

SPDX-License-Identifier: MIT
"""
from pytest import approx
from scipy.optimize import brentq

from tespy.tools.fluid_properties import CoolPropWrapper
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.fluid_properties import v_mix_pT
from tespy.tools.fluid_properties.mixtures import cond_check


def make_p_sat_and_pp_water_equal(x, fluid_data):

    _, _, _, _, p_sat, pp_water = cond_check(1e5, x, fluid_data, "water")
    return p_sat - pp_water


def test_humid_air_at_saturation():
    fluids_mass_fractions = {
        "air": 0.8,
        "water": 0.2
    }

    fluid_data = {
        k: {"mass_fraction": v, "wrapper": CoolPropWrapper(k, "HEOS")}
        for k, v in fluids_mass_fractions.items()
    }

    temperature = brentq(make_p_sat_and_pp_water_equal, 330, 350, fluid_data)

    assert approx(temperature) == 341.2075
    p = 1e5

    h = h_mix_pT(p, temperature, fluid_data, "ideal-cond")
    temp = T_mix_ph(p, h, fluid_data, "ideal-cond")

    assert approx(temp) == temperature


def test_v_mix_ideal_cond_liquid_fraction_negligible_volume():
    """Specific volume with a small liquid-water fraction must stay close to
    the dry-gas value.

    With x_liq ≈ 1 % liquid water, the liquid volume contribution is tiny
    (ρ_liq ≈ 1000 kg/m³ vs ρ_gas ≈ 1.2 kg/m³), so v_mix should be within a
    few percent of the dry-gas specific volume.

    The previous (buggy) implementation computed
      d = ρ_liq·x_liq + ρ_gas·(1−x_liq),  v = 1/d
    which produced a value ~10× too small because it averaged densities
    instead of specific volumes.
    """
    p = 1e5
    T = 285.0  # K — below dew point so water condenses

    fluids_mass_fractions = {"air": 0.99, "water": 0.01}
    fluid_data = {
        k: {"mass_fraction": v, "wrapper": CoolPropWrapper(k, "HEOS")}
        for k, v in fluids_mass_fractions.items()
    }

    # Reference: dry air specific volume at same T and p
    fluid_data_dry = {
        "air": {"mass_fraction": 1.0, "wrapper": CoolPropWrapper("air", "HEOS")}
    }
    v_dry = v_mix_pT(p, T, fluid_data_dry, "ideal")

    v_cond = v_mix_pT(p, T, fluid_data, "ideal-cond")

    # With only 1 % liquid water the specific volume must be within 5 % of
    # the dry-gas value, not an order of magnitude smaller.
    assert v_cond == approx(v_dry, rel=0.05)
