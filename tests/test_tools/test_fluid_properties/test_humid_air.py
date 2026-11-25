# -*- coding: utf-8

"""Module for testing the humid air property functions.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_tools/test_fluid_properties/test_humid_air.py

SPDX-License-Identifier: MIT
"""
import pytest
from CoolProp.CoolProp import HAPropsSI

from tespy.tools.fluid_properties.mixtures import h_mix_pT_humidair
from tespy.tools.fluid_properties.mixtures import s_mix_pT_humidair
from tespy.tools.fluid_properties.mixtures import v_mix_pT_humidair
from tespy.tools.fluid_properties.mixtures import viscosity_mix_pT_humidair

FORWARD_FUNCTIONS = [
    (h_mix_pT_humidair, "H"),
    (s_mix_pT_humidair, "S"),
    (v_mix_pT_humidair, "V"),
    (viscosity_mix_pT_humidair, "Visc"),
]


def _fluid_data(w):
    return {
        "air": {"mass_fraction": 1 / (1 + w)},
        "water": {"mass_fraction": w / (1 + w)},
    }


@pytest.mark.parametrize("function, key", FORWARD_FUNCTIONS)
def test_water_above_saturation_is_not_part_of_humid_air(function, key):
    p, T = 101325, 293.15
    w_sat = HAPropsSI("W", "P", p, "T", T, "R", 1)
    assert function(p, T, _fluid_data(0.05)) == pytest.approx(
        HAPropsSI(key, "P", p, "T", T, "W", w_sat)
    )


@pytest.mark.parametrize("function, key", FORWARD_FUNCTIONS)
def test_unsaturated_humid_air(function, key):
    p, T, w = 101325, 293.15, 0.005
    assert function(p, T, _fluid_data(w)) == pytest.approx(
        HAPropsSI(key, "P", p, "T", T, "W", w)
    )
