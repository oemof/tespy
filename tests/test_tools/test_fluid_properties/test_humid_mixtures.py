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
