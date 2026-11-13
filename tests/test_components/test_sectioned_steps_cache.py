# -*- coding: utf-8

"""Tests of the per iteration cache of the phase boundary steps and of the
saturation step search of the sectioned heat exchangers.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_sectioned_steps_cache.py

SPDX-License-Identifier: MIT
"""
import numpy as np
from pytest import approx
from scipy.optimize import brentq

from tespy.components import MovingBoundaryHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.components.heat_exchangers import sectioned
from tespy.components.heat_exchangers.sectioned import _isclose
from tespy.components.heat_exchangers.sectioned import _step_at_saturation
from tespy.components.heat_exchangers.sectioned import identify_step_at_saturation
from tespy.connections import Connection
from tespy.networks import Network
from tespy.tools.fluid_properties.wrappers import CoolPropWrapper
from tespy.tools.fluid_properties.wrappers.base import next_cache_generation


def _condenser():
    nw = Network(iterinfo=False)
    nw.units.set_defaults(temperature="degC", pressure="bar")
    hx = MovingBoundaryHeatExchanger("condenser")
    c1 = Connection(Source("hot in"), "out1", hx, "in1")
    c2 = Connection(hx, "out1", Sink("hot out"), "in1")
    c3 = Connection(Source("cold in"), "out1", hx, "in2")
    c4 = Connection(hx, "out2", Sink("cold out"), "in1")
    nw.add_conns(c1, c2, c3, c4)
    c1.set_attr(fluid={"NH3": 1}, m=1, td_dew=60, T=120)
    c2.set_attr(td_bubble=5)
    c3.set_attr(fluid={"water": 1}, p=1, T=50)
    c4.set_attr(T=60)
    hx.set_attr(dp1=0.5, dp2=0.1)
    nw.solve("design")
    nw.assert_convergence()
    return nw, hx


def test_moving_steps_cached_within_a_generation(monkeypatch):
    _, hx = _condenser()
    calls = {"n": 0}
    original = sectioned.SectionedHeatExchanger._get_moving_steps

    def counting(c1, c2):
        calls["n"] += 1
        return original(c1, c2)

    monkeypatch.setattr(
        sectioned.SectionedHeatExchanger, "_get_moving_steps",
        staticmethod(counting)
    )
    next_cache_generation()
    first = hx._moving_steps(0)
    hx._moving_steps(0)
    hx._moving_steps(1)
    hx._moving_steps(1)
    assert calls["n"] == 2
    assert hx._moving_steps(0) == first

    # a changed state of one side only recomputes that side
    hx.inl[0].h.val_SI += 1
    hx._moving_steps(0)
    hx._moving_steps(1)
    assert calls["n"] == 3

    # a new iteration starts from an empty cache
    next_cache_generation()
    hx._moving_steps(1)
    assert calls["n"] == 4


def test_step_at_saturation_matches_brentq():
    wrapper = CoolPropWrapper("R290")
    fluid_data = {"R290": {"wrapper": wrapper, "mass_fraction": 1.0}}
    p_in = 20e5
    h_in = wrapper.h_pQ(p_in, 0) - 5e4
    delta_h = 4e5
    for delta_p in (-1.5e5, -2e4, 5e4):
        for Q in (0, 1):
            h_sat = wrapper.h_pQ(p_in, Q)
            args = (p_in, h_in, delta_p, delta_h, Q, fluid_data)
            x0 = (h_sat - h_in) / delta_h
            x = _step_at_saturation(0, 1, x0, *args)
            reference = brentq(
                identify_step_at_saturation, 0, 1, args=args, xtol=1e-14
            )
            assert x == approx(reference, abs=1e-9)
            assert 0 <= x <= 1


def test_isclose_matches_numpy():
    for a, b in ((1.0, 1.0 + 1e-6), (1.0, 1.0 + 1e-4), (0.0, 1e-9), (3e5, 3e5 + 10)):
        assert _isclose(a, b) == bool(np.isclose(a, b))
