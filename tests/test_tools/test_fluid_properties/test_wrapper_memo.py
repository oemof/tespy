# -*- coding: utf-8

"""Tests of the per iteration property memo of the CoolProp wrapper.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_tools/test_fluid_properties/test_wrapper_memo.py

SPDX-License-Identifier: MIT
"""
import numpy as np
import pytest

from tespy.components import Compressor
from tespy.components import CycleCloser
from tespy.components import SimpleHeatExchanger
from tespy.components import Valve
from tespy.connections import Connection
from tespy.networks import Network
from tespy.tools.fluid_properties.wrappers import CoolPropWrapper
from tespy.tools.fluid_properties.wrappers.base import cache_generation
from tespy.tools.fluid_properties.wrappers.base import next_cache_generation
from tespy.tools.fluid_properties.wrappers.base import property_memo

STATES = [(p, h) for p in (2e5, 8e5, 20e5) for h in (2.5e5, 4.5e5, 6.5e5)]


def _count_updates(monkeypatch):
    calls = {"n": 0}
    original = CoolPropWrapper._update

    def counting(self, input_pair, a, b):
        calls["n"] += 1
        return original(self, input_pair, a, b)

    monkeypatch.setattr(CoolPropWrapper, "_update", counting)
    return calls


def test_memo_values_match_direct_evaluation(monkeypatch):
    wrapper = CoolPropWrapper("R290")
    methods = ["T_ph", "s_ph", "d_ph", "Q_ph", "phase_ph"]

    monkeypatch.setitem(property_memo, "enabled", False)
    direct = [[getattr(wrapper, m)(p, h) for m in methods] for p, h in STATES]

    monkeypatch.setitem(property_memo, "enabled", True)
    next_cache_generation()
    memo = [[getattr(wrapper, m)(p, h) for m in methods] for p, h in STATES]
    # a second pass is served from the memo
    repeated = [[getattr(wrapper, m)(p, h) for m in methods] for p, h in STATES]

    assert memo == direct
    assert repeated == direct


def test_memo_answers_repeats_without_flash(monkeypatch):
    wrapper = CoolPropWrapper("R290")
    next_cache_generation()
    calls = _count_updates(monkeypatch)

    wrapper.T_ph(10e5, 5e5)
    wrapper.s_ph(10e5, 5e5)
    wrapper.d_ph(10e5, 5e5)
    wrapper.T_ph(10e5, 5e5)

    assert calls["n"] == 1


def test_memo_cleared_with_the_generation(monkeypatch):
    wrapper = CoolPropWrapper("R290")
    next_cache_generation()
    wrapper.T_ph(10e5, 5e5)
    wrapper.T_ph(12e5, 5e5)
    assert len(wrapper._memo["data"]) == 2

    next_cache_generation()
    calls = _count_updates(monkeypatch)
    wrapper.T_ph(10e5, 5e5)

    assert calls["n"] == 1
    assert len(wrapper._memo["data"]) == 1
    assert wrapper._memo["generation"] == cache_generation.value


def test_memo_shared_between_instances_of_a_fluid(monkeypatch):
    first = CoolPropWrapper("water")
    second = CoolPropWrapper("water")
    other = CoolPropWrapper("R290")
    assert first._memo is second._memo
    assert first._memo is not other._memo

    next_cache_generation()
    calls = _count_updates(monkeypatch)
    first.T_ph(1e5, 3e5)
    second.T_ph(1e5, 3e5)

    assert calls["n"] == 1


def test_memo_size_cap_evicts_oldest(monkeypatch):
    monkeypatch.setitem(property_memo, "max_states", 10)
    wrapper = CoolPropWrapper("R290")
    next_cache_generation()
    enthalpies = np.linspace(3e5, 6e5, 12)
    for h in enthalpies:
        wrapper.T_ph(10e5, h)

    keys = list(wrapper._memo["data"])
    assert len(keys) <= 10
    # the most recent states stay, the oldest ones are gone
    assert keys[-1][1] == enthalpies[-1]
    assert all(key[1] != enthalpies[0] for key in keys)
    # values stay correct after the memo was cleared by the cap
    monkeypatch.setitem(property_memo, "enabled", False)
    direct = wrapper.T_ph(10e5, 3e5)
    monkeypatch.setitem(property_memo, "enabled", True)
    assert wrapper.T_ph(10e5, 3e5) == direct


def test_memo_disabled(monkeypatch):
    monkeypatch.setitem(property_memo, "enabled", False)
    wrapper = CoolPropWrapper("R290")
    next_cache_generation()
    wrapper._memo["data"].clear()
    T = wrapper.T_ph(10e5, 5e5)

    assert len(wrapper._memo["data"]) == 0
    assert T > 0


def test_incompressible_back_end_through_the_memo():
    incomp = CoolPropWrapper("Water", "INCOMP")
    next_cache_generation()
    assert incomp.T_ph(1e5, 1e5) > 273
    assert incomp.phase_ph(1e5, 1e5) == "l"


def _heat_pump():
    nw = Network(iterinfo=False)
    nw.units.set_defaults(temperature="degC", pressure="bar", heat="kW")
    cc = CycleCloser("cycle closer")
    cp = Compressor("compressor")
    cd = SimpleHeatExchanger("condenser")
    va = Valve("valve")
    ev = SimpleHeatExchanger("evaporator")
    c1 = Connection(cc, "out1", cp, "in1", label="1")
    c2 = Connection(cp, "out1", cd, "in1", label="2")
    c3 = Connection(cd, "out1", va, "in1", label="3")
    c4 = Connection(va, "out1", ev, "in1", label="4")
    c0 = Connection(ev, "out1", cc, "in1", label="0")
    nw.add_conns(c1, c2, c3, c4, c0)
    cd.set_attr(pr=0.98, Q=-1000)
    ev.set_attr(pr=0.98)
    cp.set_attr(eta_s=0.85)
    c1.set_attr(T=20, x=1, fluid={"R290": 1})
    c3.set_attr(T=80, x=0)
    return nw


def test_solve_results_independent_of_memo(monkeypatch):
    results = {}
    for enabled in (False, True):
        monkeypatch.setitem(property_memo, "enabled", enabled)
        nw = _heat_pump()
        nw.solve("design")
        nw.assert_convergence()
        results[enabled] = [
            value for c in nw.conns["object"]
            for value in (c.m.val_SI, c.p.val_SI, c.h.val_SI)
        ]

    assert results[True] == pytest.approx(results[False], rel=1e-9)


def test_set_property_memo(monkeypatch):
    from tespy import get_property_memo
    from tespy import set_property_memo
    monkeypatch.setitem(property_memo, "enabled", True)
    monkeypatch.setitem(property_memo, "max_states", 20000)
    set_property_memo(enabled=False, max_states=100)
    assert get_property_memo() == {"enabled": False, "max_states": 100}
    with pytest.raises(ValueError):
        set_property_memo(max_states=1)
