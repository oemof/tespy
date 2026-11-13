# -*- coding: utf-8

"""Module for testing the SolutionConnection.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_connections/test_solution_connection.py

SPDX-License-Identifier: MIT
"""
import importlib.util

from CoolProp.CoolProp import PropsSI
from pytest import approx
from pytest import fixture
from pytest import mark
from pytest import raises

from tespy.components import Sink
from tespy.components import Source
from tespy.connections import SolutionConnection
from tespy.networks import Network
from tespy.tools.fluid_properties.solutions import base
from tespy.tools.fluid_properties.solutions import get_solution_pair
from tespy.tools.fluid_properties.solutions.libr_water import LiBrWater

needs_teqp = mark.skipif(
    importlib.util.find_spec("teqp") is None, reason="teqp is not installed"
)


@fixture(autouse=True)
def builtin_pairs_only(monkeypatch):
    # user code, e.g. the tutorials, may have registered further pairs with
    # the same key component, which disables the inference from fluid names
    pairs = {
        name: pair for name, pair in base._SOLUTION_PAIRS.items()
        if name in ("libr-water", "nh3-water")
    }
    monkeypatch.setattr(base, "_SOLUTION_PAIRS", pairs)


def _connection(**kwargs):
    return SolutionConnection(
        Source("source"), "out1", Sink("sink"), "in1", **kwargs
    )


def _solve(c):
    nw = Network(iterinfo=False)
    nw.add_conns(c)
    nw.solve("design")
    nw.assert_convergence()
    return nw


def test_pair_inferred_from_fluid_names():
    c = _connection()
    c.set_attr(fluid={"LiBr": 0.5, "H2O": 0.5})
    assert c.mixing_rule == "libr-water"


@needs_teqp
def test_pair_inferred_from_ammonia():
    c = _connection()
    c.set_attr(fluid={"NH3": 0.5, "H2O": 0.5})
    assert c.mixing_rule == "nh3-water"


def test_fallback_to_libr_water_without_any_pair_information():
    c = _connection()
    c.set_attr(xi=0.55, p=1e4, T=330, m=1)
    assert c.mixing_rule is None
    _solve(c)
    assert c.mixing_rule == "libr-water"
    assert c.xi.val_SI == approx(0.55)


def test_conflicting_pair_raises():
    c = _connection()
    c.set_attr(fluid={"LiBr": 0.5, "H2O": 0.5})
    with raises(ValueError, match="cannot be changed"):
        c.set_attr(mixing_rule="nh3-water")


def test_mixing_rule_must_be_a_pair():
    c = _connection()
    with raises(ValueError, match="working pair"):
        c.set_attr(mixing_rule="ideal")


def test_custom_pair_with_same_key_is_not_inferred():
    class LimitedLiBrWater(LiBrWater):
        name = "libr-water-limited"
        w_max = 0.7

    base._SOLUTION_PAIRS[LimitedLiBrWater.name] = LimitedLiBrWater()
    c = _connection()
    c.set_attr(fluid={"LiBr": 0.5, "H2O": 0.5})
    assert c.mixing_rule is None
    c.set_attr(mixing_rule="libr-water-limited")
    assert c.mixing_rule == "libr-water-limited"


@needs_teqp
def test_pure_ammonia_uses_the_reference_state_of_the_pair():
    pair = get_solution_pair("nh3-water")
    c = _connection(mixing_rule="nh3-water")
    c.set_attr(xi=1, T=313.15, x=0, m=1)
    _solve(c)
    assert c.xi.val_SI == approx(1, abs=2e-5)
    assert c.h.val_SI == approx(pair.h_pQ(c.p.val_SI, 0, c.xi.val_SI), abs=1e-3)
    # the reference states of G4-01 and of CoolProp's ammonia differ by a few
    # kJ/kg, the connection must stay on the basis of the pair
    h_coolprop = PropsSI("H", "T", 313.15, "Q", 0, "ammonia")
    assert abs(c.h.val_SI - h_coolprop) > 1e3


@needs_teqp
def test_pressure_from_temperature_and_vapor_mass_fraction():
    # CSL-0249 state 13: evaporator outlet
    c = _connection(mixing_rule="nh3-water")
    c.set_attr(xi=0.9995, T=-8.800141259714 + 273.15, x=0.99, m=1)
    _solve(c)
    assert c.p.val_SI == approx(2.905681496458e+05, rel=1e-6)
    assert c.h.val_SI == approx(1.583403098026e+06, abs=1)


def test_libr_water_pressure_from_temperature_and_vapor_mass_fraction():
    pair = get_solution_pair("libr-water")
    c = _connection(mixing_rule="libr-water")
    c.set_attr(xi=0.55, T=330, x=0.1, m=1)
    _solve(c)
    assert c.p.val_SI == approx(pair.p_bubble(330, 0.55 / 0.9), rel=1e-8)
    assert c.x.val_SI == approx(0.1, abs=1e-8)
