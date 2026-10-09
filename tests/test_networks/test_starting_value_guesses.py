# -*- coding: utf-8

"""Module for testing temperature and quality starting value guesses.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_starting_value_guesses.py

SPDX-License-Identifier: MIT
"""
import numpy as np
from CoolProp.CoolProp import PropsSI
from pytest import approx

from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.networks import Network


def _single_connection_network():
    nw = Network(iterinfo=False)
    nw.units.set_defaults(temperature="degC", pressure="bar")
    c = Connection(Source("source"), "out1", Sink("sink"), "in1", label="c1")
    nw.add_conns(c)
    return nw, c


def test_temperature_guess_generates_enthalpy():
    nw, c = _single_connection_network()
    c.set_attr(fluid={"R290": 1}, p0=5, T0=20)
    nw.solve("design", init_only=True)
    expected = PropsSI("H", "P", 5e5, "T", 293.15, "R290")
    assert c.h.val_SI == approx(expected, rel=1e-3)


def test_quality_guess_generates_enthalpy():
    nw, c = _single_connection_network()
    c.set_attr(fluid={"R290": 1}, p0=5, x0=0.5)
    nw.solve("design", init_only=True)
    expected = PropsSI("H", "P", 5e5, "Q", 0.5, "R290")
    assert c.h.val_SI == approx(expected, rel=1e-3)


def test_quality_and_temperature_guess_generate_pressure():
    nw, c = _single_connection_network()
    c.set_attr(fluid={"R290": 1}, T0=20, x0=0.5)
    nw.solve("design", init_only=True)
    expected = PropsSI("P", "T", 293.15, "Q", 0.5, "R290")
    assert c.p.val_SI == approx(expected, rel=1e-3)
    expected_h = PropsSI("H", "P", expected, "Q", 0.5, "R290")
    assert c.h.val_SI == approx(expected_h, rel=1e-3)


def test_quality_guess_at_supercritical_pressure_is_dropped():
    nw, c = _single_connection_network()
    c.set_attr(fluid={"R290": 1}, p=60, x0=0.5)
    nw.solve("design", init_only=True)
    assert c.p.val_SI == approx(60e5)

    nw, c = _single_connection_network()
    c.set_attr(fluid={"R290": 1}, p0=60, x0=0.5)
    nw.solve("design", init_only=True)
    assert c.p.val_SI == approx(60e5)


def _solved_heater_network():
    nw = Network(iterinfo=False)
    nw.units.set_defaults(temperature="degC", pressure="bar")
    she = SimpleHeatExchanger("heater")
    c1 = Connection(Source("source"), "out1", she, "in1", label="c1")
    c2 = Connection(she, "out1", Sink("sink"), "in1", label="c2")
    nw.add_conns(c1, c2)
    c1.set_attr(fluid={"R290": 1}, m=1, p=5, T=20)
    she.set_attr(pr=1, Q=1e4)
    nw.solve("design")
    nw.assert_convergence()
    return nw, c2


def test_guess_overrides_previous_solution():
    nw, c2 = _solved_heater_network()
    c2.set_attr(x0=0.5)
    nw.solve("design", init_only=True)
    expected = PropsSI("H", "P", 5e5, "Q", 0.5, "R290")
    assert c2.h.val_SI == approx(expected, rel=1e-3)


def test_guess_is_cleared_by_successful_solve():
    nw, c2 = _solved_heater_network()
    c2.set_attr(x0=0.5)
    nw.solve("design")
    nw.assert_convergence()
    h2 = c2.h.val_SI
    assert np.isnan(c2.x.val0)

    nw.solve("design")
    assert c2.h.val_SI == approx(h2)


def test_temperature_guess_beats_generic_values_in_solve():
    nw, c = _single_connection_network()
    c.set_attr(fluid={"R290": 1}, m=1, p=5, T0=20)
    nw.presolve("design", check=False)
    expected = PropsSI("H", "P", 5e5, "T", 293.15, "R290")
    values = {
        (label, prop): value
        for _, represented in nw.get_variable_values().items()
        for label, prop, value in represented
    }
    assert values[("c1", "h")] == approx(expected, rel=1e-3)
