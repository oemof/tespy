# -*- coding: utf-8

"""Module for testing the UserDefinedVariable feature.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_user_defined_variable.py

SPDX-License-Identifier: MIT
"""
from CoolProp.CoolProp import PropsSI
from pytest import approx
from pytest import raises

from tespy.components import Sink
from tespy.components import Source
from tespy.components import Valve
from tespy.connections import Connection
from tespy.networks import Network
from tespy.tools import UserDefinedEquation
from tespy.tools import UserDefinedVariable
from tespy.tools.helpers import TESPyNetworkError


def kv_relation(ude):
    c_in, c_out = ude.conns
    dp = c_in.p.val_SI - c_out.p.val_SI
    rho = 1 / c_in.calc_vol()
    return c_in.m.val_SI - ude.params["kv"].val_SI * (dp * rho) ** 0.5


def kv_dependents(ude):
    c_in, c_out = ude.conns
    return [c_in.m, c_in.p, c_in.h, c_out.p, ude.params["kv"].variable]


def _kv_model():
    nw = Network(iterinfo=False)
    nw.units.set_defaults(temperature="degC", pressure="bar")
    va = Valve("valve")
    c1 = Connection(Source("source"), "out1", va, "in1", label="c1")
    c2 = Connection(va, "out1", Sink("sink"), "in1", label="c2")
    nw.add_conns(c1, c2)
    c1.set_attr(fluid={"water": 1}, T=20, p=10, m=5)
    c2.set_attr(p=8)

    kv = UserDefinedVariable("kv", val0=1e-3)
    ude = UserDefinedEquation(
        "kv relation", kv_relation, kv_dependents,
        conns=[c1, c2], params={"kv": kv}
    )
    return nw, kv, ude


def _expected_kv():
    rho = PropsSI("D", "P", 10e5, "T", 293.15, "water")
    return 5 / (2e5 * rho) ** 0.5


def test_variable_without_closing_equation_raises():
    nw, kv, _ = _kv_model()
    nw.add_udv(kv)
    with raises(TESPyNetworkError):
        nw.solve("design")


def test_deactivated_variable_with_equation_raises():
    nw, kv, ude = _kv_model()
    nw.add_udv(kv)
    nw.add_ude(ude)
    kv.set_attr(val=1e-3, is_var=False)
    with raises(TESPyNetworkError):
        nw.solve("design")


def test_duplicate_label_raises():
    nw, kv, _ = _kv_model()
    nw.add_udv(kv)
    duplicate = UserDefinedVariable("kv", val0=1.0)
    with raises(Exception):
        nw.add_udv(duplicate)


def test_kv_identification():
    nw, kv, ude = _kv_model()
    nw.add_udv(kv)
    nw.add_ude(ude)
    nw.solve("design")
    nw.assert_convergence()
    assert kv.val_SI == approx(_expected_kv(), rel=1e-4)


def test_kv_identification_simultaneous():
    nw, kv, ude = _kv_model()
    nw.add_udv(kv)
    nw.add_ude(ude)
    nw.solve("design", block_solve=False)
    nw.assert_convergence()
    assert kv.val_SI == approx(_expected_kv(), rel=1e-4)


def test_prediction_with_fixed_variable():
    nw, kv, ude = _kv_model()
    nw.add_udv(kv)
    nw.add_ude(ude)
    nw.solve("design")
    nw.assert_convergence()
    kv_value = kv.val_SI

    kv.set_attr(is_var=False)
    c1 = nw.get_conn("c1")
    c1.set_attr(m=None)
    nw.get_conn("c2").set_attr(p=6)
    nw.solve("design")
    nw.assert_convergence()
    rho = PropsSI("D", "P", 10e5, "T", 293.15, "water")
    expected = kv_value * (4e5 * rho) ** 0.5
    assert c1.m.val_SI == approx(expected, rel=1e-4)


def test_bounds_hold_solution():
    nw, _, ude = _kv_model()
    kv = UserDefinedVariable("kv", val0=1e-3, min_val=1e-6, max_val=1.0)
    ude.params["kv"] = kv
    nw.add_udv(kv)
    nw.add_ude(ude)
    nw.solve("design")
    nw.assert_convergence()
    assert kv.val_SI == approx(_expected_kv(), rel=1e-4)
    assert 1e-6 <= kv.val_SI <= 1.0


def test_bounds_prevent_false_convergence():
    nw, _, ude = _kv_model()
    upper = _expected_kv() / 2
    kv = UserDefinedVariable("kv", val0=1e-4, min_val=0, max_val=upper)
    ude.params["kv"] = kv
    nw.add_udv(kv)
    nw.add_ude(ude)
    nw.solve("design")
    assert nw.status != 0


def test_staged_modification():
    nw, kv, ude = _kv_model()
    nw.add_udv(kv)
    nw.add_ude(ude)
    nw.presolve("design")
    nw.set_variable_value("kv", "variable", 5e-4)
    values = {
        (label, prop): value
        for _, represented in nw.get_variable_values().items()
        for label, prop, value in represented
    }
    assert values[("kv", "variable")] == approx(5e-4)
    nw.solve_continue()
    nw.assert_convergence()
    assert kv.val_SI == approx(_expected_kv(), rel=1e-4)


def test_variable_listed_with_label():
    nw, kv, ude = _kv_model()
    nw.add_udv(kv)
    nw.add_ude(ude)
    nw.presolve("design")
    represented = [
        (label, prop)
        for _, entries in nw.get_variable_values().items()
        for label, prop, _ in entries
    ]
    assert ("kv", "variable") in represented


def test_get_and_del_udv():
    nw, kv, ude = _kv_model()
    nw.add_udv(kv)
    nw.add_ude(ude)
    assert nw.get_udv("kv") is kv
    nw.del_udv(kv)
    with raises(TESPyNetworkError):
        nw.solve("design")
