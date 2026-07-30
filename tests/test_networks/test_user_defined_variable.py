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


def test_variable_shared_by_two_equations():
    """Two valves in series share one flow coefficient: neither equation
    alone determines the variable, it couples with the intermediate
    pressure into one square block."""
    nw = Network(iterinfo=False)
    nw.units.set_defaults(temperature="degC", pressure="bar")
    va1 = Valve("valve 1")
    va2 = Valve("valve 2")
    c1 = Connection(Source("source"), "out1", va1, "in1", label="c1")
    c2 = Connection(va1, "out1", va2, "in1", label="c2")
    c3 = Connection(va2, "out1", Sink("sink"), "in1", label="c3")
    nw.add_conns(c1, c2, c3)
    c1.set_attr(fluid={"water": 1}, T=20, p=10, m=5)
    c3.set_attr(p=6)

    coefficient = UserDefinedVariable("C", val0=1e-3)

    def valve_flow(ude):
        c_in, c_out = ude.conns
        dp = c_in.p.val_SI - c_out.p.val_SI
        rho = 1 / c_in.calc_vol()
        return c_in.m.val_SI - ude.params["C"].val_SI * (dp * rho) ** 0.5

    def valve_flow_dependents(ude):
        c_in, c_out = ude.conns
        # the variable object is listed directly, like a connection
        # property container
        return [c_in.m, c_in.p, c_in.h, c_out.p, ude.params["C"]]

    nw.add_udv(coefficient)
    nw.add_ude(UserDefinedEquation(
        "valve 1 flow", valve_flow, valve_flow_dependents,
        conns=[c1, c2], params={"C": coefficient}
    ))
    nw.add_ude(UserDefinedEquation(
        "valve 2 flow", valve_flow, valve_flow_dependents,
        conns=[c2, c3], params={"C": coefficient}
    ))
    nw.solve("design")
    nw.assert_convergence()

    # identical coefficients split the pressure difference about evenly
    assert nw.get_conn("c2").p.val_SI == approx(8e5, rel=1e-3)
    rho = PropsSI("D", "P", 10e5, "T", 293.15, "water")
    assert coefficient.val_SI == approx(5 / (2e5 * rho) ** 0.5, rel=1e-3)

    blocks = nw.problem.decompose().blocks
    udv_blocks = [
        b for b in blocks
        if any("valve" in str(label) for label, _ in b.equation_labels)
    ]
    assert len(udv_blocks) == 1
    assert udv_blocks[0].kind == "square"
    assert len(udv_blocks[0].equations) == 2


def test_two_variables_chained_equations():
    """The first equation determines the first variable, the second
    equation relates the second variable to the first: two chained scalar
    blocks, the relating equation touches no connection at all."""
    nw, kv, kv_ude = _kv_model()
    scaled = UserDefinedVariable("kv scaled", val0=1.0)

    def relation_func(ude):
        return ude.params["scaled"].val_SI - ude.params["kv"].val_SI * 3600

    def relation_dependents(ude):
        return [ude.params["kv"], ude.params["scaled"]]

    relation_ude = UserDefinedEquation(
        "kv scaling", relation_func, relation_dependents,
        params={"kv": kv, "scaled": scaled}
    )
    nw.add_udv(kv, scaled)
    nw.add_ude(kv_ude)
    nw.add_ude(relation_ude)
    nw.solve("design")
    nw.assert_convergence()

    expected = _expected_kv()
    assert kv.val_SI == approx(expected, rel=1e-4)
    assert scaled.val_SI == approx(expected * 3600, rel=1e-4)

    decomposition = nw.problem.decompose()
    block_of = {}
    for block in decomposition.blocks:
        for label, _ in block.equation_labels:
            block_of[label] = block
    kv_block = block_of["kv relation"]
    relation_block = block_of["kv scaling"]
    assert relation_block.kind == "scalar"
    assert kv_block.id in decomposition.precedence[relation_block.id]


def test_get_and_del_udv():
    nw, kv, ude = _kv_model()
    nw.add_udv(kv)
    nw.add_ude(ude)
    assert nw.get_udv("kv") is kv
    nw.del_udv(kv)
    with raises(TESPyNetworkError, match="not part of the network"):
        nw.solve("design")


def test_never_added_variable_raises():
    """An equation referencing a variable that was never registered on the
    network must fail with a descriptive error instead of treating the
    value as a constant silently."""
    nw, kv, ude = _kv_model()
    nw.add_ude(ude)
    with raises(TESPyNetworkError, match="not part of the network"):
        nw.solve("design")


def test_non_variable_property_dependent_raises():
    """Listing a property that can never be a solver variable, e.g. a
    temperature container, must fail with a descriptive error instead of
    silently leaving the equation with incomplete dependencies."""
    nw, kv, _ = _kv_model()

    def faulty_dependents(ude):
        c_in, c_out = ude.conns
        return [c_in.m, c_in.T, c_out.p, ude.params["kv"].variable]

    ude = UserDefinedEquation(
        "kv relation", kv_relation, faulty_dependents,
        conns=[nw.get_conn("c1"), nw.get_conn("c2")], params={"kv": kv}
    )
    nw.add_udv(kv)
    nw.add_ude(ude)
    with raises(TESPyNetworkError, match="never a variable of the solver"):
        nw.solve("design")
