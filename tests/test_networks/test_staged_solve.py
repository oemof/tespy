# -*- coding: utf-8

"""Module for testing the staged solve API.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_staged_solve.py

SPDX-License-Identifier: MIT
"""
from pytest import approx
from pytest import raises

from tespy.components import Compressor
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.networks import Network
from tespy.tools.helpers import TESPyNetworkError


def _compressor_network():
    nw = Network(iterinfo=False)
    so = Source("source")
    compressor = Compressor("compressor")
    si = Sink("sink")
    c1 = Connection(so, "out1", compressor, "in1", label="c1")
    c2 = Connection(compressor, "out1", si, "in1", label="c2")
    nw.add_conns(c1, c2)
    c1.set_attr(fluid={"air": 1}, m=1, p=1e5, T=293.15)
    compressor.set_attr(pr=3, eta_s=0.8)
    return nw, c2


def test_staged_solve_matches_one_shot():
    reference, _ = _compressor_network()
    reference.solve("design")
    reference.assert_convergence()

    nw, c2 = _compressor_network()
    nw.presolve("design")
    nw.solve_continue()
    nw.assert_convergence()

    assert c2.h.val_SI == approx(
        reference.get_conn("c2").h.val_SI, rel=1e-6
    )


def test_variable_modification_between_stages():
    nw, c2 = _compressor_network()
    nw.presolve("design")

    # the outlet enthalpy is the only remaining variable, modify its
    # starting value through the loaded variable space
    variables = nw.problem.variables_dict
    assert len(variables) == 1
    container = variables[0]["obj"]
    container.val_SI = 6e5

    nw.solve_continue()
    nw.assert_convergence()
    # the solver must have moved the variable to the solution
    assert c2.h.val_SI != approx(6e5)


def test_solve_continue_requires_presolve():
    nw, _ = _compressor_network()
    with raises(TESPyNetworkError):
        nw.solve_continue()


def test_solve_continue_consumes_preparation():
    nw, _ = _compressor_network()
    nw.presolve("design")
    nw.solve_continue()
    nw.assert_convergence()
    with raises(TESPyNetworkError):
        nw.solve_continue()


def test_solve_continue_after_init_only():
    nw, _ = _compressor_network()
    nw.solve("design", init_only=True)
    nw.solve_continue()
    nw.assert_convergence()


def _open_mass_flow_network():
    nw = Network(iterinfo=False)
    so = Source("source")
    compressor = Compressor("compressor")
    si = Sink("sink")
    c1 = Connection(so, "out1", compressor, "in1", label="c1")
    c2 = Connection(compressor, "out1", si, "in1", label="c2")
    nw.add_conns(c1, c2)
    # mass flow is not specified: it remains a variable representing the
    # mass flow of both connections, solved by the power specification
    c1.set_attr(fluid={"air": 1}, p=1e5, T=293.15)
    c2.set_attr(T=450)
    compressor.set_attr(pr=3, P=2e5)
    return nw, c1, c2


def test_get_variable_values_lists_linear_dependents():
    nw, c1, c2 = _open_mass_flow_network()
    nw.presolve("design")

    values = nw.get_variable_values()
    mass_flow = [
        represents for (_, prop), represents in values.items() if prop == "m"
    ]
    assert len(mass_flow) == 1
    represents = mass_flow[0]
    assert {(lbl, prop) for lbl, prop, _ in represents} == {
        ("c1", "m"), ("c2", "m")
    }
    values = {value for _, _, value in represents}
    assert len(values) == 1

    nw.print_variable_values()


def test_set_variable_value_propagates_to_dependents():
    nw, c1, c2 = _open_mass_flow_network()
    nw.presolve("design")

    nw.set_variable_value("c2", "m", 5)
    assert c1.m.val_SI == approx(5)
    assert c2.m.val_SI == approx(5)

    nw.solve_continue()
    nw.assert_convergence()


def test_set_variable_value_rejects_non_variables():
    nw, c1, c2 = _open_mass_flow_network()
    nw.presolve("design")

    with raises(TESPyNetworkError):
        nw.set_variable_value("c1", "p", 2e5)
    with raises(KeyError):
        nw.set_variable_value("does not exist", "m", 1)
    with raises(KeyError):
        nw.set_variable_value("compressor", "igva", 1)


def test_presolve_checks_determination():
    nw, c2 = _compressor_network()
    c2.set_attr(T=450)
    with raises(TESPyNetworkError):
        nw.presolve("design")
    # inspection of the ill determined problem is possible without check
    nw.solve("design", init_only=True)
    assert len(nw.problem.variables_dict) == 0
