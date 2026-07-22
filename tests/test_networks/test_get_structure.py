# -*- coding: utf-8

"""Module for testing the serializable structure description.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_get_structure.py

SPDX-License-Identifier: MIT
"""
import json
import warnings

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
    return nw


def _get_structure(nw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        return nw.get_structure()


def test_get_structure_requires_preprocessing():
    nw = _compressor_network()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        with raises(TESPyNetworkError):
            nw.get_structure()


def test_get_structure_warns_about_instability():
    nw = _compressor_network()
    nw.solve("design", init_only=True)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        nw.get_structure()
    assert any(
        issubclass(warning.category, FutureWarning) for warning in caught
    )


def test_get_structure_is_json_serializable():
    nw = _compressor_network()
    nw.solve("design", init_only=True)
    structure = _get_structure(nw)
    assert json.loads(json.dumps(structure)) is not None


def test_get_structure_joins_with_topology_identifiers():
    nw = _compressor_network()
    nw.solve("design", init_only=True)
    structure = _get_structure(nw)

    component_labels = {c["label"] for c in structure["components"]}
    connection_labels = {c["label"] for c in structure["connections"]}
    object_labels = component_labels | connection_labels

    assert component_labels == {"source", "compressor", "sink"}
    assert connection_labels == {"c1", "c2"}
    for variable in structure["variables"]:
        assert variable["object"] in object_labels
    for equation in structure["equations"]:
        assert equation["object"] in object_labels

    for connection in structure["connections"]:
        assert connection["source"] in component_labels
        assert connection["target"] in component_labels


def test_get_structure_states_and_solver_space():
    nw = _compressor_network()
    nw.solve("design", init_only=True)
    structure = _get_structure(nw)

    by_key = {
        (v["object"], v["property"]): v for v in structure["variables"]
    }
    # user specified values
    assert by_key[("c1", "m")]["state"] == "specified"
    assert by_key[("c1", "p")]["state"] == "specified"
    # h of c1 presolved from p and T, p of c2 presolved through pr
    assert by_key[("c1", "h")]["state"] == "presolved"
    assert by_key[("c2", "p")]["state"] == "presolved"
    assert by_key[("c2", "m")]["state"] == "presolved"
    # outlet enthalpy is solved by the eta_s equation
    assert by_key[("c2", "h")]["state"] == "variable"

    # the affine relation of the pr link
    p2 = by_key[("c2", "p")]
    assert p2["reference"] != p2["id"] or p2["factor"] != 1.0

    # scalar solver indices and fluid solver indices cover the variable space
    solver_indices = set()
    for variable in structure["variables"]:
        if variable["solver_index"] is None:
            continue
        if isinstance(variable["solver_index"], dict):
            solver_indices |= set(variable["solver_index"].values())
        else:
            solver_indices.add(variable["solver_index"])
    assert solver_indices == set(range(nw.problem.variable_counter))


def test_get_structure_equation_kinds_and_origins():
    nw = _compressor_network()
    nw.solve("design", init_only=True)
    structure = _get_structure(nw)

    by_key = {(e["object"], e["name"]): e for e in structure["equations"]}

    pr = by_key[("compressor", "pr")]
    assert pr["kind"] == "affine"
    assert pr["origin"] == "specification"
    assert pr["state"] == "consumed"

    mass_flow = by_key[("compressor", "mass_flow_constraints")]
    assert mass_flow["kind"] == "affine"
    assert mass_flow["origin"] == "topology"
    assert mass_flow["state"] == "consumed"

    eta_s = by_key[("compressor", "eta_s")]
    assert eta_s["kind"] == "nonlinear"
    assert eta_s["origin"] == "specification"
    assert eta_s["state"] == "active"
    assert len(eta_s["solver_indices"]) == 1
    assert len(eta_s["variables"]) > 0

    # T of c1 was consumed to presolve h
    temperature = by_key[("c1", "T")]
    assert temperature["state"] == "consumed"
    assert temperature["kind"] == "nonlinear"
