# -*- coding: utf-8

"""Module for testing the structure graph of the equation system.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_structure_graph.py

SPDX-License-Identifier: MIT
"""
from tespy.components import Compressor
from tespy.components import Merge
from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Splitter
from tespy.connections import Connection
from tespy.connections import Ref
from tespy.networks import Network


def _linear_network():
    nw = Network(iterinfo=False)
    so = Source("source")
    compressor = Compressor("compressor")
    si = Sink("sink")
    c1 = Connection(so, "out1", compressor, "in1", label="c1")
    c2 = Connection(compressor, "out1", si, "in1", label="c2")
    nw.add_conns(c1, c2)
    c1.set_attr(fluid={"air": 1}, m=1, p=1e5, T=293.15)
    compressor.set_attr(pr=3, eta_s=0.8)
    return nw, c1, c2, compressor


def _splitter_merge_network():
    nw = Network(iterinfo=False)
    so = Source("source")
    splitter = Splitter("splitter")
    heater1 = SimpleHeatExchanger("heater 1")
    heater2 = SimpleHeatExchanger("heater 2")
    merge = Merge("merge")
    si = Sink("sink")
    c1 = Connection(so, "out1", splitter, "in1", label="c1")
    c2 = Connection(splitter, "out1", heater1, "in1", label="c2")
    c3 = Connection(heater1, "out1", merge, "in1", label="c3")
    c4 = Connection(splitter, "out2", heater2, "in1", label="c4")
    c5 = Connection(heater2, "out1", merge, "in2", label="c5")
    c6 = Connection(merge, "out1", si, "in1", label="c6")
    nw.add_conns(c1, c2, c3, c4, c5, c6)
    c1.set_attr(fluid={"air": 1}, m=1, p=1e5, T=293.15)
    c2.set_attr(m=0.4)
    # no pr on heater 2: a second affine pressure path from splitter to merge
    # would form a cycle
    heater1.set_attr(pr=0.98, Q=1e4)
    heater2.set_attr(Q=2e4)
    return nw


def test_affine_edges_and_no_linear_rows_on_linear_branch():
    nw, c1, c2, compressor = _linear_network()
    nw.solve("design", init_only=True)
    graph = nw.problem.structure_graph

    # mass flow equality and pr are 2-entry rows, the compressor has no
    # splitting or merging, so no n-ary rows exist
    assert len(graph.linear_rows) == 0
    assert len(graph.edges_with_factors) > 0

    m_cols = (
        nw.problem._object_to_variable_lookup[c1]["m"],
        nw.problem._object_to_variable_lookup[c2]["m"],
    )
    p_cols = (
        nw.problem._object_to_variable_lookup[c1]["p"],
        nw.problem._object_to_variable_lookup[c2]["p"],
    )
    edges = {(col1, col2): factor for col1, col2, factor in graph.edges_with_factors}
    assert tuple(sorted(m_cols)) in edges
    assert tuple(sorted(p_cols)) in edges
    # the pr edge carries the pressure ratio as factor
    assert abs(edges[tuple(sorted(p_cols))]) == 3

    assert graph.find_cycle() is None


def test_nary_balances_appear_in_nonlinear_incidence():
    nw = _splitter_merge_network()
    nw.solve("design", init_only=True)
    graph = nw.problem.structure_graph

    # splitter and merge mass balances are declared as residual functions,
    # not as structure matrix rows, so they are part of the solver incidence
    # with three or more structural variables each
    assert len(graph.linear_rows) == 0
    nary_equations = [
        cols for cols in graph.nonlinear_incidence.values() if len(cols) > 2
    ]
    assert len(nary_equations) >= 2
    # rows recorded as n-ary must not be consumed as affine edges
    assert set(graph.linear_rows) & set(graph.edge_eq_idx.values()) == set()


def test_nonlinear_incidence_attached_after_prepare():
    nw, c1, c2, compressor = _linear_network()
    nw.solve("design", init_only=True)
    graph = nw.problem.structure_graph

    assert len(graph.nonlinear_incidence) == len(
        nw.problem._incidence_matrix
    )
    # incidence refers to structural variable numbers, all of which must be
    # known to the variable lookup
    for cols in graph.nonlinear_incidence.values():
        for col in cols:
            assert col in nw.problem._variable_lookup


def test_equation_origin_classification():
    nw, c1, c2, compressor = _linear_network()
    nw.solve("design", init_only=True)
    graph = nw.problem.structure_graph
    lookup = nw.problem._object_to_variable_lookup

    m_edge = tuple(sorted((lookup[c1]["m"], lookup[c2]["m"])))
    p_edge = tuple(sorted((lookup[c1]["p"], lookup[c2]["p"])))

    origins = nw.problem._equation_set_origin
    # mass flow equality is a mandatory constraint, pr is user imposed
    assert origins[graph.edge_eq_idx[m_edge]] == "topology"
    assert origins[graph.edge_eq_idx[p_edge]] == "specification"


def test_mass_flow_ref_edge_is_specification_origin():
    nw = Network(iterinfo=False)
    so_a = Source("source a")
    si_a = Sink("sink a")
    so_b = Source("source b")
    si_b = Sink("sink b")
    ca = Connection(so_a, "out1", si_a, "in1", label="ca")
    cb = Connection(so_b, "out1", si_b, "in1", label="cb")
    nw.add_conns(ca, cb)
    ca.set_attr(fluid={"air": 1}, m=1, p=1e5, T=293.15)
    cb.set_attr(fluid={"water": 1}, m=Ref(ca, 0.5, 0), p=1e5, T=293.15)
    nw.solve("design", init_only=True)

    # the referenced mass flow links both variables in the affine graph,
    # but as a specification: consumers distinguishing physical topology
    # from user imposed links rely on this classification
    lookup = nw.problem._object_to_variable_lookup
    ref_edge = tuple(sorted((lookup[ca]["m"], lookup[cb]["m"])))
    graph = nw.problem.structure_graph
    assert ref_edge in graph.edge_eq_idx
    origin = nw.problem._equation_set_origin[graph.edge_eq_idx[ref_edge]]
    assert origin == "specification"


def test_structural_key_stable_and_spec_sensitive():
    nw, c1, c2, compressor = _linear_network()
    nw.solve("design", init_only=True)
    key_first = nw.problem.structure_graph.structural_key()

    nw.solve("design", init_only=True)
    key_second = nw.problem.structure_graph.structural_key()
    assert key_first == key_second

    # changing a specified value must not change the key
    compressor.set_attr(eta_s=0.85)
    nw.solve("design", init_only=True)
    assert nw.problem.structure_graph.structural_key() == key_first

    # moving a specification changes the key
    compressor.set_attr(eta_s=None)
    c2.set_attr(T=450)
    nw.solve("design", init_only=True)
    assert nw.problem.structure_graph.structural_key() != key_first
