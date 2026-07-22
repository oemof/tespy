# -*- coding: utf-8

"""Module for testing block-wise solving of the equation system.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_block_solve.py

SPDX-License-Identifier: MIT
"""
from pytest import approx
from pytest import raises

from tespy.components import Compressor
from tespy.components import Merge
from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Splitter
from tespy.connections import Connection
from tespy.networks import Network
from tespy.tools.helpers import TESPyNetworkError


def _compressor_network(eta_s_and_T=False):
    nw = Network(iterinfo=False)
    so = Source("source")
    compressor = Compressor("compressor")
    si = Sink("sink")
    c1 = Connection(so, "out1", compressor, "in1", label="c1")
    c2 = Connection(compressor, "out1", si, "in1", label="c2")
    nw.add_conns(c1, c2)
    c1.set_attr(fluid={"air": 1}, m=1, p=1e5, T=293.15)
    if eta_s_and_T:
        # pressure ratio free: eta_s and T couple outlet p and h in a
        # single square block
        compressor.set_attr(eta_s=0.8)
        c2.set_attr(T=450)
    else:
        compressor.set_attr(pr=3, eta_s=0.8)
    return nw


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
    heater1.set_attr(pr=0.98, Q=1e4)
    heater2.set_attr(Q=2e4)
    return nw


def _assert_same_results(build_network, **kwargs):
    reference = build_network(**kwargs)
    reference.solve("design")
    reference.assert_convergence()

    nw = build_network(**kwargs)
    nw.solve("design", block_solve=True)
    nw.assert_convergence()

    for c_ref in reference.conns["object"]:
        c = nw.get_conn(c_ref.label)
        for prop in ("m", "p", "h"):
            assert c.get_attr(prop).val_SI == approx(
                c_ref.get_attr(prop).val_SI, rel=1e-5
            ), f"{c_ref.label} ({prop})"
    return nw


def test_block_solve_scalar_chain():
    nw = _assert_same_results(_compressor_network)
    blocks = nw.problem.decompose().blocks
    assert [block.kind for block in blocks] == ["scalar"]
    assert all(block.status == 0 for block in blocks)


def test_block_solve_square_block():
    nw = _assert_same_results(_compressor_network, eta_s_and_T=True)
    blocks = nw.problem.decompose().blocks
    assert [block.kind for block in blocks] == ["square"]
    assert blocks[0].variables == [0, 1]


def test_block_inspection_methods():
    nw = _splitter_merge_network()
    nw.presolve("design")

    blocks = nw.get_blocks()
    assert all(block["kind"] == "scalar" for block in blocks)
    assert all(block["status"] is None for block in blocks)
    for block in blocks:
        assert all(dep < block["id"] for dep in block["prerequisites"])

    nw.print_blocks()
    nw.print_incidence_matrix(block_order=True)

    nw.solve_continue()
    nw.assert_convergence()
    assert all(block["status"] == 0 for block in nw.get_blocks())


def _compressor_heater_network():
    nw = Network(iterinfo=False)
    so = Source("source")
    compressor = Compressor("compressor")
    heater = SimpleHeatExchanger("heater")
    si = Sink("sink")
    c1 = Connection(so, "out1", compressor, "in1", label="c1")
    c2 = Connection(compressor, "out1", heater, "in1", label="c2")
    c3 = Connection(heater, "out1", si, "in1", label="c3")
    nw.add_conns(c1, c2, c3)
    c1.set_attr(fluid={"air": 1}, m=1, p=1e5, T=293.15)
    # pressure ratio free: eta_s and T couple outlet p and h in a square
    # block, the heater outlet enthalpy follows as a scalar block
    compressor.set_attr(eta_s=0.8)
    c2.set_attr(T=450)
    heater.set_attr(pr=0.98, Q=-1e4)
    return nw


def test_pause_on_block_failure():
    # the square compressor block needs more than two newton iterations,
    # so it deterministically fails with max_iter=2 and pauses the solver
    nw = _compressor_heater_network()
    nw.solve(
        "design", max_iter=2, pause_on_block_failure=True
    )
    assert nw.status == 20
    assert not nw.converged
    assert nw.problem.paused

    failed = nw.problem._paused_driver.paused_block
    variables = nw.get_variable_values(block=failed.id)
    assert set(variables) < set(nw.get_variable_values())
    nw.print_variable_values(block=failed.id)
    with raises(KeyError):
        nw.get_variable_values(block=-1)

    # the linear system of the failing iteration is recorded together with
    # the declared incidence, the variable values and the states of the
    # involved connections
    data = nw.get_block_jacobian(failed.id)
    assert data["jacobian"].shape == (
        len(data["equations"]), len(data["variables"])
    )
    assert data["residual"].shape == (len(data["equations"]),)
    assert data["expected"].shape == data["jacobian"].shape
    assert len(data["values"]) == len(data["variables"])
    nw.print_block_jacobian(failed.id)
    # current states show the entry state while paused, the failure states
    # cover the same connections
    states = nw.get_block_states(failed.id)
    labels = [state["label"] for state in states]
    assert "c2" in labels
    failure_states = nw.get_block_states(failed.id, at="failure")
    assert [state["label"] for state in failure_states] == labels
    nw.print_block_states(failed.id)
    nw.print_block_states(failed.id, at="failure")
    with raises(ValueError):
        nw.get_block_states(failed.id, at="entry")
    # the later block was never attempted, so there is nothing recorded
    with raises(TESPyNetworkError):
        nw.get_block_jacobian(failed.id + 1)
    with raises(TESPyNetworkError):
        nw.get_block_states(failed.id + 1, at="failure")
    # the current states of an unattempted block are available
    assert nw.get_block_states(failed.id + 1)

    # only the failed block's variables are open for modification
    with raises(TESPyNetworkError):
        nw.set_variable_value("c3", "h", 4e5)
    nw.set_variable_value("c2", "h", 4.5e5)

    # the retry runs with a sufficient iteration budget and continues
    # through the remaining blocks
    nw.solve_continue()
    assert not nw.problem.paused
    nw.assert_convergence()


def test_pause_then_escalate():
    nw = _compressor_heater_network()
    nw.solve("design", max_iter=2, pause_on_block_failure=True)
    assert nw.status == 20

    # turning the pause off hands the still failing block over to the
    # standard escalation stages
    nw.solve_continue(max_iter=2, pause_on_block_failure=False)
    assert nw.status != 20
    assert not nw.problem.paused


def test_pause_released_by_fresh_solve():
    nw = _compressor_heater_network()
    nw.solve("design", max_iter=2, pause_on_block_failure=True)
    assert nw.status == 20

    nw.solve("design")
    nw.assert_convergence()


def test_block_solve_sequential_precedence():
    nw = _assert_same_results(_splitter_merge_network)
    decomposition = nw.problem.decompose()
    assert all(block.kind == "scalar" for block in decomposition.blocks)
    assert all(block.status == 0 for block in decomposition.blocks)
    # at least one block must depend on an earlier one (heater 2 outlet
    # enthalpy needs the mass flow of its branch first)
    assert any(deps for deps in decomposition.precedence.values())
    for block_id, prerequisites in decomposition.precedence.items():
        assert all(dep < block_id for dep in prerequisites)
