# -*- coding: utf-8

"""Tests of the relative increment filter of the newton iteration.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_increment_filter.py

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
from tespy.tools.global_vars import ERR


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


@pytest.fixture
def prepared_problem():
    nw = _heat_pump()
    nw.solve("design", init_only=True)
    p = nw._problem
    n = p.variable_counter
    p.increment = np.full(n, np.inf)
    p.variable_weights = np.array([1e6, 1e5, 10.0, 1.0][:n] + [1.0] * max(0, n - 4))
    return p


def test_seed_is_never_filtered(prepared_problem):
    p = prepared_problem
    assert not p._increment_filter().any()


def test_zero_increments_always_filtered(prepared_problem):
    p = prepared_problem
    p.increment = np.zeros(p.variable_counter)
    p.increment[0] = 0.5 * p.variable_weights[0]
    mask = p._increment_filter()
    assert not mask[0]
    assert mask[1:].all()


def test_converged_variable_kept_only_while_contracting(prepared_problem):
    p = prepared_problem
    n = p.variable_counter
    # one variable barely moves, another one moves a lot: no filtering
    p.increment = np.full(n, 1e-8) * p.variable_weights
    p.increment[-1] = 0.1 * p.variable_weights[-1]
    assert not p._increment_filter().any()
    # everything moves little: the converged ones are filtered
    p.increment[-1] = 1e-4 * p.variable_weights[-1]
    mask = p._increment_filter()
    assert mask[:-1].all()
    assert not mask[-1]


def test_absolute_filter_when_disabled(prepared_problem):
    p = prepared_problem
    p.increment_filter_rtol = None
    p.increment = np.full(p.variable_counter, 1e-8) * p.variable_weights
    assert not p._increment_filter().any()
    p.increment = np.full(p.variable_counter, 0.5 * ERR ** 2)
    assert p._increment_filter().all()


def test_solve_results_independent_of_filter():
    results = {}
    iterations = {}
    for rtol in (None, 1e-6):
        nw = _heat_pump()
        nw.solve("design", increment_filter=rtol)
        nw.assert_convergence()
        results[rtol] = [
            value for c in nw.conns["object"]
            for value in (c.m.val_SI, c.p.val_SI, c.h.val_SI)
        ]
        iterations[rtol] = nw._problem.iter

    assert results[1e-6] == pytest.approx(results[None], rel=1e-8)
    assert iterations[1e-6] == iterations[None]


def test_invalid_increment_filter():
    # value errors of the solver options are reported as a crashed solve,
    # like an invalid line_search option
    nw = _heat_pump()
    nw.solve("design", increment_filter=2)
    assert nw.status == 99
