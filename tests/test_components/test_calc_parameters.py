# -*- coding: utf-8

"""Tests for the calc= parameter dispatch and _topological_sort.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_calc_parameters.py

SPDX-License-Identifier: MIT
"""

import pytest

from tespy.components import Sink
from tespy.components import Source
from tespy.components.basics.subsystem_interface import SubsystemInterface
from tespy.components.component import _topological_sort
from tespy.connections import Connection
from tespy.networks import Network
from tespy.tools.data_containers import ComponentProperties as dc_cp

# ---------------------------------------------------------------------------
# Minimal stub for unit tests – only needs the calc_deps attribute
# ---------------------------------------------------------------------------

class _DC:
    def __init__(self, calc_deps=None):
        self.calc_deps = calc_deps or []


# ---------------------------------------------------------------------------
# Unit tests for _topological_sort
# ---------------------------------------------------------------------------

def test_topological_sort_dag():
    """A chain a→b→c is returned with a before b before c."""
    items = {
        "a": _DC(),
        "b": _DC(calc_deps=["a"]),
        "c": _DC(calc_deps=["b"]),
    }
    order = _topological_sort(items)
    assert set(order) == {"a", "b", "c"}
    assert order.index("a") < order.index("b")
    assert order.index("b") < order.index("c")


def test_topological_sort_cycle_raises():
    """A mutual dependency a↔b is detected and raises ValueError."""
    items = {
        "a": _DC(calc_deps=["b"]),
        "b": _DC(calc_deps=["a"]),
    }
    with pytest.raises(ValueError, match="Cycle detected"):
        _topological_sort(items)


# ---------------------------------------------------------------------------
# Integration test: SubsystemInterface subclass with fake calc= parameters
# ---------------------------------------------------------------------------

class _FakeComponent(SubsystemInterface):
    """Pass-through component carrying three chained result parameters."""

    def get_parameters(self):
        params = super().get_parameters()
        params["base"] = dc_cp(
            is_result=True,
            calc=lambda: 3.0,
        )
        params["derived"] = dc_cp(
            is_result=True,
            calc=lambda: self.base.val_SI * 2.0,
            calc_deps=["base"],
        )
        params["leaf"] = dc_cp(
            is_result=True,
            calc=lambda: self.derived.val_SI + self.base.val_SI,
            calc_deps=["base", "derived"],
        )
        return params


def test_calc_parameters_respects_dependency_order():
    """calc= parameters are evaluated in dependency order after solving."""
    nw = Network(iterinfo=False)
    nw.units.set_defaults(pressure="bar", pressure_difference="bar", temperature="degC")

    so = Source("so")
    si = Sink("si")
    comp = _FakeComponent("comp")

    c1 = Connection(so, "out1", comp, "in1")
    c2 = Connection(comp, "out1", si, "in1")
    nw.add_conns(c1, c2)

    c1.set_attr(fluid={"water": 1}, T=20, p=1, m=1)
    nw.solve("design")
    nw.assert_convergence()

    assert comp.base.val_SI == pytest.approx(3.0)
    assert comp.derived.val_SI == pytest.approx(6.0)   # base * 2
    assert comp.leaf.val_SI == pytest.approx(9.0)      # derived + base
