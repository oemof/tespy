# -*- coding: utf-8

"""Module for testing the structural decomposition of the equation system.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_decomposition.py

SPDX-License-Identifier: MIT
"""
from pytest import raises

from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.networks import Network
from tespy.solver import dulmage_mendelsohn
from tespy.tools.helpers import TESPyNetworkError


def test_diagonal_system_yields_scalar_blocks():
    decomposition = dulmage_mendelsohn({0: [0], 1: [1]}, 2)
    assert [block.kind for block in decomposition.blocks] == ["scalar"] * 2
    assert decomposition.precedence == {0: [], 1: []}
    assert decomposition.defective_blocks == []


def test_sequential_system_yields_precedence():
    decomposition = dulmage_mendelsohn({0: [0], 1: [0, 1]}, 2)
    assert [block.kind for block in decomposition.blocks] == ["scalar"] * 2
    first, second = decomposition.blocks
    assert first.equations == [0]
    assert second.equations == [1]
    assert decomposition.precedence[second.id] == [first.id]


def test_coupled_system_yields_square_block():
    decomposition = dulmage_mendelsohn({0: [0, 1], 1: [0, 1]}, 2)
    assert len(decomposition.blocks) == 1
    block = decomposition.blocks[0]
    assert block.kind == "square"
    assert block.equations == [0, 1]
    assert block.variables == [0, 1]


def test_independent_blocks_are_separated():
    decomposition = dulmage_mendelsohn(
        {0: [0], 1: [1], 2: [2, 3], 3: [2, 3]}, 4
    )
    kinds = sorted(block.kind for block in decomposition.blocks)
    assert kinds == ["scalar", "scalar", "square"]
    assert all(deps == [] for deps in decomposition.precedence.values())


def test_underdetermined_part_detected():
    decomposition = dulmage_mendelsohn({0: [0, 1]}, 2)
    assert len(decomposition.defective_blocks) == 1
    block = decomposition.defective_blocks[0]
    assert block.kind == "underdetermined"
    assert block.equations == [0]
    assert block.variables == [0, 1]


def test_overdetermined_part_detected():
    decomposition = dulmage_mendelsohn({0: [0], 1: [0]}, 1)
    assert len(decomposition.defective_blocks) == 1
    block = decomposition.defective_blocks[0]
    assert block.kind == "overdetermined"
    assert block.equations == [0, 1]
    assert block.variables == [0]


def test_mixed_defects_do_not_contaminate_square_part():
    # eq 0 solves var 0, eqs 1 and 2 both constrain var 1 (overdetermined),
    # var 2 has no equation (underdetermined)
    decomposition = dulmage_mendelsohn({0: [0], 1: [1], 2: [1]}, 3)
    kinds = {block.kind for block in decomposition.blocks}
    assert kinds == {"scalar", "overdetermined", "underdetermined"}
    square = [b for b in decomposition.blocks if b.kind == "scalar"][0]
    assert square.equations == [0]
    assert square.variables == [0]


def test_underdetermination_error_names_variables():
    nw = Network(iterinfo=False)
    so = Source("source")
    si = Sink("sink")
    c = Connection(so, "out1", si, "in1", label="c1")
    nw.add_conns(c)
    c.set_attr(fluid={"water": 1}, m=1, p=1e5)

    with raises(TESPyNetworkError) as error:
        nw.solve("design")
    assert "Structural analysis" in str(error.value)
    assert "The original problem holds" in str(error.value)
