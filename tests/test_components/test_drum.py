# -*- coding: utf-8

"""Module for testing components of type merge.
This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_merge.py
SPDX-License-Identifier: MIT
"""
from pytest import approx
from pytest import fixture

from tespy.components import Drum
from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Splitter
from tespy.connections import Connection
from tespy.networks import Network


@fixture()
def drum_network_setup():
    nw = Network(T_unit="C", p_unit="bar")
    dr = Drum("drum")
    so = Source("liquid")
    si = Sink("vapor")
    c1 = Connection(so, "out1", dr, "in1", label="1")
    c2 = Connection(dr, "out2", si, "in1", label="2")
    nw.add_conns(c1, c2)
    c1.set_attr(fluid={"R290": 1}, m=10, Td_bp=-10, T=50)
    yield nw


def test_drum_with_blowdown(drum_network_setup):
    nw = drum_network_setup
    sp = Splitter("blowdown splitter")
    si = Sink("blowdown sink")
    eva = SimpleHeatExchanger("evaporator")
    dr = nw.get_comp("drum")
    c3 = Connection(dr, "out1", sp, "in1", label="3")
    c4 = Connection(sp, "out1", si, "in1", label="4")
    c5 = Connection(sp, "out2", eva, "in1", label="5")
    c6 = Connection(eva, "out1", dr, "in2", label="6")

    c6.set_attr(x=0.7, m=15)

    nw.add_conns(c3, c4, c5, c6)

    nw.solve("design")
    assert nw.status == 0
    assert 0.72728 == approx(c4.m.val_SI)
