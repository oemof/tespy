# -*- coding: utf-8

"""Module for testing components of type merge.
This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_merge.py
SPDX-License-Identifier: MIT
"""

from numpy.testing import assert_almost_equal
from pytest import approx
from pytest import fixture

from tespy.components import Merge
from tespy.components import Node
from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Splitter
from tespy.connections import Connection
from tespy.networks import Network


@fixture
def node_network():
    nw = Network()

    so1 = Source("source1")
    so2 = Source("source2")
    si1 = Sink("sink1")
    si2 = Sink("sink2")

    node = Node("node")
    node.set_attr(num_in=2, num_out=2)

    c1 = Connection(so1, "out1", node, "in1", label="i1")
    c2 = Connection(so2, "out1", node, "in2", label="i2")
    c3 = Connection(node, "out1", si1, "in1", label="o1")
    c4 = Connection(node, "out2", si2, "in1", label="o2")

    nw.add_conns(c1, c2, c3, c4)

    return nw


@fixture
def merge_splitter_network():
    nw = Network()

    so1 = Source("source1")
    so2 = Source("source2")
    si1 = Sink("sink1")
    si2 = Sink("sink2")

    merge = Merge("merge")
    splitter = Splitter("splitter")
    merge.set_attr(num_in=2)
    splitter.set_attr(num_out=2)

    c1 = Connection(so1, "out1", merge, "in1", label="i1")
    c2 = Connection(so2, "out1", merge, "in2", label="i2")
    c0 = Connection(merge, "out1", splitter, "in1", label="c1")
    c3 = Connection(splitter, "out1", si1, "in1", label="o1")
    c4 = Connection(splitter, "out2", si2, "in1", label="o2")

    nw.add_conns(c1, c2, c0, c3, c4)

    return nw


def test_node_equals_merge_splitter_model(node_network, merge_splitter_network):

    c1, c2, c3 = node_network.get_conn(["i1", "i2", "o1"])
    c1.set_attr(fluid={"O2": 1}, m=5, p=1e5, T=300)
    c2.set_attr(fluid={"N2": 1}, m=5, T=350)
    c3.set_attr(m=2)

    node_network.solve("design")

    c1, c2, c3 = merge_splitter_network.get_conn(["i1", "i2", "o1"])
    c1.set_attr(fluid={"O2": 1}, m=5, p=1e5, T=300)
    c2.set_attr(fluid={"N2": 1}, m=5, T=350)
    c3.set_attr(m=2)

    merge_splitter_network.solve("design")

    # the node network has the index subset
    index = node_network.results["Connection"].index
    number_cols = node_network.results["Connection"].dtypes == float
    # this transposes the mask into the actual column names so the order
    # in both networks is the same
    number_cols = number_cols[number_cols].index

    assert_almost_equal(
        node_network.results["Connection"].loc[index, number_cols].values,
        merge_splitter_network.results["Connection"].loc[index, number_cols].values
    )


def test_node_with_set_outflow_composition(node_network):
    c1, c2, c3, c4 = node_network.get_conn(["i1", "i2", "o1", "o2"])

    c1.set_attr(fluid={"O2": 1}, m=5, p=1e5, T=300)
    c2.set_attr(fluid={"N2": 1}, T=350)
    c3.set_attr(m=3)
    c4.set_attr(fluid={"O2": 0.4, "N2": 0.6})

    node_network.solve("design")

    assert approx(c3.m.val_SI + c4.m.val_SI) == c1.m.val_SI + c2.m.val_SI
    assert approx(c3.fluid.val["O2"]) == c4.fluid.val["O2"]
    assert (
        approx((c3.m.val_SI  + c4.m.val_SI) * c4.fluid.val["O2"])
        == c1.fluid.val["O2"] * c1.m.val_SI
    )
    assert approx(c3.h.val_SI) == c4.h.val_SI
    assert (
        approx(c3.h.val_SI * (c1.m.val_SI + c2.m.val_SI))
        == (c1.h.val_SI * c1.m.val_SI + c2.h.val_SI * c2.m.val_SI)
    )

def test_node_with_set_outflow_composition_and_temperature(node_network):
    c1, c2, c3, c4 = node_network.get_conn(["i1", "i2", "o1", "o2"])

    c1.set_attr(fluid={"O2": 1}, m=5, p=1e5, T=300)
    c2.set_attr(fluid={"N2": 1})
    c3.set_attr(m=3)
    c4.set_attr(fluid={"O2": 0.4, "N2": 0.6}, T=340)

    node_network.solve("design")

    assert approx(c3.m.val_SI + c4.m.val_SI) == c1.m.val_SI + c2.m.val_SI
    assert approx(c3.fluid.val["O2"]) == c4.fluid.val["O2"]
    assert (
        approx((c3.m.val_SI  + c4.m.val_SI) * c4.fluid.val["O2"])
        == c1.fluid.val["O2"] * c1.m.val_SI
    )
    assert approx(c3.h.val_SI) == c4.h.val_SI
    assert (
        approx(c3.h.val_SI * (c1.m.val_SI + c2.m.val_SI))
        == (c1.h.val_SI * c1.m.val_SI + c2.h.val_SI * c2.m.val_SI)
    )


def test_node_with_set_pressure_at_outlet(node_network):
    c1, c2, c3, c4 = node_network.get_conn(["i1", "i2", "o1", "o2"])

    c1.set_attr(fluid={"O2": 1}, m=5, T=300)
    c2.set_attr(fluid={"N2": 1})
    c3.set_attr(m=3, p=1e5)
    c4.set_attr(m=12, T=340)

    node_network.solve("design")

    assert approx(c3.m.val_SI + c4.m.val_SI) == c1.m.val_SI + c2.m.val_SI
    assert approx(c3.fluid.val["O2"]) == c4.fluid.val["O2"]
    assert (
        approx((c3.m.val_SI  + c4.m.val_SI) * c4.fluid.val["O2"])
        == c1.fluid.val["O2"] * c1.m.val_SI
    )
    assert approx(c3.h.val_SI) == c4.h.val_SI
    assert (
        approx(c3.h.val_SI * (c1.m.val_SI + c2.m.val_SI))
        == (c1.h.val_SI * c1.m.val_SI + c2.h.val_SI * c2.m.val_SI)
    )


def test_splitter_outlet_order_beyond_nine_ports():
    nw = Network()
    nw.units.set_defaults(pressure="bar", enthalpy="kJ/kg")
    num_out = 11
    sp = Splitter("splitter", num_out=num_out)

    conns = []
    for i in range(1, num_out + 1):
        c = Connection(sp, f"out{i}", Sink(f"sink {i}"), "in1", label=f"out{i}")
        c.set_attr(m=i)
        conns.append(c)

    c_in = Connection(Source("source"), "out1", sp, "in1", label="in")
    nw.add_conns(c_in, *conns)
    c_in.set_attr(fluid={"water": 1}, p=1, h=100)

    nw.solve("design")
    nw.assert_convergence()

    assert [c.source_id for c in sp.outl] == [f"out{i}" for i in range(1, num_out + 1)]


def test_node_closing_recirculation_loop():
    nw = Network()
    nw.units.set_defaults(pressure="bar", temperature="degC")
    node = Node("node")
    heater = SimpleHeatExchanger("heater")
    c1 = Connection(Source("source"), "out1", node, "in1", label="1")
    c2 = Connection(node, "out1", Sink("sink"), "in1", label="2")
    c3 = Connection(node, "out2", heater, "in1", label="3")
    c4 = Connection(heater, "out1", node, "in2", label="4")
    nw.add_conns(c1, c2, c3, c4)

    c1.set_attr(fluid={"water": 1}, m=1, p=1, T=20)
    c3.set_attr(m=0.5)
    heater.set_attr(Q=1e4)

    nw.solve("design")
    nw.assert_convergence()

    assert c2.m.val_SI == approx(1)
    assert c2.h.val_SI == approx(c3.h.val_SI)
