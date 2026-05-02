# -*- coding: utf-8

"""Tests for HeatConnection, HeatBus, HeatSource and HeatSink components.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tests/test_networks/test_heat_connections.py

SPDX-License-Identifier: MIT
"""

from pytest import approx

from tespy.components import HeatBus
from tespy.components import HeatSink
from tespy.components import HeatSource
from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.connections import HeatConnection
from tespy.networks import Network


def _fluid_network():
    """Return a minimal network with water flowing through a SimpleHeatExchanger."""
    nw = Network(iterinfo=False)
    nw.units.set_defaults(pressure="bar", temperature="degC", heat="kW")

    src = Source("fluid source")
    si = Sink("fluid sink")
    hx = SimpleHeatExchanger("heat exchanger", pr=1)

    c1 = Connection(src, "out1", hx, "in1", label="c1")
    c2 = Connection(hx, "out1", si, "in1", label="c2")
    nw.add_conns(c1, c2)

    c1.set_attr(fluid={"water": 1}, T=20, p=1, m=1)

    return nw, hx, c1, c2


class TestHeatConnectionOutlet:
    """SimpleHeatExchanger dissipates heat via a HeatConnection on the outlet."""

    def test_energy_matches_component_Q(self):
        nw, hx, _, _ = _fluid_network()

        ambient = HeatSink("ambient")
        h1 = HeatConnection(hx, "heat", ambient, "heat", label="h1")
        nw.add_conns(h1)

        hx.set_attr(Q=-10)
        nw.solve("design")
        nw.assert_convergence()

        assert h1.E.val_SI == approx(-hx.Q.val_SI, rel=1e-4)

    def test_energy_preserved_in_design_state(self):
        nw, hx, _, _ = _fluid_network()

        ambient = HeatSink("ambient")
        h1 = HeatConnection(hx, "heat", ambient, "heat", label="h1")
        nw.add_conns(h1)

        hx.set_attr(Q=-10)
        nw.solve("design")
        state = nw.save(as_dict=True)

        assert "HeatConnection" in state["Connection"]
        assert "h1" in state["Connection"]["HeatConnection"]
        assert state["Connection"]["HeatConnection"]["h1"]["E"] == approx(
            h1.E.val, rel=1e-4
        )


class TestHeatConnectionInlet:
    """SimpleHeatExchanger receives heat via a HeatConnection on the inlet."""

    def test_energy_matches_component_Q(self):
        nw, hx, c1, c2 = _fluid_network()

        source = HeatSource("heat source")
        h1 = HeatConnection(source, "heat", hx, "heat", label="h1")
        nw.add_conns(h1)

        h1.set_attr(E=10)
        nw.solve("design")
        nw.assert_convergence()

        assert hx.Q.val_SI == approx(h1.E.val_SI, rel=1e-4)

    def test_component_Q_drives_connection_energy(self):
        nw, hx, c1, c2 = _fluid_network()

        source = HeatSource("heat source")
        h1 = HeatConnection(source, "heat", hx, "heat", label="h1")
        nw.add_conns(h1)

        hx.set_attr(Q=10)
        nw.solve("design")
        nw.assert_convergence()

        assert h1.E.val_SI == approx(hx.Q.val_SI, rel=1e-4)


class TestHeatBus:
    """HeatBus distributes multiple heat flows while enforcing energy balance."""

    def test_single_in_single_out_balance(self):
        nw = Network(iterinfo=False)
        nw.units.set_defaults(heat="kW")

        src = HeatSource("source")
        bus = HeatBus("bus", num_in=1, num_out=1)
        si = HeatSink("sink")

        h1 = HeatConnection(src, "heat", bus, "heat_in1", label="h1")
        h2 = HeatConnection(bus, "heat_out1", si, "heat", label="h2")
        nw.add_conns(h1, h2)

        h1.set_attr(E=50)
        nw.solve("design")
        nw.assert_convergence()

        assert h2.E.val == approx(h1.E.val, rel=1e-4)

    def test_two_in_one_out_balance(self):
        nw = Network(iterinfo=False)
        nw.units.set_defaults(heat="kW")

        src1 = HeatSource("source 1")
        src2 = HeatSource("source 2")
        bus = HeatBus("bus", num_in=2, num_out=1)
        si = HeatSink("sink")

        h1 = HeatConnection(src1, "heat", bus, "heat_in1", label="h1")
        h2 = HeatConnection(src2, "heat", bus, "heat_in2", label="h2")
        h3 = HeatConnection(bus, "heat_out1", si, "heat", label="h3")
        nw.add_conns(h1, h2, h3)

        h1.set_attr(E=30)
        h2.set_attr(E=20)
        nw.solve("design")
        nw.assert_convergence()

        assert h3.E.val == approx(h1.E.val + h2.E.val, rel=1e-4)

    def test_one_in_two_out_balance(self):
        nw = Network(iterinfo=False)
        nw.units.set_defaults(heat="kW")

        src = HeatSource("source")
        bus = HeatBus("bus", num_in=1, num_out=2)
        si1 = HeatSink("sink 1")
        si2 = HeatSink("sink 2")

        h1 = HeatConnection(src, "heat", bus, "heat_in1", label="h1")
        h2 = HeatConnection(bus, "heat_out1", si1, "heat", label="h2")
        h3 = HeatConnection(bus, "heat_out2", si2, "heat", label="h3")
        nw.add_conns(h1, h2, h3)

        h2.set_attr(E=15)
        h3.set_attr(E=35)
        nw.solve("design")
        nw.assert_convergence()

        assert h1.E.val == approx(h2.E.val + h3.E.val, rel=1e-4)
