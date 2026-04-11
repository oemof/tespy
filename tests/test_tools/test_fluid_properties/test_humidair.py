# -*- coding: utf-8

"""Module for testing fluid properties of HAConnection.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_tools/test_fluid_properties/test_humidair.py

SPDX-License-Identifier: MIT
"""
from CoolProp.CoolProp import HAPropsSI
from pytest import fixture

from tespy.components import MovingBoundaryHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.connections import HAConnection
from tespy.networks import Network


@fixture
def humidair_network():
    nw = Network()
    nw.units.set_defaults(
        temperature="°C",
        pressure="bar",
        heat="kW"
    )


    so = Source("source")
    hex = MovingBoundaryHeatExchanger("heat exchanger")

    si = Sink("sink")

    so2 = Source("refrigerant source")
    si2 = Sink("refrigerant sink")

    a1 = Connection(so2, "out1", hex, "in2", label="a1")
    a2 = Connection(hex, "out2", si2, "in1", label="a2")

    c1 = HAConnection(so, "out1", hex, "in1", label="c1")
    c2 = HAConnection(hex, "out1", si, "in1", label="c2")

    nw.add_conns(a1, a2, c1, c2)

    return nw


def test_specification_humidity_ratio(humidair_network):

    nw = humidair_network
    a1, a2, c1, c2 = nw.get_conn(["a1", "a2", "c1", "c2"])
    hex = nw.get_comp("heat exchanger")

    a1.set_attr(fluid={"NH3": 1}, m=2, x=0.3)
    a2.set_attr(td_dew=5, T_dew=-20)

    c1.set_attr(
        p=1,
        T=15,
        w=0.004 #-> use this to directly impose fluid composition, then do not specify r and fluid balance for this example
    )
    c2.set_attr(T=5)

    hex.set_attr(dp1=0, dp2=0)

    nw.solve("design")
    nw.assert_convergence()


def test_specification_relative_humidity(humidair_network):

    nw = humidair_network
    a1, a2, c1, c2 = nw.get_conn(["a1", "a2", "c1", "c2"])
    hex = nw.get_comp("heat exchanger")

    a1.set_attr(fluid={"NH3": 1}, m=2, x=0.3)
    a2.set_attr(td_dew=5, T_dew=-20)

    c1.set_attr(
        p=1,
        T=15,
        r=0.9,  # relative humidity
        fluid0={"air": 0.99, "water": 0.01}, # there must be some kind of composition information present
        fluid_balance=True  # fluid balance is not automatically closed
    )
    c2.set_attr(T=5)

    hex.set_attr(dp1=0, dp2=0)

    nw.solve("design")
    nw.assert_convergence()
