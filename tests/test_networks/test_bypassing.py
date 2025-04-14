# -*- coding: utf-8

"""Module for testing network properties.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_network.py

SPDX-License-Identifier: MIT
"""
from tespy.components import HeatExchanger
from tespy.components import Merge
from tespy.components import Pump
from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Splitter
from tespy.components import Turbine
from tespy.components import Valve
from tespy.connections import Connection
from tespy.networks import Network


def test_bypass_system():
    nw = Network(p_unit="bar", T_unit="C", h_unit="kJ / kg")

    source = Source("In")
    sink = Sink("Out")
    heater = SimpleHeatExchanger("Heater")

    c1 = Connection(source, "out1", heater, "in1", "1")
    c2 = Connection(heater, "out1", sink, "in1", "2")

    nw.add_conns(c1, c2)

    c1.set_attr(T=100, p=2, m=1, fluid={"water":1})
    heater.set_attr(Q=2e6, pr=0.9)

    nw.solve("design")
    nw._convergence_check()

    msg = "In non-bypass mode, pressure ratio must be the specified value"
    assert round(c2.p.val, 4) == round(c1.p.val * 0.9, 4), msg
    msg = "In non-bypass mode, temperature must result from heat input"
    assert round(c2.T.val, 4) == 116.9113, msg

    heater.set_attr(bypass=True)

    nw.solve("design")
    nw._convergence_check()
    msg = "With bypass enabled, temperature an pressure should not change."
    assert round(c2.T.val, 4) == round(c1.T.val, 4), msg


def test_bypass_regenerative_preheater():
    nw = Network(p_unit='bar', T_unit='C', h_unit='kJ / kg')

    source = Source("In")
    turbine_hp = Turbine("high pressure turbine")
    extraction = Splitter("steam extraction")
    turbine_lp = Turbine("low pressure turbine")
    merge = Merge("condensate merge")
    main_condenser = SimpleHeatExchanger("main condenser")
    pump = Pump("feed water pump")

    preheater = HeatExchanger("preheater")
    valve = Valve("condensate throttle")

    sink = Sink("Out")


    c1 = Connection(source, "out1", turbine_hp, "in1", "1")
    c2 = Connection(turbine_hp, "out1", extraction, "in1", "2")
    c3 = Connection(extraction, "out1", turbine_lp, "in1", "3")
    c4 = Connection(turbine_lp, "out1", merge, "in1", "4")
    c5 = Connection(merge, "out1", main_condenser, "in1", "5")
    c6 = Connection(main_condenser, "out1", pump, "in1", "6")
    c7 = Connection(pump, "out1", preheater, "in2", "7")
    c8 = Connection(preheater, "out2", sink, "in1", "8")
    nw.add_conns(c1, c2, c3, c4, c5, c6, c7, c8)

    c11 = Connection(extraction, "out2", preheater, "in1", "11")
    c12 = Connection(preheater, "out1", valve, "in1", "12")
    c13 = Connection(valve, "out1", merge, "in2", "13")
    nw.add_conns(c11, c12, c13)

    c1.set_attr(T=500, p=100, m=10, fluid={"water":1})
    c2.set_attr(p=7)

    c4.set_attr(p=0.5)
    c6.set_attr(p=0.5, x=0)
    c7.set_attr(p=100)

    c11.set_attr(m=1)
    c12.set_attr(x=0)

    turbine_hp.set_attr(eta_s=0.9)
    turbine_lp.set_attr(eta_s=0.88)
    pump.set_attr(eta_s=0.75)
    dp1 = 0.2
    dp2 = 1
    preheater.set_attr(dp1=dp1, dp2=dp2)

    nw.solve("design")
    nw._convergence_check()

    assert round(c12.p.val, 2) == round(c11.p.val - dp1, 2)
    assert round(c8.p.val, 2) == round(c7.p.val - dp2, 2)
    assert round(c11.h.val, 2) > round(c12.h.val, 2)
    assert round(c8.h.val, 2) > round(c7.h.val, 2)

    preheater.set_attr(bypass=True)
    c12.set_attr(p=None, x=None)
    c8.set_attr(p=None)
    nw.solve("design")
    nw._convergence_check()
    assert round(c12.p.val, 2) == round(c11.p.val, 2)
    assert round(c8.p.val, 2) == round(c7.p.val, 2)
    assert round(c11.h.val, 2) == round(c12.h.val, 2)
    assert round(c8.h.val, 2) == round(c7.h.val, 2)
