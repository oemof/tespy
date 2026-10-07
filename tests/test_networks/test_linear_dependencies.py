# -*- coding: utf-8

"""Module for testing some Networks that in older versions would have caused
linear dependencies.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_linear_dependencies.py

SPDX-License-Identifier: MIT
"""
from pytest import fixture

from tespy.components import Condenser
from tespy.components import CycleCloser
from tespy.components import HeatExchanger
from tespy.components import Merge
from tespy.components import PowerBus
from tespy.components import Pump
from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Splitter
from tespy.components import SteamTurbine
from tespy.components import Valve
from tespy.connections import Connection
from tespy.connections import PowerConnection
from tespy.networks import Network


@fixture
def nw():
    nw = Network()
    nw.units.set_defaults(
        temperature="degC", pressure="bar", pressure_difference="bar"
    )
    return nw


def test_preheater_evaporator_pinch_without_starting_values(nw):
    heat_source_in = Source("heat source inlet")
    heat_source_out = Sink("heat source outlet")

    wf_source = Source("working fluid source")
    preheater = HeatExchanger("preheater")
    evaporator = HeatExchanger("evaporator")
    wf_sink = Sink("working fluid sink")

    a1 = Connection(heat_source_in, "out1", evaporator, "in1", label="a1")
    a2 = Connection(evaporator, "out1", preheater, "in1", label="a2")
    a3 = Connection(preheater, "out1", heat_source_out, "in1", label="a3")

    b1 = Connection(wf_source, "out1", preheater, "in2", label="b4")
    b2 = Connection(preheater, "out2", evaporator, "in2", label="b5")
    b3 = Connection(evaporator, "out2", wf_sink, "in1", label="b0")

    nw.add_conns(a1, a2, a3)
    nw.add_conns(b1, b2, b3)

    working_fluid = "Isopentane"

    a1.set_attr(fluid={"water": 1}, p=10, T=160, m=100)

    b1.set_attr(td_bubble=60)
    b2.set_attr(td_bubble=4)
    b3.set_attr(fluid={working_fluid: 1}, x=1, T=140)

    preheater.set_attr(pr1=1, pr2=1)
    evaporator.set_attr(pr1=1, pr2=1, ttd_l=5)

    nw.solve("design")
    nw.assert_convergence()


def test_network_from_issue_922(nw):
    # https://github.com/oemof/tespy/discussions/922

    cw_source=Source("Cw_in")
    da=Merge("da",num_in=3)
    cond_mix=Merge("cond_mix",num_in=3)
    cep=Pump("cep")
    bfp=Pump("bfp")
    hpt1 = SteamTurbine("High Pressure Turbine1")
    ipt1_0 = SteamTurbine("IP Turbine1")
    ipt1_1 = SteamTurbine("IP Turbine2")
    lpt1_0 = SteamTurbine("LP Turbine1")
    lpt1_1 = SteamTurbine("LP Turbine2")
    v1=Valve("lphdrip")
    v2=Valve("hphdrip")
    v3=Valve("cep_disch")
    #Introduce LP Turbine and LP Heater
    bfpt=SteamTurbine("bfpt")
    cond1= Condenser("cond1")
    cw_sink = Sink("cw_sink")
    s1=Splitter("s1",num_out=3) #HPT Outlet
    shaft2=PowerBus("shaft2", num_in=1, num_out=1)
    #for introducing HP Heater
    s2=Splitter("s2",num_out=2)

    # for LP Trubine
    s3=Splitter("s3",num_out=2)

    sg1 = SimpleHeatExchanger("steam generator1")
    sg2 = SimpleHeatExchanger("steam generator2")
    cc = CycleCloser("cycle closer")

    c1_0 = Connection(cc, "out1", hpt1, "in1", label="c1")
    c1 = Connection(sg1, "out1", cc, "in1")

    #introducing HP Heater from CRH
    hph1= Condenser("hph1")
    lph1= Condenser("lph1")

    #Rheater Introduced
    c2 = Connection(s2, "out1", sg2, "in1")
    c2_5 = Connection(hpt1, "out1", s2, "in1")
    c2_6 = Connection(s2, "out2", hph1, "in1")
    c2_8 = Connection(hph1, "out1", v2, "in1")
    c2_7 = Connection(v2, "out1", da, "in3")
    c2_4 = Connection(sg2, "out1", ipt1_0, "in1")
    c2_0 = Connection(ipt1_0, "out1", s1, "in1")

    c2_1 = Connection(s1, "out1", ipt1_1, "in1")
    c2_2 = Connection(s1, "out2", da, "in2")
    c2_3 = Connection(s1, "out3", bfpt, "in1")

    #cooling water connections
    c3 = Connection(cw_source, "out1", cond1, "in2")
    c4 = Connection(cond1, "out2", cw_sink, "in1")
    c5=Connection(cond1, "out1", cep, "in1")

    #introduce LPH1
    c6_1=Connection(cep, "out1", lph1, "in2")
    c6_4=Connection(lph1, "out2", v3, "in1")
    c6_3=Connection(v3, "out1", da, "in1")
    c6_2=Connection(da, "out1", bfp, "in1")
    #hph1
    c7=Connection(bfp, "out1", hph1, "in2")
    c7_1=Connection(hph1, "out2", sg1, "in1")

    # introduce connections for LP Turbine and lph
    c9=Connection(ipt1_1, "out1", lpt1_0, "in1")
    c10=Connection(lpt1_0, "out1", s3, "in1")
    c10_1=Connection(s3, "out1", lpt1_1, "in1")
    c10_2=Connection(s3, "out2", lph1, "in1")
    c10_3=Connection(lph1, "out1", v1, "in1")
    c10_4=Connection(v1, "out1", cond_mix, "in3")
    c8_1=Connection(lpt1_1, "out1", cond_mix, "in1")
    c8_2=Connection(bfpt, "out1", cond_mix, "in2")
    c8_3=Connection(cond_mix, "out1", cond1, "in1")
    e8=PowerConnection(bfpt, "power", shaft2, "power_in1", label="e8")
    e9=PowerConnection(shaft2, "power_out1", bfp, "power", label="e9")

    nw.add_conns(
        c1, c1_0, c2, c2_0, c2_1, c2_2, c2_3, c2_4, c2_5, c2_6, c2_7, c2_8,
        c3, c4, c5, c6_1, c6_2, c6_3, c6_4, c7, c7_1, c8_1, c8_2, c8_3, c9,
        c10, c10_1, c10_2, c10_3, c10_4, e8, e9
    )

    cep.set_attr(eta_s=0.9)
    bfp.set_attr(eta_s=0.9)
    # 4. Set Parameters
    c1_0.set_attr(p=256, T=600,m=1827, fluid={"Water": 1})
    c2.set_attr(p=40)
    c2_0.set_attr(p=12)
    c2_4.set_attr(T=610)
    cond1.set_attr(ttd_u=2.5)
    c3.set_attr(p=2, T=33, fluid={"Water": 1})
    c4.set_attr( p=1.5,T=42)
    cond1.set_attr(pr1=1)
    c6_2.set_attr(x=0)
    sg1.set_attr(pr=0.9)
    sg2.set_attr(pr=0.9)
    c9.set_attr(p=5.5)
    c10.set_attr(p=2)
    hpt1.set_attr(eta_s=0.9)
    ipt1_0.set_attr(eta_s=0.9)
    ipt1_1.set_attr(eta_s=0.9)
    lpt1_0.set_attr(eta_s=0.9)
    lpt1_1.set_attr(eta_s=0.9)
    bfpt.set_attr(eta_s=0.9)
    hph1.set_attr(ttd_u=5)
    hph1.set_attr(pr2=.95)
    lph1.set_attr(ttd_u=4.31)
    lph1.set_attr(pr2=0.95)
    v1.set_attr(pr=0.95)
    v2.set_attr(pr=0.95)
    v3.set_attr(pr=0.5)
    nw.solve("design")
    nw.assert_convergence()
