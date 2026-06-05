
from pytest import mark

from tespy.components import Compressor
from tespy.components import HeatExchanger
from tespy.components import Merge
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Splitter
from tespy.components import Turbine
from tespy.components import Valve
from tespy.connections import Connection
from tespy.connections import Ref
from tespy.networks import Network


def test_R601_converges():
    orc = Network()

    caso = Source('cooling air source')
    casi = Sink('cooling air sink')

    turb = Turbine('turbine')
    cond = HeatExchanger('condenser')

    c1 = Connection(caso, 'out1', cond, 'in2', 'c1')
    c2 = Connection(cond, 'out2', casi, 'in1', 'c2')

    b1 = Connection(Source("source"), 'out1', turb, 'in1', 'b1')
    b2 = Connection(turb, 'out1', cond, 'in1', 'b2')
    b3 = Connection(cond, 'out1', Sink("sink"), 'in1', 'b3')

    orc.add_conns(c1, c2, b1, b2, b3)

    # boundary conditions
    air = {"Air":1}
    orcfluid = {"R601": 1}

    lsT = 140 + 273.15
    lsx = 1

    coT = 20 + 273.15
    cop = 1e5
    pr_lossless = 1

    eta_t = 0.85
    ttd_l_c = 15
    t_rise = 10

    c1.set_attr(fluid=air, p=cop, T=coT)
    c2.set_attr(T=Ref(c1, 1, t_rise))
    b1.set_attr(T=lsT, x=lsx, fluid=orcfluid, m=1)
    b3.set_attr(x=0)

    turb.set_attr(eta_s = eta_t)
    cond.set_attr(pr1=pr_lossless, pr2=pr_lossless)
    cond.set_attr(ttd_l=ttd_l_c)

    orc.solve('design')
    assert orc.status == 2


def test_two_compressors_side_stream_bypass():
    nw = Network()

    inlet = Source("Inlet")
    side_stream = Source("side stream source")
    outlet = Sink("outlet")
    imerge = Merge("Inlet Merge")
    smerge = Merge("Side Stream Merge")
    comp1 = Compressor("Comp1")
    comp2 = Compressor("Comp2")
    splitter = Splitter("Splitter")
    bypass = Valve("Bypass")

    c01 = Connection(inlet, "out1", imerge, "in1", label="Inlet")
    c02 = Connection(imerge, "out1", comp1, "in1", label="Compressor Inlet")
    c03 = Connection(comp1, "out1", smerge, "in1", label="Compressor 1 Outlet")
    s01 = Connection(side_stream, "out1", smerge, "in2", label="Side Stream")
    c04 = Connection(smerge, "out1", comp2, "in1", label="Compressor 2 Inlet")
    c05 = Connection(comp2, "out1", splitter, "in1", label="Compressor 2 Outlet")
    c06 = Connection(splitter, "out2", bypass, "in1", label="Bypass Inlet")
    c07 = Connection(bypass, "out1", imerge, "in2", label="Bypass Outlet")
    c08 = Connection(splitter, "out1", outlet, "in1", label="Outlet")
    nw.add_conns(c01, c02, c03, c04, c05, c06, c07, c08, s01)

    c01.set_attr(fluid={"Water": 1}, x=1, p=1e5, m=1)
    s01.set_attr(fluid={"Water": 1}, x=1, p=2e5, m=1)
    c06.set_attr(m=0)
    c08.set_attr(p=3e5)
    comp1.set_attr(eta_s=0.85)
    comp2.set_attr(eta_s=0.85)

    nw.solve("design")
    assert nw.status == 0


def _make_single_compressor_bypass():
    """Return (nw, c01, c02, c03, c04, c05, c06, comp) with base specs set."""
    nw = Network()
    inlet = Source("Inlet")
    merge = Merge("Merge")
    comp = Compressor("Compressor")
    splitter = Splitter("Splitter")
    bypass = Valve("Bypass")
    outlet = Sink("outlet")

    c01 = Connection(inlet, "out1", merge, "in1", label="Inlet")
    c02 = Connection(merge, "out1", comp, "in1", label="Compressor Inlet")
    c03 = Connection(comp, "out1", splitter, "in1", label="Compressor Outlet")
    c04 = Connection(splitter, "out2", bypass, "in1", label="Bypass Inlet")
    c05 = Connection(bypass, "out1", merge, "in2", label="Bypass Outlet")
    c06 = Connection(splitter, "out1", outlet, "in1", label="Outlet")
    nw.add_conns(c01, c02, c03, c04, c05, c06)

    c01.set_attr(fluid={"Water": 1}, x=1, p=1e5, m=1)
    comp.set_attr(eta_s=0.85, pr=1.5)
    return nw, c01, c02, c03, c04, c05, c06, comp


def test_single_compressor_bypass_ref_bypass_to_outlet():
    nw, c01, c02, c03, c04, c05, c06, comp = _make_single_compressor_bypass()
    c04.set_attr(m=Ref(c03, 0.01, 0))
    nw.solve("design")
    assert nw.status == 0


def test_single_compressor_bypass_ref_bypass_to_outlet_with_m0():
    nw, c01, c02, c03, c04, c05, c06, comp = _make_single_compressor_bypass()
    c04.set_attr(m=Ref(c03, 0.01, 0), m0=0.01)
    nw.solve("design")
    assert nw.status == 0


def test_single_compressor_bypass_ref_outlet_to_comp_outlet():
    nw, c01, c02, c03, c04, c05, c06, comp = _make_single_compressor_bypass()
    c06.set_attr(m=Ref(c03, 0.1, 0))
    nw.solve("design")
    assert nw.status == 0


def test_single_compressor_bypass_ref_outlet_to_comp_outlet_with_m0():
    nw, c01, c02, c03, c04, c05, c06, comp = _make_single_compressor_bypass()
    c06.set_attr(m=Ref(c03, 0.1, 0), m0=1)
    nw.solve("design")
    assert nw.status == 0


def test_single_compressor_bypass_ref_comp_outlet_to_outlet():
    nw, c01, c02, c03, c04, c05, c06, comp = _make_single_compressor_bypass()
    c03.set_attr(m=Ref(c06, 1 / 0.1, 0))
    nw.solve("design")
    assert nw.status == 0
