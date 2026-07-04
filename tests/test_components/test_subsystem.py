from pytest import approx
from pytest import raises

from tespy.components import Compressor
from tespy.components import HeatSink
from tespy.components import HeatSource
from tespy.components import PowerSource
from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Splitter
from tespy.components import Subsystem
from tespy.connections import Connection
from tespy.connections import HeatConnection
from tespy.connections import PowerConnection
from tespy.networks import Network
from tespy.tools.helpers import TESPyComponentError


class MySubsystem(Subsystem):
    def __init__(self, label):
        self.num_in = 1
        self.num_out = 2
        super().__init__(label)

    def create_network(self):
        splitter = Splitter("splitter", num_out=2)

        c = Connection(self.inlet, "out1", splitter, "in1")
        self.add_conns(c)
        c = Connection(splitter, "out1", self.outlet, "in1")
        self.add_conns(c)
        c = Connection(splitter, "out2", self.outlet, "in2")
        self.add_conns(c)


def test_subsystem_add_connections_individually():
    sub = MySubsystem("my subsystem")
    assert len(sub.comps) == 3


class SubsystemWithDuplicateConnectionLabel(Subsystem):
    def __init__(self, label):
        self.num_in = 2
        self.num_out = 2
        super().__init__(label)

    def create_network(self):
        c1 = Connection(self.inlet, "out1", self.outlet, "in1", label="c1")
        c2 = Connection(self.inlet, "out2", self.outlet, "in2", label="c1")
        self.add_conns(c1, c2)


def test_subsystem_duplicated_connection_label():
    with raises(TESPyComponentError):
        SubsystemWithDuplicateConnectionLabel("test")


class SubsystemWithDuplicateComponentLabel(Subsystem):
    def __init__(self, label):
        self.num_in = 2
        self.num_out = 2
        super().__init__(label)

    def create_network(self):

        heatex1 = SimpleHeatExchanger("heatex1")
        heatex2 = SimpleHeatExchanger("heatex1")
        c1 = Connection(self.inlet, "out1", heatex1, "in1", label="c1")
        c2 = Connection(heatex1, "out1", self.outlet, "in1", label="c2")

        c3 = Connection(self.inlet, "out2", heatex2, "in1", label="c3")
        c4 = Connection(heatex2, "out1", self.outlet, "in2", label="c4")

        self.add_conns(c1, c2, c3, c4)


def test_subsystem_duplicated_component_label():
    with raises(TESPyComponentError):
        SubsystemWithDuplicateComponentLabel("test")


def test_power_interface():
    """Subsystem with a power inlet passes E unchanged through SubsystemInterface."""

    class CompressorSubsystem(Subsystem):
        def __init__(self, label):
            self.num_in = 1
            self.num_out = 1
            self.num_power_in = 1
            super().__init__(label)

        def create_network(self):
            comp = Compressor("compressor")
            c1 = Connection(self.inlet, "out1", comp, "in1", label="c1")
            c2 = Connection(comp, "out1", self.outlet, "in1", label="c2")
            p1 = PowerConnection(
                self.inlet, "power_out1", comp, "power", label="p1"
            )
            self.add_conns(c1, c2, p1)


    nw = Network()
    nw.units.set_defaults(
        temperature="°C", pressure="bar"
    )
    source = Source("source")
    sink = Sink("sink")
    ps = PowerSource("power_source")
    sub = CompressorSubsystem("sub")

    c_in = Connection(source, "out1", sub.inlet, "in1", label="c_in")
    c_out = Connection(sub.outlet, "out1", sink, "in1", label="c_out")
    p_in = PowerConnection(ps, "power", sub.inlet, "power_in1", label="p_in")

    nw.add_conns(c_in, c_out, p_in)
    nw.add_subsystems(sub)

    c_in.set_attr(fluid={"air": 1}, T=20, p=1, m=1)
    sub.get_comp("compressor").set_attr(eta_s=0.8)
    c_out.set_attr(p=5)

    nw.solve("design")
    nw.assert_convergence()

    assert sub.inlet.power_inl[0].E.val_SI == approx(sub.inlet.power_outl[0].E.val_SI)


def test_heat_interface():
    """Subsystem with a heat inlet passes E unchanged through SubsystemInterface."""

    class HeatedFluidSubsystem(Subsystem):
        def __init__(self, label):
            self.num_in = 1
            self.num_out = 1
            self.num_heat_in = 1
            super().__init__(label)

        def create_network(self):
            hx = SimpleHeatExchanger("heat_exchanger")
            c1 = Connection(self.inlet, "out1", hx, "in1", label="c1")
            c2 = Connection(hx, "out1", self.outlet, "in1", label="c2")
            h1 = HeatConnection(
                self.inlet, "heat_out1", hx, "heat", label="h1"
            )
            self.add_conns(c1, c2, h1)

    nw = Network(iterinfo=False)
    nw.units.set_defaults(pressure="bar", temperature="degC")
    source = Source("source")
    sink = Sink("sink")
    hs = HeatSource("heat_source")
    sub = HeatedFluidSubsystem("sub")

    c_in = Connection(source, "out1", sub.inlet, "in1", label="c_in")
    c_out = Connection(sub.outlet, "out1", sink, "in1", label="c_out")
    h_in = HeatConnection(hs, "heat", sub.inlet, "heat_in1", label="h_in")

    nw.add_conns(c_in, c_out, h_in)
    nw.add_subsystems(sub)

    c_in.set_attr(fluid={"water": 1}, T=20, p=1, m=1)
    sub.get_comp("heat_exchanger").set_attr(dp=0)
    h_in.set_attr(E=10000)

    nw.solve("design")
    nw.assert_convergence()

    assert sub.inlet.heat_inl[0].E.val_SI == approx(sub.inlet.heat_outl[0].E.val_SI)


def test_del_subsystems():
    nw = Network()
    source = Source("source")
    sink1 = Sink("sink1")
    sink2 = Sink("sink2")

    sub = MySubsystem("my subsystem")

    c_in = Connection(source, "out1", sub.inlet, "in1", label="c_in")
    c_out1 = Connection(sub.outlet, "out1", sink1, "in1", label="c_out1")
    c_out2 = Connection(sub.outlet, "out2", sink2, "in1", label="c_out2")
    nw.add_conns(c_in, c_out1, c_out2)
    nw.add_subsystems(sub)

    assert "my subsystem" in nw.subsystems
    assert len(nw.conns) == 6  # 3 subsystem-internal + 3 external

    nw.del_subsystems(sub)

    assert "my subsystem" not in nw.subsystems
    # the three internal connections are gone; the three external ones remain
    assert len(nw.conns) == 3
    assert "c_in" in nw.conns.index
    assert "c_out1" in nw.conns.index
    assert "c_out2" in nw.conns.index
