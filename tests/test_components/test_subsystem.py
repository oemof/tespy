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


def _compressor_network():
    nw = Network(iterinfo=False)
    nw.units.set_defaults(pressure="bar", temperature="degC")
    so = Source("air inlet")
    cp = Compressor("compressor")
    si = Sink("air outlet")
    ps = PowerSource("grid")

    c1 = Connection(so, "out1", cp, "in1", label="1")
    c2 = Connection(cp, "out1", si, "in1", label="2")
    p1 = PowerConnection(ps, "power", cp, "power", label="p1")
    nw.add_conns(c1, c2, p1)

    c1.set_attr(fluid={"air": 1}, T=20, p=1, m=1)
    c2.set_attr(p=5)
    cp.set_attr(eta_s=0.8)
    return nw


def test_subsystem_from_network():
    nw = _compressor_network()
    nw.solve("design")
    nw.assert_convergence()

    sub = Subsystem.from_network("compressor block", nw)

    assert sub.interface_map == {
        "in1": "air inlet", "out1": "air outlet", "power_in1": "grid"
    }
    assert sub.get_comp("inlet").label == "compressor block_inlet"
    assert sub.get_comp("compressor").label == "compressor block_compressor"
    assert sub.get_conn("1").label == "compressor block_1"
    # all specifications are carried over unchanged
    assert sub.get_conn("1").T.is_set
    assert sub.get_conn("1").fluid.is_set
    assert sub.get_conn("2").p.is_set
    assert sub.get_comp("compressor").eta_s.is_set
    # the source network is untouched
    assert nw.get_conn("1").T.is_set
    assert "air inlet" in nw.comps.index

    # the boundary conditions travel with the subsystem, so the external
    # connections stay unspecified
    nw2 = Network(iterinfo=False)
    nw2.units.set_defaults(pressure="bar", temperature="degC")
    so = Source("source")
    si = Sink("sink")
    ps = PowerSource("power source")

    c_in = Connection(so, "out1", sub.inlet, "in1", label="c_in")
    c_out = Connection(sub.outlet, "out1", si, "in1", label="c_out")
    p_in = PowerConnection(ps, "power", sub.inlet, "power_in1", label="p_in")
    nw2.add_conns(c_in, c_out, p_in)
    nw2.add_subsystems(sub)

    nw2.solve("design")
    nw2.assert_convergence()

    assert p_in.E.val_SI == approx(nw.get_conn("p1").E.val_SI)


def test_subsystem_from_network_structure_only():
    """A never solved, unparametrized network converts to a subsystem."""
    nw = Network(iterinfo=False)
    so = Source("air inlet")
    cp = Compressor("compressor")
    si = Sink("air outlet")
    ps = PowerSource("grid")

    c1 = Connection(so, "out1", cp, "in1", label="1")
    c2 = Connection(cp, "out1", si, "in1", label="2")
    p1 = PowerConnection(ps, "power", cp, "power", label="p1")
    nw.add_conns(c1, c2, p1)

    sub = Subsystem.from_network("block", nw)

    assert sub.interface_map == {
        "in1": "air inlet", "out1": "air outlet", "power_in1": "grid"
    }

    nw2 = Network(iterinfo=False)
    nw2.units.set_defaults(pressure="bar", temperature="degC")
    so = Source("source")
    si = Sink("sink")
    ps = PowerSource("power source")

    c_in = Connection(so, "out1", sub.inlet, "in1", label="c_in")
    c_out = Connection(sub.outlet, "out1", si, "in1", label="c_out")
    p_in = PowerConnection(ps, "power", sub.inlet, "power_in1", label="p_in")
    nw2.add_conns(c_in, c_out, p_in)
    nw2.add_subsystems(sub)

    c_in.set_attr(fluid={"air": 1}, T=20, p=1, m=1)
    c_out.set_attr(p=5)
    sub.get_comp("compressor").set_attr(eta_s=0.8)

    nw2.solve("design")
    nw2.assert_convergence()

    assert p_in.E.val_SI > 0


def test_subsystem_from_network_interface_exceptions():
    nw = _compressor_network()
    sub = Subsystem.from_network(
        "sub", nw, interface_exceptions=["air inlet", "grid"]
    )

    assert sub.num_in == 0
    assert sub.num_out == 1
    assert sub.num_power_in == 0
    assert sub.interface_map == {"out1": "air outlet"}
    assert sub.get_comp("air inlet").label == "sub_air inlet"
    assert sub.get_comp("grid").label == "sub_grid"


def test_subsystem_from_network_repeatable():
    nw = _compressor_network()
    sub1 = Subsystem.from_network("block1", nw)

    nw_dict = nw.export()
    sub2 = Subsystem.from_network("block2", nw)

    assert sub1.get_comp("compressor").label == "block1_compressor"
    assert sub2.get_comp("compressor").label == "block2_compressor"
    assert nw_dict["Component"]["Source"] == nw.export()["Component"]["Source"]


def test_subsystem_from_network_rewiring_covers_registry():
    """Every registered connection class must be handled by from_network."""
    from tespy.components.subsystem import _BOUNDARY_REWIRING
    from tespy.connections.connection import connection_registry

    covered = {
        conn_key
        for conn_keys, _, _ in _BOUNDARY_REWIRING
        for conn_key in conn_keys
    }
    missing = set(connection_registry.items) - covered
    msg = (
        f"The connection classes {missing} are not handled by "
        "Subsystem.from_network. Please add them to the _BOUNDARY_REWIRING "
        "table in tespy/components/subsystem.py."
    )
    assert not missing, msg


def test_subsystem_from_network_reserved_label():
    nw = Network(iterinfo=False)
    so = Source("inlet")
    si = Sink("sink")
    c1 = Connection(so, "out1", si, "in1", label="1")
    nw.add_conns(c1)

    with raises(TESPyComponentError):
        Subsystem.from_network("sub", nw)


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
