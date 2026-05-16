from pytest import raises

from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Splitter
from tespy.components import Subsystem
from tespy.connections import Connection
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
