from pytest import raises

from tespy.components import Sink
from tespy.components import Source
from tespy.components import Splitter
from tespy.components import Subsystem
from tespy.connections import Connection
from tespy.networks import Network


class SplitterSubsystem(Subsystem):
    def __init__(self, label):
        self.num_in = 1
        self.num_out = 2
        super().__init__(label)

    def create_network(self):
        splitter = Splitter("splitter", num_out=2)
        self.add_conns(
            Connection(self.inlet, "out1", splitter, "in1", label="feed"),
            Connection(splitter, "out1", self.outlet, "in1", label="first"),
            Connection(splitter, "out2", self.outlet, "in2", label="second"),
        )


def _network_with_splitter_subsystem():
    nw = Network(iterinfo=False)
    sub = SplitterSubsystem("sub")
    c_in = Connection(Source("source"), "out1", sub, "in1", label="c_in")
    c_out1 = Connection(sub, "out1", Sink("sink1"), "in1", label="c_out1")
    c_out2 = Connection(sub, "out2", Sink("sink2"), "in1", label="c_out2")
    nw.add_subsystems(sub)
    nw.add_conns(c_in, c_out1, c_out2)
    return nw, sub


def test_get_conn_miss_raises():
    sub = SplitterSubsystem("sub")
    with raises(KeyError, match="sub_feed"):
        sub.get_conn("sub_feed")


def test_get_comp_miss_raises():
    sub = SplitterSubsystem("sub")
    with raises(KeyError, match="sub_splitter"):
        sub.get_comp("sub_splitter")


def test_del_subsystems_with_boundary_leaves_valid_network():
    nw, sub = _network_with_splitter_subsystem()
    source = Source("other source")
    sink = Sink("other sink")
    nw.add_conns(Connection(source, "out1", sink, "in1", label="other"))

    nw.del_subsystems(sub, with_boundary=True)

    assert list(nw.conns.index) == ["other"]
    assert list(nw.comps.index) == ["other sink", "other source"]
    nw.check_topology()
