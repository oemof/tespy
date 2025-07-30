from tespy.components import Splitter
from tespy.components import Subsystem
from tespy.connections import Connection


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
