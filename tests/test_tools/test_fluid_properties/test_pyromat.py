from tespy.networks import Network
from tespy.connections import Connection
from tespy.components import Source
from tespy.components import Sink
from tespy.components import Turbine
from tespy.tools.fluid_properties.wrappers import PyromatWrapper


class TestPyromat:

    def setup_method(self):
        self.nwk = Network()

        so = Source("Source")
        tu = Turbine("Pump")
        si = Sink("Sink")

        c1 = Connection(so, "out1", tu, "in1", label="1")
        c2 = Connection(tu, "out1", si, "in1", label="2")

        self.nwk.add_conns(c1, c2)

        tu.set_attr(eta_s=0.9)

        c1.set_attr(v=1, p=1e5, T=500, fluid={"H2O": 1})
        c2.set_attr(p=1e4)

    def test_pyromat(self):
        c1, c2 = self.nwk.get_conn(["1", "2"])

        self.nwk.solve("design")
        self.nwk._convergence_check()

        h_out_ref = round(c2.h.val_SI / 1e3)
        T_out_ref = round(c2.T.val_SI)
        x_out_ref = round(c2.x.val_SI, 3)

        self.setup_method()
        c1, c2 = self.nwk.get_conn(["1", "2"])
        c1.set_attr(fluid={"mp::H2O": 1}, fluid_engines={"H2O": PyromatWrapper})

        self.nwk.solve("design")
        self.nwk._convergence_check()

        assert h_out_ref == round(c2.h.val_SI / 1e3)
        assert T_out_ref == round(c2.T.val_SI)
        assert x_out_ref == round(c2.x.val_SI, 3)
