# -*- coding: utf-8

"""Module for testing fluid properties of gas mixtures.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_tools/test_fluid_properties/test_iapws.py

SPDX-License-Identifier: MIT
"""
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Turbine
from tespy.connections import Connection
from tespy.networks import Network
from tespy.tools.fluid_properties.wrappers import IAPWSWrapper


class TestIAPWS:

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

    def test_iapws_95(self):
        c1, c2 = self.nwk.get_conn(["1", "2"])

        self.nwk.solve("design")
        self.nwk.assert_convergence()

        h_out_ref = round(c2.h.val_SI, 3)
        T_out_ref = round(c2.T.val_SI, 3)
        x_out_ref = round(c2.x.val_SI, 3)

        self.setup_method()
        c1, c2 = self.nwk.get_conn(["1", "2"])
        c1.set_attr(fluid={"IF95::H2O": 1}, fluid_engines={"H2O": IAPWSWrapper})

        self.nwk.solve("design")
        self.nwk.assert_convergence()

        assert h_out_ref == round(c2.h.val_SI, 3)
        assert T_out_ref == round(c2.T.val_SI, 3)
        assert x_out_ref == round(c2.x.val_SI, 3)

    def test_iapws_97(self):
        c1, c2 = self.nwk.get_conn(["1", "2"])
        c1.set_attr(fluid={"IF97::H2O": 1})

        self.nwk.solve("design")
        self.nwk.assert_convergence()

        h_out_ref = round(c2.h.val_SI / 1000)
        T_out_ref = round(c2.T.val_SI)
        x_out_ref = round(c2.x.val_SI, 3)

        self.setup_method()
        c1, c2 = self.nwk.get_conn(["1", "2"])
        c1.set_attr(fluid={"IF97::H2O": 1}, fluid_engines={"H2O": IAPWSWrapper})

        self.nwk.solve("design")
        self.nwk.assert_convergence()

        assert h_out_ref == round(c2.h.val_SI / 1000)
        assert T_out_ref == round(c2.T.val_SI)
        assert x_out_ref == round(c2.x.val_SI, 3)
