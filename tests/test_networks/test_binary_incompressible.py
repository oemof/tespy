# -*- coding: utf-8

"""Module for testing mixing rule propagation in networks.
This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_binary_incompressible.py
SPDX-License-Identifier: MIT
"""

from tespy.components import HeatExchanger
from tespy.components import Pump
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.networks import Network


class TestBinaryIncompressibles:


    def setup_method(self):
        self.nw = Network()

        so1 = Source("Source1")
        so2 = Source("Source2")
        pu = Pump("Pump")
        he = HeatExchanger("Heat")
        si1 = Sink("Sink1")
        si2 = Sink("Sink2")

        c1 = Connection(so1, "out1", pu, "in1", label="1")
        c2 = Connection(pu, "out1", he, "in2", label="2")
        c3 = Connection(he, "out2", si1, "in1", label="3")

        c4 = Connection(so2, "out1", he, "in1", label="4")
        c5 = Connection(he, "out1", si2, "in1", label="5")

        self.nw.add_conns(c1, c2, c3, c4, c5)

        c1.set_attr(v=1, p=1e5, T=300, fluid={"INCOMP::MPG[0.2]": 1})
        c2.set_attr(h0=2e4)
        c3.set_attr(p=1e5, T=320, h0=1e5)

        he.set_attr(pr1=0.98, pr2=0.98, ttd_l=10)
        pu.set_attr(eta_s=0.7)

        c4.set_attr(p=1e5, T=350, fluid={"INCOMP::MEG[0.2]": 1})
        c5.set_attr(h0=1e5)

        self.nw.solve("design")

    def test_binaries(self):
        self.nw._convergence_check()
