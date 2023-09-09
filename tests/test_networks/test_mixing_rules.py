# -*- coding: utf-8

"""Module for testing mixing rule propagation in networks.
This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_mixing_rules.py
SPDX-License-Identifier: MIT
"""

from pytest import approx

from tespy.components import Compressor
from tespy.components import HeatExchangerSimple
from tespy.components import Merge
from tespy.components import Pipe
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Splitter
from tespy.components import Turbine
from tespy.connections import Connection
from tespy.networks import Network


class TestMixingRules:

    def setup_method(self):

        self.nwk = Network(T_unit="C", p_unit="bar", h_unit="kJ / kg")
        self.nwk.set_attr(iterinfo=False)

        so1 = Source("air")
        so2 = Source("Other gases")
        m = Merge("gas mixing")
        p = Pipe("Pipe", pr=1, Q=0)
        sp = Splitter("Splitter")
        t = Turbine("Turbine", pr=.1, eta_s=.8)
        cp = Compressor("Compressor", pr=10, eta_s=.8)
        si1 = Sink("Sink1")
        si2 = Sink("Sink2")

        c1 = Connection(so1, "out1", m, "in1", label="1")
        c2 = Connection(so2, "out1", m, "in2", label="2")
        c3 = Connection(m, "out1", p, "in1", label="3")
        c4 = Connection(p, "out1", sp, "in1", label="4")
        c5 = Connection(sp, "out1", t, "in1", label="5")
        c6 = Connection(t, "out1", si1, "in1", label="6")
        c7 = Connection(sp, "out2", cp, "in1", label="7")
        c8 = Connection(cp, "out1", si2, "in1", label="8")

        self.nwk.add_conns(c1, c2, c3, c4, c5, c6, c7, c8)

    def test_ideal_ideal_cond(self):

        c1, c2, c6, c7 = self.nwk.get_conn(["1", "2", "6", "7"])
        c1.set_attr(fluid={"N2": 0.76, "O2": 0.23, "Ar": 0.01}, m=10, T=400, p=1, mixing_rule="ideal-cond")
        c2.set_attr(fluid={"H2O": 1}, m=.5, T=400)
        c7.set_attr(m=4)

        p, t, cp = self.nwk.get_comp(["Pipe", "Turbine", "Compressor"])
        p.set_attr(pr=1, Q=0)
        t.set_attr(pr=.1, eta_s=.8)
        cp.set_attr(pr=10, eta_s=.8)

        self.nwk.solve("design")
        self.nwk._convergence_check()

        target = c2.m.val_SI / (c1.m.val_SI + c2.m.val_SI)
        msg = f"The H2O mass fraction in connection 7 must be {target}"
        assert c7.fluid.val["H2O"] == approx(target), msg

        h_ideal_cond = c6.h.val_SI
        for c in self.nwk.conns["object"]:
            c.mixing_rule = "ideal"

        self.nwk.solve("design")
        self.nwk._convergence_check()

        target = h_ideal_cond
        msg = f"The enthalpy at connection 6 must be equal to {target}"
        assert c6.h.val_SI == approx(target), msg

        c1.set_attr(T=200)
        c2.set_attr(T=200)

        self.nwk.solve("design")
        self.nwk._convergence_check()

        h_ideal = c6.h.val_SI
        for c in self.nwk.conns["object"]:
            c.mixing_rule = "ideal-cond"

        self.nwk.solve("design")
        self.nwk._convergence_check()

        target = h_ideal
        msg = (
            "Using ideal-cond mixing, the enthalpy at connection 6 must be "
            f"larger than using ideal mixing rule ({target})"
        )
        assert c6.h.val_SI > target, msg
