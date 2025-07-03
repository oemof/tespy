# -*- coding: utf-8

"""Module for testing components of type merge.
This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_merge.py
SPDX-License-Identifier: MIT
"""

from pytest import approx

from tespy.components import Merge
from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Splitter
from tespy.connections import Connection
from tespy.networks import Network


class TestMerge:

    def setup_method(self):
        self.nwk = Network(T_unit="C", p_unit="bar", h_unit="kJ / kg")

        so1 = Source("Source1")
        so2 = Source("Source2")
        me = Merge("Merge")
        si = Sink("Sink")

        c1 = Connection(so1, "out1", me, "in1", label="1")
        c2 = Connection(so2, "out1", me, "in2", label="2")
        c3 = Connection(me, "out1", si, "in1", label="3")

        self.nwk.add_conns(c1, c2, c3)

    def test_single_fluid_at_outlet(self):

        c1, c2, c3 = self.nwk.get_conn(["1", "2", "3"])
        c1.set_attr(m=5, p=10, h=200)
        c2.set_attr(m=5, h=200)
        c3.set_attr(fluid={"water": 1})

        self.nwk.solve("design")
        self.nwk.assert_convergence()
        assert self.nwk.status == 0

        target = c1.m.val_SI + c2.m.val_SI
        msg = f"Target value for mass flow at connection 3 must be {target}."
        assert c3.m.val_SI == approx(target), msg

    def test_massflows_from_two_fluid_fractions(self):

        c1, c2, c3 = self.nwk.get_conn(["1", "2", "3"])
        c1.set_attr(m=5, p=10, h=200, fluid={"N2": 1})
        c2.set_attr(h=200, fluid={"O2": 1})
        c3.set_attr(fluid={"N2": 0.3, "O2": 0.7})

        self.nwk.solve("design")
        self.nwk.assert_convergence()

        target = c1.m.val / c3.fluid.val["N2"]
        msg = f"Target value for mass flow at connection 3 is {target}"
        assert c3.m.val_SI == approx(target), msg


class TestCyclicMerging:
    """
    Testing issue raised in https://github.com/oemof/tespy/issues/424
    """

    def setup_method(self):

        self.nwk = Network(T_unit="C", p_unit="bar", h_unit="kJ / kg")

        source = Source("source1")
        merge = Merge("merge")
        component1 = SimpleHeatExchanger("comp1", pr=1)
        splitter = Splitter("splitter")
        component2 = SimpleHeatExchanger("comp2")
        sink = Sink("sink")

        c1 = Connection(source, "out1", merge, "in1", label="1")
        c2 = Connection(merge, "out1", component1, "in1", label="2")
        c3 = Connection(component1, "out1", splitter, "in1", label="3")
        c4 = Connection(splitter, "out1", component2, "in1", label="4")
        c5 = Connection(component2, "out1", merge, "in2", label="5")
        c6 = Connection(splitter, "out2", sink, "in1", label="6")

        self.nwk.add_conns(c1, c2, c3, c4, c5, c6)

    def test_single_fluid_setup(self):

        c1, c3, c4, c5, c6 = self.nwk.get_conn(["1", "3", "4", "5", "6"])

        c1.set_attr(p=1, h=200, m=10, fluid={"R134a": 1})
        c3.set_attr(h=180)
        c4.set_attr(m=1)
        c5.set_attr(h=170)

        self.nwk.solve("design")
        self.nwk.assert_convergence()

        target = c1.m.val_SI
        msg = f"Target value for mass flow at connection 3 is {target}"
        assert c6.m.val_SI == approx(target), msg

    def test_two_fluid_setup(self):
        c1, c3, c4, c5, c6 = self.nwk.get_conn(["1", "3", "4", "5", "6"])

        c1.set_attr(p=1, h=200, m=10, fluid={"R134a": 1})
        c3.set_attr(h=180)
        c4.set_attr(m=1)
        c5.set_attr(h=170)

        self.nwk.solve("design")
        self.nwk.assert_convergence()

        target = c1.m.val_SI
        msg = f"Target value for mass flow at connection 3 is {target}"
        assert c6.m.val_SI == approx(target), msg
