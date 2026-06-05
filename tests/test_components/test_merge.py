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
from tespy.connections import Ref
from tespy.networks import Network


class TestMerge:

    def setup_method(self):
        self.nwk = Network()
        self.nwk.units.set_defaults(**{
            "pressure": "bar", "pressure_difference": "bar",
            "temperature": "degC", "enthalpy": "kJ/kg"
        })

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


class TestMergeNumEq:
    """Cover the three branches of Merge._update_num_eq."""

    def setup_method(self):
        self.nwk = Network(iterinfo=False)
        so1 = Source("so1")
        so2 = Source("so2")
        self.merge = Merge("merge")
        si = Sink("si")
        self.c1 = Connection(so1, "out1", self.merge, "in1", label="1")
        self.c2 = Connection(so2, "out1", self.merge, "in2", label="2")
        self.c3 = Connection(self.merge, "out1", si, "in1", label="3")
        self.nwk.add_conns(self.c1, self.c2, self.c3)

    def test_branch1_fluid_propagated_not_set_on_merge_connections(self):
        """Branch 1: fluid propagated from upstream; merge connections carry no
        user-specified fractions and no variable fractions.

        The fluid is set only on the source-to-SHE connections. The SHE passes
        it through, but the resulting SHE-to-merge connections carry the fluid
        purely via propagation: is_set and is_var remain empty on all three of
        the merge's connections. The merge therefore contributes only the mass
        flow balance equation (fluid_eq=0, mass_flow num_eq=1).
        """
        nwk = Network(iterinfo=False)
        so1, so2, si = Source("so1"), Source("so2"), Sink("si")
        she1 = SimpleHeatExchanger("she1", pr=1, Q=0)
        she2 = SimpleHeatExchanger("she2", pr=1, Q=0)
        merge = Merge("merge")

        c1 = Connection(so1, "out1", she1, "in1", label="1")
        c2 = Connection(she1, "out1", merge, "in1", label="2")
        c3 = Connection(so2, "out1", she2, "in1", label="3")
        c4 = Connection(she2, "out1", merge, "in2", label="4")
        c5 = Connection(merge, "out1", si, "in1", label="5")
        nwk.add_conns(c1, c2, c3, c4, c5)

        c1.set_attr(fluid={"Air": 1}, p=1e5, T=300, m=1)
        c3.set_attr(fluid={"Air": 1}, T=350, m=0.5)

        nwk.solve("design")
        nwk.assert_convergence()

        assert all(len(c.fluid.is_set) == 0 for c in merge.inl + merge.outl)
        assert merge.variable_fluids == set()
        assert merge.constraints["fluid_constraints"].num_eq == 0
        assert merge.constraints["mass_flow_constraints"].num_eq == 1
        assert c5.m.val_SI == approx(1.5)

    def test_branch2_all_fluids_set(self):
        """Branch 2: all fluid fractions specified directly on merge connections,
        none variable.

        With a single pure fluid fixed everywhere the fluid balance equation
        already implies mass conservation, so mass_flow_constraints is
        suppressed (num_eq=0) and fluid_constraints provides one equation per
        fluid.
        """
        self.c1.set_attr(fluid={"Air": 1}, p=1e5, T=300, m=1)
        self.c2.set_attr(fluid={"Air": 1}, T=350, m=0.5)

        self.nwk.solve("design")
        self.nwk.assert_convergence()

        assert self.merge.variable_fluids == set()
        assert self.merge.constraints["mass_flow_constraints"].num_eq == 0
        assert self.merge.constraints["fluid_constraints"].num_eq == 1
        assert self.c3.m.val_SI == approx(1.5)

    def test_branch3_variable_fluid_fractions(self):
        """Branch 3: outlet fluid fractions are solver variables.

        Merging two different pure fluids means the outlet fractions are
        unknown. fluid_constraints provides one equation per variable fluid;
        mass_flow_constraints stays at its default of 1.
        """
        self.c1.set_attr(fluid={"O2": 1}, p=1e5, T=300, m=1)
        self.c2.set_attr(fluid={"N2": 1}, T=300, m=0.5)

        self.nwk.solve("design")
        self.nwk.assert_convergence()

        assert len(self.merge.variable_fluids) == 2
        assert self.merge.constraints["fluid_constraints"].num_eq == 2
        assert self.merge.constraints["mass_flow_constraints"].num_eq == 1
        assert self.c3.m.val_SI == approx(1.5)
        assert self.c3.fluid.val["O2"] == approx(1 / 1.5, rel=1e-4)
        assert self.c3.fluid.val["N2"] == approx(0.5 / 1.5, rel=1e-4)

    def test_branch3_ref_in2_to_in1(self):
        """Branch 3: N2 inlet (in2) mass flow referenced to O2 inlet (in1).

        SHE on the O2 path fixes the absolute O2 mass flow via Q and T_out.
        N2 mass flow = 0.5 * O2 mass flow → outlet fractions O2:N2 = 2:1.
        """
        nwk = Network(iterinfo=False)
        so_o2, so_n2, si = Source("so_O2"), Source("so_N2"), Sink("si")
        she = SimpleHeatExchanger("she", pr=1)
        merge = Merge("merge")
        c1 = Connection(so_o2, "out1", she, "in1", label="1")
        c2 = Connection(she, "out1", merge, "in1", label="2")
        c3 = Connection(so_n2, "out1", merge, "in2", label="3")
        c4 = Connection(merge, "out1", si, "in1", label="4")
        nwk.add_conns(c1, c2, c3, c4)
        c1.set_attr(fluid={"O2": 1}, p=1e5, T=300)
        c2.set_attr(T=320)
        she.set_attr(Q=5000)
        c3.set_attr(fluid={"N2": 1}, T=300, m=Ref(c2, 0.5, 0))
        nwk.solve("design")
        nwk.assert_convergence()
        assert c3.m.val_SI == approx(c2.m.val_SI * 0.5, rel=1e-4)
        assert c4.fluid.val["O2"] == approx(2 / 3, rel=1e-3)
        assert c4.fluid.val["N2"] == approx(1 / 3, rel=1e-3)
        assert she.Q.val_SI == approx(c2.m.val_SI * (c2.h.val_SI - c1.h.val_SI), rel=1e-4)

    def test_branch3_ref_in1_to_in2(self):
        """Branch 3: O2 inlet (in1) mass flow referenced to N2 inlet (in2).

        SHE on the N2 path fixes the absolute N2 mass flow via Q and T_out.
        O2 mass flow = 2 * N2 mass flow → outlet fractions O2:N2 = 2:1.
        """
        nwk = Network(iterinfo=False)
        so_o2, so_n2, si = Source("so_O2"), Source("so_N2"), Sink("si")
        she = SimpleHeatExchanger("she", pr=1)
        merge = Merge("merge")
        c1 = Connection(so_n2, "out1", she, "in1", label="1")
        c2 = Connection(so_o2, "out1", merge, "in1", label="2")
        c3 = Connection(she, "out1", merge, "in2", label="3")
        c4 = Connection(merge, "out1", si, "in1", label="4")
        nwk.add_conns(c1, c2, c3, c4)
        c1.set_attr(fluid={"N2": 1}, p=1e5, T=300)
        c3.set_attr(T=320)
        she.set_attr(Q=5000)
        c2.set_attr(fluid={"O2": 1}, T=300, m=Ref(c3, 2.0, 0))
        nwk.solve("design")
        nwk.assert_convergence()
        assert c2.m.val_SI == approx(c3.m.val_SI * 2.0, rel=1e-4)
        assert c4.fluid.val["O2"] == approx(2 / 3, rel=1e-3)
        assert c4.fluid.val["N2"] == approx(1 / 3, rel=1e-3)
        assert she.Q.val_SI == approx(c3.m.val_SI * (c3.h.val_SI - c1.h.val_SI), rel=1e-4)

    def test_branch3_ref_in1_to_out(self):
        """Branch 3: O2 inlet (in1) mass flow referenced to the outlet.

        SHE on the N2 path fixes the absolute N2 mass flow.
        O2 mass flow = (2/3) * outlet mass flow; from the mass balance this
        implies O2 = 2 * N2 → outlet fractions O2:N2 = 2:1.
        """
        nwk = Network(iterinfo=False)
        so_o2, so_n2, si = Source("so_O2"), Source("so_N2"), Sink("si")
        she = SimpleHeatExchanger("she", pr=1)
        merge = Merge("merge")
        c1 = Connection(so_n2, "out1", she, "in1", label="1")
        c2 = Connection(so_o2, "out1", merge, "in1", label="2")
        c3 = Connection(she, "out1", merge, "in2", label="3")
        c4 = Connection(merge, "out1", si, "in1", label="4")
        nwk.add_conns(c1, c2, c3, c4)
        c1.set_attr(fluid={"N2": 1}, p=1e5, T=300)
        c3.set_attr(T=320)
        she.set_attr(Q=5000)
        c2.set_attr(fluid={"O2": 1}, T=300, m=Ref(c4, 2 / 3, 0))
        nwk.solve("design")
        nwk.assert_convergence()
        assert c2.m.val_SI == approx(c4.m.val_SI * 2 / 3, rel=1e-4)
        assert c4.fluid.val["O2"] == approx(2 / 3, rel=1e-3)
        assert c4.fluid.val["N2"] == approx(1 / 3, rel=1e-3)

    def test_branch3_ref_in2_to_out(self):
        """Branch 3: N2 inlet (in2) mass flow referenced to the outlet.

        SHE on the O2 path fixes the absolute O2 mass flow.
        N2 mass flow = (1/3) * outlet mass flow; from the mass balance this
        implies N2 = 0.5 * O2 → outlet fractions O2:N2 = 2:1.
        """
        nwk = Network(iterinfo=False)
        so_o2, so_n2, si = Source("so_O2"), Source("so_N2"), Sink("si")
        she = SimpleHeatExchanger("she", pr=1)
        merge = Merge("merge")
        c1 = Connection(so_o2, "out1", she, "in1", label="1")
        c2 = Connection(she, "out1", merge, "in1", label="2")
        c3 = Connection(so_n2, "out1", merge, "in2", label="3")
        c4 = Connection(merge, "out1", si, "in1", label="4")
        nwk.add_conns(c1, c2, c3, c4)
        c1.set_attr(fluid={"O2": 1}, p=1e5, T=300)
        c2.set_attr(T=320)
        she.set_attr(Q=5000)
        c3.set_attr(fluid={"N2": 1}, T=300, m=Ref(c4, 1 / 3, 0))
        nwk.solve("design")
        nwk.assert_convergence()
        assert c3.m.val_SI == approx(c4.m.val_SI * 1 / 3, rel=1e-4)
        assert c4.fluid.val["O2"] == approx(2 / 3, rel=1e-3)
        assert c4.fluid.val["N2"] == approx(1 / 3, rel=1e-3)


class TestCyclicMerging:
    """
    Testing issue raised in https://github.com/oemof/tespy/issues/424
    """

    def setup_method(self):

        self.nwk = Network()
        self.nwk.units.set_defaults(**{
            "pressure": "bar", "pressure_difference": "bar",
            "temperature": "degC", "enthalpy": "kJ/kg"
        })

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
