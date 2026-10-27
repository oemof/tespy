# -*- coding: utf-8

"""Tests of the Connection class: references, presolving and fluid handling.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_connections/test_connection.py

SPDX-License-Identifier: MIT
"""
from CoolProp.CoolProp import get_global_param_string
from pytest import approx
from pytest import mark

from tespy.components import Merge
from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.connections import Ref
from tespy.networks import Network
from tespy.tools.fluid_properties.functions import p_sat_TQ
from tespy.tools.fluid_properties.wrappers import FluidPropertyWrapper
from tespy.tools.units import SI_UNITS


class TestConnections:

    def setup_method(self):
        """Set up the model."""
        self.nw = Network()
        self.nw.units.set_defaults(**{
            "pressure": "bar",
            "pressure_difference": "bar",
            "temperature": "°C",
            "volumetric_flow": "l/s",
            "mass_flow": "t/h"
        })

        so1 = Source("source 1")
        so2 = Source("source 2")
        si1 = Sink("sink 1")
        si2 = Sink("sink 2")

        c1 = Connection(so1, "out1", si1, "in1", label="Some example label")
        c2 = Connection(so2, "out1", si2, "in1")

        self.nw.add_conns(c1, c2)

        c1.set_attr(m=1, p=1, T=25, fluid={"Air": 1})
        c2.set_attr(m=0.5, p=10, T=25, fluid={"Air": 1})

        self.nw.solve("design")

    def test_volumetric_flow_reference(self):
        """Test the referenced volumetric flow."""
        c1, c2 = self.nw.get_conn(
            ["Some example label", "source 2:out1_sink 2:in1"]
        )
        c2.set_attr(m=None, v=Ref(c1, 1, 0))
        self.nw.solve("design")

        m_expected = round(c1.m.val * c1.vol.val / c2.vol.val, 4)
        m_is = round(c2.m.val, 4)
        msg = (
            "The mass flow of the connection 2 should be equal to "
            f"{m_expected} kg/s, but is {m_is} kg/s"
        )
        assert m_is == m_expected, msg

        c2.set_attr(v=Ref(c1, 2, 10))
        self.nw.solve("design")

        v_expected = round(c1.v.val * 2 + 10, 4)
        v_is = round(c2.v.val, 4)
        msg = (
            "The mass flow of the connection 2 should be equal to "
            f"{v_expected} l/s, but is {v_is} l/s"
        )
        assert v_is == v_expected, msg

    def test_temperature_reference(self):
        """Test the referenced temperature."""
        c1, c2 = self.nw.get_conn(
            ["Some example label", "source 2:out1_sink 2:in1"]
        )
        c2.set_attr(T=None)
        c2.set_attr(T=Ref(c1, 1, 0))

        self.nw.solve("design")

        T_expected = round(c1.T.val, 4)
        T_is = round(c2.T.val, 4)
        msg = (
            "The temperature of the connection 2 should be equal to "
            f"{T_expected} C, but is {T_is} C"
        )
        assert T_is == T_expected, msg

        delta = -75
        c2.set_attr(T=Ref(c1, 1.5, delta))
        self.nw.solve("design")

        delta_SI = self.nw.units.ureg.Quantity(
            delta, self.nw.units.default["temperature_difference"]
        ).m_as(SI_UNITS["temperature_difference"])
        assert round(delta_SI, 4) == round(c2.T.reference.ref.delta_SI, 4)

        T_expected = round(c1.T.val_SI * 1.5 + delta_SI, 4)
        T_is = round(c2.T.val_SI, 4)
        msg = (
            "The temperature of the connection 2 should be equal to "
            f"{T_expected} C, but is {T_is} C"
        )
        assert T_is == T_expected, msg

    def test_primary_reference(self):
        """Test referenced primary variable."""
        c1, c2 = self.nw.get_conn(
            ["Some example label", "source 2:out1_sink 2:in1"]
        )
        c2.set_attr(m=None)
        c2.set_attr(m=Ref(c1, 1, 0))

        self.nw.solve("design")

        m_expected = round(c1.m.val, 4)
        m_is = round(c2.m.val, 4)
        msg = (
            "The mass flow of the connection 2 should be equal to "
            f"{m_expected} kg/s, but is {m_is} kg/s"
        )
        assert m_is == m_expected, msg

        delta = -0.5
        c2.set_attr(m=Ref(c1, 2, delta))
        self.nw.solve("design")
        delta_SI = self.nw.units.ureg.Quantity(
            delta, self.nw.units.default["mass_flow"]
        ).m_as(SI_UNITS["mass_flow"])
        assert round(delta_SI, 4) == round(c2.m.reference.ref.delta_SI, 4)

        m_expected = round(c1.m.val_SI * 2 + delta_SI, 4)
        m_is = round(c2.m.val_SI, 4)
        msg = (
            "The mass flow of the connection 2 should be equal to "
            f"{m_expected} kg/s, but is {m_is} kg/s"
        )
        assert m_is == m_expected, msg

    def test_ref_and_numerical_value_are_mutually_exclusive(self):
        """Setting a Ref clears any numerical value and vice versa."""
        c1, c2 = self.nw.get_conn(
            ["Some example label", "source 2:out1_sink 2:in1"]
        )
        # numerical value → Ref: numerical value must be unset
        c2.set_attr(m=0.5)
        c2.set_attr(m=Ref(c1, 1, 0))
        assert not c2.m.is_set
        assert c2.m.reference.is_set
        self.nw.solve("design")
        self.nw.assert_convergence()
        assert round(c2.m.val_SI, 4) == round(c1.m.val_SI, 4)

        # Ref → numerical value: Ref must be unset
        c2.set_attr(m=0.3)
        assert not c2.m.reference.is_set
        assert c2.m.is_set
        self.nw.solve("design")
        self.nw.assert_convergence()
        assert round(c2.m.val, 4) == round(0.3, 4)



_skipif_no_refprop = mark.skipif(
    get_global_param_string("REFPROP_version") == "n/a",
    reason="This test requires REFPROP.",
)

_ZEOTROPIC_FLUID = "ISOBUTAN[0.5]&IPENTANE[0.5]|mass"


@_skipif_no_refprop
class TestZeotropicConnectionPresolve:
    """Regression tests for the T+x presolve fix.

    Before the fix, _presolve used p_sat_T (Q=0.5) for all quality values.
    For zeotropic mixtures this gave the wrong pressure when x=0 or x=1.
    """

    def setup_method(self):
        self.nw = Network()
        self.nw.units.set_defaults(temperature="°C", pressure="bar", pressure_difference="bar")
        self.nw.iterinfo = False
        src = Source("src")
        snk = Sink("snk")
        pipe = SimpleHeatExchanger("pipe")
        self.c_in = Connection(src, "out1", pipe, "in1", label="c_in")
        self.c_out = Connection(pipe, "out1", snk, "in1", label="c_out")
        self.nw.add_conns(self.c_in, self.c_out)
        pipe.set_attr(dp=0, Q=0)

    def _solve(self, T_degC, x):
        self.c_in.set_attr(
            fluid={f"REFPROP::{_ZEOTROPIC_FLUID}": 1}, T=T_degC, x=x, m=1
        )
        self.nw.solve("design")
        self.nw.assert_convergence()

    def test_presolve_x1_gives_dew_pressure(self):
        self._solve(150, x=1)
        T_K = self.c_in.T.val_SI
        p_dew = p_sat_TQ(T_K, 1, self.c_in.fluid_data)
        p_mid = p_sat_TQ(T_K, 0.5, self.c_in.fluid_data)
        assert self.c_in.p.val_SI == approx(p_dew, rel=1e-4)
        assert self.c_in.p.val_SI != approx(p_mid, rel=1e-3)

    def test_presolve_x0_gives_bubble_pressure(self):
        self._solve(100, x=0)
        T_K = self.c_in.T.val_SI
        p_bubble = p_sat_TQ(T_K, 0, self.c_in.fluid_data)
        p_mid = p_sat_TQ(T_K, 0.5, self.c_in.fluid_data)
        assert self.c_in.p.val_SI == approx(p_bubble, rel=1e-4)
        assert self.c_in.p.val_SI != approx(p_mid, rel=1e-3)

    def test_presolve_intermediate_x_between_bubble_and_dew(self):
        self._solve(120, x=0.5)
        T_K = self.c_in.T.val_SI
        p_bubble = p_sat_TQ(T_K, 0, self.c_in.fluid_data)
        p_dew = p_sat_TQ(T_K, 1, self.c_in.fluid_data)
        p_lo, p_hi = sorted([p_bubble, p_dew])
        assert p_lo < self.c_in.p.val_SI < p_hi

    def test_calc_T_dew_exceeds_calc_T_bubble(self):
        self._solve(120, x=0.5)
        assert self.c_in.calc_T_dew() > self.c_in.calc_T_bubble()


def test_wrapper_kwargs_injection():
    so = Source("source")
    si = Sink("sink")
    c = Connection(so, "out1", si, "in1", label="c")
    c.set_attr(
        fluid={"H2O": 1},
        fluid_wrapper_kwargs={"H2O": {"testkeyword": "data"}},
        fluid_engines={"H2O": FluidPropertyWrapper}
    )
    # if kwargs could not be passed to the Wrapper instantiation this would
    # raise an error
    c._create_fluid_wrapper()


def test_fluid_balance_converges_with_variable_fractions():
    nw = Network()
    nw.units.set_defaults(pressure="bar", pressure_difference="bar", temperature="degC")
    merge = Merge("merge", num_in=3)
    c1 = Connection(Source("air"), "out1", merge, "in1", label="1")
    c2 = Connection(Source("carbon dioxide"), "out1", merge, "in2", label="2")
    c3 = Connection(Source("nitrogen"), "out1", merge, "in3", label="3")
    c4 = Connection(merge, "out1", Sink("sink"), "in1", label="4")
    nw.add_conns(c1, c2, c3, c4)

    c1.set_attr(
        fluid={"O2": 0.23, "N2": 0.76, "Ar": 0.01, "CO2": 0}, p=1, T=20, m=5
    )
    c2.set_attr(fluid={"CO2": 1}, T=20, m=5)
    c3.set_attr(fluid={"N2": 1}, T=20)
    c4.set_attr(fluid={"N2": 0.4}, fluid_balance=True)

    nw.solve("design")
    nw.assert_convergence()

    assert sum(c4.fluid.val.values()) == approx(1)
    assert c3.m.val_SI == approx(1 / 3)
