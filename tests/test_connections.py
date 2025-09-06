# -*- coding: utf-8

"""Module for testing busses.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_connections.py

SPDX-License-Identifier: MIT
"""
from CoolProp.CoolProp import get_global_param_string
from pytest import approx
from pytest import fixture
from pytest import mark

from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.connections import Ref
from tespy.networks import Network
from tespy.tools.fluid_properties.functions import T_bubble_p
from tespy.tools.fluid_properties.functions import T_dew_p
from tespy.tools.units import SI_UNITS


class TestConnections:

    def setup_method(self):
        """Set up the model."""
        self.nw = Network()
        self.nw.units.set_defaults(**{
            "pressure": "bar",
            # TODO: replace C with degC in next major version
            "temperature": "C",
            "volumetric_flow": "l/s",
            "mass_flow": "t/h"
        })

        so1 = Source('source 1')
        so2 = Source('source 2')
        si1 = Sink('sink 1')
        si2 = Sink('sink 2')

        c1 = Connection(so1, 'out1', si1, 'in1', label='Some example label')
        c2 = Connection(so2, 'out1', si2, 'in1')

        self.nw.add_conns(c1, c2)

        c1.set_attr(m=1, p=1, T=25, fluid={'Air': 1})
        c2.set_attr(m=0.5, p=10, T=25, fluid={'Air': 1})

        self.nw.solve('design')

    def test_volumetric_flow_reference(self):
        """Test the referenced volumetric flow."""
        c1, c2 = self.nw.get_conn(
            ['Some example label', 'source 2:out1_sink 2:in1']
        )
        c2.set_attr(m=None, v=Ref(c1, 1, 0))
        self.nw.solve('design')

        m_expected = round(c1.m.val * c1.vol.val / c2.vol.val, 4)
        m_is = round(c2.m.val, 4)
        msg = (
            'The mass flow of the connection 2 should be equal to '
            f'{m_expected} kg/s, but is {m_is} kg/s'
        )
        assert m_is == m_expected, msg

        c2.set_attr(v=Ref(c1, 2, 10))
        self.nw.solve('design')

        v_expected = round(c1.v.val * 2 + 10, 4)
        v_is = round(c2.v.val, 4)
        msg = (
            'The mass flow of the connection 2 should be equal to '
            f'{v_expected} l/s, but is {v_is} l/s'
        )
        assert v_is == v_expected, msg

    def test_temperature_reference(self):
        """Test the referenced temperature."""
        c1, c2 = self.nw.get_conn(
            ['Some example label', 'source 2:out1_sink 2:in1']
        )
        c2.set_attr(T=None)
        c2.set_attr(T=Ref(c1, 1, 0))

        self.nw.solve('design')

        T_expected = round(c1.T.val, 4)
        T_is = round(c2.T.val, 4)
        msg = (
            'The temperature of the connection 2 should be equal to '
            f'{T_expected} C, but is {T_is} C'
        )
        assert T_is == T_expected, msg

        delta = -75
        c2.set_attr(T=Ref(c1, 1.5, delta))
        self.nw.solve('design')

        delta_SI = self.nw.units.ureg.Quantity(
            delta, self.nw.units.default["temperature_difference"]
        ).to(SI_UNITS["temperature_difference"]).magnitude
        assert round(delta_SI, 4) == round(c2.T_ref.ref.delta_SI, 4)

        T_expected = round(c1.T.val_SI * 1.5 + delta_SI, 4)
        T_is = round(c2.T.val_SI, 4)
        msg = (
            'The temperature of the connection 2 should be equal to '
            f'{T_expected} C, but is {T_is} C'
        )
        assert T_is == T_expected, msg

    def test_primary_reference(self):
        """Test referenced primary variable."""
        c1, c2 = self.nw.get_conn(
            ['Some example label', 'source 2:out1_sink 2:in1']
        )
        c2.set_attr(m=None)
        c2.set_attr(m=Ref(c1, 1, 0))

        self.nw.solve('design')

        m_expected = round(c1.m.val, 4)
        m_is = round(c2.m.val, 4)
        msg = (
            'The mass flow of the connection 2 should be equal to '
            f'{m_expected} kg/s, but is {m_is} kg/s'
        )
        assert m_is == m_expected, msg

        delta = -0.5
        c2.set_attr(m=Ref(c1, 2, delta))
        self.nw.solve('design')
        delta_SI = self.nw.units.ureg.Quantity(
            delta, self.nw.units.default["mass_flow"]
        ).to(SI_UNITS["mass_flow"]).magnitude
        assert round(delta_SI, 4) == round(c2.m_ref.ref.delta_SI, 4)

        m_expected = round(c1.m.val_SI * 2 + delta_SI, 4)
        m_is = round(c2.m.val_SI, 4)
        msg = (
            'The mass flow of the connection 2 should be equal to '
            f'{m_expected} kg/s, but is {m_is} kg/s'
        )
        assert m_is == m_expected, msg


@fixture
def simple_test_network():
    nw = Network()
    nw.units.set_defaults(
        temperature="degC",
        pressure="bar"
    )

    so = Source("source")
    si = Sink("sink")

    heatexchanger = SimpleHeatExchanger("heatexchanger")

    c1 = Connection(so, "out1", heatexchanger, "in1", label="c1")
    c2 = Connection(heatexchanger, "out1", si, "in1", label="c2")

    nw.add_conns(c1, c2)
    return nw


@mark.skipif(
    get_global_param_string("REFPROP_version") == "n/a",
    reason='This test requires REFPROP, dependency is missing.'
)
def test_td_bubble_and_td_dew_in_iterations(simple_test_network):
    nw = simple_test_network
    c1, c2 = nw.get_conn(["c1", "c2"])
    heatexchanger = nw.get_comp("heatexchanger")
    # R513A
    c1.set_attr(m=1, p=10, fluid={"REFPROP::R134A[0.44]&R1234yf[0.56]|mass": 1})
    delta_T = 10
    c2.set_attr(td_bubble=delta_T)

    # settings to prevent preprocessing of temperatures
    heatexchanger.set_attr(Q=1e5, zeta=0)

    nw.solve("design")

    nw.assert_convergence()

    T_with_td_bubble_set = c2.T.val_SI

    c2.set_attr(td_dew=-delta_T, td_bubble=None)
    nw.solve("design")

    nw.assert_convergence()

    T_with_td_dew_set = c2.T.val_SI
    T_dew = T_dew_p(c2.p.val_SI, c2.fluid_data)
    T_bubble = T_bubble_p(c2.p.val_SI, c2.fluid_data)

    assert approx(T_with_td_bubble_set) == T_bubble - delta_T
    # - delta_T because dew line temperature is set negative here!
    assert approx(T_with_td_dew_set) == T_dew - delta_T
    assert T_with_td_bubble_set != T_with_td_dew_set


def test_td_bubble_larger_0(simple_test_network):
    nw = simple_test_network
    c1, c2 = nw.get_conn(["c1", "c2"])
    heatexchanger = nw.get_comp("heatexchanger")

    c1.set_attr(m=1, p=10, fluid={"R290": 1})
    delta_T = 10
    c2.set_attr(td_bubble=delta_T)

    # settings to prevent preprocessing of temperatures
    heatexchanger.set_attr(Q=1e5, zeta=0)
    nw.solve("design")
    nw.assert_convergence()

    assert approx(c2.T.val_SI + delta_T) == T_bubble_p(c2.p.val_SI, c2.fluid_data)


def test_td_bubble_equals_0(simple_test_network):
    nw = simple_test_network
    c1, c2 = nw.get_conn(["c1", "c2"])
    heatexchanger = nw.get_comp("heatexchanger")

    c1.set_attr(m=1, p=10, fluid={"R290": 1})
    delta_T = 0
    c2.set_attr(td_bubble=delta_T)

    # settings to prevent preprocessing of temperatures
    heatexchanger.set_attr(Q=1e5, zeta=0)
    nw.solve("design")
    nw.assert_convergence()

    assert approx(c2.T.val_SI + delta_T) == T_bubble_p(c2.p.val_SI, c2.fluid_data)


def test_td_bubble_in_preprocessing1(simple_test_network):
    nw = simple_test_network
    c1, c2 = nw.get_conn(["c1", "c2"])
    heatexchanger = nw.get_comp("heatexchanger")

    c1.set_attr(m=1, p=10, fluid={"R290": 1})
    delta_T = 10
    c2.set_attr(td_bubble=delta_T)

    # settings to prevent preprocessing of temperatures
    heatexchanger.set_attr(Q=1e5, dp=0)
    nw.solve("design")
    nw.assert_convergence()

    assert approx(c2.T.val_SI + delta_T) == T_bubble_p(c2.p.val_SI, c2.fluid_data)


def test_td_bubble_in_preprocessing2(simple_test_network):
    nw = simple_test_network
    c1, c2 = nw.get_conn(["c1", "c2"])
    heatexchanger = nw.get_comp("heatexchanger")

    c1.set_attr(m=1, fluid={"R290": 1})
    delta_T = 10
    c2.set_attr(T=50, td_bubble=delta_T)

    # settings to prevent preprocessing of temperatures
    heatexchanger.set_attr(Q=1e5, dp=0)
    nw.solve("design")
    nw.assert_convergence()

    assert approx(c2.T.val_SI + delta_T) == T_bubble_p(c2.p.val_SI, c2.fluid_data)


def test_td_bubble_equals_0_in_preprocessing1(simple_test_network):
    nw = simple_test_network
    c1, c2 = nw.get_conn(["c1", "c2"])
    heatexchanger = nw.get_comp("heatexchanger")

    c1.set_attr(m=1, p=10, fluid={"R290": 1})
    delta_T = 0
    c2.set_attr(td_bubble=delta_T)

    # settings to prevent preprocessing of temperatures
    heatexchanger.set_attr(Q=1e5, dp=0)
    nw.solve("design")
    nw.assert_convergence()

    assert approx(c2.T.val_SI + delta_T) == T_bubble_p(c2.p.val_SI, c2.fluid_data)


def test_td_bubble_equals_0_in_preprocessing2(simple_test_network):
    nw = simple_test_network
    c1, c2 = nw.get_conn(["c1", "c2"])
    heatexchanger = nw.get_comp("heatexchanger")

    c1.set_attr(m=1, fluid={"R290": 1})
    delta_T = 0
    c2.set_attr(T=50, td_bubble=delta_T)

    # settings to prevent preprocessing of temperatures
    heatexchanger.set_attr(Q=1e5, dp=0)
    nw.solve("design")
    nw.assert_convergence()

    assert approx(c2.T.val_SI + delta_T) == T_bubble_p(c2.p.val_SI, c2.fluid_data)


def test_td_dew_larger_0(simple_test_network):
    nw = simple_test_network
    c1, c2 = nw.get_conn(["c1", "c2"])
    heatexchanger = nw.get_comp("heatexchanger")

    c1.set_attr(m=1, p=10, fluid={"R290": 1})
    delta_T = 10.0
    # setting zeta does not work here, not sure why
    c2.set_attr(td_dew=delta_T, v=Ref(c1, 1, 0))

    # settings to prevent preprocessing of temperatures
    heatexchanger.set_attr(Q=0)
    nw.solve("design")
    nw.assert_convergence()

    # minus delta T, because td_dew is higher than T_dew
    assert approx(c2.T.val_SI - delta_T) == T_dew_p(c2.p.val_SI, c2.fluid_data)


def test_td_dew_equals_0(simple_test_network):
    nw = simple_test_network
    c1, c2 = nw.get_conn(["c1", "c2"])
    heatexchanger = nw.get_comp("heatexchanger")

    c1.set_attr(m=1, p=10, fluid={"R290": 1})
    delta_T = 0.0
    # setting zeta does not work here, not sure why
    c2.set_attr(td_dew=delta_T, v=Ref(c1, 1, 0))

    # settings to prevent preprocessing of temperatures
    heatexchanger.set_attr(Q=0)
    nw.solve("design")
    nw.assert_convergence()

    # minus delta T, because td_dew is higher than T_dew
    assert approx(c2.T.val_SI - delta_T) == T_dew_p(c2.p.val_SI, c2.fluid_data)


def test_td_dew_in_preprocessing1(simple_test_network):
    nw = simple_test_network
    c1, c2 = nw.get_conn(["c1", "c2"])
    heatexchanger = nw.get_comp("heatexchanger")

    c1.set_attr(m=1, p=10, fluid={"R290": 1})
    delta_T = 10
    c2.set_attr(td_dew=delta_T)

    # settings to prevent preprocessing of temperatures
    heatexchanger.set_attr(Q=1e5, dp=0)
    nw.solve("design")
    nw.assert_convergence()

    # minus delta T, because td_dew is higher than T_dew
    assert approx(c2.T.val_SI - delta_T) == T_dew_p(c2.p.val_SI, c2.fluid_data)


def test_td_dew_in_preprocessing2(simple_test_network):
    nw = simple_test_network
    c1, c2 = nw.get_conn(["c1", "c2"])
    heatexchanger = nw.get_comp("heatexchanger")

    c1.set_attr(m=1, fluid={"R290": 1})
    delta_T = 10
    c2.set_attr(T=50, td_dew=delta_T)

    # settings to prevent preprocessing of temperatures
    heatexchanger.set_attr(Q=1e5, dp=0)
    nw.solve("design")
    nw.assert_convergence()

    # minus delta T, because td_dew is higher than T_dew
    assert approx(c2.T.val_SI - delta_T) == T_dew_p(c2.p.val_SI, c2.fluid_data)


def test_td_dew_equals_0_in_preprocessing1(simple_test_network):
    nw = simple_test_network
    c1, c2 = nw.get_conn(["c1", "c2"])
    heatexchanger = nw.get_comp("heatexchanger")

    c1.set_attr(m=1, p=10, fluid={"R290": 1})
    delta_T = 0
    c2.set_attr(td_dew=delta_T)

    # settings to prevent preprocessing of temperatures
    heatexchanger.set_attr(Q=1e5, dp=0)
    nw.solve("design")
    nw.assert_convergence()

    # minus delta T, because td_dew is higher than T_dew
    assert approx(c2.T.val_SI - delta_T) == T_dew_p(c2.p.val_SI, c2.fluid_data)


def test_td_dew_equals_0_in_preprocessing2(simple_test_network):
    nw = simple_test_network
    c1, c2 = nw.get_conn(["c1", "c2"])
    heatexchanger = nw.get_comp("heatexchanger")

    c1.set_attr(m=1, fluid={"R290": 1})
    delta_T = 0
    c2.set_attr(T=50, td_dew=delta_T)

    # settings to prevent preprocessing of temperatures
    heatexchanger.set_attr(Q=1e5, dp=0)
    nw.solve("design")
    nw.assert_convergence()

    # minus delta T, because td_dew is higher than T_dew
    assert approx(c2.T.val_SI - delta_T) == T_dew_p(c2.p.val_SI, c2.fluid_data)
