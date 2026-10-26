# -*- coding: utf-8

"""Tests of the saturation related specifications of connections.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_connections/test_connection_saturation.py

SPDX-License-Identifier: MIT
"""
import pytest
from CoolProp.CoolProp import get_global_param_string
from pytest import approx
from pytest import fixture
from pytest import mark

from tespy.components import Merge
from tespy.components import Pipe
from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.connections import Ref
from tespy.networks import Network
from tespy.tools.fluid_properties.functions import T_bubble_p
from tespy.tools.fluid_properties.functions import T_dew_p
from tespy.tools.fluid_properties.functions import p_bubble_T
from tespy.tools.fluid_properties.functions import p_dew_T
from tespy.tools.helpers import TESPyNetworkError


@fixture
def simple_test_network():
    nw = Network()
    nw.units.set_defaults(
        temperature="degC",
        pressure="bar",
        pressure_difference="bar"
    )

    so = Source("source")
    si = Sink("sink")

    heatexchanger = SimpleHeatExchanger("heatexchanger")

    c1 = Connection(so, "out1", heatexchanger, "in1", label="c1")
    c2 = Connection(heatexchanger, "out1", si, "in1", label="c2")

    nw.add_conns(c1, c2)
    return nw


def test_td_dew_convergence_helper(simple_test_network):
    """The old convergence helper for td_dew led to temporarily oscillating
    residuals"""
    nw = simple_test_network

    c1, c2 = nw.get_conn(["c1", "c2"])
    heatexchanger = nw.get_comp("heatexchanger")

    c1.set_attr(m=1, p=10, fluid={"water": 1})
    c2.set_attr(td_bubble=0)

    # settings to prevent preprocessing of temperatures
    heatexchanger.set_attr(Q=1e5, zeta_d4=0)

    nw.solve("design")
    nw.assert_convergence()
    assert nw.problem.iter < 10


def test_td_bubble_convergence_helper(simple_test_network):
    """The old convergence helper for td_dew led to temporarily oscillating
    residuals"""
    nw = simple_test_network

    c1, c2 = nw.get_conn(["c1", "c2"])
    heatexchanger = nw.get_comp("heatexchanger")

    c1.set_attr(m=1, p=10, fluid={"water": 1})
    c2.set_attr(td_bubble=0)

    # settings to prevent preprocessing of temperatures
    heatexchanger.set_attr(Q=1e5, zeta_d4=0)

    nw.solve("design")
    nw.assert_convergence()
    assert nw.problem.iter < 10


@mark.skipif(
    get_global_param_string("REFPROP_version") == "n/a",
    reason="This test requires REFPROP, dependency is missing."
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
    heatexchanger.set_attr(Q=1e5, zeta_d4=0)

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
    heatexchanger.set_attr(Q=1e5, zeta_d4=0)
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
    heatexchanger.set_attr(Q=1e5, zeta_d4=0)
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


def test_td_bubble_in_postprocessing(simple_test_network):
    nw = simple_test_network
    c1, c2 = nw.get_conn(["c1", "c2"])
    heatexchanger = nw.get_comp("heatexchanger")

    c1.set_attr(m=1, fluid={"R290": 1})
    c2.set_attr(T=50, p=20)

    # settings to prevent preprocessing of temperatures
    heatexchanger.set_attr(Q=1e5, dp=0)
    nw.solve("design")
    nw.assert_convergence()

    # minus delta T, because td_dew is higher than T_dew
    assert approx(c2.td_bubble.val) == 7.26229


def test_T_bubble_specified(simple_test_network):
    nw = simple_test_network
    c1, c2 = nw.get_conn(["c1", "c2"])
    heatexchanger = nw.get_comp("heatexchanger")

    c1.set_attr(m=1, fluid={"R290": 1})
    c2.set_attr(T=50, T_bubble=60)

    # settings to prevent preprocessing of temperatures
    heatexchanger.set_attr(Q=1e5, dp=0)
    nw.solve("design")
    nw.assert_convergence()

    # minus delta T, because td_dew is higher than T_dew
    assert approx(c2.td_bubble.val) == 10


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


def test_td_dew_in_postprocessing(simple_test_network):
    nw = simple_test_network
    c1, c2 = nw.get_conn(["c1", "c2"])
    heatexchanger = nw.get_comp("heatexchanger")

    c1.set_attr(m=1, fluid={"R290": 1})
    c2.set_attr(T=50, p=10)

    # settings to prevent preprocessing of temperatures
    heatexchanger.set_attr(Q=1e5, dp=0)
    nw.solve("design")
    nw.assert_convergence()

    # minus delta T, because td_dew is higher than T_dew
    assert approx(c2.td_dew.val) == 23.05767


def test_T_dew_specified(simple_test_network):
    nw = simple_test_network
    c1, c2 = nw.get_conn(["c1", "c2"])
    heatexchanger = nw.get_comp("heatexchanger")

    c1.set_attr(m=1, fluid={"R290": 1})
    c2.set_attr(T=70, T_dew=60)

    # settings to prevent preprocessing of temperatures
    heatexchanger.set_attr(Q=1e5, dp=0)
    nw.solve("design")
    nw.assert_convergence()

    # minus delta T, because td_dew is higher than T_dew
    assert approx(c2.td_dew.val) == 10


def _pipe_network():
    nw = Network(iterinfo=False)
    nw.units.set_defaults(pressure="bar", temperature="degC")
    pipe = Pipe("pipe")
    c1 = Connection(Source("source"), "out1", pipe, "in1", label="c1")
    c2 = Connection(pipe, "out1", Sink("sink"), "in1", label="c2")
    nw.add_conns(c1, c2)
    pipe.set_attr(pr=1, Q=0)
    return nw, c1, c2


def test_dew_temperature_determines_free_pressure():
    nw, c1, c2 = _pipe_network()
    c1.set_attr(fluid={"water": 1}, m=1, T=20, T_dew=150)
    nw.solve("design")
    nw.assert_convergence()
    assert c1.p.val_SI == pytest.approx(p_dew_T(150 + 273.15, c1.fluid_data))


def test_bubble_temperature_determines_free_pressure():
    nw, c1, c2 = _pipe_network()
    c1.set_attr(fluid={"water": 1}, m=1, T=20, T_bubble=120)
    nw.solve("design")
    nw.assert_convergence()
    assert c1.p.val_SI == pytest.approx(p_bubble_T(120 + 273.15, c1.fluid_data))


@pytest.mark.parametrize("parameter", ["T_dew", "T_bubble"])
def test_saturation_temperature_on_linked_pressure_is_overdetermined(parameter):
    """The outlet pressure follows the inlet pressure through the pipe, a
    saturation temperature on the outlet competes with it."""
    nw, c1, c2 = _pipe_network()
    c1.set_attr(fluid={"water": 1}, m=1, p=1, T=20)
    c2.set_attr(**{parameter: 150})
    with pytest.raises(TESPyNetworkError, match=f"{parameter} specification overdetermines"):
        nw.solve("design")


@pytest.mark.parametrize(
    "parameter", ["T_dew", "T_bubble", "x", "td_dew", "td_bubble"]
)
def test_saturation_specification_with_variable_composition_is_rejected(parameter):
    """Saturation properties are not available for a variable composition,
    the specification must not be dropped silently."""
    nw = Network(iterinfo=False)
    nw.units.set_defaults(pressure="bar", temperature="degC")
    merge = Merge("merge")
    c1 = Connection(Source("source 1"), "out1", merge, "in1", label="c1")
    c2 = Connection(Source("source 2"), "out1", merge, "in2", label="c2")
    c3 = Connection(merge, "out1", Sink("sink"), "in1", label="c3")
    nw.add_conns(c1, c2, c3)
    c1.set_attr(fluid={"water": 1}, m=1, T=20)
    c2.set_attr(fluid={"water": 0.5, "N2": 0.5}, m=1, T=20)
    c3.set_attr(**{parameter: 0.5})
    with pytest.raises(TESPyNetworkError, match=f"{parameter}.*variable fluid composition"):
        nw.solve("design")


@pytest.mark.parametrize("parameter", ["T_dew", "T_bubble"])
def test_saturation_temperature_with_pressure_is_rejected(parameter):
    """A saturation temperature together with the pressure on the same
    connection is a direct overdetermination."""
    nw, c1, c2 = _pipe_network()
    c1.set_attr(fluid={"water": 1}, m=1, p=5, **{parameter: 150})
    with pytest.raises(TESPyNetworkError, match="simultaneously specify pressure"):
        nw.solve("design")
