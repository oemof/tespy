# -*- coding: utf-8

"""Module for testing helper functions.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tests/test_tools/test_units.py

SPDX-License-Identifier: MIT
"""
import warnings

from pint import UnitRegistry
from pytest import approx
from pytest import fixture
from pytest import raises
from pytest import warns

from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.networks import Network
from tespy.tools.units import Units


def test_units_set_defaults_incompatible():
    units = Units()
    with raises(ValueError):
        units.set_defaults(temperature="meter")


def test_units_set_defaults_temperature_difference_incompatible():
    units = Units()
    with raises(ValueError):
        units.set_defaults(temperature_difference="degC")


def test_units_set_ureg():
    units = Units()
    ureg = units.ureg
    ureg_custom = UnitRegistry()
    units.set_ureg(ureg_custom)
    assert units.ureg is not ureg
    assert units.ureg is ureg_custom


@fixture
def pipe_network():
    """Return (nw, c_in, c_out, pipe) for a Source->SimpleHeatExchanger->Sink."""
    nw = Network()
    source = Source("source")
    sink = Sink("sink")
    pipe = SimpleHeatExchanger("pipe")

    c1 = Connection(source, "out1", pipe, "in1", label="c1")
    c2 = Connection(pipe, "out1", sink, "in1", label="c2")
    nw.add_conns(c1, c2)
    return nw


def test_pressure_unit_propagates_to_pressure_difference(pipe_network):
    """Setting only 'pressure' must also apply to dp (pressure_difference).

    When the user sets pressure='bar' without explicitly setting
    pressure_difference, the backwards-compatibility path should carry 'bar'
    over to pressure_difference as well, so dp=1 means 1 bar = 100 000 Pa.
    """
    nw = pipe_network
    c1, c2 = nw.get_conn(["c1", "c2"])
    pipe = nw.get_comp("pipe")

    nw.units.set_defaults(
        pressure="bar", temperature="degC"
    )

    c1.set_attr(fluid={"water": 1}, m=1, T=20, p=5)
    pipe.set_attr(Q=0, dp=1)

    nw.solve("design")

    assert approx(c1.p.val - c2.p.val) == pipe.dp.val
    assert approx(c1.p.val_SI - c2.p.val_SI) == 1e5


def test_pressure_difference_unit_independent_from_pressure_unit(pipe_network):
    """pressure_difference unit must be independent of the pressure unit.

    When pressure='bar' and pressure_difference='Pa', dp=100000 (Pa) must
    produce the same 1 bar pressure drop as in the previous test.  The .val
    of dp (100 000) will differ from the .val difference of inlet and outlet
    pressures (1.0 bar).
    """
    nw = pipe_network
    c1, c2 = nw.get_conn(["c1", "c2"])
    pipe = nw.get_comp("pipe")

    nw.units.set_defaults(
        pressure="bar", pressure_difference="Pa", temperature="degC"
    )

    c1.set_attr(fluid={"water": 1}, m=1, T=20, p=5)
    pipe.set_attr(Q=0, dp=1e5)

    nw.solve("design")

    # dp.val is in Pa, calculated difference from pressure values in bar
    p_diff_bar = c1.p.val - c2.p.val
    assert p_diff_bar != approx(pipe.dp.val)
    # but the underlying SI values must agree: 1 bar = 100 000 Pa
    assert approx(c1.p.val_SI - c2.p.val_SI) == pipe.dp.val_SI


def _solve_pipe_network_in_degC(nw):
    c1 = nw.get_conn("c1")
    pipe = nw.get_comp("pipe")
    nw.units.set_defaults(
        temperature="degC", pressure="bar", pressure_difference="bar"
    )
    c1.set_attr(fluid={"water": 1}, m=1, T=20, p=5)
    pipe.set_attr(Q=0, dp=1)
    nw.solve("design")


def test_set_defaults_no_warning_before_solve(pipe_network):
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        pipe_network.units.set_defaults(temperature="degF")


def test_set_defaults_after_solve_warns(pipe_network):
    nw = pipe_network
    _solve_pipe_network_in_degC(nw)
    with warns(UserWarning, match="changing default units"):
        nw.units.set_defaults(temperature="degF")


def test_default_change_preserves_physical_values(pipe_network):
    nw = pipe_network
    _solve_pipe_network_in_degC(nw)
    c1, c2 = nw.get_conn(["c1", "c2"])
    T_spec_SI = c1.T.val_SI
    T_result_SI = c2.T.val_SI

    with warns(UserWarning):
        nw.units.set_defaults(temperature="degF")
    nw.solve("design")

    assert c1.T.unit == "degree_Fahrenheit"
    assert c1.T.val == approx(68.0)
    assert c1.T.val_SI == approx(T_spec_SI)
    assert c2.T.unit == "degree_Fahrenheit"
    assert c2.T.val_SI == approx(T_result_SI)


def test_respecified_uses_new_defaults(pipe_network):
    nw = pipe_network
    _solve_pipe_network_in_degC(nw)
    c1 = nw.get_conn("c1")

    with warns(UserWarning):
        nw.units.set_defaults(temperature="degF")
    c1.set_attr(T=68)
    nw.solve("design")

    assert c1.T.unit == "degree_Fahrenheit"
    assert c1.T.val_SI == approx(293.15)


def test_specified_quantity_keeps_unit_on_default_change(pipe_network):
    nw = pipe_network
    _solve_pipe_network_in_degC(nw)
    c1, c2 = nw.get_conn(["c1", "c2"])
    c1.set_attr(T=nw.units.ureg.Quantity(20, "degC"))
    nw.solve("design")

    with warns(UserWarning):
        nw.units.set_defaults(temperature="degF")
    nw.solve("design")

    assert c1.T.unit == "degree_Celsius"
    assert c1.T.val == approx(20.0)
    assert c1.T.val_SI == approx(293.15)
    assert c2.T.unit == "degree_Fahrenheit"
