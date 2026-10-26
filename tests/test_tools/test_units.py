# -*- coding: utf-8

"""Module for testing helper functions.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tests/test_tools/test_units.py

SPDX-License-Identifier: MIT
"""
import json
import warnings

from pint import UnitRegistry
from pytest import approx
from pytest import fixture
from pytest import raises
from pytest import warns

from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Subsystem
from tespy.connections import Connection
from tespy.connections import Ref
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


def test_pressure_difference_unit_in_bar(pipe_network):
    """With pressure and pressure_difference in bar, dp=1 means 1 bar."""
    nw = pipe_network
    c1, c2 = nw.get_conn(["c1", "c2"])
    pipe = nw.get_comp("pipe")

    nw.units.set_defaults(
        pressure="bar", pressure_difference="bar", temperature="degC"
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


def _solve_pipe_network_with_ref_delta(nw, delta):
    c1, c2 = nw.get_conn(["c1", "c2"])
    pipe = nw.get_comp("pipe")
    nw.units.set_defaults(
        temperature="degC", pressure="bar", pressure_difference="bar"
    )
    c1.set_attr(fluid={"water": 1}, m=1, T=20, p=5)
    pipe.set_attr(Q=0)
    c2.set_attr(p=Ref(c1, 1, delta))
    nw.solve("design")


def test_ref_delta_preserved_on_default_change(pipe_network):
    nw = pipe_network
    _solve_pipe_network_with_ref_delta(nw, -2)
    c2 = nw.get_conn("c2")
    assert c2.p.reference.ref.delta_SI == approx(-2e5)
    assert c2.p.val_SI == approx(3e5)

    with warns(UserWarning):
        nw.units.set_defaults(pressure_difference="psi")
    nw.solve("design")

    assert c2.p.reference.ref.delta_SI == approx(-2e5)
    assert c2.p.reference.ref.delta == approx(-29.00755, rel=1e-5)
    assert c2.p.val_SI == approx(3e5)


def test_respecified_ref_delta_uses_new_defaults(pipe_network):
    nw = pipe_network
    _solve_pipe_network_with_ref_delta(nw, -2)
    c1, c2 = nw.get_conn(["c1", "c2"])

    with warns(UserWarning):
        nw.units.set_defaults(pressure_difference="psi")
    c2.set_attr(p=Ref(c1, 1, -1))
    nw.solve("design")

    assert c2.p.reference.ref.delta_SI == approx(-6894.757, rel=1e-5)
    assert c2.p.val_SI == approx(5e5 - 6894.757, rel=1e-5)


def test_ref_quantity_delta_pinned_on_default_change(pipe_network):
    nw = pipe_network
    delta = nw.units.ureg.Quantity(-2, "bar")
    _solve_pipe_network_with_ref_delta(nw, delta)
    c2 = nw.get_conn("c2")
    assert c2.p.reference.ref.delta_SI == approx(-2e5)
    assert c2.p.val_SI == approx(3e5)

    with warns(UserWarning):
        nw.units.set_defaults(pressure_difference="psi")
    nw.solve("design")

    assert c2.p.reference.ref.delta_SI == approx(-2e5)
    assert str(c2.p.reference.ref.delta.units) == "bar"
    assert c2.p.val_SI == approx(3e5)


def test_ref_quantity_delta_export_import(pipe_network):
    nw = pipe_network
    delta = nw.units.ureg.Quantity(-2, "bar")
    _solve_pipe_network_with_ref_delta(nw, delta)

    data = nw.export()
    nw2 = Network.from_dict(data)

    ref = nw2.get_conn("c2").p.reference.ref
    assert ref.delta.magnitude == approx(-2)
    assert str(ref.delta.units) == "bar"

    nw2.solve("design")
    nw2.assert_convergence()
    assert nw2.get_conn("c2").p.val_SI == approx(3e5)


def test_ref_quantity_delta_save_restore_specifications(pipe_network):
    nw = pipe_network
    delta = nw.units.ureg.Quantity(-2, "bar")
    _solve_pipe_network_with_ref_delta(nw, delta)

    data = nw.save_specifications()

    nw.clear_specifications()
    nw.restore_specifications(data)

    ref = nw.get_conn("c2").p.reference.ref
    assert ref.delta.magnitude == approx(-2)
    assert str(ref.delta.units) == "bar"

    nw.solve("design")
    nw.assert_convergence()
    assert nw.get_conn("c2").p.val_SI == approx(3e5)


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


def _specify_pipe_network(nw):
    c1 = nw.get_conn("c1")
    pipe = nw.get_comp("pipe")
    c1.set_attr(fluid={"water": 1}, m=1, T=293.15, p=5e5)
    pipe.set_attr(Q=0, dp=1e5)


def test_export_plain_values_without_unit(pipe_network):
    """Only specified values carry a unit in the export, unspecified values
    are interpreted in the default units of the network they are loaded to"""
    nw = pipe_network
    nw.units.set_defaults(
        temperature="degC", pressure="bar", pressure_difference="bar"
    )
    c1, c2 = nw.get_conn(["c1", "c2"])
    pipe = nw.get_comp("pipe")
    c1.set_attr(fluid={"water": 1}, m=1, T=20, p=5)
    pipe.set_attr(Q=0, dp=1)

    data = nw.export()
    c1_data = data["Connection"]["Connection"]["c1"]
    c2_data = data["Connection"]["Connection"]["c2"]
    pipe_data = data["Component"]["SimpleHeatExchanger"]["pipe"]

    assert c1_data["T"]["unit"] == "degree_Celsius"
    assert c1_data["p"]["unit"] == "bar"
    assert c2_data["T"]["unit"] is None
    assert c2_data["p"]["unit"] is None
    assert pipe_data["pr"]["unit"] is None

    # calculated values carry the unit they are reported in
    nw.solve("design")
    data = nw.export()
    c2_data = data["Connection"]["Connection"]["c2"]
    assert c2_data["T"]["unit"] == "degree_Celsius"
    assert c2_data["T"]["val"] == approx(c2.T.val)
    assert data["Component"]["SimpleHeatExchanger"]["pipe"]["pr"]["unit"] == "dimensionless"


def test_subsystem_from_network_uses_target_network_units(pipe_network):
    """A subsystem created from a network reports in the default units of the
    network it is added to, not in the units of the source network"""
    source = pipe_network
    source.units.set_defaults(
        temperature="degF", pressure="MPa", pressure_difference="kPa"
    )
    c1 = source.get_conn("c1")
    c1.set_attr(fluid={"water": 1}, m=1, T=68, p=0.5)
    source.get_comp("pipe").set_attr(Q=0, dp=100)
    sub = Subsystem.from_network("sub", source)

    nw = Network()
    nw.units.set_defaults(
        temperature="degC", pressure="bar", pressure_difference="bar"
    )
    c_in = Connection(Source("source"), "out1", sub.inlet, "in1", label="in")
    c_out = Connection(sub.outlet, "out1", Sink("sink"), "in1", label="out")
    nw.add_conns(c_in, c_out)
    nw.add_subsystems(sub)

    nw.solve("design")
    nw.assert_convergence()

    c2 = sub.get_conn("c2")
    assert c2.p.unit == "bar"
    assert c2.p.val == approx(4)
    assert c2.T.unit == "degree_Celsius"
    assert c2.T.val == approx(20, abs=0.1)
    assert nw.results["Connection"].loc["sub_c2", "p_unit"] == "bar"
    assert nw.results["Connection"].loc["sub_c2", "T_unit"] == "degree_Celsius"


def test_subsystem_from_network_keeps_starting_value_units(pipe_network):
    """Starting values taken from a solved network keep the units of that
    network, independent of the defaults of the target network"""
    nw = pipe_network
    nw.units.set_defaults(
        temperature="degC", pressure="bar", pressure_difference="bar",
        enthalpy="kJ/kg", mass_flow="t/h"
    )
    c1 = nw.get_conn("c1")
    pipe = nw.get_comp("pipe")
    c1.set_attr(fluid={"water": 1}, m=3.6, T=20, p=5)
    pipe.set_attr(Q=0, dp=1)
    nw.solve("design")
    source_state = {
        key: nw.get_conn("c2").get_attr(key).val_SI for key in ("m", "p", "h")
    }
    sub = Subsystem.from_network("sub", nw, keep="starting_values")

    nw2 = Network()
    c_in = Connection(Source("source"), "out1", sub.inlet, "in1", label="in")
    c_out = Connection(sub.outlet, "out1", Sink("sink"), "in1", label="out")
    nw2.add_conns(c_in, c_out)
    nw2.add_subsystems(sub)
    c_in.set_attr(fluid={"water": 1}, m=1, T=293.15, p=5e5)
    sub.get_comp("pipe").set_attr(Q=0, dp=1e5)

    # after preprocessing the starting values hold the state of the source
    # network in SI, the starting value units are the ones of the source
    nw2.solve("design", init_only=True)
    c2 = sub.get_conn("c2")
    for key, value_SI in source_state.items():
        assert c2.get_attr(key).val_SI == approx(value_SI), key
    # the enthalpy is the only solver variable of c2, mass flow and pressure
    # are determined by the presolving
    assert str(c2.h.val0_with_unit.units) == "kilojoule / kilogram"
    assert c2.h.val0_with_unit.m_as("J/kg") == approx(source_state["h"])

    nw2.solve("design")
    nw2.assert_convergence()
    assert str(c2.p.val0_with_unit.units) == "bar"
    assert c2.p.val0_with_unit.m_as("Pa") == approx(4e5)
    assert str(c2.m.val0_with_unit.units) == "metric_ton / hour"
    assert c2.p.unit == "pascal"
    assert c2.h.unit == "joule / kilogram"


def test_import_ignores_file_units_of_unspecified_values(pipe_network):
    """Files of version 0.11.2 and earlier carry SI units on unspecified
    values, they must not pin the display unit of the imported values"""
    _specify_pipe_network(pipe_network)
    data = pipe_network.export()
    data["Connection"]["Connection"]["c2"]["T"]["unit"] = "K"
    data["Connection"]["Connection"]["c2"]["p"]["unit"] = "Pa"
    data["Component"]["SimpleHeatExchanger"]["pipe"]["Q"]["unit"] = "W"

    nw = Network.from_dict(data)
    nw.units.set_defaults(
        temperature="degC", pressure="bar", pressure_difference="bar",
        heat="kW"
    )
    nw.get_comp("pipe").set_attr(Q=None)
    nw.get_conn("c2").set_attr(T=20)
    nw.solve("design")
    nw.assert_convergence()

    c2 = nw.get_conn("c2")
    assert c2.p.unit == "bar"
    assert c2.p.val == approx(4)
    assert c2.T.unit == "degree_Celsius"
    assert nw.get_comp("pipe").Q.unit == "kilowatt"


def test_quantity_from_foreign_registry(pipe_network):
    """Quantities from any pint registry are accepted: the conversion cache
    must not compare units of different registries"""
    nw = pipe_network
    nw.units.set_defaults(
        temperature="degC", pressure="bar", pressure_difference="bar"
    )
    c1, c2 = nw.get_conn(["c1", "c2"])
    pipe = nw.get_comp("pipe")
    ureg = UnitRegistry()
    c1.set_attr(fluid={"water": 1}, m=1, T=ureg.Quantity(20, "degC"), p=5)
    c2.set_attr(p=ureg.Quantity(4, "bar"))
    pipe.set_attr(Q=0)

    nw.solve("design")
    nw.assert_convergence()

    assert c1.T.val_SI == approx(293.15)
    assert c2.p.val_SI == approx(4e5)
    assert pipe.dp.val == approx(1)


def test_incompatible_units_rejected(pipe_network):
    """Units are checked against the quantity of the container when they are
    set, no network is required"""
    c1 = pipe_network.get_conn("c1")
    with raises(ValueError, match="not compatible with temperature"):
        c1.set_attr(T=UnitRegistry().Quantity(5, "bar"))
    with raises(ValueError, match="not compatible with pressure"):
        c1.p.set_attr(val=5, unit="degC", is_set=True)
    with raises(ValueError, match="not compatible with temperature_difference"):
        c1.set_attr(td_dew=UnitRegistry().Quantity(5, "degC"))
    c1.set_attr(td_dew=UnitRegistry().Quantity(5, "kelvin"))
    assert c1.td_dew.unit == "kelvin"


def test_incompatible_file_unit_rejected_on_import(pipe_network):
    _specify_pipe_network(pipe_network)
    data = pipe_network.export()
    data["Connection"]["Connection"]["c1"]["T"]["unit"] = "bar"
    with raises(ValueError, match="not compatible with temperature"):
        Network.from_dict(data)
