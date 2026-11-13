# -*- coding: utf-8

"""Tests of the exergy results of connections.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_connections/test_exergy_results.py

SPDX-License-Identifier: MIT
"""
import numpy as np
import pint
import pytest

from tespy.components import Compressor
from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.connections import HAConnection
from tespy.networks import Network
from tespy.tools.data_containers import OptionalResult
from tespy.tools.fluid_properties.functions import calc_physical_exergy

EXERGY_RESULTS = ["e_T", "e_M", "e_PH", "E_T", "E_M", "E_PH"]


@pytest.fixture
def network():
    nw = Network(iterinfo=False)
    nw.units.set_defaults(
        pressure="bar", temperature="degC", specific_energy="kJ/kg", power="kW"
    )
    comp = Compressor("compressor")
    c1 = Connection(Source("source"), "out1", comp, "in1", label="1")
    c2 = Connection(comp, "out1", Sink("sink"), "in1", label="2")
    nw.add_conns(c1, c2)
    c1.set_attr(fluid={"air": 1}, p=1, T=20, m=1)
    c2.set_attr(p=5)
    comp.set_attr(eta_s=0.8)
    nw.solve("design")
    assert nw.status == 0
    return nw


def test_exergy_results_are_optional_results(network):
    c = network.get_conn("2")
    assert list(c.optional_results()) == EXERGY_RESULTS
    for key in EXERGY_RESULTS:
        data = c.get_attr(key)
        assert isinstance(data, OptionalResult)
        assert data.result_only
        assert not data.is_calculated
    with pytest.raises(TypeError):
        c.set_attr(e_PH=1)


def test_physical_exergy_values_and_units(network):
    c = network.get_conn("2")
    c._get_physical_exergy(1e5, 293.15)
    e_T, e_M = calc_physical_exergy(
        c.h.val_SI, c.s.val_SI, c.p.val_SI, 1e5, 293.15, c.fluid_data,
        c.mixing_rule, c.T.val_SI
    )
    assert c.e_T.val_SI == pytest.approx(e_T)
    assert c.e_M.val_SI == pytest.approx(e_M)
    assert c.e_PH.val_SI == pytest.approx(e_T + e_M)
    assert c.E_PH.val_SI == pytest.approx((e_T + e_M) * c.m.val_SI)
    assert c.e_PH.unit == network.units.canonical("kJ/kg")
    assert c.e_PH.val == pytest.approx(c.e_PH.val_SI / 1e3)
    assert c.E_PH.unit == network.units.canonical("kW")
    assert c.E_PH.val == pytest.approx(c.E_PH.val_SI / 1e3)


def test_optional_results_in_results_dataframe(network):
    exergy_cols = [col for key in EXERGY_RESULTS for col in (key, f"{key}_unit")]

    results = network.results["Connection"]
    base_cols = list(results.columns)
    assert not set(exergy_cols) & set(base_cols)

    c = network.get_conn("2")
    c._get_physical_exergy(1e5, 293.15)
    results = network.results["Connection"]
    assert list(results.columns) == base_cols + exergy_cols
    for key in EXERGY_RESULTS:
        data = c.get_attr(key)
        assert results.loc["2", key] == pytest.approx(data.val)
        assert results.loc["2", f"{key}_unit"] == data.unit
        assert np.isnan(results.loc["1", key])
    assert results.loc["2", "p"] == pytest.approx(c.p.val)

    network.get_conn("1")._get_physical_exergy(1e5, 293.15)
    results = network.results["Connection"]
    assert results["e_PH"].notna().all()

    network.solve("design")
    assert not c.e_PH.is_calculated
    assert list(network.results["Connection"].columns) == base_cols


def test_ambient_state_as_quantity(network):
    c = network.get_conn("2")
    c._get_physical_exergy(1e5, 293.15)
    reference = c.e_PH.val_SI
    ureg = pint.UnitRegistry()
    c._get_physical_exergy(1 * ureg.bar, ureg.Quantity(20, "degC"))
    assert c.e_PH.val_SI == pytest.approx(reference)


def test_deprecated_attributes_return_SI_values(network):
    c = network.get_conn("2")
    c._get_physical_exergy(1e5, 293.15)
    for old, new in (
        ("ex_therm", "e_T"), ("ex_mech", "e_M"), ("ex_physical", "e_PH"),
        ("Ex_therm", "E_T"), ("Ex_mech", "E_M"), ("Ex_physical", "E_PH")
    ):
        with pytest.warns(FutureWarning, match=f"'{new}.val_SI'"):
            assert getattr(c, old) == c.get_attr(new).val_SI


def test_humid_air_connection_exergy_results():
    nw = Network(iterinfo=False)
    nw.units.set_defaults(pressure="bar", temperature="degC")
    he = SimpleHeatExchanger("heater")
    c1 = HAConnection(Source("source"), "out1", he, "in1", label="1")
    c2 = HAConnection(he, "out1", Sink("sink"), "in1", label="2")
    nw.add_conns(c1, c2)
    c1.set_attr(p=1, T=15, w=0.004, m=1)
    c2.set_attr(T=35)
    he.set_attr(dp=0)
    nw.solve("design")
    assert nw.status == 0
    assert list(c2.optional_results()) == EXERGY_RESULTS
    c2._get_physical_exergy(1e5, 293.15)
    assert c2.e_PH.val_SI > 0
    assert c2.E_PH.val_SI == pytest.approx(c2.e_PH.val_SI * c2.m.val_SI)
    assert nw.results["HAConnection"].loc["2", "e_PH"] == pytest.approx(c2.e_PH.val)
    nw.solve("design")
    assert not c2.e_PH.is_calculated
