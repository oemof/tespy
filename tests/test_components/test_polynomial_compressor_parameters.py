# -*- coding: utf-8

"""Tests of the swept volume and polynomial parameters of the
PolynomialCompressor.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_polynomial_compressor_parameters.py

SPDX-License-Identifier: MIT
"""
import numpy as np
import pandas as pd
import pytest

from tespy.components import PolynomialCompressor
from tespy.components import Sink
from tespy.components import Source
from tespy.components.displacementmachinery.polynomial_compressor import calc_EN12900
from tespy.components.displacementmachinery.polynomial_compressor import generate_eta_polys_from_data
from tespy.components.displacementmachinery.polynomial_compressor import swept_volume_from_displacement
from tespy.connections import Connection
from tespy.networks import Network
from tespy.tools.characteristics import EN12900Polynomial
from tespy.tools.data_containers import CharacteristicData
from tespy.tools.data_containers import ComponentParameter


@pytest.fixture
def reference_state():
    return {
        "T_sh": 20,
        "T_sc": 0,
        "frequency_poly": 50.0,
        "swept_volume": swept_volume_from_displacement(214, 20.0),
    }


@pytest.fixture
def eta_polys(reference_state):
    power = pd.DataFrame(
        columns=[10, 7.5, 5, 0, -5, -10], index=[30, 40, 50], dtype=float
    )
    cooling = power.copy()
    cooling.loc[30] = [465600, 424100, 385500, 316700, 257900, 208000]
    cooling.loc[40] = [418900, 380400, 344800, 281400, 227400, 181600]
    cooling.loc[50] = [365900, 331300, 299200, 242100, 193700, 152900]
    power.loc[30] = [62.0, 61.8, 61.8, 61.8, 61.7, 61.3]
    power.loc[40] = [78.0, 78.0, 78.0, 78.0, 77.7, 76.8]
    power.loc[50] = [99.2, 99.2, 99.2, 98.9, 98.1, 96.5]
    return generate_eta_polys_from_data(
        power * 1000, cooling, "R134a", T_sh=reference_state["T_sh"],
        T_sc=reference_state["T_sc"], frequency=reference_state["frequency_poly"],
        swept_volume=reference_state["swept_volume"]
    )


def _network(instance, **defaults):
    nw = Network(iterinfo=False)
    units = {
        "pressure": "bar", "pressure_difference": "bar", "temperature": "degC"
    }
    units.update(defaults)
    nw.units.set_defaults(**units)
    source = Source("source")
    sink = Sink("sink")
    c1 = Connection(source, "out1", instance, "in1")
    c2 = Connection(instance, "out1", sink, "in1")
    nw.add_conns(c1, c2)
    return nw, c1, c2


def test_parameter_types():
    compressor = PolynomialCompressor("compressor")
    assert isinstance(compressor.swept_volume, ComponentParameter)
    assert compressor.swept_volume.quantity == "volume"
    assert isinstance(compressor.eta_s_poly, CharacteristicData)
    assert compressor.eta_s_poly.num_inputs == 2
    assert compressor.eta_s_poly.x_quantity == "temperature"
    assert compressor.eta_s_poly.y_quantity == "temperature"


def test_reference_state_sets_swept_volume():
    compressor = PolynomialCompressor("compressor")
    with pytest.warns(FutureWarning):
        compressor.set_attr(reference_state={
            "displacement": 214, "frequency_displacement": 20.0
        })
    assert compressor.swept_volume.is_set
    assert compressor.swept_volume.val == pytest.approx(
        swept_volume_from_displacement(214, 20.0)
    )


def test_reference_state_dict_in_helpers_is_deprecated(reference_state):
    power = pd.DataFrame(
        columns=[10, 0, -10], index=[30, 40, 50], dtype=float
    )
    cooling = power.copy()
    power.loc[:] = [[62e3, 61.8e3, 61.3e3]] * 3
    cooling.loc[:] = [[465600, 316700, 208000]] * 3
    with pytest.warns(FutureWarning, match="dictionary is deprecated"):
        legacy = generate_eta_polys_from_data(power, cooling, "R134a", reference_state)
    new = generate_eta_polys_from_data(
        power, cooling, "R134a", T_sh=20, T_sc=0, frequency=50.0,
        swept_volume=reference_state["swept_volume"]
    )
    np.testing.assert_allclose(legacy[0], new[0])
    with pytest.raises(ValueError, match="swept_volume"):
        generate_eta_polys_from_data(power, cooling, "R134a", T_sh=20, T_sc=0, frequency=50.0)


def test_generating_functions_follow_unit_system(reference_state):
    """A table in the units of a network gives the standard form polynomial
    of the same table in °C and W."""
    from tespy.tools.units import Units
    power = pd.DataFrame(
        columns=[10, 0, -10], index=[30, 40, 50], dtype=float
    )
    cooling = power.copy()
    power.loc[:] = [[62e3, 61.8e3, 61.3e3]] * 3
    cooling.loc[:] = [[465600, 316700, 208000]] * 3
    standard = generate_eta_polys_from_data(
        power, cooling, "R134a", T_sh=20, T_sc=0, frequency=50.0,
        swept_volume=reference_state["swept_volume"]
    )

    units = Units()
    units.set_defaults(
        temperature="K", power="kW", heat="kW", frequency="1/min",
        volume="cm**3"
    )
    power_K = power.copy() / 1000
    power_K.columns = power.columns + 273.15
    power_K.index = power.index + 273.15
    cooling_K = cooling.copy() / 1000
    cooling_K.columns = power_K.columns
    cooling_K.index = power_K.index
    in_units = generate_eta_polys_from_data(
        power_K, cooling_K, "R134a", T_sh=20, T_sc=0, frequency=3000,
        swept_volume=reference_state["swept_volume"] * 1e6, units=units
    )
    np.testing.assert_allclose(in_units[0], standard[0], rtol=1e-6)
    np.testing.assert_allclose(in_units[1], standard[1], rtol=1e-6)

    Q_ = units.ureg.Quantity
    quantities = generate_eta_polys_from_data(
        power, cooling, "R134a", T_sh=Q_(20, "delta_degC"), T_sc=0,
        frequency=Q_(50, "Hz"), swept_volume=Q_(reference_state["swept_volume"], "m**3")
    )
    np.testing.assert_allclose(quantities[0], standard[0], rtol=1e-6)


def test_polynomial_carries_degC(eta_polys, reference_state):
    """A compressor whose polynomial parameters declare no units reads the
    °C the EN 12900 polynomial carries."""
    from tespy.tools.data_containers import CharacteristicData
    eta_s_poly, eta_vol_poly = eta_polys
    compressor = PolynomialCompressor("compressor")
    nw, c1, c2 = _network(compressor, temperature="K")
    compressor.parameters["eta_s_poly"] = CharacteristicData(
        num_inputs=2, x_quantity="temperature", y_quantity="temperature"
    )
    compressor.eta_s_poly = compressor.parameters["eta_s_poly"]
    c1.set_attr(fluid={"R134a": 1}, T=293.15, td_dew=reference_state["T_sh"])
    c2.set_attr(T_dew=323.15)
    compressor.set_attr(
        eta_s_poly=EN12900Polynomial(eta_s_poly),
        eta_vol_poly=EN12900Polynomial(eta_vol_poly),
        swept_volume=reference_state["swept_volume"],
        frequency=reference_state["frequency_poly"],
        dissipation_ratio=0
    )
    nw.solve("design")
    nw.assert_convergence()
    assert compressor.eta_s.val == pytest.approx(calc_EN12900(eta_s_poly, 0, 50), rel=1e-3)


def test_network_file_of_version_0_11_loads(reference_state):
    """Polynomial coefficients stored as val and the reference_state dictionary
    of files written before 0.12 are translated on import."""
    coefficients = list(np.linspace(0.1, 1, 10) * 1e-2)
    compressor = PolynomialCompressor("compressor")
    nw, c1, c2 = _network(compressor)
    c1.set_attr(fluid={"R134a": 1})
    data = nw.export()
    legacy = data["Component"]["PolynomialCompressor"]["compressor"]
    legacy["eta_s_poly"] = {"val": coefficients, "is_set": True}
    legacy["eta_vol_poly"] = {"val": float("nan"), "is_set": False}
    legacy["reference_state"] = {
        "val": {"displacement": 214, "frequency_displacement": 20.0,
                "T_sh": 20, "T_sc": 0, "frequency_poly": 50.0},
        "is_set": True
    }
    legacy.pop("swept_volume")

    with pytest.warns(FutureWarning):
        restored = Network.from_dict(data).get_comp("compressor")
    assert isinstance(restored.eta_s_poly.char_func, EN12900Polynomial)
    np.testing.assert_allclose(restored.eta_s_poly.val, coefficients)
    assert restored.eta_s_poly.is_set
    assert restored.eta_vol_poly.char_func is None
    assert restored.swept_volume.is_set
    assert restored.swept_volume.val == pytest.approx(
        swept_volume_from_displacement(214, 20.0)
    )


def test_coefficient_list_is_wrapped():
    compressor = PolynomialCompressor("compressor")
    with pytest.warns(FutureWarning):
        compressor.set_attr(eta_s_poly=list(range(10)))
    assert isinstance(compressor.eta_s_poly.char_func, EN12900Polynomial)
    with pytest.raises(ValueError):
        compressor.set_attr(eta_s_poly=EN12900Polynomial([1, 2, 3]))


@pytest.mark.parametrize("temperature_unit", ["degC", "K"])
def test_polynomials_evaluated_in_degC(eta_polys, reference_state, temperature_unit):
    eta_s_poly, eta_vol_poly = eta_polys
    compressor = PolynomialCompressor("compressor")
    nw, c1, c2 = _network(compressor, temperature=temperature_unit)
    offset = 273.15 if temperature_unit == "K" else 0
    c1.set_attr(fluid={"R134a": 1}, T=20 + offset, td_dew=reference_state["T_sh"])
    c2.set_attr(T_dew=50 + offset)
    compressor.set_attr(
        eta_s_poly=EN12900Polynomial(eta_s_poly),
        eta_vol_poly=EN12900Polynomial(eta_vol_poly),
        swept_volume=reference_state["swept_volume"],
        frequency=reference_state["frequency_poly"],
        dissipation_ratio=0
    )
    nw.solve("design")
    nw.assert_convergence()

    assert compressor.eta_s.val == pytest.approx(calc_EN12900(eta_s_poly, 0, 50), rel=1e-3)
    assert compressor.eta_vol.val == pytest.approx(calc_EN12900(eta_vol_poly, 0, 50), rel=1e-3)
    assert c1.v.val_SI / compressor.eta_vol.val == pytest.approx(
        reference_state["swept_volume"] * reference_state["frequency_poly"], rel=1e-3
    )


def test_swept_volume_in_network_units(reference_state):
    compressor = PolynomialCompressor("compressor")
    nw, c1, c2 = _network(compressor, volume="cm**3")
    c1.set_attr(fluid={"R134a": 1}, x=1, T=5)
    c2.set_attr(p=6, T=150)
    compressor.set_attr(
        eta_vol=0.8, frequency=25.0, dissipation_ratio=0,
        swept_volume=reference_state["swept_volume"] * 1e6
    )
    nw.solve("design")
    nw.assert_convergence()
    assert compressor.swept_volume.val_SI == pytest.approx(reference_state["swept_volume"])
    assert c1.v.val_SI == pytest.approx(
        reference_state["swept_volume"] * 25.0 * 0.8, rel=1e-3
    )
