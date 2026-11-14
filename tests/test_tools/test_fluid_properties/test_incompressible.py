# -*- coding: utf-8

"""Module for testing fluid properties of incompressibles.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_tools/test_fluid_properties/test_incompressible.py

SPDX-License-Identifier: MIT
"""
import numpy as np
from pytest import approx
from pytest import fixture
from pytest import mark

from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.fluid_properties import isentropic
from tespy.tools.fluid_properties import s_mix_ph
from tespy.tools.fluid_properties import v_mix_pT
from tespy.tools.fluid_properties.wrappers import IncompressibleFluidWrapper


@fixture
def property_data():
    """Heat transfer oil data in SI units.

    Source: ISTpure Data Sheet - Standard Heat Transfer Oil, page 3 (metric
    units), https://de.scribd.com/document/553659821/ISTpure-Data-Sheet-Standard-Heat-Transfer-Oil
    """
    return {
        "temperature_data": np.array([
            0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140,
            150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270,
            280, 290, 300, 310, 315
        ]) + 273.15,
        "heat_capacity_data": np.array([
            1.894, 1.930, 1.967, 2.003, 2.040, 2.076, 2.113, 2.150, 2.186,
            2.223, 2.259, 2.296, 2.332, 2.369, 2.405, 2.442, 2.478, 2.515,
            2.552, 2.588, 2.625, 2.661, 2.698, 2.734, 2.771, 2.807, 2.844,
            2.880, 2.917, 2.954, 2.990, 3.027, 3.045
        ]) * 1e3,
        "density_data": np.array([
            890, 884, 877, 870, 863, 857, 850, 843, 836, 830, 823, 816, 810,
            803, 796, 789, 783, 776, 769, 763, 756, 749, 742, 736, 729, 722,
            715, 709, 702, 695, 689, 682, 678
        ]),
        "viscosity_data": np.array([
            160, 70.0, 37.0, 23.0, 15.0, 11.0, 8.10, 6.30, 5.00, 4.10, 3.40,
            2.90, 2.40, 2.20, 1.90, 1.70, 1.50, 1.30, 1.20, 1.10, 1.00, 0.88,
            0.81, 0.75, 0.69, 0.64, 0.59, 0.55, 0.51, 0.48, 0.45, 0.42, 0.41
        ]) * 1e-3,
        "conductivity_data": np.array([
            0.1361, 0.1354, 0.1347, 0.1340, 0.1332, 0.1325, 0.1318, 0.1310,
            0.1303, 0.1295, 0.1287, 0.1280, 0.1272, 0.1264, 0.1256, 0.1248,
            0.1240, 0.1232, 0.1224, 0.1216, 0.1208, 0.1200, 0.1191, 0.1183,
            0.1174, 0.1165, 0.1157, 0.1149, 0.1140, 0.1131, 0.1123, 0.1114,
            0.1109
        ])
    }


@fixture
def wrapper_instance(property_data):
    return IncompressibleFluidWrapper(
        "fluid name",
        None,
        **property_data
    )


def test_setup_wrapper(wrapper_instance):
    assert wrapper_instance is not None


def test_enthalpy_forwards_backwards(property_data, wrapper_instance):
    temperature_data = property_data["temperature_data"]
    temperature_check = []

    for temperature in temperature_data:
        h = wrapper_instance.h_pT(1e5, temperature)
        temperature_check.append(wrapper_instance.T_ph(1e5, h))

    np.testing.assert_allclose(temperature_check, temperature_data)


def test_enthalpy_forwards_backwards_constant_heat_capacity():
    # constant data produce a heat capacity polynomial slope of (nearly)
    # zero, the standard quadratic formula in the enthalpy inversion loses
    # all precision in this case
    wrapper = IncompressibleFluidWrapper(
        "fluid name",
        None,
        temperature_data=np.array([273.15, 673.15]),
        heat_capacity_data=np.array([1000.0, 1000.0]),
        density_data=np.array([0.9, 0.9]),
        viscosity_data=np.array([3e-5, 3e-5])
    )
    temperature_data = np.array([273.15, 373.15, 573.15, 673.15])
    temperature_check = [
        wrapper.T_ph(1e5, wrapper.h_pT(1e5, temperature))
        for temperature in temperature_data
    ]

    np.testing.assert_allclose(temperature_check, temperature_data)


def test_entropy_forwards_backwards(property_data, wrapper_instance):
    temperature_data = property_data["temperature_data"]
    temperature_check = []

    for temperature in temperature_data:
        s = wrapper_instance.s_pT(1e5, temperature)
        temperature_check.append(wrapper_instance.T_ps(1e5, s))

    np.testing.assert_allclose(temperature_check, temperature_data)


def test_entropy_enthalpy_roundtrip(property_data, wrapper_instance):
    temperature_data = property_data["temperature_data"]
    temperature_check = []

    for temperature in temperature_data:
        s = wrapper_instance.s_pT(1e5, temperature)
        h = wrapper_instance.h_ps(1e5, s)
        temperature_check.append(wrapper_instance.T_ph(1e5, h))

    np.testing.assert_allclose(temperature_check, temperature_data)


def test_isentropic(property_data, wrapper_instance):
    temperature_data = property_data["temperature_data"]

    enthalpy_inflow = wrapper_instance.h_pT(1e5, temperature_data)
    enthalpy_outflow = wrapper_instance.isentropic(1e5, enthalpy_inflow, 2e5)
    rho = wrapper_instance.d_pT(None, temperature_data)

    np.testing.assert_allclose(
        enthalpy_outflow - enthalpy_inflow, 1e5 / rho
    )


def test_density(property_data, wrapper_instance):

    density_data = property_data["density_data"]
    temperature_data = property_data["temperature_data"]

    density = wrapper_instance.d_pT(None, temperature_data)
    np.testing.assert_allclose(density, density_data, rtol=1e-3)


def test_heat_capacity(property_data, wrapper_instance):

    capacity_data = property_data["heat_capacity_data"]
    temperature_data = property_data["temperature_data"]
    # capacity is not implemented as method in wrapper, using enthalpy over
    # temperature differences
    d = 1e-3
    capacity = [
        (h_u - h_l) / (2 * d)
        for h_l, h_u in
        zip(
            wrapper_instance.h_pT(1e5, temperature_data - d),
            wrapper_instance.h_pT(1e5, temperature_data + d)
        )
    ]
    np.testing.assert_allclose(capacity, capacity_data, rtol=1e-3)


def test_viscosity(property_data, wrapper_instance):

    viscosity_data = property_data["viscosity_data"]
    temperature_data = property_data["temperature_data"]

    viscosity = wrapper_instance.viscosity_pT(None, temperature_data)
    # the data span more than two orders of magnitude, the largest deviation
    # of the fit is at the lowest temperature
    np.testing.assert_allclose(viscosity, viscosity_data, rtol=1e-1)


def test_conductivity(property_data, wrapper_instance):

    conductivity_data = property_data["conductivity_data"]
    temperature_data = property_data["temperature_data"]

    conductivity = wrapper_instance.conductivity_pT(None, temperature_data)
    # the data are slightly curved, the fit is linear
    np.testing.assert_allclose(conductivity, conductivity_data, rtol=5e-3)


@mark.parametrize("pressure", [1e6, 1e7])
def test_entropy_forwards_backwards_off_reference_pressure(
        property_data, wrapper_instance, pressure
):
    temperature_data = property_data["temperature_data"]
    temperature_check = [
        wrapper_instance.T_ps(pressure, wrapper_instance.s_pT(pressure, T))
        for T in temperature_data
    ]

    np.testing.assert_allclose(temperature_check, temperature_data)


@mark.parametrize("pressure", [1e6, 1e7])
def test_enthalpy_forwards_backwards_off_reference_pressure(
        property_data, wrapper_instance, pressure
):
    temperature_data = property_data["temperature_data"]
    temperature_check = [
        wrapper_instance.T_ph(pressure, wrapper_instance.h_pT(pressure, T))
        for T in temperature_data
    ]

    np.testing.assert_allclose(temperature_check, temperature_data)


def test_entropy_pressure_derivative(property_data, wrapper_instance):
    # Maxwell relation (ds/dp)_T = -(dv/dT)_p = (drho/dT)_p / rho ** 2
    temperature_data = property_data["temperature_data"][1:-1]
    p_1, p_2 = 1e5, 1e6
    d = 1e-3

    slope = (
        wrapper_instance.s_pT(p_2, temperature_data)
        - wrapper_instance.s_pT(p_1, temperature_data)
    ) / (p_2 - p_1)

    rho = wrapper_instance.d_pT(None, temperature_data)
    drho_dT = (
        wrapper_instance.d_pT(None, temperature_data + d)
        - wrapper_instance.d_pT(None, temperature_data - d)
    ) / (2 * d)

    np.testing.assert_allclose(slope, drho_dT / rho ** 2, rtol=1e-6)


def test_enthalpy_pressure_derivative(property_data, wrapper_instance):
    # (dh/dp)_T = v - T * (dv/dT)_p from dh = T ds + v dp and the Maxwell
    # relation
    temperature_data = property_data["temperature_data"][1:-1]
    p_1, p_2 = 1e5, 1e6
    d = 1e-3

    slope = (
        wrapper_instance.h_pT(p_2, temperature_data)
        - wrapper_instance.h_pT(p_1, temperature_data)
    ) / (p_2 - p_1)

    v = 1 / wrapper_instance.d_pT(None, temperature_data)
    dv_dT = (
        1 / wrapper_instance.d_pT(None, temperature_data + d)
        - 1 / wrapper_instance.d_pT(None, temperature_data - d)
    ) / (2 * d)

    np.testing.assert_allclose(slope, v - temperature_data * dv_dT, rtol=1e-6)


def test_isentropic_enthalpy_change_via_entropy(property_data, wrapper_instance):
    # along an isentrope dh = v dp, v barely changes for a liquid
    p_1, p_2 = 1e5, 1e7

    for temperature in property_data["temperature_data"][1:-1]:
        h_1 = wrapper_instance.h_pT(p_1, temperature)
        h_2 = wrapper_instance.h_ps(p_2, wrapper_instance.s_pT(p_1, temperature))
        v = 1 / wrapper_instance.d_pT(None, temperature)

        assert h_2 - h_1 == approx(v * (p_2 - p_1), rel=1e-2)


def test_isentropic_mixture(property_data, wrapper_instance):
    heavy_data = property_data.copy()
    heavy_data["density_data"] = property_data["density_data"] * 1.1
    heavy_wrapper = IncompressibleFluidWrapper("heavy fluid", None, **heavy_data)
    fluid_data = {
        "fluid name": {"wrapper": wrapper_instance, "mass_fraction": 0.4},
        "heavy fluid": {"wrapper": heavy_wrapper, "mass_fraction": 0.6},
    }
    p_1, p_2 = 1e5, 1e7

    for temperature in property_data["temperature_data"][1:-1]:
        h_1 = h_mix_pT(p_1, temperature, fluid_data, "incompressible")
        h_2 = isentropic(p_1, h_1, p_2, fluid_data, "incompressible")
        v = v_mix_pT(p_1, temperature, fluid_data, "incompressible")

        s_1 = s_mix_ph(p_1, h_1, fluid_data, "incompressible")
        s_2 = s_mix_ph(p_2, h_2, fluid_data, "incompressible")

        assert s_2 == approx(s_1, abs=1e-6)
        assert h_2 - h_1 == approx(v * (p_2 - p_1), rel=1e-2)
