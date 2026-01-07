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

from tespy.tools.fluid_properties.wrappers import IncompressibleFluidWrapper


@fixture
def property_data():
    # sample data from a publicly available datasheet
    # https://petrocanadalubricants.com/api/sitecore/lubesapi/downloadresource?docID=IM-7852E&type=TechData&lang=english&name=CALFLO%20AF
    return {
        "temperature_data": np.array([292.647, 310.808, 366.241, 421.673, 477.108, 532.542, 588.826, 618.580]),
        "heat_capacity_data": np.array([1901.775, 1961.529, 2143.908, 2326.287, 2508.674, 2691.060, 2876.242, 2974.135]) * 1000,
        "density_data": np.array([863.811, 852.596, 818.368, 784.139, 749.909, 715.678, 680.924, 662.551]),
        "viscosity_data": np.array([0.050335, 0.028525, 0.007075, 0.002500, 0.00111, 0.000579, 0.000334, 0.000259])
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
        h = wrapper_instance.h_pT(None, temperature)
        temperature_check.append(wrapper_instance.T_ph(None, h))

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

    enthalpy_inflow = wrapper_instance.h_pT(None, temperature_data)
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
            wrapper_instance.h_pT(None, temperature_data - d),
            wrapper_instance.h_pT(None, temperature_data + d)
        )
    ]
    np.testing.assert_allclose(capacity, capacity_data, rtol=1e-3)


def test_viscosity(property_data, wrapper_instance):

    viscosity_data = property_data["viscosity_data"]
    temperature_data = property_data["temperature_data"]

    viscosity = wrapper_instance.viscosity_pT(None, temperature_data)
    # allow higher tolerance for viscosity
    np.testing.assert_allclose(viscosity, viscosity_data, rtol=1e-2)
