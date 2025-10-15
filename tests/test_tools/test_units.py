# -*- coding: utf-8

"""Module for testing helper functions.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_tools/test_units.py

SPDX-License-Identifier: MIT
"""
from pint import UnitRegistry
from pytest import raises

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
