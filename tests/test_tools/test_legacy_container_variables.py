# -*- coding: utf-8

"""Module for testing variable declarations with the deprecated containers.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_tools/test_legacy_container_variables.py

SPDX-License-Identifier: MIT
"""
import pytest

from tespy.tools.data_containers import ComponentParameter
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import FluidProperties as dc_prop


def test_potential_var_declares_variable_parameter():
    with pytest.warns(FutureWarning, match="ComponentProperties"):
        parameter = dc_cp(min_val=15, _potential_var=True)

    assert isinstance(parameter, ComponentParameter)
    assert parameter._allows_var

    parameter.accept("var")
    assert parameter.is_set
    assert parameter.is_var


def test_parameter_without_potential_var_cannot_be_variable():
    with pytest.warns(FutureWarning, match="ComponentProperties"):
        parameter = dc_cp(min_val=15)

    assert not parameter._allows_var
    with pytest.raises(TypeError, match="result-only"):
        parameter.accept("var")


def test_potential_var_kept_on_connection_properties():
    with pytest.warns(FutureWarning, match="FluidProperties"):
        parameter = dc_prop(_potential_var=True)

    assert parameter._potential_var
