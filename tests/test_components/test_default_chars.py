# -*- coding: utf-8

"""Tests for default characteristic line and map loading with inheritance.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_default_chars.py

SPDX-License-Identifier: MIT
"""
import json
import os

import numpy as np
import pytest

from tespy import __datapath__
from tespy.components.component import component_registry
from tespy.tools.data_containers import ComponentCharacteristicMaps as dc_cm
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc


def _build_cases():
    with open(os.path.join(__datapath__, 'char_lines.json')) as f:
        char_lines = json.load(f)
    with open(os.path.join(__datapath__, 'char_maps.json')) as f:
        char_maps = json.load(f)

    cases = []
    for comp_name, cls in sorted(component_registry.items.items()):
        inst = cls(comp_name)
        for param, data in inst.parameters.items():
            if isinstance(data, dc_cc):
                data_dict = char_lines
                get_last = lambda entry: entry['y'][-1]
                fallback = 1.0
            elif isinstance(data, dc_cm):
                data_dict = char_maps
                get_last = lambda entry: np.asarray(entry['z']).ravel()[-1]
                fallback = 1.0
            else:
                continue

            expected = fallback
            for mro_cls in cls.__mro__:
                if mro_cls.__name__ in data_dict and param in data_dict[mro_cls.__name__]:
                    expected = get_last(data_dict[mro_cls.__name__][param]['DEFAULT'])
                    break

            cases.append(
                pytest.param(cls, param, expected, id=f"{comp_name}.{param}")
            )
    return cases


_CASES = _build_cases()


def _load_default_chars(comp):
    comp.num_vars = 0
    comp.prop_specifications = {}
    comp.var_specifications = {}
    comp.group_specifications = {}
    comp.char_specifications = {}
    comp._rhs = {}
    comp._equation_set_lookup = {}
    comp._mode = 'design'
    comp._setup_user_imposed_constraints(0, 0)


@pytest.mark.parametrize("cls, param, expected_last", _CASES)
def test_default_char_loaded(cls, param, expected_last):
    inst = cls(cls.__name__)
    _load_default_chars(inst)
    data = inst.parameters[param]
    if isinstance(data, dc_cc):
        actual = data.char_func.y[-1]
    else:
        actual = data.char_func.z.ravel()[-1]
    assert actual == pytest.approx(expected_last)
