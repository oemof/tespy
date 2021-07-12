# -*- coding: utf-8

"""Module for testing reading charactersitic lines.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_tools/test_characteristics.py

SPDX-License-Identifier: MIT
"""
import json
import os
import shutil

import numpy as np
from pkg_resources import resource_filename

from tespy.tools.characteristics import CharLine
from tespy.tools.characteristics import CharMap
from tespy.tools.characteristics import load_custom_char
from tespy.tools.characteristics import load_default_char
from tespy.tools.helpers import extend_basic_path


def test_custom_CharLine_import():
    """Test importing a custom characteristc lines."""

    # we need to write some data to the path first, using defaults
    data_path = resource_filename('tespy.data', 'char_lines.json')
    path = extend_basic_path('data')
    tmp_path = extend_basic_path('tmp_dir_for_testing')

    if os.path.exists(path):
        for f in os.listdir(path):
            shutil.copy(src=path + '/' + f, dst=tmp_path)

    with open(data_path) as f:
        raw_data = json.loads(f.read())

    data = raw_data['heat exchanger']['kA_char2']
    with open(os.path.join(path, 'char_lines.json'), 'w') as outfile:
        json.dump(data, outfile)

    char_original = load_default_char('heat exchanger', 'kA_char2',
                                      'EVAPORATING FLUID', CharLine)
    char_custom = load_custom_char('EVAPORATING FLUID', CharLine)

    shutil.rmtree(path, ignore_errors=True)

    if os.path.exists(tmp_path):
        path = extend_basic_path('data')
        for f in os.listdir(tmp_path):
            shutil.copy(src=tmp_path + '/' + f, dst=path)

        shutil.rmtree(tmp_path, ignore_errors=True)

    x_cond = np.array_equal(char_original.x, char_custom.x)
    y_cond = np.array_equal(char_original.y, char_custom.y)

    msg = ('The x values from the custom characteristic line ' +
           str(char_custom.x) + ' must be identical to the x values from '
           'the default characteristic line ' + str(char_original.x) + ' '
           'as these have been duplicated before load.')
    assert x_cond, msg

    msg = ('The y values from the custom characteristic line ' +
           str(char_custom.y) + ' must be identical to the y values from '
           'the default characteristic line ' + str(char_original.y) + ' '
           'as these have been duplicated before load.')
    assert y_cond, msg


def test_custom_CharMap_import():
    """Test importing a custom characteristc map."""

    # we need to write some data to the path first, using defaults
    data_path = resource_filename('tespy.data', 'char_maps.json')
    path = extend_basic_path('data')
    tmp_path = extend_basic_path('tmp_dir_for_testing')

    if os.path.exists(path):
        for f in os.listdir(path):
            shutil.copy(src=path + '/' + f, dst=tmp_path)

    with open(data_path) as f:
        raw_data = json.loads(f.read())

    data = raw_data['compressor']['char_map_pr']
    with open(os.path.join(path, 'char_maps.json'), 'w') as outfile:
        json.dump(data, outfile)

    char_original = load_default_char(
        'compressor', 'char_map_pr', 'DEFAULT', CharMap)
    char_custom = load_custom_char('DEFAULT', CharMap)

    x_cond = np.array_equal(char_original.x, char_custom.x)
    y_cond = np.array_equal(char_original.y, char_custom.y)
    z_cond = np.array_equal(char_original.z, char_custom.z)

    shutil.rmtree(path, ignore_errors=True)

    if os.path.exists(tmp_path):
        path = extend_basic_path('data')
        for f in os.listdir(tmp_path):
            shutil.copy(src=tmp_path + '/' + f, dst=path)

        shutil.rmtree(tmp_path, ignore_errors=True)

    msg = ('The x values from the custom characteristic line ' +
           str(char_custom.x) + ' must be identical to the x values from '
           'the default characteristic line ' + str(char_original.x) + ' '
           'as these have been duplicated before load.')
    assert x_cond, msg

    msg = ('The y values from the custom characteristic line ' +
           str(char_custom.y) + ' must be identical to the y values from '
           'the default characteristic line ' + str(char_original.y) + ' '
           'as these have been duplicated before load.')
    assert y_cond, msg

    msg = ('The z values from the custom characteristic line ' +
           str(char_custom.z) + ' must be identical to the z values from '
           'the default characteristic line ' + str(char_original.z) + ' '
           'as these have been duplicated before load.')
    assert z_cond, msg


def test_CharLine_evaluation():
    """Test the characteristc line evaluation."""

    # create a characteristc line with values of y=(x-2)^2
    line = CharLine(x=[0, 1, 2, 3, 4], y=[4, 1, 0, 1, 4])

    # test evaluation at given x value (x=3, y=1)
    x = 3
    y = line.evaluate(x)
    msg = ('The evaluation of x=' + str(x) + ' must be 1.0, but is ' +
           str(y) + '.')
    assert y == 1.0, msg

    # test evaluation at x=0.5 to force interpolation, result: y=2.5
    x = 0.5
    y = line.evaluate(x)
    msg = ('The evaluation of x=' + str(x) + ' must be 2.5, but is ' +
           str(y) + '.')
    assert y == 2.5, msg

    # test evaluation at x=-1 to check lower limits, result: y=4
    x = -1
    y = line.evaluate(x)
    msg = ('The evaluation of x=' + str(x) + ' must be 4, but is ' +
           str(y) + '.')
    assert y == 4.0, msg

    # test evaluation at x=5 to check upper limits, result: y=4
    x = 5
    y = line.evaluate(x)
    msg = ('The evaluation of x=' + str(x) + ' must be 4, but is ' +
           str(y) + '.')
    assert y == 4.0, msg


def test_CharLine_extrapolation():
    """Test the characteristc line with extrapolation."""

    # create a characteristc line with values of y=(x-2)^2
    line = CharLine(x=[0, 1, 2, 3, 4], y=[4, 1, 0, 1, 4], extrapolate=True)

    # test evaluation at x=-1 to check lower limits, result: y=7
    x = -1
    y = line.evaluate(x)
    msg = ('The evaluation of x=' + str(x) + ' must be 7, but is ' +
           str(y) + '.')
    assert y == 7.0, msg

    # test evaluation at x=5 to check upper limits, result: y=7
    x = 5
    y = line.evaluate(x)
    msg = ('The evaluation of x=' + str(x) + ' must be 7, but is ' +
           str(y) + '.')
    assert y == 7.0, msg


def test_CharMap_evaluation():
    """Test the characteristc map evaluation."""

    # create a characteristc line with values of y=(x-2)^2
    x = [1, 2, 3]
    y = np.array([[1, 2, 3], [2, 3, 4], [3, 4, 5]])
    z = y ** 0.5
    map = CharMap(x=x, y=y, z=z)

    # test evaluation at x=2 and y=3, result: z1=1.73, z2=9
    x = 2
    y = 3
    z = map.evaluate(x=x, y=y)
    msg = ('The evaluation of x=' + str(x) + ' and y=' + str(y) + ' for z '
           'must be 1.73, but is ' + str(round(z, 2)) + '.')
    assert round(z, 2) == 1.73, msg

    # test evaluation at x=0 and y=0 for lower value range limit,
    # result: z=1
    x = 0
    y = 0
    z = map.evaluate(x=x, y=y)
    msg = ('The evaluation of x=' + str(x) + ' and y=' + str(y) + ' for z '
           'must be 1.0, but is ' + str(round(z, 1)) + '.')
    assert round(z, 1) == 1.0, msg

    # test evaluation at x=4 and y=6 for upper value range limit,
    # result: z=2.24
    x = 4
    y = 6
    z = map.evaluate(x=x, y=y)
    msg = ('The evaluation of x=' + str(x) + ' and y=' + str(y) + ' for z '
           'must be 2.24, but is ' + str(round(z, 2)) + '.')
    assert round(z, 2) == 2.24, msg

    # check, if bound errors go through
    map.get_domain_errors(x, y, 'Componentlabel')
