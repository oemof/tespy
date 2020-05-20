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

from tespy.tools.characteristics import char_line
from tespy.tools.characteristics import char_map
from tespy.tools.characteristics import compressor_map
from tespy.tools.characteristics import load_custom_char
from tespy.tools.characteristics import load_default_char
from tespy.tools.helpers import extend_basic_path


def test_custom_char_line_import():
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
                                      'EVAPORATING FLUID', char_line)
    char_custom = load_custom_char('EVAPORATING FLUID', char_line)

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
    assert x_cond is True, msg

    msg = ('The y values from the custom characteristic line ' +
           str(char_custom.y) + ' must be identical to the y values from '
           'the default characteristic line ' + str(char_original.y) + ' '
           'as these have been duplicated before load.')
    assert y_cond is True, msg


def test_custom_char_map_import():
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

    data = raw_data['compressor']['char_map']
    with open(os.path.join(path, 'char_maps.json'), 'w') as outfile:
        json.dump(data, outfile)

    char_original = load_default_char('compressor', 'char_map',
                                      'DEFAULT', compressor_map)
    char_custom = load_custom_char('DEFAULT', compressor_map)

    x_cond = np.array_equal(char_original.x, char_custom.x)
    y_cond = np.array_equal(char_original.y, char_custom.y)
    z1_cond = np.array_equal(char_original.z1, char_custom.z1)
    z2_cond = np.array_equal(char_original.z2, char_custom.z2)

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
    assert x_cond is True, msg

    msg = ('The y values from the custom characteristic line ' +
           str(char_custom.y) + ' must be identical to the y values from '
           'the default characteristic line ' + str(char_original.y) + ' '
           'as these have been duplicated before load.')
    assert y_cond is True, msg

    msg = ('The z1 values from the custom characteristic line ' +
           str(char_custom.z1) + ' must be identical to the z1 values from '
           'the default characteristic line ' + str(char_original.z1) + ' '
           'as these have been duplicated before load.')
    assert z1_cond is True, msg

    msg = ('The z2 values from the custom characteristic line ' +
           str(char_custom.z2) + ' must be identical to the z2 values from '
           'the default characteristic line ' + str(char_original.z2) + ' '
           'as these have been duplicated before load.')
    assert z2_cond is True, msg


def test_char_line_evaluation():
    """Test the characteristc line evaluation."""

    # create a characteristc line with values of y=(x-2)^2
    line = char_line(x=[0, 1, 2, 3, 4], y=[4, 1, 0, 1, 4])

    # test evaluation at given x value (x=3, y=1)
    x = 3
    y = line.evaluate(x)
    msg = ("The evaluation of x=" + str(x) + " must be 1.0, but is " +
           str(y) + ".")
    assert y == 1.0, msg

    # test evaluation at x=0.5 to force interpolation, result: y=2.5
    x = 0.5
    y = line.evaluate(x)
    msg = ("The evaluation of x=" + str(x) + " must be 2.5, but is " +
           str(y) + ".")
    assert y == 2.5, msg

    # test evaluation at x=-1 to check lower limits, result: y=4
    x = -1
    y = line.evaluate(x)
    msg = ("The evaluation of x=" + str(x) + " must be 4, but is " +
           str(y) + ".")
    assert y == 4.0, msg

    # test evaluation at x=5 to check upper limits, result: y=4
    x = 5
    y = line.evaluate(x)
    msg = ("The evaluation of x=" + str(x) + " must be 4, but is " +
           str(y) + ".")
    assert y == 4.0, msg


def test_char_line_extrapolation():
    """Test the characteristc line with extrapolation."""

    # create a characteristc line with values of y=(x-2)^2
    line = char_line(x=[0, 1, 2, 3, 4], y=[4, 1, 0, 1, 4], extrapolate=True)

    # test evaluation at x=-1 to check lower limits, result: y=7
    x = -1
    y = line.evaluate(x)
    msg = ("The evaluation of x=" + str(x) + " must be 7, but is " +
           str(y) + ".")
    assert y == 7.0, msg

    # test evaluation at x=5 to check upper limits, result: y=7
    x = 5
    y = line.evaluate(x)
    msg = ("The evaluation of x=" + str(x) + " must be 7, but is " +
           str(y) + ".")
    assert y == 7.0, msg


def test_char_map_evaluation():
    """Test the characteristc map evaluation."""

    # create a characteristc line with values of y=(x-2)^2
    x = [1, 2, 3]
    y = np.array([[1, 2, 3], [2, 3, 4], [3, 4, 5]])
    z1 = y ** 0.5
    z2 = y ** 2
    map = char_map(x=x, y=y, z1=z1, z2=z2)

    # test evaluation at x=2 and y=3, result: z1=1.73, z2=9
    x = 2
    y = 3
    z1, z2 = map.evaluate(x=x, y=y)
    msg = ("The evaluation of x=" + str(x) + " and y=" + str(y) + " for z1 "
           "must be 1.73, but is " + str(round(z1, 2)) + ".")
    assert round(z1, 2) == 1.73, msg

    msg = ("The evaluation of x=" + str(x) + " and y=" + str(y) + " for z2 "
           "must be 9.0, but is " + str(round(z2, 1)) + ".")
    assert round(z2, 1) == 9.0, msg

    # test evaluation at x=0 and y=0 for lower value range limit,
    # result: z1=1, z2=1
    x = 0
    y = 0
    z1, z2 = map.evaluate(x=x, y=y)
    msg = ("The evaluation of x=" + str(x) + " and y=" + str(y) + " for z1 "
           "must be 1.0, but is " + str(round(z1, 1)) + ".")
    assert round(z1, 1) == 1.0, msg

    msg = ("The evaluation of x=" + str(x) + " and y=" + str(y) + " for z2 "
           "must be 1.0, but is " + str(round(z2, 1)) + ".")
    assert round(z2, 1) == 1.0, msg

    # test evaluation at x=4 and y=6 for upper value range limit,
    # result: z1=2.24, z2=25
    x = 4
    y = 6
    z1, z2 = map.evaluate(x=x, y=y)
    msg = ("The evaluation of x=" + str(x) + " and y=" + str(y) + " for z1 "
           "must be 2.24, but is " + str(round(z1, 2)) + ".")
    assert round(z1, 2) == 2.24, msg

    msg = ("The evaluation of x=" + str(x) + " and y=" + str(y) + " for z2 "
           "must be 25.0, but is " + str(round(z2, 1)) + ".")
    assert round(z2, 1) == 25.0, msg

    # check, if bound errors go through
    map.get_bound_errors(x, y, 'Componentlabel')


def test_compressor_map_evaluation():
    """Test the characteristc compressor map evaluation."""

    # create a characteristc line with values of y=(x-2)^2
    x = [1, 2, 3]
    y = np.array([[1, 2, 3], [2, 3, 4], [3, 4, 5]])
    z1 = y ** 0.5
    z2 = y ** 2
    map = compressor_map(x=x, y=y, z1=z1, z2=z2)
    igva = 20

    # test evaluation at x=2 and y=3, result: z1=1.55, z2=13.7
    # f_igva1  = (1 - igva / 100)
    # -->      = 0.8
    # f_igva2  = (1 - igva ** 2 / 10000)
    # -->      = 0.96
    # --> yarr = yarr * f_igva
    # -->      = [1.6, 2.4, 3.2] with y=3
    # --> z1, z2 between second and third value at 0.75 (3-2.4)/(3.2-2.4)
    # --> z1 0.75 * (2 - 3 ** 0.5) * 0.8 + 3 ** 0.5 * 0.8
    # --> z2 0.75 * (16 - 9) * 0.96 + 9 * 0.96
    x = 2
    y = 3
    z1, z2 = map.evaluate(x=x, y=y, igva=igva)
    msg = ("The evaluation of x=" + str(x) + " and y=" + str(y) + " for z1 "
           "must be 1.55, but is " + str(round(z1, 2)) + ".")
    assert round(z1, 2) == 1.55, msg

    msg = ("The evaluation of x=" + str(x) + " and y=" + str(y) + " for z2 "
           "must be 13.7, but is " + str(round(z2, 1)) + ".")
    assert round(z2, 2) == 13.68, msg

    # test evaluation at x=0 and y=0 for lower value range limit,
    # result: z1=0.8, z2=0.96
    x = 0
    y = 0
    z1, z2 = map.evaluate(x=x, y=y, igva=igva)
    msg = ("The evaluation of x=" + str(x) + " and y=" + str(y) + " for z1 "
           "must be 0.8, but is " + str(round(z1, 1)) + ".")
    assert round(z1, 1) == 0.8, msg

    msg = ("The evaluation of x=" + str(x) + " and y=" + str(y) + " for z2 "
           "must be 0.96, but is " + str(round(z2, 1)) + ".")
    assert round(z2, 2) == 0.96, msg

    # test evaluation at x=4 and y=6 for upper value range limit,
    # result: z1=1.79, z2=24
    x = 4
    y = 6
    z1, z2 = map.evaluate(x=x, y=y, igva=igva)
    msg = ("The evaluation of x=" + str(x) + " and y=" + str(y) + " for z1 "
           "must be 1.79, but is " + str(round(z1, 2)) + ".")
    assert round(z1, 2) == 1.79, msg

    msg = ("The evaluation of x=" + str(x) + " and y=" + str(y) + " for z2 "
           "must be 24.0, but is " + str(round(z2, 1)) + ".")
    assert round(z2, 1) == 24.0, msg
