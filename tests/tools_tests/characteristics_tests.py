# -*- coding: utf-8

"""Module for testing reading charactersitic lines.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/tools_tests/characteristics_tests.py

SPDX-License-Identifier: MIT
"""

from nose.tools import eq_

from tespy.tools.characteristics import (char_line, compressor_map,
                                         load_default_char, load_custom_char)

from tespy.tools.helpers import extend_basic_path
from pkg_resources import resource_filename
import os
import json
import numpy as np
import shutil


class characteristics_loader_tests:

    def setup(self):
        # create data path and write json files into path
        self.path = extend_basic_path('data')

    def test_custom_char_line_import(self):
        """Test importing a custom characteristc lines."""

        # we need to write some data to the path first, using defaults
        data_path = resource_filename('tespy.data', 'char_lines.json')

        with open(data_path) as f:
            raw_data = json.loads(f.read())

        data = raw_data['heat exchanger']['kA_char2']
        with open(os.path.join(self.path, 'char_lines.json'), 'w') as outfile:
            json.dump(data, outfile)

        char_original = load_default_char('heat exchanger', 'kA_char2',
                                          'EVAPORATING FLUID', char_line)
        char_custom = load_custom_char('EVAPORATING FLUID', char_line)

        shutil.rmtree(self.path, ignore_errors=True)

        x_cond = np.array_equal(char_original.x, char_custom.x)
        y_cond = np.array_equal(char_original.y, char_custom.y)

        msg = ('The x values from the custom characteristic line ' +
               str(char_custom.x) + ' must be identical to the x values from '
               'the default characteristic line ' + str(char_original.x) + ' '
               'as these have been duplicated before load.')
        eq_(True, x_cond, msg)

        msg = ('The y values from the custom characteristic line ' +
               str(char_custom.y) + ' must be identical to the y values from '
               'the default characteristic line ' + str(char_original.y) + ' '
               'as these have been duplicated before load.')
        eq_(True, y_cond, msg)

    def test_custom_char_map_import(self):
        """Test importing a custom characteristc map."""

        # we need to write some data to the path first, using defaults
        data_path = resource_filename('tespy.data', 'char_maps.json')

        with open(data_path) as f:
            raw_data = json.loads(f.read())

        data = raw_data['compressor']['char_map']
        with open(os.path.join(self.path, 'char_maps.json'), 'w') as outfile:
            json.dump(data, outfile)

        char_original = load_default_char('compressor', 'char_map',
                                          'DEFAULT', compressor_map)
        char_custom = load_custom_char('DEFAULT', compressor_map)

        x_cond = np.array_equal(char_original.x, char_custom.x)
        y_cond = np.array_equal(char_original.y, char_custom.y)
        z1_cond = np.array_equal(char_original.z1, char_custom.z1)
        z2_cond = np.array_equal(char_original.z2, char_custom.z2)

        shutil.rmtree(self.path, ignore_errors=True)

        msg = ('The x values from the custom characteristic line ' +
               str(char_custom.x) + ' must be identical to the x values from '
               'the default characteristic line ' + str(char_original.x) + ' '
               'as these have been duplicated before load.')
        eq_(True, x_cond, msg)

        msg = ('The y values from the custom characteristic line ' +
               str(char_custom.y) + ' must be identical to the y values from '
               'the default characteristic line ' + str(char_original.y) + ' '
               'as these have been duplicated before load.')
        eq_(True, y_cond, msg)

        msg = ('The z1 values from the custom characteristic line ' +
               str(char_custom.z1) + ' must be identical to the z1 values from '
               'the default characteristic line ' + str(char_original.z1) + ' '
               'as these have been duplicated before load.')
        eq_(True, z1_cond, msg)

        msg = ('The z2 values from the custom characteristic line ' +
               str(char_custom.z2) + ' must be identical to the z2 values from '
               'the default characteristic line ' + str(char_original.z2) + ' '
               'as these have been duplicated before load.')
        eq_(True, z2_cond, msg)


def test_characteristic_line():
    """Test the characteristc line without extrapolation."""

    # create a characteristc line with values of y=(x-2)^2
    line = char_line(x=[0, 1, 2, 3, 4], y=[4, 1, 0, 1, 4])

    # test evaluation at given x value (x=3, y=1)
    x = 3
    y = line.evaluate(x)
    msg = ("The evaluation of x=" + str(x) + " must be 1.0, but is" +
           str(y) + ".")
    eq_(1.0, round(y, 1), msg)

    # test evaluation at x=0.5 to force interpolation, result: y=2.5
    x = 0.5
    y = line.evaluate(x)
    msg = ("The evaluation of x=" + str(x) + " must be 2.5, but is" +
           str(y) + ".")
    eq_(2.5, round(y, 1), msg)

    # test evaluation at x=-1 to check lower limits, result: y=4
    x = -1
    y = line.evaluate(x)
    msg = ("The evaluation of x=" + str(x) + " must be 4, but is" +
           str(y) + ".")
    eq_(4, round(y, 1), msg)

    # test evaluation at x=5 to check upper limits, result: y=4
    x = 5
    y = line.evaluate(x)
    msg = ("The evaluation of x=" + str(x) + " must be 4, but is" +
           str(y) + ".")
    eq_(4, round(y, 1), msg)


def test_extrapolation():
    """Test the characteristc line with extrapolation."""

    # create a characteristc line with values of y=(x-2)^2
    line = char_line(x=[0, 1, 2, 3, 4], y=[4, 1, 0, 1, 4], extrapolate=True)

    # test evaluation at given x value (x=3, y=1)
    x = 1
    y = line.evaluate(x)
    msg = ("The evaluation of x=" + str(x) + " must be 1.0, but is" +
           str(y) + ".")
    eq_(1.0, round(y, 1), msg)

    # test evaluation at x=0.5 to force interpolation, result: y=2.5
    x = 0.5
    y = line.evaluate(x)
    msg = ("The evaluation of x=" + str(x) + " must be 2.5, but is" +
           str(y) + ".")
    eq_(2.5, round(y, 1), msg)

    # test evaluation at x=-1 to check lower limits, result: y=7
    x = -1
    y = line.evaluate(x)
    msg = ("The evaluation of x=" + str(x) + " must be 7, but is" +
           str(y) + ".")
    eq_(7.0, round(y, 1), msg)

    # test evaluation at x=5 to check upper limits, result: y=7
    x = 5
    y = line.evaluate(x)
    msg = ("The evaluation of x=" + str(x) + " must be 7, but is" +
           str(y) + ".")
    eq_(7.0, round(y, 1), msg)
