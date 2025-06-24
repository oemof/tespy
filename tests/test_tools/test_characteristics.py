# -*- coding: utf-8

"""Module for testing reading charactersitic lines.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_tools/test_characteristics.py

SPDX-License-Identifier: MIT
"""
import numpy as np

from tespy import __datapath__
from tespy.tools.characteristics import CharLine
from tespy.tools.characteristics import CharMap


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
