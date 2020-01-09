# -*- coding: utf-8

from nose.tools import eq_

from tespy.tools.helpers import newton


def func(params, x):
    return x ** 2 + x - 20


def deriv(params, x):
    return 2 * x + 1


def test_newton_bounds():
    """
    Test newton algorithm value limit handling.

    Try to calculate a zero crossing of a quadratic function in three
    tries.

    - zero crossing within limits, starting value near 4
    - zero crossing within limits, starting value near -5

    - zero crossing below minimum
    - zero crossing above maximum

    The function is x^2 + x - 20, there crossings are -5 and 4.
    """
    result = newton(func, deriv, [], 0, valmin=-10, valmax=10, val0=0)
    msg = ('The newton algorithm should find the zero crossing at 4.0. ' +
           str(round(result, 1)) + ' was found instead.')
    eq_(4.0, result, msg)

    result = newton(func, deriv, [], 0, valmin=-10, valmax=10, val0=-10)
    msg = ('The newton algorithm should find the zero crossing at -5.0. ' +
           str(round(result, 1)) + ' was found instead.')
    eq_(-5.0, result, msg)

    result = newton(func, deriv, [], 0, valmin=-4, valmax=-2, val0=-3)
    msg = ('The newton algorithm should not be able to find a zero crossing. '
           'The value ' + str(round(result, 1)) + ' was found, but the '
           'algorithm should have found the lower boundary of -4.0.')
    eq_(-4.0, result, msg)

    result = newton(func, deriv, [], 0, valmin=-20, valmax=-10, val0=-10)
    msg = ('The newton algorithm should not be able to find a zero crossing. '
           'The value ' + str(round(result, 1)) + ' was found, but the '
           'algorithm should have found the upper boundary of -10.0.')
    eq_(-10.0, result, msg)
