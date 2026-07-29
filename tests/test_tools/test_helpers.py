# -*- coding: utf-8

"""Module for testing helper functions.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_tools/test_helpers.py

SPDX-License-Identifier: MIT
"""
from pytest import approx

from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import ScalarVariable as dc_scavar
from tespy.tools.helpers import _get_dependents
from tespy.tools.helpers import newton_with_kwargs


def func(x, **kwargs):
    return x ** 2 + x - 20


def deriv(x, **kwargs):
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
    kwargs = {"function": func, "parameter": "x"}
    result = newton_with_kwargs(deriv, 0, valmin=-10, valmax=10, val0=0, **kwargs)
    msg = ('The newton algorithm should find the zero crossing at 4.0. ' +
           str(round(result, 1)) + ' was found instead.')
    assert 4.0 == approx(result), msg

    result = newton_with_kwargs(deriv, 0, valmin=-10, valmax=10, val0=-10, **kwargs)
    msg = ('The newton algorithm should find the zero crossing at -5.0. ' +
           str(round(result, 1)) + ' was found instead.')
    assert -5.0 == approx(result), msg

    result = newton_with_kwargs(deriv, 0, valmin=-4, valmax=-2, val0=-3, **kwargs)
    msg = ('The newton algorithm should not be able to find a zero crossing. '
           'The value ' + str(round(result, 1)) + ' was found, but the '
           'algorithm should have found the lower boundary of -4.0.')
    assert -4.0 == approx(result), msg

    result = newton_with_kwargs(deriv, 0, valmin=-20, valmax=-10, val0=-10, **kwargs)
    msg = ('The newton algorithm should not be able to find a zero crossing. '
           'The value ' + str(round(result, 1)) + ' was found, but the '
           'algorithm should have found the upper boundary of -10.0.')
    assert -10.0 == approx(result), msg


def test_get_dependents_excludes_dangling_containers():
    """A container with the variable flag but without a reference container
    belongs to an object outside of the network, e.g. a UserDefinedVariable
    that was removed while an equation still lists it as dependent. It must
    not enter the incidence."""
    reference = dc_scavar(_is_var=True)
    attached = dc_cp(_is_var=True)
    attached._reference_container = reference
    dangling = dc_cp(_is_var=True)
    not_a_variable = dc_cp(_is_var=False)

    flat = _get_dependents([attached, dangling, not_a_variable])
    assert flat == [{reference}]

    nested = _get_dependents([[attached], [dangling, not_a_variable]])
    assert nested == [{reference}, set()]
