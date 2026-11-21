# -*- coding: utf-8

"""Module for testing helper functions.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_tools/test_helpers.py

SPDX-License-Identifier: MIT
"""
from pytest import approx
from pytest import mark

from tespy.tools.global_vars import FLUID_ALIASES
from tespy.tools.helpers import fluidalias_in_list
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
    msg = (
        "The newton algorithm should find the zero crossing at 4.0. "
        f"{round(result, 1)} was found instead."
    )
    assert 4.0 == approx(result), msg

    result = newton_with_kwargs(deriv, 0, valmin=-10, valmax=10, val0=-10, **kwargs)
    msg = (
        "The newton algorithm should find the zero crossing at -5.0. "
        f"{round(result, 1)} was found instead."
    )
    assert -5.0 == approx(result), msg

    result = newton_with_kwargs(deriv, 0, valmin=-4, valmax=-2, val0=-3, **kwargs)
    msg = (
        "The newton algorithm should not be able to find a zero crossing. The "
        f"value {round(result, 1)} was found, but the algorithm should "
        "have found the lower boundary of -4.0."
    )
    assert -4.0 == approx(result), msg

    result = newton_with_kwargs(deriv, 0, valmin=-20, valmax=-10, val0=-10, **kwargs)
    msg = (
        "The newton algorithm should not be able to find a zero crossing. The "
        f"value {round(result, 1)} was found, but the algorithm should "
        "have found the upper boundary of -10.0."
    )
    assert -10.0 == approx(result), msg


@mark.parametrize("names", [
    ["n-Hexane", "nHexane", "Hexane"],
    ["Methane", "CH4", "methane"],
    ["Water", "H2O", "water"],
    ["n-Butane", "butane", "R600"],
])
def test_fluid_aliases_are_the_same_for_every_name(names):
    """The CoolProp primary name is an alias like every other name."""
    aliases = [FLUID_ALIASES.get_fluid(name) for name in names]
    assert all(alias == aliases[0] for alias in aliases)
    assert set(names) <= aliases[0]


@mark.parametrize("fluid,fluid_list", [
    ("n-Hexane", ["n-Hexane"]),
    ("nHexane", ["n-Hexane"]),
    ("n-Hexane", ["nHexane"]),
    ("Methane", ["Methane"]),
])
def test_fluidalias_in_list_with_primary_name(fluid, fluid_list):
    assert fluidalias_in_list(fluid, fluid_list) == set(fluid_list)


def test_fluid_aliases_of_unknown_fluid():
    assert FLUID_ALIASES.get_fluid("not a CoolProp fluid") == {"not a CoolProp fluid"}
