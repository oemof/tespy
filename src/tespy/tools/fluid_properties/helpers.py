# -*- coding: utf-8

"""Module for fluid property helper functions.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/tools/fluid_properties/helpers.py

SPDX-License-Identifier: MIT
"""

import math

import CoolProp.CoolProp as CP
import numpy as np

from tespy.tools.global_vars import ERR
from tespy.tools.helpers import central_difference
from tespy.tools.helpers import newton_with_kwargs
from tespy.tools.logger import logger


def _is_larger_than_precision(value):
    return value > ERR


def _check_mixing_rule(mixing_rule, mixing_functions, propertyfunction):
    if mixing_rule not in mixing_functions:
        msg = (
            f"The mixing rule '{mixing_rule}' is not available for "
            f"the fluid property functions for {propertyfunction}. Available "
            f"rules are '" + "', '".join(mixing_functions.keys()) + "'."
        )
        logger.exception(msg)
        raise KeyError(msg)


def get_number_of_fluids(fluid_data):
    return sum([1 for f in fluid_data.values() if _is_larger_than_precision(f["mass_fraction"])])


def get_pure_fluid(fluid_data):
    for f in fluid_data.values():
        if _is_larger_than_precision(f["mass_fraction"]):
            return f


def single_fluid(fluid_data):
    r"""Return the name of the pure fluid in a fluid vector."""
    if get_number_of_fluids(fluid_data) > 1:
        return None
    else:
        for fluid, data in fluid_data.items():
            if _is_larger_than_precision(data["mass_fraction"]):
                return fluid


def get_molar_fractions(fluid_data):
    molarflow = {
        key: value["mass_fraction"] / value["wrapper"]._molar_mass
        for key, value in fluid_data.items()
    }
    molarflow_sum = sum(molarflow.values())
    return {key: value / molarflow_sum for key, value in molarflow.items()}


def inverse_temperature_mixture(p=None, target_value=None, fluid_data=None, T0=None, f=None):
    # calculate the fluid properties for fluid mixtures
    valmin, valmax = get_mixture_temperature_range(fluid_data)
    if T0 is None or T0 == 0 or np.isnan(T0):
        T0 = (valmin + valmax) / 2.0

    valmax *= 2

    function_kwargs = {
        "p": p, "fluid_data": fluid_data, "T": T0,
        "function": f, "parameter": "T" , "delta": 0.01
    }
    return newton_with_kwargs(
        central_difference,
        target_value,
        val0=T0,
        valmin=valmin,
        valmax=valmax,
        **function_kwargs
    )


def get_mixture_temperature_range(fluid_data):
    valmin = max(
        [v["wrapper"]._T_min for v in fluid_data.values() if _is_larger_than_precision(v["mass_fraction"])]
    ) + 0.1
    valmax = min(
        [v["wrapper"]._T_max for v in fluid_data.values() if _is_larger_than_precision(v["mass_fraction"])]
    ) - 0.1
    return valmin, valmax


def calc_molar_mass_mixture(fluid_data, molar_fractions):
    return sum([x * fluid_data[fluid]["wrapper"]._molar_mass for fluid, x in molar_fractions.items()])


def fluid_structure(fluid):
    r"""
    Return the checmical formula of fluid.

    Parameters
    ----------
    fluid : str
        Name of the fluid.

    Returns
    -------
    parts : dict
        Dictionary of the chemical base elements as keys and the number of
        atoms in a molecule as values.

    Example
    -------
    Get the chemical formula of methane.

    >>> from tespy.tools.fluid_properties.helpers import fluid_structure
    >>> elements = fluid_structure('methane')
    >>> elements['C'], elements['H']
    (1, 4)
    """
    parts = {}
    for element in CP.get_fluid_param_string(
            fluid, 'formula').split('}'):
        if element != '':
            el = element.split('_{')
            parts[el[0]] = int(el[1])

    return parts


def darcy_friction_factor(re, ks, d):
    r"""
    Calculate the Darcy friction factor.

    Parameters
    ----------
    re : float
        Reynolds number re / 1.

    ks : float
        Pipe roughness ks / m.

    d : float
        Pipe diameter/characteristic lenght d / m.

    Returns
    -------
    darcy_friction_factor : float
        Darcy friction factor :math:`\lambda` / 1

    Note
    ----
    **Laminar flow** (:math:`re \leq 2320`)

    .. math::

        \lambda = \frac{64}{re}

    **turbulent flow** (:math:`re > 2320`)

    *hydraulically smooth:* :math:`\frac{re \cdot k_{s}}{d} < 65`

    .. math::

        \lambda = \begin{cases}
        0.03164 \cdot re^{-0.25} & re \leq 10^4\\
        \left(1.8 \cdot \log \left(re\right) -1.5 \right)^{-2} &
        10^4 < re < 10^6\\
        solve \left(0 = 2 \cdot \log\left(re \cdot \sqrt{\lambda} \right) -0.8
        - \frac{1}{\sqrt{\lambda}}\right) & re \geq 10^6\\
        \end{cases}

    *transition zone and hydraulically rough:*

    .. math::

        \lambda = solve \left( 0 = 2 \cdot \log \left( \frac{2.51}{re \cdot
        \sqrt{\lambda}} + \frac{k_{s}}{d \cdot 3.71} \right) -
        \frac{1}{\sqrt{\lambda}} \right)

    Reference: :cite:`Nirschl2018`.

    Example
    -------
    Calculate the Darcy friction factor at different hydraulic states.

    >>> from tespy.tools.fluid_properties.helpers import darcy_friction_factor
    >>> ks = 5e-5
    >>> d = 0.05
    >>> re_laminar = 2000
    >>> re_turb_smooth = 5000
    >>> re_turb_trans = 70000
    >>> re_high = 1000000
    >>> d_high = 0.8
    >>> re_very_high = 6000000
    >>> d_very_high = 1
    >>> ks_low = 1e-5
    >>> ks_rough = 1e-3
    >>> darcy_friction_factor(re_laminar, ks, d)
    0.032
    >>> round(darcy_friction_factor(re_turb_smooth, ks, d), 3)
    0.038
    >>> round(darcy_friction_factor(re_turb_trans, ks, d), 3)
    0.023
    >>> round(darcy_friction_factor(re_turb_trans, ks_rough, d), 3)
    0.049
    >>> round(darcy_friction_factor(re_high, ks, d_high), 3)
    0.012
    >>> round(darcy_friction_factor(re_very_high, ks_low, d_very_high), 3)
    0.009
    """
    if re <= 2320:
        return 64 / re
    else:
        if re * ks / d < 65:
            if re <= 1e4:
                return blasius(re)
            elif re < 1e6:
                return hanakov(re)
            else:
                l0 = 0.02
                function_kwargs = {
                    "function": prandtl_karman,
                    "parameter": "darcy_friction_factor",
                    "reynolds": re

                }
                return newton_with_kwargs(
                    prandtl_karman_derivative,
                    0,
                    val0=l0,
                    valmin=0.00001,
                    valmax=0.2,
                    **function_kwargs
                )

        else:
            l0 = 0.002
            function_kwargs = {
                "function": colebrook,
                "parameter": "darcy_friction_factor",
                "reynolds": re,
                "ks": ks,
                "diameter": d,
                "delta": 0.001
            }
            return newton_with_kwargs(
                central_difference,
                0,
                val0=l0,
                valmin=0.0001,
                valmax=0.2,
                **function_kwargs
            )


def blasius(re):
    """
    Calculate friction coefficient according to Blasius.

    Parameters
    ----------
    re : float
        Reynolds number.

    Returns
    -------
    darcy_friction_factor : float
        Darcy friction factor.
    """
    return 0.3164 * re ** (-0.25)


def hanakov(re):
    """
    Calculate friction coefficient according to Hanakov.

    Parameters
    ----------
    re : float
        Reynolds number.

    Returns
    -------
    darcy_friction_factor : float
        Darcy friction factor.
    """
    return (1.8 * math.log10(re) - 1.5) ** (-2)


def prandtl_karman(reynolds, darcy_friction_factor, **kwargs):
    """
    Calculate friction coefficient according to Prandtl and v. Kármán.

    Applied in smooth conditions.

    Parameters
    ----------
    re : float
        Reynolds number.

    darcy_friction_factor : float
        Darcy friction factor.

    Returns
    -------
    darcy_friction_factor : float
        Darcy friction factor.
    """
    return (
        2 * math.log10(reynolds * darcy_friction_factor ** 0.5)
        - 0.8 - 1 / darcy_friction_factor ** 0.5
    )


def prandtl_karman_derivative(reynolds, darcy_friction_factor, **kwargs):
    """Calculate derivative for Prandtl and v. Kármán equation."""
    return (
        1 / (darcy_friction_factor * math.log(10))
        + 0.5 * darcy_friction_factor ** (-1.5)
    )


def colebrook(reynolds, ks, diameter, darcy_friction_factor, **kwargs):
    """
    Calculate friction coefficient accroding to Colebrook-White equation.

    Applied in transition zone and rough conditions.

    Parameters
    ----------
    re : float
        Reynolds number.

    ks : float
        Equivalent sand roughness.

    d : float
        Pipe's diameter.

    darcy_friction_factor : float
        Darcy friction factor.

    Returns
    -------
    darcy_friction_factor : float
        Darcy friction factor.
    """
    return (
        2 * math.log10(
            2.51 / (reynolds * darcy_friction_factor ** 0.5) + ks
            / (3.71 * diameter)
        ) + 1 / darcy_friction_factor ** 0.5
    )
