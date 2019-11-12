# -*- coding: utf-8

"""
.. module:: helpers
    :synopsis: helpers for frequently used functionalities

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>
"""

import CoolProp as CP
from CoolProp.CoolProp import PropsSI as CPPSI

import math
import numpy as np
from scipy import interpolate
import pandas as pd
import os
import collections

import logging

import warnings
warnings.simplefilter("ignore", RuntimeWarning)

global err
err = 1e-6
global molar_masses
molar_masses = {}
global gas_constants
gas_constants = {}
gas_constants['uni'] = 8.3144598

# %%


def molar_mass_flow(flow):
    r"""
    Calculates molar mass flow.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

    Returns
    -------
    mm : float
        Molar mass flow mm / (mol/s).

        .. math::

            mm = \sum_{i} \left( \frac{x_{i}}{M_{i}} \right)
    """
    mm = 0
    for fluid, x in flow.items():
        if x > err:
            mm += x / molar_masses[fluid]
    return mm

# %%


def num_fluids(fluids):
    r"""
    Returns number of fluids in fluid mixture.

    Parameters
    ----------
    fluids : dict
        Fluid mass fractions.

    Returns
    -------
    n : int
        Number of fluids in fluid mixture n / 1.

        .. math::

            n = \sum_{i} \left( \begin{cases}
            0 & x_{i} < \epsilon \\
            1 & x_{i} \geq \epsilon
            \end{cases} \right)\;
            \forall i \in \text{network fluids}
    """
    n = 0
    for fluid, x in fluids.items():
        if x > err:
            n += 1

    return n

# %%


def single_fluid(fluids):
    r"""
    Returns the name of the pure fluid in a fluid vector.

    Parameters
    ----------
    fluids : dict
        Fluid mass fractions.

    Returns
    -------
    fluid : str
        Name of the single fluid.
    """
    if num_fluids(fluids) == 1:
        for fluid, x in fluids.items():
            if x > err:
                return fluid
    else:
        return []

# %%


def fluid_structure(fluid):
    r"""
    Returns the checmical formula of fluid.

    Parameters
    ----------
    fluid : str
        Name of the fluid.

    Returns
    -------
    parts : dict
        Dictionary of the chemical base elements as keys and the number of
        atoms in a molecule as values.
    """
    parts = {}
    for element in CP.CoolProp.get_fluid_param_string(
            fluid, 'formula').split('}'):
        if element != '':
            el = element.split('_{')
            parts[el[0]] = int(el[1])

    return parts

# %%


def lamb(re, ks, d):
    r"""
    Calculates the darcy friction factor from the moody diagram.

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
    lamb : float
        Darcy friction factor lamb / 1

    Note
    ----
    **Laminar flow** (:math:`re \leq 2320`)

    .. math::
        \lambda = \frac{64}{re}

    **turbulent flow** (:math:`re > 2320`)

    *hydraulically smooth:* :math:`\frac{re \cdot k_{s}}{d} < 65`

    .. math::
        \lambda = \begin{cases}
        0.03164 \cdot re^{-0.25} & re \leq 10^5\\
        0.0032 + 0.221 \cdot re^{-0.237} & 10^5 < re < 5 \cdot 10^6\\
        solve \left(0 = 2 \cdot \log\left(re \cdot \sqrt{\lambda} \right) -0.8
        - \frac{1}{\sqrt{\lambda}}\right) & re \geq 5 \cdot 10^6 \\
        \end{cases}

    *transition zone:* :math:`65 \leq \frac{re \cdot k_{s}}{d} \leq 1300`

    .. math::
        \lambda = solve \left( 0 = 2 \cdot \log \left( \frac{2.51}{re \cdot
        \sqrt{\lambda}} + \frac{k_{s}}{d} \cdot 0.269 \right) -
        \frac{1}{\sqrt{\lambda}} \right)

    *hydraulically rough:* :math:`\frac{re \cdot k_{s}}{d} > 1300`

    .. math::
        \lambda = \frac{1}{\left( 2\cdot \log \left( \frac{3.71 \cdot d}{k_{s}}
        \right) \right)}

    Example
    -------
    >>> from tespy import hlp
    >>> ks = 5e-5
    >>> d = 0.05
    >>> re_laminar = 2000
    >>> re_turb_smooth = 20000
    >>> re_turb_trans = 70000
    >>> re_high = 1000000
    >>> d_high = 0.8
    >>> re_very_high = 6000000
    >>> d_very_high = 1
    >>> ks_low = 1e-5
    >>> ks_rough = 1e-3
    >>> hlp.lamb(re_laminar, ks, d)
    0.032
    >>> round(hlp.lamb(re_turb_smooth, ks, d), 3)
    0.027
    >>> round(hlp.lamb(re_turb_trans, ks, d), 3)
    0.023
    >>> round(hlp.lamb(re_turb_trans, ks_rough, d), 3)
    0.049
    >>> round(hlp.lamb(re_high, ks, d_high), 3)
    0.012
    >>> round(hlp.lamb(re_very_high, ks_low, d_very_high), 3)
    0.009
    """
    if re <= 2320:
        return 64 / re
    else:
        if re * ks / d < 65:
            if re <= 1e5:
                return 0.3164 * re ** (-0.25)
            elif re > 1e5 and re < 5e6:
                return 0.0032 + 0.221 * re ** (-0.237)
            else:
                l0 = 0.0001
                return newton(lamb_smooth, dlamb_smooth_dlamb, [re],
                              0, val0=l0, valmin=0.00001, valmax=0.2)

        elif re * ks / d > 1300:
            return 1 / (2 * math.log(3.71 * d / ks, 10)) ** 2

        else:
            l0 = 0.002
            return newton(lamb_trans, dlamb_trans_dlamb, [re, ks, d], 0,
                          val0=l0, valmin=0.0001, valmax=0.2)


def lamb_smooth(params, lamb):
    """
    TODO: DOCS!
    """
    re = params[0]
    return 2 * math.log(re * math.sqrt(lamb), 10) - 0.8 - 1 / math.sqrt(lamb)


def dlamb_smooth_dlamb(params, lamb):
    """
    TODO: DOCS!
    """
    return 1 / (lamb * math.log(10)) + 1 / 2 * lamb ** (-1.5)


def lamb_trans(params, lamb):
    """
    TODO: DOCS!
    """
    re, ks, d = params[0], params[1], params[2]
    return (2 * math.log(2.51 / (re * math.sqrt(lamb)) + ks / d * 0.269, 10) +
            1 / math.sqrt(lamb))


def dlamb_trans_dlamb(params, lamb):
    """
    TODO: DOCS!
    """
    d = 0.001
    return (lamb_trans(params, lamb + d) -
            lamb_trans(params, lamb - d)) / (2 * d)

# %%


def modify_path_os(path):
    """
    Modify the path according the os. Also detects weather the path
    specification is absolute or relative and adjusts the path respectively.

    Parameters
    ----------
    path : str
        Path to modify.

    Returns
    -------
    path : str
        Modified path.
    """

    if os.name == 'nt':
        # windows
        path = path.replace('/', '\\')
        if path[0] != '\\' and path[1:2] != ':' and path[0] != '.':
            # relative path
            path = '.\\' + path
    elif os.name == 'posix':
        # linux, max
        path = path.replace('\\', '/')
        if path[0] != '/' and path[0] != '.':
            # absolute path
            path = './' + path
    else:
        # unkown os
        msg = 'Unknown operating system, using posix pathing logic.'
        logging.warning(msg)

    return path
