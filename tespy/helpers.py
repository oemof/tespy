"""
.. module:: helpers
    :platforms: all
    :synopsis: helpers for frequently used functionalities

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>
"""

import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI as CPPSI

import math
import numpy as np
import pandas as pd
import sys

from scipy.optimize import fsolve

global err
err = 1e-6
global molar_masses
molar_masses = {}
global gas_constants
gas_constants = {}


class memorise:

    def __init__(self, num_fl):
        memorise.T_ph = np.empty((0, num_fl + 3), float)
        memorise.T_ps = np.empty((0, num_fl + 4), float)
        memorise.v_ph = np.empty((0, num_fl + 3), float)
        memorise.visc_ph = np.empty((0, num_fl + 3), float)
        memorise.s_ph = np.empty((0, num_fl + 3), float)


class MyNetworkError(Exception):
    pass


class MyConnectionError(Exception):
    pass


class MyComponentError(Exception):
    pass


class MyConvergenceError(Exception):
    pass


def query_yes_no(question, default='yes'):
    """
    in prompt query

    :param question: question to ask in prompt
    :type question: str
    :param default: default answer
    :type default: str
    :returns: bool
    """
    valid = {'yes': True,
             'y': True,
             'ye': True,
             'no': False,
             'n': False}
    if default is None:
        prompt = '[y / n]'
    elif default == 'yes':
        prompt = '[Y / n]'
    elif default == 'no':
        prompt = '[y / N]'

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write('Please respond with \'yes\' or \'no\' '
                             '(or \'y\' or \'n\').\n')


def newton(func, deriv, flow, k):
    r"""
    find zero crossings of function func with 1-D newton algorithm,
    required for reverse functions of fluid mixtures

    :param func: function to find zero crossing in
    :type func: function
    :param deriv: derivative of the function
    :type deriv: function
    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :param k: target value for function func
    :type k: numeric
    :returns: val (float) - val, so that func(flow, val) = k

    .. math::

        x_{i+1} = x_{i} - \frac{f(x_{i})}{\frac{df}{dx}(x_{i})}\\
        f(x_{n}) \leq \epsilon, \; n < 10\\
        n: \text{number of iterations}

    """
    res = 1
    val = 300
    i = 0
    while abs(res) >= err:
        res = k - func(flow, val)
        val += res / deriv(flow, val)

        if val < 69:
            val = 69
        i += 1

        if i > 10:
            raise ValueError('No value found.')

    return val


def T_mix_ph(flow):
    r"""
    calculates the temperature from pressure and enthalpy,
    uses CoolProp reverse functions for pure fluids, newton for mixtures
        - check if property has already been memorised
        - calculate property otherwise

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :returns: T (float) - temperature in K

    **fluid mixtures**

    .. math::

        T_{mix}\left(p,h\right) = T_{i}\left(pp_{i},h_{i}\right)\;
        \forall i \in \text{fluid components}\\

        h_{i} = h \left(pp_{i}, T_{mix} \right)\\
        pp: \text{partial pressure}

    """
    a = memorise.T_ph[:, 0:-1]
    b = np.array([flow[1], flow[2]] + list(flow[3].values()))
    ix = np.where(np.all(abs(a - b) <= err**2, axis=1))[0]
    if ix.size == 1:
        return memorise.T_ph[ix, -1][0]
    else:
        if num_fluids(flow[3]) > 1:
            val = newton(h_mix_pT, dh_mix_pdT, flow, flow[2])
            new = np.array([[flow[1], flow[2]] + list(flow[3].values()) +
                            [val]])
            memorise.T_ph = np.append(memorise.T_ph, new, axis=0)
            return val
        else:
            for fluid, x in flow[3].items():
                if x > err:
                    val = CPPSI('T', 'H', flow[2], 'P', flow[1], fluid)
                    return val


def T_mix_ps(flow, s):
    r"""
    calculates the temperature from pressure and entropy,
    uses CoolProp reverse functions for pure fluids, newton for mixtures

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :param s: entropy in J / (kg * K)
    :type s: numeric
    :returns: T (float) - temperature in K

    **fluid mixtures**

    .. math::

        T_{mix}\left(p,s\right) = T_{i}\left(pp_{i},s_{i}\right)\;
        \forall i \in \text{fluid components}\\

        s_{i} = s \left(pp_{i}, T_{mix} \right)\\
        pp: \text{partial pressure}
    """
    a = memorise.T_ps[:, 0:-1]
    b = np.array([flow[1], flow[2]] + list(flow[3].values()) + [s])
    ix = np.where(np.all(abs(a - b) <= err**2, axis=1))[0]
    if ix.size == 1:
        return memorise.T_ps[ix, -1][0]
    else:
        if num_fluids(flow[3]) > 1:
            val = newton(s_mix_pT, ds_mix_pdT, flow, s)
            new = np.array([[flow[1], flow[2]] + list(flow[3].values()) +
                            [s, val]])
            memorise.T_ps = np.append(memorise.T_ps, new, axis=0)
            return val
        else:
            for fluid, x in flow[3].items():
                if x > err:
                    val = CPPSI('T', 'S', s, 'P', flow[1], fluid)
                    new = np.array([[flow[1], flow[2]] +
                                    list(flow[3].values()) + [s, val]])
                    memorise.T_ps = np.append(memorise.T_ps, new, axis=0)
                    return val


def dT_mix_dph(flow):
    r"""
    calculates partial derivate of temperature to pressure at
    constant enthalpy and fluid composition

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :returns: dT / dp (float) - derivative in K / Pa

    .. math::

        \frac{\partial T_{mix}}{\partial p} = \frac{T_{mix}(p+d,h)-
        T_{mix}(p-d,h)}{2 \cdot d}
    """
    d = 1
    u = flow.copy()
    l = flow.copy()
    u[1] += d
    l[1] -= d
    return (T_mix_ph(u) - T_mix_ph(l)) / (2 * d)


def dT_mix_pdh(flow):
    r"""
    method to calculate partial derivate of temperature to enthalpy at
    constant pressure and fluid composition

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :returns: dT / dh (float) - derivative in (K * kg) / J

    .. math::

        \frac{\partial T_{mix}}{\partial h} = \frac{T_{mix}(p,h+d)-
        T_{mix}(p,h-d)}{2 \cdot d}
    """
    d = 1
    u = flow.copy()
    l = flow.copy()
    u[2] += d
    l[2] -= d
    return (T_mix_ph(u) - T_mix_ph(l)) / (2 * d)


def dT_mix_ph_dfluid(flow):
    r"""
    calculates partial derivates of temperature to fluid composition at
    constant pressure and enthalpy

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :returns: dT / dfluid (np.array of floats) - derivatives in K

    .. math::

        \frac{\partial T_{mix}}{\partial fluid_{i}} =
        \frac{T_{mix}(p,h,fluid_{i}+d)-
        T_{mix}(p,h,fluid_{i}-d)}{2 \cdot d}
    """
    d = 1e-5
    u = flow.copy()
    l = flow.copy()
    vec_deriv = []
    for fluid, x in flow[3].items():
        if x > err:
            u[3][fluid] += d
            l[3][fluid] -= d
            vec_deriv += [(T_mix_ph(u) - T_mix_ph(l)) / (2 * d)]
            u[3][fluid] -= d
            l[3][fluid] += d
        else:
            vec_deriv += [0]

    return np.asarray(vec_deriv)


def h_mix_pT(flow, T):
    r"""
    calculates enthalpy from pressure and temperature

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :param T: temperature in K
    :type T: numeric
    :returns: h (float) - enthalpy in J / kg

    .. math::
        h_{mix}(p,T)=\sum_{i} h(pp_{i},T,fluid_{i})\;
        \forall i \in \text{fluid components}\\
        pp: \text{partial pressure}

    """

    n = molar_massflow(flow[3])

    h = 0
    for fluid, x in flow[3].items():
        if x > err:
            ni = x / molar_masses[fluid]
            h += CPPSI('H', 'P', flow[1] * ni / n, 'T', T, fluid) * x

    return h


def dh_mix_pdT(flow, T):
    r"""
    calculates partial derivate of enthalpy to temperature at constant pressure

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :param T: temperature in K
    :type T: numeric
    :returns: dh / dT (float) - derivative in J / (kg * K)

    .. math::

        \frac{\partial h_{mix}}{\partial T} =
        \frac{h_{mix}(p,T+d)-h_{mix}(p,T-d)}{2 \cdot d}
    """
    n = molar_massflow(flow[3])
    d = 2
    h_u = 0
    h_l = 0
    for fluid, x in flow[3].items():
        if x > err:
            pp = flow[1] * x / (molar_masses[fluid] * n)
            h_u += CPPSI('H', 'P', pp, 'T', T + d, fluid) * x
            h_l += CPPSI('H', 'P', pp, 'T', T - d, fluid) * x

    return (h_u - h_l) / (2 * d)


def h_mix_pQ(flow, Q):
    """
    calculates enthalpy from pressure and quality

    .. note::

       This function works for pure fluids only!

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :param Q: fraction of vapour mass to total mass in 1
    :type Q: numeric
    :returns: h (float) - enthalpy in J / kg
    """
    n = molar_massflow(flow[3])

    h = 0
    for fluid, x in flow[3].items():
        if x > err:
            pp = flow[1] * x / (molar_masses[fluid] * n)
            pcrit = CPPSI('Pcrit', fluid)
            while pp > pcrit:
                flow[1] = flow[1] * 0.95

            h += CPPSI('H', 'P', pp, 'Q', Q, fluid) * x

    return h


def dh_mix_dpQ(flow, Q):
    r"""
    calculates partial derivative of enthalpy to pressure at constant quality

    .. note::

       This function works for pure fluids only!

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :param Q: fraction of vapour mass to total mass in 1
    :type Q: numeric
    :returns: dh / dp (float) - derivative in J / (kg * Pa)

    .. math::

        \frac{\partial h_{mix}}{\partial p} =
        \frac{h_{mix}(p+d,Q)-h_{mix}(p-d,Q)}{2 \cdot d}\\
        Q: \text{vapour mass fraction}
    """
    d = 1
    u = flow.copy()
    l = flow.copy()
    u[1] += d
    l[1] -= d
    return (h_mix_pQ(u, Q) - h_mix_pQ(l, Q)) / (2 * d)


def v_mix_ph(flow):
    r"""
    calculates specific volume from pressure and enthalpy
    uses CoolProp reverse functions for pure fluids, newton for mixtures

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :returns: v (float) - specific volume in kg / m :sup:`3`

    **fluid mixtures**

    .. math::

        v_{mix}\left(p,h\right) = v\left(p,T_{mix}(p,h)\right)
    """
    a = memorise.v_ph[:, 0:-1]
    b = np.array([flow[1], flow[2]] + list(flow[3].values()))
    ix = np.where(np.all(abs(a - b) <= err**2, axis=1))[0]
    if ix.size == 1:
        return memorise.v_ph[ix, -1][0]
    else:
        if num_fluids(flow[3]) > 1:
            val = v_mix_pT(flow, T_mix_ph(flow))
            new = np.array([[flow[1], flow[2]] + list(flow[3].values()) +
                            [val]])
            memorise.v_ph = np.append(memorise.v_ph, new, axis=0)
            return val
        else:
            for fluid, x in flow[3].items():
                if x > err:
                    val = 1 / CPPSI('D', 'P', flow[1], 'H', flow[2], fluid)
                    new = np.array([[flow[1], flow[2]] +
                                    list(flow[3].values()) + [val]])
                    memorise.v_ph = np.append(memorise.v_ph, new, axis=0)
                    return val


def v_mix_pT(flow, T):
    r"""
    calculates specific volume from pressure and temperature

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :param T: temperature in K
    :type T: numeric
    :returns: v (float) - specific volume in kg / m :sup:`3`

    .. math::
        v_{mix}(p,T)=\sum_{i} v(pp_{i},T,fluid_{i})\;
        \forall i \in \text{fluid components}\\
        pp: \text{partial pressure}
    """
    n = molar_massflow(flow[3])

    d = 0
    for fluid, x in flow[3].items():
        if x > err:
            pp = flow[1] * x / (molar_masses[fluid] * n)
            d += CPPSI('D', 'P', pp, 'T', T, fluid) * x

    return 1 / d


def visc_mix_ph(flow):
    r"""
    calculates dynamic viscosity from pressure and enthalpy,
    uses CoolProp reverse functions for pure fluids, newton for mixtures

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :returns: v (float) - specific volume in Pa s

    **fluid mixtures**

    .. math::

        \eta_{mix}\left(p,h\right) = \eta\left(p,T_{mix}(p,h)\right)
    """
    if num_fluids(flow[3]) > 1:
        return visc_mix_pT(flow, T_mix_ph(flow))
    else:
        for fluid, x in flow[3].items():
            if x > err:
                return CPPSI('V', 'P', flow[1], 'H', flow[2], fluid)


def visc_mix_pT(flow, T):
    r"""
    calculates dynamic viscosity from pressure and temperature

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :param T: temperature in K
    :type T: numeric
    :returns: v (float) - specific volume in kg / m :sup:`3`

    .. math::
        \eta_{mix}(p,T)=\frac{\sum_{i} \left( \eta(p,T,fluid_{i}) \cdot y_{i}
        \cdot M_{i} \right)}
        {\sum_{i} \left(y_{i} \cdot M_{i} \right)}\;
        \forall i \in \text{fluid components}\\
        y: \text{volume fraction}\\
        M: \text{molar mass}
    """
    n = molar_massflow(flow[3])

    a = 0
    b = 0
    for fluid, x in flow[3].items():
        if x > err:
            bi = x * math.sqrt(molar_masses[fluid]) / (molar_masses[fluid] * n)
            b += bi
            a += bi * CPPSI('V', 'P', flow[1], 'T', T, fluid)

    return a / b


def s_mix_ph(flow):
    r"""
    calculates entropy from pressure and enthalpy
    uses CoolProp reverse functions for pure fluids, newton for mixtures

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :returns: s (float) - entropy in J / (kg * K)

    **fluid mixtures**

    .. math::

        s_{mix}\left(p,h\right) = s\left(p,T_{mix}(p,h)\right)
    """
    a = memorise.s_ph[:, 0:-1]
    b = np.array([flow[1], flow[2]] + list(flow[3].values()))
    ix = np.where(np.all(abs(a - b) <= err**2, axis=1))[0]
    if ix.size == 1:
        return memorise.s_ph[ix, -1][0]
    else:
        if num_fluids(flow[3]) > 1:
            val = s_mix_pT(flow, T_mix_ph(flow))
            new = np.array([[flow[1], flow[2]] + list(flow[3].values()) +
                            [val]])
            memorise.s_ph = np.append(memorise.s_ph, new, axis=0)
            return val
        else:
            for fluid, x in flow[3].items():
                if x > err:
                    val = CPPSI('S', 'P', flow[1], 'H', flow[2], fluid)
                    new = np.array([[flow[1], flow[2]] +
                                    list(flow[3].values()) + [val]])
                    memorise.s_ph = np.append(memorise.s_ph, new, axis=0)
                    return val


def s_mix_pT(flow, T):
    r"""
    calculates entropy from pressure and temperature
    uses CoolProp reverse functions for pure fluids, newton for mixtures

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :param T: temperature in K
    :type T: numeric
    :returns: s (float) - entropy in J / (kg * K)

    .. math::
        s_{mix}(p,T)=\sum_{i} x_{i} \cdot s(pp_{i},T,fluid_{i})-
        \sum_{i} x_{i} \cdot R_{i} \cdot \ln \frac{pp_{i}}{p}\;
        \forall i \in \text{fluid components}\\
        pp: \text{partial pressure}\\
        R: \text{gas constant}
    """
    n = molar_massflow(flow[3])

    s = 0
    for fluid, x in flow[3].items():
        if x > err:
            pp = flow[1] * x / (molar_masses[fluid] * n)
            s += CPPSI('S', 'P', pp, 'T', T, fluid) * x
            s -= (x * gas_constants[fluid] / molar_masses[fluid] *
                  math.log(pp / flow[1]))

    return s


def ds_mix_pdT(flow, T):
    r"""
    calculates partial derivate of entropy to temperature at constant pressure

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :param T: temperature in K
    :type T: numeric
    :returns: ds / dT (float) - derivative in J / (kg * K :sup:`2`)

    .. math::

        \frac{\partial s_{mix}}{\partial T} =
        \frac{s_{mix}(p,T+d)-s_{mix}(p,T-d)}{2 \cdot d}
    """
    n = molar_massflow(flow[3])
    d = 2
    s_u = 0
    s_l = 0
    for fluid, x in flow[3].items():
        if x > err:
            pp = flow[1] * x / (molar_masses[fluid] * n)
            s_u += CPPSI('S', 'P', pp, 'T', T + d, fluid) * x
            s_l += CPPSI('S', 'P', pp, 'T', T - d, fluid) * x

    return (s_u - s_l) / (2 * d)


def molar_massflow(flow):
    r"""
    calculates molar massflow

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :returns: mm (float) - molar massflow in mol / s

    .. math::
        mm = \sum_{i} \left( \frac{x_{i}}{M_{i}} \right)
    """
    mm = 0
    for fluid, x in flow.items():
        if x > err:
            try:
                mm += x / molar_masses[fluid]
            except:
                mm += x / CPPSI('molar_mass', fluid)

    return mm


def num_fluids(fluids):
    r"""
    calculates number of fluids in fluid vector

    :param fluids: fluid vector {fluid: mass fraction}
    :type fluids: dict
    :returns: n (int) - number of fluids in fluid vector in 1

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


def fluid_structure(fluid):
    """
    gets the chemical formular of a fluid

    :param fluid: alias of the fluid
    :type fluid: str
    :returns: parts (dict) - returns the elements of the fluid {'element': n}
    """
    parts = {}
    for element in CP.get_fluid_param_string(fluid, 'formula').split('}'):
        if element != '':
            el = element.split('_{')
            parts[el[0]] = int(el[1])

    return parts

def lamb(re, ks, d):
    r"""
    calculates darcy friction factor

    :param re: reynolds number in 1
    :type re: numeric
    :param ks: roughness in m
    :type ks: numeric
    :param d: pipe diameter in m
    :type d: numeric
    :returns: lambda (float) - darcy friction factor in 1

    **laminar flow:** :math:`re \leq 2320`

    .. math::
        \lambda = \frac{64}{re}

    **turbulent flow:** :math:`re > 2320`

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
        \lambda = \frac{1}{\left( 2\cdot \log \left( \frac{3.71 \cdot d}{k_{s}} \right) \right)}


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
                func = lambda l: (2 * math.log(re * math.sqrt(l), 10) -
                                  0.8 - 1 / math.sqrt(l))
                return fsolve(func, l0)
        elif re * ks / d > 1300:
            return 1 / (2 * math.log(3.71 * d / ks, 10)) ** 2
        else:
            l0 = 0.002
            func = lambda l: (2 * math.log(2.51 / (re * math.sqrt(l)) +
                              ks / d * 0.269, 10) + 1 / math.sqrt(l))
            return fsolve(func, l0)
