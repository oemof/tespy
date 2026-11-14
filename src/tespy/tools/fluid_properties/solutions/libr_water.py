# -*- coding: utf-8

"""Module for the LiBr-water working pair.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/tools/fluid_properties/solutions/libr_water.py

SPDX-License-Identifier: MIT
"""
import math
from functools import lru_cache

import CoolProp as CP
from scipy.optimize import brentq

from .base import SolutionPair

M_H2O = 0.018015268
M_LIBR = 0.08685
T_C = 647.096
P_C = 22.064e6
RHO_C = 17873
H_C = 37548.5
S_C = 79.3933
T_T = 273.16
CP_T = 76.0226
T_0 = 221
X_MAX = 0.4

# coefficients of Patek and Klomfar (2006), equations (1) to (5) as tuples of
# (m, n, t, a)
_P_COEFFICIENTS = tuple(zip(
    (3, 4, 4, 8, 1, 1, 4, 6),
    (0, 5, 6, 3, 0, 2, 6, 0),
    (0, 0, 0, 0, 1, 1, 1, 1),
    (
        -2.41303e2, 1.91750e7, -1.75521e8, 3.25430e7, 3.92571e2, -2.12626e3,
        1.85127e8, 1.91216e3
    ),
))
_RHO_COEFFICIENTS = tuple(zip((1, 1), (0, 6), (1.746, 4.709)))
_H_COEFFICIENTS = tuple(zip(
    (
        1, 1, 2, 3, 6, 1, 3, 5, 4, 5, 5, 6, 6, 1, 2, 2, 2, 5, 6, 7, 1, 1, 2,
        2, 2, 3, 1, 1, 1, 1
    ),
    (
        0, 1, 6, 6, 2, 0, 0, 4, 0, 4, 5, 5, 6, 0, 3, 5, 7, 0, 3, 1, 0, 4, 2,
        6, 7, 0, 0, 1, 2, 3
    ),
    (
        0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4,
        4, 4, 4, 5, 5, 5, 5
    ),
    (
        2.27431, -7.99511, 3.85239e2, -1.63940e4, -4.22562e2, 1.13314e-1,
        -8.33474, -1.73833e4, 6.49763, 3.24552e3, -1.34643e4, 3.99322e4,
        -2.58877e5, -1.93046e-3, 2.80616, -4.04479e1, 1.45342e2, -2.74873,
        -4.49743e2, -1.21794e1, -5.83739e-3, 2.33910e-1, 3.41888e-1,
        8.85259, -1.78731e1, 7.35179e-2, -1.79430e-4, 1.84261e-3,
        -6.24282e-3, 6.84765e-3
    ),
))
_S_COEFFICIENTS = tuple(zip(
    (
        1, 1, 2, 3, 6, 1, 3, 5, 1, 2, 2, 4, 5, 5, 6, 6, 1, 3, 5, 7, 1, 1, 1,
        2, 3, 1, 1, 1, 1
    ),
    (
        0, 1, 6, 6, 2, 0, 0, 4, 0, 0, 4, 0, 4, 5, 2, 5, 0, 4, 0, 1, 0, 2, 4,
        7, 1, 0, 1, 2, 3
    ),
    (
        0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4,
        4, 4, 5, 5, 5, 5
    ),
    (
        1.53091, -4.52564, 6.98302e+2, -2.16664e+4, -1.47533e+3, 8.47012e-2,
        -6.59523, -2.95331e+4, 9.56314e-3, -1.88679e-1, 9.31752, 5.78104,
        1.38931e+4, -1.71762e+4, 4.15108e+2, -5.55647e+4, -4.23409e-3,
        3.05242e+1, -1.67620, 1.48283e+1, 3.03055e-3, -4.01810e-2,
        1.49252e-1, 2.59240, -1.77421e-1, -6.99650e-5, 6.05007e-4,
        -1.65228e-3, 1.22966e-3
    ),
))
_PS_H2O = (-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502)
_PS_H2O_EXPONENTS = (1, 1.5, 3, 3.5, 4, 7.5)
_RHO_H2O = (1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352, -6.7469445e5)
_RHO_H2O_EXPONENTS = (1 / 3, 2 / 3, 5 / 3, 16 / 3, 43 / 3, 110 / 3)
_H_H2O = (-4.37196e-1, 3.03440e-1, -1.29582, -1.76410e-1)
_H_H2O_EXPONENTS = (1 / 3, 2 / 3, 5 / 6, 21 / 6)
_S_H2O = (-3.34112e-1, -8.47987e-1, -9.11980e-1, -1.64046)
_S_H2O_EXPONENTS = (1 / 3, 1, 8 / 3, 8)


def _mole_fraction(w):
    return (w / M_LIBR) / (w / M_LIBR + (1 - w) / M_H2O)


def _mass_fraction(x):
    return x * M_LIBR / (x * M_LIBR + (1 - x) * M_H2O)


def _molar_mass(x):
    return x * M_LIBR + (1 - x) * M_H2O


def _reduced_distance_to_critical(T):
    # the correlations of saturated water end at the critical point, the
    # clamp keeps the solver iterations beyond it real valued
    return max(1 - T / T_C, 1e-12)


def _ln_ps_h2o(T):
    tau = _reduced_distance_to_critical(T)
    return math.log(P_C) + T_C / T * sum(
        a * tau ** e for a, e in zip(_PS_H2O, _PS_H2O_EXPONENTS)
    )


def _ps_h2o(T):
    return math.exp(_ln_ps_h2o(T))


def _T_sat_h2o(p):
    ln_p = math.log(min(p, P_C * (1 - 1e-9)))
    T = 1 / (1 / 373.124 - math.log(p / 101325) / 4900)
    T = min(max(T, 200), T_C - 1)
    for _ in range(50):
        tau = _reduced_distance_to_critical(T)
        f = sum(a * tau ** e for a, e in zip(_PS_H2O, _PS_H2O_EXPONENTS))
        df = sum(
            a * e * tau ** (e - 1)
            for a, e in zip(_PS_H2O, _PS_H2O_EXPONENTS)
        )
        residual = math.log(P_C) + T_C / T * f - ln_p
        derivative = -T_C / T ** 2 * f - df / T
        step = residual / derivative
        T = min(T - step, T_C * (1 - 1e-12))
        if abs(step) < 1e-11 * T:
            return T
    msg = f"The water saturation temperature at p={p} Pa did not converge."
    raise ValueError(msg)


def _rho_h2o(T):
    theta = _reduced_distance_to_critical(T)
    return RHO_C * (1 + sum(
        b * theta ** e for b, e in zip(_RHO_H2O, _RHO_H2O_EXPONENTS)
    ))


def _h_h2o(T):
    theta = _reduced_distance_to_critical(T)
    return H_C * (1 + sum(
        a * theta ** e for a, e in zip(_H_H2O, _H_H2O_EXPONENTS)
    ))


def _s_h2o(T):
    theta = _reduced_distance_to_critical(T)
    return S_C * (1 + sum(
        a * theta ** e for a, e in zip(_S_H2O, _S_H2O_EXPONENTS)
    ))


def _boiling_point_shift(x):
    r"""Return :math:`A_0(x)` and :math:`A_1(x)` of the sum in equation (1),
    which is linear in the reduced temperature:
    :math:`\sum = A_0 + A_1 T / T_c`."""
    A0 = 0
    A1 = 0
    for m, n, t, a in _P_COEFFICIENTS:
        term = a * x ** m * (X_MAX - x) ** n
        if t == 0:
            A0 += term
        else:
            A1 += term
    return A0, A1


def _excess(coefficients, T, x):
    tau = T_C / (T - T_0)
    return sum(
        a * x ** m * (X_MAX - x) ** n * tau ** t
        for m, n, t, a in coefficients
    )


def _p_bubble_molar(T, x):
    A0, A1 = _boiling_point_shift(x)
    return _ps_h2o(T - A0 - A1 * T / T_C)


def _h_saturated(T, x):
    """Molar enthalpy of the saturated solution, equation (4)."""
    return (1 - x) * _h_h2o(T) + H_C * _excess(_H_COEFFICIENTS, T, x)


def _s_saturated(T, x):
    """Molar entropy of the saturated solution, equation (5)."""
    return (1 - x) * _s_h2o(T) + S_C * _excess(_S_COEFFICIENTS, T, x)


def _rho_saturated(T, x):
    """Molar density of the saturated solution, equation (2)."""
    return (1 - x) * _rho_h2o(T) + RHO_C * sum(
        a * x ** m * (T / T_C) ** n for m, n, a in _RHO_COEFFICIENTS
    )


class _WaterVapor:
    """Superheated water vapor from CoolProp's IAPWS-95."""

    def __init__(self):
        self.AS = CP.AbstractState("HEOS", "Water")

    def state(self, p, T):
        """Return enthalpy, entropy and specific volume."""
        self.AS.specify_phase(CP.iphase_gas)
        try:
            self.AS.update(CP.PT_INPUTS, p, T)
        finally:
            self.AS.unspecify_phase()
        return self.AS.hmass(), self.AS.smass(), 1 / self.AS.rhomass()


class LiBrWater(SolutionPair):
    r"""
    LiBr-water working pair after Pátek and Klomfar (2006).

    The liquid solution follows the explicit formulation of
    :cite:`Patek2006`: vapor pressure, density, enthalpy and entropy of the
    saturated solution as function of temperature and composition, with the
    reference state of IAPWS-95 for water. The liquid is incompressible, its
    enthalpy rises by :math:`v \left(p - p_\text{bubble}\right)` above the
    bubble pressure, its entropy does not depend on pressure. The vapor is
    pure water, evaluated with IAPWS-95 at the actual pressure and
    temperature.

    Below the bubble pressure the state splits into the saturated liquid of
    the composition :math:`w_\text{l} > w` in equilibrium at the given
    pressure and temperature and superheated water vapor. The vapor mass
    fraction follows from the LiBr balance :math:`Q = 1 - w / w_\text{l}`.

    The mass fraction :code:`w` is the LiBr mass fraction. The formulation is
    valid from 273 K to 500 K over the full composition range up to the
    crystallization line; outside of that range all functions extrapolate
    smoothly to support the solver iterations.
    """

    name = "libr-water"
    refrigerant = "H2O"
    absorbent = "LiBr"
    key = "LiBr"
    volatile_absorbent = False
    molar_mass = {"H2O": M_H2O, "LiBr": M_LIBR}
    T_min = 273.15
    T_max = 500
    w_max = _mass_fraction(0.395)

    def __init__(self):
        self._vapor = _WaterVapor()

    def p_bubble(self, T, w):
        return _p_bubble_molar(T, _mole_fraction(w))

    def T_bubble(self, p, w):
        A0, A1 = _boiling_point_shift(_mole_fraction(w))
        return (_T_sat_h2o(p) + A0) / (1 - A1 / T_C)

    def w_bubble(self, p, T):
        theta = _T_sat_h2o(p)

        def residual(x):
            A0, A1 = _boiling_point_shift(x)
            return T - A0 - A1 * T / T_C - theta

        x_max = _mole_fraction(self.w_max)
        if residual(0) <= 0:
            return 0.0
        if residual(x_max) > 0:
            msg = (
                f"No LiBr-water solution with a LiBr mass fraction below "
                f"{self.w_max:.4f} is in equilibrium at p={p} Pa and T={T} K."
            )
            raise ValueError(msg)
        return _mass_fraction(brentq(residual, 0, x_max, xtol=1e-15))

    def y_eq(self, p, w):
        return 0.0

    def _liquid(self, p, T, w):
        x = _mole_fraction(w)
        M = _molar_mass(x)
        v = 1 / (M * _rho_saturated(T, x))
        dp = p - _p_bubble_molar(T, x)
        return (
            _h_saturated(T, x) / M + v * dp,
            _s_saturated(T, x) / M,
            v,
        )

    def _two_phase(self, p, T, w, w_liquid):
        Q = 1 - w / w_liquid
        liquid = self._liquid(p, T, w_liquid)
        vapor = self._vapor.state(p, T)
        return tuple((1 - Q) * l + Q * v for l, v in zip(liquid, vapor))

    @lru_cache(maxsize=1024)
    def _state_pT(self, p, T, w):
        if p >= self.p_bubble(T, w):
            return self._liquid(p, T, w)
        try:
            w_liquid = self.w_bubble(p, T)
        except ValueError:
            w_liquid = self.w_max
        return self._two_phase(p, T, w, max(w_liquid, w))

    def h_pT(self, p, T, w):
        return self._state_pT(p, T, w)[0]

    def s_pT(self, p, T, w):
        return self._state_pT(p, T, w)[1]

    def v_pT(self, p, T, w):
        return self._state_pT(p, T, w)[2]

    @lru_cache(maxsize=1024)
    def _invert(self, p, value, w, index):
        T_bubble = self.T_bubble(p, w)
        value_bubble = self._liquid(p, T_bubble, w)[index]
        if value <= value_bubble:
            return self._invert_liquid(p, value, w, index, T_bubble)

        T_max = self.T_bubble(p, self.w_max)
        if value <= self._two_phase(p, T_max, w, self.w_max)[index]:
            def residual(w_liquid):
                T = self.T_bubble(p, w_liquid)
                return self._two_phase(p, T, w, w_liquid)[index] - value

            # on the bubble line the residual at w is round-off with either
            # sign, brentq evaluates it as python float
            w = float(w)
            if residual(w) >= 0:
                return T_bubble
            w_liquid = brentq(residual, w, self.w_max, xtol=1e-14)
            return self.T_bubble(p, w_liquid)

        def residual(T):
            return self._two_phase(p, T, w, self.w_max)[index] - value

        for T_upper in (T_max + 10, T_max + 100, 1000, 2000):
            if residual(T_upper) >= 0:
                return brentq(residual, T_max, T_upper, xtol=1e-12)
        # linear continuation far beyond any meaningful state, reached only
        # by solver trial points
        d = 1e-3
        slope = (residual(T_upper) - residual(T_upper - d)) / d
        return T_upper - residual(T_upper) / slope

    def _invert_liquid(self, p, value, w, index, T_bubble):
        def residual(T):
            return self._liquid(p, T, w)[index] - value

        d = 1e-3
        T_lower = T_0 + 20
        residual_lower = residual(T_lower)
        if residual_lower >= 0:
            # linear continuation below the range of the formulation
            slope = (residual(T_lower + d) - residual_lower) / d
            return T_lower - residual_lower / slope

        lower, upper = T_lower, T_bubble
        T = upper - 10 if upper - 10 > lower else (lower + upper) / 2
        for _ in range(100):
            r = residual(T)
            if r < 0:
                lower = T
            else:
                upper = T
            slope = (residual(T + d) - residual(T - d)) / (2 * d)
            T_new = T - r / slope if slope > 0 else (lower + upper) / 2
            if not lower <= T_new <= upper:
                T_new = (lower + upper) / 2
            if abs(T_new - T) < 1e-10:
                return T_new
            T = T_new
        msg = (
            f"The liquid LiBr-water temperature at p={p} Pa and w={w} did not "
            "converge."
        )
        raise ValueError(msg)

    def T_ph(self, p, h, w, T0=None):
        return self._invert(p, h, w, 0)

    def T_ps(self, p, s, w, T0=None):
        return self._invert(p, s, w, 1)

    def _liquid_fraction_ph(self, p, h, w):
        T = self.T_ph(p, h, w)
        if p >= self.p_bubble(T, w):
            return w
        try:
            return max(self.w_bubble(p, T), w)
        except ValueError:
            return self.w_max

    def Q_ph(self, p, h, w):
        return 1 - w / self._liquid_fraction_ph(p, h, w)

    def phase_ph(self, p, h, w):
        if h <= self.h_pQ(p, 0, w):
            return "l"
        return "tp"

    def p_TQ(self, T, Q, w):
        if not 0 <= Q < 1:
            msg = (
                "The vapor mass fraction of a LiBr-water solution must be "
                f"within [0, 1), {Q} was given."
            )
            raise ValueError(msg)
        return self.p_bubble(T, w / (1 - Q))

    def h_pQ(self, p, Q, w):
        if not 0 <= Q < 1:
            msg = (
                "The vapor mass fraction of a LiBr-water solution must be "
                f"within [0, 1), {Q} was given."
            )
            raise ValueError(msg)
        w_liquid = w / (1 - Q)
        T = self.T_bubble(p, w_liquid)
        return self._two_phase(p, T, w, w_liquid)[0]
