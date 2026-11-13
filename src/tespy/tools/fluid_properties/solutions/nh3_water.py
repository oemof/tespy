# -*- coding: utf-8

"""Module for the NH3-water working pair.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/tools/fluid_properties/solutions/nh3_water.py

SPDX-License-Identifier: MIT
"""
import math
from functools import lru_cache

import CoolProp as CP
import numpy as np
from scipy.optimize import brentq

from .base import SolutionPair

M_NH3 = 17.03026e-3
M_H2O = 18.015268e-3
R = 8.314471
T_REF = 500.0
RHO_REF = 15000.0
X_LIMIT = 1e-10

# ideal gas part of IAPWS G4-01, Table 2: water a1 to a8 with the Einstein
# terms 4 to 8, ammonia a9 to a14
_A = (
    0, -7.720435, 8.649358, 3.006320, 0.012436, 0.97315, 1.279500,
    0.969560, 0.248730, -16.444285, 4.036946, -1.0, 10.69955, -1.775436,
    0.82374034
)
_THETA = (0, 0, 0, 0, 1.666, 4.578, 10.018, 11.964, 35.600)
_AMMONIA_EXPONENTS = {12: 1 / 3, 13: -1.5, 14: -1.75}


def _mole_fraction(w):
    return (w / M_NH3) / (w / M_NH3 + (1 - w) / M_H2O)


def _mass_fraction(x):
    return x * M_NH3 / (x * M_NH3 + (1 - x) * M_H2O)


def _molar_mass(x):
    return x * M_NH3 + (1 - x) * M_H2O


def _composition(x):
    return np.array([x, 1 - x])


def _clip(x):
    return min(max(x, X_LIMIT), 1 - X_LIMIT)


def _ideal(T, rho, x):
    r"""Ideal part :math:`\alpha^0` and :math:`\tau \partial \alpha^0 /
    \partial \tau` of the reduced Helmholtz energy, IAPWS G4-01 equation (4)
    to (6), :math:`\tau = T_\text{ref} / T`."""
    tau = T_REF / T
    water = _A[1] + _A[2] * tau + _A[3] * math.log(tau)
    water_tau = _A[2] * tau + _A[3]
    for i in range(4, 9):
        water += _A[i] * math.log(1 - math.exp(-_THETA[i] * tau))
        water_tau += _A[i] * _THETA[i] * tau / (math.exp(_THETA[i] * tau) - 1)
    ammonia = _A[9] + _A[10] * tau + _A[11] * math.log(tau)
    ammonia_tau = _A[10] * tau + _A[11]
    for i, exponent in _AMMONIA_EXPONENTS.items():
        ammonia += _A[i] * tau ** exponent
        ammonia_tau += _A[i] * exponent * tau ** exponent
    mixing = x * math.log(x) + (1 - x) * math.log(1 - x)
    alpha = (
        math.log(rho / RHO_REF) + (1 - x) * water + x * ammonia + mixing
    )
    return alpha, (1 - x) * water_tau + x * ammonia_tau


def _increasing_root(f, x0, step, xtol, lower_limit=-np.inf, upper_limit=np.inf):
    r"""Root of an increasing function :math:`f`, bracketed from :math:`x_0`
    in steps growing from :code:`step`.

    An evaluation that fails above a value with a negative residual counts
    as beyond the root, e.g. beyond the critical line of the mixture where no
    phase equilibrium exists.
    """
    def evaluate(x):
        try:
            return f(x)
        except ValueError:
            return None

    x0 = min(max(x0, lower_limit), upper_limit)
    value = evaluate(x0)
    step = abs(step)
    if value is None or value > 0:
        upper = x0
        lower = x0
        for _ in range(60):
            lower = max(lower - step, lower_limit)
            value = evaluate(lower)
            if value is not None and value < 0:
                break
            if value is not None:
                upper = lower
            if lower == lower_limit:
                raise ValueError("The root could not be bracketed.")
            step *= 2
        else:
            raise ValueError("The root could not be bracketed.")
    else:
        lower = x0
        upper = x0
        for _ in range(60):
            candidate = min(upper + step, upper_limit)
            value = evaluate(candidate)
            if value is None:
                # beyond the existence of the equilibrium: refine towards
                # the last valid point
                step /= 2
                if step < xtol:
                    raise ValueError("The root could not be bracketed.")
                continue
            if value > 0:
                upper = candidate
                break
            lower = candidate
            if candidate == upper_limit:
                raise ValueError("The root could not be bracketed.")
            step *= 2
        else:
            raise ValueError("The root could not be bracketed.")
    return brentq(f, lower, upper, xtol=xtol, rtol=1e-14)


class _PureSaturation:
    """Saturation of the pure components from CoolProp, extrapolated above
    the critical point, as starting values of the phase equilibria."""

    def __init__(self, fluid):
        self.AS = CP.AbstractState("HEOS", fluid)
        self.T_crit = self.AS.T_critical()
        self.p_crit = self.AS.p_critical()
        self.T_triple = self.AS.Ttriple()

    def p(self, T):
        if T >= self.T_crit:
            return self.p_crit * math.exp(5.4 * (1 - self.T_crit / T))
        self.AS.update(CP.QT_INPUTS, 0, max(T, self.T_triple))
        return self.AS.p()

    def T(self, p):
        if p >= self.p_crit:
            return self.T_crit / (1 - math.log(p / self.p_crit) / 5.4)
        self.AS.update(CP.PQ_INPUTS, p, 0)
        return self.AS.T()


class AmmoniaWater(SolutionPair):
    r"""
    NH3-water working pair after IAPWS G4-01 (Tillner-Roth and Friend 1998).

    The residual Helmholtz energy, its derivatives, the fugacity coefficients
    and the Newton solvers of the phase equilibrium come from the compiled
    kernel of :code:`teqp` (:cite:`Bell2022`, model
    :code:`AmmoniaWaterTillnerRoth`), which is an optional dependency of
    TESPy. The ideal part of the guideline, the density solver, the bubble
    and dew points, the flash at given pressure and temperature and the
    inversions for temperature from enthalpy and entropy are implemented
    here. The phase equilibria start from successive substitution on the
    fugacity coefficients and are polished by the Newton solvers of teqp.

    The mass fraction :code:`w` is the NH3 mass fraction. The vapor mass
    fraction outside of the two-phase region is 0 for the subcooled liquid
    and 1 for the superheated vapor.
    """

    name = "nh3-water"
    refrigerant = "NH3"
    absorbent = "H2O"
    key = "NH3"
    volatile_absorbent = True
    molar_mass = {"NH3": M_NH3, "H2O": M_H2O}
    T_min = 200
    T_max = 600

    def __init__(self):
        self._teqp_model = None
        self._saturation = None
        self._guess = {}

    @property
    def _model(self):
        if self._teqp_model is None:
            try:
                import teqp
            except ImportError as e:
                msg = (
                    "The NH3-water working pair requires the package teqp, "
                    "install it with 'pip install teqp'."
                )
                raise ImportError(msg) from e
            self._teqp_model = teqp.AmmoniaWaterTillnerRoth()
        return self._teqp_model

    @property
    def _pure(self):
        """Saturation of pure ammonia and pure water."""
        if self._saturation is None:
            self._saturation = (
                _PureSaturation("Ammonia"), _PureSaturation("Water")
            )
        return self._saturation

    def _pressure(self, T, rho, z):
        return rho * R * T * (1 + self._model.get_Ar01(T, rho, z))

    def _density(self, T, p, x, phase):
        z = _composition(x)
        model = self._model
        if phase == "l":
            rho = 1 / (x / 4e4 + (1 - x) / 5.5e4)
        else:
            rho = p / (R * T)
        for _ in range(100):
            Ar01 = model.get_Ar01(T, rho, z)
            Ar02 = model.get_Ar02(T, rho, z)
            derivative = R * T * (1 + 2 * Ar01 + Ar02)
            if derivative <= 0:
                rho *= 1.1 if phase == "l" else 0.9
                continue
            step = (rho * R * T * (1 + Ar01) - p) / derivative
            rho_new = rho - step
            if rho_new <= 0:
                rho_new = rho / 2
            if abs(rho_new - rho) < 1e-12 * rho:
                return rho_new
            rho = rho_new
        msg = (
            f"The NH3-water density at T={T} K, p={p} Pa and x={x} did not "
            "converge."
        )
        raise ValueError(msg)

    def _caloric(self, T, rho, x):
        """Mass specific enthalpy, entropy and specific volume."""
        z = _composition(x)
        model = self._model
        alpha0, alpha0_tau = _ideal(T, rho, x)
        Ar00 = model.get_Ar00(T, rho, z)
        Ar10 = model.get_Ar10(T, rho, z)
        Ar01 = model.get_Ar01(T, rho, z)
        M = _molar_mass(x)
        h = R * T * (1 + alpha0_tau + Ar10 + Ar01) / M
        s = R * (alpha0_tau + Ar10 - alpha0 - Ar00) / M
        return h, s, 1 / (rho * M)

    def _ln_phi(self, T, p, x, phase):
        rho = self._density(T, p, x, phase)
        z = _composition(x)
        return np.log(self._model.get_fugacity_coefficients(T, rho * z)), rho

    def _equilibrium(self, kind, T, spec):
        r"""
        Phase equilibrium at temperature *T* with the composition *spec* of
        one phase given: the liquid for :code:`kind="bubble"`, the vapor for
        :code:`kind="dew"`.

        Returns the pressure, the mole fraction of the other phase and the
        molar density vectors of the liquid and the vapor.
        """
        guess = self._guess.get(kind)
        if guess is not None:
            result = self._polish(kind, T, spec, *guess)
            if result is not None:
                return result
        result = self._polish(kind, T, spec, *self._substitution(kind, T, spec))
        if result is None:
            msg = (
                f"The NH3-water {kind} point at T={T} K and x={spec} did not "
                "converge."
            )
            raise ValueError(msg)
        return result

    def _substitution(self, kind, T, spec):
        """Successive substitution on the fugacity coefficients."""
        z = _composition(spec)
        ammonia, water = self._pure
        pn = ammonia.p(T)
        pw = water.p(T)
        if kind == "bubble":
            p = spec * pn + (1 - spec) * pw
            other = np.array([spec * pn, (1 - spec) * pw]) / p
        else:
            p = 1 / (spec / pn + (1 - spec) / pw)
            other = np.array([spec * p / pn, (1 - spec) * p / pw])
            other /= other.sum()
        other_phase = "g" if kind == "bubble" else "l"
        spec_phase = "l" if kind == "bubble" else "g"
        for _ in range(200):
            try:
                ln_phi_spec, rho_spec = self._ln_phi(T, p, spec, spec_phase)
                ln_phi_other, rho_other = self._ln_phi(
                    T, p, _clip(other[0]), other_phase
                )
            except ValueError:
                p *= 0.7 if kind == "bubble" else 1.3
                continue
            rho_liquid = rho_spec if kind == "bubble" else rho_other
            rho_vapor = rho_other if kind == "bubble" else rho_spec
            if rho_liquid < 3 * rho_vapor:
                p *= 1.3 if kind == "bubble" else 0.7
                continue
            if kind == "bubble":
                K = np.exp(ln_phi_spec - ln_phi_other)
                total = (K * z).sum()
                other_new = K * z / total
                p_new = p * total
            else:
                K = np.exp(ln_phi_other - ln_phi_spec)
                total = (z / K).sum()
                other_new = z / K / total
                p_new = p / total
            converged = (
                abs(p_new / p - 1) < 1e-6
                and abs(other_new[0] - other[0]) < 1e-7
            )
            p, other = 0.5 * (p + p_new), other_new
            if converged:
                break
        x_other = _clip(other[0])
        rho_spec = self._density(T, p, spec, spec_phase) * z
        rho_other = (
            self._density(T, p, x_other, other_phase) * _composition(x_other)
        )
        if kind == "bubble":
            return rho_spec, rho_other
        return rho_other, rho_spec

    def _polish(self, kind, T, spec, rho_liquid, rho_vapor):
        z = _composition(spec)
        if kind == "bubble":
            rho_liquid = rho_liquid.sum() * z
            first, second = rho_liquid, rho_vapor
        else:
            rho_vapor = rho_vapor.sum() * z
            first, second = rho_vapor, rho_liquid
        try:
            code, rho_first, rho_second = self._model.mix_VLE_Tx(
                T, first, second, z, 1e-12, 1e-12, 1e-12, 1e-12, 50
            )
        except (ValueError, RuntimeError):
            return None
        if kind == "bubble":
            rho_liquid, rho_vapor = rho_first, rho_second
        else:
            rho_liquid, rho_vapor = rho_second, rho_first
        if not self._valid(T, rho_liquid, rho_vapor):
            return None
        self._guess[kind] = (rho_liquid, rho_vapor)
        x_liquid = rho_liquid[0] / rho_liquid.sum()
        x_vapor = rho_vapor[0] / rho_vapor.sum()
        p = self._pressure(T, rho_liquid.sum(), _composition(x_liquid))
        other = x_vapor if kind == "bubble" else x_liquid
        return p, other, rho_liquid, rho_vapor

    def _valid(self, T, rho_liquid, rho_vapor, p=None):
        """Check a phase equilibrium by its own residuals."""
        if not (
                np.all(np.isfinite(rho_liquid)) and np.all(np.isfinite(rho_vapor))
                and np.all(rho_liquid >= 0) and np.all(rho_vapor >= 0)
                and rho_liquid.sum() > 1.5 * rho_vapor.sum()
            ):
            return False
        model = self._model
        x = rho_liquid / rho_liquid.sum()
        y = rho_vapor / rho_vapor.sum()
        p_liquid = self._pressure(T, rho_liquid.sum(), x)
        p_vapor = self._pressure(T, rho_vapor.sum(), y)
        if p is None:
            p = p_liquid
        if not (abs(p_liquid / p - 1) < 1e-8 and abs(p_vapor / p - 1) < 1e-8):
            return False
        with np.errstate(divide="ignore", invalid="ignore"):
            ln_f_liquid = np.log(
                model.get_fugacity_coefficients(T, rho_liquid) * x
            )
            ln_f_vapor = np.log(
                model.get_fugacity_coefficients(T, rho_vapor) * y
            )
        return bool(np.all(np.abs(ln_f_liquid - ln_f_vapor) < 1e-8))

    def _raoult_temperature(self, kind, p, spec):
        """Bubble or dew point temperature of an ideal solution."""
        ammonia, water = self._pure

        def residual(T):
            if kind == "bubble":
                p_ideal = spec * ammonia.p(T) + (1 - spec) * water.p(T)
            else:
                p_ideal = 1 / (spec / ammonia.p(T) + (1 - spec) / water.p(T))
            return math.log(p_ideal / p)

        return brentq(residual, 150, 1500, xtol=1e-6)

    def _saturation_temperature(self, kind, p, spec):
        r"""Bubble or dew point temperature at *p*, bracketed from the ideal
        solution estimate."""
        def residual(T):
            return math.log(self._equilibrium(kind, T, spec)[0] / p)

        return _increasing_root(
            residual, self._raoult_temperature(kind, p, spec), 5, 1e-12
        )

    @lru_cache(maxsize=1024)
    def _bubble_point(self, p, x):
        T = self._saturation_temperature("bubble", p, x)
        return T, self._equilibrium("bubble", T, x)[1]

    @lru_cache(maxsize=1024)
    def _dew_point(self, p, y):
        T = self._saturation_temperature("dew", p, y)
        return T, self._equilibrium("dew", T, y)[1]

    def p_bubble(self, T, w):
        return self._equilibrium("bubble", T, _clip(_mole_fraction(w)))[0]

    def p_dew(self, T, w):
        return self._equilibrium("dew", T, _clip(_mole_fraction(w)))[0]

    def T_bubble(self, p, w):
        return self._bubble_point(p, _clip(_mole_fraction(w)))[0]

    def T_dew(self, p, w):
        return self._dew_point(p, _clip(_mole_fraction(w)))[0]

    def y_eq(self, p, w):
        return _mass_fraction(self._bubble_point(p, _clip(_mole_fraction(w)))[1])

    def w_eq(self, p, y):
        return _mass_fraction(self._dew_point(p, _clip(_mole_fraction(y)))[1])

    @lru_cache(maxsize=1024)
    def w_bubble(self, p, T):
        def residual(x):
            return math.log(self._equilibrium("bubble", T, x)[0] / p)

        ammonia, water = self._pure
        x0 = _clip((p - water.p(T)) / (ammonia.p(T) - water.p(T)))
        return _mass_fraction(
            _increasing_root(residual, x0, 0.05, 1e-14, X_LIMIT, 1 - X_LIMIT)
        )

    def _flash(self, p, T, x):
        """Liquid and vapor mole fractions and molar density vectors of the
        two-phase state at *p* and *T* of the overall mole fraction *x*."""
        guesses = [self._guess.get("flash")]
        for kind in ("bubble", "dew"):
            try:
                guesses.append(self._equilibrium(kind, T, x)[2:])
            except ValueError:
                pass
        for guess in guesses:
            if guess is None:
                continue
            try:
                result = self._model.mix_VLE_Tp(T, p, *guess)
            except (ValueError, RuntimeError):
                continue
            rho_liquid, rho_vapor = result.rhovecL, result.rhovecV
            if self._valid(T, rho_liquid, rho_vapor, p):
                self._guess["flash"] = (rho_liquid, rho_vapor)
                return (
                    rho_liquid[0] / rho_liquid.sum(),
                    rho_vapor[0] / rho_vapor.sum(),
                    rho_liquid, rho_vapor
                )

        # the liquid of the two-phase state is poorer in ammonia than the
        # overall composition
        def residual(x_liquid):
            return math.log(self._equilibrium("bubble", T, x_liquid)[0] / p)

        x_liquid = _increasing_root(residual, x, -x / 2, 1e-14, X_LIMIT, x)
        _, y, rho_liquid, rho_vapor = self._equilibrium("bubble", T, x_liquid)
        self._guess["flash"] = (rho_liquid, rho_vapor)
        return x_liquid, y, rho_liquid, rho_vapor

    def _single_phase(self, p, T, x, phase):
        return self._caloric(T, self._density(T, p, x, phase), x)

    def _boundaries(self, p, x):
        try:
            return self._bubble_point(p, x)[0], self._dew_point(p, x)[0]
        except ValueError:
            return None, None

    @lru_cache(maxsize=4096)
    def _state(self, p, T, w):
        """Enthalpy, entropy, specific volume and vapor mass fraction."""
        x = _clip(_mole_fraction(w))
        T_bubble, T_dew = self._boundaries(p, x)
        if T_bubble is None:
            return (*self._single_phase(p, T, x, "l"), 0)
        if T <= T_bubble:
            return (*self._single_phase(p, T, x, "l"), 0)
        if T >= T_dew:
            return (*self._single_phase(p, T, x, "g"), 1)
        x_liquid, x_vapor, rho_liquid, rho_vapor = self._flash(p, T, x)
        liquid = self._caloric(T, rho_liquid.sum(), x_liquid)
        vapor = self._caloric(T, rho_vapor.sum(), x_vapor)
        w_liquid = _mass_fraction(x_liquid)
        w_vapor = _mass_fraction(x_vapor)
        Q = (w - w_liquid) / (w_vapor - w_liquid)
        return (*((1 - Q) * l + Q * v for l, v in zip(liquid, vapor)), Q)

    def h_pT(self, p, T, w):
        return self._state(p, T, w)[0]

    def s_pT(self, p, T, w):
        return self._state(p, T, w)[1]

    def v_pT(self, p, T, w):
        return self._state(p, T, w)[2]

    @lru_cache(maxsize=4096)
    def _invert(self, p, value, w, index):
        x = _clip(_mole_fraction(w))
        T_bubble, T_dew = self._boundaries(p, x)
        if T_bubble is None:
            return self._invert_single(p, value, x, index, "l", 230, 700)
        value_bubble = self._single_phase(p, T_bubble, x, "l")[index]
        if value <= value_bubble:
            return self._invert_single(
                p, value, x, index, "l", self.T_min, T_bubble
            )
        value_dew = self._single_phase(p, T_dew, x, "g")[index]
        if value >= value_dew:
            return self._invert_single(
                p, value, x, index, "g", T_dew, 2 * self.T_max
            )

        def residual(T):
            return self._state(p, T, w)[index] - value

        return brentq(residual, T_bubble, T_dew, xtol=1e-12, rtol=1e-14)

    def _invert_single(self, p, value, x, index, phase, lower, upper):
        def residual(T):
            return self._single_phase(p, T, x, phase)[index] - value

        d = 1e-3
        for _ in range(60):
            try:
                residual_lower = residual(lower)
                break
            except ValueError:
                # no density root of the phase, e.g. far below the bubble
                # point of a water-rich liquid
                lower = 0.5 * (lower + upper)
        if residual_lower > 0:
            # linear continuation outside of the temperature range, reached
            # only by solver trial points
            slope = (residual(lower + d) - residual_lower) / d
            return lower - residual_lower / slope
        residual_upper = residual(upper)
        if residual_upper < 0:
            slope = (residual_upper - residual(upper - d)) / d
            return upper - residual_upper / slope
        return brentq(residual, lower, upper, xtol=1e-12, rtol=1e-14)

    def T_ph(self, p, h, w, T0=None):
        return self._invert(p, h, w, 0)

    def T_ps(self, p, s, w, T0=None):
        return self._invert(p, s, w, 1)

    def Q_ph(self, p, h, w):
        return self._state(p, self.T_ph(p, h, w), w)[3]

    def phase_ph(self, p, h, w):
        x = _clip(_mole_fraction(w))
        T_bubble, T_dew = self._boundaries(p, x)
        if T_bubble is None:
            return "l"
        if h <= self._single_phase(p, T_bubble, x, "l")[0]:
            return "l"
        if h >= self._single_phase(p, T_dew, x, "g")[0]:
            return "g"
        return "tp"

    def h_pQ(self, p, Q, w):
        x = _clip(_mole_fraction(w))
        T_bubble, T_dew = self._boundaries(p, x)
        if T_bubble is None:
            msg = f"No two-phase NH3-water state at p={p} Pa and w={w}."
            raise ValueError(msg)
        if Q <= 0:
            return self._single_phase(p, T_bubble, x, "l")[0]
        if Q >= 1:
            return self._single_phase(p, T_dew, x, "g")[0]

        def residual(T):
            return self._state(p, T, w)[3] - Q

        T = brentq(residual, T_bubble, T_dew, xtol=1e-12, rtol=1e-14)
        return self._state(p, T, w)[0]
