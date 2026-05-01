# -*- coding: utf-8

"""Module for fluid property mixture routines.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/tools/fluid_properties/mixtures.py

SPDX-License-Identifier: MIT
"""

import math

from CoolProp.CoolProp import HAPropsSI

from tespy.tools.global_vars import FLUID_ALIASES
from tespy.tools.global_vars import gas_constants
from tespy.tools.logger import logger

from .helpers import _is_larger_than_precision
from .helpers import calc_molar_mass_mixture
from .helpers import get_molar_fractions


def h_mix_pT_ideal(p=None, T=None, fluid_data=None, **kwargs):
    molar_fractions = get_molar_fractions(fluid_data)

    h = 0
    for fluid, data in fluid_data.items():

        if _is_larger_than_precision(data["mass_fraction"]):
            pp = p * molar_fractions[fluid]
            h += data["wrapper"].h_pT(pp, T) * data["mass_fraction"]

    return h


def h_mix_pT_ideal_cond(p=None, T=None, fluid_data=None, **kwargs):

    water_alias = _get_fluid_alias("H2O", fluid_data)
    if water_alias:
        water_alias = next(iter(water_alias))
        mass_fractions_gas, molar_fraction_gas, mass_liquid, _, p_sat, pp_water = cond_check(p, T, fluid_data, water_alias)
        # at saturation liquid mass may be zero, but we cannot calculate water properties with pT
        if not _is_larger_than_precision(mass_liquid) and abs(pp_water - p_sat) / p_sat > 1e-6:
            return h_mix_pT_ideal(p, T, fluid_data, **kwargs)
        h = 0
        for fluid, data in fluid_data.items():
            if _is_larger_than_precision(data["mass_fraction"]):
                if fluid == water_alias:
                    # at saturation liquid mass may be zero, but we cannot calculate water properties with pT
                    if mass_liquid > 0:
                        h += fluid_data[water_alias]["wrapper"].h_QT(0, T) * mass_liquid
                    h += fluid_data[water_alias]["wrapper"].h_QT(1, T) * mass_fractions_gas[fluid] * (1 - mass_liquid)
                else:
                    pp = p * molar_fraction_gas[fluid]
                    h += data["wrapper"].h_pT(pp, T) * mass_fractions_gas[fluid] * (1 - mass_liquid)
        return h
    else:
        return h_mix_pT_ideal(p, T, fluid_data, **kwargs)


def h_mix_pT_forced_gas(p, T, fluid_data, **kwargs):
    molar_fractions = get_molar_fractions(fluid_data)

    h = 0
    for fluid, data in fluid_data.items():

        if _is_larger_than_precision(data["mass_fraction"]):
            pp = p * molar_fractions[fluid]
            if fluid == "H2O" and pp >= data["wrapper"]._p_min:
                if T <= data["wrapper"].T_sat(pp):
                    h += data["wrapper"].h_QT(1, T) * data["mass_fraction"]
                else:
                    h += data["wrapper"].h_pT(pp, T) * data["mass_fraction"]
            else:
                h += data["wrapper"].h_pT(pp, T) * data["mass_fraction"]

    return h


def h_mix_pT_incompressible(p, T, fluid_data, **kwargs):

    h = 0
    for data in fluid_data.values():
        if _is_larger_than_precision(data["mass_fraction"]):
            h += data["wrapper"].h_pT(p, T) * data["mass_fraction"]

    return h


def w_mix_fluid_data(fluid_data):

    water_alias = _get_fluid_alias("H2O", fluid_data)
    water_alias = next(iter(water_alias))

    air_alias = _get_fluid_alias("air", fluid_data)
    air_alias = next(iter(air_alias))

    return (
        fluid_data[water_alias]["mass_fraction"]
        / fluid_data[air_alias]["mass_fraction"]
    )

def w_mix_pTrh_humidair(p, T, rh):
    return HAPropsSI("W", "P", p, "T", T, "RH", rh)  # kg water/kg dry air

def w_mix_pT_humidair(p, T, fluid_data, **kwargs):
    w_def = w_mix_fluid_data(fluid_data)
    w_max = w_mix_pTrh_humidair(p, T, 1.0)
    if w_def > w_max:
        _msg = f"Humidity ratio {w_def:.4f} exceeds maximum value of {w_max:.4f} for given p and T. Check fluid composition."
        return w_max
    return w_def

def h_mix_pT_humidair(p, T, fluid_data, **kwargs):
    w = w_mix_pT_humidair(p, T, fluid_data, **kwargs)
    return HAPropsSI("H", "P", p, "T", T, "W", w)

def w_mix_phrh_humidair(p, h, rh):
    return HAPropsSI("W", "P", p, "H", h, "RH", rh)  # kg water/kg dry air

def w_mix_ph_humidair(p, h, fluid_data, **kwargs):
    w_def = w_mix_fluid_data(fluid_data)
    w_max = w_mix_phrh_humidair(p, h, 1.0)
    if w_def > w_max:
        _msg = f"Humidity ratio {w_def:.4f} exceeds maximum value of {w_max:.4f} for given p and T. Check fluid composition."
        return w_max
    return w_def

def w_mix_psrh_humidair(p, s, rh):
    return HAPropsSI("W", "P", p, "S", s, "RH", rh)  # kg water/kg dry air

def w_mix_ps_humidair(p, s, fluid_data, **kwargs):
    w_def = w_mix_fluid_data(fluid_data)
    w_max = w_mix_psrh_humidair(p, s, 1.0)
    if w_def > w_max:
        _msg = f"Humidity ratio {w_def:.4f} exceeds maximum value of {w_max:.4f} for given p and T. Check fluid composition."
        return w_max
    return w_def

def s_mix_pT_ideal(p=None, T=None, fluid_data=None, **kwargs):
    molar_fractions = get_molar_fractions(fluid_data)

    s = 0
    for fluid, data in fluid_data.items():

        if _is_larger_than_precision(data["mass_fraction"]):
            pp = p * molar_fractions[fluid]
            s += data["wrapper"].s_pT(pp, T) * data["mass_fraction"]

    return s


def s_mix_pT_ideal_cond(p=None, T=None, fluid_data=None, **kwargs):

    water_alias = _get_fluid_alias("H2O", fluid_data)
    if water_alias:
        water_alias = next(iter(water_alias))
        mass_fractions_gas, molar_fraction_gas, mass_liquid, _, p_sat, pp_water = cond_check(p, T, fluid_data, water_alias)
        # at saturation liquid mass may be zero, but we cannot calculate water properties with pT
        if not _is_larger_than_precision(mass_liquid) and abs(pp_water - p_sat) / p_sat > 1e-6:
            return s_mix_pT_ideal(p, T, fluid_data, **kwargs)
        s = 0
        for fluid, data in fluid_data.items():
            if _is_larger_than_precision(data["mass_fraction"]):
                if fluid == water_alias:
                    if mass_liquid > 0:
                        s += fluid_data[water_alias]["wrapper"].s_QT(0, T) * mass_liquid
                    s += fluid_data[water_alias]["wrapper"].s_QT(1, T) * mass_fractions_gas[fluid] * (1 - mass_liquid)
                else:
                    pp = p * molar_fraction_gas[fluid]
                    s += data["wrapper"].s_pT(pp, T) * mass_fractions_gas[fluid] * (1 - mass_liquid)
        return s
    else:
        return s_mix_pT_ideal(p, T, fluid_data, **kwargs)


def s_mix_pT_incompressible(p=None, T=None, fluid_data=None, **kwargs):

    s = 0
    for data in fluid_data.values():

        if _is_larger_than_precision(data["mass_fraction"]):
            s += data["wrapper"].s_pT(p, T) * data["mass_fraction"]

    return s


def s_mix_pT_humidair(p, T, fluid_data, **kwargs):
    w = w_mix_fluid_data(fluid_data)
    return HAPropsSI("S", "P", p, "T", T, "W", w)


def v_mix_pT_ideal(p=None, T=None, fluid_data=None, **kwargs):
    molar_fractions = get_molar_fractions(fluid_data)

    d = 0
    for fluid, data in fluid_data.items():

        if _is_larger_than_precision(data["mass_fraction"]):
            pp = p * molar_fractions[fluid]
            d += data["wrapper"].d_pT(pp, T)

    return 1 / d


def v_mix_pT_ideal_cond(p=None, T=None, fluid_data=None, **kwargs):

    water_alias = _get_fluid_alias("H2O", fluid_data)
    if water_alias:
        water_alias = next(iter(water_alias))
        _, molar_fraction_gas, mass_liquid, _, p_sat, pp_water = cond_check(p, T, fluid_data, water_alias)
        # at saturation liquid mass may be zero, but we cannot calculate water properties with pT
        if not _is_larger_than_precision(mass_liquid) and abs(pp_water - p_sat) / p_sat > 1e-6:
            return v_mix_pT_ideal(p, T, fluid_data, **kwargs)
        d = 0
        for fluid, data in fluid_data.items():
            if _is_larger_than_precision(data["mass_fraction"]):
                if fluid == water_alias:
                    if mass_liquid > 0:
                        d += fluid_data[water_alias]["wrapper"].d_QT(0, T) * mass_liquid
                    d += fluid_data[water_alias]["wrapper"].d_QT(1, T) * (1 - mass_liquid)
                else:
                    pp = p * molar_fraction_gas[fluid]
                    d += data["wrapper"].d_pT(pp, T) * (1 - mass_liquid)
        return 1 / d
    else:
        return v_mix_pT_ideal(p, T, fluid_data, **kwargs)


def v_mix_pT_incompressible(p=None, T=None, fluid_data=None, **kwargs):

    v = 0
    for data in fluid_data.values():
        if _is_larger_than_precision(data["mass_fraction"]):
            v += 1 / data["wrapper"].d_pT(p, T) * data["mass_fraction"]

    return v


def v_mix_pT_humidair(p, T, fluid_data, **kwargs):
    w = w_mix_fluid_data(fluid_data)
    return HAPropsSI("V", "P", p, "T", T, "W", w)


def viscosity_mix_pT_ideal(p=None, T=None, fluid_data=None, **kwargs):
    r"""
    Calculate dynamic viscosity from pressure and temperature.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

    T : float
        Temperature T / K.

    Returns
    -------
    visc : float
        Dynamic viscosity visc / Pa s.

    Note
    ----
    Calculation for fluid mixtures.

    .. math::

        \eta_{mix}(p,T)=\frac{\sum_{i} \left( \eta(p,T,fluid_{i}) \cdot y_{i}
        \cdot \sqrt{M_{i}} \right)}
        {\sum_{i} \left(y_{i} \cdot \sqrt{M_{i}} \right)}\;
        \forall i \in \text{fluid components}\\
        y: \text{volume fraction}\\
        M: \text{molar mass}

    Reference: :cite:`Herning1936`.
    """
    molar_fractions = get_molar_fractions(fluid_data)

    a = 0
    b = 0
    for fluid, data in fluid_data.items():
        if _is_larger_than_precision(data["mass_fraction"]):
            bi = molar_fractions[fluid] * data["wrapper"]._molar_mass ** 0.5
            b += bi
            a += bi * data["wrapper"].viscosity_pT(p, T)

    return a / b


def viscosity_mix_pT_incompressible(p=None, T=None, fluid_data=None, **kwargs):

    viscosity = 0
    for data in fluid_data.values():
        if _is_larger_than_precision(data["mass_fraction"]):
            viscosity += data["wrapper"].viscosity_pT(p, T) * data["mass_fraction"]

    return viscosity


def viscosity_mix_pT_humidair(p, T, fluid_data, **kwargs):
    w = w_mix_fluid_data(fluid_data)
    return HAPropsSI("Visc", "P", p, "T", T, "W", w)


def exergy_chemical_ideal_cond(pamb, Tamb, fluid_data, Chem_Ex):

    molar_fractions = get_molar_fractions(fluid_data)
    water_alias = _get_fluid_alias("H2O", fluid_data)
    if water_alias:
        water_alias = next(iter(water_alias))
        _, molar_fractions_gas, _, molar_liquid, _, _ = cond_check(
            pamb, Tamb, fluid_data, water_alias
        )
    else:
        molar_fractions_gas = molar_fractions
        molar_liquid = 0

    ex_cond = 0
    ex_dry = 0
    for fluid, x in molar_fractions_gas.items():
        if x == 0:
            continue

        fluid_aliases = FLUID_ALIASES.get_fluid(fluid)

        if molar_liquid > 0 and "water" in fluid_aliases:
            y = [
                Chem_Ex[k][2] for k in fluid_aliases if k in Chem_Ex
            ]
            ex_cond += molar_liquid * y[0]

        y = [Chem_Ex[k][3] for k in fluid_aliases if k in Chem_Ex]
        ex_dry += x * y[0] + Tamb * gas_constants['uni'] * 1e-3 * x * math.log(x)

    ex_chemical = ex_cond + ex_dry * (1 - molar_liquid)
    ex_chemical *= 1 / calc_molar_mass_mixture(
        fluid_data, molar_fractions
    )

    return ex_chemical * 1e3  # Data from Chem_Ex are in kJ / mol


def _get_fluid_alias(fluid, fluid_data):
    return (
        FLUID_ALIASES.get_fluid(fluid)
        & set([
            f for f in fluid_data
            if _is_larger_than_precision(fluid_data[f]["mass_fraction"])
        ])
    )


def cond_check(p, T, fluid_data, water_alias):
    """Check if water is partially condensing in gaseous mixture.

    Parameters
    ----------
    p : float
        pressure of mixture
    T : float
        temperature of mixture
    fluid_data : dict
        Dictionary of fluid data:
        :code:`{fluid_name: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`
    water_alias : str
        label of the water in the fluid_data dictionary

    Returns
    -------
    tuple
        Tuple containing gas phase mass specific and molar specific
        compositions and overall liquid water mass fraction, as well as
        saturation pressure of water and partial pressure of water in gas phase
    """
    molar_fractions = get_molar_fractions(fluid_data)
    molar_fractions_gas = molar_fractions
    mass_fractions_gas = {f: v["mass_fraction"] for f, v in fluid_data.items()}
    water_mass_liquid = 0
    water_molar_liquid = 0
    # make something up here for now
    p_sat = p / 2
    pp_water = p / 3

    if fluid_data[water_alias]["wrapper"]._is_below_T_critical(T):
        p_sat = fluid_data[water_alias]["wrapper"].p_sat(T)
        pp_water = p * molar_fractions[water_alias]

        if p_sat < pp_water:
            water_molar_gas = (1 - molar_fractions[water_alias]) / (p / p_sat - 1)
            water_molar_liquid = molar_fractions[water_alias] - water_molar_gas
            x_gas_sum = 1 - water_molar_liquid

            molar_fractions_gas = {f: x / x_gas_sum for f, x in molar_fractions.items()}
            molar_fractions_gas[water_alias] = water_molar_gas / x_gas_sum

            water_mass_liquid = (
                water_molar_liquid
                * fluid_data[water_alias]["wrapper"]._molar_mass
                / calc_molar_mass_mixture(fluid_data, molar_fractions)
            )

            molar_mass_mixture = calc_molar_mass_mixture(fluid_data, molar_fractions_gas)
            mass_fractions_gas = {
                fluid: (
                    x / molar_mass_mixture
                    * fluid_data[fluid]["wrapper"]._molar_mass
                )
                for fluid, x in molar_fractions_gas.items()
            }

    return mass_fractions_gas, molar_fractions_gas, water_mass_liquid, water_molar_liquid, p_sat, pp_water


class MixingRuleRegistry:
    """Registry of mixing rule functions for multi-component fluid properties.

    Use :meth:`register` to add custom mixing rules at runtime.

    Example
    -------
    Register a custom enthalpy mixing rule:

    >>> def my_h_pT(p, T, fluid_data, **kwargs): ...
    >>> MIXING_RULES.register("my-rule", h_pT=my_h_pT)
    """

    def __init__(self):
        self._h_pT = {}
        self._s_pT = {}
        self._v_pT = {}
        self._viscosity_pT = {}
        self._T_ph = {}
        self._T_ps = {}
        self._exergy_chemical = {}

    def register(
        self, name, *, h_pT=None, s_pT=None, v_pT=None,
        viscosity_pT=None, exergy_chemical=None,
        T_ph_inversion=True, T_ps_inversion=True,
    ):
        """Register a mixing rule.

        Parameters
        ----------
        name : str
            Mixing rule identifier used as the *mixing_rule* argument.
        h_pT : callable, optional
            :code:`h(p, T, fluid_data, **kwargs) -> float`
        s_pT : callable, optional
            :code:`s(p, T, fluid_data, **kwargs) -> float`
        v_pT : callable, optional
            :code:`v(p, T, fluid_data, **kwargs) -> float`
        viscosity_pT : callable, optional
            :code:`visc(p, T, fluid_data, **kwargs) -> float`
        exergy_chemical : callable, optional
            :code:`ex(pamb, Tamb, fluid_data, Chem_Ex) -> float`
        T_ph_inversion : bool
            When *True* (default), *h_pT* is also registered as the Newton
            residual for :code:`T(p, h)` inversion.  Set to *False* when the
            function is not monotonic in T (e.g. :code:`"forced-gas"`).
        T_ps_inversion : bool
            Analogous flag for :code:`T(p, s)` inversion via *s_pT*.
        """
        if h_pT is not None:
            self._h_pT[name] = h_pT
            if T_ph_inversion:
                self._T_ph[name] = h_pT
        if s_pT is not None:
            self._s_pT[name] = s_pT
            if T_ps_inversion:
                self._T_ps[name] = s_pT
        if v_pT is not None:
            self._v_pT[name] = v_pT
        if viscosity_pT is not None:
            self._viscosity_pT[name] = viscosity_pT
        if exergy_chemical is not None:
            self._exergy_chemical[name] = exergy_chemical

    def _get(self, registry, name, label):
        if name not in registry:
            available = sorted(registry.keys())
            msg = (
                f"The mixing rule '{name}' is not available for the fluid "
                f"property function for {label}. Available rules are '"
                + "', '".join(available) + "'."
            )
            logger.exception(msg)
            raise KeyError(msg)
        return registry[name]

    def h_pT(self, name):
        return self._get(self._h_pT, name, "enthalpy")

    def s_pT(self, name):
        return self._get(self._s_pT, name, "entropy")

    def v_pT(self, name):
        return self._get(self._v_pT, name, "specific volume")

    def viscosity_pT(self, name):
        return self._get(self._viscosity_pT, name, "viscosity")

    def T_ph(self, name):
        return self._get(self._T_ph, name, "temperature (from enthalpy)")

    def T_ps(self, name):
        return self._get(self._T_ps, name, "temperature (from entropy)")

    def exergy_chemical(self, name):
        return self._get(self._exergy_chemical, name, "chemical exergy")


MIXING_RULES = MixingRuleRegistry()

MIXING_RULES.register(
    "ideal",
    h_pT=h_mix_pT_ideal,
    s_pT=s_mix_pT_ideal,
    v_pT=v_mix_pT_ideal,
    viscosity_pT=viscosity_mix_pT_ideal,
)
MIXING_RULES.register(
    "ideal-cond",
    h_pT=h_mix_pT_ideal_cond,
    s_pT=s_mix_pT_ideal_cond,
    v_pT=v_mix_pT_ideal_cond,
    viscosity_pT=viscosity_mix_pT_ideal,
    exergy_chemical=exergy_chemical_ideal_cond,
)
MIXING_RULES.register(
    "incompressible",
    h_pT=h_mix_pT_incompressible,
    s_pT=s_mix_pT_incompressible,
    v_pT=v_mix_pT_incompressible,
    viscosity_pT=viscosity_mix_pT_incompressible,
)
MIXING_RULES.register(
    "forced-gas",
    h_pT=h_mix_pT_forced_gas,
    T_ph_inversion=False,
)
MIXING_RULES.register(
    "humidair",
    h_pT=h_mix_pT_humidair,
    s_pT=s_mix_pT_humidair,
    v_pT=v_mix_pT_humidair,
    viscosity_pT=viscosity_mix_pT_humidair,
)
