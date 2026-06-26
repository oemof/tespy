# -*- coding: utf-8

"""Module for fluid property mixture routines.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/tools/fluid_properties/mixtures.py

SPDX-License-Identifier: MIT
"""

import CoolProp as _CP
from CoolProp.CoolProp import HAPropsSI

from tespy.tools.global_vars import FLUID_ALIASES
from tespy.tools.logger import logger

from .helpers import _is_larger_than_precision
from .helpers import calc_molar_mass_mixture
from .helpers import get_molar_fractions


def h_mix_pT_ideal(p=None, T=None, fluid_data=None, **kwargs):
    r"""
    Calculate specific enthalpy of an ideal gas mixture.

    Applies Dalton's law to assign partial pressures and sums the
    mass-fraction-weighted component enthalpies.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    **kwargs
        Ignored; present for interface compatibility.

    Returns
    -------
    float
        Specific enthalpy in J/kg.

    Notes
    -----
    Molar fractions are derived from mass fractions:

    .. math::

        y_i = \frac{x_i / M_i}{\sum_j x_j / M_j}

    Each component is evaluated at its partial pressure
    :math:`p_i = p \cdot y_i`. The mixture enthalpy is:

    .. math::

        h_\text{mix}(p, T) = \sum_i x_i \cdot h_i(p_i, T)
    """
    molar_fractions = get_molar_fractions(fluid_data)

    h = 0
    for fluid, data in fluid_data.items():

        if _is_larger_than_precision(data["mass_fraction"]):
            pp = p * molar_fractions[fluid]
            h += data["wrapper"].h_pT(pp, T) * data["mass_fraction"]

    return h


def h_mix_pT_ideal_cond(p=None, T=None, fluid_data=None, **kwargs):
    r"""
    Calculate specific enthalpy of an ideal gas mixture with water condensation.

    Extends :func:`h_mix_pT_ideal` by checking whether the water vapour
    partial pressure exceeds the saturation pressure at *T*. If condensation
    occurs the water content is split between a saturated gas phase and a
    liquid phase, and the enthalpy is computed accordingly.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    **kwargs
        Ignored; present for interface compatibility.

    Returns
    -------
    float
        Specific enthalpy in J/kg.

    Notes
    -----
    If no H2O is present, or the water partial pressure is below saturation,
    the result equals :func:`h_mix_pT_ideal`.

    When :math:`p_\text{sat}(T) < p \cdot y_\text{H2O}`, the liquid water
    mass fraction :math:`x_\text{liq}` is determined by :func:`cond_check`.
    The enthalpy is then:

    .. math::

        h_\text{mix} =
            x_\text{liq} \cdot h_\text{H2O}(Q{=}0,\,T)
            + (1 - x_\text{liq}) \left[
                x_\text{H2O}^\text{gas} \cdot h_\text{H2O}(Q{=}1,\,T)
                + \sum_{i \neq \text{H2O}}
                    x_i^\text{gas} \cdot h_i\!\left(p_i^\text{gas},\,T\right)
            \right]

    where :math:`x_i^\text{gas}` and :math:`p_i^\text{gas}` are the gas-phase
    mass fractions and partial pressures after condensate removal.
    """

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
    r"""
    Calculate specific enthalpy treating all H2O as gas phase.

    Like :func:`h_mix_pT_ideal`, but forces the water component to be
    evaluated as saturated vapour whenever the temperature is at or below the
    saturation point at the water partial pressure.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    **kwargs
        Ignored; present for interface compatibility.

    Returns
    -------
    float
        Specific enthalpy in J/kg.

    Notes
    -----
    For all components other than H2O the partial-pressure enthalpy is used:

    .. math::

        h_i = h_i(p_i, T), \quad p_i = p \cdot y_i

    For H2O, when :math:`T \le T_\text{sat}(p_\text{H2O})`, the saturated
    vapour enthalpy is substituted:

    .. math::

        h_\text{H2O} = h_\text{H2O}(Q{=}1,\,T)

    The mixture enthalpy is:

    .. math::

        h_\text{mix}(p, T) = \sum_i x_i \cdot h_i
    """
    molar_fractions = get_molar_fractions(fluid_data)
    water_aliases = _get_fluid_alias("H2O", fluid_data)

    h = 0
    for fluid, data in fluid_data.items():

        if _is_larger_than_precision(data["mass_fraction"]):
            pp = p * molar_fractions[fluid]
            if fluid in water_aliases and pp >= data["wrapper"]._p_min:
                if T <= data["wrapper"].T_sat(pp):
                    h += data["wrapper"].h_QT(1, T) * data["mass_fraction"]
                else:
                    h += data["wrapper"].h_pT(pp, T) * data["mass_fraction"]
            else:
                h += data["wrapper"].h_pT(pp, T) * data["mass_fraction"]

    return h


def h_mix_pT_incompressible(p, T, fluid_data, **kwargs):
    r"""
    Calculate specific enthalpy of an incompressible mixture.

    Each component is evaluated at the full mixture pressure (no partial
    pressure splitting) and contributions are summed by mass fraction.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    **kwargs
        Ignored; present for interface compatibility.

    Returns
    -------
    float
        Specific enthalpy in J/kg.

    Notes
    -----
    .. math::

        h_\text{mix}(p, T) = \sum_i x_i \cdot h_i(p, T)
    """

    h = 0
    for data in fluid_data.values():
        if _is_larger_than_precision(data["mass_fraction"]):
            h += data["wrapper"].h_pT(p, T) * data["mass_fraction"]

    return h


def w_mix_fluid_data(fluid_data):
    r"""
    Return the humidity ratio of a moist-air mixture from fluid composition.

    Parameters
    ----------
    fluid_data : dict
        Fluid property data containing H2O and air components.

    Returns
    -------
    float
        Humidity ratio in kg\ :sub:`water` / kg\ :sub:`dry air`.

    Notes
    -----
    .. math::

        w = \frac{x_\text{H2O}}{x_\text{air}}
    """

    water_alias = _get_fluid_alias("H2O", fluid_data)
    water_alias = next(iter(water_alias))

    air_alias = _get_fluid_alias("air", fluid_data)
    air_alias = next(iter(air_alias))

    return (
        fluid_data[water_alias]["mass_fraction"]
        / fluid_data[air_alias]["mass_fraction"]
    )

def w_mix_pTrh_humidair(p, T, rh):
    r"""
    Return humidity ratio at given pressure, temperature and relative humidity.

    Thin wrapper around :func:`CoolProp.CoolProp.HAPropsSI`.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    T : float
        Temperature in K.
    rh : float
        Relative humidity (0–1).

    Returns
    -------
    float
        Humidity ratio in kg\ :sub:`water` / kg\ :sub:`dry air`.
    """
    return HAPropsSI("W", "P", p, "T", T, "RH", rh)  # kg water/kg dry air

def w_mix_pT_humidair(p, T, fluid_data, **kwargs):
    r"""
    Return effective humidity ratio for a humid-air mixture at given *p* and *T*.

    Compares the mixture humidity ratio derived from *fluid_data* with the
    saturation limit at RH = 1 and caps it if condensation would occur.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data containing H2O and air components.
    **kwargs
        Ignored; present for interface compatibility.

    Returns
    -------
    float
        Humidity ratio in kg\ :sub:`water` / kg\ :sub:`dry air`.

    Notes
    -----
    .. math::

        w = \min\!\left(w_\text{fluid\_data},\; w_\text{sat}(p, T)\right)

    where :math:`w_\text{sat}(p, T) = w(p, T, \text{RH}{=}1)`.
    """
    w_def = w_mix_fluid_data(fluid_data)
    w_max = w_mix_pTrh_humidair(p, T, 1.0)
    if w_def > w_max:
        _msg = f"Humidity ratio {w_def:.4f} exceeds maximum value of {w_max:.4f} for given p and T. Check fluid composition."
        return w_max
    return w_def

def h_mix_pT_humidair(p, T, fluid_data, **kwargs):
    r"""
    Calculate specific enthalpy of a humid-air mixture.

    Delegates to :func:`CoolProp.CoolProp.HAPropsSI` using the humidity ratio
    obtained from :func:`w_mix_pT_humidair`.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data containing H2O and air components.
    **kwargs
        Forwarded to :func:`w_mix_pT_humidair`.

    Returns
    -------
    float
        Specific enthalpy in J/kg.
    """
    w = w_mix_pT_humidair(p, T, fluid_data, **kwargs)
    return HAPropsSI("H", "P", p, "T", T, "W", w)

def w_mix_phrh_humidair(p, h, rh):
    r"""
    Return humidity ratio at given pressure, enthalpy and relative humidity.

    Thin wrapper around :func:`CoolProp.CoolProp.HAPropsSI`.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    h : float
        Specific enthalpy in J/kg.
    rh : float
        Relative humidity (0–1).

    Returns
    -------
    float
        Humidity ratio in kg\ :sub:`water` / kg\ :sub:`dry air`.
    """
    return HAPropsSI("W", "P", p, "H", h, "RH", rh)  # kg water/kg dry air

def w_mix_ph_humidair(p, h, fluid_data, **kwargs):
    r"""
    Return effective humidity ratio for a humid-air mixture at given *p* and *h*.

    Compares the mixture humidity ratio derived from *fluid_data* with the
    saturation limit at RH = 1 and caps it if condensation would occur.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    h : float
        Specific enthalpy in J/kg.
    fluid_data : dict
        Fluid property data containing H2O and air components.
    **kwargs
        Ignored; present for interface compatibility.

    Returns
    -------
    float
        Humidity ratio in kg\ :sub:`water` / kg\ :sub:`dry air`.

    Notes
    -----
    .. math::

        w = \min\!\left(w_\text{fluid\_data},\; w_\text{sat}(p, h)\right)

    where :math:`w_\text{sat}(p, h) = w(p, h, \text{RH}{=}1)`.
    """
    w_def = w_mix_fluid_data(fluid_data)
    w_max = w_mix_phrh_humidair(p, h, 1.0)
    if w_def > w_max:
        _msg = f"Humidity ratio {w_def:.4f} exceeds maximum value of {w_max:.4f} for given p and T. Check fluid composition."
        return w_max
    return w_def

def w_mix_psrh_humidair(p, s, rh):
    r"""
    Return humidity ratio at given pressure, entropy and relative humidity.

    Thin wrapper around :func:`CoolProp.CoolProp.HAPropsSI`.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    s : float
        Specific entropy in J/(kg K).
    rh : float
        Relative humidity (0–1).

    Returns
    -------
    float
        Humidity ratio in kg\ :sub:`water` / kg\ :sub:`dry air`.
    """
    return HAPropsSI("W", "P", p, "S", s, "RH", rh)  # kg water/kg dry air

def w_mix_ps_humidair(p, s, fluid_data, **kwargs):
    r"""
    Return effective humidity ratio for a humid-air mixture at given *p* and *s*.

    Compares the mixture humidity ratio derived from *fluid_data* with the
    saturation limit at RH = 1 and caps it if condensation would occur.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    s : float
        Specific entropy in J/(kg K).
    fluid_data : dict
        Fluid property data containing H2O and air components.
    **kwargs
        Ignored; present for interface compatibility.

    Returns
    -------
    float
        Humidity ratio in kg\ :sub:`water` / kg\ :sub:`dry air`.

    Notes
    -----
    .. math::

        w = \min\!\left(w_\text{fluid\_data},\; w_\text{sat}(p, s)\right)

    where :math:`w_\text{sat}(p, s) = w(p, s, \text{RH}{=}1)`.
    """
    w_def = w_mix_fluid_data(fluid_data)
    w_max = w_mix_psrh_humidair(p, s, 1.0)
    if w_def > w_max:
        _msg = f"Humidity ratio {w_def:.4f} exceeds maximum value of {w_max:.4f} for given p and T. Check fluid composition."
        return w_max
    return w_def

def s_mix_pT_ideal(p=None, T=None, fluid_data=None, **kwargs):
    r"""
    Calculate specific entropy of an ideal gas mixture.

    Applies Dalton's law to assign partial pressures and sums the
    mass-fraction-weighted component entropies.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    **kwargs
        Ignored; present for interface compatibility.

    Returns
    -------
    float
        Specific entropy in J/(kg K).

    Notes
    -----
    Molar fractions are derived from mass fractions:

    .. math::

        y_i = \frac{x_i / M_i}{\sum_j x_j / M_j}

    Each component is evaluated at its partial pressure
    :math:`p_i = p \cdot y_i`. The mixture entropy is:

    .. math::

        s_\text{mix}(p, T) = \sum_i x_i \cdot s_i(p_i, T)
    """
    molar_fractions = get_molar_fractions(fluid_data)

    s = 0
    for fluid, data in fluid_data.items():

        if _is_larger_than_precision(data["mass_fraction"]):
            pp = p * molar_fractions[fluid]
            s += data["wrapper"].s_pT(pp, T) * data["mass_fraction"]

    return s


def s_mix_pT_ideal_cond(p=None, T=None, fluid_data=None, **kwargs):
    r"""
    Calculate specific entropy of an ideal gas mixture with water condensation.

    Extends :func:`s_mix_pT_ideal` by checking whether the water vapour
    partial pressure exceeds the saturation pressure at *T*. If condensation
    occurs the water content is split between a saturated gas phase and a
    liquid phase, and the entropy is computed accordingly.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    **kwargs
        Ignored; present for interface compatibility.

    Returns
    -------
    float
        Specific entropy in J/(kg K).

    Notes
    -----
    If no H2O is present, or the water partial pressure is below saturation,
    the result equals :func:`s_mix_pT_ideal`.

    When :math:`p_\text{sat}(T) < p \cdot y_\text{H2O}`, the liquid water
    mass fraction :math:`x_\text{liq}` is determined by :func:`cond_check`.
    The entropy is then:

    .. math::

        s_\text{mix} =
            x_\text{liq} \cdot s_\text{H2O}(Q{=}0,\,T)
            + (1 - x_\text{liq}) \left[
                x_\text{H2O}^\text{gas} \cdot s_\text{H2O}(Q{=}1,\,T)
                + \sum_{i \neq \text{H2O}}
                    x_i^\text{gas} \cdot s_i\!\left(p_i^\text{gas},\,T\right)
            \right]

    where :math:`x_i^\text{gas}` and :math:`p_i^\text{gas}` are the gas-phase
    mass fractions and partial pressures after condensate removal.
    """

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
    r"""
    Calculate specific entropy of an incompressible mixture.

    Each component is evaluated at the full mixture pressure (no partial
    pressure splitting) and contributions are summed by mass fraction.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    **kwargs
        Ignored; present for interface compatibility.

    Returns
    -------
    float
        Specific entropy in J/(kg K).

    Notes
    -----
    .. math::

        s_\text{mix}(p, T) = \sum_i x_i \cdot s_i(p, T)
    """

    s = 0
    for data in fluid_data.values():

        if _is_larger_than_precision(data["mass_fraction"]):
            s += data["wrapper"].s_pT(p, T) * data["mass_fraction"]

    return s


def s_mix_pT_humidair(p, T, fluid_data, **kwargs):
    r"""
    Calculate specific entropy of a humid-air mixture.

    Delegates to :func:`CoolProp.CoolProp.HAPropsSI` using the humidity ratio
    obtained from :func:`w_mix_fluid_data`.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data containing H2O and air components.
    **kwargs
        Ignored; present for interface compatibility.

    Returns
    -------
    float
        Specific entropy in J/(kg K).
    """
    w = w_mix_fluid_data(fluid_data)
    return HAPropsSI("S", "P", p, "T", T, "W", w)


def v_mix_pT_ideal(p=None, T=None, fluid_data=None, **kwargs):
    r"""
    Calculate specific volume of an ideal gas mixture.

    Sums the partial densities evaluated at each component's partial pressure
    and returns the reciprocal.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    **kwargs
        Ignored; present for interface compatibility.

    Returns
    -------
    float
        Specific volume in m³/kg.

    Notes
    -----
    Using Dalton's law, each component contributes its density at its partial
    pressure :math:`p_i = p \cdot y_i`:

    .. math::

        \rho_\text{mix}(p, T) = \sum_i \rho_i(p_i, T)

        v_\text{mix} = \frac{1}{\rho_\text{mix}}
    """
    molar_fractions = get_molar_fractions(fluid_data)

    d = 0
    for fluid, data in fluid_data.items():

        if _is_larger_than_precision(data["mass_fraction"]):
            pp = p * molar_fractions[fluid]
            d += data["wrapper"].d_pT(pp, T)

    return 1 / d


def v_mix_pT_ideal_cond(p=None, T=None, fluid_data=None, **kwargs):
    r"""
    Calculate specific volume of an ideal gas mixture with water condensation.

    Extends :func:`v_mix_pT_ideal` by checking whether the water vapour
    partial pressure exceeds the saturation pressure at *T*. If condensation
    occurs the water content is split between a saturated gas phase and a
    liquid phase, and the specific volume is computed accordingly.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    **kwargs
        Ignored; present for interface compatibility.

    Returns
    -------
    float
        Specific volume in m³/kg.

    Notes
    -----
    If no H2O is present, or the water partial pressure is below saturation,
    the result equals :func:`v_mix_pT_ideal`.

    When :math:`p_\text{sat}(T) < p \cdot y_\text{H2O}`, the liquid water
    mass fraction :math:`x_\text{liq}` is determined by :func:`cond_check`.
    The gas-phase density is the sum of component mass densities at their
    gas-phase partial pressures (same as :func:`v_mix_pT_ideal`):

    .. math::

        \rho_\text{gas} =
            \rho_\text{H2O}(Q{=}1,\,T)
            + \sum_{i \neq \text{H2O}} \rho_i\!\left(p_i^\text{gas},\,T\right)

    The mixture specific volume is the mass-fraction-weighted sum of phase
    specific volumes:

    .. math::

        v_\text{mix} =
            x_\text{liq} \cdot v_\text{H2O}(Q{=}0,\,T)
            + (1 - x_\text{liq}) \cdot \frac{1}{\rho_\text{gas}}
    """

    water_alias = _get_fluid_alias("H2O", fluid_data)
    if water_alias:
        water_alias = next(iter(water_alias))
        _, molar_fraction_gas, mass_liquid, _, p_sat, pp_water = cond_check(p, T, fluid_data, water_alias)
        # at saturation liquid mass may be zero, but we cannot calculate water properties with pT
        if not _is_larger_than_precision(mass_liquid) and abs(pp_water - p_sat) / p_sat > 1e-6:
            return v_mix_pT_ideal(p, T, fluid_data, **kwargs)
        d_gas = fluid_data[water_alias]["wrapper"].d_QT(1, T)
        for fluid, data in fluid_data.items():
            if _is_larger_than_precision(data["mass_fraction"]) and fluid != water_alias:
                pp = p * molar_fraction_gas[fluid]
                d_gas += data["wrapper"].d_pT(pp, T)

        v_liq = 1 / fluid_data[water_alias]["wrapper"].d_QT(0, T)
        v_gas = 1 / d_gas
        return mass_liquid * v_liq + (1 - mass_liquid) * v_gas
    else:
        return v_mix_pT_ideal(p, T, fluid_data, **kwargs)


def v_mix_pT_incompressible(p=None, T=None, fluid_data=None, **kwargs):
    r"""
    Calculate specific volume of an incompressible mixture.

    Applies volume additivity: each component's specific volume is weighted
    by its mass fraction and contributions are summed.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    **kwargs
        Ignored; present for interface compatibility.

    Returns
    -------
    float
        Specific volume in m³/kg.

    Notes
    -----
    .. math::

        v_\text{mix}(p, T) = \sum_i x_i \cdot v_i(p, T)
    """

    v = 0
    for data in fluid_data.values():
        if _is_larger_than_precision(data["mass_fraction"]):
            v += 1 / data["wrapper"].d_pT(p, T) * data["mass_fraction"]

    return v


def v_mix_pT_humidair(p, T, fluid_data, **kwargs):
    r"""
    Calculate specific volume of a humid-air mixture.

    Delegates to :func:`CoolProp.CoolProp.HAPropsSI` using the humidity ratio
    obtained from :func:`w_mix_fluid_data`.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data containing H2O and air components.
    **kwargs
        Ignored; present for interface compatibility.

    Returns
    -------
    float
        Specific volume in m³/kg.
    """
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
    r"""
    Calculate dynamic viscosity of an incompressible mixture.

    Applies a simple mass-fraction-weighted average of the component
    viscosities at the full mixture pressure.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    **kwargs
        Ignored; present for interface compatibility.

    Returns
    -------
    float
        Dynamic viscosity in Pa s.

    Notes
    -----
    .. math::

        \eta_\text{mix}(p, T) = \sum_i x_i \cdot \eta_i(p, T)
    """

    viscosity = 0
    for data in fluid_data.values():
        if _is_larger_than_precision(data["mass_fraction"]):
            viscosity += data["wrapper"].viscosity_pT(p, T) * data["mass_fraction"]

    return viscosity


def viscosity_mix_pT_humidair(p, T, fluid_data, **kwargs):
    r"""
    Calculate dynamic viscosity of a humid-air mixture.

    Delegates to :func:`CoolProp.CoolProp.HAPropsSI` using the humidity ratio
    obtained from :func:`w_mix_fluid_data`.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data containing H2O and air components.
    **kwargs
        Ignored; present for interface compatibility.

    Returns
    -------
    float
        Dynamic viscosity in Pa s.
    """
    w = w_mix_fluid_data(fluid_data)
    return HAPropsSI("Visc", "P", p, "T", T, "W", w)


def _get_fluid_alias(fluid, fluid_data):
    r"""
    Return the set of aliases present in *fluid_data* for a given fluid name.

    Parameters
    ----------
    fluid : str
        Canonical fluid name (e.g. :code:`"H2O"`, :code:`"air"`).
    fluid_data : dict
        Fluid property data.

    Returns
    -------
    set
        Intersection of known aliases for *fluid* with the keys in
        *fluid_data* that carry a non-negligible mass fraction.
    """
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

    def register(
        self, name, *, h_pT=None, s_pT=None, v_pT=None,
        viscosity_pT=None,
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

_LIBR_AS = _CP.AbstractState("INCOMP", "LiBr")
_LIBR_PROPS_CACHE = (None, None, None, None)  # (p, T, xi, (h, s, v))
_LIBR_T_EPS = 0.01


def _libr_T_limits(fluid_data):
    water = _get_fluid_alias("H2O", fluid_data)
    for fluid, data in fluid_data.items():
        if fluid not in water and _is_larger_than_precision(data["mass_fraction"]):
            return data["wrapper"]._T_min, data["wrapper"]._T_max
    return 273.0, 500.0


def _xi_libr(fluid_data):
    water = _get_fluid_alias("H2O", fluid_data)
    if water:
        return 1.0 - fluid_data[next(iter(water))]["mass_fraction"]
    return sum(d["mass_fraction"] for d in fluid_data.values())


def _xi_sat_libr(p, T):
    r"""LiBr mass fraction at which :math:`p_\text{sat}(T, \xi) = p`.

    Uses :func:`scipy.optimize.brentq` over :math:`\xi \in [0.001, 0.749]`.
    :math:`p_\text{sat}` is monotonically decreasing in :math:`\xi` (higher
    LiBr concentration lowers the water vapour pressure).
    """
    from scipy.optimize import brentq

    def residual(xi):
        _LIBR_AS.set_mass_fractions([xi])
        _LIBR_AS.update(_CP.QT_INPUTS, 0, T)
        return _LIBR_AS.p() - p

    r_lo = residual(0.001)
    if r_lo <= 0:
        return 0.001
    r_hi = residual(0.749)
    if r_hi >= 0:
        return 0.749
    return brentq(residual, 0.001, 0.749, xtol=1e-6)


def _T_sat_libr(p, xi):
    r"""Saturation temperature of LiBr-H2O at pressure *p* and mass fraction *xi*.

    Inverts :func:`_xi_sat_libr` over temperature using
    :func:`scipy.optimize.brentq`.  :math:`p_\text{sat}(T, \xi)` is
    monotonically increasing in *T*.
    """
    from scipy.optimize import brentq

    xi_safe = max(0.001, min(0.749, xi))

    def residual(T):
        _LIBR_AS.set_mass_fractions([xi_safe])
        _LIBR_AS.update(_CP.QT_INPUTS, 0, T)
        return _LIBR_AS.p() - p

    T_min = 273.02
    T_max = 499.98
    if residual(T_min) >= 0:
        return T_min
    if residual(T_max) <= 0:
        return T_max
    return brentq(residual, T_min, T_max, xtol=1e-6)


def _p_sat_libr(T, xi):
    r"""Vapour pressure of LiBr-H2O at temperature *T* and mass fraction *xi*."""
    xi_safe = max(0.001, min(0.749, xi))
    _LIBR_AS.set_mass_fractions([xi_safe])
    _LIBR_AS.update(_CP.QT_INPUTS, 0, T)
    return _LIBR_AS.p()


def _libr_props_compute(p, T, xi, fluid_data):
    r"""Compute :math:`(h, s, v)` for LiBr-H2O at *(p, T, xi)*.

    When :math:`p \geq p_\text{sat}(T, \xi)` the solution is fully liquid and
    properties are read directly from CoolProp's :code:`INCOMP::LiBr` backend.

    When :math:`p < p_\text{sat}(T, \xi)` the solution is in equilibrium with
    water vapour.  The saturated LiBr fraction :math:`\xi_\text{sat}` is found
    via :func:`_xi_sat_libr` such that :math:`p_\text{sat}(T, \xi_\text{sat}) = p`.
    Per-unit-mass balances give liquid fraction
    :math:`m_l = \xi / \xi_\text{sat}` and vapour fraction
    :math:`m_v = 1 - m_l`.  Properties are the weighted sum of the saturated
    liquid solution and pure water vapour.

    """
    _LIBR_AS.set_mass_fractions([xi])
    _LIBR_AS.update(_CP.QT_INPUTS, 0, T)
    p_sat = _LIBR_AS.p()

    if p >= p_sat:
        _LIBR_AS.update(_CP.PT_INPUTS, p, T)
        return _LIBR_AS.hmass(), _LIBR_AS.smass(), 1.0 / _LIBR_AS.rhomass()

    xi_sat = _xi_sat_libr(p, T)
    m_l = xi / xi_sat
    m_v = 1.0 - m_l

    _LIBR_AS.set_mass_fractions([xi_sat])
    _LIBR_AS.update(_CP.QT_INPUTS, 0, T)
    p_sat_new = _LIBR_AS.p()
    _LIBR_AS.update(_CP.PT_INPUTS, p_sat_new + 1.0, T)
    h_l = _LIBR_AS.hmass()
    s_l = _LIBR_AS.smass()
    v_l = 1.0 / _LIBR_AS.rhomass()

    water = _get_fluid_alias("H2O", fluid_data)
    w = fluid_data[next(iter(water))]["wrapper"]
    h_v = w.h_QT(1, T)
    s_v = w.s_QT(1, T)
    v_v = 1.0 / w.d_QT(1, T)

    h0 = m_l * h_l + m_v * h_v
    s0 = m_l * s_l + m_v * s_v
    v0 = m_l * v_l + m_v * v_v
    return h0, s0, v0


def _get_libr_props(p, T, fluid_data):
    global _LIBR_PROPS_CACHE
    xi = _xi_libr(fluid_data)
    cp, cT, cxi, cached = _LIBR_PROPS_CACHE
    if p == cp and T == cT and xi == cxi:
        return cached
    props = _libr_props_compute(p, T, xi, fluid_data)
    _LIBR_PROPS_CACHE = (p, T, xi, props)
    return props


def h_mix_pT_libr_water(p, T, fluid_data, **kwargs):
    r"""
    Calculate specific enthalpy of a LiBr-water solution.

    Uses CoolProp's :code:`INCOMP::LiBr` backend (Patek-Klomfar correlations).
    Handles both the sub-saturated liquid and the two-phase (solution + water
    vapour) region; see :func:`_libr_props_compute` for details.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data containing LiBr and H2O components.
    **kwargs
        Ignored; present for interface compatibility.

    Returns
    -------
    float
        Specific enthalpy in J/kg.
    """
    return _get_libr_props(p, T, fluid_data)[0]


def s_mix_pT_libr_water(p, T, fluid_data, **kwargs):
    r"""
    Calculate specific entropy of a LiBr-water solution.

    Uses CoolProp's :code:`INCOMP::LiBr` backend.  Handles both the
    sub-saturated liquid and the two-phase region; see
    :func:`_libr_props_compute` for details.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data containing LiBr and H2O components.
    **kwargs
        Ignored; present for interface compatibility.

    Returns
    -------
    float
        Specific entropy in J/(kg K).
    """
    return _get_libr_props(p, T, fluid_data)[1]


def v_mix_pT_libr_water(p, T, fluid_data, **kwargs):
    r"""
    Calculate specific volume of a LiBr-water solution.

    Uses CoolProp's :code:`INCOMP::LiBr` backend.  Handles both the
    sub-saturated liquid and the two-phase region; see
    :func:`_libr_props_compute` for details.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data containing LiBr and H2O components.
    **kwargs
        Ignored; present for interface compatibility.

    Returns
    -------
    float
        Specific volume in m³/kg.
    """
    return _get_libr_props(p, T, fluid_data)[2]


def phase_mix_ph_libr_water(p, h, fluid_data):
    r"""Return the phase of a LiBr-water stream given *(p, h)*.

    Returns :code:`"l"` when the solution is fully liquid
    (:math:`p \geq p_\text{sat}(T, \xi)`) and :code:`"tp"` when the system
    is in the two-phase (solution + water vapour) region.
    """
    T = T_mix_ph_libr_water(p, h, fluid_data)
    xi = _xi_libr(fluid_data)
    _LIBR_AS.set_mass_fractions([xi])
    _LIBR_AS.update(_CP.QT_INPUTS, 0, T)
    return "l" if p >= _LIBR_AS.p() else "tp"


def T_mix_ph_libr_water(p, h, fluid_data, T0=None):
    r"""
    Invert :func:`h_mix_pT_libr_water` to recover temperature.

    Uses :func:`scipy.optimize.brentq` over the valid temperature range of
    the LiBr wrapper, inset by :code:`2 * _LIBR_T_EPS` from each boundary.
    This avoids Newton convergence failures near the liquid-to-two-phase
    boundary where :math:`dh/dT` changes sharply.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    h : float
        Specific enthalpy in J/kg.
    fluid_data : dict
        Fluid property data containing LiBr and H2O components.
    T0 : float, optional
        Ignored; present for interface compatibility.

    Returns
    -------
    float
        Temperature in K.
    """
    from scipy.optimize import brentq

    xi = _xi_libr(fluid_data)
    T_lim_min, T_lim_max = _libr_T_limits(fluid_data)
    T_lo = T_lim_min + 2 * _LIBR_T_EPS
    T_hi = T_lim_max - 2 * _LIBR_T_EPS

    def residual(T):
        return _libr_props_compute(p, T, xi, fluid_data)[0] - h

    # clamp value to be within T_lo/T_hi to prevent inter-iteration errors
    h_lo = _libr_props_compute(p, T_lo, xi, fluid_data)[0]
    if h <= h_lo:
        return T_lo
    h_hi = _libr_props_compute(p, T_hi, xi, fluid_data)[0]
    if h >= h_hi:
        return T_hi

    return brentq(residual, T_lo, T_hi, xtol=1e-6)


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
MIXING_RULES.register(
    "libr_water",
    h_pT=h_mix_pT_libr_water,
    s_pT=s_mix_pT_libr_water,
    v_pT=v_mix_pT_libr_water,
)
