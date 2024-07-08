# -*- coding: utf-8

"""Module for fluid property mixture routines.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/tools/fluid_properties/mixtures.py

SPDX-License-Identifier: MIT
"""

import math

import CoolProp as CP

from tespy.tools.global_vars import gas_constants

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

    water_alias = _water_in_mixture(fluid_data)
    if water_alias:
        water_alias = next(iter(water_alias))
        mass_fractions_gas, molar_fraction_gas, mass_liquid, molar_liquid = cond_check(p, T, fluid_data, water_alias)
        if not _is_larger_than_precision(mass_liquid):
            return h_mix_pT_ideal(p, T, fluid_data, **kwargs)
        h = 0
        for fluid, data in fluid_data.items():
            if _is_larger_than_precision(data["mass_fraction"]):
                if fluid == water_alias:
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


def s_mix_pT_ideal(p=None, T=None, fluid_data=None, **kwargs):
    molar_fractions = get_molar_fractions(fluid_data)

    s = 0
    for fluid, data in fluid_data.items():

        if _is_larger_than_precision(data["mass_fraction"]):
            pp = p * molar_fractions[fluid]
            s += data["wrapper"].s_pT(pp, T) * data["mass_fraction"]

    return s


def s_mix_pT_ideal_cond(p=None, T=None, fluid_data=None, **kwargs):

    water_alias = _water_in_mixture(fluid_data)
    if water_alias:
        water_alias = next(iter(water_alias))
        mass_fractions_gas, molar_fraction_gas, mass_liquid, molar_liquid = cond_check(p, T, fluid_data, water_alias)
        if mass_liquid == 0:
            return s_mix_pT_ideal(p, T, fluid_data, **kwargs)
        s = 0
        for fluid, data in fluid_data.items():
            if _is_larger_than_precision(data["mass_fraction"]):
                if fluid == water_alias:
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


def v_mix_pT_ideal(p=None, T=None, fluid_data=None, **kwargs):
    molar_fractions = get_molar_fractions(fluid_data)

    d = 0
    for fluid, data in fluid_data.items():

        if _is_larger_than_precision(data["mass_fraction"]):
            pp = p * molar_fractions[fluid]
            d += data["wrapper"].d_pT(pp, T)

    return 1 / d


def v_mix_pT_ideal_cond(p=None, T=None, fluid_data=None, **kwargs):

    water_alias = _water_in_mixture(fluid_data)
    if water_alias:
        water_alias = next(iter(water_alias))
        mass_fractions_gas, molar_fraction_gas, mass_liquid, molar_liquid = cond_check(p, T, fluid_data, water_alias)
        if mass_liquid == 0:
            return v_mix_pT_ideal(p, T, fluid_data, **kwargs)
        d = 0
        for fluid, data in fluid_data.items():
            if _is_larger_than_precision(data["mass_fraction"]):
                if fluid == water_alias:
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


def exergy_chemical_ideal_cond(pamb, Tamb, fluid_data, Chem_Ex):

    molar_fractions = get_molar_fractions(fluid_data)
    water_alias = _water_in_mixture(fluid_data)
    if water_alias:
        water_alias = next(iter(water_alias))
        _, molar_fractions_gas, _, molar_liquid = cond_check(
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

        fluid_aliases = fluid_data[fluid]["wrapper"]._aliases

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


def _water_in_mixture(fluid_data):
    water_aliases = set(CP.CoolProp.get_aliases("H2O"))
    return water_aliases & set([f for f in fluid_data if _is_larger_than_precision(fluid_data[f]["mass_fraction"])])


def cond_check(p, T, fluid_data, water_alias):
    """Check if water is partially condensing in gaseous mixture.

    Parameters
    ----------
    y_i : dict
        Mass specific fluid composition.
    x_i : dict
        Mole specific fluid composition.
    p : float
        Pressure of mass flow.
    n : float
        Molar mass flow.
    T : float
        Temperature of mass flow.

    Returns
    -------
    tuple
        Tuple containing gas phase mass specific and molar specific
        compositions and overall liquid water mass fraction.
    """
    molar_fractions = get_molar_fractions(fluid_data)
    molar_fractions_gas = molar_fractions
    mass_fractions_gas = {f: v["mass_fraction"] for f, v in fluid_data.items()}
    water_mass_liquid = 0
    water_molar_liquid = 0

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

    return mass_fractions_gas, molar_fractions_gas, water_mass_liquid, water_molar_liquid


T_MIX_PH_REVERSE = {
    "ideal": h_mix_pT_ideal,
    "ideal-cond": h_mix_pT_ideal_cond,
    "incompressible": h_mix_pT_incompressible
}


T_MIX_PS_REVERSE = {
    "ideal": s_mix_pT_ideal,
    "ideal-cond": s_mix_pT_ideal_cond,
    "incompressible": s_mix_pT_incompressible
}


H_MIX_PT_DIRECT = {
    "ideal": h_mix_pT_ideal,
    "ideal-cond": h_mix_pT_ideal_cond,
    "incompressible": h_mix_pT_incompressible,
    "forced-gas": h_mix_pT_forced_gas
}


S_MIX_PT_DIRECT = {
    "ideal": s_mix_pT_ideal,
    "ideal-cond": s_mix_pT_ideal_cond,
    "incompressible": s_mix_pT_incompressible
}


V_MIX_PT_DIRECT = {
    "ideal": v_mix_pT_ideal,
    "ideal-cond": v_mix_pT_ideal_cond,
    "incompressible": v_mix_pT_incompressible
}


VISCOSITY_MIX_PT_DIRECT = {
    "ideal": viscosity_mix_pT_ideal,
    "ideal-cond": viscosity_mix_pT_ideal,
    "incompressible": viscosity_mix_pT_incompressible
}

EXERGY_CHEMICAL = {
    "ideal-cond": exergy_chemical_ideal_cond,
}
