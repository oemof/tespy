# -*- coding: utf-8

"""Module for fluid property functions.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/tools/fluid_properties/functions.py

SPDX-License-Identifier: MIT
"""

from .helpers import _check_mixing_rule
from .helpers import get_number_of_fluids
from .helpers import get_pure_fluid
from .helpers import inverse_temperature_mixture
from .mixtures import EXERGY_CHEMICAL
from .mixtures import H_MIX_PT_DIRECT
from .mixtures import S_MIX_PT_DIRECT
from .mixtures import T_MIX_PH_REVERSE
from .mixtures import T_MIX_PS_REVERSE
from .mixtures import V_MIX_PT_DIRECT
from .mixtures import VISCOSITY_MIX_PT_DIRECT


def isentropic(p_1, h_1, p_2, fluid_data, mixing_rule=None, T0=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["wrapper"].isentropic(p_1, h_1, p_2)
    else:
        s_1 = s_mix_ph(p_1, h_1, fluid_data, mixing_rule)
        T_2 = T_mix_ps(p_2, s_1, fluid_data, mixing_rule)
        return h_mix_pT(p_2, T_2, fluid_data, mixing_rule)


def calc_physical_exergy(h, s, p, pamb, Tamb, fluid_data, mixing_rule=None, T0=None):
    r"""
    Calculate specific physical exergy.

    Physical exergy is allocated to a thermal and a mechanical share according
    to :cite:`Morosuk2019`.

    Parameters
    ----------
    pamb : float
        Ambient pressure p0 / Pa.

    Tamb : float
        Ambient temperature T0 / K.

    Returns
    -------
    e_ph : tuple
        Specific thermal and mechanical exergy
        (:math:`e^\mathrm{T}`, :math:`e^\mathrm{M}`) in J / kg.

        .. math::

            e^\mathrm{T} = \left( h - h \left( p, T_0 \right) \right) -
            T_0 \cdot \left(s - s\left(p, T_0\right)\right)

            e^\mathrm{M}=\left(h\left(p,T_0\right)-h\left(p_0,T_0\right)\right)
            -T_0\cdot\left(s\left(p, T_0\right)-s\left(p_0,T_0\right)\right)

            e^\mathrm{PH} = e^\mathrm{T} + e^\mathrm{M}
    """
    h_T0_p = h_mix_pT(p, Tamb, fluid_data, mixing_rule)
    s_T0_p = s_mix_pT(p, Tamb, fluid_data, mixing_rule)
    ex_therm = (h - h_T0_p) - Tamb * (s - s_T0_p)
    h0 = h_mix_pT(pamb, Tamb, fluid_data, mixing_rule)
    s0 = s_mix_pT(pamb, Tamb, fluid_data, mixing_rule)
    ex_mech = (h_T0_p - h0) - Tamb * (s_T0_p - s0)
    return ex_therm, ex_mech


def calc_chemical_exergy(pamb, Tamb, fluid_data, Chem_Ex, mixing_rule=None, T0=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        fluid_aliases = pure_fluid["wrapper"]._aliases
        y = [Chem_Ex[k][Chem_Ex[k][4]] for k in fluid_aliases if k in Chem_Ex]
        return y[0] / pure_fluid["wrapper"]._molar_mass * 1e3
    else:
        _check_mixing_rule(mixing_rule, EXERGY_CHEMICAL, "chemical exergy")
        return EXERGY_CHEMICAL[mixing_rule](pamb, Tamb, fluid_data, Chem_Ex)


def T_mix_ph(p, h, fluid_data, mixing_rule=None, T0=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["wrapper"].T_ph(p, h)
    else:
        _check_mixing_rule(mixing_rule, T_MIX_PH_REVERSE, "temperature (from enthalpy)")
        kwargs = {
            "p": p, "target_value": h, "fluid_data": fluid_data, "T0": T0,
            "f": T_MIX_PH_REVERSE[mixing_rule]
        }
        return inverse_temperature_mixture(**kwargs)


def dT_mix_pdh(p, h, fluid_data, mixing_rule=None, T0=None):
    d = 1e-1
    upper = T_mix_ph(p, h + d, fluid_data, mixing_rule=mixing_rule, T0=T0)
    lower = T_mix_ph(p, h - d, fluid_data, mixing_rule=mixing_rule, T0=upper)
    return (upper - lower) / (2 * d)


def dT_mix_dph(p, h, fluid_data, mixing_rule=None, T0=None):
    d = 1e-1
    upper = T_mix_ph(p + d, h, fluid_data, mixing_rule=mixing_rule, T0=T0)
    lower = T_mix_ph(p - d, h, fluid_data, mixing_rule=mixing_rule, T0=upper)
    return (upper - lower) / (2 * d)


def dT_mix_ph_dfluid(p, h, fluid, fluid_data, mixing_rule=None, T0=None):
    d = 1e-5
    fluid_data[fluid]["mass_fraction"] += d
    upper = T_mix_ph(p, h, fluid_data, mixing_rule=mixing_rule, T0=T0)
    fluid_data[fluid]["mass_fraction"] -= 2 * d
    lower = T_mix_ph(p, h, fluid_data, mixing_rule=mixing_rule, T0=upper)
    fluid_data[fluid]["mass_fraction"] += d
    return (upper - lower) / (2 * d)


def h_mix_pT(p, T, fluid_data, mixing_rule=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["wrapper"].h_pT(p, T)
    else:
        _check_mixing_rule(mixing_rule, H_MIX_PT_DIRECT, "enthalpy")
        return H_MIX_PT_DIRECT[mixing_rule](p, T, fluid_data)


def h_mix_pQ(p, Q, fluid_data, mixing_rule=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["wrapper"].h_pQ(p, Q)
    else:
        msg = "Saturation function cannot be called on mixtures."
        raise ValueError(msg)


def dh_mix_dpQ(p, Q, fluid_data, mixing_rule=None):
    d = 0.1
    upper = h_mix_pQ(p + d, Q, fluid_data)
    lower = h_mix_pQ(p - d, Q, fluid_data)
    return (upper - lower) / (2 * d)


def Q_mix_ph(p, h, fluid_data, mixing_rule=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["wrapper"].Q_ph(p, h)
    else:
        msg = "Saturation function cannot be called on mixtures."
        raise ValueError(msg)

def phase_mix_ph(p, h, fluid_data, mixing_rule=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["wrapper"].phase_ph(p, h)
    else:
        msg = "State function cannot be called on mixtures."
        raise ValueError(msg)


def p_sat_T(T, fluid_data, mixing_rule=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["wrapper"].p_sat(T)
    else:
        msg = "Saturation function cannot be called on mixtures."
        raise ValueError(msg)


def T_sat_p(p, fluid_data, mixing_rule=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["wrapper"].T_sat(p)
    else:
        msg = "Saturation function cannot be called on mixtures."
        raise ValueError(msg)


def dT_sat_dp(p, fluid_data, mixing_rule=None):
    d = 0.01
    upper = T_sat_p(p + d, fluid_data)
    lower = T_sat_p(p - d, fluid_data)
    return (upper - lower) / (2 * d)


def s_mix_ph(p, h, fluid_data, mixing_rule=None, T0=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["wrapper"].s_ph(p, h)
    else:
        T = T_mix_ph(p, h , fluid_data, mixing_rule, T0)
        return s_mix_pT(p, T, fluid_data, mixing_rule)



def s_mix_pT(p, T, fluid_data, mixing_rule=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["wrapper"].s_pT(p, T)
    else:
        _check_mixing_rule(mixing_rule, S_MIX_PT_DIRECT, "entropy")
        return S_MIX_PT_DIRECT[mixing_rule](p, T, fluid_data)


def T_mix_ps(p, s, fluid_data, mixing_rule=None, T0=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["wrapper"].T_ps(p, s)
    else:
        _check_mixing_rule(mixing_rule, T_MIX_PS_REVERSE, "temperature (from entropy)")
        kwargs = {
            "p": p, "target_value": s, "fluid_data": fluid_data, "T0": T0,
            "f": T_MIX_PS_REVERSE[mixing_rule]
        }
        return inverse_temperature_mixture(**kwargs)


def v_mix_ph(p, h, fluid_data, mixing_rule=None, T0=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return 1 / pure_fluid["wrapper"].d_ph(p, h)
    else:
        T = T_mix_ph(p, h , fluid_data, mixing_rule, T0)
        return v_mix_pT(p, T, fluid_data, mixing_rule)


def dv_mix_dph(p, h, fluid_data, mixing_rule=None, T0=None):
    d = 1e-1
    upper = v_mix_ph(p + d, h, fluid_data, mixing_rule=mixing_rule, T0=T0)
    lower = v_mix_ph(p - d, h, fluid_data, mixing_rule=mixing_rule, T0=upper)
    return (upper - lower) / (2 * d)


def dv_mix_pdh(p, h, fluid_data, mixing_rule=None, T0=None):
    d = 1e-1
    upper = v_mix_ph(p, h + d, fluid_data, mixing_rule=mixing_rule, T0=T0)
    lower = v_mix_ph(p, h - d, fluid_data, mixing_rule=mixing_rule, T0=upper)
    return (upper - lower) / (2 * d)


def v_mix_pT(p, T, fluid_data, mixing_rule=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return 1 / pure_fluid["wrapper"].d_pT(p, T)
    else:
        _check_mixing_rule(mixing_rule, V_MIX_PT_DIRECT, "specific volume")
        return V_MIX_PT_DIRECT[mixing_rule](p, T, fluid_data)


def viscosity_mix_ph(p, h, fluid_data, mixing_rule=None, T0=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["wrapper"].viscosity_ph(p, h)
    else:
        T = T_mix_ph(p, h , fluid_data, mixing_rule, T0)
        return viscosity_mix_pT(p, T, fluid_data, mixing_rule)


def viscosity_mix_pT(p, T, fluid_data, mixing_rule=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["wrapper"].viscosity_pT(p, T)
    else:
        _check_mixing_rule(mixing_rule, V_MIX_PT_DIRECT, "viscosity")
        return VISCOSITY_MIX_PT_DIRECT[mixing_rule](p, T, fluid_data)
