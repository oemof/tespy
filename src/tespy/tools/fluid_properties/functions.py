from .helpers import _check_mixing_rule
from .helpers import get_pure_fluid
from .helpers import get_number_of_fluids
from .helpers import inverse_temperature_mixture
from .mixtures import T_MIX_PH_REVERSE
from .mixtures import T_MIX_PS_REVERSE
from .mixtures import H_MIX_PT_DIRECT
from .mixtures import S_MIX_PT_DIRECT
from .mixtures import V_MIX_PT_DIRECT
from .mixtures import VISCOSITY_MIX_PT_DIRECT


def isentropic(p_1, h_1, p_2, fluid_data, mixing_rule=None, T0=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["wrapper"].isentropic(p_1, h_1, p_2)
    else:
        # print(mixing_rule, fluid_data)
        # print(T0)
        s_1 = s_mix_ph(p_1, h_1, fluid_data, mixing_rule)
        T_2 = T_mix_ps(p_2, s_1, fluid_data, mixing_rule)
        return h_mix_pT(p_2, T_2, fluid_data, mixing_rule)


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
        # print(p, h, fluid_data)
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
        # print(T)
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
