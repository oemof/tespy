from .helpers import get_pure_fluid
from .helpers import get_number_of_fluids
from .helpers import inverse_temperature_mixture
from .mixtures import T_MIX_PH_REVERSE
from .mixtures import T_MIX_PS_REVERSE
from .mixtures import H_MIX_PT_DIRECT
from .mixtures import S_MIX_PT_DIRECT
from .mixtures import V_MIX_PT_DIRECT
from .mixtures import VISCOSITY_MIX_PT_DIRECT


def T_mix_ph(p, h, fluid_data, mixing_rule=None, T0=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["property_object"].T_ph(p, h)
    else:
        kwargs = {
            "p": p, "target_value": h, "fluid_data": fluid_data, "T0": T0,
            "f": T_MIX_PH_REVERSE[mixing_rule]
        }
        return inverse_temperature_mixture(**kwargs)


def h_mix_pT(p, T, fluid_data, mixing_rule=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["property_object"].h_pT(p, T)
    else:
        return H_MIX_PT_DIRECT[mixing_rule](p, T, fluid_data)


def s_mix_pT(p, T, fluid_data, mixing_rule=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["property_object"].s_pT(p, T)
    else:
        return S_MIX_PT_DIRECT[mixing_rule](p, T, fluid_data)


def T_mix_ps(p, s, fluid_data, mixing_rule=None, T0=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["property_object"].T_ps(p, s)
    else:
        kwargs = {
            "p": p, "target_value": s, "fluid_data": fluid_data, "T0": T0,
            "f": T_MIX_PS_REVERSE[mixing_rule]
        }
        return inverse_temperature_mixture(**kwargs)


def v_mix_ph(p, h, fluid_data, mixing_rule=None, T0=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return 1 / pure_fluid["property_object"].d_ph(p, h)
    else:
        T = T_mix_ph(p, h , fluid_data, mixing_rule, T0)
        return v_mix_pT(p, T, fluid_data, mixing_rule)


def v_mix_pT(p, T, fluid_data, mixing_rule=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return 1 / pure_fluid["property_object"].d_pT(p, T)
    else:
        return V_MIX_PT_DIRECT[mixing_rule](p, T, fluid_data)


def viscosity_mix_ph(p, h, fluid_data, mixing_rule=None, T0=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["property_object"].viscosity_ph(p, h)
    else:
        T = T_mix_ph(p, h , fluid_data, mixing_rule, T0)
        return viscosity_mix_pT(p, T, fluid_data, mixing_rule)


def viscosity_mix_pT(p, T, fluid_data, mixing_rule=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["property_object"].viscosity_pT(p, T)
    else:
        return VISCOSITY_MIX_PT_DIRECT[mixing_rule](p, T, fluid_data)
