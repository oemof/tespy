from tespy.tools.fluid_properties.functions import calc_physical_exergy
from tespy.tools.fluid_properties.functions import h_mix_pT
from tespy.tools.fluid_properties.functions import s_mix_pT
from tespy.tools.fluid_properties.wrappers import CoolPropWrapper


def test_physical_exergy_at_T_smaller_T_min():

    wrapper = CoolPropWrapper("NaK", "INCOMP")
    fluid_data = {
        "NaK": {
            "mass_fraction": 1,
            "wrapper": wrapper
        }
    }

    p = 2e5
    pamb = 1e5
    temp = 600
    h = h_mix_pT(p, temp, fluid_data)
    s = s_mix_pT(p, temp, fluid_data)
    Tamb = 100

    ex_th, ex_mech = calc_physical_exergy(h, s, p, pamb, Tamb, fluid_data)
    # ex_mech == 0 is usually only the case at p = pamb
    # here all the exergy is assigend to the thermal one due to the minimum
    # temperature violation
    assert ex_mech == 0


def test_physical_exergy_at_T_larger_T_min():

    wrapper = CoolPropWrapper("NaK", "INCOMP")
    fluid_data = {
        "NaK": {
            "mass_fraction": 1,
            "wrapper": wrapper
        }
    }

    p = 2e5
    pamb = 1e5
    temp = 800
    h = h_mix_pT(p, temp, fluid_data)
    s = s_mix_pT(p, temp, fluid_data)
    Tamb = 600

    ex_th, ex_mech = calc_physical_exergy(h, s, p, pamb, Tamb, fluid_data)
    # ex_mech == 0 is usually only the case at p = pamb
    # here all the exergy is assigend to the thermal one due to the minimum
    # temperature violation
    assert ex_mech > 0
    assert ex_th > 0