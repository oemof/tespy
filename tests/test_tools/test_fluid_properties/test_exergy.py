from tespy.tools.fluid_properties.wrappers import CoolPropWrapper
from tespy.tools.fluid_properties.functions import calc_physical_exergy
from tespy.tools.fluid_properties.functions import h_mix_pT
from tespy.tools.fluid_properties.functions import s_mix_pT


def test_physical_exergy_at_T_smaller_T_min():

    wrapper = CoolPropWrapper("NaK", "INCOMP")
    fluid_data = {
        "NaK": {
            "mass_fraction": 1,
            "wrapper": wrapper
        }
    }

    p = 1e5
    pamb = 1e5
    temp = 600
    h = h_mix_pT(p, temp, fluid_data)
    s = s_mix_pT(p, temp, fluid_data)
    Tamb = 100

    calc_physical_exergy(h, s, p, pamb, Tamb, fluid_data)