
from CoolProp.CoolProp import AbstractState
import CoolProp as CP
import numpy as np


libr = AbstractState("INCOMP", "LiBr")

libr.set_mass_fractions([0.5])
libr.update(CP.QT_INPUTS, 0, 300)
libr.p()
libr.set_mass_fractions([0.2])
libr.update(CP.QT_INPUTS, 0, 300)
libr.p()

# saturation pressure by temperature and libr mass fraction x

def psat_TX(T, X):
    libr.set_mass_fractions([X])
    libr.update(CP.QT_INPUTS, 0, T)
    return libr.p()


def Xsat_Tp(T, p):
    X_min = 0
    X_max = 0.75
    X = .3
    d = 1e-5
    while True:
        res = psat_TX(T, X) - p
        deriv = (psat_TX(T, X + d) - psat_TX(T, X - d)) / (2 * d)
        X -= res / deriv
        # print(X)

        if abs(res) < 1e-6:
            if X < X_min:
                return X_min
            elif X > X_max:
                return X_max
            break

        if X >= X_max:
            X = X_max - d
        elif X <= X_min:
            X = X_min + d

    psat_TX(T, X)
    return X


def Tsat_pX(p, X):
    return newton(psat_TX, p, X)


def newton(func, p, X):
    d = 1e-3
    # T = 300
    T_min = libr.trivial_keyed_output(CP.iT_min)
    T_max = libr.trivial_keyed_output(CP.iT_max)
    T = (T_min + T_max) / 2
    while True:
        res = func(T, X) - p
        deriv = (func(T + d, X) - func(T -d, X)) / (2 * d)
        T -= res / deriv

        if T >= T_max:
            T = T_max - d
        elif T <= T_min:
            T = T_min + d

        if abs(res) < 1e-6:
            return T

# watercycle
T_evap = 275.15
p_evap = CP.CoolProp.PropsSI("P", "Q", 0, "T", T_evap, "water")

T_cond = 309.15
p_cond = CP.CoolProp.PropsSI("P", "Q", 0, "T", T_cond, "water")

# heat source temperature
T_sol_des_out = 350

# heat sink temperature
T_sol_abs_out = 306.15

x_rich = Xsat_Tp(T_sol_abs_out, p_evap)

x_water_rich = 1 - x_rich
T_water_abs_in = T_evap
h_water_abs_in = CP.CoolProp.PropsSI("H", "Q", 1, "T", T_water_abs_in, "water")

h_sol_abs_out = libr.hmass()
s_sol_abs_out = libr.smass()

eta_s_pump = 0.9
libr.update(CP.PSmass_INPUTS, p_cond, s_sol_abs_out)
h_pump_out_s = libr.hmass()
h_pump_out = h_sol_abs_out + (h_pump_out_s - h_sol_abs_out) / eta_s_pump

x_poor = Xsat_Tp(T_sol_des_out, p_cond)
h_poor = libr.hmass()#update(CP.PT_INPUTS, p_cond, T_sol_des_out)

delta_T = 0
T_water_des_out = T_sol_des_out - delta_T
h_water_des_out = CP.CoolProp.PropsSI("H", "Q", 1, "T", T_water_des_out, "water")

x_water_poor = 1 - x_poor

m_water = 1
m_rich = m_water * (1 - x_water_poor) / (x_water_rich - x_water_poor)

m_poor = m_rich - m_water

delta_H_abs = -m_water * h_water_abs_in - m_poor * h_poor + m_rich * h_sol_abs_out
delta_H_des = m_water * h_water_des_out + m_poor * h_poor - m_rich * h_pump_out
