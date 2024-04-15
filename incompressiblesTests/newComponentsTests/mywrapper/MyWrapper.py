import CoolProp.CoolProp as CP
from tespy.tools.fluid_properties.wrappers import FluidPropertyWrapper, CoolPropWrapper
import numpy as np
import matplotlib.pyplot as plt

# # coefficients   a      b       c    d        
# COEF = {
#    "protein": {
#        "unit" : "C",
#        "cp": [2008.2,     1.2089, -1.3129*1e-3,    0.0],
#        "d" : [1329.9,    -0.5184,          0.0,    0.0],
#     }
# }

class MyWrapper(FluidPropertyWrapper):
    def __init__(self, fluid, back_end=None, Tref = 293.15, coefs=[]) -> None:
        super().__init__(fluid, back_end)
        if self.fluid not in coefs:
            msg = "Fluid not available in database"
            raise KeyError(msg)

        # get coefs (converted to kelvin) and calculate reference
        self.T0 = Tref
        self.get_coefs(coefs)
        
        self._molar_mass = 1
        self._T_min = 100
        self._T_max = 2000
        self._p_min = 1000
        self._p_max = 10000000

    def get_coefs(self, coefs):
        if coefs[self.fluid]["unit"] == "C":
            self.C_c = coefs[self.fluid]["cp"]
            self.C_d = coefs[self.fluid]["d"]
            # convert coefficients
            T_C = np.linspace(1,50)
            cp = self.cp_pT(None,T_C)
            d = self.d_pT(None,T_C)
            T_K = np.linspace(1+273.15,50+273.15)
            self.C_c = list(np.polyfit(T_K, cp, len(coefs[self.fluid]["cp"])-1))
            self.C_c = self.C_c[::-1]
            self.C_d = list(np.polyfit(T_K, d, len(coefs[self.fluid]["d"])-1))
            self.C_d = self.C_d[::-1]
        elif coefs[self.fluid]["unit"] == "K":
            self.C_c = coefs[self.fluid]["cp"]
            self.C_d = coefs[self.fluid]["d"]
        else:
            ValueError("unit is not C or K")

    def cp_pT(self, p, T):
        return np.sum([self.C_c[i] * T**i for i in range(len(self.C_c))], axis=0)
   
    def d_pT(self, p, T, **kwargs):
        return np.sum([self.C_d[i] * T**i for i in range(len(self.C_d))], axis=0)

    def u_pT(self, p, T):
        integral = 0
        for i in range(len(self.C_c)):
            integral += (1 / (i + 1)) * self.C_c[i] * (T**(i + 1) - self.T0**(i + 1))
        return integral 

    def h_pT(self, p, T, **kwargs):
        u = self.u_pT(p, T)
        d = self.d_pT(p, T)
        return u - p/d

    def s_pT(self, p, T, **kwargs):
        integral = self.C_c[0] * np.log(T / self.T0)
        for i in range(len(self.C_c) - 1):
            integral += (1 / (i + 1)) * self.C_c[i + 1] * (T**(i + 1) - self.T0**(i + 1))
        return integral 

    def T_ph(self, p, h):
        return self.newton(self.h_pT, self.cp_pT, h, p)   

    def T_ps(self, p, s):
        return self.newton(self.s_pT, self.dsdT, s, p)
    def dsdT(self, p, T):
        return self.cp_pT(p, T)/T
    def h_ps(self, p, s):
        T = self.T_ps(p, s)
        return self.h_pT(p, T)

    def s_ph(self, p, h):
        T = self.T_ph(p, h)
        return self.s_pT(p, T)

    def isentropic(self, p_1, h_1, p_2):
        return self.h_ps(p_2, self.s_ph(p_1, h_1))

    def newton(self, func, deriv, val, p):
        # default valaues
        T = 300
        valmin = -1000
        valmax = 3000
        max_iter = 10
        tol_rel = 1e-6
        # start newton loop
        expr = True
        i = 0
        while expr:
            # calculate function residual and new value
            res = val - func(p, T)
            T += res / deriv(p, T)
            # check for value ranges
            if T < valmin:
                T = valmin
            if T > valmax:
                T = valmax
            i += 1
            if i > max_iter:
                break
            expr = abs(res / val) >= tol_rel
        return T    



if __name__ == "__main__":

    fluidwrap = MyWrapper("protein") 

    T = 300
    p = 1e5 
    u = fluidwrap.u_pT(p, T)
    d = fluidwrap.d_pT(p, T)
    h = fluidwrap.h_pT(p, T)
    s = fluidwrap.s_pT(p, T)
    print(f"u = {u}    d = {d}    h = {h}    s = {s}")

    T = 273.15
    u = fluidwrap.u_pT(p, T)
    d = fluidwrap.d_pT(p, T)
    h = fluidwrap.h_pT(p, T)
    s = fluidwrap.s_pT(p, T)
    print(f"u = {u}    d = {d}    h = {h}    s = {s}")

    T = 373.15
    u = fluidwrap.u_pT(p, T)
    d = fluidwrap.d_pT(p, T)
    h = fluidwrap.h_pT(p, T)
    s = fluidwrap.s_pT(p, T)
    print(f"u = {u}    d = {d}    h = {h}    s = {s}")
    T = fluidwrap.T_ph(p,h)
    s = fluidwrap.s_ph(p,h)
    print(f"recalc: T = {T}    s = {s}")
    T = fluidwrap.T_ps(p,s)
    h = fluidwrap.h_ps(p,s)
    print(f"recalc: T = {T}    h = {h}")

    CP_cp = []
    CP_k  = []
    CP_d  = []
    CP_h  = []
    CP_s  = []

    wrap_cp = []
    wrap_d = []
    wrap_h = []
    wrap_s = []

    p = 101325 * 5

    #Specific heat, kJ/(kg·K) 
    Tplt = np.linspace(273.15,373.15)
    for T in Tplt:
        CP_cp      += [CP.PropsSI('C','T',T,'P',p,'INCOMP::FoodProtein')]
        CP_k       += [CP.PropsSI('L','T',T,'P',p,'INCOMP::FoodProtein')]
        CP_d       += [CP.PropsSI('D','T',T,'P',p,'INCOMP::FoodProtein')]
        CP_h       += [CP.PropsSI('H','T',T,'P',p,'INCOMP::FoodProtein')]
        CP_s       += [CP.PropsSI('S','T',T,'P',p,'INCOMP::FoodProtein')]
        wrap_cp    += [fluidwrap.cp_pT(p, T)]
        wrap_d     += [fluidwrap.d_pT(p, T)]
        wrap_h    += [fluidwrap.h_pT(p, T)]
        wrap_s     += [fluidwrap.s_pT(p, T)]



    fig, ax = plt.subplots(2, 2, figsize=(16, 8))
    ax = ax.flatten()
    [a.grid() for a in ax]
    [a.set_xlabel('temperature, K') for a in ax]

    ax[0].plot(Tplt, wrap_cp)
    ax[1].plot(Tplt, wrap_d)
    ax[2].plot(Tplt, wrap_h)
    ax[3].plot(Tplt, wrap_s)

    ax[0].scatter(Tplt, CP_cp)
    ax[1].scatter(Tplt, CP_d)
    ax[2].scatter(Tplt, CP_h)
    ax[3].scatter(Tplt, CP_s)

    ax[0].set_ylabel('cp')
    ax[1].set_ylabel('d')
    ax[2].set_ylabel('h')
    ax[3].set_ylabel('s')

    plt.show()

    print("hey")