from CoolProp.CoolProp import AbstractState
import CoolProp as CP
import numpy as np

from tespy.components.component import Component


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


# x_range = np.linspace(0, 0.75, 76)
# for x in x_range:
#     psat_TX(300, x)
#     Tsat_pX(1e4, x)


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
    return X

# X = Xsat_Tp(350, 1e4)
# libr.set_mass_fractions([X])
# print(X)

# p_max = psat_TX(400, 0)
# print(p_max, Xsat_Tp(400, p_max))
# for p in np.linspace(1e4, p_max):
#     # print(p)
#     print(p, Xsat_Tp(400, p))

# libr.update(CP.PT_INPUTS, 1e4, 400)

# Tsat_pX(1e4, .3)
# libr
# T = libr.s
# libr.update(CP.QT_INPUTS, 0, T)

# print(libr.trivial_keyed_output(CP.iT_min))
# print(libr.trivial_keyed_output(CP.iT_freeze))
# print(libr.trivial_keyed_output(CP.iP))

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


from tespy.tools.fluid_properties.mixtures import xsat_pT_incomp_solution


class Absorber(Component):

    def __init__(self, label, **kwargs):
        super().__init__(label, **kwargs)

    @staticmethod
    def inlets():
        return ["in1", "in2"]

    @staticmethod
    def outlets():
        return ["out1"]

    def preprocess(self, num_vars):
        self.h2o = "water"
        super().preprocess(num_vars)

    def get_parameters(self):
        return {}

    def get_mandatory_constraints(self):
        constraints = {
            'mass_flow_constraints': {
                'func': self.mass_flow_func, 'deriv': self.mass_flow_deriv,
                'constant_deriv': True,# 'latex': self.mass_flow_func_doc,
                'num_eq': 1},
            'fluid_constraints': {
                'func': self.fluid_func, 'deriv': self.fluid_deriv,
                'constant_deriv': False,# 'latex': self.fluid_func_doc,
                'num_eq': 1},
            'pressure_constraints': {
                'func': self.pressure_equality_func,
                'deriv': self.pressure_equality_deriv,
                'constant_deriv': True,
                'latex': self.pressure_equality_func_doc,
                'num_eq': 2},
            "saturation_constraints_libr": {
                "func": self.saturated_solution_libr_func,
                "deriv": self.saturated_solution_libr_deriv,
                "constant_deriv": False,
                "num_eq": 1
            },
        }
        if "LiBr" in self.outl[0].fluid.is_var:
            constraints["saturation_constraints_water"] = {
                "func": self.saturated_solution_water_func,
                "deriv": self.saturated_solution_water_deriv,
                "constant_deriv": False,
                "num_eq": 1
            }
        return constraints

    def mass_flow_func(self):
        return self.inl[0].m.val_SI + self.inl[1].m.val_SI - self.outl[0].m.val_SI

    def mass_flow_deriv(self, k):
        for c in self.inl:
            if c.m.is_var:
                self.jacobian[k, c.m.J_col] = 1

        if self.outl[0].m.is_var:
            self.jacobian[k, self.outl[0].m.J_col] = -1

    def fluid_func(self):
        return (
            self.inl[0].m.val_SI
            + self.inl[1].m.val_SI * self.inl[1].fluid.val[self.h2o]
            - self.outl[0].m.val_SI * self.outl[0].fluid.val[self.h2o]
        )

    def fluid_deriv(self, increment_filter, k):
        if self.inl[0].m.is_var:
            self.jacobian[k, self.inl[0].m.J_col] = 1
        if self.inl[1].m.is_var:
            self.jacobian[k, self.inl[1].m.J_col] = self.inl[1].fluid.val[self.h2o]
        if self.outl[0].m.is_var:
            self.jacobian[k, self.outl[0].m.J_col] = -self.outl[0].fluid.val[self.h2o]

    def pressure_equality_func(self):
        residual = []
        for c in self.inl :
            residual += [c.p.val_SI - self.outl[0].p.val_SI]
        return residual

    def pressure_equality_deriv(self, k):
        r"""
        Calculate partial derivatives for all pressure equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the fluid equations.
        """
        for c in self.inl:
            if c.p.is_var:
                self.jacobian[k, c.p.J_col] = 1
            if self.outl[0].p.is_var:
                self.jacobian[k, self.outl[0].p.J_col] = -1
            k += 1

    def saturated_solution_water_func(self):
        outl = self.outl[0]
        return 1 - outl.fluid.val[outl.solvent] - outl.fluid.val[self.h2o]


    def saturated_solution_water_deriv(self, increment_filter, k):
        outl = self.outl[0]
        if self.h2o in outl.fluid.is_var:
            self.jacobian[k, outl.fluid.J_col[self.h2o]] = -1
        if outl.solvent in outl.fluid.is_var:
            self.jacobian[k, outl.fluid.J_col[outl.solvent]] = -1
        # pass

    def saturated_solution_libr_func(self):
        outl = self.outl[0]
        x_previous = outl.fluid.val[outl.solvent]
        T = outl.calc_T()
        x_libr = xsat_pT_incomp_solution(outl.p.val_SI, T, outl.fluid_data, solvent=outl.solvent)
        outl.fluid_data[outl.solvent]["wrapper"].AS.set_mass_fractions([x_previous])
        return x_libr - outl.fluid.val[outl.solvent]

    def saturated_solution_libr_deriv(self, increment_filter, k):
        outl = self.outl[0]
        if outl.p.is_var:
            deriv = self.numeric_deriv(self.saturated_solution_libr_func, "p", outl)
            self.jacobian[k, outl.p.J_col] = deriv
        if outl.h.is_var:
            deriv = self.numeric_deriv(self.saturated_solution_libr_func, "h", outl)
            self.jacobian[k, outl.h.J_col] = deriv
        if outl.solvent in outl.fluid.is_var:
            self.jacobian[k, outl.fluid.J_col[outl.solvent]] = -1

    @staticmethod
    def is_branch_source():
        return True

    def start_branch(self):
        outconn = self.outl[0]
        branch = {
            "connections": [outconn],
            "components": [self, outconn.target],
            "subbranches": {}
        }
        outconn.target.propagate_to_target(branch)

        return {outconn.label: branch}

    def propagate_to_target(self, branch):
        return

    def propagate_wrapper_to_target(self, branch):
        if branch["connections"][-1] == self.inl[0]:
            return

        if self in branch["components"]:
            return

        outconn = self.outl[0]
        branch["connections"] += [outconn]
        outconn.target.propagate_wrapper_to_target(branch)

    def new_constraints(self):
        constraints = {
            'mass_flow_constraints': {
                'func': self.mass_flow_func, 'deriv': self.mass_flow_deriv,
                'constant_deriv': True,# 'latex': self.mass_flow_func_doc,
                'num_eq': 1},
            'fluid_constraints': {
                'func': self.fluid_func, 'deriv': self.fluid_deriv,
                'constant_deriv': False,# 'latex': self.fluid_func_doc,
                'num_eq': 1},
            'pressure_constraints': {
                'func': self.pressure_equality_func,
                'deriv': self.pressure_equality_deriv,
                'constant_deriv': True,
                'latex': self.pressure_equality_func_doc,
                'num_eq': 2},
            "saturation_constraints_libr": {
                "func": self.saturated_solution_libr_func,
                "deriv": self.saturated_solution_libr_deriv,
                "constant_deriv": False,
                "num_eq": 1
            },
        }
        if "LiBr" in self.outl[0].fluid.is_var:
            constraints["saturation_constraints_water"] = {
                "func": self.saturated_solution_water_func,
                "deriv": self.saturated_solution_water_deriv,
                "constant_deriv": False,
                "num_eq": 1
            }
        return constraints


from tespy.components import Source, Sink
from tespy.networks import Network
from tespy.connections import Connection


nw = Network()

water = Source("water source")
rich = Sink("rich solution")
poor = Source("poor solution")

absorber = Absorber("absorber")

c1 = Connection(water, "out1", absorber, "in1", label="1")
c2 = Connection(poor, "out1", absorber, "in2", label="2")
c3 = Connection(absorber, "out1", rich, "in1", label="3")

nw.add_conns(c1, c2, c3)

c1.set_attr(fluid={"water": 1}, p=p_evap, m=1, x=1)
c2.set_attr(fluid={"water": x_water_poor}, h=h_poor, mixing_rule="incomp-solution", solvent="LiBr")
c3.set_attr(fluid={"INCOMP::LiBr": x_rich}, h=h_sol_abs_out, p0=p_evap)
nw.set_attr()
nw.solve("design")

# how does the reference temperature affect the energy balance
p_ref = 1e5
for T_ref in range(274, 400):
    h_water_ref = CP.CoolProp.PropsSI("H", "P", p_ref, "T", T_ref, "water")
    h_poor_ref = CP.CoolProp.PropsSI("H", "P", p_ref, "T", T_ref, f"INCOMP::LiBr[{x_poor}]")
    h_rich_ref = CP.CoolProp.PropsSI("H", "P", p_ref, "T", T_ref, f"INCOMP::LiBr[{x_rich}]")
    print((c1.h.val_SI - h_water_ref) * c1.m.val_SI + (c2.h.val_SI - h_poor_ref) * c2.m.val_SI - (c3.h.val_SI - h_rich_ref) * c3.m.val_SI)

# some checks
assert round(c1.m.val_SI * c1.fluid.val["water"] + c2.m.val_SI * c2.fluid.val["water"], 4) == round(c3.m.val_SI * c3.fluid.val["water"], 4)
assert round(c2.m.val_SI * c2.fluid.val["LiBr"], 4) == round(c3.m.val_SI * c3.fluid.val["LiBr"], 4)

# replace initial draft of constraints -> maybe add saturation constraint as parameter
absorber.get_mandatory_constraints = absorber.new_constraints
# now there is two additional equations in case both fluids at outlet 1 are variable
# and one additional equation if none of the fluids is variale
# fluid composition as function of saturation temperature and pressure
c3.fluid.is_set = set()
c1.set_attr(p=None)
c3.set_attr(h=None, T=306.15, p=p_evap)
c3.h.val0 = c3.h.val_SI
nw.solve("design")
# check results, should be identical to previous ones
assert round(c1.m.val_SI * c1.fluid.val["water"] + c2.m.val_SI * c2.fluid.val["water"], 4) == round(c3.m.val_SI * c3.fluid.val["water"], 4)
assert round(c2.m.val_SI * c2.fluid.val["LiBr"], 4) == round(c3.m.val_SI * c3.fluid.val["LiBr"], 4)

# fix the either the water or the libr mass fraction and the temperature or pressure
# gives us the other respective value, pressure in this case
c3.set_attr(fluid={"INCOMP::LiBr": 0.45}, T=295, p=None)
nw.solve("design")
assert round(c1.m.val_SI * c1.fluid.val["water"] + c2.m.val_SI * c2.fluid.val["water"], 4) == round(c3.m.val_SI * c3.fluid.val["water"], 4)
assert round(c2.m.val_SI * c2.fluid.val["LiBr"], 4) == round(c3.m.val_SI * c3.fluid.val["LiBr"], 4)

print(f"T3: {c3.T.val_SI}, p3: {c3.p.val_SI}, fluid3: {c3.fluid.val}")
