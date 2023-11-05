
from tespy.components.component import Component

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
        if self.outl[0].solvent in self.outl[0].fluid.is_var:
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
            self.inl[1].m.val_SI * self.inl[1].fluid.val["LiBr"]
            - self.outl[0].m.val_SI * self.outl[0].fluid.val["LiBr"]
        )

    def fluid_deriv(self, increment_filter, k):
        outl = self.outl[0]
        inl = self.inl[1]

        if inl.m.is_var:
            self.jacobian[k, inl.m.J_col] = inl.fluid.val[inl.solvent]
        if inl.solvent in inl.fluid.is_var:
            self.jacobian[k, inl.fluid.J_col[inl.solvent]] = inl.m.val_SI
        if outl.m.is_var:
            self.jacobian[k, outl.m.J_col] = -outl.fluid.val[outl.solvent]
        if outl.solvent in outl.fluid.is_var:
            self.jacobian[k, outl.fluid.J_col[outl.solvent]] = -outl.m.val_SI

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

    def saturated_solution_libr_func(self):
        outl = self.outl[0]
        x_previous = outl.fluid.val[outl.solvent]
        T = outl.calc_T()
        x_libr = xsat_pT_incomp_solution(outl.p.val_SI, T, outl.fluid_data, solvent=outl.solvent, x0=x_previous)
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
            self.jacobian[k, outl.fluid.J_col[outl.solvent]] = self.numeric_deriv(self.saturated_solution_libr_func, outl.solvent, outl)

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
            branch["components"] += [self]
            return

        if self in branch["components"]:
            return

        outconn = self.outl[0]
        branch["connections"] += [outconn]
        branch["components"] += [self]
        outconn.target.propagate_wrapper_to_target(branch)


if __name__ == "__main__":

    from tespy.components import Source, Sink
    from tespy.networks import Network
    from tespy.connections import Connection

    # this runs the references calculations and imports all variables from there
    from sorption_reference import *

    nw = Network()

    water = Source("water source")
    rich = Sink("rich solution")
    poor = Source("poor solution")

    absorber = Absorber("absorber")

    c1 = Connection(water, "out1", absorber, "in1", label="1")
    c2 = Connection(poor, "out1", absorber, "in2", label="2")
    c3 = Connection(absorber, "out1", rich, "in1", label="3")

    nw.add_conns(c1, c2, c3)

    c1.set_attr(fluid={"water": 1}, p=0.01e5, m=1, x=1)
    c2.set_attr(fluid={"water": x_water_poor, "INCOMP::LiBr": 1 - x_water_poor}, h=h_poor, mixing_rule="incomp-solution", solvent="LiBr")
    c3.set_attr(h=h_sol_abs_out, p0=0.01e5, m0=m_rich)

    nw.solve("design")
    nw.print_results()
    print(c3.fluid.val)

    c2.set_attr(m=10)
    c3.set_attr(h=None)

    nw.solve("design")
    nw.print_results()

    print(c1.m.val_SI * c1.h.val_SI + c2.m.val_SI * c2.h.val_SI - c3.m.val_SI * c3.h.val_SI)

    # how does the reference temperature affect the energy balance
    # eb_prev = 0
    # p_ref = 1e5
    # for T_ref in range(274, 350):
    #     h_water_ref = CP.CoolProp.PropsSI("H", "P", p_ref, "T", T_ref, "water")
    #     h_poor_ref = CP.CoolProp.PropsSI("H", "P", p_ref, "T", T_ref, f"INCOMP::LiBr[{c2.fluid.val['LiBr']}]")
    #     h_rich_ref = CP.CoolProp.PropsSI("H", "P", p_ref, "T", T_ref, f"INCOMP::LiBr[{c3.fluid.val['LiBr']}]")
    #     (CP.CoolProp.PropsSI("Cpmass", "P", p_ref, "T", T_ref, "water"))
    #     eb = (c1.h.val_SI - h_water_ref) * c1.m.val_SI + (c2.h.val_SI - h_poor_ref) * c2.m.val_SI - (c3.h.val_SI - h_rich_ref) * c3.m.val_SI
    #     eb_prev = eb

    # s

    # some checks
    # nw.print_results()
    print(m_rich, c3.m.val_SI)
    assert round(c1.m.val_SI * c1.fluid.val["water"] + c2.m.val_SI * c2.fluid.val["water"], 4) == round(c3.m.val_SI * c3.fluid.val["water"], 4)
    assert round(c2.m.val_SI * c2.fluid.val["LiBr"], 4) == round(c3.m.val_SI * c3.fluid.val["LiBr"], 4)

    # use temperature instead of enthalpy
    c2.set_attr(m=None)
    c3.set_attr(T=306.15)
    nw.solve("design")
    nw.print_results()
    # check results, should be identical to previous ones
    assert round(c1.m.val_SI * c1.fluid.val["water"] + c2.m.val_SI * c2.fluid.val["water"], 4) == round(c3.m.val_SI * c3.fluid.val["water"], 4)
    assert round(c2.m.val_SI * c2.fluid.val["LiBr"], 4) == round(c3.m.val_SI * c3.fluid.val["LiBr"], 4)

    # fix the either the water or the libr mass fraction and the temperature or pressure
    # gives us the other respective value, pressure in this case
    c1.set_attr(p=None)
    c3.set_attr(fluid={"INCOMP::LiBr": 0.50}, T=305)
    nw.solve("design")
    assert round(c1.m.val_SI * c1.fluid.val["water"] + c2.m.val_SI * c2.fluid.val["water"], 4) == round(c3.m.val_SI * c3.fluid.val["water"], 4)
    assert round(c2.m.val_SI * c2.fluid.val["LiBr"], 4) == round(c3.m.val_SI * c3.fluid.val["LiBr"], 4)

    print(f"T3: {c3.T.val_SI}, p3: {c3.p.val_SI}, fluid3: {c3.fluid.val}")

    nw.print_results()
