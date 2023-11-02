
from tespy.components.component import Component

from tespy.tools.fluid_properties.mixtures import xsat_pT_incomp_solution


class Desorber(Component):

    def __init__(self, label, **kwargs):
        super().__init__(label, **kwargs)

    @staticmethod
    def inlets():
        return ["in1"]

    @staticmethod
    def outlets():
        return ["out1", "out2"]

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
        if self.outl[1].solvent in self.outl[1].fluid.is_var:
            constraints["saturation_constraints_water"] = {
                "func": self.saturated_solution_water_func,
                "deriv": self.saturated_solution_water_deriv,
                "constant_deriv": False,
                "num_eq": 1
            }
        return constraints

    def mass_flow_func(self):
        return self.inl[0].m.val_SI - self.outl[0].m.val_SI - self.outl[1].m.val_SI

    def mass_flow_deriv(self, k):
        inl = self.inl[0]
        if inl.m.is_var:
            self.jacobian[k, inl.m.J_col] = 1

        for outl in self.outl:
            if outl.m.is_var:
                self.jacobian[k, outl.m.J_col] = -1

    def fluid_func(self):
        inl = self.inl[0]
        outl = self.outl[1]
        return (
            inl.m.val_SI * inl.fluid.val[inl.solvent]
            - outl.m.val_SI * outl.fluid.val[outl.solvent]
        )

    def fluid_deriv(self, increment_filter, k):
        outl = self.outl[1]
        inl = self.inl[0]

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
        for c in self.outl :
            residual += [c.p.val_SI - self.inl[0].p.val_SI]
        return residual

    def pressure_equality_deriv(self, k):
        r"""
        Calculate partial derivatives for all pressure equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the fluid equations.
        """
        for c in self.outl:
            if c.p.is_var:
                self.jacobian[k, c.p.J_col] = 1
            if self.inl[0].p.is_var:
                self.jacobian[k, self.inl[0].p.J_col] = -1
            k += 1

    def saturated_solution_water_func(self):
        outl = self.outl[1]
        return 1 - outl.fluid.val[outl.solvent] - outl.fluid.val[self.h2o]

    def saturated_solution_water_deriv(self, increment_filter, k):
        outl = self.outl[1]
        if self.h2o in outl.fluid.is_var:
            self.jacobian[k, outl.fluid.J_col[self.h2o]] = -1
        if outl.solvent in outl.fluid.is_var:
            self.jacobian[k, outl.fluid.J_col[outl.solvent]] = -1

    def saturated_solution_libr_func(self):
        outl = self.outl[1]
        x_previous = outl.fluid.val[outl.solvent]
        T = outl.calc_T()
        x_libr = xsat_pT_incomp_solution(outl.p.val_SI, T, outl.fluid_data, solvent=outl.solvent, x0=x_previous)
        outl.fluid_data[outl.solvent]["wrapper"].AS.set_mass_fractions([x_previous])
        return x_libr - outl.fluid.val[outl.solvent]

    def saturated_solution_libr_deriv(self, increment_filter, k):
        outl = self.outl[1]
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

    @staticmethod
    def is_wrapper_branch_source():
        return True

    def start_branch(self):
        branches = {}
        for outconn in self.outl:
            branch = {
                "connections": [outconn],
                "components": [self, outconn.target],
                "subbranches": {}
            }
            outconn.target.propagate_to_target(branch)

            branches[outconn.label] = branch
        return branches

    def start_fluid_wrapper_branch(self):
        outconn = self.outl[0]
        branch = {
            "connections": [outconn],
            "components": [self]
        }
        outconn.target.propagate_wrapper_to_target(branch)

        return {outconn.label: branch}

    def propagate_to_target(self, branch):
        return

    def propagate_wrapper_to_target(self, branch):
        if self in branch["components"]:
            return

        outconn = self.outl[1]
        branch["connections"] += [outconn]
        branch["components"] += [self]
        outconn.target.propagate_wrapper_to_target(branch)


if __name__ == "__main__":

    from tespy.components import Source, Sink, SimpleHeatExchanger
    from tespy.networks import Network
    from tespy.connections import Connection

    # this runs the references calculations and imports all variables from there
    from sorption_reference import *


    nw = Network()

    water = Sink("water source")
    rich = Source("rich solution")
    poor = Sink("poor solution")

    desorber = Desorber("desorber")
    valve = SimpleHeatExchanger("valve")

    c5 = Connection(rich, "out1", desorber, "in1", label="5")
    c10 = Connection(desorber, "out1", water, "in1", label="10")
    c6 = Connection(desorber, "out2", valve, "in1", label="6")
    c7 = Connection(valve, "out1", poor, "in1", label="7")

    nw.add_conns(c5, c10, c6, c7)

    c5.set_attr(fluid={"water": x_water_rich, "INCOMP::LiBr": 1 - x_water_rich}, p=p_cond, m=10, h=h_pump_out, mixing_rule="incomp-solution", solvent="LiBr")
    c10.set_attr(fluid={"water": 1}, x=1)
    c6.set_attr(h=h_poor, p0=0.1e5, m0=m_rich)
    c7.set_attr(p=p_cond, h=h_poor * 1.2)

    nw.solve("design")
    print(c6.fluid.val)
    print(c10.h.val_SI * c10.m.val_SI + c6.h.val_SI * c6.m.val_SI - c5.m.val_SI * c5.h.val_SI)
    print(c7.h.val_SI)
    print(c7.T.val_SI)

    c10.set_attr(x=None, T=320)
    nw.solve("design")
    print(c6.fluid.val)
    print(c10.h.val_SI * c10.m.val_SI + c6.h.val_SI * c6.m.val_SI - c5.m.val_SI * c5.h.val_SI)
    exit()
    nw.print_results()

    print(T_sol_des_out)
    print(h_poor)
    print(c6.h.val_SI)
    print(c6.T.val_SI)
    print(c10.T.val_SI)
    print(c6.fluid.val)
    c6.set_attr(m=9)
    c6.set_attr(h=None)

    nw.solve("design")
    nw.print_results()
    print(c6.fluid.val)

    c6.set_attr(m=None, T=350)
    nw.solve("design")
    nw.print_results()
    print(c6.fluid.val)

    c6.set_attr(fluid={"INCOMP::LiBr": 0.7})
    c5.set_attr(p=None)

    nw.solve("design")
    nw.print_results()
    print(c6.fluid.val)
