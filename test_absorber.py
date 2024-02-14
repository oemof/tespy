
from tespy.tools.data_containers import ComponentProperties as dc_cp
from sorption import Sorption


class Absorber(Sorption):

    @staticmethod
    def inlets():
        return ["in1", "in2"]

    @staticmethod
    def outlets():
        return ["out1"]

    def get_parameters(self):
        return {
            'Q': dc_cp(
                max_val=0,
                func=self.heat_func,
                num_eq=1,
                deriv=self.heat_deriv
            )
        }

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

    def propagate_wrapper_to_target(self, branch):
        if branch["connections"][-1] == self.inl[1]:
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

    c1 = Connection(water, "out1", absorber, "in2", label="1")
    c2 = Connection(poor, "out1", absorber, "in1", label="2")
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
