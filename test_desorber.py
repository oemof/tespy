from tespy.tools.data_containers import ComponentProperties as dc_cp
from sorption import Sorption


class Desorber(Sorption):

    @staticmethod
    def inlets():
        return ["in1"]

    @staticmethod
    def outlets():
        return ["out1", "out2"]

    def get_parameters(self):
        return {
            'Q': dc_cp(
                min_val=0,
                func=self.heat_func,
                num_eq=1,
                deriv=self.heat_deriv
            )
        }

    def pressure_equality_func(self):
        residual = []
        for c in self.outl:
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
        branches = {}
        for outconn in self.outl:
            branch = {
                "connections": [outconn],
                "components": [self]
            }
            outconn.target.propagate_wrapper_to_target(branch)
            branches[outconn.label] = branch

        return branches

    def propagate_wrapper_to_target(self, branch):
        if self in branch["components"]:
            return

        outconn = self.outl[0]
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
    c10 = Connection(desorber, "out2", water, "in1", label="10")
    c6 = Connection(desorber, "out1", valve, "in1", label="6")
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
