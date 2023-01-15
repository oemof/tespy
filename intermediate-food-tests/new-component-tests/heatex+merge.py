import logging


from tespy.components import HeatExchangerSimple, Source, Sink, Merge
from tespy.connections import Connection
from tespy.networks import Network
import numpy as np

from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp


class DiabaticSimpleHeatExchanger(HeatExchangerSimple):

    @staticmethod
    def component():
        return 'diabatic simple heat exchanger'

    def get_variables(self):
        variables = super().get_variables()
        variables["eta"] = dc_cp(min_val=1e-5, val=1, max_val=1)
        variables["Q_loss"] = dc_cp(max_val=0, val=0, is_result=True)
        variables["Q_total"] = dc_cp(is_result=True)
        variables["energy_group"] = dc_gcp(
            elements=['Q_total', 'eta'],
            num_eq=1,
            latex=self.energy_balance_func_doc,
            func=self.energy_balance2_func, deriv=self.energy_balance2_deriv
        )

        return variables

    def energy_balance2_func(self):
        r"""
        Equation for pressure drop calculation.

        Returns
        -------
        residual : float
            Residual value of equation:

            .. math::

                0 =\dot{m}_{in}\cdot\left( h_{out}-h_{in}\right) -\dot{Q}
        """
        if self.Q_total.val < 0:
            return self.inl[0].m.val_SI * (
                self.outl[0].h.val_SI - self.inl[0].h.val_SI
            ) * self.eta.val - self.Q_total.val
        else:
            return self.inl[0].m.val_SI * (
                self.outl[0].h.val_SI - self.inl[0].h.val_SI
            ) - self.Q_total.val * self.eta.val

    def energy_balance2_deriv(self, increment_filter, k):
        r"""
        Calculate partial derivatives of energy balance.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        self.jacobian[k, 0, 0] = (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI)
        self.jacobian[k, 0, 2] = -self.inl[0].m.val_SI
        self.jacobian[k, 1, 2] = self.inl[0].m.val_SI
        # custom variable Q
        if self.Q_total.is_var:
            if self.Q_total.val < 0:
                self.jacobian[k, 2 + self.Q.var_pos, 0] = -1
            else:
                self.jacobian[k, 2 + self.Q.var_pos, 0] = -self.eta.val

        if self.eta.is_var:
            if self.Q_total.val < 0:
                self.jacobian[k, 2 + self.eta.var_pos, 0] = self.inl[0].m.val_SI * (
                self.outl[0].h.val_SI - self.inl[0].h.val_SI
                )
            else:
                self.jacobian[k, 2 + self.eta.var_pos, 0] = -self.Q_total.val

    def calc_parameters(self):
        super().calc_parameters()

        if self.eta.is_set:
            if self.Q.val < 0:
                self.Q_loss.val = self.Q.val * (1 - self.eta.val)
            else:
                self.Q_loss.val = -self.Q.val * (1 / self.eta.val - 1)

            self.Q_total.val = self.Q.val - self.Q_loss.val


class MergeWithPressureLoss(Merge):

    @staticmethod
    def component():
        return 'merge with pressure losses'

    def get_variables(self):
        variables = super().get_variables()
        variables["pr"] = dc_cp(
            min_val=0,
            deriv=self.pr_deriv,
            func=self.pr_func,
            latex=self.pr_func_doc,
            num_eq=1
        )
        return variables

    def get_mandatory_constraints(self):
        return {
            'mass_flow_constraints': {
                'func': self.mass_flow_func, 'deriv': self.mass_flow_deriv,
                'constant_deriv': True, 'latex': self.mass_flow_func_doc,
                'num_eq': 1},
            'fluid_constraints': {
                'func': self.fluid_func, 'deriv': self.fluid_deriv,
                'constant_deriv': False, 'latex': self.fluid_func_doc,
                'num_eq': self.num_nw_fluids},
            'energy_balance_constraints': {
                'func': self.energy_balance_func,
                'deriv': self.energy_balance_deriv,
                'constant_deriv': False, 'latex': self.energy_balance_func_doc,
                'num_eq': 1}
        }

    def pr_func(self):
        r"""
        Equation for pressure drop.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = p_\mathrm{in,1} \cdot pr - p_\mathrm{out,1}
        """
        return self.inl[0].p.val_SI * self.pr.val - self.outl[0].p.val_SI

    def pr_deriv(self, increment_filter, k):
        r"""
        Calculate the partial derivatives for combustion pressure ratio.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of equation in Jacobian matrix.
        """
        self.jacobian[k, 0, 1] = self.pr.val
        self.jacobian[k, self.num_i, 1] = -1

    def calc_parameters(self):

        self.pr.val = self.outl[0].p.val_SI / self.inl[0].p.val_SI
        for i in range(self.num_i):
            if self.inl[i].p.val < self.outl[0].p.val:
                msg = (
                    f"The pressure at inlet {i + 1} is lower than the pressure "
                    f"at the outlet of component {self.label}."
                )
                logging.warning(msg)


# caution, must write "Water" (capital W) in INCOMP backend -> CoolProp bug? Intentional?
fluids = ["INCOMP::Water", "INCOMP::T66"]


nw = Network(fluids=fluids, p_unit="bar", T_unit="C")

so = Source("Source")
so2 = Source("Source2")

#  Variant 2: Q is m (h_2 - h_1), Q_total is taking efficiency into account and represents the heat transfer over system
# boundary. For heat transfer into the system: Q = Q_total * eta, for heat transfer from the system: Q_total = Q * eta
he = DiabaticSimpleHeatExchanger("Heater")
me = MergeWithPressureLoss("Merge")
si = Sink("Sink")

c1 = Connection(so, "out1", he, "in1", label="1")
c2 = Connection(he, "out1", me, "in1", label="2")
c3 = Connection(so2, "out1", me, "in2", label="3")
c4 = Connection(me, "out1", si, "in1", label="4")

nw.add_conns(c1, c2, c3, c4)

# set some generic data for starting values
c1.set_attr(m=1, p=1.2, h=0.5e5, fluid={"Water": 0.9, "T66": 0.1})
c2.set_attr(h=2.2e5)
# mix with pure water
c3.set_attr(m=0.05, p=1.1, h=0.5e5, fluid={"Water": 1, "T66": 0})

# set pressure ratios of heater and merge
he.set_attr(pr=1)
me.set_attr(pr=0.9)

nw.solve("design")

# use temperature to make it relatable
c1.set_attr(h=None, T=30)
c2.set_attr(h=None, T=50)

nw.solve("design")

# add some heat
c2.set_attr(T=None)
# efficiency is used for postprocessing here
he.set_attr(Q=1e5, eta=0.9)

nw.solve("design")
nw.print_results()

c2.set_attr(T=50)

# impose over system boundary heat transfer (cannot be lower than actual heat transfer, efficiency value cannot be > 1!)
# In this case, efficiency decreases
he.set_attr(Q=None, Q_total=1.5e5, eta="var")

nw.solve("design")
nw.print_results()

# with set efficiency, temperature cannot be set anymore
c2.set_attr(T=None)
he.set_attr(Q_total=1.5e5, eta=.5)

nw.solve("design")
nw.print_results()

# now cooling instead of heating, CoolProp or TESPy have issues with freezing temperatures, so > 0°C
c2.set_attr(T=5)
he.set_attr(Q_total=None, eta=None)

nw.solve("design")
nw.print_results()

he.set_attr(Q_total=-.6e5, eta="var")

nw.solve("design")
nw.print_results()
