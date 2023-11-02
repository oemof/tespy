import logging

from tespy.components import Merge, Splitter
from tespy.tools.data_containers import ComponentProperties as dc_cp


# Fictious Energy Supply models (energy flows modelled as mass flows)
# No real use for tespy I guess

class MassFactorEnergySupply(Splitter):

    @staticmethod
    def component():
        return 'mass factor vapor compression cycle using COP for converting electricity to heat and cooling (energy flows modelled using tespy mass balances)'

    def get_parameters(self):
        variables = super().get_parameters()
        variables["COP"] = dc_cp(
            min_val=0,
            deriv=self.COP_deriv,
            func=self.COP_func,
            latex=self.mass_flow_func_doc,
            num_eq=1
        )
        return variables

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        del constraints['pressure_constraints']
        del constraints['energy_balance_constraints']
        #del constraints['mass_flow_constraints']
        return constraints   

    def COP_func(self):
        r"""
        Equation for COP.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = p_\mathrm{in,1} \cdot pr - p_\mathrm{out,1}
        """
        return self.inl[0].m.val_SI * self.COP.val - self.outl[0].m.val_SI

    def COP_deriv(self, increment_filter, k):
        r"""
        Calculate the partial derivatives for combustion pressure ratio.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of equation in Jacobian matrix.
        """
        inl = self.inl[0]
        outl = self.outl[0]
        if inl.m.is_var:
            self.jacobian[k, inl.m.J_col] = self.COP.val
        if outl.m.is_var:
            self.jacobian[k, outl.m.J_col] = -1

    def calc_parameters(self):
        super().calc_parameters()
        self.COP.val = self.outl[0].m.val_SI / (self.outl[0].m.val_SI - (-self.outl[1].m.val_SI))


class MassLossEnergySupply(Splitter):

    @staticmethod
    def component():
        return 'mass loss model for splitting energy flows (modelled using tespy mass balances)'

    def get_parameters(self):
        variables = super().get_parameters()
        variables["Loss"] = dc_cp(
            min_val=0,
            deriv=self.Loss_deriv,
            func=self.Loss_func,
            latex=self.mass_flow_func_doc,
            num_eq=1
        )
        return variables

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        del constraints['pressure_constraints']
        del constraints['energy_balance_constraints']
        return constraints   

    def Loss_func(self):
        return self.inl[0].m.val_SI * (1-self.Loss.val) - self.outl[0].m.val_SI

    def Loss_deriv(self, increment_filter, k):
        inl = self.inl[0]
        outl = self.outl[0]
        if inl.m.is_var:
            self.jacobian[k, inl.m.J_col] = (1-self.Loss.val)
        if outl.m.is_var:
            self.jacobian[k, outl.m.J_col] = -1

    def calc_parameters(self):
        super().calc_parameters()
        self.Loss.val = (self.inl[0].m.val_SI - self.outl[0].m.val_SI)/self.inl[0].m.val_SI

class MergeEnergySupply(Merge):

    @staticmethod
    def component():
        return 'merge without pressure/energy constraints'

    def get_parameters(self):
        variables = super().get_parameters()
        return variables

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        del constraints['pressure_constraints']
        del constraints['energy_balance_constraints']
        return constraints

class SplitterEnergySupply(Splitter):

    @staticmethod
    def component():
        return 'Splitter without pressure/energy constraints'

    def get_parameters(self):
        variables = super().get_parameters()
        return variables

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        del constraints['pressure_constraints']
        del constraints['energy_balance_constraints']
        return constraints   