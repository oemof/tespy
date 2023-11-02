import logging

from tespy.components import Merge, Splitter
from tespy.tools.data_containers import ComponentProperties as dc_cp


# Fictious Energy Supply models (energy flows modelled as mass flows)
# No real use for tespy I guess

class MassFactorVCC(Splitter):

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



class MassFactorVCCWithPressureLoss(MassFactorVCC):

    @staticmethod
    def component():
        return 'mass factor vapor compression cycle using COP for converting electricity to heat and cooling (energy flows modelled using tespy mass balances)'

    def get_parameters(self):
        variables = super().get_parameters()
        variables["pr"] = dc_cp(
            min_val=0,
            deriv=self.pr_deriv,
            func=self.pr_func,
            latex=self.pr_func_doc,
            num_eq=1
        )
        return variables

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

    def get_mandatory_constraints(self):
        return {
            'mass_flow_constraints': {
                'func': self.mass_flow_func, 'deriv': self.mass_flow_deriv,
                'constant_deriv': True, 'latex': self.mass_flow_func_doc,
                'num_eq': 1},
            'fluid_constraints': {
                'func': self.fluid_func, 'deriv': self.fluid_deriv,
                'constant_deriv': True, 'latex': self.fluid_func_doc,
                'num_eq': self.num_o * self.num_nw_fluids},
            'energy_balance_constraints': {
                'func': self.energy_balance_func,
                'deriv': self.energy_balance_deriv,
                'constant_deriv': True, 'latex': self.energy_balance_func_doc,
                'num_eq': self.num_o},
        }

    def calc_parameters(self):
        super().calc_parameters()
        self.pr.val = self.outl[0].p.val_SI / self.inl[0].p.val_SI
        for i in range(self.num_i):
            if self.inl[i].p.val < self.outl[0].p.val:
                msg = (
                    f"The pressure at inlet {i + 1} is lower than the pressure "
                    f"at the outlet of component {self.label}."
                )
                logging.warning(msg)




class MassFactorLossModel(Splitter):

    @staticmethod
    def component():
        return 'mass factor loss model for splitting energy flows (modelled using tespy mass balances)'

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


class MassFactorLossModelWithPressureLoss(MassFactorLossModel):

    @staticmethod
    def component():
        return 'mass factor loss model for splitting energy flows (modelled using tespy mass balances)'

    def get_parameters(self):
        variables = super().get_parameters()
        variables["pr"] = dc_cp(
            min_val=0,
            deriv=self.pr_deriv,
            func=self.pr_func,
            latex=self.pr_func_doc,
            num_eq=1
        )
        return variables

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

    def get_mandatory_constraints(self):
        return {
            'mass_flow_constraints': {
                'func': self.mass_flow_func, 'deriv': self.mass_flow_deriv,
                'constant_deriv': True, 'latex': self.mass_flow_func_doc,
                'num_eq': 1},
            'fluid_constraints': {
                'func': self.fluid_func, 'deriv': self.fluid_deriv,
                'constant_deriv': True, 'latex': self.fluid_func_doc,
                'num_eq': self.num_o * self.num_nw_fluids},
            'energy_balance_constraints': {
                'func': self.energy_balance_func,
                'deriv': self.energy_balance_deriv,
                'constant_deriv': True, 'latex': self.energy_balance_func_doc,
                'num_eq': self.num_o},
        }

    def calc_parameters(self):
        super().calc_parameters()
        self.pr.val = self.outl[0].p.val_SI / self.inl[0].p.val_SI
        for i in range(self.num_i):
            if self.inl[i].p.val < self.outl[0].p.val:
                msg = (
                    f"The pressure at inlet {i + 1} is lower than the pressure "
                    f"at the outlet of component {self.label}."
                )
                logging.warning(msg)

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

class MassFactorVCCEnergySupply(MassFactorVCC):

    @staticmethod
    def component():
        return 'mass factor vapor compression cycle using COP for converting electricity to heat and cooling (energy flows modelled using tespy mass balances, without pressure/enthalpy constraints)'

    def get_parameters(self):
        variables = super().get_parameters()
        return variables

class MassFactorLossModelEnergySupply(MassFactorLossModel):

    @staticmethod
    def component():
        return 'mass factor loss model for splitting energy flows (modelled using tespy mass balances, without pressure/enthalpy constraints)'

    def get_parameters(self):
        variables = super().get_parameters()
        return variables
    