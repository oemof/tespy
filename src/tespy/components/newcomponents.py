import logging

from tespy.components import SimpleHeatExchanger, Merge, Separator, Splitter, HeatExchanger
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.fluid_properties import T_mix_ph

from tespy.components.component import Component


import warnings

import numpy as np

class DiabaticSimpleHeatExchanger(SimpleHeatExchanger):

    @staticmethod
    def component():
        return 'diabatic simple heat exchanger'

    def get_parameters(self):
        variables = super().get_parameters()
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


        i = self.inl[0]
        o = self.outl[0]
        if i.m.is_var:
            self.jacobian[k, i.m.J_col] = (o.h.val_SI - i.h.val_SI)
        if i.h.is_var:
            self.jacobian[k, i.h.J_col] = -i.m.val_SI
        if o.h.is_var:
            self.jacobian[k, o.h.J_col] = i.m.val_SI
        # custom variable Q
        if self.Q_total.is_var:
            if self.Q_total.val < 0:
                self.jacobian[k, self.Q_total.J_col] = -1
            else:
                self.jacobian[k, self.Q_total.J_col] = -self.eta.val

        if self.eta.is_var:
            if self.Q_total.val < 0:
                self.jacobian[k, self.eta.J_col] = self.inl[0].m.val_SI * (
                self.outl[0].h.val_SI - self.inl[0].h.val_SI
                )
            else:
                self.jacobian[k, self.eta.J_col] = -self.Q_total.val

    def calc_parameters(self):
        super().calc_parameters()

        if self.eta.is_set:
            if self.Q.val < 0:
                self.Q_loss.val = self.Q.val * (1 - self.eta.val)
            else:
                self.Q_loss.val = -self.Q.val * (1 / self.eta.val - 1)

            self.Q_total.val = self.Q.val - self.Q_loss.val


class SimpleHeatExchangerDeltaP(SimpleHeatExchanger):

    @staticmethod
    def component():
        return 'simple heat exchanger with pressure drop'

    def get_parameters(self):
        variables = super().get_parameters()
        variables["deltaP"] = dc_cp(
            min_val=0,
            deriv=self.deltaP_deriv,
            func=self.deltaP_func,
            latex=self.pr_func_doc,
            num_eq=1,
        )
        return variables

    def deltaP_func(self):
        r"""
        Equation for pressure drop.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = p_\mathrm{in,1} \cdot pr - p_\mathrm{out,1}
        """

        return self.inl[0].p.val_SI - self.deltaP.val*1e5 - self.outl[0].p.val_SI

    def deltaP_deriv(self, increment_filter, k, pr='', inconn=0, outconn=0):
        r"""
        Calculate the partial derivatives for pressure drop.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of equation in Jacobian matrix.

        pr : str
            Component parameter to evaluate the pr_func on, e.g.
            :code:`pr1`.

        inconn : int
            Connection index of inlet.

        outconn : int
            Connection index of outlet.
        """
        
        deltaP = self.get_attr("deltaP")
        i = self.inl[inconn]
        o = self.outl[inconn]
        if i.p.is_var:
            self.jacobian[k, i.p.J_col] = 1
        if o.p.is_var:
            self.jacobian[k, o.p.J_col] = -1
        if deltaP.is_var:
            self.jacobian[k, self.pr.J_col] = 1


    def calc_parameters(self):
        super().calc_parameters()
        self.deltaP.val = (self.inl[0].p.val_SI - self.outl[0].p.val_SI)/1e5


class SimpleHeatExchangerDeltaPLossFactor(SimpleHeatExchangerDeltaP):

    @staticmethod
    def component():
        return 'diabatic simple heat exchanger'

    def get_parameters(self):
        variables = super().get_parameters()
        variables["LF"] = dc_cp(min_val=0, val=0, max_val=1, is_result=True)
        variables["Q_loss"] = dc_cp(is_result=True)
        variables["Q_total"] = dc_cp(is_result=True)       
        variables["energy_group1"] = dc_gcp(
                elements=['LF', 'Q_total'],
                func=self.Q_total_func,
                deriv=self.Q_total_deriv,
                latex=self.energy_balance_func_doc, num_eq=1)
        variables["energy_group2"] = dc_gcp(
                elements=['Q_loss', 'Q_total'],
                func=self.Q_total_func,
                deriv=self.Q_total_deriv,
                latex=self.energy_balance_func_doc, num_eq=1)
        variables["energy_group3"] = dc_gcp(
                elements=['Q_loss', 'LF'],
                func=self.Q_total_func,
                deriv=self.Q_total_deriv,
                latex=self.energy_balance_func_doc, num_eq=1)                
        return variables

    def Q_total_func(self):
        r"""
        Equation for total heat flow rate

        """
        # self.Q_loss.val is negative and Q_total is positive (and vice versa)

        if self.energy_group2.is_set:
            self.LF.val = -self.Q_loss.val/(self.inl[0].m.val_SI * (self.outl[0].h.val_SI - self.inl[0].h.val_SI))
        if self.energy_group3.is_set:
            self.Q_total.val = -self.Q_loss.val*(1+self.LF.val)/self.LF.val 

        return self.inl[0].m.val_SI * (self.outl[0].h.val_SI - self.inl[0].h.val_SI)*(1+self.LF.val) - self.Q_total.val
      

    def Q_total_deriv(self, increment_filter, k):
        r"""
        Calculate partial derivatives of Q_total

        """

        i = self.inl[0]
        o = self.outl[0]
        if i.m.is_var:
            self.jacobian[k, i.m.J_col] = (o.h.val_SI - i.h.val_SI)*(1+self.LF.val)
        if i.h.is_var:
            self.jacobian[k, i.h.J_col] = -i.m.val_SI*(1+self.LF.val)
        if o.h.is_var:
            self.jacobian[k, o.h.J_col] = i.m.val_SI*(1+self.LF.val)
        if self.Q_total.is_var:
            self.jacobian[k, self.Q_total.J_col] = -1
        if self.LF.is_var: 
            self.jacobian[k, self.LF.J_col] = self.inl[0].m.val_SI * (self.outl[0].h.val_SI - self.inl[0].h.val_SI)
        if self.Q_loss.is_var:
            self.jacobian[k, self.Q_loss.J_col] = -(1+self.LF.val)/self.LF.val

    def calc_parameters(self):
        super().calc_parameters()

        # repeat calculations to ensure variables are assigned
        if self.Q_total.is_set:
            self.Q_loss.val = self.Q.val-self.Q_total.val
            self.LF.val = -self.Q_loss.val / self.Q.val
        elif self.LF.is_set:
            self.Q_total.val = self.Q.val * (1+self.LF.val)
            self.Q_loss.val = self.Q.val-self.Q_total.val
        else:
            self.Q_total.val = self.Q.val-self.Q_loss.val
            self.LF.val = -self.Q_loss.val/self.Q.val        


class SimpleHeatExchangerDeltaPLfKpi(SimpleHeatExchangerDeltaP):

    @staticmethod
    def component():
        return 'simple heat exchanger with loss factor and KPI'

    def get_parameters(self):
        variables = super().get_parameters()
        variables["LF"] = dc_cp(min_val=0, val=0, max_val=1, is_result=True)
        variables["Q_loss"] = dc_cp(is_result=True)
        variables["KPI"] = dc_cp(
            deriv=self.KPI_deriv,
            func=self.KPI_func,
            latex=self.pr_func_doc,
            num_eq=1)
        return variables

    def energy_balance_func(self):
        r"""
        Equation for total heat flow rate

        """

        return self.inl[0].m.val_SI * (self.outl[0].h.val_SI - self.inl[0].h.val_SI)*(1+self.LF.val) - self.Q.val

    def energy_balance_deriv(self, increment_filter, k):
        r"""
        Calculate partial derivatives of Q_total

        """

        i = self.inl[0]
        o = self.outl[0]
        if i.m.is_var:
            self.jacobian[k, i.m.J_col] = (o.h.val_SI - i.h.val_SI)*(1+self.LF.val) 
        if i.h.is_var:
            self.jacobian[k, i.h.J_col] = -i.m.val_SI*(1+self.LF.val) 
        if o.h.is_var:
            self.jacobian[k, o.h.J_col] = i.m.val_SI*(1+self.LF.val) 
        if self.Q.is_var:
            self.jacobian[k, self.Q.J_col] = -1
        if self.LF.is_var: 
            self.jacobian[k, self.LF.J_col] = self.inl[0].m.val_SI * (self.outl[0].h.val_SI - self.inl[0].h.val_SI)

    def KPI_func(self):
        r"""
        Equation for total heat flow rate

        """
        return self.inl[0].m.val_SI * (self.outl[0].h.val_SI - self.inl[0].h.val_SI)*(1+self.LF.val) - self.KPI.val * self.inl[0].m.val_SI 

    def KPI_deriv(self, increment_filter, k):
        r"""
        Calculate partial derivatives of Q_total

        """
        i = self.inl[0]
        o = self.outl[0]
        if i.m.is_var:
            self.jacobian[k, i.m.J_col] = (o.h.val_SI - i.h.val_SI)*(1+self.LF.val) - self.KPI.val
        if i.h.is_var:
            self.jacobian[k, i.h.J_col] = -i.m.val_SI*(1+self.LF.val) 
        if o.h.is_var:
            self.jacobian[k, o.h.J_col] = i.m.val_SI*(1+self.LF.val) 
        if self.LF.is_var: 
            self.jacobian[k, self.LF.J_col] = self.inl[0].m.val_SI * (self.outl[0].h.val_SI - self.inl[0].h.val_SI)
        if self.KPI.is_var:
            self.jacobian[k, self.Q_loss.J_col] = -self.inl[0].m.val_SI 

    def calc_parameters(self):
        super().calc_parameters()
        self.Q.val = self.inl[0].m.val_SI * (self.outl[0].h.val_SI - self.inl[0].h.val_SI)*(1+self.LF.val)
        # repeat calculations to ensure variables are assigned
        if self.KPI.is_set:
            self.Q.val = self.KPI.val * self.inl[0].m.val_SI 
        else:
            self.KPI.val = self.Q.val / self.inl[0].m.val_SI 
        self.Q_loss.val = - self.LF.val * self.Q.val


class TwoStreamHeatExchanger(HeatExchanger):

    @staticmethod
    def component():
        return 'two stream heat exchanger with min ttd (pinch)'

    def get_parameters(self):
        variables = super().get_parameters()
        variables['ttd_min'] = dc_cp(
                min_val=0, num_eq=1, func=self.ttd_min_func,
                deriv=self.ttd_min_deriv, latex=self.ttd_u_func_doc)
        return variables

    def _calc_dTs(self):
        i1 = self.inl[0]
        o1 = self.outl[0]
        i2 = self.inl[1]
        o2 = self.outl[1]

        T_i1 = i1.calc_T(T0=i1.T.val_SI)
        T_o1 = o1.calc_T(T0=o1.T.val_SI)
        T_i2 = i2.calc_T(T0=i2.T.val_SI)
        T_o2 = o2.calc_T(T0=o2.T.val_SI)

        dTa = abs(T_i1-T_o2)
        dTb = abs(T_i2-T_o1)
        return dTa,dTb

    def ttd_min_func(self):
        r"""
        Equation for minimum terminal temperature difference.
        """

        dTa,dTb = self._calc_dTs()

        if dTa < dTb:
            return self.ttd_min.val - dTa
        else:
            return self.ttd_min.val - dTb

        # T_o2 = o.calc_T(T0=o.T.val_SI)
        # return self.ttd_u.val - T_i1 + T_o2


    def ttd_min_deriv(self, increment_filter, k):
        """
        Calculate partial derivates for minimum terminal temperature difference..

        """
        f = self.ttd_min_func
        for c in [self.inl[0], self.inl[1], self.outl[0], self.outl[1]]:
            if self.is_variable(c.p): #, increment_filter): increment filter may detect no change on the wrong end 
                self.jacobian[k, c.p.J_col] = self.numeric_deriv(f, 'p', c)
            if self.is_variable(c.h): #, increment_filter):
                self.jacobian[k, c.h.J_col] = self.numeric_deriv(f, 'h', c)

    def calc_parameters(self):
        super().calc_parameters()
        if not self.ttd_min.is_set:
            self.ttd_min.val = min(self._calc_dTs())




class MergeDeltaP(Merge):

    @staticmethod
    def component():
        return 'merge with pressure losses'

    def get_parameters(self):
        variables = super().get_parameters()
        variables["deltaP"] = dc_cp(
            min_val=0,
            deriv=self.deltaP_deriv,
            func=self.deltaP_func,
            latex=self.pr_func_doc,
            num_eq=1,
        )
        return variables

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        del constraints['pressure_constraints']
        return constraints

    def deltaP_func(self):
        r"""
        Equation for pressure drop.

        """
        p_in_min = min([i.p.val_SI for i in self.inl])
        return p_in_min - self.deltaP.val*1e5 - self.outl[0].p.val_SI

    def deltaP_deriv(self, increment_filter, k):
        r"""
        Calculate the partial derivatives for pressure drop.

        """
        p_in = [i.p.val_SI for i in self.inl]
        p_min_index = p_in.index(min(p_in))

        if self.inl[p_min_index].p.is_var:
            self.jacobian[k, self.inl[p_min_index].p.J_col] = 1 #self.pr.val
        if self.outl[0].p.is_var:
            self.jacobian[k, self.outl[0].p.J_col] = -1

    def calc_parameters(self):
        super().calc_parameters()
        Pmin = min([i.p.val_SI for i in self.inl])
        Pmax = max([i.p.val_SI for i in self.inl])
        if abs(self.outl[0].p.val_SI - Pmin) >= abs(self.outl[0].p.val_SI - Pmax):
            self.deltaP.val = (Pmin - self.outl[0].p.val_SI)/1e5
        else:
            self.deltaP.val = (Pmax - self.outl[0].p.val_SI)/1e5


class SeparatorWithSpeciesSplits(Separator):

    def __init__(self, label, **kwargs):
        #self.set_attr(**kwargs)
        # need to assign the number of outlets before the variables are set
        self.num_out = 2 # default
        for key in kwargs:
            if key == 'num_out':
                self.num_out=kwargs[key]
        super().__init__(label, **kwargs)


    @staticmethod
    def component():
        return 'separator with species flow splits'

    def get_parameters(self):
        variables = super().get_parameters()
        variables["SFS"] = dc_cp_SFS(
            min_val=0,
            deriv=self.SFS_deriv,
            func=self.SFS_func,
            latex=self.pr_func_doc,
            num_eq=1,
        )
        return variables

    def SFS_func(self):
        r"""
        Equation for SFS.

        """

        fluid = self.SFS.split_fluid
        out_i = int(self.SFS.split_outlet[3:]) - 1
        i = self.inl[0]
        o = self.outl[out_i]

        res = i.fluid.val[fluid] * i.m.val_SI * self.SFS.val \
            - o.fluid.val[fluid] * o.m.val_SI

        return res

    def SFS_deriv(self, increment_filter, k):
        r"""
        Calculate the partial derivatives for SFS.

        """

        fluid = self.SFS.split_fluid
        out_i = int(self.SFS.split_outlet[3:]) - 1

        i = self.inl[0]
        o = self.outl[out_i]
        if i.m.is_var:
            self.jacobian[k, i.m.J_col] = i.fluid.val[fluid] * self.SFS.val
        if fluid in i.fluid.is_var:
            self.jacobian[k, i.fluid.J_col[fluid]] = i.m.val_SI * self.SFS.val
        if o.m.is_var:
            self.jacobian[k, o.m.J_col] = -o.fluid.val[fluid]
        if fluid in o.fluid.is_var:
            self.jacobian[k, o.fluid.J_col[fluid]] = -o.m.val_SI


class SeparatorWithSpeciesSplitsDeltaT(SeparatorWithSpeciesSplits):

    @staticmethod
    def component():
        return 'separator with species flow splits and dT on outlets'

    def get_parameters(self):
        variables = super().get_parameters()
        variables["deltaT"] = dc_cp(
            deriv=self.energy_balance_deriv, # same as before
            func=self.energy_balance_deltaT_func,
            latex=self.pr_func_doc,
            num_eq=self.num_out
        )
        variables["Q"] = dc_cp(is_result=True)
        #variables["Qout"] = dc_cpa()
        return variables

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        self.variable_fluids = self.variable_fluids = set(self.inl[0].fluid.back_end.keys()) 
        num_fluid_eq = len(self.variable_fluids)
        constraints['fluid_constraints'] = {
            'func': self.fluid_func, 'deriv': self.fluid_deriv,
            'constant_deriv': False, 'latex': self.fluid_func_doc,
            'num_eq': num_fluid_eq}
        del constraints['energy_balance_constraints']
        return constraints

    def energy_balance_deltaT_func(self):
        r"""
        Calculate deltaT residuals.

        """
        i = self.inl[0]
        if i.T.is_set:
            T_in = i.T.val_SI
        else:
            # calculate T_in
            if i.T.val0 > 0:
                T_in = T_mix_ph(i.p.val_SI,i.h.val_SI,i.fluid_data,i.mixing_rule,i.T.val0) 
            else:
                T_in = T_mix_ph(i.p.val_SI,i.h.val_SI,i.fluid_data,i.mixing_rule) 
        
        residual = []
        for o in self.outl:
            residual += [T_in - self.deltaT.val - T_mix_ph(o.p.val_SI,o.h.val_SI,o.fluid_data,o.mixing_rule, T0=T_in)] # use T_in as guess
        return residual

    def calc_parameters(self):
        super().calc_parameters()
        i = self.inl[0]
        self.Q.val = np.sum([o.m.val_SI * (o.h.val_SI - i.h.val_SI) for o in self.outl])

        Tmin = min([o.T.val_SI for o in self.outl])
        Tmax = max([o.T.val_SI for o in self.outl])
        if abs(i.T.val_SI - Tmin) >= abs(i.T.val_SI - Tmax):
            self.deltaT.val = i.T.val_SI - Tmin
        else:
            self.deltaT.val = i.T.val_SI - Tmax

class SeparatorWithSpeciesSplitsDeltaP(SeparatorWithSpeciesSplits):

    @staticmethod
    def component():
        return 'separator with species flow splits and dT and Pr on outlets'

    def get_parameters(self):
        variables = super().get_parameters()
        variables["deltaP"] = dc_cp(
            min_val=0,
            deriv=self.deltaP_deriv,
            func=self.deltaP_func,
            latex=self.pr_func_doc,
            num_eq=self.num_out,
        )
        return variables

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        self.variable_fluids = self.variable_fluids = set(self.inl[0].fluid.back_end.keys()) 
        num_fluid_eq = len(self.variable_fluids)
        constraints['fluid_constraints'] = {
            'func': self.fluid_func, 'deriv': self.fluid_deriv,
            'constant_deriv': False, 'latex': self.fluid_func_doc,
            'num_eq': num_fluid_eq}        
        del constraints['pressure_constraints']
        return constraints   

    def deltaP_func(self):
        r"""
        Equation for pressure drop.

        """
        residual = []
        p_in = self.inl[0].p.val_SI
        for o in self.outl:
            residual += [p_in - self.deltaP.val*1e5 - o.p.val_SI]
        return residual

    def deltaP_deriv(self, increment_filter, k):
        r"""
        Calculate the partial derivatives for pressure drop

        """

        i = self.inl[0]
        for o in self.outl:
            if i.p.is_var:
                self.jacobian[k, i.p.J_col] = 1
            if o.p.is_var:                
                self.jacobian[k, o.p.J_col] = -1
            k += 1            

    def calc_parameters(self):
        super().calc_parameters()

        Pmin = min([i.p.val_SI for i in self.outl])
        Pmax = max([i.p.val_SI for i in self.outl])
        if abs(self.inl[0].p.val_SI - Pmin) >= abs(self.inl[0].p.val_SI - Pmax):
            self.deltaP.val = (self.inl[0].p.val_SI - Pmin)/1e5
        else:
            self.deltaP.val = (self.inl[0].p.val_SI - Pmax)/1e5


class SeparatorWithSpeciesSplitsDeltaTDeltaP(SeparatorWithSpeciesSplitsDeltaT, SeparatorWithSpeciesSplitsDeltaP):

    @staticmethod
    def component():
        return 'separator with species flow splits and dT and Pr on outlets'

    def get_parameters(self):
        variables = super().get_parameters()
        return variables

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        #del constraints['pressure_constraints']
        #del constraints['energy_balance_constraints']
        return constraints


class SeparatorWithSpeciesSplitsDeltaTDeltaPBus(SeparatorWithSpeciesSplitsDeltaTDeltaP):

    @staticmethod
    def component():
        return 'separator with species flow splits and dT and Pr on outlets + Bus connection on Q'

    def get_parameters(self):
        variables = super().get_parameters()
        return variables

    def bus_func(self, bus):
        r"""
        Calculate the value of the bus function.

        """
        return np.sum([o.m.val_SI * (o.h.val_SI - self.inl[0].h.val_SI) for o in self.outl])

    def bus_func_doc(self, bus):
        r"""
        Return LaTeX string of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus.Bus
            TESPy bus object.

        Returns
        -------
        latex : str
            LaTeX string of bus function.
        """
        return (
            r'\dot{m}_\mathrm{in} \cdot \left(h_\mathrm{out} - '
            r'h_\mathrm{in} \right)')

    def bus_deriv(self, bus):
        r"""
        Calculate partial derivatives of the bus function.

        """

        f = self.calc_bus_value
        if self.inl[0].m.is_var:
            if self.inl[0].m.J_col not in bus.jacobian:
                bus.jacobian[self.inl[0].m.J_col] = 0
            bus.jacobian[self.inl[0].m.J_col] -= self.numeric_deriv(f, 'm', self.inl[0], bus=bus)

        if self.inl[0].h.is_var:
            if self.inl[0].h.J_col not in bus.jacobian:
                bus.jacobian[self.inl[0].h.J_col] = 0
            bus.jacobian[self.inl[0].h.J_col] -= self.numeric_deriv(f, 'h', self.inl[0], bus=bus)

        for o in self.outl:
            if o.h.is_var:
                if o.h.J_col not in bus.jacobian:
                    bus.jacobian[o.h.J_col] = 0
                bus.jacobian[o.h.J_col] -= self.numeric_deriv(f, 'h', o, bus=bus)        
            if o.m.is_var:
                if o.m.J_col not in bus.jacobian:
                    bus.jacobian[o.m.J_col] = 0
                bus.jacobian[o.m.J_col] -= self.numeric_deriv(f, 'm', o, bus=bus)        


class SplitterDeltaP(Splitter):

    def __init__(self, label, **kwargs):
        #self.set_attr(**kwargs)
        # need to assign the number of outlets before the variables are set
        for key in kwargs:
            if key == 'num_out':
                self.num_out=kwargs[key]
        super().__init__(label, **kwargs)    

    @staticmethod
    def component():
        return 'Splitter with pressure losses'

    def get_parameters(self):
        variables = super().get_parameters()
        variables["deltaP"] = dc_cp(
            min_val=0,
            deriv=self.deltaP_deriv,
            func=self.deltaP_func,
            latex=self.pr_func_doc,
            num_eq=self.num_out,
        )
        return variables

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        del constraints['pressure_constraints']
        return constraints

    def deltaP_func(self):
        r"""
        Equation for pressure drop.

        """
        #return self.inl[0].p.val_SI * self.pr.val - self.outl[0].p.val_SI
        residual = []
        p_in = self.inl[0].p.val_SI
        for o in self.outl:
            residual += [p_in - self.deltaP.val*1e5 - o.p.val_SI]
        return residual

    def deltaP_deriv(self, increment_filter, k):
        r"""
        Calculate the partial derivatives for combustion pressure ratio.

        """

        i = self.inl[0]
        for o in self.outl:
            if i.p.is_var:
                self.jacobian[k, i.p.J_col] = 1
            if o.p.is_var:                
                self.jacobian[k, o.p.J_col] = -1
            k += 1

    def calc_parameters(self):
        super().calc_parameters()

        Pmin = min([i.p.val_SI for i in self.outl])
        Pmax = max([i.p.val_SI for i in self.outl])
        if abs(self.inl[0].p.val_SI - Pmin) >= abs(self.inl[0].p.val_SI - Pmax):
            self.deltaP.val = (self.inl[0].p.val_SI - Pmin)/1e5
        else:
            self.deltaP.val = (self.inl[0].p.val_SI - Pmax)/1e5

class SplitterWithFlowSplitter(Splitter):

    @staticmethod
    def component():
        return 'splitter with flow split ratios'

    def get_parameters(self):
        variables = super().get_parameters()
        variables["FS"] = dc_cp_FS(
            min_val=0,
            deriv=self.FS_deriv,
            func=self.FS_func,
            latex=self.pr_func_doc,
            num_eq=1,
        )
        return variables

    def FS_func(self):
        r"""
        Equation for flow split.

        """

        out_i = int(self.FS.split_outlet[3:]) - 1
        res = self.inl[0].m.val_SI * self.FS.val - self.outl[out_i].m.val_SI
        return res

    def FS_deriv(self, increment_filter, k):
        r"""
        Calculate the partial derivatives for flow split

        """

        out_i = int(self.FS.split_outlet[3:]) - 1

        i = self.inl[0]
        o = self.outl[out_i]
        if i.m.is_var:
            self.jacobian[k, i.m.J_col]     = self.FS.val
        if o.m.is_var:
            self.jacobian[k, o.m.J_col]     = -1


class SplitterWithFlowSplitterDeltaP(SplitterWithFlowSplitter, SplitterDeltaP):

    @staticmethod
    def component():
        return 'splitter with flow split ratios and pressure drop'

    def get_parameters(self):
        variables = super().get_parameters()
        return variables


#%% Class containers

class dc_cp_SFS(dc_cp):
    """
    Data container for simple properties.
    + SFS_fluid
    + SFS_outlet
    """
    @staticmethod
    def attr():
        attributes = dc_cp.attr()
        attributes.update({'split_fluid' : None, 'split_outlet' : None})
        return attributes
    
    @staticmethod
    def _serializable_keys():
        keys = dc_cp._serializable_keys()
        keys.append("split_fluid")
        keys.append("split_outlet")
        return keys

class dc_cp_FS(dc_cp):
    """
    Data container for component properties.
    + FS_outlet
    """
    @staticmethod
    def attr():
        attributes = dc_cp.attr()
        attributes.update({'split_outlet' : None})
        return attributes

    @staticmethod
    def _serializable_keys():
        keys = dc_cp._serializable_keys()
        keys.append("split_outlet")
        return keys


# class MergeWithPressureLoss(MergeDeltaP):

#     def __init__(self, label, **kwargs):
#         super().__init__(label, **kwargs)
#         msg = (
#             "The API for the component MergeWithPressureLoss will change with "
#             "the next major release, please import MergeDeltaP instead."
#         )
#         warnings.warn(msg, FutureWarning)

# class SeparatorWithSpeciesSplitsAndDeltaT(SeparatorWithSpeciesSplitsDeltaT):

#     def __init__(self, label, **kwargs):
#         super().__init__(label, **kwargs)
#         msg = (
#             "The API for the component SeparatorWithSpeciesSplitsAndDeltaT will change with "
#             "the next major release, please import SeparatorWithSpeciesSplitsDeltaT instead."
#         )
#         warnings.warn(msg, FutureWarning)        

# class SeparatorWithSpeciesSplitsAndDeltaTAndPr(SeparatorWithSpeciesSplitsDeltaTDeltaP):

#     def __init__(self, label, **kwargs):
#         super().__init__(label, **kwargs)
#         msg = (
#             "The API for the component SeparatorWithSpeciesSplitsAndDeltaTAndPr will change with "
#             "the next major release, please import SeparatorWithSpeciesSplitsDeltaTDeltaP instead."
#         )
#         warnings.warn(msg, FutureWarning)   

# class SeparatorWithSpeciesSplitsAndDeltaTAndPrAndBus(SeparatorWithSpeciesSplitsDeltaTDeltaPBus):

#     def __init__(self, label, **kwargs):
#         super().__init__(label, **kwargs)
#         msg = (
#             "The API for the component SeparatorWithSpeciesSplitsAndDeltaTAndPrAndBus will change with "
#             "the next major release, please import SeparatorWithSpeciesSplitsDeltaTDeltaPBus instead."
#         )
#         warnings.warn(msg, FutureWarning)   
        

# class SplitterWithPressureLoss(SplitterDeltaP):

#     def __init__(self, label, **kwargs):
#         super().__init__(label, **kwargs)
#         msg = (
#             "The API for the component SeparatorWithSpeciesSplitsAndDeltaTAndPr will change with "
#             "the next major release, please import SeparatorWithSpeciesSplitsDeltaTDeltaP instead."
#         )
#         warnings.warn(msg, FutureWarning)   
        
# class SplitterWithPressureLoss(SplitterDeltaP):

#     def __init__(self, label, **kwargs):
#         super().__init__(label, **kwargs)
#         msg = (
#             "The API for the component SplitterWithPressureLoss will change with "
#             "the next major release, please import SplitterDeltaP instead."
#         )
#         warnings.warn(msg, FutureWarning)   


      

# class SplitWithFlowSplitter(SplitterWithFlowSplitter):

#     def __init__(self, label, **kwargs):
#         super().__init__(label, **kwargs)
#         msg = (
#             "The API for the component SplitWithFlowSplitter will change with "
#             "the next major release, please import SplitterWithFlowSplitter instead."
#         )
#         warnings.warn(msg, FutureWarning)   


# class SplitWithFlowSplitterDeltaP(SplitterWithFlowSplitterDeltaP):

#     def __init__(self, label, **kwargs):
#         super().__init__(label, **kwargs)
#         msg = (
#             "The API for the component SplitWithFlowSplitterDeltaP will change with "
#             "the next major release, please import SplitterWithFlowSplitterDeltaP instead."
#         )
#         warnings.warn(msg, FutureWarning)   
