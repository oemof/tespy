import logging

from tespy.components import SimpleHeatExchanger, Merge, Separator, Splitter
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.fluid_properties import T_mix_ph

from tespy.components.component import Component

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


class MergeWithPressureLoss(Merge):

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


class SplitterWithPressureLoss(Splitter):

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


class SeparatorWithSpeciesSplits(Separator):

    def __init__(self, label, **kwargs):
        #self.set_attr(**kwargs)
        # need to assign the number of outlets before the variables are set
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
        # residual = []
        # for fluid, x in self.inl[0].fluid.val.items():
        #     res = x * self.inl[0].m.val_SI
        #     for o in self.outl:
        #         res -= o.fluid.val[fluid] * o.m.val_SI
        #     residual += [res]
        # return residual

        fluid = self.SFS.split_fluid
        out_i = int(self.SFS.split_outlet[3:]) - 1
        inl = self.inl[0]
        outl = self.outl[out_i]

        res = inl.fluid.val[fluid] * inl.m.val_SI * self.SFS.val \
            - outl.fluid.val[fluid] * outl.m.val_SI

        #print(res)
        return res

    def SFS_deriv(self, increment_filter, k):
        r"""
        Calculate the partial derivatives for SFS.

        """

        # j=0
        # self.jacobian[k, j, 0] = self.inl[j].fluid.val[self.split_fluid] * self.TS.val
        # self.jacobian[k, j, i + 3] = self.inl[j].m.val_SI * self.TS.val

        # i = 0
        # for fluid, x in self.outl[0].fluid.val.items():
        #     j = 0
        #     for inl in self.inl:
        #         self.jacobian[k, j, 0] = inl.fluid.val[fluid]
        #         self.jacobian[k, j, i + 3] = inl.m.val_SI
        #         j += 1
        #     self.jacobian[k, j, 0] = -x
        #     self.jacobian[k, j, i + 3] = -self.outl[0].m.val_SI
        #     i += 1
        #     k += 1

        fluid_index = list(self.inl[0].fluid.val.keys()).index(self.SFS.split_fluid)
        fluid = self.SFS.split_fluid
        out_i = int(self.SFS.split_outlet[3:]) - 1

        i = fluid_index
        j = 0
        inl = self.inl[0]
        outl = self.outl[out_i]
        if inl.m.is_var:
            self.jacobian[k, inl.m.J_col] = inl.fluid.val[fluid] * self.SFS.val

        if fluid in inl.fluid.is_var:
            self.jacobian[k, inl.fluid.J_col[fluid]] = inl.m.val_SI * self.SFS.val

        if outl.m.is_var:
            self.jacobian[k, outl.m.J_col] = -outl.fluid.val[fluid]
        if fluid in outl.fluid.is_var:
            self.jacobian[k, outl.fluid.J_col[fluid]] = -outl.m.val_SI

        #print(self.jacobian)
        #print(self.jacobian[k,:,:])




class SeparatorWithSpeciesSplitsAndDeltaT(SeparatorWithSpeciesSplits):



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
        del constraints['energy_balance_constraints']
        return constraints

    def energy_balance_deltaT_func(self):
        r"""
        Calculate delta T derivatives.

        """
        residual = []
        T_in = T_mix_ph(self.inl[0].get_flow(), T0=self.inl[0].T.val_SI)
        i=0
        for o in self.outl:
            residual += [T_in - self.deltaT.val - T_mix_ph(o.get_flow(), T0=o.T.val_SI)]
            i+=1
        return residual

    def calc_parameters(self):
        super().calc_parameters()
        self.Q.val = np.sum([o.m.val_SI * (o.h.val_SI - self.inl[0].h.val_SI) for o in self.outl])

        Tmin = min([i.T.val_SI for i in self.outl])
        Tmax = max([i.T.val_SI for i in self.outl])
        if abs(self.inl[0].T.val_SI - Tmin) >= abs(self.inl[0].T.val_SI - Tmax):
            self.deltaT.val = self.inl[0].T.val_SI - Tmin
        else:
            self.deltaT.val = self.inl[0].T.val_SI - Tmax
        # self.inl[0].T.val_SI - min([i.T.val_SI for i in self.outl])

class SeparatorWithSpeciesSplitsAndPr(SeparatorWithSpeciesSplits):

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
        j = 0
        for c in self.outl:
            self.jacobian[k, 0, 1] = 1
            self.jacobian[k, j + 1, 1] = -1
            j += 1
            k += 1

    def calc_parameters(self):
        super().calc_parameters()

        Pmin = min([i.p.val_SI for i in self.outl])
        Pmax = max([i.p.val_SI for i in self.outl])
        if abs(self.inl[0].p.val_SI - Pmin) >= abs(self.inl[0].p.val_SI - Pmax):
            self.deltaP.val = (self.inl[0].p.val_SI - Pmin)/1e5
        else:
            self.deltaP.val = (self.inl[0].p.val_SI - Pmax)/1e5


class SeparatorWithSpeciesSplitsAndDeltaTAndPr(SeparatorWithSpeciesSplitsAndDeltaT, SeparatorWithSpeciesSplitsAndPr):

    @staticmethod
    def component():
        return 'separator with species flow splits and dT and Pr on outlets'

    def get_parameters(self):
        variables = super().get_parameters()
        return variables

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        del constraints['pressure_constraints']
        del constraints['energy_balance_constraints']
        return constraints


class SeparatorWithSpeciesSplitsAndDeltaTAndPrAndBus(SeparatorWithSpeciesSplitsAndDeltaTAndPr):

    @staticmethod
    def component():
        return 'separator with species flow splits and dT and Pr on outlets + Bus connection on Q'

    def get_parameters(self):
        variables = super().get_parameters()
        return variables

    def bus_func(self, bus):
        r"""
        Calculate the value of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus.Bus
            TESPy bus object.

        Returns
        -------
        val : float
            Value of energy transfer :math:`\dot{E}`. This value is passed to
            :py:meth:`tespy.components.component.Component.calc_bus_value`
            for value manipulation according to the specified characteristic
            line of the bus.

            .. math::

                \dot{E} = \dot{m}_{in} \cdot \left( h_{out} - h_{in} \right)
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

        Parameters
        ----------
        bus : tespy.connections.bus.Bus
            TESPy bus object.

        Returns
        -------
        deriv : ndarray
            Matrix of partial derivatives.
        """
#        for o in self.outl:
#            self.Qout.val += [o.m.val_SI * (o.h.val_SI - self.inl[0].h.val_SI)]
#        return np.sum(self.Qout.val)

        deriv = np.zeros((1, len(self.outl)+1, self.num_nw_vars))
        f = self.calc_bus_value
        deriv[0, 0, 2] = self.numeric_deriv(f, 'h', 0, bus=bus)
        i = 0
        for o in self.outl:
            i = i+1
            deriv[0, i, 0] = self.numeric_deriv(f, 'm', i, bus=bus)
            deriv[0, i, 2] = self.numeric_deriv(f, 'h', i, bus=bus)
        return deriv




class SplitWithFlowSplitter(Splitter):

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
        Equation for pressure drop.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = p_\mathrm{in,1} \cdot pr - p_\mathrm{out,1}
        """

        out_i = int(self.FS.split_outlet[3:]) - 1
        res = self.inl[0].m.val_SI * self.FS.val - self.outl[out_i].m.val_SI

        #print(res)
        return res

    def FS_deriv(self, increment_filter, k):
        r"""
        Calculate the partial derivatives for combustion pressure ratio.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of equation in Jacobian matrix.
        """

        out_i = int(self.FS.split_outlet[3:]) - 1

        inl = self.inl[0]
        outl = self.outl[out_i]
        j = 0
        if inl.m.is_var:
            self.jacobian[k, inl.m.J_col]     = self.FS.val
        if outl.m.is_var:
            self.jacobian[k, outl.m.J_col]     = -1

        #print(self.jacobian)
        #print(self.jacobian[k,:,:])


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
        attributes.update({'split_fluid' : str, 'split_outlet' : str})
        return attributes

class dc_cp_FS(dc_cp):
    """
    Data container for component properties.
    + FS_outlet
    """
    @staticmethod
    def attr():
        attributes = dc_cp.attr()
        attributes.update({'split_outlet' : str})
        return attributes



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
        return {
            'energy_balance_constraints': {
                'func': self.energy_balance_func,
                'deriv': self.energy_balance_deriv,
                'constant_deriv': True, 'latex': self.energy_balance_func_doc,
                'num_eq': self.num_o},
            'pressure_constraints': {
                'func': self.pressure_equality_func,
                'deriv': self.pressure_equality_deriv,
                'constant_deriv': True,
                'latex': self.pressure_equality_func_doc,
                'num_eq': self.num_i + self.num_o - 1}
        }


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
            'pressure_constraints': {
                'func': self.pressure_equality_func,
                'deriv': self.pressure_equality_deriv,
                'constant_deriv': True,
                'latex': self.pressure_equality_func_doc,
                'num_eq': self.num_i + self.num_o - 1}
        }


    def Loss_func(self):
        return self.inl[0].m.val_SI * (1-self.Loss.val) - self.outl[0].m.val_SI

    def Loss_deriv(self, increment_filter, k):
        self.jacobian[k  ,            0, 0] = (1-self.Loss.val)
        self.jacobian[k  ,   self.num_i, 0] = -1

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
        return {
            'mass_flow_constraints': {
                'func': self.mass_flow_func, 'deriv': self.mass_flow_deriv,
                'constant_deriv': True, 'latex': self.mass_flow_func_doc,
                'num_eq': 1},
            'fluid_constraints': {
                'func': self.fluid_func, 'deriv': self.fluid_deriv,
                'constant_deriv': False, 'latex': self.fluid_func_doc,
                'num_eq': self.num_nw_fluids},
        }

class SplitterEnergySupply(Splitter):

    @staticmethod
    def component():
        return 'Splitter without pressure/energy constraints'

    def get_parameters(self):
        variables = super().get_parameters()
        return variables

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
        }



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
        }
