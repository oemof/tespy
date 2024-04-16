import logging

from tespy.components import Merge, Splitter
from tespy.components import Sink

from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import ComponentPropertiesArray as dc_cpa

from tespy.components.nodes.base import NodeBase
from tespy.tools.data_containers import SimpleDataContainer as dc_simple


# Fictious Energy Supply models (energy flows modelled as mass flows)
# No real use for tespy I guess

class HeatPumpEnergySupply(Splitter):

    @staticmethod
    def component():
        return 'mass factor vapor compression cycle using COP for converting electricity to heat and cooling (energy flows modelled using tespy mass balances)'

    def outlets(self):
        if self.num_out.is_set:
            return ['out' + str(i + 1) for i in range(self.num_out.val)]
        else:
            self.set_attr(num_out=1)
            return self.outlets()    

    def get_parameters(self):
        variables = super().get_parameters()
        variables["COP"] = dc_cp(
            min_val=0,
            deriv=self.COP_deriv,
            func=self.COP_func,
            latex=self.mass_flow_func_doc,
            num_eq=1
        )
        variables['Heating'] = dc_cp(min_val=0, is_result=True)
        variables['Cooling'] = dc_cp(min_val=0, is_result=True)        
        return variables

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        del constraints['pressure_constraints']
        del constraints['energy_balance_constraints']
        del constraints['mass_flow_constraints']
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
        self.Heating.val = self.outl[0].m.val_SI
        self.Cooling.val = -(self.outl[0].m.val_SI-self.inl[0].m.val_SI)



class HeatPumpUsefullLossEnergySupply(Splitter):

    """
    COP sets self.outl[0].m.val_SI
    UsefullCoolingRatio sets self.outl[1].m.val_SI
    """

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
        variables["UsefullCoolingRatio"] = dc_cp(
            min_val=0,
            deriv=self.usefull_deriv,
            func=self.usefull_func,
            latex=self.mass_flow_func_doc,
            num_eq=1
        )
        variables["NonUseFullCooling"] = dc_cp(min_val=0, is_result=True)
        variables['Heating'] = dc_cp(min_val=0, is_result=True)
        variables['Cooling'] = dc_cp(min_val=0, is_result=True)
        return variables

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        del constraints['pressure_constraints']
        del constraints['energy_balance_constraints']
        del constraints['mass_flow_constraints']
        return constraints   

    def COP_func(self):
        return self.inl[0].m.val_SI * self.COP.val - self.outl[0].m.val_SI

    def COP_deriv(self, increment_filter, k):
        inl = self.inl[0]
        outl = self.outl[0]
        if inl.m.is_var:
            self.jacobian[k, inl.m.J_col] = self.COP.val
        if outl.m.is_var:
            self.jacobian[k, outl.m.J_col] = -1

    def usefull_func(self):
        return self.inl[0].m.val_SI * (self.COP.val-1) * self.UsefullCoolingRatio.val + self.outl[1].m.val_SI # cooling negative

    def usefull_deriv(self, increment_filter, k):
        inl = self.inl[0]
        outl = self.outl[1]
        if inl.m.is_var:
            self.jacobian[k, inl.m.J_col] = (self.COP.val-1) * self.UsefullCoolingRatio.val 
        if outl.m.is_var:
            self.jacobian[k, outl.m.J_col] = 1       

    def calc_parameters(self):
        super().calc_parameters()
        # if not self.COP.is_set:
        #     self.COP.val = self.outl[0].m.val_SI / (self.outl[0].m.val_SI - (-self.outl[1].m.val_SI))
        self.Heating.val = self.outl[0].m.val_SI
        self.Cooling.val = self.outl[1].m.val_SI
        self.NonUseFullCooling.val = self.inl[0].m.val_SI * (self.COP.val - 1) * (1 - self.UsefullCoolingRatio.val)


class RefUnitEnergySupply(Splitter):

    @staticmethod
    def component():
        return 'mass factor vapor compression cycle using COP for converting electricity to heat and cooling (energy flows modelled using tespy mass balances)'

    def outlets(self):
        if self.num_out.is_set:
            return ['out' + str(i + 1) for i in range(self.num_out.val)]
        else:
            self.set_attr(num_out=1)
            return self.outlets()    

    def get_parameters(self):
        variables = super().get_parameters()
        variables["COP"] = dc_cp(
            min_val=0,
            deriv=self.COP_deriv,
            func=self.COP_func,
            latex=self.mass_flow_func_doc,
            num_eq=1
        )
        variables['Heating'] = dc_cp(min_val=0, is_result=True)
        variables['Cooling'] = dc_cp(min_val=0, is_result=True)        
        return variables

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        del constraints['pressure_constraints']
        del constraints['energy_balance_constraints']
        del constraints['mass_flow_constraints']
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
        return self.inl[0].m.val_SI * self.COP.val + self.outl[0].m.val_SI

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
            self.jacobian[k, outl.m.J_col] = 1

    def calc_parameters(self):
        super().calc_parameters()
        self.Cooling.val = self.outl[0].m.val_SI
        self.Heating.val = -self.outl[0].m.val_SI + self.inl[0].m.val_SI


class RefUnitUsefullLossEnergySupply(Splitter):

    """
    COP sets self.outl[1].m.val_SI
    UsefullHeatingRatio sets self.outl[0].m.val_SI
    """


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
        variables["UsefullHeatingRatio"] = dc_cp(
            min_val=0,
            deriv=self.usefull_deriv,
            func=self.usefull_func,
            latex=self.mass_flow_func_doc,
            num_eq=1
        )
        variables["NonUseFullHeating"] = dc_cp(min_val=0, is_result=True)        
        variables['Heating'] = dc_cp(min_val=0, is_result=True)
        variables['Cooling'] = dc_cp(min_val=0, is_result=True)
        return variables

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        del constraints['pressure_constraints']
        del constraints['energy_balance_constraints']
        del constraints['mass_flow_constraints']
        return constraints   

    def COP_func(self):
        return self.inl[0].m.val_SI * self.COP.val + self.outl[1].m.val_SI # cooling negative

    def COP_deriv(self, increment_filter, k):
        inl = self.inl[0]
        outl = self.outl[1]
        if inl.m.is_var:
            self.jacobian[k, inl.m.J_col] = self.COP.val
        if outl.m.is_var:
            self.jacobian[k, outl.m.J_col] = 1

    def usefull_func(self):
        return self.inl[0].m.val_SI * (self.COP.val + 1) * self.UsefullHeatingRatio.val - self.outl[0].m.val_SI

    def usefull_deriv(self, increment_filter, k):
        inl = self.inl[0]
        outl = self.outl[0]
        if inl.m.is_var:
            self.jacobian[k, inl.m.J_col] =  (self.COP.val + 1) * self.UsefullHeatingRatio.val
        if outl.m.is_var:
            self.jacobian[k, outl.m.J_col] = -1       

    def calc_parameters(self):
        super().calc_parameters()
        # if not self.COP.is_set:
        #     self.COP.val = self.outl[0].m.val_SI / (self.outl[0].m.val_SI - (-self.outl[1].m.val_SI))
        self.Heating.val = self.outl[0].m.val_SI
        self.Cooling.val = self.outl[1].m.val_SI
        self.NonUseFullHeating.val = self.inl[0].m.val_SI * (self.COP.val + 1) * (1 - self.UsefullHeatingRatio.val)



class MassLossEnergySupply(Splitter):

    @staticmethod
    def component():
        return 'mass loss model for splitting energy flows (modelled using tespy mass balances)'

    def get_parameters(self):
        variables = super().get_parameters()
        variables["LossRatio"] = dc_cp(
            min_val=0,
            deriv=self.Loss_deriv,
            func=self.Loss_func,
            latex=self.mass_flow_func_doc,
            num_eq=1
        )
        variables['Energy'] = dc_cp(min_val=0, is_result=True)           
        variables['EnergyLoss'] = dc_cp(min_val=0, is_result=True)           
        return variables
    
    def outlets(self):
        if self.num_out.is_set:
            return ['out' + str(i + 1) for i in range(self.num_out.val)]
        else:
            self.set_attr(num_out=1)
            return self.outlets()    

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        del constraints['mass_flow_constraints']
        del constraints['pressure_constraints']
        del constraints['energy_balance_constraints']
        return constraints   

    def Loss_func(self):
        return self.inl[0].m.val_SI * (1-self.LossRatio.val) - sum([o.m.val_SI for o in self.outl])

    def Loss_deriv(self, increment_filter, k):
        inl = self.inl[0]
        if inl.m.is_var:
            self.jacobian[k, inl.m.J_col] = (1-self.LossRatio.val)
        for o in self.outl:
            if o.m.is_var:
                self.jacobian[k, o.m.J_col] = -1      

    def calc_parameters(self):
        super().calc_parameters()
        mout = sum([o.m.val_SI for o in self.outl])
        if not self.LossRatio.is_set:
            self.LossRatio.val = (self.inl[0].m.val_SI - mout)/self.inl[0].m.val_SI
        self.EnergyLoss.val = self.inl[0].m.val_SI - mout
        self.Energy.val = mout

class BoilerEffEnergySupply(Splitter):

    @staticmethod
    def component():
        return 'mass efficiency model for splitting energy flows (modelled using tespy mass balances)'

    def get_parameters(self):
        variables = super().get_parameters()
        variables["BoilerEff"] = dc_cp(
            min_val=0,
            deriv=self.Eff_deriv,
            func=self.Eff_func,
            latex=self.mass_flow_func_doc,
            num_eq=1
        )
        variables['Energy'] = dc_cp(min_val=0, is_result=True)           
        variables['EnergyLoss'] = dc_cp(min_val=0, is_result=True)           
        return variables
    
    def outlets(self):
        if self.num_out.is_set:
            return ['out' + str(i + 1) for i in range(self.num_out.val)]
        else:
            self.set_attr(num_out=1)
            return self.outlets()    

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        del constraints['mass_flow_constraints']
        del constraints['pressure_constraints']
        del constraints['energy_balance_constraints']
        return constraints   

    def Eff_func(self):
        return self.inl[0].m.val_SI * self.BoilerEff.val - sum([o.m.val_SI for o in self.outl])

    def Eff_deriv(self, increment_filter, k):
        inl = self.inl[0]
        if inl.m.is_var:
            self.jacobian[k, inl.m.J_col] = self.BoilerEff.val
        for o in self.outl:
            if o.m.is_var:
                self.jacobian[k, o.m.J_col] = -1

    def calc_parameters(self):
        super().calc_parameters()
        mout = sum([o.m.val_SI for o in self.outl])
        if not self.BoilerEff.is_set:
            self.BoilerEff.val = mout/self.inl[0].m.val_SI
        self.EnergyLoss.val = self.inl[0].m.val_SI - mout
        self.Energy.val = mout



class BoilerEffUsefullLossEnergySupply(Splitter):

    """
    BoilerEff defines self.outl[0].m.val_SI
    UsefullLossRatio defines self.outl[1].m.val_SI
    """

    @staticmethod
    def component():
        return 'mass efficiency model for splitting energy flows (modelled using tespy mass balances)'

    def get_parameters(self):
        variables = super().get_parameters()
        variables["BoilerEff"] = dc_cp(
            min_val=0,
            deriv=self.Eff_deriv,
            func=self.Eff_func,
            latex=self.mass_flow_func_doc,
            num_eq=1
        )
        variables["UsefullLossRatio"] = dc_cp(
            min_val=0,
            deriv=self.usefull_deriv,
            func=self.usefull_func,
            latex=self.mass_flow_func_doc,
            num_eq=1
        )        
        variables['Energy'] = dc_cp(min_val=0, is_result=True)           
        variables['EnergyLoss'] = dc_cp(min_val=0, is_result=True)           
        variables['EnergyLossUsefull'] = dc_cp(min_val=0, is_result=True)           
        return variables
    
    def outlets(self):
        if self.num_out.is_set:
            return ['out' + str(i + 1) for i in range(self.num_out.val)]
        else:
            self.set_attr(num_out=2)
            return self.outlets() 

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        del constraints['mass_flow_constraints']
        del constraints['pressure_constraints']
        del constraints['energy_balance_constraints']
        return constraints   

    def Eff_func(self):
        return self.inl[0].m.val_SI * self.BoilerEff.val - self.outl[0].m.val_SI

    def Eff_deriv(self, increment_filter, k):
        inl = self.inl[0]
        outl = self.outl[0]
        if inl.m.is_var:
            self.jacobian[k, inl.m.J_col] = self.BoilerEff.val
        if outl.m.is_var:
            self.jacobian[k, outl.m.J_col] = -1

    def usefull_func(self):
        return self.inl[0].m.val_SI * (1-self.BoilerEff.val) * self.UsefullLossRatio.val - self.outl[1].m.val_SI

    def usefull_deriv(self, increment_filter, k):
        inl = self.inl[0]
        outl = self.outl[1]
        if inl.m.is_var:
            self.jacobian[k, inl.m.J_col] = (1-self.BoilerEff.val) * self.UsefullLossRatio.val
        if outl.m.is_var:
            self.jacobian[k, outl.m.J_col] = -1

    def calc_parameters(self):
        super().calc_parameters()
        # mout = sum([o.m.val_SI for o in self.outl])
        # if not self.BoilerEff.is_set:
        #     self.BoilerEff.val = mout/self.inl[0].m.val_SI
        self.EnergyLossUsefull.val = self.outl[1].m.val_SI
        self.EnergyLoss.val = self.inl[0].m.val_SI - sum([o.m.val_SI for o in self.outl])
        self.Energy.val = self.outl[0].m.val_SI



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
    


# class SourceEnergy(NodeBase):

#     def __init__(self, label, **kwargs):
#         #self.set_attr(**kwargs)
#         # need to assign the number of outlets before the variables are set
#         for key in kwargs:
#             if key == 'num_out':
#                 self.num_out=kwargs[key]
#         super().__init__(label, **kwargs)    

#     @staticmethod
#     def component():
#         return 'Source'

#     def get_parameters(self):
#         variables = super().get_parameters()
#         variables['num_out'] = dc_simple()
#         variables["Energy"] = dc_cp(
#             min_val=0,
#             deriv=self.mass_flow_deriv,
#             func=self.mass_flow_func,
#             latex=self.mass_flow_func_doc,
#             num_eq=1
#         )        
#         return variables

#     def get_mandatory_constraints(self):
#         return {}

#     def outlets(self):
#         if self.num_out.is_set:
#             return ['out' + str(i + 1) for i in range(self.num_out.val)]
#         else:
#             self.set_attr(num_out=2)
#             return self.outlets()

#     @staticmethod
#     def is_branch_source():
#         return True

#     def start_branch(self):
#         branches = {}
#         for outconn in self.outl:
#             branch = {
#                 "connections": [outconn],
#                 "components": [self, outconn.target],
#                 "subbranches": {}
#             }
#             outconn.target.propagate_to_target(branch)
#             branches[outconn.label] = branch
#         return branches

#     def start_fluid_wrapper_branch(self):
#         branches = {}
#         for outconn in self.outl:
#             branch = {
#                 "connections": [outconn],
#                 "components": [self]
#             }
#             outconn.target.propagate_wrapper_to_target(branch)
#             branches[outconn.label] = branch
#         return branches

#     def mass_flow_func(self):
#         r"""
#         Calculate the residual value for mass flow balance equation.

#         Returns
#         -------
#         res : float
#             Residual value of equation.

#             .. math::

#                 0 = \sum \dot{m}_{in,i} - \sum \dot{m}_{out,j} \;
#                 \forall i \in inlets, \forall j \in outlets
#         """
#         res = self.Energy.val
#         for o in self.outl:
#             res -= o.m.val_SI
#         return res

#     def mass_flow_deriv(self, increment_filter, k):
#         r"""
#         Calculate partial derivatives for mass flow equation.

#         Returns
#         -------
#         deriv : list
#             Matrix with partial derivatives for the fluid equations.
#         """
#         for o in self.outl:
#             if o.m.is_var:
#                 self.jacobian[k, o.m.J_col] = -1

#     def calc_parameters(self):
#         super().calc_parameters()
#         if not self.Energy.is_set:
#             self.Energy.val = sum([o.m.val_SI for o in self.outl])


class SourceEnergy(NodeBase):

    def __init__(self, label, **kwargs):
        #self.set_attr(**kwargs)
        # need to assign the number of outlets before the variables are set
        for key in kwargs:
            if key == 'num_out':
                self.num_out=kwargs[key]
        super().__init__(label, **kwargs)    

    @staticmethod
    def component():
        return 'Source'

    def get_parameters(self):
        variables = super().get_parameters()
        variables['num_out'] = dc_simple()
        variables["Energy"] = dc_cp(
            min_val=0,
            deriv=self.mass_flow_deriv,
            func=self.mass_flow_func,
            latex=self.mass_flow_func_doc,
            num_eq=1
        )        
        variables["EnergyArray"] = dc_cpa(
            min_val=0,
            deriv=self.mass_flow_array_deriv,
            func=self.mass_flow_array_func,
            latex=self.mass_flow_func_doc,
            #num_eq=self.num_out
        )        
                
        return variables

    def get_mandatory_constraints(self):
        return {}

    def outlets(self):
        if self.num_out.is_set:
            return ['out' + str(i + 1) for i in range(self.num_out.val)]
        else:
            self.set_attr(num_out=2)
            return self.outlets()

    @staticmethod
    def is_branch_source():
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

    def mass_flow_func(self):
        res = self.Energy.val
        for o in self.outl:
            res -= o.m.val_SI
        return res

    def mass_flow_deriv(self, increment_filter, k):
        for o in self.outl:
            if o.m.is_var:
                self.jacobian[k, o.m.J_col] = -1

    def mass_flow_array_func(self):
        residual  = []
        for i,is_set in enumerate(self.EnergyArray.is_set):
            if is_set:
                o = self.outl[i]
                residual += [self.EnergyArray.val[i] - o.m.val_SI]
        return residual

    def mass_flow_array_deriv(self, increment_filter, k):
        m=0
        for i,is_set in enumerate(self.EnergyArray.is_set):
            if is_set:
                o = self.outl[i]
                self.jacobian[k + m, o.m.J_col] = -1
                m=m+1

    def calc_parameters(self):
        super().calc_parameters()
        if not self.Energy.is_set:
            self.Energy.val = sum([o.m.val_SI for o in self.outl])
        for i,o in enumerate(self.outl):
            self.EnergyArray.val[i] = o.m.val_SI


class SinkEnergy(Sink):

    @staticmethod
    def component():
        return 'sink with energy '

    def get_parameters(self):
        variables = super().get_parameters()
        variables["Energy"] = dc_cp(
            min_val=0,
            deriv=self.mass_flow_deriv,
            func=self.mass_flow_func,
            latex=self.mass_flow_func_doc,
            num_eq=1
        )        
        return variables        

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        return constraints


    def mass_flow_func(self):
        r"""
        Calculate the residual value for mass flow balance equation.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                0 = \sum \dot{m}_{in,i} - \sum \dot{m}_{out,j} \;
                \forall i \in inlets, \forall j \in outlets
        """
        res = self.Energy.val
        for i in self.inl:
            res -= i.m.val_SI
        return res
    
    def mass_flow_func_doc(self, label):
        pass

    def mass_flow_deriv(self, increment_filter, k):
        r"""
        Calculate partial derivatives for mass flow equation.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        for i in self.inl:
            if i.m.is_var:
                self.jacobian[k, i.m.J_col] = -1

    def calc_parameters(self):
        super().calc_parameters()
        if not self.Energy.is_set:
            self.Energy.val = self.inl[0].m.val_SI


class FlowEnergy(Splitter):

    @staticmethod
    def component():
        return 'flow with energy '

    def get_parameters(self):
        variables = super().get_parameters()
        variables["Energy"] = dc_cp(
            min_val=0,
            deriv=self.Energy_mass_flow_deriv,
            func=self.Energy_mass_flow_func,
            latex=self.Energy_mass_flow_func_doc,
            num_eq=1
        )        
        return variables        

    def outlets(self):
        self.set_attr(num_out=1)
        return ['out1']

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        #del constraints['mass_flow_constraints']
        del constraints['pressure_constraints']
        del constraints['energy_balance_constraints']
        return constraints       

    def Energy_mass_flow_func(self):
        res = self.Energy.val
        for i in self.inl:
            res -= i.m.val_SI
        return res
    
    def Energy_mass_flow_deriv(self, increment_filter, k):
        for i in self.inl:
            if i.m.is_var:
                self.jacobian[k, i.m.J_col] = -1    
    
    def Energy_mass_flow_func_doc(self, label):
        pass

    def calc_parameters(self):
        super().calc_parameters()
        if not self.Energy.is_set:
            self.Energy.val = self.inl[0].m.val_SI        