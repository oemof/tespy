# -*- coding: utf-8

"""Module of class Pipe.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/piping/pipe.py

SPDX-License-Identifier: MIT
"""
import math
#import warnings

#import numpy as np

from tespy.components.heat_exchangers.simple import SimpleHeatExchanger
from tespy.components.component import component_registry
#from tespy.tools import logger
#from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
from tespy.tools.document_models import generate_latex_eq
#from tespy.tools.fluid_properties import s_mix_ph
#from tespy.tools.fluid_properties.helpers import darcy_friction_factor as dff
#from tespy.tools.helpers import convert_to_SI
from tespy.tools.fluid_properties.wrappers import CoolPropWrapper

@component_registry
class Pipeline(SimpleHeatExchanger):
    r"""
    A basic heat exchanger representing a heat source or heat sink.

    The component SimpleHeatExchanger is the parent class for the components:

    - :py:class:`tespy.components.heat_exchangers.solar_collector.SolarCollector`
    - :py:class:`tespy.components.heat_exchangers.parabolic_trough.ParabolicTrough`
    - :py:class:`tespy.components.piping.pipe.Pipe`

    **Mandatory Equations**

    - :py:meth:`tespy.components.component.Component.fluid_func`
    - :py:meth:`tespy.components.component.Component.mass_flow_func`

    **Optional Equations**

    - :py:meth:`tespy.components.component.Component.pr_func`
    - :py:meth:`tespy.components.component.Component.zeta_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.energy_balance_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.darcy_group_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.hw_group_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.kA_group_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.kA_char_group_func`

    Inlets/Outlets

    - in1
    - out1

    Image

    .. image:: /api/_images/Pipe.svg
       :alt: flowsheet of the simple heat exchanger
       :align: center
       :class: only-light

    .. image:: /api/_images/Pipe_darkmode.svg
       :alt: flowsheet of the simple heat exchanger
       :align: center
       :class: only-dark

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    design_path : str
        Path to the components design case.

    local_offdesign : boolean
        Treat this component in offdesign mode in a design calculation.

    local_design : boolean
        Treat this component in design mode in an offdesign calculation.

    char_warnings : boolean
        Ignore warnings on default characteristics usage for this component.

    printout : boolean
        Include this component in the network's results printout.

    Q : float, dict, :code:`"var"`
        Heat transfer, :math:`Q/\text{W}`.

    pr : float, dict, :code:`"var"`
        Outlet to inlet pressure ratio, :math:`pr/1`.

    zeta : float, dict, :code:`"var"`
        Geometry independent friction coefficient,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    D : float, dict, :code:`"var"`
        Diameter of the pipes, :math:`D/\text{m}`.

    L : float, dict, :code:`"var"`
        Length of the pipes, :math:`L/\text{m}`.

    ks : float, dict, :code:`"var"`
        Pipe's roughness, :math:`ks/\text{m}`.

    darcy_group : str, dict
        Parametergroup for pressure drop calculation based on pipes dimensions
        using darcy weissbach equation.

    ks_HW : float, dict, :code:`"var"`
        Pipe's roughness, :math:`ks/\text{1}`.

    hw_group : str, dict
        Parametergroup for pressure drop calculation based on pipes dimensions
        using hazen williams equation.

    kA : float, dict, :code:`"var"`
        Area independent heat transfer coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char : tespy.tools.characteristics.CharLine, dict
        Characteristic line for heat transfer coefficient.

    Tamb : float, dict
        Ambient temperature, provide parameter in network's temperature unit.

    kA_group : str, dict
        Parametergroup for heat transfer calculation from ambient temperature
        and area independent heat transfer coefficient kA.

    insulation_m: float
        thickness of insulation

    insulation_tc: float 
        thermal conductivity insulation
        
    material: str, float
        material of pipe: 'Steel', 'Carbon Steel', 'Cast Iron', 'Stainless Steel', 'PVC', 'CommercialCopper'
        or heat conductivity of material: float

    pipe_thickness: float
        thickness of pipe
    
    environment_media: str
        environment media around the pipe: air, 'gravel' , 'stones' ,'dry soil', 'moist soil'
        
    wind_velocity: float
        Mean velocity of the wind #in m/s
        
    pipe_depth: float
        pipe depth in thew ground #in m

    

    Example
    -------
    >>> from tespy.components import Sink, Source,PipeNet
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> nw = Network()
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg', m_unit = 'kg / s',iterinfo=False)
    >>> so = Source('source 1')
    >>> si = Sink('sink 1')
    >>> pi = PipeNet('pipeline')
    >>> inc = Connection(so, 'out1', pi, 'in1')
    >>> outg = Connection(pi, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)
    >>> inc.set_attr(fluid={'water': 1}, m=100, T=300, p=40)
    >>> pi.set_attr(pr=0.95, Tamb = 20, L=1000, D='var',  ks=4.57e-5,
    >>>                     insulation_m=0.1 ,insulation_tc= 0.035, pipe_thickness=0.003,material='Steel', 
    >>>                     environment_media= 'air', wind_velocity = 0.1,
    >>>                     #environment_media= 'dry soil', pipe_depth = 5,
    >>>             ) 
    >>> nw.solve('design')
    >>> nw.save('tmp')
    """

    @staticmethod
    def component():
        return 'Pipeline'

    def get_parameters(self):
        parameters=super().get_parameters()
        parameters['Q_ohc_group']=dc_gcp(
                elements=['insulation_m','insulation_tc', 'Tamb', 'material',"pipe_thickness"], num_eq=1,
                latex=self.ohc_group_func_doc,
                func=self.ohc_group_func, 
                deriv=self.ohc_group_deriv)
        parameters['insulation_m']=dc_cp(min_val=1e-3, max_val=1e1)
        parameters['insulation_tc']=dc_cp(min_val=1e-3, max_val=1e2) #thermal conductivity insulation
        parameters['material']=dc_simple(val='Steel')
        parameters['pipe_thickness']=dc_cp(min_val=0, max_val=1e-1)
        parameters['environment_media']=dc_simple(val='air') 
        parameters['wind_velocity']=dc_cp(min_val=1e-3, max_val=10) #in m/s
        parameters['pipe_depth']= dc_cp(min_val=1e-2, max_val=1e2) #in m
        return(parameters)
        
    def ohc_group_func(self): 
        T_in = self.inl[0].calc_T() 
        T_out = self.outl[0].calc_T() 
        
        Diameters= [self.D.val, 
            self.D.val+2*self.pipe_thickness.val, 
            self.D.val+2*self.pipe_thickness.val +2* self.insulation_m.val]
        air= CoolPropWrapper('air')

        A= self.L.val * math.pi *(Diameters[2] ) #outer surface area per definition
        dTm= (T_out - T_in) / math.log((self.Tamb.val - T_in) / (self.Tamb.val - T_out))
        
        #heat transfer resistance
        R_sum=[]
        
        '''
        inner heat transfer resistance neglected yet
        R_int = 1/alpha_i *Diameters[2]/ Diameters[0]        
        R_sum.append(R_int)
        '''

        # pipe wall heat transfer resistance
        if Diameters[1] > Diameters[0]: 
            if isinstance(self.material.val, str): 
                R_sum.append(Diameters[1]/self.pipe_tc(self.material.val, (T_in -T_out)/2) *math.log( Diameters[1]/ Diameters[0])/2) 
            if isinstance(self.material.val, float):
                R_sum.append(Diameters[1]/self.material.val *math.log( Diameters[1]/ Diameters[0])/2)
            #potential bug in DWSIM: factor 2 is missing
        
        # insulation heat transfer resistance
        if self.insulation_m.val != 0: 
            R_sum.append(Diameters[2]/self.insulation_tc.val *math.log(Diameters[2]/Diameters[1])/2) 
            #potential bug in DWSIM: factor 2 is missing
        # external heat transfer resistance (to environment)
        lambda_ground ={'gravel':1.1, 'stones':1.95, 'dry soil':0.5, 'moist soil':2.2}
        if self.environment_media.val == 'air':
            '''
            Gnielinski, V.: Berechnung mittlerer Wärme- und Stoffübergangskoeffizienten an laminar 
            und turbulent überströmten Einzelkörpern mithilfe einer einheitlichen Gleichung. 
            Forsch. Ing.-Wes. 41(5), 145–153 (1975)
            '''
            Re= self.wind_velocity.val *math.pi/2*(Diameters[1] +self.insulation_m.val *2)/ air.viscosity_pT(101300,self.Tamb.val+273)*air.d_pT(101300,self.Tamb.val+273)
            Pr=air.AS.Prandtl()
            Nu_lam =0.664 *Re**0.5 *Pr**(1/3) 
            Nu_turb = 0.037* Re**(0.8) * Pr /(1+ 2.443 * Re**(-0.1)* (Pr**(2/3)-1))
            Nu_ext =0.3 +(Nu_lam**2 +Nu_turb**2)**0.5
            alpha_ext = Nu_ext/(math.pi/2*(Diameters[1] +self.insulation_m.val *2) )  * air.AS.conductivity() #W/m²/K
            R_sum.append(1/alpha_ext)
            
        elif self.environment_media.val in lambda_ground.keys():
            '''
            First order approximation of multipole method for a single pipe in the ground.
            No surface resistance
            Source:
            Wallentén P. Steady-state heat loss from insulated pipes. Licentiate Thesis. Division of Building Physics;1991
            '''
            Beta= lambda_ground[self.environment_media.val]/self.insulation_tc.val * math.log(Diameters[2]/ Diameters[0])
            _h = math.log(2*self.pipe_depth.val/Diameters[0])+ Beta + 1/ (1- (2*self.pipe_depth.val/Diameters[0])**2 *(1+Beta)/(1-Beta)) 
            R_soil = _h / (2* math.pi*lambda_ground[self.environment_media.val])
            return self.inl[0].m.val_SI * (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI
            )- dTm  / R_soil
        
        else:
            raise Exception('No valid environmental media.')
        if len(R_sum)==0:
            raise Exception("No heat transfer resistance. Check input values.") 
                
        return self.inl[0].m.val_SI * (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI
        )- dTm * A / sum(R_sum)
    
    def pipe_tc(self, material: str, T: float) -> float:
        
        '''
        Function from DWSIM:
        https://github.com/DanWBR/dwsim/blob/windows/DWSIM.UnitOperations/UnitOperations/Pipe.vb#L1370-L1405

        '''
        
        lamda = 0.0
        if material in ["AoComum", "Steel"]:
            lamda = (-0.000000004 * T ** 3) - (0.00002 * T ** 2) + (0.021 * T) + 33.743
        elif material in ["AoCarbono", "CarbonSteel", "Carbon Steel"]:
            lamda = (0.000000007 * T ** 3) - (0.00002 * T ** 2) - (0.0291 * T) + 70.765
        elif material in ["FerroBottomido", "CastIron", "Cast Iron"]:
            lamda = (-0.00000008 * T ** 3) + (0.0002 * T ** 2) - (0.211 * T) + 127.99
        elif material in ["AoInoxidvel", "StainlessSteel", "Stainless Steel"]:
            lamda = 14.6 + 0.0127 * (T - 273.15)
        elif material in ["PVC", "PVC+PFRV"]:
            lamda = 0.16
        elif material in ["CommercialCopper", "CommercialCopper"]:
            lamda = 420.75 - 0.068493 * T

        return lamda  # W/m.K

    def ohc_group_deriv(self, increment_filter, k): 
        i = self.inl[0]
        o = self.outl[0]
        func= self.ohc_group_func
        if i.m.is_var:
            self.jacobian[k, i.m.J_col] = o.h.val_SI - i.h.val_SI
        if i.h.is_var:
            self.jacobian[k, i.h.J_col] = self.numeric_deriv(func, 'h', i)
        if o.h.is_var:
            self.jacobian[k, o.h.J_col] = self.numeric_deriv(func, 'h', o)
        if i.p.is_var:
            self.jacobian[k, i.p.J_col] = self.numeric_deriv(func, 'p', i)
        if o.p.is_var:
            self.jacobian[k, o.p.J_col] = self.numeric_deriv(func, 'p', o)
        
    def ohc_group_func_doc(self, label):
        r"""
        Equation for ohc calculation.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        latex = (
            r'0 = \dot{m}_\mathrm{in} \cdot \left(h_\mathrm{out} - '
            r'h_\mathrm{in} \right) -\dot{Q}'
        )
        return generate_latex_eq(self, latex, label)