# -*- coding: utf-8

"""Module of class Pipe.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/piping/pipeline.py

SPDX-License-Identifier: MIT
"""
import math

from tespy.components.heat_exchangers.simple import SimpleHeatExchanger
from tespy.components.component import component_registry
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.data_containers import SimpleDataContainer as dc_simple

from tespy.tools.helpers import convert_to_SI
from tespy.tools.fluid_properties.wrappers import CoolPropWrapper


@component_registry
class Pipeline(SimpleHeatExchanger):
    r"""
    Pipeline for modeling steam networks.

    **Mandatory Equations**

    - :py:meth:`tespy.components.component.Component.fluid_func`
    - :py:meth:`tespy.components.component.Component.mass_flow_func`

    **Optional Equations**

    - :py:meth:`tespy.components.component.Component.pr_func`
    - :py:meth:`tespy.components.component.Component.zeta_func`
    - :py:meth:`tespy.components.component.Component.dp_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.energy_balance_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.darcy_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.hazen_williams_func`
    - :py:meth:`tespy.components.piping.pipeline.Pipeline.ohc_group_func`

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

    Tamb : float, dict
        Ambient temperature, provide parameter in network's temperature unit.

    insulation_thickness: float
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
    >>> from tespy.components import Sink, Source, Pipeline
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> nw = Network()
    >>> nw.set_attr(
    ...    p_unit='bar', T_unit='C', h_unit='kJ / kg', m_unit = 'kg / s',
    ...    iterinfo=False
    ... )
    >>> so = Source('source 1')
    >>> si = Sink('sink 1')
    >>> pi = Pipeline('pipeline')
    >>> inc = Connection(so, 'out1', pi, 'in1')
    >>> outg = Connection(pi, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)
    >>> inc.set_attr(fluid={'water': 1}, m=100, T=300, p=40)
    >>> pi.set_attr(
    ...     pr=0.95, Tamb=20, L=1000, D='var', ks=4.57e-5,
    ...     insulation_thickness=0.1, insulation_tc=0.035, pipe_thickness=0.003,
    ...     material='Steel', environment_media='air', wind_velocity=0.1
    ... )
    >>> nw.solve('design')
    >>> pi.D.val
    """

    @staticmethod
    def component():
        return 'Pipeline'

    def preprocess(self, num_nw_vars):
        self.air= CoolPropWrapper('air')

        super().preprocess(num_nw_vars)

    def get_parameters(self):
        parameters=super().get_parameters()
        # remove unused parameters
        for k in [ "kA_char_group", "kA_char"]:
            del parameters[k]

        parameters['Q_ohc_group_surface']=dc_gcp(
            elements=['insulation_thickness', 'insulation_tc', 'Tamb', 'material', 'pipe_thickness', 'environment_media','wind_velocity'],
            num_eq=1,
            func=self.ohc_surface_group_func,
            deriv=self.ohc_surface_group_deriv
        )
        parameters['Q_ohc_group_subsurface']=dc_gcp(
            elements=['insulation_thickness', 'insulation_tc', 'Tamb', 'material', 'pipe_thickness', 'environment_media','pipe_depth'],
            num_eq=1,
            func=self.ohc_subsurface_group_func,
            deriv=self.ohc_subsurface_group_deriv
        )
        parameters['insulation_thickness']=dc_cp(min_val=1e-3, max_val=1e1)
        parameters['insulation_tc']=dc_cp(min_val=1e-3, max_val=1e2)
        parameters['material']=dc_simple(val='Steel')
        parameters['pipe_thickness']=dc_cp(min_val=0, max_val=1)
        parameters['environment_media']=dc_simple(val='air')
        parameters['wind_velocity']=dc_cp(min_val=1e-3, max_val=10)
        parameters['pipe_depth']= dc_cp(min_val=1e-2, max_val=1e2)
        return parameters

    def ohc_surface_group_func(self):
        r"""Heat transfer calculation based on pipe material, insulation and
        surrounding ambient conditions fur surface pipes.

        Returns
        -------
        float
            Residual value of equation

            .. math::

                0 = \dot m \cdot \left(h_\text{out}-h_\text{in}\right)-
                \Delta T_\text{log} \cdot A \cdot U

                U = \frac{1}{\frac{1}{\alpha_\text{inner}} +
                R_\text{conductance} + \frac{1}{\alpha_\text{outer}}}

                \Delta T_{log} =
                \frac{T_{in}-T_{out}}{\ln{\frac{T_{in}-T_{amb}}
                {T_{out}-T_{amb}}}}
        
        Reference: :cite:`gnielinski1975`
        """

        diameters= [
            self.D.val,
            self.D.val + 2 * self.pipe_thickness.val,
            self.D.val + 2 * self.pipe_thickness.val + 2 * self.insulation_thickness.val
        ]

        # outer surface area per definition
        area = self.L.val * math.pi * diameters[2]

        # heat transfer resistance
        R_sum = []

        '''
        inner heat transfer resistance neglected yet
        R_int = 1/alpha_i *Diameters[2]/ Diameters[0]
        R_sum.append(R_int)
        '''

        # pipe wall heat transfer resistance
        pipe_tc ={'Steel':46.5, 'Carbon Steel':46, 'Cast Iron':48.8, 'Stainless Steel':21, 'PVC':0.23, 'Copper': 380}
        if diameters[1] > diameters[0]:
            if isinstance(self.material.val, str):
                wall_conductivity = pipe_tc[self.material.val]
            else:
                wall_conductivity = self.material.val
            R_sum.append(
                diameters[1] / wall_conductivity
                * math.log(diameters[1] / diameters[0]) / 2
            )

        # insulation heat transfer resistance
        if self.insulation_thickness.val != 0:
            R_sum.append(
                diameters[2] / self.insulation_tc.val
                * math.log(diameters[2] / diameters[1]) / 2
            )
        # external heat transfer resistance (to environment)
        Re = (
            self.wind_velocity.val * math.pi / 2
            * (diameters[1] + self.insulation_thickness.val * 2)
            / self.air.viscosity_pT(101300, self.Tamb.val_SI)
            * self.air.d_pT(101300, self.Tamb.val_SI)
        )
        Pr = self.air.AS.Prandtl()
        Nu_lam = 0.664 * Re ** 0.5 *Pr ** (1 / 3)
        Nu_turb = (
            0.037 * Re ** 0.8 * Pr
            / (1+ 2.443 * Re** (-0.1) * (Pr ** (2 / 3) - 1))
        )
        Nu_ext = 0.3 + (Nu_lam ** 2 + Nu_turb ** 2) ** 0.5
        alpha_ext = (
            Nu_ext
            / (math.pi / 2 * (diameters[1] + self.insulation_thickness.val *2))
            * self.air.AS.conductivity()
        ) #W/m²/K
        R_sum.append(1 / alpha_ext)

        if len(R_sum) == 0:
            raise ValueError("No heat transfer resistance. Check input values.")
      
        i = self.inl[0]
        o = self.outl[0]

        return (i.m.val_SI * (o.h.val_SI - i.h.val_SI) 
                + area / sum(R_sum) * self._deltaT_log()
        )

    def ohc_subsurface_group_func(self):
        r"""Heat transfer calculation based on pipe material, insulation and
        surrounding ambient conditions for subsurface pipes.

        Returns
        -------
        float
            Residual value of equation

            .. math::

                0 = \dot m \cdot \left(h_\text{out}-h_\text{in}\right)-
                \Delta T_\text{log} \cdot A \cdot U

                U = \frac{1}{\frac{1}{\alpha_\text{inner}} +
                R_\text{conductance} + \frac{1}{\alpha_\text{outer}}}

                \Delta T_{log} =
                \frac{T_{in}-T_{out}}{\ln{\frac{T_{in}-T_{amb}}
                {T_{out}-T_{amb}}}}
        
                First order approximation of multipole method for a single pipe in the ground.
        
        Assume no surface resistance.

        Reference: :cite:`wallenten1991`
        """

        diameters= [
            self.D.val,
            self.D.val + 2 * self.pipe_thickness.val,
            self.D.val + 2 * self.pipe_thickness.val + 2 * self.insulation_thickness.val
        ]

        '''
        inner heat transfer resistance neglected yet
        R_int = 1/alpha_i *Diameters[2]/ Diameters[0]
        R_sum.append(R_int)
        '''

        # external heat transfer resistance (to environment)
        ground_conductivity ={
            'gravel': 1.1, 'stones': 1.95, 'dry soil': 0.5, 'moist soil': 2.2
        }

        Beta = (
            ground_conductivity[self.environment_media.val]
            / self.insulation_tc.val * math.log(diameters[2] / diameters[0])
        )
        _h = (
            math.log(2 * self.pipe_depth.val / diameters[0]) + Beta
            + 1 / (
                1 - (2 * self.pipe_depth.val / diameters[0]) ** 2
                * (1 + Beta) / (1 - Beta)
            )
        )
        R_soil = (
            _h / (2 * math.pi * ground_conductivity[self.environment_media.val])
        )
        i = self.inl[0]
        o = self.outl[0]

        return (
            i.m.val_SI * (o.h.val_SI - i.h.val_SI) + 
            1 / R_soil * self._deltaT_log()
        )

    def ohc_subsurface_group_deriv(self, increment_filter, k):
        """Calculate the partial derivatives of the ohc equation

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        i = self.inl[0]
        o = self.outl[0]
        func= self.ohc_subsurface_group_func
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
        if self.D.is_var:
            self.jacobian[k, self.D.J_col] = self.numeric_deriv(func, 'D', None)

    def ohc_surface_group_deriv(self, increment_filter, k):
        """Calculate the partial derivatives of the ohc equation

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        i = self.inl[0]
        o = self.outl[0]
        func= self.ohc_surface_group_func
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
        if self.D.is_var:
            self.jacobian[k, self.D.J_col] = self.numeric_deriv(func, 'D', None)

