# -*- coding: utf-8

"""Module of class Pipe.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/piping/pipe.py

SPDX-License-Identifier: MIT
"""

import math

from tespy.components.heat_exchangers.simple import SimpleHeatExchanger
from tespy.components.component import component_registry
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.data_containers import SimpleDataContainer as dc_simple

from tespy.tools.fluid_properties.wrappers import CoolPropWrapper


@component_registry
class Pipe(SimpleHeatExchanger):
    r"""
    The Pipe is a subclass of a SimpleHeatExchanger.

    There are two different types of pipes available: An at the surface and
    a subsurface buried pipe. The implementation is based on
    :cite:`gnielinski1975` (surface) and :cite:`wallenten1991` (subsurface).

    **Mandatory Equations**

    - :py:meth:`tespy.components.component.Component.fluid_func`
    - :py:meth:`tespy.components.component.Component.mass_flow_func`

    **Optional Equations**

    - :py:meth:`tespy.components.component.Component.pr_func`
    - :py:meth:`tespy.components.component.Component.zeta_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.energy_balance_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.darcy_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.hazen_williams_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.kA_group_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.kA_char_group_func`
    - :py:meth:`tespy.components.heat_exchangers.piping.Pipe.ohc_surface_group_func`
    - :py:meth:`tespy.components.heat_exchangers.piping.Pipe.ohc_subsurface_group_func`

    Inlets/Outlets

    - in1
    - out1

    Optional inlets/outlets

    - heat

    Image

    .. image:: /api/_images/Pipe.svg
       :alt: flowsheet of the pipe
       :align: center
       :class: only-light

    .. image:: /api/_images/Pipe_darkmode.svg
       :alt: flowsheet of the pipe
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

    Q : float, dict
        Heat transfer, :math:`Q/\text{W}`.

    pr : float, dict
        Outlet to inlet pressure ratio, :math:`pr/1`.

    zeta : float, dict
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
        Ambient temperature, provide parameter in network's temperature
        unit, :math:`Tamb/\text{K}`.

    kA_group : str, dict
        Parametergroup for heat transfer calculation from ambient temperature
        and area independent heat transfer coefficient kA.

    insulation_thickness: float
        thickness of insulation, :math:`insulation_thickness/\text{m}`.

    insulation_tc: float
        thermal conductivity insulation,
        :math:`insulation_tc/\frac{\text{W}}{\text{m}\text{K}}`.

    material: str, float
        material of pipe: 'Steel', 'Carbon Steel', 'Cast Iron',
        'Stainless Steel', 'PVC', 'CommercialCopper' or user-specified heat
        conductivity of material: float

    pipe_thickness: float
        thickness of pipe, :math:`pipe_thickness/\text{m}`.

    environment_media: str
        environment media around the pipe: 'air', 'gravel', 'stones',
        'dry soil', 'moist soil'.

    wind_velocity: float
        Mean velocity of the wind. Needs to be greater than zero,
        :math:`wind_velocity/\frac{\text{m}}{\text{s}}`.

    pipe_depth: float
        pipe depth in the ground, :math:`pipe_depth/\text{m}`

    Example
    -------
    A mass flow of 10 kg/s hot ethanol is transported in a pipeline. The pipe
    is considered adiabatic, in the first approach and has a length of 100
    meters. We can calculate the diameter required at a given pressure loss of
    2.5 %. After we determined the required diameter, we can predict pressure
    loss at a different mass flow through the pipeline. Afterwards heat losses
    can be calculated by defining insulation and environment parameters. The
    heat losses of a subsurface pipe can be compared to heat losses of a
    surface pipe.

    >>> from tespy.components import Sink, Source, Pipe
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> import os
    >>> nw = Network()
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg', iterinfo=False)
    >>> so = Source('source 1')
    >>> si = Sink('sink 1')
    >>> pi = Pipe('pipeline')
    >>> pi.set_attr(pr=0.975, Q=0, L=100, D='var', ks=5e-5)
    >>> inc = Connection(so, 'out1', pi, 'in1')
    >>> outg = Connection(pi, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)
    >>> inc.set_attr(fluid={'ethanol': 1}, m=10, T=30, p=3)
    >>> nw.solve('design')
    >>> round(pi.D.val, 3)
    0.119
    >>> outg.p.val / inc.p.val == pi.pr.val
    True
    >>> inc.set_attr(m=15)
    >>> pi.set_attr(pr=None)
    >>> pi.set_attr(D=pi.D.val)
    >>> nw.solve('design')
    >>> round(pi.pr.val, 2)
    0.94

    In the second section the example shows how to calcualte the heat losses of
    the pipe to the ambient considering insulation. For this, we will look at a
    pipe transporting hot water. Since we change the fluid, we should also give
    a reasonable guess value for the outflow connection of the pipe as the
    initial guess originates from the previos calculation using ethanol as
    fluid.

    >>> inc.set_attr(fluid={'water': 1, 'ethanol': 0}, T=100)
    >>> outg.set_attr(h0=300)
    >>> pi.set_attr(
    ...     D='var', Q=None, pr=0.975,
    ...     Tamb=0, environment_media='dry soil', pipe_depth=5,
    ...     insulation_thickness=0.1, insulation_tc=0.035,
    ...     pipe_thickness=0.003, material='Steel'
    ... )
    >>> nw.solve('design')
    >>> round(pi.Q.val, 2)
    -17.81

    We can reuse many of the given parameters of the pipe. By unsetting the
    pipe's depth and setting the environment media and wind velocity instead
    the analogous method for surface pipes is applied. Observe, how the
    overall heat loss increases.

    >>> pi.set_attr(
    ...     pipe_depth=None, environment_media='air', wind_velocity=2
    ... )
    >>> nw.solve('design')
    >>> round(pi.Q.val, 2)
    -2434.12
    """

    def _preprocess(self, row_idx):
        self.air = CoolPropWrapper('air')

        super()._preprocess(row_idx)

    def get_parameters(self):
        parameters=super().get_parameters()

        parameters['Q_ohc_group_surface']=dc_gcp(
            elements=[
                'insulation_thickness', 'insulation_tc', 'Tamb', 'material',
                'pipe_thickness', 'environment_media', 'wind_velocity'
            ],
            num_eq_sets=1,
            func=self.ohc_surface_group_func,
            dependents=self.ohc_surface_group_dependents
        )
        parameters['Q_ohc_group_subsurface']=dc_gcp(
            elements=[
                'insulation_thickness', 'insulation_tc', 'Tamb', 'material',
                'pipe_thickness', 'environment_media','pipe_depth'
            ],
            num_eq_sets=1,
            func=self.ohc_subsurface_group_func,
            dependents=self.ohc_subsurface_group_dependents
        )
        parameters['insulation_thickness']=dc_cp(min_val=1e-3, max_val=1e1)
        parameters['insulation_tc']=dc_cp(min_val=1e-3, max_val=1e2)
        parameters['material']=dc_simple(val='Steel')
        parameters['pipe_thickness']=dc_cp(min_val=0, max_val=1)
        parameters['environment_media']=dc_simple(val='soil')
        parameters['wind_velocity']=dc_cp(min_val=1e-10, max_val=20)
        parameters['pipe_depth']= dc_cp(min_val=1e-2, max_val=1e2)
        return parameters

    def ohc_surface_group_func(self):
        r"""Heat transfer calculation based on pipe material, insulation and
        surrounding ambient conditions fur surface pipes.
        Valid for forced convection.

        Returns
        -------
        float
            Residual value of equation

            .. math::

                0 = \dot m \cdot \left(h_\text{out}-h_\text{in}\right)-
                \Delta T_\text{log} \cdot A \cdot U

                U = \frac{1}{\frac{1}{\alpha_\text{inner}} +
                R_\text{conductance} + \frac{1}{\alpha_\text{outer}}}

                \alpha_\text{outer} = \frac{Nu_\text{l} \cdot \lambda}{l}

                Nu_\text{l}= 0.3 + \sqrt{Nu_\text{l, lam}^{2} +
                Nu_\text{l, turb}^{2}}

                Nu_\text{l, turb} = \frac{0.037 Re_l^{0.8} \cdot
                Pr}{1+2.443 \cdot Re_l^{-0.1}\cdot (Pr^{2/3}-1)}

                Nu_\text{l, lam} = 0.664 \sqrt{Re_l}\cdot \sqrt[3]{Pr}

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
        pipe_tc ={
            'Steel':46.5, 'Carbon Steel':46, 'Cast Iron':48.8,
            'Stainless Steel':21, 'PVC':0.23, 'Copper': 380
        }
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

        return (
            i.m.val_SI * (o.h.val_SI - i.h.val_SI)
            + area / sum(R_sum) * self._calculate_td_log()
        )

    def ohc_surface_group_dependents(self):
        return (
            [self.inl[0].m]
            + [var for c in self.inl + self.outl for var in [c.p, c.h]]
            + [self.D]
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
            i.m.val_SI * (o.h.val_SI - i.h.val_SI)
            + 1 / R_soil * self._calculate_td_log()
        )

    def ohc_subsurface_group_dependents(self):
        return (
            [self.inl[0].m]
            + [var for c in self.inl + self.outl for var in [c.p, c.h]]
            + [self.D]
        )
