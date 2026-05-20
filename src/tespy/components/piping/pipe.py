# -*- coding: utf-8

"""Module of class Pipe.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/piping/pipe.py

SPDX-License-Identifier: MIT
"""

import math

from tespy.components.component import component_registry
from tespy.components.heat_exchangers.simple import SimpleHeatExchanger
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
from tespy.tools.fluid_properties.wrappers import CoolPropWrapper
from tespy.tools.logger import logger


@component_registry
class Pipe(SimpleHeatExchanger):
    r"""
    The Pipe is a subclass of a SimpleHeatExchanger.

    There are two different types of pipes available: An at the surface and
    a subsurface buried pipe. The implementation is based on
    :cite:`gnielinski1975` (surface) and :cite:`wallenten1991` (subsurface).

    .. image:: /api/_images/components/Pipe.svg
       :alt: flowsheet of the pipe
       :align: center
       :class: only-light

    .. image:: /api/_images/components/Pipe_darkmode.svg
       :alt: flowsheet of the pipe
       :align: center
       :class: only-dark

    Ports
    -----

    Fluid inlets: in1

    Fluid outlets: out1

    Power inlets: heat

    Power outlets: heat

    Heat inlets: heat

    Heat outlets: heat

    Mandatory Equations
    -------------------

    - mass flow equality constraint(s): :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
    - fluid composition equality constraint(s): :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`

    When a power or heat connector is attached:

    - energy_connector_balance: :py:meth:`energy_connector_balance_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.energy_connector_balance_func>`

    Parameters
    ----------

    char_warnings : bool
        Ignore warnings on default characteristics usage for this component.

    D : float, dict, :code:`"var"`
        Diameter of channel. Quantity: :code:`length`. Can be set as a system
        variable by passing :code:`"var"` as its value.

    darcy_group : GroupedComponentProperties
        Darcy-Weißbach equation for pressure loss. Elements: :code:`L`,
        :code:`ks`, :code:`D`.
        Equation: :py:meth:`darcy_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.darcy_func>`.

    design : list
        List containing design parameters (stated as String).

    design_path : str
        Path to the components design case.

    dissipative : bool


    dp : float, dict
        Inlet to outlet absolute pressure change. Quantity:
        :code:`pressure_difference`.
        Equation: :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`.

    environment_media : str


    flow_speed : float, dict
        Flow speed at inlet of pipe. Quantity: :code:`speed`.

    flow_speed_group : GroupedComponentProperties
        Equation connecting volumetric flow, flow speed and diameter of pipe.
        Elements: :code:`D`, :code:`flow_speed`.
        Equation: :py:meth:`flow_speed_func <tespy.components.piping.pipe.Pipe.flow_speed_func>`.

    hw_group : GroupedComponentProperties
        Hazen-Williams equation for pressure loss. Elements: :code:`L`,
        :code:`ks_HW`, :code:`D`.
        Equation: :py:meth:`hazen_williams_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.hazen_williams_func>`.

    insulation_tc : float, dict
        Thermal conductivity of insulation. Quantity:
        :code:`thermal_conductivity`.

    insulation_thickness : float, dict
        Thickness of pipe insulation. Quantity: :code:`length`.

    kA : float, dict, :code:`"var"`
        Heat transfer coefficient considering ambient temperature. Quantity:
        :code:`heat_transfer_coefficient`. Can be set as a system variable by
        passing :code:`"var"` as its value.

    kA_char : tespy.tools.characteristics.CharLine, dict
        Heat transfer coefficient lookup table for offdesign.

    kA_char_group : GroupedComponentProperties
        Heat transfer from design heat transfer coefficient, modifier lookup
        table and ambient temperature. Elements: :code:`kA_char`, :code:`Tamb`.
        Equation: :py:meth:`kA_char_group_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.kA_char_group_func>`.

    kA_group : GroupedComponentProperties
        Equation for heat transfer based on ambient temperature and heat
        transfer coefficient. Elements: :code:`kA`, :code:`Tamb`.
        Equation: :py:meth:`kA_group_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.kA_group_func>`.

    ks : float, dict, :code:`"var"`
        Roughness of wall material. Quantity: :code:`length`. Can be set as a
        system variable by passing :code:`"var"` as its value.

    ks_HW : float, dict, :code:`"var"`
        Hazen-Williams roughness. Can be set as a system variable by passing
        :code:`"var"` as its value.

    L : float, dict, :code:`"var"`
        Length of channel. Quantity: :code:`length`. Can be set as a system
        variable by passing :code:`"var"` as its value.

    label : str
        The label of the component.

    local_design : bool
        Treat this component in design mode in an offdesign calculation.

    local_offdesign : bool
        Treat this component in offdesign mode in a design calculation.

    material : str


    offdesign : list
        List containing offdesign parameters (stated as String).

    pipe_depth : float, dict
        Depth of buried pipe. Quantity: :code:`length`.

    pipe_thickness : float, dict
        Wall thickness of pipe. Quantity: :code:`length`.

    power_connector_location : str


    pr : float, dict
        Outlet to inlet pressure ratio. Quantity: :code:`ratio`.
        Equation: :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`.

    printout : bool
        Include this component in the network's results printout.

    Q : float, dict
        Heat transfer. Quantity: :code:`heat`.
        Equation: :py:meth:`energy_balance_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.energy_balance_func>`.

    Q_ohc_group_subsurface : GroupedComponentProperties
        Equation for heat loss of buried pipes. Elements:
        :code:`insulation_thickness`, :code:`insulation_tc`, :code:`Tamb`,
        :code:`material`, :code:`pipe_thickness`, :code:`environment_media`,
        :code:`pipe_depth`.
        Equation: :py:meth:`ohc_subsurface_group_func <tespy.components.piping.pipe.Pipe.ohc_subsurface_group_func>`.

    Q_ohc_group_surface : GroupedComponentProperties
        Equation for heat loss of surface pipes. Elements:
        :code:`insulation_thickness`, :code:`insulation_tc`, :code:`Tamb`,
        :code:`material`, :code:`pipe_thickness`, :code:`environment_media`,
        :code:`wind_velocity`.
        Equation: :py:meth:`ohc_surface_group_func <tespy.components.piping.pipe.Pipe.ohc_surface_group_func>`.

    Tamb : float, dict
        Ambient temperature. Quantity: :code:`temperature`.

    wind_velocity : float, dict
        Velocity of wind at insulation surface. Quantity: :code:`speed`.

    zeta : float, dict
        Non-dimensional friction coefficient for pressure loss calculation.
        Equation: :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`.

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
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "pressure_difference": "bar",
    ...     "temperature": "degC", "enthalpy": "kJ/kg"
    ... })
    >>> so = Source("source 1")
    >>> si = Sink("sink 1")
    >>> pi = Pipe("pipeline")
    >>> pi.set_attr(pr=0.975, Q=0, L=100, D="var", ks=5e-5)
    >>> inc = Connection(so, "out1", pi, "in1")
    >>> outg = Connection(pi, "out1", si, "in1")
    >>> nw.add_conns(inc, outg)
    >>> inc.set_attr(fluid={"ethanol": 1}, m=10, T=30, p=3)
    >>> nw.solve("design")
    >>> round(pi.D.val, 3)
    0.119
    >>> round(outg.p.val / inc.p.val, 3) == round(pi.pr.val, 3)
    True
    >>> inc.set_attr(m=15)
    >>> pi.set_attr(pr=None)
    >>> pi.set_attr(D=pi.D.val)
    >>> nw.solve("design")
    >>> round(pi.pr.val, 2)
    0.94

    In the second section the example shows how to calculate the heat losses of
    the pipe to the ambient considering insulation. For this, we will look at a
    pipe transporting hot water. Since we change the fluid, we should also give
    a reasonable guess value for the outflow connection of the pipe as the
    initial guess originates from the previous calculation using ethanol as
    fluid.

    >>> inc.set_attr(fluid={"water": 1, "ethanol": 0}, T=100)
    >>> outg.set_attr(h0=300)
    >>> pi.set_attr(
    ...     D="var", Q=None, pr=0.975,
    ...     Tamb=0, environment_media="dry soil", pipe_depth=5,
    ...     insulation_thickness=0.1, insulation_tc=0.035,
    ...     pipe_thickness=0.003, material="Steel"
    ... )
    >>> nw.solve("design")
    >>> round(pi.Q.val, 2)
    -1780.74

    We can reuse many of the given parameters of the pipe. By unsetting the
    pipe"s depth and setting the environment media and wind velocity instead
    the analogous method for surface pipes is applied. Observe, how the
    overall heat loss increases.

    >>> pi.Q_ohc_group_subsurface.is_set = False
    >>> pi.set_attr(
    ...     pipe_depth=None, environment_media="air", wind_velocity=2.0
    ... )
    >>> nw.solve("design")
    >>> round(pi.Q.val, 2)
    -2434.12
    """

    def _preprocess(self, row_idx):
        self.air = CoolPropWrapper("air")
        if self.wind_velocity.is_set:
            if self.wind_velocity.val < self.wind_velocity.min_val:
                msg = (
                    f"Minimum wind velocity is {self.wind_velocity.min_val} "
                    "for numerical reasons. The value is changed to the "
                    "specified minimum."
                )
                logger.debug(msg)
                self.wind_velocity.val = self.wind_velocity.min_val

        super()._preprocess(row_idx)

    def get_parameters(self):
        parameters=super().get_parameters()

        parameters["Q_ohc_group_surface"]=dc_gcp(
            elements=[
                "insulation_thickness", "insulation_tc", "Tamb", "material",
                "pipe_thickness", "environment_media", "wind_velocity"
            ],
            num_eq_sets=1,
            func=self.ohc_surface_group_func,
            dependents=self.ohc_surface_group_dependents,
            description="equation for heat loss of surface pipes"
        )
        parameters["Q_ohc_group_subsurface"]=dc_gcp(
            elements=[
                "insulation_thickness", "insulation_tc", "Tamb", "material",
                "pipe_thickness", "environment_media","pipe_depth"
            ],
            num_eq_sets=1,
            func=self.ohc_subsurface_group_func,
            dependents=self.ohc_subsurface_group_dependents,
            description="equation for heat loss of buried pipes"
        )
        parameters["insulation_thickness"]=dc_cp(
            min_val=1e-3, max_val=1e1, quantity="length",
            description="thickness of pipe insulation"
        )
        parameters["insulation_tc"]=dc_cp(
            min_val=1e-3, max_val=1e2, quantity="thermal_conductivity",
            description="thermal conductivity of insulation"
        )
        parameters["material"]=dc_simple(val="Steel", dtype="str")
        parameters["pipe_thickness"]=dc_cp(
            min_val=0, max_val=1, quantity="length",
            description="wall thickness of pipe"
        )
        parameters["environment_media"]=dc_simple(val="soil", dtype="str")
        parameters["wind_velocity"]=dc_cp(
            min_val=1e-6, max_val=20, quantity="speed",
            description="velocity of wind at insulation surface"
        )
        parameters["pipe_depth"]= dc_cp(
            min_val=1e-2, max_val=1e2, quantity="length",
            description="depth of buried pipe"
        )
        parameters["flow_speed"]= dc_cp(
            min_val=1e-2, max_val=1e2, quantity="speed",
            description="flow speed at inlet of pipe",
            calc=self._calc_flow_speed
        )
        parameters["flow_speed_group"]= dc_gcp(
            elements=["D", "flow_speed"],
            num_eq_sets=1,
            func=self.flow_speed_func,
            dependents=self.flow_speed_dependents,
            description="equation connecting volumetric flow, flow speed and diameter of pipe"
        )
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

                U = R_\text{conductance} + \frac{1}{\alpha_\text{outer}}

                \alpha_\text{outer} = \frac{Nu_\text{l} \cdot \lambda}{l}

                Nu_\text{l}= 0.3 + \sqrt{Nu_\text{l, lam}^{2} +
                Nu_\text{l, turb}^{2}}

                Nu_\text{l, turb} = \frac{0.037 Re_l^{0.8} \cdot
                Pr}{1+2.443 \cdot Re_l^{-0.1}\cdot (Pr^{2/3}-1)}

                Nu_\text{l, lam} = 0.664 \sqrt{Re_l}\cdot \sqrt[3]{Pr}

        Reference: :cite:`gnielinski1975`
        """

        diameters= [
            self.D.val_SI,
            self.D.val_SI + 2 * self.pipe_thickness.val_SI,
            self.D.val_SI + 2 * self.pipe_thickness.val_SI + 2 * self.insulation_thickness.val_SI
        ]

        # outer surface area per definition
        area = self.L.val_SI * math.pi * diameters[2]

        # heat transfer resistance
        R_sum = []

        """
        inner heat transfer resistance neglected yet
        R_int = 1/alpha_i *Diameters[2]/ Diameters[0]
        R_sum.append(R_int)
        """

        # pipe wall heat transfer resistance
        pipe_tc ={
            "Steel": 46.5, "Carbon Steel": 46, "Cast Iron": 48.8,
            "Stainless Steel": 21, "PVC": 0.23, "Copper": 380
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
        if self.insulation_thickness.val_SI != 0:
            R_sum.append(
                diameters[2] / self.insulation_tc.val_SI
                * math.log(diameters[2] / diameters[1]) / 2
            )
        # external heat transfer resistance (to environment)
        Re = (
            self.wind_velocity.val_SI * math.pi / 2
            * (diameters[1] + self.insulation_thickness.val_SI * 2)
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
            / (math.pi / 2 * (diameters[1] + self.insulation_thickness.val_SI *2))
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
            + [self.D, self.L]
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

            First order approximation of multipole method for a single pipe in
            the ground.


        Reference: :cite:`wallenten1991`
        """

        diameters= [
            self.D.val_SI,
            self.D.val_SI + 2 * self.pipe_thickness.val_SI,
            self.D.val_SI + 2 * self.pipe_thickness.val_SI
            + 2 * self.insulation_thickness.val_SI
        ]

        """
        inner heat transfer resistance neglected yet
        R_int = 1/alpha_i *Diameters[2]/ Diameters[0]
        R_sum.append(R_int)
        """

        # external heat transfer resistance (to environment)
        ground_conductivity ={
            "gravel": 1.1, "stones": 1.95, "dry soil": 0.5, "moist soil": 2.2
        }
        # conductivity of the pipe neglected according to the original publication
        Beta = (
            ground_conductivity[self.environment_media.val]
            / self.insulation_tc.val_SI * math.log(diameters[2] / diameters[0])
        )
        _h = (
            math.log(2 * self.pipe_depth.val_SI / diameters[0]) + Beta
            + 1 / (
                1 - (2 * self.pipe_depth.val_SI / diameters[0]) ** 2
                * (1 + Beta) / (1 - Beta)
            )
        )
        R_soil = (
            _h / (2 * math.pi * ground_conductivity[self.environment_media.val])
        )
        i = self.inl[0]
        o = self.outl[0]

        # here the resistance is a per meter of pipe resistance, therefore
        # we only multiply be length
        return (
            i.m.val_SI * (o.h.val_SI - i.h.val_SI)
            + 1 / R_soil * self._calculate_td_log() * self.L.val_SI
        )

    def ohc_subsurface_group_dependents(self):
        return (
            [self.inl[0].m]
            + [var for c in self.inl + self.outl for var in [c.p, c.h]]
            + [self.D, self.L]
        )

    def flow_speed_func(self):
        r"""Heat transfer calculation based on pipe material, insulation and
        surrounding ambient conditions for subsurface pipes.

        Returns
        -------
        float
            Residual value of equation

            .. math::

                0 = c \cdot \pi \cdot D ^ 2 - \dot m \cdot v_\text{in} * 4
        """
        return (
            self.flow_speed.val_SI * math.pi * self.D.val_SI ** 2
            - self.inl[0].m.val_SI * self.inl[0].calc_vol() * 4
        )

    def flow_speed_dependents(self):
        return [self.inl[0].m, self.inl[0].p, self.inl[0].h, self.D]

    def _calc_flow_speed(self):
        if not (self.D.is_set or self.D.is_var):
            return self.flow_speed.val_SI
        return self.inl[0].v.val_SI * 4 / (math.pi * self.D.val_SI ** 2)
