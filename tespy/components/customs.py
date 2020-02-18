# -*- coding: utf-8

"""Module for custom components.

Components in this module:

    - :func:`tespy.components.customs.orc_evaporator`

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/customs.py

SPDX-License-Identifier: MIT
"""

import logging

import numpy as np

from tespy.components.components import component

from tespy.tools.data_containers import dc_cc, dc_cp, dc_simple, dc_gcp
from tespy.tools.fluid_properties import (
        h_mix_pT, s_mix_ph, v_mix_ph, visc_mix_ph, T_mix_ph,
        dh_mix_dpQ, h_mix_pQ, T_bp_p, memorise
        )
from tespy.tools.helpers import lamb, single_fluid

# %%


class orc_evaporator(component):
    r"""
    Class orc_evaporator is the evaporator component in
    the Organic Rankine Cycle (ORC). Generally, the hot side
    of the geo-fluid from the geothermal wells keeps 2-phase.
    In order to fully use the energy in the geo-fluid,
    there are 2 inlets at the hot side.

    The ORC evaporator represents counter current evaporators. Both, 2 hot
    and 1 cold side of the evaporator, are simulated.

    Equations

        **mandatory equations**

        - :func:`tespy.components.customs.orc_evaporator.fluid_func`
        - :func:`tespy.components.customs.orc_evaporator.mass_flow_func`

        - :func:`tespy.components.customs.orc_evaporator.energy_func`

        .. math::

            0 = p_{1,in} \cdot pr1 - p_{1,out}\\
            0 = p_{2,in} \cdot pr2 - p_{2,out}\\
            0 = p_{3,in} \cdot pr3 - p_{3,out}

        - hot side steam :func:`tespy.components.components.component.zeta_func`
        - hot side brine :func:`tespy.components.components.component.zeta_func`
        - worling fluid :func:`tespy.components.components.component.zeta_func`

        **mandatory equations at outlet of the steam
        from geothermal heat source side**

        .. math::

            0 = h_{1,out} - h\left(p, x=0 \right)\\
            x: \text{vapour mass fraction}

        **mandatory equations at outlet of the working fluid
        of being evaporated side**

        .. math::

            0 = h_{3,out} - h\left(p, x=1 \right)\\
            x: \text{vapour mass fraction}

    Inlets/Outlets

        - in1, in2, in3 (index 1: steam from geothermal heat source,
        index 2: brine from geothermal heat source,
        index 3: working fluid of being evaporated)
        - out1, out2, out3 (index 1: steam from geothermal heat source,
        index 2: brine from geothermal heat source,
        index 3: working fluid of being evaporated)

    Image

        .. image:: _images/orc_evaporator.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    design_path: str
        Path to the components design case.

    local_offdesign : boolean
        Treat this component in offdesign mode in a design calculation.

    local_design : boolean
        Treat this component in design mode in an offdesign calculation.

    char_warnings: boolean
        Ignore warnings on default characteristics usage for this component.

    printout: boolean
        Include this component in the network's results printout.

    Q : String/float/tespy.helpers.dc_cp
        Heat transfer, :math:`Q/\text{W}`.

    pr1 : String/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio at hot side 1 (steam),
        :math:`pr/1`.

    pr2 : String/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio at hot side 2 (brine),
        :math:`pr/1`.

    pr3 : String/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio at cold side (working fluid),
        :math:`pr/1`.

    zeta1 : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient at hot side 1 (steam),
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    zeta2 : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient at hot side 2 (brine),
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    zeta3 : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient at cold side (working fluid),
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    subcooling : bool
        Enable/disable subcooling at oulet of the hot side 1,
        default value: disabled.

    overheating : bool
        Enable/disable overheating at oulet of the cold side,
        default value: disabled.

    Note
    ----
    The ORC evaporator has an additional equation for enthalpy
    at outlet of the steam from geothermal heat source side:
    The fluid leaves the component in saturated liquid state.
    If subcooling is activated, it possible to specify
    the enthalpy at the outgoing connection manually.

    It also has an another additional equation for enthalpy
    at outlet of the working fluid of being evaporated:
    The fluid leaves the component in saturated gas state.
    If overheating is activated, it possible to specify
    the enthalpy at the outgoing connection manually.

    Example
    -------
    A 2-phase geo-fluid is used as the heat source for evaporating
    the working fluid. The evaporator is designed for calculate the
    mass flow rate of the working fluid with known steam and
    brine mass flow rate. From this, it is possible to calculate
    the mass flow rate of the working fluid that is fully evaporated
    through the ORC evaporator.

    >>> from tespy.connections import connection
    >>> from tespy.networks import network
    >>> from tespy.components import source, sink
    >>> from tespy.components.customs import orc_evaporator
    >>> fluids = ['water', 'Isopentane']
    >>> nw = network(fluids=fluids, iterinfo=False)
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')
    >>> evaporator = orc_evaporator('geothermal orc evaporator')
    >>> evaporator.component()
    'orc_evaporator'
    >>> source_wf = source('working fluid source')
    >>> sink_wf = sink('working fluid sink')
    >>> source_s = source('steam source')
    >>> source_b = source('brine source')
    >>> sink_s = sink('steam sink')
    >>> sink_b = sink('brine sink')
    >>> eva_wf_in = connection(source_wf, 'out1', evaporator, 'in3')
    >>> eva_wf_out = connection(evaporator, 'out3', sink_wf, 'in1')
    >>> eva_steam_in = connection(source_s, 'out1', evaporator, 'in1')
    >>> eva_sink_s = connection(evaporator, 'out1', sink_s, 'in1')
    >>> eva_brine_in = connection(source_b, 'out1', evaporator, 'in2')
    >>> eva_sink_b = connection(evaporator, 'out2', sink_b, 'in1')
    >>> nw.add_conns(eva_wf_in, eva_wf_out)
    >>> nw.add_conns(eva_steam_in, eva_sink_s)
    >>> nw.add_conns(eva_brine_in, eva_sink_b)

    The orc working fluids leaves the evaporator in saturated steam state, the
    geothermal steam leaves the component in staturated liquid state. We imply
    the state of geothermal steam and brine with the corresponding mass flow as
    well as the working fluid's state at the evaporator inlet. Additionaly, the
    pressure ratios for all three streams are specified.

    >>> evaporator.set_attr(pr1=0.95, pr2=0.98, pr3=0.99)
    >>> eva_wf_in.set_attr(T=111, p=11,
    ... fluid={'water': 0, 'Isopentane': 1})
    >>> eva_steam_in.set_attr(T=147, p=4.3, m=20,
    ... fluid={'water': 1, 'Isopentane': 0})
    >>> eva_brine_in.set_attr(T=147, p=10.2, m=190,
    ... fluid={'water': 1, 'Isopentane': 0})
    >>> eva_sink_b.set_attr(T=117)
    >>> nw.solve(mode='design')

    Check the state of the steam and working fluid outlet:

    >>> eva_wf_out.x.val
    1.0
    >>> eva_sink_s.x.val
    0.0
    """

    @staticmethod
    def component():
        return 'orc_evaporator'

    @staticmethod
    def attr():
        return {'Q': dc_cp(max_val=0),
                'pr1': dc_cp(max_val=1), 'pr2': dc_cp(max_val=1),
                'pr3': dc_cp(max_val=1),
                'zeta1': dc_cp(min_val=0), 'zeta2': dc_cp(min_val=0),
                'zeta3': dc_cp(min_val=0),
                'subcooling': dc_simple(val=False),
                'overheating': dc_simple(val=False),
                'SQ1': dc_simple(), 'SQ2': dc_simple(), 'SQ3': dc_simple(),
                'Sirr': dc_simple()}

    @staticmethod
    def inlets():
        return ['in1', 'in2', 'in3']

    @staticmethod
    def outlets():
        return ['out1', 'out2', 'out3']

    def comp_init(self, nw):

        component.comp_init(self, nw)

        # number of mandatroy equations for
        # fluid balance: num_fl * 3
        # mass flow: 3
        # energy balance: 1
        self.num_eq = self.num_nw_fluids * 3 + 3 + 1
        # enthalpy hot side 1 outlet (if not subcooling): 1
        if self.subcooling.val is False:
            self.num_eq += 1
        # enthalpy cold side outlet (if not overheating): 1
        if self.overheating.val is False:
            self.num_eq += 1
        for var in [self.Q, self.pr1, self.pr2, self.pr3,
                    self.zeta1, self.zeta2, self.zeta3, ]:
            if var.is_set is True:
                self.num_eq += 1

        self.mat_deriv = np.zeros((
            self.num_eq,
            self.num_i + self.num_o + self.num_vars,
            self.num_nw_vars))

        self.vec_res = np.zeros(self.num_eq)
        pos = self.num_nw_fluids * 3
        self.mat_deriv[0:pos] = self.fluid_deriv()
        self.mat_deriv[pos:pos + 3] = self.mass_flow_deriv()

    def equations(self):
        r"""
        Calculates vector vec_res with results of equations for this component.

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        k = 0

        ######################################################################
        # equations for fluid balance
        self.vec_res[k:k + self.num_nw_fluids * 3] = self.fluid_func()
        k += self.num_nw_fluids * 3

        ######################################################################
        # equations for mass flow balance
        self.vec_res[k:k + 3] = self.mass_flow_func()
        k += 3

        ######################################################################
        # equations for energy balance
        self.vec_res[k] = self.energy_func()
        k += 1

        ######################################################################
        # equations for specified heat transfer
        if self.Q.is_set:
            self.vec_res[k] = (
                self.inl[2].m.val_SI * (
                    self.outl[2].h.val_SI - self.inl[2].h.val_SI) - self.Q.val)
            k += 1

        ######################################################################
        # equations for specified pressure ratio at hot side 1
        if self.pr1.is_set:
            self.vec_res[k] = (
                    self.pr1.val * self.inl[0].p.val_SI -
                    self.outl[0].p.val_SI)
            k += 1

        ######################################################################
        # equations for specified pressure ratio at hot side 2
        if self.pr2.is_set:
            self.vec_res[k] = (
                    self.pr2.val * self.inl[1].p.val_SI -
                    self.outl[1].p.val_SI)
            k += 1

        ######################################################################
        # equations for specified pressure ratio at cold side
        if self.pr3.is_set:
            self.vec_res[k] = (
                    self.pr3.val * self.inl[2].p.val_SI -
                    self.outl[2].p.val_SI)
            k += 1

        ######################################################################
        # equations for specified zeta at hot side 1
        if self.zeta1.is_set:
            self.vec_res[k] = self.zeta_func(zeta='zeta1', conn=0)
            k += 1

        ######################################################################
        # equations for specified zeta at hot side 2
        if self.zeta2.is_set:
            self.vec_res[k] = self.zeta_func(zeta='zeta2', conn=0)
            k += 1

        ######################################################################
        # equations for specified zeta at cold side
        if self.zeta3.is_set:
            self.vec_res[k] = self.zeta_func(zeta='zeta3', conn=0)
            k += 1

        ######################################################################
        # equation for saturated liquid at hot side 1 outlet
        if self.subcooling.val is False:
            o1 = self.outl[0].to_flow()
            self.vec_res[k] = o1[2] - h_mix_pQ(o1, 0)
            k += 1

        ######################################################################
        # equation for saturated gas at cold side outlet
        if self.overheating.val is False:
            o3 = self.outl[2].to_flow()
            self.vec_res[k] = o3[2] - h_mix_pQ(o3, 1)
            k += 1

    def derivatives(self, vec_z):
        r"""
        Calculates matrix of partial derivatives for given equations.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """

        ######################################################################
        # derivatives fluid and mass balance are static
        k = self.num_nw_fluids * 3 + 3

        ######################################################################
        # derivatives for energy balance equation
        # mat_deriv += self.energy_deriv()
        for i in range(3):
            self.mat_deriv[k, i, 0] = (
                    self.outl[i].h.val_SI - self.inl[i].h.val_SI)
            self.mat_deriv[k, i, 2] = -self.inl[i].m.val_SI

        self.mat_deriv[k, 3, 2] = self.inl[0].m.val_SI
        self.mat_deriv[k, 4, 2] = self.inl[1].m.val_SI
        self.mat_deriv[k, 5, 2] = self.inl[2].m.val_SI
        k += 1

        ######################################################################
        # derivatives for specified heat transfer
        if self.Q.is_set:
            self.mat_deriv[k, 2, 0] = (
                    self.outl[2].h.val_SI - self.inl[2].h.val_SI)
            self.mat_deriv[k, 2, 2] = -self.inl[2].m.val_SI
            self.mat_deriv[k, 5, 2] = self.inl[2].m.val_SI
            k += 1

        ######################################################################
        # derivatives for specified pressure ratio at hot side 1
        if self.pr1.is_set:
            self.mat_deriv[k, 0, 1] = self.pr1.val
            self.mat_deriv[k, 3, 1] = -1
            k += 1

        ######################################################################
        # derivatives for specified pressure ratio at hot side 2
        if self.pr2.is_set:
            self.mat_deriv[k, 1, 1] = self.pr2.val
            self.mat_deriv[k, 4, 1] = -1
            k += 1

        ######################################################################
        # derivatives for specified pressure ratio at cold side
        if self.pr3.is_set:
            self.mat_deriv[k, 2, 1] = self.pr3.val
            self.mat_deriv[k, 5, 1] = -1
            k += 1

        ######################################################################
        # derivatives for specified zeta at hot side 1
        if self.zeta1.is_set:
            f = self.zeta_func
            if not vec_z[0, 0]:
                self.mat_deriv[k, 0, 0] = self.numeric_deriv(
                    f, 'm', 0, zeta='zeta1', conn=0)
            if not vec_z[0, 1]:
                self.mat_deriv[k, 0, 1] = self.numeric_deriv(
                    f, 'p', 0, zeta='zeta1', conn=0)
            if not vec_z[0, 2]:
                self.mat_deriv[k, 0, 2] = self.numeric_deriv(
                    f, 'h', 0, zeta='zeta1', conn=0)
            if not vec_z[3, 1]:
                self.mat_deriv[k, 3, 1] = self.numeric_deriv(
                    f, 'p', 3, zeta='zeta1', conn=0)
            if not vec_z[3, 2]:
                self.mat_deriv[k, 3, 2] = self.numeric_deriv(
                    f, 'h', 3, zeta='zeta1', conn=0)
            k += 1

        ######################################################################
        # derivatives for specified zeta at hot side 2
        if self.zeta2.is_set:
            f = self.zeta_func
            if not vec_z[1, 0]:
                self.mat_deriv[k, 1, 0] = self.numeric_deriv(
                    f, 'm', 1, zeta='zeta2', conn=1)
            if not vec_z[1, 1]:
                self.mat_deriv[k, 1, 1] = self.numeric_deriv(
                    f, 'p', 1, zeta='zeta2', conn=1)
            if not vec_z[1, 2]:
                self.mat_deriv[k, 1, 2] = self.numeric_deriv(
                    f, 'h', 1, zeta='zeta2', conn=1)
            if not vec_z[4, 1]:
                self.mat_deriv[k, 4, 1] = self.numeric_deriv(
                    f, 'p', 4, zeta='zeta2', conn=1)
            if not vec_z[4, 2]:
                self.mat_deriv[k, 4, 2] = self.numeric_deriv(
                    f, 'h', 4, zeta='zeta2', conn=1)
            k += 1

        ######################################################################
        # derivatives for specified zeta at cold side
        if self.zeta3.is_set:
            f = self.zeta_func
            if not vec_z[2, 0]:
                self.mat_deriv[k, 2, 0] = self.numeric_deriv(
                    f, 'm', 2, zeta='zeta3', conn=2)
            if not vec_z[2, 1]:
                self.mat_deriv[k, 2, 1] = self.numeric_deriv(
                    f, 'p', 2, zeta='zeta3', conn=2)
            if not vec_z[2, 2]:
                self.mat_deriv[k, 2, 2] = self.numeric_deriv(
                    f, 'h', 2, zeta='zeta3', conn=2)
            if not vec_z[5, 1]:
                self.mat_deriv[k, 5, 1] = self.numeric_deriv(
                    f, 'p', 5, zeta='zeta3', conn=2)
            if not vec_z[5, 2]:
                self.mat_deriv[k, 5, 2] = self.numeric_deriv(
                    f, 'h', 5, zeta='zeta3', conn=2)
            k += 1

        ######################################################################
        # derivatives for saturated liquid at hot side 1 outlet equation
        if self.subcooling.val is False:
            o1 = self.outl[0].to_flow()
            self.mat_deriv[k, 3, 1] = -dh_mix_dpQ(o1, 0)
            self.mat_deriv[k, 3, 2] = 1
            k += 1

        ######################################################################
        # derivatives for saturated gas at cold side outlet 3 equation
        if self.overheating.val is False:
            o3 = self.outl[2].to_flow()
            self.mat_deriv[k, 5, 1] = -dh_mix_dpQ(o3, 1)
            self.mat_deriv[k, 5, 2] = 1
            k += 1

    def fluid_func(self):
        r"""
        Calculates the vector of residual values for component's fluid balance
        equations.

        Returns
        -------
        vec_res : list
            Vector of residual values for component's fluid balance.

            .. math::

                0 = fluid_{i,in_{j}} - fluid_{i,out_{j}} \;
                \forall i \in \mathrm{fluid}, \; \forall j \in inlets/outlets
        """
        vec_res = []

        for i in range(self.num_i):
            for fluid, x in self.inl[i].fluid.val.items():
                vec_res += [x - self.outl[i].fluid.val[fluid]]
        return vec_res

    def mass_flow_func(self):
        r"""
        Calculates the residual value for component's mass flow balance
        equation.

        Returns
        -------
        vec_res : list
            Vector with residual value for component's mass flow balance.

            .. math::

                0 = \dot{m}_{in,i} - \dot{m}_{out,i} \;
                \forall i \in inlets/outlets
        """
        vec_res = []
        for i in range(self.num_i):
            vec_res += [self.inl[i].m.val_SI - self.outl[i].m.val_SI]
        return vec_res

    def fluid_deriv(self):
        r"""
        Calculates the partial derivatives for all fluid balance equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((self.num_nw_fluids * 3,
                          6 + self.num_vars,
                          self.num_nw_vars))
        deriv = np.zeros((self.num_nw_fluids * self.num_i,
                          2 * self.num_i,
                          self.num_nw_vars))
        for i in range(self.num_i):
            for j in range(self.num_nw_fluids):
                deriv[i * self.num_nw_fluids + j, i, j + 3] = 1
                deriv[i * self.num_nw_fluids + j, self.num_i + i, j + 3] = -1
        return deriv

    def mass_flow_deriv(self):
        r"""
        Calculates the partial derivatives for all mass flow balance equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the mass flow balance
            equations.
        """
        deriv = np.zeros((self.num_i, 2 * self.num_i, self.num_nw_vars))
        for i in range(self.num_i):
            deriv[i, i, 0] = 1
        for j in range(self.num_i):
            deriv[j, j + self.num_i, 0] = -1
        return deriv

    def energy_func(self):
        r"""
        Equation for heat exchanger energy balance.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_{1,in} \cdot \left(h_{1,out} - h_{1,in} \right) +
                \dot{m}_{2,in} \cdot \left(h_{2,out} - h_{2,in} \right) +
                \dot{m}_{3,in} \cdot \left(h_{3,out} - h_{3,in} \right)
        """

        return (self.inl[0].m.val_SI * (self.outl[0].h.val_SI -
                                        self.inl[0].h.val_SI) +
                self.inl[1].m.val_SI * (self.outl[1].h.val_SI -
                                        self.inl[1].h.val_SI) +
                self.inl[2].m.val_SI * (self.outl[2].h.val_SI -
                                        self.inl[2].h.val_SI))

    def bus_func(self, bus):
        r"""
        Calculates the residual value of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus
            TESPy bus object.

        Returns
        -------
        val : float
            Residual value of equation.

            .. math::

                val = P \cdot f\left( \frac{P}{P_{ref}}\right)

                P = \dot{m}_{3,in} \cdot \left( h_{3,out} - h_{3,in} \right)
        """
        i = self.inl[2].to_flow()
        o = self.outl[2].to_flow()

        val = i[0] * (o[2] - i[2])
        if np.isnan(bus.P_ref):
            expr = 1
        else:
            expr = abs(val / bus.P_ref)
        return val * bus.char.evaluate(expr)

    def bus_deriv(self, bus):
        r"""
        Calculates the matrix of partial derivatives of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus
            TESPy bus object.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 6, self.num_nw_vars))
        deriv[0, 2, 0] = self.numeric_deriv(self.bus_func, 'm', 2, bus=bus)
        deriv[0, 2, 2] = self.numeric_deriv(self.bus_func, 'h', 2, bus=bus)
        deriv[0, 5, 2] = self.numeric_deriv(self.bus_func, 'h', 5, bus=bus)
        return deriv

    def initialise_source(self, c, key):
        r"""
        Returns a starting value for pressure and enthalpy at component's
        outlet.

        Parameters
        ----------
        c : tespy.connections.connection
            Connection to perform initialisation on.

        key : str
            Fluid property to retrieve.

        Returns
        -------
        val : float
            Starting value for pressure/enthalpy in SI units.

            .. math::

                val = \begin{cases}
                4 \cdot 10^5 & \text{key = 'p'}\\
                h\left(p, 200 \text{K} \right) & \text{key = 'h' at outlet 1}\\
                h\left(p, 250 \text{K} \right) & \text{key = 'h' at outlet 2}
                \end{cases}
        """
        if key == 'p':
            return 50e5
        elif key == 'h':
            flow = [c.m.val0, c.p.val_SI, c.h.val_SI, c.fluid.val]
            if c.s_id == 'out1':
                T = 200 + 273.15
                return h_mix_pT(flow, T)
            elif c.s_id == 'out2':
                T = 200 + 273.15
                return h_mix_pT(flow, T)
            else:
                T = 250 + 273.15
                return h_mix_pT(flow, T)

    def initialise_target(self, c, key):
        r"""
        Returns a starting value for pressure and enthalpy at component's
        inlet.

        Parameters
        ----------
        c : tespy.connections.connection
            Connection to perform initialisation on.

        key : str
            Fluid property to retrieve.

        Returns
        -------
        val : float
            Starting value for pressure/enthalpy in SI units.

            .. math::

                val = \begin{cases}
                4 \cdot 10^5 & \text{key = 'p'}\\
                h\left(p, 300 \text{K} \right) & \text{key = 'h' at inlet 1}\\
                h\left(p, 220 \text{K} \right) & \text{key = 'h' at outlet 2}
                \end{cases}
        """
        if key == 'p':
            return 50e5
        elif key == 'h':
            flow = [c.m.val0, c.p.val_SI, c.h.val_SI, c.fluid.val]
            if c.t_id == 'in1':
                T = 300 + 273.15
                return h_mix_pT(flow, T)
            elif c.t_id == 'in2':
                T = 300 + 273.15
                return h_mix_pT(flow, T)
            else:
                T = 220 + 273.15
                return h_mix_pT(flow, T)

    def calc_parameters(self):
        r"""
        Postprocessing parameter calculation.
        """
        # connection information
        i1 = self.inl[0].to_flow()
        i2 = self.inl[1].to_flow()
        i3 = self.inl[2].to_flow()
        o1 = self.outl[0].to_flow()
        o2 = self.outl[1].to_flow()
        o3 = self.outl[2].to_flow()

        # temperatures
        if isinstance(self, orc_evaporator):
            T_i1 = T_bp_p(i1)
        else:
            T_i1 = T_mix_ph(i1, T0=self.inl[0].T.val_SI)
        T_i2 = T_mix_ph(i2, T0=self.inl[1].T.val_SI)
        T_i3 = T_mix_ph(i3, T0=self.inl[2].T.val_SI)
        T_o1 = T_mix_ph(o1, T0=self.outl[0].T.val_SI)
        T_o2 = T_mix_ph(o2, T0=self.outl[1].T.val_SI)
        T_o3 = T_mix_ph(o3, T0=self.outl[2].T.val_SI)

        # specific volume
        v_i1 = v_mix_ph(i1, T0=T_i1)
        v_i2 = v_mix_ph(i2, T0=T_i2)
        v_i3 = v_mix_ph(i3, T0=T_i3)
        v_o1 = v_mix_ph(o1, T0=T_o1)
        v_o2 = v_mix_ph(o2, T0=T_o2)
        v_o3 = v_mix_ph(o3, T0=T_o3)

        # specific entropy
        s_i1 = s_mix_ph(i1, T0=T_i1)
        s_i2 = s_mix_ph(i2, T0=T_i2)
        s_i3 = s_mix_ph(i3, T0=T_i3)
        s_o1 = s_mix_ph(o1, T0=T_o1)
        s_o2 = s_mix_ph(o2, T0=T_o2)
        s_o3 = s_mix_ph(o3, T0=T_o3)

        # component parameters
        self.Q.val = -i3[0] * (o3[2] - i3[2])

        self.pr1.val = o1[1] / i1[1]
        self.pr2.val = o2[1] / i2[1]
        self.pr3.val = o3[1] / i3[1]
        self.zeta1.val = ((i1[1] - o1[1]) * np.pi ** 2 /
                          (8 * i1[0] ** 2 * (v_i1 + v_o1) / 2))
        self.zeta2.val = ((i2[1] - o2[1]) * np.pi ** 2 /
                          (8 * i2[0] ** 2 * (v_i2 + v_o2) / 2))
        self.zeta3.val = ((i3[1] - o3[1]) * np.pi ** 2 /
                          (8 * i3[0] ** 2 * (v_i3 + v_o3) / 2))

        self.SQ1.val = self.inl[0].m.val_SI * (s_o1 - s_i1)
        self.SQ2.val = self.inl[1].m.val_SI * (s_o2 - s_i2)
        self.SQ3.val = self.inl[2].m.val_SI * (s_o3 - s_i3)
        self.Sirr.val = self.SQ1.val + self.SQ2.val + self.SQ3.val

        self.check_parameter_bounds()
