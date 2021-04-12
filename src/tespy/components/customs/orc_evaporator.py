# -*- coding: utf-8

"""Module of class ORCEvaporator.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/customs/orc_evaporator.py

SPDX-License-Identifier: MIT
"""

import numpy as np

from tespy.components.component import Component
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import DataContainerSimple as dc_simple
from tespy.tools.document_models import generate_latex_eq
from tespy.tools.fluid_properties import dh_mix_dpQ
from tespy.tools.fluid_properties import h_mix_pQ
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.fluid_properties import s_mix_ph


class ORCEvaporator(Component):
    r"""
    Evaporator of the geothermal Organic Rankine Cycle (ORC).

    Generally, the hot side of the geo-fluid from the geothermal wells deliver
    two phases: steam and brine. In order to fully use the energy of the
    geo-fluid, there are 2 inlets at the hot side.

    The ORC evaporator represents counter current evaporators. Both, two hot
    and one cold side of the evaporator, are simulated.

    **Mandatory Equations**

    - :py:meth:`tespy.components.component.Component.fluid_func`
    - :py:meth:`tespy.components.component.Component.mass_flow_func`
    - :py:meth:`tespy.components.customs.orc_evaporator.ORCEvaporator.energy_balance_func`
    - steam side outlet state, function can be disabled by specifying
      :code:`set_attr(subcooling=True)`
      :py:meth:`tespy.components.customs.orc_evaporator.ORCEvaporator.subcooling_func`
    - working fluid outlet state, function can be disabled by specifying
      :code:`set_attr(overheating=True)`
      :py:meth:`tespy.components.customs.orc_evaporator.ORCEvaporator.overheating_func`

    **Optional Equations**

    - :py:meth:`tespy.components.customs.orc_evaporator.ORCEvaporator.energy_balance_cold_func`
    - hot side steam :py:meth:`tespy.components.component.Component.pr_func`
    - hot side brine :py:meth:`tespy.components.component.Component.pr_func`
    - worling fluid :py:meth:`tespy.components.component.Component.pr_func`
    - hot side steam :py:meth:`tespy.components.component.Component.zeta_func`
    - hot side brine :py:meth:`tespy.components.component.Component.zeta_func`
    - worling fluid :py:meth:`tespy.components.component.Component.zeta_func`


    Inlets/Outlets

    - in1, in2, in3 (index 1: steam from geothermal heat source,
      index 2: brine from geothermal heat source,
      index 3: working fluid of being evaporated)
    - out1, out2, out3 (index 1: steam from geothermal heat source,
      index 2: brine from geothermal heat source,
      index 3: working fluid of being evaporated)

    Image

    .. image:: _images/ORCEvaporator.svg
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

    pr1 : float, dict, :code:`"var"`
        Outlet to inlet pressure ratio at hot side 1 (steam),
        :math:`pr/1`.

    pr2 : float, dict, :code:`"var"`
        Outlet to inlet pressure ratio at hot side 2 (brine),
        :math:`pr/1`.

    pr3 : float, dict, :code:`"var"`
        Outlet to inlet pressure ratio at cold side (working fluid),
        :math:`pr/1`.

    zeta1 : float, dict, :code:`"var"`
        Geometry independent friction coefficient at hot side 1 (steam),
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    zeta2 : float, dict, :code:`"var"`
        Geometry independent friction coefficient at hot side 2 (brine),
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    zeta3 : float, dict, :code:`"var"`
        Geometry independent friction coefficient at cold side (working fluid),
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    subcooling : boolean
        Enable/disable subcooling at oulet of the hot side 1,
        default value: disabled (False).

    overheating : boolean
        Enable/disable overheating at oulet of the cold side,
        default value: disabled (False).

    Note
    ----
    The ORC evaporator has an additional equation for enthalpy at the outlet of
    the geothermal steam: The fluid leaves the component in saturated liquid
    state. If code:`subcooling` is activated (:code:`True`), it is possible to
    specify the enthalpy at the outgoing connection manually.

    Additionally, an equation for enthalpy at the outlet of the working fluid
    is imposed: It leaves the component in saturated gas state. If
    :code:`overheating` is enabled (:code:`True`), it is possible to specify
    the enthalpy at the outgoing connection manually.

    Example
    -------
    A two-phase geo-fluid is used as the heat source for evaporating the
    working fluid. We calculate the mass flow of the working fluid with known
    steam and brine mass flow.

    >>> from tespy.components import Source, Sink, ORCEvaporator
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> fluids = ['water', 'Isopentane']
    >>> nw = Network(fluids=fluids, iterinfo=False)
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')
    >>> evaporator = ORCEvaporator('geothermal orc evaporator')
    >>> evaporator.component()
    'orc evaporator'
    >>> source_wf = Source('working fluid source')
    >>> sink_wf = Sink('working fluid sink')
    >>> source_s = Source('steam source')
    >>> source_b = Source('brine source')
    >>> sink_s = Sink('steam sink')
    >>> sink_b = Sink('brine sink')
    >>> eva_wf_in = Connection(source_wf, 'out1', evaporator, 'in3')
    >>> eva_wf_out = Connection(evaporator, 'out3', sink_wf, 'in1')
    >>> eva_steam_in = Connection(source_s, 'out1', evaporator, 'in1')
    >>> eva_sink_s = Connection(evaporator, 'out1', sink_s, 'in1')
    >>> eva_brine_in = Connection(source_b, 'out1', evaporator, 'in2')
    >>> eva_sink_b = Connection(evaporator, 'out2', sink_b, 'in1')
    >>> nw.add_conns(eva_wf_in, eva_wf_out)
    >>> nw.add_conns(eva_steam_in, eva_sink_s)
    >>> nw.add_conns(eva_brine_in, eva_sink_b)

    The orc working fluids leaves the evaporator in saturated steam state, the
    geothermal steam leaves the component in staturated liquid state. We imply
    the state of geothermal steam and brine with the corresponding mass flow as
    well as the working fluid's state at the evaporator inlet. The pressure
    ratio is specified for each of the three streams.

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
        return 'orc evaporator'

    def get_variables(self):
        return {
            'Q': dc_cp(
                max_val=0, num_eq=1, latex=self.energy_balance_cold_func_doc,
                func=self.energy_balance_cold_func,
                deriv=self.energy_balance_cold_deriv),
            'pr1': dc_cp(
                min_val=1e-4, max_val=1, num_eq=1, deriv=self.pr_deriv,
                latex=self.pr_func_doc,
                func=self.pr_func, func_params={'pr': 'pr1'}),
            'pr2': dc_cp(
                min_val=1e-4, max_val=1, num_eq=1, latex=self.pr_func_doc,
                deriv=self.pr_deriv, func=self.pr_func,
                func_params={'pr': 'pr2', 'inconn': 1, 'outconn': 1}),
            'pr3': dc_cp(
                min_val=1e-4, max_val=1, num_eq=1, latex=self.pr_func_doc,
                deriv=self.pr_deriv, func=self.pr_func,
                func_params={'pr': 'pr3', 'inconn': 2, 'outconn': 2}),
            'zeta1': dc_cp(
                min_val=0, max_val=1e15, num_eq=1, latex=self.zeta_func_doc,
                deriv=self.zeta_deriv, func=self.zeta_func,
                func_params={'zeta': 'zeta1'}),
            'zeta2': dc_cp(
                min_val=0, max_val=1e15, num_eq=1, latex=self.zeta_func_doc,
                deriv=self.zeta_deriv, func=self.zeta_func,
                func_params={'zeta': 'zeta2', 'inconn': 1, 'outconn': 1}),
            'zeta3': dc_cp(
                min_val=0, max_val=1e15, num_eq=1, latex=self.zeta_func_doc,
                deriv=self.zeta_deriv, func=self.zeta_func,
                func_params={'zeta': 'zeta3', 'inconn': 2, 'outconn': 2}),
            'subcooling': dc_simple(
                val=False, num_eq=1, latex=self.subcooling_func_doc,
                deriv=self.subcooling_deriv, func=self.subcooling_func),
            'overheating': dc_simple(
                val=False, num_eq=1, latex=self.overheating_func_doc,
                deriv=self.overheating_deriv, func=self.overheating_func)
        }

    def get_mandatory_constraints(self):
        return {
            'mass_flow_constraints': {
                'func': self.mass_flow_func, 'deriv': self.mass_flow_deriv,
                'constant_deriv': True, 'latex': self.mass_flow_func_doc,
                'num_eq': 3},
            'fluid_constraints': {
                'func': self.fluid_func, 'deriv': self.fluid_deriv,
                'constant_deriv': True, 'latex': self.fluid_func_doc,
                'num_eq': self.num_nw_fluids * 3},
            'energy_balance_constraints': {
                'func': self.energy_balance_func,
                'deriv': self.energy_balance_deriv,
                'constant_deriv': False, 'latex': self.energy_balance_func_doc,
                'num_eq': 1}
        }

    @staticmethod
    def inlets():
        return ['in1', 'in2', 'in3']

    @staticmethod
    def outlets():
        return ['out1', 'out2', 'out3']

    def comp_init(self, nw):

        self.overheating.is_set = not self.overheating.val
        self.subcooling.is_set = not self.subcooling.val
        Component.comp_init(self, nw)

    def energy_balance_func(self):
        r"""
        Equation for heat exchanger energy balance.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                \begin{split}
                0 = &
                \dot{m}_{in,1} \cdot \left(h_{out,1} - h_{in,1} \right) \\
                &+ \dot{m}_{in,2} \cdot \left(h_{out,2} - h_{in,2} \right) \\
                &+ \dot{m}_{in,3} \cdot \left(h_{out,3} - h_{in,3} \right)
                \end{split}
        """
        return (
            self.inl[0].m.val_SI * (
                self.outl[0].h.val_SI - self.inl[0].h.val_SI) +
            self.inl[1].m.val_SI * (
                self.outl[1].h.val_SI - self.inl[1].h.val_SI) +
            self.inl[2].m.val_SI * (
                self.outl[2].h.val_SI - self.inl[2].h.val_SI))

    def energy_balance_func_doc(self, label):
        r"""
        Equation for heat exchanger energy balance.

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
            r'\begin{split}' + '\n'
            r'0 = &' + '\n'
            r'\dot{m}_\mathrm{in,1}\cdot\left(h_\mathrm{out,1}-'
            r'h_\mathrm{in,1}\right) \\' + '\n'
            r'&+ \dot{m}_\mathrm{in,2} \cdot \left(h_\mathrm{out,2} - '
            r'h_\mathrm{in,2} \right)\\' + '\n'
            r'&+ \dot{m}_\mathrm{in,3} \cdot \left(h_\mathrm{out,3} - '
            r'h_\mathrm{in,3} \right)' + '\n'
            r'\end{split}')
        return generate_latex_eq(self, latex, latex)

    def energy_balance_deriv(self, increment_filter, k):
        """
        Calculate partial derivatives of energy balance function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        for i in range(3):
            self.jacobian[k, i, 0] = (
                self.outl[i].h.val_SI - self.inl[i].h.val_SI)
            self.jacobian[k, i, 2] = -self.inl[i].m.val_SI
            self.jacobian[k, i + 3, 2] = self.inl[i].m.val_SI
        k += 1

    def energy_balance_cold_func(self):
        r"""
        Equation for cold side heat exchanger energy balance.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 =\dot{m}_{in,3} \cdot \left(h_{out,3}-h_{in,3}\right)+\dot{Q}
        """
        return self.inl[2].m.val_SI * (
            self.outl[2].h.val_SI - self.inl[2].h.val_SI) + self.Q.val

    def energy_balance_cold_func_doc(self, label):
        r"""
        Equation for cold side heat exchanger energy balance.

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
            r'0 =\dot{m}_{in,3} \cdot \left(h_{out,3}-'
            r'h_{in,3}\right)+\dot{Q}')
        return [generate_latex_eq(self, latex, label)]

    def energy_balance_cold_deriv(self, increment_filter, k):
        """
        Partial derivatives for cold side energy balance.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        self.jacobian[k, 2, 0] = self.outl[2].h.val_SI - self.inl[2].h.val_SI
        self.jacobian[k, 2, 2] = -self.inl[2].m.val_SI
        self.jacobian[k, 5, 2] = self.inl[2].m.val_SI

    def subcooling_func(self):
        r"""
        Equation for steam side outlet state.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0=h_{out,1} -h\left(p_{out,1}, x=0 \right)

        Note
        ----
        This equation is applied in case subcooling is False!
        """
        return self.outl[0].h.val_SI - h_mix_pQ(self.outl[0].get_flow(), 0)

    def subcooling_func_doc(self, label):
        r"""
        Equation for steam side outlet state.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        latex = r'0=h_\mathrm{out,1} -h\left(p_\mathrm{out,1}, x=0 \right)'
        return generate_latex_eq(self, latex, label)

    def subcooling_deriv(self, increment_filter, k):
        """
        Calculate partial derivatives for steam side outlet state.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        self.jacobian[k, 3, 1] = -dh_mix_dpQ(self.outl[0].get_flow(), 0)
        self.jacobian[k, 3, 2] = 1

    def overheating_func(self):
        r"""
        Equation for cold side outlet state.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0=h_{out,3} -h\left(p_{out,3}, x=1 \right)

        Note
        ----
        This equation is applied in case overheating is False!
        """
        return self.outl[2].h.val_SI - h_mix_pQ(self.outl[2].get_flow(), 1)

    def overheating_func_doc(self, label):
        r"""
        Equation for cold side outlet state.

        Parameters
        ----------
        label : str
            Label for equation.
        """
        latex = r'0=h_\mathrm{out,3} -h\left(p_\mathrm{out,3}, x=1 \right)'
        return generate_latex_eq(self, latex, label)

    def overheating_deriv(self, increment_filter, k):
        """
        Calculate partial derivatives for cold side outlet state.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        self.jacobian[k, 5, 1] = -dh_mix_dpQ(self.outl[0].get_flow(), 0)
        self.jacobian[k, 5, 2] = 1

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

                \dot{E} = -\dot{m}_{in,3} \cdot \left(
                h_{out,3} - h_{in,3} \right)
        """
        return -self.inl[2].m.val_SI * (
            self.outl[2].h.val_SI - self.inl[2].h.val_SI)

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
            r'-\dot{m}_\mathrm{in,3} \cdot \left(h_\mathrm{out,3} - '
            r'h_\mathrm{in,3} \right)')

    def bus_deriv(self, bus):
        r"""
        Calculate the matrix of partial derivatives of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus.Bus
            TESPy bus object.

        Returns
        -------
        deriv : ndarray
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 6, self.num_nw_vars))
        f = self.calc_bus_value
        deriv[0, 2, 0] = self.numeric_deriv(f, 'm', 2, bus=bus)
        deriv[0, 2, 2] = self.numeric_deriv(f, 'h', 2, bus=bus)
        deriv[0, 5, 2] = self.numeric_deriv(f, 'h', 5, bus=bus)
        return deriv

    def initialise_source(self, c, key):
        r"""
        Return a starting value for pressure and enthalpy at outlet.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            Connection to perform initialisation on.

        key : str
            Fluid property to retrieve.

        Returns
        -------
        val : float
            Starting value for pressure/enthalpy in SI units.

            .. math::

                val = \begin{cases}
                10 \cdot 10^5 & \text{key = 'p'}\\
                h\left(p, 473.15 \text{K} \right) &
                \text{key = 'h' at outlet 1}\\
                h\left(p, 473.15 \text{K} \right) &
                \text{key = 'h' at outlet 2}\\
                h\left(p, 523.15 \text{K} \right) &
                \text{key = 'h' at outlet 3}
                \end{cases}
        """
        if key == 'p':
            return 10e5
        elif key == 'h':
            if c.source_id == 'out1':
                T = 200 + 273.15
                return h_mix_pT(c.get_flow(), T)
            elif c.source_id == 'out2':
                T = 200 + 273.15
                return h_mix_pT(c.get_flow(), T)
            else:
                T = 250 + 273.15
                return h_mix_pT(c.get_flow(), T)

    def initialise_target(self, c, key):
        r"""
        Return a starting value for pressure and enthalpy at inlet.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            Connection to perform initialisation on.

        key : str
            Fluid property to retrieve.

        Returns
        -------
        val : float
            Starting value for pressure/enthalpy in SI units.

            .. math::

                val = \begin{cases}
                10 \cdot 10^5 & \text{key = 'p'}\\
                h\left(p, 573.15 \text{K} \right) &
                \text{key = 'h' at inlet 1}\\
                h\left(p, 573.15 \text{K} \right) &
                \text{key = 'h' at inlet 2}\\
                h\left(p, 493.15 \text{K} \right) &
                \text{key = 'h' at inlet 3}
                \end{cases}
        """
        if key == 'p':
            return 10e5
        elif key == 'h':
            if c.target_id == 'in1':
                T = 300 + 273.15
                return h_mix_pT(c.get_flow(), T)
            elif c.target_id == 'in2':
                T = 300 + 273.15
                return h_mix_pT(c.get_flow(), T)
            else:
                T = 220 + 273.15
                return h_mix_pT(c.get_flow(), T)

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        # component parameters
        self.Q.val = -self.inl[2].m.val_SI * (
            self.outl[2].h.val_SI - self.inl[2].h.val_SI)
        # pressure ratios and zeta values
        for i in range(3):
            self.get_attr('pr' + str(i + 1)).val = (
                self.outl[i].p.val_SI / self.inl[i].p.val_SI)
            self.get_attr('zeta' + str(i + 1)).val = (
                (self.inl[i].p.val_SI - self.outl[i].p.val_SI) * np.pi ** 2 / (
                    4 * self.inl[i].m.val_SI ** 2 *
                    (self.inl[i].vol.val_SI + self.outl[i].vol.val_SI)
                ))

    def entropy_balance(self):
        r"""
        Calculate entropy balance of the two-phase orc evaporator.

        The allocation of the entropy streams due to heat exchanged and due to
        irreversibility is performed by solving for T on all sides of the heat
        exchanger:

        .. math::

            h_\mathrm{out} - h_\mathrm{in} = \int_\mathrm{in}^\mathrm{out} v
            \cdot dp - \int_\mathrm{in}^\mathrm{out} T \cdot ds

        As solving :math:`\int_\mathrm{in}^\mathrm{out} v \cdot dp` for non
        isobaric processes would require perfect process knowledge (the path)
        on how specific volume and pressure change throught the component, the
        heat transfer is splitted into three separate virtual processes:

        - in->in*: decrease pressure to
          :math:`p_\mathrm{in*}=p_\mathrm{in}\cdot\sqrt{\frac{p_\mathrm{out}}{p_\mathrm{in}}}`
          without changing enthalpy.
        - in*->out* transfer heat without changing pressure.
          :math:`h_\mathrm{out*}-h_\mathrm{in*}=h_\mathrm{out}-h_\mathrm{in}`
        - out*->out decrease pressure to outlet pressure :math:`p_\mathrm{out}`
          without changing enthalpy.

        Note
        ----
        The entropy balance makes the follwing parameter available:

        .. math::

            \text{S\_Q1}=\dot{m} \cdot \left(s_\mathrm{out*,1}-s_\mathrm{in*,1}
            \right)\\
            \text{S\_Q2}=\dot{m} \cdot \left(s_\mathrm{out*,2}-s_\mathrm{in*,2}
            \right)\\
            \text{S\_Q3}=\dot{m} \cdot \left(s_\mathrm{out*,3}-s_\mathrm{in*,3}
            \right)\\
            \text{S\_Qirr}=\text{S\_Q3} - \text{S\_Q1} - \text{S\_Q2}\\
            \text{S\_irr1}=\dot{m} \cdot \left(s_\mathrm{out,1}-s_\mathrm{in,1}
            \right) - \text{S\_Q1}\\
            \text{S\_irr2}=\dot{m} \cdot \left(s_\mathrm{out,2}-s_\mathrm{in,2}
            \right) - \text{S\_Q2}\\
            \text{S\_irr3}=\dot{m} \cdot \left(s_\mathrm{out,3}-s_\mathrm{in,3}
            \right) - \text{S\_Q3}\\
            \text{S\_irr}=\sum \dot{S}_\mathrm{irr}\\
            \text{T\_mQ1}=\frac{\dot{Q}_1}{\text{S\_Q1}}\\
            \text{T\_mQ2}=\frac{\dot{Q}_2}{\text{S\_Q2}}\\
            \text{T\_mQ3}=\frac{\dot{Q}_1 + \dot{Q}_2}{\text{S\_Q3}}
        """
        self.S_irr = 0
        for i in range(3):
            inl = self.inl[i]
            out = self.outl[i]
            p_star = inl.p.val_SI * (
                self.get_attr('pr' + str(i + 1)).val) ** 0.5
            s_i_star = s_mix_ph(
                [0, p_star, inl.h.val_SI, inl.fluid.val], T0=inl.T.val_SI)
            s_o_star = s_mix_ph(
                [0, p_star, out.h.val_SI, out.fluid.val], T0=out.T.val_SI)

            setattr(self, 'S_Q' + str(i + 1),
                    inl.m.val_SI * (s_o_star - s_i_star))
            S_Q = self.get_attr('S_Q' + str(i + 1))
            setattr(self, 'S_irr' + str(i + 1),
                    inl.m.val_SI * (out.s.val_SI - inl.s.val_SI) - S_Q)
            setattr(self, 'T_mQ' + str(i + 1),
                    inl.m.val_SI * (out.h.val_SI - inl.h.val_SI) / S_Q)

            self.S_irr += self.get_attr('S_irr' + str(i + 1))

        self.S_irr += self.S_Q1 + self.S_Q2 + self.S_Q3

    def get_plotting_data(self):
        """Generate a dictionary containing FluProDia plotting information.

        Returns
        -------
        data : dict
            A nested dictionary containing the keywords required by the
            :code:`calc_individual_isoline` method of the
            :code:`FluidPropertyDiagram` class. First level keys are the
            connection index ('in1' -> 'out1', therefore :code:`1` etc.).
        """
        return {
            i + 1: {
                'isoline_property': 'p',
                'isoline_value': self.inl[i].p.val,
                'isoline_value_end': self.outl[i].p.val,
                'starting_point_property': 'v',
                'starting_point_value': self.inl[i].vol.val,
                'ending_point_property': 'v',
                'ending_point_value': self.outl[i].vol.val
            } for i in range(3)}
