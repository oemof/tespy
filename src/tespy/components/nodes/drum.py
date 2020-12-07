# -*- coding: utf-8

"""Module of class Drum.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/nodes/drum.py

SPDX-License-Identifier: MIT
"""

import numpy as np

from tespy.components.component import Component
from tespy.tools.fluid_properties import dh_mix_dpQ
from tespy.tools.fluid_properties import h_mix_pQ


class Drum(Component):
    r"""
    A drum separates saturated gas from saturated liquid.

    Equations

        **mandatory equations**

        - :py:meth:`tespy.components.nodes.drum.Drum.fluid_func`
        - :py:meth:`tespy.components.component.Component.mass_flow_func`

        .. math::

            0 = \sum_i \left(\dot{m}_{i,in} \cdot h_{i,in} \right) -
            \sum_j \left(\dot{m}_{j,out} \cdot h_{j,out} \right)\\
            \forall i \in inlets, \; \forall j \in outlet

            0 = p_{in,1} - p_{out,i}\\
            \forall i \in \mathrm{outlets}

            0 = h_{1,out} - h\left(p, x=0 \right)

            0 = h_{2,out} - h\left(p, x=1 \right)\\
            x: \text{vapour mass fraction}

    Inlets/Outlets

        - in1, in2 (index 1: from economiser, index 2: from evaporator)
        - out1, out2 (index 1: to evaporator, index 2: to superheater)

    Image

        .. image:: _images/Drum.svg
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

    Note
    ----
    If you are using a drum in a network with multiple fluids, it is likely
    the fluid propagation causes trouble. If this is the case, try to
    specify the fluid composition at another connection of your network.

    This component assumes, that the fluid composition between outlet 1 and
    inlet 2 does not change, thus there is no equation for the fluid mass
    fraction at the inlet 2!

    Example
    -------
    The drum separates saturated gas from saturated liquid. The liquid phase is
    transported to an evaporator, the staturated gas phase is extracted from
    the drum. In this example ammonia is evaporated using ambient air. A
    characteristic function is applied for the heat transfer coefficient of the
    evaporator. It is possible to load the CharLine with the function
    :code:`load_default_char` from the default lines. We want to use the
    'EVAPORATING FLUID' lines of the heat exchanger.

    >>> from tespy.components import Sink, Source, Drum, Pump, HeatExchanger
    >>> from tespy.connections import Connection, Ref
    >>> from tespy.networks import Network
    >>> from tespy.tools.characteristics import CharLine
    >>> from tespy.tools.characteristics import load_default_char as ldc
    >>> import shutil
    >>> import numpy as np
    >>> nw = Network(fluids=['NH3', 'air'], T_unit='C', p_unit='bar',
    ...     h_unit='kJ / kg', iterinfo=False)
    >>> fa = Source('feed ammonia')
    >>> amb_in = Source('air inlet')
    >>> amb_out = Sink('air outlet')
    >>> s = Sink('steam')
    >>> dr = Drum('drum')
    >>> dr.component()
    'drum'
    >>> ev = HeatExchanger('evaporator')
    >>> erp = Pump('evaporator reciculation pump')
    >>> f_dr = Connection(fa, 'out1', dr, 'in1')
    >>> dr_erp = Connection(dr, 'out1', erp, 'in1')
    >>> erp_ev = Connection(erp, 'out1', ev, 'in2')
    >>> ev_dr = Connection(ev, 'out2', dr, 'in2')
    >>> dr_s = Connection(dr, 'out2', s, 'in1')
    >>> nw.add_conns(f_dr, dr_erp, erp_ev, ev_dr, dr_s)
    >>> amb_ev = Connection(amb_in, 'out1', ev, 'in1')
    >>> ev_amb = Connection(ev, 'out1', amb_out, 'in1')
    >>> nw.add_conns(amb_ev, ev_amb)

    The ambient air enters the evaporator at 30 °C. The pinch point temperature
    difference (ttd_l) of the evaporator is at 5 K, and 1 MW of heat should be
    transferred. State of ammonia at the inlet is at -5 °C and 5 bar. From this
    design it is possible to calculate offdesign performance at 75 % part load.

    >>> char1 = ldc('heat exchanger', 'kA_char1', 'DEFAULT',
    ... CharLine)
    >>> char2 = ldc('heat exchanger', 'kA_char2', 'EVAPORATING FLUID',
    ... CharLine)
    >>> ev.set_attr(pr1=0.999, pr2=0.99, ttd_l=5, kA_char1=char1,
    ...     kA_char2=char2, design=['pr1', 'ttd_l'],
    ...     offdesign=['zeta1', 'kA_char'])
    >>> ev.set_attr(Q=-1e6)
    >>> erp.set_attr(eta_s=0.8)
    >>> f_dr.set_attr(p=5, T=-5)
    >>> erp_ev.set_attr(m=Ref(f_dr, 4, 0), fluid={'air': 0, 'NH3': 1})
    >>> amb_ev.set_attr(fluid={'air': 1, 'NH3': 0}, T=30)
    >>> ev_amb.set_attr(p=1)
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(ev_amb.T.val - erp_ev.T.val ,1)
    5.0
    >>> round(f_dr.h.val, 1)
    322.7
    >>> round(dr_erp.h.val, 1)
    364.9
    >>> round(ev_dr.h.val, 1)
    687.2
    >>> round(f_dr.m.val, 2)
    0.78
    >>> ev.set_attr(Q=-0.75e6)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(f_dr.m.val, 2)
    0.58
    >>> round(ev_amb.T.val - erp_ev.T.val ,1)
    3.0
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    @staticmethod
    def component():
        return 'drum'

    @staticmethod
    def inlets():
        return ['in1', 'in2']

    @staticmethod
    def outlets():
        return ['out1', 'out2']

    def comp_init(self, nw):

        Component.comp_init(self, nw)

        # number of mandatroy equations for
        # fluid balance: num_fl * 2
        # mass flow: 1
        # pressure: 3
        # enthalpy: 1
        # saturated liquid outlet: 1
        # saturated gas outlet: 1
        self.num_eq = self.num_nw_fluids * 2 + 7

        self.jacobian = np.zeros((
            self.num_eq,
            self.num_i + self.num_o + self.num_vars,
            self.num_nw_vars))

        self.residual = np.zeros(self.num_eq)
        pos = self.num_nw_fluids * 2
        self.jacobian[0:pos] = self.fluid_deriv()
        self.jacobian[pos:pos + 1] = self.mass_flow_deriv()
        self.jacobian[pos + 1:pos + 4] = self.pressure_deriv()

    def equations(self):
        r"""Calculate residual vector with results of equations."""
        k = 0
        ######################################################################
        # eqations for fluid balance
        self.residual[k:k + self.num_nw_fluids * 2] = self.fluid_func()
        k += self.num_nw_fluids * 2

        ######################################################################
        # eqations for mass flow balance
        self.residual[k] = self.mass_flow_func()
        k += 1

        ######################################################################
        # eqations for pressure
        p = self.inl[0].p.val_SI
        for c in [self.inl[1]] + self.outl:
            self.residual[k] = p - c.p.val_SI
            k += 1

        ######################################################################
        # eqations for enthalpy
        val = 0
        for i in self.inl:
            val += i.m.val_SI * i.h.val_SI
        for o in self.outl:
            val -= o.m.val_SI * o.h.val_SI
        self.residual[k] = val
        k += 1

        ######################################################################
        # eqations for staturated fluid state at outlets
        self.residual[k] = h_mix_pQ(
            self.outl[0].to_flow(), 0) - self.outl[0].h.val_SI
        k += 1
        self.residual[k] = h_mix_pQ(
            self.outl[1].to_flow(), 1) - self.outl[1].h.val_SI
        k += 1

    def derivatives(self, increment_filter):
        r"""Calculate partial derivatives for given equations."""
        ######################################################################
        # derivatives fluid, mass flow and pressure balance are static
        k = self.num_nw_fluids * 2 + 4

        ######################################################################
        # derivatives for energy balance equation
        i = 0
        for inl in self.inl:
            self.jacobian[k, i, 0] = inl.h.val_SI
            self.jacobian[k, i, 2] = inl.m.val_SI
            i += 1
        j = 0
        for outl in self.outl:
            self.jacobian[k, j + i, 0] = -outl.h.val_SI
            self.jacobian[k, j + i, 2] = -outl.m.val_SI
            j += 1
        k += 1

        ######################################################################
        # derivatives of equations for saturated states at outlets
        self.jacobian[k, 2, 1] = dh_mix_dpQ(self.outl[0].to_flow(), 0)
        self.jacobian[k, 2, 2] = -1
        k += 1
        self.jacobian[k, 3, 1] = dh_mix_dpQ(self.outl[1].to_flow(), 1)
        self.jacobian[k, 3, 2] = -1
        k += 1

    def fluid_func(self):
        r"""
        Calculate the vector of residual values for fluid balance equations.

        Returns
        -------
        residual : list
            Vector of residual values for component's fluid balance.

            .. math::

                0 = fluid_{i,in_1} - fluid_{i,out_{j}}\\
                \forall i \in \mathrm{fluid}, \; \forall j \in inlets

        """
        residual = []

        for o in self.outl:
            for fluid, x in self.inl[0].fluid.val.items():
                residual += [x - o.fluid.val[fluid]]
        return residual

    def fluid_deriv(self):
        r"""
        Calculate partial derivatives for all fluid balance equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((2 * self.num_nw_fluids, 4, self.num_nw_vars))
        for k in range(2):
            for i in range(self.num_nw_fluids):
                deriv[i + k * self.num_nw_fluids, 0, i + 3] = 1
                deriv[i + k * self.num_nw_fluids, k + 2, i + 3] = -1
        return deriv

    def pressure_deriv(self):
        r"""
        Calculate partial derivatives for pressure equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the pressure equations.
        """
        deriv = np.zeros((3, 4, self.num_nw_vars))
        for k in range(3):
            deriv[k, 0, 1] = 1
            deriv[k, k + 1, 1] = -1
        return deriv

    @staticmethod
    def initialise_source(c, key):
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
                10^6 & \text{key = 'p'}\\
                h\left(p, x=0 \right) & \text{key = 'h' at outlet 1}\\
                h\left(p, x=1 \right) & \text{key = 'h' at outlet 2}
                \end{cases}
        """
        if key == 'p':
            return 10e5
        elif key == 'h':
            if c.source_id == 'out1':
                return h_mix_pQ(c.to_flow(), 0)
            else:
                return h_mix_pQ(c.to_flow(), 1)

    @staticmethod
    def initialise_target(c, key):
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
                10^6 & \text{key = 'p'}\\
                h\left(p, x=0 \right) & \text{key = 'h' at inlet 1}\\
                h\left(p, x=0.7 \right) & \text{key = 'h' at inlet 2}
                \end{cases}
        """
        if key == 'p':
            return 10e5
        elif key == 'h':
            if c.target_id == 'in1':
                return h_mix_pQ(c.to_flow(), 0)
            else:
                return h_mix_pQ(c.to_flow(), 0.7)

    def propagate_fluid_to_target(self, inconn, start):
        r"""
        Propagate the fluids towards connection's target in recursion.

        Parameters
        ----------
        inconn : tespy.connections.connection.Connection
            Connection to initialise.

        start : tespy.components.component.Component
            This component is the fluid propagation starting point.
            The starting component is saved to prevent infinite looping.
        """
        if self != start:
            start = self
            for outconn in self.outl:
                for fluid, x in inconn.fluid.val.items():
                    if (outconn.fluid.val_set[fluid] is False and
                            outconn.good_starting_values is False):
                        outconn.fluid.val[fluid] = x

                outconn.target.propagate_fluid_to_target(outconn, start)

    def propagate_fluid_to_source(self, outconn, start):
        r"""
        Propagate the fluids towards connection's source in recursion.

        Parameters
        ----------
        outconn : tespy.connections.connection.Connection
            Connection to initialise.

        start : tespy.components.component.Component
            This component is the fluid propagation starting point.
            The starting component is saved to prevent infinite looping.
        """
        if self != start:
            start = self
            for inconn in self.inl:
                for fluid, x in outconn.fluid.val.items():
                    if (inconn.fluid.val_set[fluid] is False and
                            inconn.good_starting_values is False):
                        inconn.fluid.val[fluid] = x

                inconn.source.propagate_fluid_to_source(inconn, start)

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
            1: {
                'isoline_property': 'p',
                'isoline_value': self.outl[0].p.val,
                'isoline_value_end': self.outl[1].p.val,
                'starting_point_property': 'v',
                'starting_point_value': self.outl[0].vol.val,
                'ending_point_property': 'v',
                'ending_point_value': self.outl[1].vol.val
            }}
