# -*- coding: utf-8

"""Module of class Turbomachine.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/turbomachinery/turbomachine.py

SPDX-License-Identifier: MIT
"""

import numpy as np

from tespy.components.component import Component
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import DataContainerSimple as dc_simple
from tespy.tools.global_vars import err


class Turbomachine(Component):
    r"""
    Parent class for compressor, pump and turbine.

    Equations

        **mandatory equations**

        - :py:meth:`tespy.components.component.Component.fluid_func`
        - :py:meth:`tespy.components.component.Component.mass_flow_func`

        **optional equations**

        .. math::

            0 = \dot{m}_{in} \cdot \left( h_{out} - h_{in} \right) - P\\
            0 = pr \cdot p_{in} - p_{out}

        isentropic efficiency equations (optional)

        - :py:meth:`tespy.components.turbomachinery.compressor.Compressor.eta_s_func`
        - :py:meth:`tespy.components.turbomachinery.pump.Pump.eta_s_func`
        - :py:meth:`tespy.components.turbomachinery.turbine.Turbine.eta_s_func`

        **additional equations**

        - :py:meth:`tespy.components.turbomachinery.compressor.Compressor.additional_equations`
        - :py:meth:`tespy.components.turbomachinery.pump.Pump.additional_equations`
        - :py:meth:`tespy.components.turbomachinery.turbine.Turbine.additional_equations`

    Inlets/Outlets

        - in1
        - out1

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

    P : float, tespy.tools.data_containers.ComponentProperties
        Power, :math:`P/\text{W}`

    pr : float, tespy.tools.data_containers.ComponentProperties
        Outlet to inlet pressure ratio, :math:`pr/1`

    Example
    -------
    For an example please refer to:

    - :class:`tespy.components.turbomachinery.compressor.Compressor`
    - :class:`tespy.components.turbomachinery.pump.Pump`
    - :class:`tespy.components.turbomachinery.turbine.Turbine`
    """

    @staticmethod
    def component():
        return 'turbomachine'

    @staticmethod
    def attr():
        return {'P': dc_cp(), 'pr': dc_cp(), 'Sirr': dc_simple()}

    @staticmethod
    def inlets():
        return ['in1']

    @staticmethod
    def outlets():
        return ['out1']

    def comp_init(self, nw):

        Component.comp_init(self, nw)

        # number of mandatroy equations for
        # fluid balance: num_fl
        # mass flow: 1
        self.num_eq = self.num_nw_fluids + 1
        for var in [self.P, self.pr]:
            if var.is_set:
                self.num_eq += 1

        self.jacobian = np.zeros((
            self.num_eq,
            self.num_i + self.num_o + self.num_vars,
            self.num_nw_vars))

        self.residual = np.zeros(self.num_eq)
        pos = self.num_nw_fluids
        self.jacobian[0:pos] = self.fluid_deriv()
        self.jacobian[pos:pos + 1] = self.mass_flow_deriv()

    def equations(self):
        r"""Calculate residual vector with results of equations."""
        k = 0
        ######################################################################
        # eqations for fluids
        if (any(np.absolute(self.residual[k:self.num_nw_fluids])) > err ** 2 or
                self.it % 4 == 0 or self.always_all_equations):
            self.residual[k:self.num_nw_fluids] = self.fluid_func()
        k += self.num_nw_fluids

        ######################################################################
        # eqations for mass flow balance
        self.residual[k] = self.mass_flow_func()
        k += 1

        ######################################################################
        # eqations for specified power
        if self.P.is_set:
            self.residual[k] = self.inl[0].m.val_SI * (
                self.outl[0].h.val_SI - self.inl[0].h.val_SI) - self.P.val

            k += 1

        ######################################################################
        # eqations for specified pressure ratio
        if self.pr.is_set:
            self.residual[k] = (
                self.pr.val * self.inl[0].p.val_SI - self.outl[0].p.val_SI)
            k += 1

        ######################################################################
        # additional equations
        self.additional_equations(k)

    def additional_equations(self, k):
        r"""Calculate results of additional equations."""
        return

    def derivatives(self, increment_filter):
        r"""Calculate partial derivatives for given equations."""
        ######################################################################
        # derivatives fluid and mass balance are static
        k = self.num_nw_fluids + 1

        ######################################################################
        # derivatives for specified power
        if self.P.is_set:
            self.jacobian[k, 0, 0] = (
                self.outl[0].h.val_SI - self.inl[0].h.val_SI)
            self.jacobian[k, 0, 2] = -self.inl[0].m.val_SI
            self.jacobian[k, 1, 2] = self.inl[0].m.val_SI
            k += 1

        ######################################################################
        # derivatives for specified pressure ratio
        if self.pr.is_set:
            self.jacobian[k, 0, 1] = self.pr.val
            self.jacobian[k, 1, 1] = -1
            k += 1

        ######################################################################
        # derivatives for additional equations
        self.additional_derivatives(increment_filter, k)

    def additional_derivatives(self, increment_filter, k):
        r"""Calculate partial derivatives for given additional equations."""
        return

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

                \dot{E} = \dot{m}_{in} \cdot \left(h_{out} - h_{in} \right)
        """
        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()
        val = i[0] * (o[2] - i[2])

        return val

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
        deriv = np.zeros((1, 2, self.num_nw_vars))
        f = self.calc_bus_value
        deriv[0, 0, 0] = self.numeric_deriv(f, 'm', 0, bus=bus)
        deriv[0, 0, 2] = self.numeric_deriv(f, 'h', 0, bus=bus)
        deriv[0, 1, 2] = self.numeric_deriv(f, 'h', 1, bus=bus)
        return deriv

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        self.P.val = self.inl[0].m.val_SI * (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI)
        self.pr.val = self.outl[0].p.val_SI / self.inl[0].p.val_SI

    def entropy_balance(self):
        r"""
        Calculate entropy balance of turbomachine.

        Note
        ----
        The entropy balance makes the follwing parameter available:

        .. math::

            \text{S\_irr}=\dot{m} \cdot \left(s_\mathrm{out}-s_\mathrm{in}
            \right)\\
        """
        self.S_irr = self.inl[0].m.val_SI * (
            self.outl[0].s.val_SI - self.inl[0].s.val_SI)

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
                'isoline_property': 's',
                'isoline_value': self.inl[0].s.val,
                'isoline_value_end': self.outl[0].s.val,
                'starting_point_property': 'v',
                'starting_point_value': self.inl[0].vol.val,
                'ending_point_property': 'v',
                'ending_point_value': self.outl[0].vol.val
            }
        }
