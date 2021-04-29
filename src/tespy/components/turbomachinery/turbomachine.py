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
from tespy.tools.document_models import generate_latex_eq


class Turbomachine(Component):
    r"""
    Parent class for compressor, pump and turbine.

    **Mandatory Equations**

    - :py:meth:`tespy.components.component.Component.fluid_func`
    - :py:meth:`tespy.components.component.Component.mass_flow_func`

    **Optional Equations**

    - :py:meth:`tespy.components.component.Component.pr_func`
    - :py:meth:`tespy.components.turbomachinery.turbomachine.Turbomachine.energy_balance_func`

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

    P : float, dict
        Power, :math:`P/\text{W}`

    pr : float, dict, :code:`"var"`
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

    def get_variables(self):
        return {
            'P': dc_cp(
                deriv=self.energy_balance_deriv, num_eq=1,
                func=self.energy_balance_func,
                latex=self.energy_balance_func_doc),
            'pr': dc_cp(
                deriv=self.pr_deriv, num_eq=1,
                func=self.pr_func, func_params={'pr': 'pr'},
                latex=self.pr_func_doc)
        }

    @staticmethod
    def inlets():
        return ['in1']

    @staticmethod
    def outlets():
        return ['out1']

    def energy_balance_func(self):
        r"""
        Calculate energy balance of a turbomachine.

        Returns
        -------
        residual : float
            Residual value of turbomachine energy balance

            .. math::

                0=\dot{m}_{in}\cdot\left(h_{out}-h_{in}\right)-P
        """
        return self.inl[0].m.val_SI * (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI) - self.P.val

    def energy_balance_func_doc(self, label):
        r"""
        Calculate energy balance of a turbomachine.

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
            r'0=\dot{m}_\mathrm{in}\cdot\left(h_\mathrm{out}-h_\mathrm{in}'
            r'\right)-P')
        return generate_latex_eq(self, latex, label)

    def energy_balance_deriv(self, increment_filter, k):
        r"""
        Calculate partial derivatives of energy balance of a turbomachine.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        self.jacobian[k, 0, 0] = self.outl[0].h.val_SI - self.inl[0].h.val_SI
        self.jacobian[k, 0, 2] = -self.inl[0].m.val_SI
        self.jacobian[k, 1, 2] = self.inl[0].m.val_SI

    def bus_func(self, bus):
        r"""
        Calculate the value of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus.Bus
            TESPy bus object.

        Returns
        -------
        residual : float
            Value of energy transfer :math:`\dot{E}`. This value is passed to
            :py:meth:`tespy.components.component.Component.calc_bus_value`
            for value manipulation according to the specified characteristic
            line of the bus.

            .. math::

                \dot{E} = \dot{m}_{in} \cdot \left(h_{out} - h_{in} \right)
        """
        return self.inl[0].m.val_SI * (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI)

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
