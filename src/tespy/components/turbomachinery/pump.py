# -*- coding: utf-8

"""Module of class Pump.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/turbomachinery/pump.py

SPDX-License-Identifier: MIT
"""

import logging

import numpy as np

from tespy.components.component import Component
from tespy.components.turbomachinery.turbomachine import Turbomachine
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import DataContainerSimple as dc_simple
from tespy.tools.fluid_properties import isentropic
from tespy.tools.fluid_properties import v_mix_ph
from tespy.tools.global_vars import err


class Pump(Turbomachine):
    r"""
    Class for axial or radial pumps.

    Equations

        **mandatory equations**

        - :py:meth:`tespy.components.component.Component.fluid_func`
        - :py:meth:`tespy.components.component.Component.mass_flow_func`

        **optional equations**

        .. math::

            0 = \dot{m}_{in} \cdot \left( h_{out} - h_{in} \right) - P\\
            0 = pr \cdot p_{in} - p_{out}

        **additional equations**

        - :py:meth:`tespy.components.turbomachinery.pump.Pump.additional_equations`

    Inlets/Outlets

        - in1
        - out1

    Image

        .. image:: _images/Pump.svg
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

    P : float, tespy.tools.data_containers.ComponentProperties
        Power, :math:`P/\text{W}`

    eta_s : float, tespy.tools.data_containers.ComponentProperties
        Isentropic efficiency, :math:`\eta_s/1`

    pr : float, tespy.tools.data_containers.ComponentProperties
        Outlet to inlet pressure ratio, :math:`pr/1`

    eta_s_char : tespy.tools.characteristics.CharLine, tespy.tools.data_containers.ComponentCharacteristics
        Characteristic curve for isentropic efficiency, provide CharLine as
        function :code:`func`.

    flow_char : tespy.tools.characteristics.CharLine, tespy.tools.data_containers.ComponentCharacteristics
        Characteristic curve for pressure rise vs. volumetric flow rate,
        provide CharLine as function :code:`func`.
        :math:`x/\frac{\text{m}^3}{\text{s}} \, y/\text{Pa}`.

    Example
    -------
    A pump with a known pump curve (difference pressure as function of
    volumetric flow) pumps 1,5 l/s of water in design conditions. E.g. for a
    given isentropic efficiency it is possible to calculate power consumption
    and pressure at the pump.

    >>> from tespy.components import Sink, Source, Pump
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> from tespy.tools.characteristics import CharLine
    >>> from tespy.tools import ComponentCharacteristics as dc_cc
    >>> import shutil
    >>> fluid_list = ['water']
    >>> nw = Network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ...     h_unit='kJ / kg', v_unit='l / s', iterinfo=False)
    >>> si = Sink('sink')
    >>> so = Source('source')
    >>> pu = Pump('pump')
    >>> pu.component()
    'pump'
    >>> inc = Connection(so, 'out1', pu, 'in1')
    >>> outg = Connection(pu, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)

    After that we calculate offdesign performance using
    the pump curve and a characteristic function for the pump efficiency. We
    can calulate the offdesign efficiency and the volumetric flow, if the
    difference pressure changed. The default characteristc lines are to be
    found in the :py:mod:`tespy.data` module. Of course, you are able to
    specify your own characteristcs.

    >>> v = np.array([0, 0.4, 0.8, 1.2, 1.6, 2]) / 1000
    >>> dp = np.array([15, 14, 12, 9, 5, 0]) * 1e5
    >>> char = CharLine(x=v, y=dp)
    >>> pu.set_attr(eta_s=0.8, flow_char=dc_cc(func=char, is_set=True),
    ... design=['eta_s'], offdesign=['eta_s_char'])
    >>> inc.set_attr(fluid={'water': 1}, p=1, T=20, v=1.5, design=['v'])
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(pu.pr.val, 0)
    7.0
    >>> round(outg.p.val - inc.p.val, 0)
    6.0
    >>> round(pu.P.val, 0)
    1125.0
    >>> outg.set_attr(p=12)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(pu.eta_s.val, 2)
    0.71
    >>> round(inc.v.val, 1)
    0.9
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    @staticmethod
    def component():
        return 'pump'

    @staticmethod
    def attr():
        return {
            'P': dc_cp(min_val=0),
            'eta_s': dc_cp(min_val=0, max_val=1), 'eta_s_char': dc_cc(),
            'pr': dc_cp(min_val=1),
            'flow_char': dc_cc(),
            'Sirr': dc_simple()
        }

    def comp_init(self, nw):

        Component.comp_init(self, nw)

        # number of mandatroy equations for
        # fluid balance: num_fl
        # mass flow: 1
        self.num_eq = self.num_nw_fluids + 1
        for var in [self.P, self.pr, self.eta_s, self.eta_s_char,
                    self.flow_char]:
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

    def additional_equations(self, k):
        r"""
        Calculate results of additional equations.

        Equations

            **optional equations**

            - :py:meth:`tespy.components.turbomachinery.pump.Pump.eta_s_func`
            - :py:meth:`tespy.components.turbomachinery.pump.Pump.eta_s_char_func`
            - :py:meth:`tespy.components.turbomachinery.pump.Pump.flow_char_func`
        """
        ######################################################################
        # eqations for specified isentropic efficiency
        if self.eta_s.is_set:
            if (np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0 or
                    self.always_all_equations):
                self.residual[k] = self.eta_s_func()
            k += 1

        ######################################################################
        # equations for specified isentropic efficiency characteristics
        if self.eta_s_char.is_set:
            if (np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0 or
                    self.always_all_equations):
                self.residual[k] = self.eta_s_char_func()
            k += 1

        ######################################################################
        # equations for specified pressure rise vs. flowrate characteristics
        if self.flow_char.is_set:
            if (np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0 or
                    self.always_all_equations):
                self.residual[k] = self.flow_char_func()
            k += 1

    def additional_derivatives(self, increment_filter, k):
        r"""Calculate partial derivatives for given additional equations."""
        ######################################################################
        # derivatives for specified isentropic efficiency
        if self.eta_s.is_set:
            f = self.eta_s_func
            if not increment_filter[0, 1]:
                self.jacobian[k, 0, 1] = self.numeric_deriv(f, 'p', 0)
            if not increment_filter[1, 1]:
                self.jacobian[k, 1, 1] = self.numeric_deriv(f, 'p', 1)
            if not increment_filter[0, 2]:
                self.jacobian[k, 0, 2] = self.numeric_deriv(f, 'h', 0)
            self.jacobian[k, 1, 2] = -self.eta_s.val
            k += 1

        ######################################################################
        # derivatives for specified isentropic efficiency characteristics
        if self.eta_s_char.is_set:
            f = self.eta_s_char_func
            if not increment_filter[0, 0]:
                self.jacobian[k, 0, 0] = self.numeric_deriv(f, 'm', 0)
            if not increment_filter[0, 1]:
                self.jacobian[k, 0, 1] = self.numeric_deriv(f, 'p', 0)
            if not increment_filter[0, 2]:
                self.jacobian[k, 0, 2] = self.numeric_deriv(f, 'h', 0)
            if not increment_filter[1, 1]:
                self.jacobian[k, 1, 1] = self.numeric_deriv(f, 'p', 1)
            if not increment_filter[1, 2]:
                self.jacobian[k, 1, 2] = self.numeric_deriv(f, 'h', 1)
            k += 1

        ######################################################################
        # derivatives for specified pressure rise vs. flowrate characteristics
        if self.flow_char.is_set:
            f = self.flow_char_func
            if not increment_filter[0, 0]:
                self.jacobian[k, 0, 0] = self.numeric_deriv(f, 'm', 0)
            if not increment_filter[0, 2]:
                self.jacobian[k, 0, 2] = self.numeric_deriv(f, 'h', 0)
            for i in range(2):
                if not increment_filter[i, 1]:
                    self.jacobian[k, i, 1] = self.numeric_deriv(f, 'p', i)
            k += 1

    def eta_s_func(self):
        r"""
        Equation for given isentropic efficiency.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                0 = -\left( h_{out} - h_{in} \right) \cdot \eta_{s,c} +
                \left( h_{out,s} - h_{in} \right)
        """
        return (
            -(self.outl[0].h.val_SI - self.inl[0].h.val_SI) * self.eta_s.val +
            (isentropic(
                self.inl[0].to_flow(), self.outl[0].to_flow(),
                T0=self.inl[0].T.val_SI) - self.inl[0].h.val_SI))

    def eta_s_char_func(self):
        r"""
        Equation for given isentropic efficiency characteristic.

        Returns
        -------
        res : ndarray
            Residual value of equation.

            .. math::

                0 = \left( h_{out} - h_{in} \right) \cdot
                \frac{\Delta h_{s,ref}}{\Delta h_{ref}}
                \cdot char\left( \frac{\dot{m}_{in} \cdot
                v_{in}}{\dot{m}_{in,ref} \cdot v_{in,ref}} \right) -
                \left( h_{out,s} - h_{in} \right)
        """
        i = self.inl[0]
        o = self.outl[0]
        v_i = v_mix_ph(i.to_flow(), T0=self.inl[0].T.val_SI)
        expr = i.m.val_SI * v_i / (i.v.design)

        return (
            (o.h.val_SI - i.h.val_SI) * self.eta_s.design *
            self.eta_s_char.func.evaluate(expr) - (
                isentropic(
                    i.to_flow(), o.to_flow(), T0=self.inl[0].T.val_SI) -
                i.h.val_SI))

    def flow_char_func(self):
        r"""
        Equation for given flow characteristic of a pump.

        Returns
        -------
        res : ndarray
            Residual value of equation.

            .. math::

                0 = p_{out} - p_{in} - char\left( \dot{m}_{in}
                \cdot v_{in} \right)
        """
        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()

        expr = i[0] * v_mix_ph(i, T0=self.inl[0].T.val_SI)

        return o[1] - i[1] - self.flow_char.func.evaluate(expr)

    def convergence_check(self):
        r"""
        Perform a convergence check.

        Note
        ----
        Manipulate enthalpies/pressure at inlet and outlet if not specified by
        user to match physically feasible constraints.
        """
        i, o = self.inl, self.outl

        if not o[0].p.val_set and o[0].p.val_SI < i[0].p.val_SI:
            o[0].p.val_SI = o[0].p.val_SI * 2
        if not i[0].p.val_set and o[0].p.val_SI < i[0].p.val_SI:
            i[0].p.val_SI = o[0].p.val_SI * 0.5

        if not o[0].h.val_set and o[0].h.val_SI < i[0].h.val_SI:
            o[0].h.val_SI = o[0].h.val_SI * 1.1
        if not i[0].h.val_set and o[0].h.val_SI < i[0].h.val_SI:
            i[0].h.val_SI = o[0].h.val_SI * 0.9

        if self.flow_char.is_set:
            expr = i[0].m.val_SI * v_mix_ph(i[0].to_flow(), T0=i[0].T.val_SI)

            if expr > self.flow_char.func.x[-1] and not i[0].m.val_set:
                i[0].m.val_SI = (self.flow_char.func.x[-1] /
                                 v_mix_ph(i[0].to_flow(), T0=i[0].T.val_SI))
            elif expr < self.flow_char.func.x[1] and not i[0].m.val_set:
                i[0].m.val_SI = (self.flow_char.func.x[0] /
                                 v_mix_ph(i[0].to_flow(), T0=i[0].T.val_SI))
            else:
                pass

    @staticmethod
    def initialise_Source(c, key):
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
                3 \cdot 10^5 & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 10e5
        elif key == 'h':
            return 3e5

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
                10^5 & \text{key = 'p'}\\
                2.9 \cdot 10^5 & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 1e5
        elif key == 'h':
            return 2.9e5

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        Turbomachine.calc_parameters(self)

        self.eta_s.val = (
            (isentropic(
                self.inl[0].to_flow(), self.outl[0].to_flow(),
                T0=self.inl[0].T.val_SI) -
             self.inl[0].h.val_SI) /
            (self.outl[0].h.val_SI - self.inl[0].h.val_SI))

        if self.eta_s_char.is_set:
            # get bound errors for isentropic efficiency characteristics
            expr = self.inl[0].v.val_SI / self.inl[0].v.design
            self.eta_s_char.func.get_bound_errors(expr, self.label)

        if self.flow_char.is_set:
            # get bound errors for flow characteristics
            expr = self.inl[0].v.val_SI
            self.flow_char.func.get_bound_errors(expr, self.label)

        self.check_parameter_bounds()

    def exergy_balance(self, T0):
        r"""
        Calculate exergy balance of a pump.

        Parameters
        ----------
        T0 : float
            Ambient temperature T0 / K.

        .. math::

            \dot{E}_\mathrm{P} =
            \begin{cases}
            \dot{E}_\mathrm{out}^\mathrm{PH} - \dot{E}_\mathrm{in}^\mathrm{PH}
            & T_\mathrm{in}, T_\mathrm{out} \geq T_0\\
            \dot{E}_\mathrm{out}^\mathrm{T} + \dot{E}_\mathrm{out}^\mathrm{M} -
            \dot{E}_\mathrm{in}^\mathrm{M}
            & T_\mathrm{out} > T_0 \leq T_\mathrm{in}\\
            \dot{E}_\mathrm{out}^\mathrm{M} - \dot{E}_\mathrm{in}^\mathrm{M}
            & T_0 \geq T_\mathrm{in}, T_\mathrm{out}\\
            \end{cases}

            \dot{E}_\mathrm{F} =
            \begin{cases}
            P & T_\mathrm{in}, T_\mathrm{out} \geq T_0\\
            P + \dot{E}_\mathrm{in}^\mathrm{T}
            & T_\mathrm{out} > T_0 \leq T_\mathrm{in}\\
            P + \dot{E}_\mathrm{in}^\mathrm{T} -\dot{E}_\mathrm{out}^\mathrm{T}
            & T_0 \geq T_\mathrm{in}, T_\mathrm{out}\\
            \end{cases}

            \dot{E}_\mathrm{bus} = P
        """
        if self.inl[0].T.val_SI >= T0 and self.outl[0].T.val_SI >= T0:
            self.E_P = self.outl[0].Ex_physical - self.inl[0].Ex_physical
            self.E_F = self.P.val
        elif self.inl[0].T.val_SI <= T0 and self.outl[0].T.val_SI > T0:
            self.E_P = self.outl[0].Ex_therm + (
                self.outl[0].Ex_mech - self.inl[0].Ex_mech)
            self.E_F = self.P.val + self.inl[0].Ex_therm
        elif self.inl[0].T.val_SI <= T0 and self.outl[0].T.val_SI <= T0:
            self.E_P = self.outl[0].Ex_mech - self.inl[0].Ex_mech
            self.E_F = self.P.val + (
                self.inl[0].Ex_therm - self.outl[0].Ex_therm)
        else:
            msg = ('Exergy balance of a pump, where outlet temperature is '
                   'smaller than inlet temperature is not implmented.')
            logging.warning(msg)
            self.E_P = np.nan
            self.E_F = np.nan

        self.E_bus = self.P.val
        self.E_D = self.E_F - self.E_P
        self.epsilon = self.E_P / self.E_F
