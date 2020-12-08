# -*- coding: utf-8

"""Module of class Compressor.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/turbomachinery/compressor.py

SPDX-License-Identifier: MIT
"""

import logging

import numpy as np

from tespy.components.component import Component
from tespy.components.turbomachinery.turbomachine import Turbomachine
from tespy.tools.characteristics import CompressorMap
from tespy.tools.characteristics import load_default_char as ldc
from tespy.tools.data_containers import ComponentCharacteristicMaps as dc_cm
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import DataContainerSimple as dc_simple
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.fluid_properties import isentropic
from tespy.tools.global_vars import err


class Compressor(Turbomachine):
    r"""
    Class for axial or radial compressor.

    Equations

        **mandatory equations**

        - :py:meth:`tespy.components.component.Component.fluid_func`
        - :py:meth:`tespy.components.component.Component.mass_flow_func`

        **optional equations**

        .. math::

            0 = \dot{m}_{in} \cdot \left( h_{out} - h_{in} \right) - P\\
            0 = pr \cdot p_{in} - p_{out}

        **additional equations**

        - :py:meth:`tespy.components.turbomachinery.compressor.Compressor.additional_equations`

    Inlets/Outlets

        - in1
        - out1

    Image

        .. image:: _images/Compressor.svg
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

    char_map : tespy.tools.characteristics.CompressorMap, tespy.tools.data_containers.ComponentCharacteristicMaps
        Characteristic map for pressure rise and isentropic efficiency vs.
        nondimensional mass flow, see
        :py:class:`tespy.tools.characteristics.CompressorMap` for further
        information. Provide a CompressorMap as function :code:`func`.

    igva : str, float, tespy.tools.data_containers.ComponentProperties
        Inlet guide vane angle, :math:`igva/^\circ`.

    Example
    -------
    Create an air compressor model and calculate the power required for
    compression of 50 l/s of ambient air to 5 bars. Using a generic compressor
    map how does the efficiency change in different operation mode (e.g. 90 %
    of nominal volumetric flow)?

    >>> from tespy.components import Sink, Source, Compressor
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> import shutil
    >>> fluid_list = ['air']
    >>> nw = Network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ... h_unit='kJ / kg', v_unit='l / s', iterinfo=False)
    >>> si = Sink('sink')
    >>> so = Source('source')
    >>> comp = Compressor('compressor')
    >>> comp.component()
    'compressor'
    >>> inc = Connection(so, 'out1', comp, 'in1')
    >>> outg = Connection(comp, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)

    Specify the compressor parameters: nominal efficiency and pressure ratio.
    For offdesign mode the characteristic map is selected instead of the
    isentropic efficiency. For offdesign, the inlet guide vane angle should be
    variable in order to maintain the same pressure ratio at a different
    volumetric flow.

    >>> comp.set_attr(pr=5, eta_s=0.8, design=['eta_s'],
    ... offdesign=['char_map'])
    >>> inc.set_attr(fluid={'air': 1}, p=1, T=20, v=50)
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(comp.P.val, 0)
    12772.0
    >>> inc.set_attr(v=45)
    >>> comp.set_attr(igva='var')
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(comp.eta_s.val, 2)
    0.77
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    @staticmethod
    def component():
        return 'compressor'

    @staticmethod
    def attr():
        return {
            'P': dc_cp(min_val=0),
            'eta_s': dc_cp(min_val=0, max_val=1),
            'eta_s_char': dc_cc(param='m'),
            'pr': dc_cp(min_val=1),
            'igva': dc_cp(min_val=-90, max_val=90, d=1e-3, val=0),
            'char_map': dc_cm(),
            'Sirr': dc_simple()
        }

    def comp_init(self, nw):

        Component.comp_init(self, nw)

        if self.char_map.func is None:
            self.char_map.func = ldc(
                self.component(), 'char_map', 'DEFAULT', CompressorMap)

            if self.char_warnings:
                msg = ('Created characteristic map for parameter char_map '
                       'at component ' + self.label + ' from default data.\n'
                       'You can specify your own data using component.char_map'
                       '.set_attr(func=custom_char).\n'
                       'If you want to disable these warnings use '
                       'component.char_warnings=False.')
                logging.warning(msg)

        # number of mandatroy equations for
        # fluid balance: num_fl
        # mass flow: 1
        self.num_eq = self.num_nw_fluids + 1
        # characteristic map delivers two equations
        for var in [self.P, self.pr, self.eta_s, self.char_map, self.char_map,
                    self.eta_s_char]:
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

            - :py:meth:`tespy.components.turbomachinery.compressor.Compressor.eta_s_func`
            - :py:meth:`tespy.components.turbomachinery.compressor.Compressor.eta_s_char_func`
            - :py:meth:`tespy.components.turbomachinery.compressor.Compressor.char_map_func`
        """
        ######################################################################
        # eqations for specified isentropic efficiency
        if self.eta_s.is_set:
            if (np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0 or
                    self.always_all_equations):
                self.residual[k] = self.eta_s_func()
            k += 1

        ######################################################################
        # equation for specified isentropic efficiency characteristics
        if self.eta_s_char.is_set:
            if (np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0 or
                    self.always_all_equations):
                self.residual[k] = self.eta_s_char_func()
            k += 1

        ######################################################################
        # equations for specified characteristic map
        if self.char_map.is_set:
            if (any(np.absolute(self.residual[k:k + 2])) > err ** 2 or
                    self.it % 4 == 0 or self.always_all_equations):
                self.residual[k:k + 2] = self.char_map_func()
            k += 2

    def additional_derivatives(self, increment_filter, k):
        r"""Calculate partial derivatives for given additional equations."""
        ######################################################################
        # derivatives for specified isentropic efficiency
        if self.eta_s.is_set:
            f = self.eta_s_func
            if not increment_filter[0, 1]:
                self.jacobian[k, 0, 1] = self.numeric_deriv(f, 'p', 0)
            if not increment_filter[0, 1]:
                self.jacobian[k, 1, 1] = self.numeric_deriv(f, 'p', 1)
            if not increment_filter[0, 1]:
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
            if not increment_filter[1, 1]:
                self.jacobian[k, 1, 1] = self.numeric_deriv(f, 'p', 1)
            if not increment_filter[0, 2]:
                self.jacobian[k, 0, 2] = self.numeric_deriv(f, 'h', 0)
            if not increment_filter[1, 2]:
                self.jacobian[k, 1, 2] = self.numeric_deriv(f, 'h', 1)
            k += 1

        ######################################################################
        # derivatives for specified characteristic map
        if self.char_map.is_set:
            f = self.char_map_func
            if not increment_filter[0, 0]:
                m11 = self.numeric_deriv(f, 'm', 0)
                self.jacobian[k, 0, 0] = m11[0]
                self.jacobian[k + 1, 0, 0] = m11[1]
            if not increment_filter[0, 1]:
                p11 = self.numeric_deriv(f, 'p', 0)
                self.jacobian[k, 0, 1] = p11[0]
                self.jacobian[k + 1, 0, 1] = p11[1]
            if not increment_filter[0, 2]:
                h11 = self.numeric_deriv(f, 'h', 0)
                self.jacobian[k, 0, 2] = h11[0]
                self.jacobian[k + 1, 0, 2] = h11[1]

            if not increment_filter[1, 1]:
                p21 = self.numeric_deriv(f, 'p', 1)
                self.jacobian[k, 1, 1] = p21[0]
                self.jacobian[k + 1, 1, 1] = p21[1]
            if not increment_filter[1, 2]:
                h21 = self.numeric_deriv(f, 'h', 1)
                self.jacobian[k, 1, 2] = h21[0]
                self.jacobian[k + 1, 1, 2] = h21[1]

            if self.igva.is_var:
                igva = self.numeric_deriv(f, 'igva', 1)
                if self.igva.is_var:
                    self.jacobian[k, 2 + self.igva.var_pos, 0] = igva[0]
                    self.jacobian[k + 1, 2 + self.igva.var_pos, 0] = igva[1]
            k += 2

    def eta_s_func(self):
        r"""
        Equation for given isentropic efficiency of a compressor.

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
                \cdot char\left( \dot{m}_{in} \cdot v_{in} \right) -
                \left( h_{out,s} - h_{in} \right)
        """
        # actual values
        i = self.inl[0]
        o = self.outl[0]
        expr = 1
        if self.eta_s_char.param == 'm':
            if not np.isnan(i.m.design):
                expr = i.m.val_SI / i.m.design
        elif self.eta_s_char.param == 'pr':
            if not np.isnan(self.pr.design):
                expr = (o.p.val_SI / i.p.val_SI) / self.pr.design
        else:
            msg = ('Must provide a parameter for eta_s_char at component ' +
                   self.label + '.')
            logging.error(msg)
            raise ValueError(msg)

        return (
            self.eta_s.design * self.eta_s_char.func.evaluate(expr) *
            (o.h.val_SI - i.h.val_SI) - (isentropic(
                i.to_flow(), o.to_flow(), T0=self.inl[0].T.val_SI) -
             i.h.val_SI))

    def char_map_func(self):
        r"""
        Equation for characteristic map of compressor.

        Returns
        -------
        res : ndarray (Z1, Z2)
            Residual values of equations.

        Note
        ----
        - X: speedline index (rotational speed is constant)
        - Y: nondimensional mass flow
        - Z1: pressure ratio equation
        - Z2: isentropic efficiency equation
        - igva: variable inlet guide vane angle (assumed 0° if not
          specified)

        .. math::

            X = \sqrt{\frac{T_\mathrm{1,ref}}{T_\mathrm{1}}}

            Y = \frac{\dot{m}_\mathrm{1} \cdot p_\mathrm{1,ref}}
            {\dot{m}_\mathrm{1,ref} \cdot p_\mathrm{1} \cdot X}

            Z1 = \frac{p_2 \cdot p_\mathrm{1,ref}}
            {p_1 \cdot p_\mathrm{2,ref}}-
            pr_{c}(char(m, igva))

            Z2 = \frac{\eta_\mathrm{s,c}}{\eta_\mathrm{s,c,ref}} -
            \eta_{s,c}(char(m, igva))
        """
        i = self.inl[0]
        o = self.outl[0]
        T_i = T_mix_ph(i.to_flow(), T0=self.inl[0].T.val_SI)

        x = np.sqrt(i.T.design / T_i)
        y = (i.m.val_SI * i.p.design) / (i.m.design * i.p.val_SI * x)

        pr, eta = self.char_map.func.evaluate(x, y, igva=self.igva.val)

        z1 = (o.p.val_SI / i.p.val_SI) / self.pr.design - pr
        z2 = ((
            isentropic(i.to_flow(), o.to_flow(), T0=self.inl[0].T.val_SI) -
            i.h.val_SI) / (o.h.val_SI - i.h.val_SI) / self.eta_s.design - eta)

        return np.array([z1, z2])

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
            o[0].p.val_SI = o[0].p.val_SI * 1.1

        if not o[0].h.val_set and o[0].h.val_SI < i[0].h.val_SI:
            o[0].h.val_SI = o[0].h.val_SI * 1.1

        if not i[0].p.val_set and o[0].p.val_SI < i[0].p.val_SI:
            i[0].p.val_SI = o[0].p.val_SI * 0.9
        if not i[0].h.val_set and o[0].h.val_SI < i[0].h.val_SI:
            i[0].h.val_SI = o[0].h.val_SI * 0.9

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
                6 \cdot 10^5 & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 10e5
        elif key == 'h':
            return 6e5

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
                4 \cdot 10^5 & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 1e5
        elif key == 'h':
            return 4e5

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        Turbomachine.calc_parameters(self)

        self.eta_s.val = (
            (isentropic(
                self.inl[0].to_flow(), self.outl[0].to_flow(),
                T0=self.inl[0].T.val_SI) - self.inl[0].h.val_SI) /
            (self.outl[0].h.val_SI - self.inl[0].h.val_SI))

        if self.char_map.is_set:
            # get bound errors for characteristic map
            x = np.sqrt(self.inl[0].T.design / self.inl[0].T.val_SI)
            y = (self.inl[0].m.val_SI * self.inl[0].p.design) / (
                self.inl[0].m.design * self.inl[0].p.val_SI * x)
            self.char_map.func.get_bound_errors(
                x, y, self.igva.val, self.label)

        if self.eta_s_char.is_set:
            # get bound errors for isentropic efficiency characteristics
            expr = 1
            if self.eta_s_char.param == 'm':
                if not np.isnan(self.inl[0].m.design):
                    expr = self.inl[0].m.val_SI / self.inl[0].m.design
            elif self.eta_s_char.param == 'pr':
                if not np.isnan(self.pr.design):
                    expr = self.pr.val / self.pr.design

            self.eta_s_char.func.get_bound_errors(expr, self.label)

        self.check_parameter_bounds()

    def exergy_balance(self, Tamb):
        r"""
        Calculate exergy balance of a compressor.

        Note
        ----
        .. math::

            \dot{E}_\mathrm{P} = \dot{m}_\mathrm{in} \cdot \left(
            e_\mathrm{ph,out} - e_\mathrm{ph,in} \right)\\
            \dot{E}_\mathrm{F} = P
        """
        self.E_P = self.inl[0].m.val_SI * (
            self.outl[0].ex_physical - self.inl[0].ex_physical)
        self.E_F = self.P.val
        self.E_D = self.E_F - self.E_P
        self.epsilon = self.E_P / self.E_F
