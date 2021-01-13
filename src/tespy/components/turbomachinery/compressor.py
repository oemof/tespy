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
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.fluid_properties import isentropic


class Compressor(Turbomachine):
    r"""
    Class for axial or radial compressor.

    **Mandatory Equations**

    - :py:meth:`tespy.components.component.Component.fluid_func`
    - :py:meth:`tespy.components.component.Component.mass_flow_func`

    **Optional Equations**

    - :py:meth:`tespy.components.component.Component.pr_func`
    - :py:meth:`tespy.components.turbomachinery.turbomachine.Turbomachine.energy_balance_func`
    - :py:meth:`tespy.components.turbomachinery.compressor.Compressor.eta_s_func`
    - :py:meth:`tespy.components.turbomachinery.compressor.Compressor.eta_s_char_func`
    - :py:meth:`tespy.components.turbomachinery.compressor.Compressor.char_map_func`

    Inlets/Outlets

    - in1
    - out1

    Image

    .. image:: _images/Compressor.svg
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

    def get_variables(self):
        return {
            'P': dc_cp(
                min_val=0, num_eq=1,
                deriv=self.energy_balance_deriv,
                func=self.energy_balance_func,
                latex=self.energy_balance_func_doc),
            'eta_s': dc_cp(
                min_val=0, max_val=1, num_eq=1,
                deriv=self.eta_s_deriv,
                func=self.eta_s_func, latex=self.eta_s_func_doc),
            'eta_s_char': dc_cc(
                param='m', num_eq=1,
                deriv=self.eta_s_char_deriv,
                func=self.eta_s_char_func, latex=self.eta_s_char_func_doc),
            'pr': dc_cp(
                min_val=1, num_eq=1,
                deriv=self.pr_deriv,
                func=self.pr_func, func_params={'pr': 'pr'},
                latex=self.pr_func_doc),
            'igva': dc_cp(min_val=-90, max_val=90, d=1e-3, val=0),
            'char_map': dc_cm(
                deriv=self.char_map_deriv, num_eq=2,
                func=self.char_map_func, latex=self.char_map_func_doc)
        }

    def comp_init(self, nw):
        Component.comp_init(self, nw)

        if self.char_map.char_func is None:
            self.char_map.char_func = ldc(
                self.component(), 'char_map', 'DEFAULT', CompressorMap)

    def eta_s_func(self, doc=False):
        r"""
        Equation for given isentropic efficiency of a compressor.

        Parameters
        ----------
        doc : boolean
            Return equation in LaTeX format instead of value.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = -\left( h_{out} - h_{in} \right) \cdot \eta_{s} +
                \left( h_{out,s} - h_{in} \right)
        """
        return (
            -(self.outl[0].h.val_SI - self.inl[0].h.val_SI) *
            self.eta_s.val + (isentropic(
                self.inl[0].to_flow(), self.outl[0].to_flow(),
                T0=self.inl[0].T.val_SI) - self.inl[0].h.val_SI))

    def eta_s_func_doc(self, label):
        r"""
        Equation for given isentropic efficiency of a compressor.

        Parameters
        ----------
        doc : boolean
            Return equation in LaTeX format instead of value.

        Returns
        -------
        residual : float
            Residual value of equation.
        """
        latex = (
            r'0 =-\left(h_\mathrm{out}-h_\mathrm{in}\right)\cdot'
            r'\eta_\mathrm{s}+\left(h_\mathrm{out,s}-h_\mathrm{in}\right)')
        return [self.generate_latex(latex, label)]

    def eta_s_deriv(self, increment_filter, k):
        r"""
        Partial derivatives for isentropic efficiency.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        f = self.eta_s_func
        if not increment_filter[0, 1]:
            self.jacobian[k, 0, 1] = self.numeric_deriv(f, 'p', 0)
        if not increment_filter[0, 1]:
            self.jacobian[k, 1, 1] = self.numeric_deriv(f, 'p', 1)
        if not increment_filter[0, 1]:
            self.jacobian[k, 0, 2] = self.numeric_deriv(f, 'h', 0)
        self.jacobian[k, 1, 2] = -self.eta_s.val

    def eta_s_char_func(self, doc=False):
        r"""
        Equation for given isentropic efficiency characteristic.

        Parameters
        ----------
        doc : boolean
            Return equation in LaTeX format instead of value.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \left(h_{out}-h_{in}\right) \cdot \eta_{s,design}
                \cdot f\left( expr \right) -\left( h_{out,s} - h_{in} \right)
        """
        p = self.eta_s_char.param
        expr = self.get_char_expr(p, **self.eta_s_char.char_params)
        if not expr:
            msg = ('Please choose a valid parameter, you want to link the '
                   'isentropic efficiency to at component ' + self.label + '.')
            logging.error(msg)
            raise ValueError(msg)

        i = self.inl[0]
        o = self.outl[0]
        return (
            self.eta_s.design * self.eta_s_char.char_func.evaluate(expr) *
            (o.h.val_SI - i.h.val_SI) - (isentropic(
                i.to_flow(), o.to_flow(), T0=self.inl[0].T.val_SI) -
             i.h.val_SI))

    def eta_s_char_func_doc(self, label):
        r"""
        Equation for given isentropic efficiency characteristic.

        Parameters
        ----------
        doc : boolean
            Return equation in LaTeX format instead of value.

        Returns
        -------
        residual : float
            Residual value of equation.
        """
        p = self.eta_s_char.param
        expr = self.get_char_expr_doc(p, **self.eta_s_char.char_params)
        if not expr:
            msg = ('Please choose a valid parameter, you want to link the '
                   'isentropic efficiency to at component ' + self.label + '.')
            logging.error(msg)
            raise ValueError(msg)

        latex = (
            r'0=\left(h_\mathrm{out}-h_\mathrm{in}\right)\cdot'
            r'\eta_\mathrm{s,design}\cdot f\left( ' + expr + r' \right)-'
            r'\left( h_{out,s} - h_{in} \right)')
        return [self.generate_latex(latex, label + '_' + p)]

    def eta_s_char_deriv(self, increment_filter, k):
        r"""
        Partial derivatives for isentropic efficiency characteristic.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
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

    def char_map_func(self):
        r"""
        Equation for characteristic map of compressor.

        Returns
        -------
        residual : ndarray
            Residual values of equations.

        Note
        ----
        - X: speedline index (rotational speed is constant)
        - Y: nondimensional mass flow
        - igva: variable inlet guide vane angle (assumed 0° if not
          specified)

        .. math::

            X = \sqrt{\frac{T_\mathrm{in,design}}{T_\mathrm{in}}}\\
            Y = \frac{\dot{m}_\mathrm{in} \cdot p_\mathrm{in,design}}
            {\dot{m}_\mathrm{in,design} \cdot p_\mathrm{in} \cdot X}\\
            0 = \frac{p_\mathrm{out} \cdot p_\mathrm{in,design}}
            {p_\mathrm{in} \cdot p_\mathrm{out,design}}-
            f\left(X,Y,igva\right)\\
            0 = \frac{\eta_\mathrm{s}}{\eta_\mathrm{s,design}} -
            f\left(X,Y,igva\right)
        """
        i = self.inl[0]
        o = self.outl[0]
        T = T_mix_ph(i.to_flow(), T0=self.inl[0].T.val_SI)

        x = np.sqrt(i.T.design / T)
        y = (i.m.val_SI * i.p.design) / (i.m.design * i.p.val_SI * x)

        pr, eta = self.char_map.char_func.evaluate(
            x, y, igva=self.igva.val)

        z1 = (o.p.val_SI / i.p.val_SI) / self.pr.design - pr
        z2 = (
            (isentropic(i.to_flow(), o.to_flow(), T0=T) -
             i.h.val_SI) / (o.h.val_SI - i.h.val_SI) / self.eta_s.design -
            eta)

        return np.array([z1, z2])

    def char_map_func_doc(self, label):
        r"""
        Equation for characteristic map of compressor.

        Parameters
        ----------
        doc : boolean
            Return equation in LaTeX format instead of value.

        Returns
        -------
        residual : ndarray
            Residual values of equations.
        """
        latex = (
            r'\begin{split}'
            r'X = &\sqrt{\frac{T_\mathrm{in,design}}{T_\mathrm{in}}}\\'
            '\n'
            r'Y = &\frac{\dot{m}_\mathrm{in} \cdot p_\mathrm{in,design}}'
            r'{\dot{m}_\mathrm{in,design} \cdot p_\mathrm{in} \cdot X}\\'
            '\n'
            r'0 = &\frac{p_\mathrm{out} \cdot p_\mathrm{in,design}}'
            r'{p_\mathrm{in} \cdot p_\mathrm{out,design}}-'
            r'f\left(X,Y,igva\right)\\' + '\n'
            r'0 = &\frac{\eta_\mathrm{s}}{\eta_\mathrm{s,design}} -'
            r'f\left(X,Y,igva\right)' + '\n'
            r'\end{split}'
        )
        return [self.generate_latex(latex, label), '']

    def char_map_deriv(self, increment_filter, k):
        r"""
        Partial derivatives for compressor map characteristic.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
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

        self.check_parameter_bounds()

    def exergy_balance(self, T0):
        r"""
        Calculate exergy balance of a compressor.

        Parameters
        ----------
        T0 : float
            Ambient temperature T0 / K.

        Note
        ----
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
            msg = ('Exergy balance of a compressor, where outlet temperature '
                   'is smaller than inlet temperature is not implmented.')
            logging.warning(msg)
            self.E_P = np.nan
            self.E_F = np.nan

        self.E_bus = self.P.val
        self.E_D = self.E_F - self.E_P
        self.epsilon = self.E_P / self.E_F
