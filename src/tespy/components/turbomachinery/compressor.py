# -*- coding: utf-8

"""Module of class Compressor.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/turbomachinery/compressor.py

SPDX-License-Identifier: MIT
"""

import numpy as np

from tespy.components.component import component_registry
from tespy.components.turbomachinery.base import Turbomachine
from tespy.tools import logger
from tespy.tools.data_containers import ComponentCharacteristicMaps as dc_cm
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.document_models import generate_latex_eq
from tespy.tools.fluid_properties import isentropic


@component_registry
class Compressor(Turbomachine):
    r"""
    Class for axial or radial compressor.

    **Mandatory Equations**

    - :py:meth:`tespy.components.component.Component.fluid_func`
    - :py:meth:`tespy.components.component.Component.mass_flow_func`

    **Optional Equations**

    - :py:meth:`tespy.components.component.Component.pr_func`
    - :py:meth:`tespy.components.turbomachinery.base.Turbomachine.energy_balance_func`
    - :py:meth:`tespy.components.turbomachinery.compressor.Compressor.eta_s_func`
    - :py:meth:`tespy.components.turbomachinery.compressor.Compressor.eta_s_char_func`
    - :py:meth:`tespy.components.turbomachinery.compressor.Compressor.char_map_eta_s_func`
    - :py:meth:`tespy.components.turbomachinery.compressor.Compressor.char_map_pr_func`

    Inlets/Outlets

    - in1
    - out1

    Image

    .. image:: /api/_images/Compressor.svg
       :alt: flowsheet of the compressor
       :align: center
       :class: only-light

    .. image:: /api/_images/Compressor_darkmode.svg
       :alt: flowsheet of the compressor
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

    P : float, dict
        Power, :math:`P/\text{W}`

    eta_s : float, dict
        Isentropic efficiency, :math:`\eta_s/1`

    pr : float, dict, :code:`"var"`
        Outlet to inlet pressure ratio, :math:`pr/1`

    eta_s_char : tespy.tools.characteristics.CharLine, dict
        Characteristic curve for isentropic efficiency, provide CharLine as
        function :code:`func`.

    char_map : tespy.tools.characteristics.CharMap, dict
        Characteristic map for pressure rise and isentropic efficiency vs.
        nondimensional mass flow, see
        :py:class:`tespy.tools.characteristics.CharMap` for further
        information. Provide a CompressorMap as function :code:`func`.

    igva : float, dict, :code:`"var"`
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
    >>> nw = Network(p_unit='bar', T_unit='C', h_unit='kJ / kg', v_unit='l / s',
    ... iterinfo=False)
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
    ... offdesign=['char_map_pr', 'char_map_eta_s'])
    >>> inc.set_attr(fluid={'air': 1}, p=1, T=20, v=50)
    >>> nw.solve('design')
    >>> nw.save('tmp.json')
    >>> round(comp.P.val, 0)
    12772.0
    >>> inc.set_attr(v=45)
    >>> comp.set_attr(igva='var')
    >>> nw.solve('offdesign', design_path='tmp.json')
    >>> round(comp.eta_s.val, 2)
    0.77
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    @staticmethod
    def component():
        return 'compressor'

    def get_parameters(self):
        parameters = super().get_parameters()
        parameters.update({
            'eta_s': dc_cp(
                min_val=0, max_val=1, num_eq=1,
                deriv=self.eta_s_deriv,
                func=self.eta_s_func, latex=self.eta_s_func_doc),
            'eta_s_char': dc_cc(
                param='m', num_eq=1,
                deriv=self.eta_s_char_deriv,
                func=self.eta_s_char_func, latex=self.eta_s_char_func_doc),
            'igva': dc_cp(min_val=-90, max_val=90, d=1e-3, val=0),
            'char_map_eta_s': dc_cm(),
            'char_map_eta_s_group': dc_gcp(
                elements=['char_map_eta_s', 'igva'], num_eq=1,
                latex=self.char_map_eta_s_func_doc,
                func=self.char_map_eta_s_func,
                deriv=self.char_map_eta_s_deriv),
            'char_map_pr': dc_cm(),
            'char_map_pr_group': dc_gcp(
                elements=['char_map_pr', 'igva'],
                deriv=self.char_map_pr_deriv, num_eq=1,
                func=self.char_map_pr_func, latex=self.char_map_pr_func_doc)
        })
        return parameters

    def eta_s_func(self):
        r"""
        Equation for given isentropic efficiency of a compressor.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = -\left( h_{out} - h_{in} \right) \cdot \eta_{s} +
                \left( h_{out,s} - h_{in} \right)
        """
        i = self.inl[0]
        o = self.outl[0]
        return (
            (o.h.val_SI - i.h.val_SI) * self.eta_s.val - (
                isentropic(
                    i.p.val_SI,
                    i.h.val_SI,
                    o.p.val_SI,
                    i.fluid_data,
                    i.mixing_rule,
                    T0=None
                ) - self.inl[0].h.val_SI
            )
        )

    def eta_s_func_doc(self, label):
        r"""
        Equation for given isentropic efficiency of a compressor.

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
            r'0 =-\left(h_\mathrm{out}-h_\mathrm{in}\right)\cdot'
            r'\eta_\mathrm{s}+\left(h_\mathrm{out,s}-h_\mathrm{in}\right)')
        return generate_latex_eq(self, latex, label)

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
        i = self.inl[0]
        o = self.outl[0]
        f = self.eta_s_func
        if self.is_variable(i.p, increment_filter):
            self.jacobian[k, i.p.J_col] = self.numeric_deriv(f, 'p', i)
        if self.is_variable(o.p, increment_filter):
            self.jacobian[k, o.p.J_col] = self.numeric_deriv(f, 'p', o)
        if self.is_variable(i.h, increment_filter):
            self.jacobian[k, i.h.J_col] = self.numeric_deriv(f, 'h', i)
        if self.is_variable(o.h, increment_filter):
            self.jacobian[k, o.h.J_col] = self.eta_s.val

    def eta_s_char_func(self):
        r"""
        Equation for given isentropic efficiency characteristic.

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
            logger.error(msg)
            raise ValueError(msg)

        i = self.inl[0]
        o = self.outl[0]
        return (
            (o.h.val_SI - i.h.val_SI)
            * self.eta_s.design * self.eta_s_char.char_func.evaluate(expr)
            - (
                isentropic(
                    i.p.val_SI,
                    i.h.val_SI,
                    o.p.val_SI,
                    i.fluid_data,
                    i.mixing_rule,
                    T0=None
                ) - i.h.val_SI
            )
        )

    def eta_s_char_func_doc(self, label):
        r"""
        Equation for given isentropic efficiency characteristic.

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
            r'0=\left(h_\mathrm{out}-h_\mathrm{in}\right)\cdot'
            r'\eta_\mathrm{s,design}\cdot f\left(X\right)-'
            r'\left( h_{out,s} - h_{in} \right)')
        return generate_latex_eq(self, latex, label)

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
        i = self.inl[0]
        o = self.outl[0]
        if self.is_variable(i.m, increment_filter):
            self.jacobian[k, i.m.J_col] = self.numeric_deriv(f, 'm', i)
        if self.is_variable(i.p, increment_filter):
            self.jacobian[k, i.p.J_col] = self.numeric_deriv(f, 'p', i)
        if self.is_variable(i.h, increment_filter):
            self.jacobian[k, i.h.J_col] = self.numeric_deriv(f, 'h', i)
        if self.is_variable(o.p, increment_filter):
            self.jacobian[k, o.p.J_col] = self.numeric_deriv(f, 'p', o)
        if self.is_variable(o.h, increment_filter):
            self.jacobian[k, o.h.J_col] = self.numeric_deriv(f, 'h', o)

    def char_map_pr_func(self):
        r"""
        Calculate pressure ratio from characteristic map.

        Returns
        -------
        residual : float
            Residual value of equations.

        Note
        ----
        - X: speedline index (rotational speed is constant)
        - Y: nondimensional mass flow
        - igva: variable inlet guide vane angle for value manipulation
          according to :cite:`GasTurb2018`.

        .. math::

            X = \sqrt{\frac{T_\mathrm{in,design}}{T_\mathrm{in}}}\\
            Y = \frac{\dot{m}_\mathrm{in} \cdot p_\mathrm{in,design}}
            {\dot{m}_\mathrm{in,design} \cdot p_\mathrm{in} \cdot X}\\
            \vec{Y} = f\left(X,Y\right)\cdot\left(1-\frac{igva}{100}\right)\\
            \vec{Z} = f\left(X,Y\right)\cdot\left(1-\frac{igva}{100}\right)\\
            0 = \frac{p_{out} \cdot p_{in,design}}
            {p_\mathrm{in} \cdot p_\mathrm{out,design}}-
            f\left(Y,\vec{Y},\vec{Z}\right)
        """
        i = self.inl[0]
        o = self.outl[0]

        x = np.sqrt(i.T.design / i.calc_T())
        y = (i.m.val_SI * i.p.design) / (i.m.design * i.p.val_SI * x)

        yarr, zarr = self.char_map_pr.char_func.evaluate_x(x)
        # value manipulation with igva
        yarr *= (1 - self.igva.val / 100)
        zarr *= (1 - self.igva.val / 100)
        pr = self.char_map_pr.char_func.evaluate_y(y, yarr, zarr)

        return (o.p.val_SI / i.p.val_SI) - pr * self.pr.design

    def char_map_pr_func_doc(self, label):
        r"""
        Get LaTeX equation for pressure ratio from characteristic map.

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
            r'X = &\sqrt{\frac{T_\mathrm{in,design}}{T_\mathrm{in}}}\\'
            '\n'
            r'Y = &\frac{\dot{m}_\mathrm{in} \cdot p_\mathrm{in,design}}'
            r'{\dot{m}_\mathrm{in,design} \cdot p_\mathrm{in} \cdot X}\\'
            '\n'
            r'\vec{Y} = &f\left(X,Y\right)\cdot\left(1-\frac{igva}{100}\right)'
            r'\\' + '\n'
            r'\vec{Z} = &'
            r'f\left(X,Y\right)\cdot\left(1-\frac{igva}{100}\right)\\' + '\n'
            r'0 = &\frac{p_\mathrm{out} \cdot p_\mathrm{in,design}}'
            r'{p_\mathrm{in} \cdot p_\mathrm{out,design}}-'
            r'f\left(Y,\vec{Y},\vec{Z}\right)\\' + '\n'
            r'\end{split}'
        )
        return generate_latex_eq(self, latex, label)

    def char_map_pr_deriv(self, increment_filter, k):
        r"""
        Partial derivatives for compressor map characteristic.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        f = self.char_map_pr_func
        i = self.inl[0]
        o = self.outl[0]
        if self.is_variable(i.m, increment_filter):
            self.jacobian[k, i.m.J_col] = self.numeric_deriv(f, 'm', i)
        if self.is_variable(i.p, increment_filter):
            self.jacobian[k, i.p.J_col] = self.numeric_deriv(f, 'p', i)
        if self.is_variable(i.h, increment_filter):
            self.jacobian[k, i.h.J_col] = self.numeric_deriv(f, 'h', i)
        if self.is_variable(o.p, increment_filter):
            self.jacobian[k, o.p.J_col] = 1 / i.p.val_SI
        if self.is_variable(o.h, increment_filter):
            self.jacobian[k, o.h.J_col] = self.numeric_deriv(f, 'h', o)

        if self.igva.is_var:
            self.jacobian[k, self.igva.J_col] = self.numeric_deriv(
                f, 'igva', None
            )

    def char_map_eta_s_func(self):
        r"""
        Calculate isentropic efficiency from characteristic map.

        Returns
        -------
        residual : float
            Residual value of equation.

        Note
        ----
        - X: speedline index (rotational speed is constant)
        - Y: nondimensional mass flow
        - igva: variable inlet guide vane angle for value manipulation
          according to :cite:`GasTurb2018`.

        .. math::

            X = \sqrt{\frac{T_\mathrm{in,design}}{T_\mathrm{in}}}\\
            Y = \frac{\dot{m}_\mathrm{in} \cdot p_\mathrm{in,design}}
            {\dot{m}_\mathrm{in,design} \cdot p_\mathrm{in} \cdot X}\\
            \vec{Y} = f\left(X,Y\right)\cdot\left(1-\frac{igva}{100}\right)\\
            \vec{Z}=f\left(X,Y\right)\cdot\left(1-\frac{igva^2}{10000}\right)\\
            0 = \frac{\eta_\mathrm{s}}{\eta_\mathrm{s,design}} -
            f\left(Y,\vec{Y},\vec{Z}\right)
        """
        i = self.inl[0]
        o = self.outl[0]

        x = np.sqrt(i.T.design / i.calc_T())
        y = (i.m.val_SI * i.p.design) / (i.m.design * i.p.val_SI * x)

        yarr, zarr = self.char_map_eta_s.char_func.evaluate_x(x)
        # value manipulation with igva
        yarr *= (1 - self.igva.val / 100)
        zarr *= (1 - self.igva.val ** 2 / 10000)
        eta = self.char_map_eta_s.char_func.evaluate_y(y, yarr, zarr)

        return (
            (
            isentropic(
                i.p.val_SI,
                i.h.val_SI,
                o.p.val_SI,
                i.fluid_data,
                i.mixing_rule,
                T0=i.T.val_SI
            ) - i.h.val_SI)
            / (o.h.val_SI - i.h.val_SI) - eta * self.eta_s.design
        )

    def char_map_eta_s_func_doc(self, label):
        r"""
        Get LaTeX equation for isentropic efficiency from characteristic map.

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
            r'\begin{split}'
            r'X = &\sqrt{\frac{T_\mathrm{in,design}}{T_\mathrm{in}}}\\'
            '\n'
            r'Y = &\frac{\dot{m}_\mathrm{in} \cdot p_\mathrm{in,design}}'
            r'{\dot{m}_\mathrm{in,design} \cdot p_\mathrm{in} \cdot X}\\'
            '\n'
            r'\vec{Y} = &f\left(X,Y\right)\cdot\left(1-\frac{igva}{100}\right)'
            r'\\' + '\n'
            r'\vec{Z} = &'
            r'f\left(X,Y\right)\cdot\left(1-\frac{igva^2}{10000}\right)\\'
            '\n'
            r'0 = &\frac{\eta_\mathrm{s}}{\eta_\mathrm{s,design}} -'
            r'f\left(Y,\vec{Y},\vec{Z}\right)' + '\n'
            r'\end{split}'
        )
        return generate_latex_eq(self, latex, label)

    def char_map_eta_s_deriv(self, increment_filter, k):
        r"""
        Partial derivatives for compressor map characteristic.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        f = self.char_map_eta_s_func
        i = self.inl[0]
        o = self.outl[0]
        if self.is_variable(i.m, increment_filter):
            self.jacobian[k, i.m.J_col] = self.numeric_deriv(f, 'm', i)
        if self.is_variable(i.p, increment_filter):
            self.jacobian[k, i.p.J_col] = self.numeric_deriv(f, 'p', i)
        if self.is_variable(i.h, increment_filter):
            self.jacobian[k, i.h.J_col] = self.numeric_deriv(f, 'h', i)
        if self.is_variable(o.p, increment_filter):
            self.jacobian[k, o.p.J_col] = self.numeric_deriv(f, 'p', o)
        if self.is_variable(o.h, increment_filter):
            self.jacobian[k, o.h.J_col] = self.numeric_deriv(f, 'h', o)

        if self.igva.is_var:
            self.jacobian[k, self.igva.J_col] = self.numeric_deriv(
                f, 'igva', None
            )

    def convergence_check(self):
        r"""
        Perform a convergence check.

        Note
        ----
        Manipulate enthalpies/pressure at inlet and outlet if not specified by
        user to match physically feasible constraints.
        """
        i, o = self.inl, self.outl

        if o[0].p.is_var and o[0].p.val_SI < i[0].p.val_SI:
            o[0].p.val_SI = i[0].p.val_SI * 1.1

        if o[0].h.is_var and o[0].h.val_SI < i[0].h.val_SI:
            o[0].h.val_SI = i[0].h.val_SI * 1.1

        if i[0].p.is_var and o[0].p.val_SI < i[0].p.val_SI:
            i[0].p.val_SI = o[0].p.val_SI * 0.9

        if i[0].h.is_var and o[0].h.val_SI < i[0].h.val_SI:
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
        super().calc_parameters()

        i = self.inl[0]
        o = self.outl[0]
        self.eta_s.val =  (
            isentropic(
                i.p.val_SI,
                i.h.val_SI,
                o.p.val_SI,
                i.fluid_data,
                i.mixing_rule,
                T0=None
            ) - self.inl[0].h.val_SI
        ) / (o.h.val_SI - i.h.val_SI)

    def check_parameter_bounds(self):
        r"""Check parameter value limits."""
        super().check_parameter_bounds()

        for data in [self.char_map_pr, self.char_map_eta_s]:
            if data.is_set:
                x = np.sqrt(self.inl[0].T.design / self.inl[0].T.val_SI)
                y = (self.inl[0].m.val_SI * self.inl[0].p.design) / (
                    self.inl[0].m.design * self.inl[0].p.val_SI * x)
                yarr = data.char_func.get_domain_errors_x(x, self.label)
                yarr *= (1 - self.igva.val / 100)
                data.char_func.get_domain_errors_y(y, yarr, self.label)

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
            logger.warning(msg)
            self.E_P = np.nan
            self.E_F = np.nan

        self.E_bus = {
            "chemical": 0, "physical": 0, "massless": self.P.val
        }
        self.E_D = self.E_F - self.E_P
        self.epsilon = self._calc_epsilon()
