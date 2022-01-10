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
from tespy.tools.data_containers import ComponentCharacteristicMaps as dc_cm
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.document_models import generate_latex_eq
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.fluid_properties import isentropic


class DiabaticCompressor(Turbomachine):
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
    - :py:meth:`tespy.components.turbomachinery.compressor.Compressor.char_map_eta_s_func`
    - :py:meth:`tespy.components.turbomachinery.compressor.Compressor.char_map_pr_func`

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
    ... offdesign=['char_map_pr', 'char_map_eta_s'])
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
            'Q': dc_cp(max_val=0, is_result=True),
            'eta_s': dc_cp(
                min_val=0, max_val=1, num_eq=1,
                deriv=self.eta_s_deriv,
                func=self.eta_s_func, latex=self.eta_s_func_doc),
            'pr': dc_cp(
                min_val=1, num_eq=1,
                deriv=self.pr_deriv,
                func=self.pr_func, func_params={'pr': 'pr'},
                latex=self.pr_func_doc),
        }

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
        return (
            -(self.outl[0].h.val_SI - self.inl[0].h.val_SI) *
            self.eta_s.val + (isentropic(
                self.inl[0].get_flow(), self.outl[0].get_flow(),
                T0=self.inl[0].T.val_SI) - self.inl[0].h.val_SI))

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
        f = self.eta_s_func
        if not increment_filter[0, 1]:
            self.jacobian[k, 0, 1] = self.numeric_deriv(f, 'p', 0)
        if not increment_filter[0, 1]:
            self.jacobian[k, 1, 1] = self.numeric_deriv(f, 'p', 1)
        if not increment_filter[0, 1]:
            self.jacobian[k, 0, 2] = self.numeric_deriv(f, 'h', 0)
        self.jacobian[k, 1, 2] = -self.eta_s.val

    def energy_balance_func(self):
        return self.inl[0].m.val_SI * (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI) - self.P.val - self.Q.val

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
                self.inl[0].get_flow(), self.outl[0].get_flow(),
                T0=self.inl[0].T.val_SI) - self.inl[0].h.val_SI) /
            (self.outl[0].h.val_SI - self.inl[0].h.val_SI))

        self.P.val = self.inl[0].m.val_SI * (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI) / 0.95
        self.Q.val = -0.05 * self.P.val

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
