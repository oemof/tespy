# -*- coding: utf-8

"""Module of class SolarCollector.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/heat_exchangers/solar_collector.py

SPDX-License-Identifier: MIT
"""

import logging

import numpy as np

from tespy.components.component import Component
from tespy.components.heat_exchangers.heat_exchanger_simple import HeatExchangerSimple
from tespy.tools.data_containers import dc_cp
from tespy.tools.data_containers import dc_gcp
from tespy.tools.data_containers import dc_simple
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.fluid_properties import s_mix_ph
from tespy.tools.fluid_properties import v_mix_ph
from tespy.tools.global_vars import err


class SolarCollector(HeatExchangerSimple):
    r"""
    The solar collector calculates heat output from irradiance.

    Equations

        **mandatory equations**

        - :py:meth:`tespy.components.component.Component.fluid_func`
        - :py:meth:`tespy.components.component.Component.mass_flow_func`

        **optional equations**

        .. math::

            0 = \dot{m}_{in} \cdot \left(h_{out} - h_{in} \right) -
            \dot{Q}

            0 = p_{in} \cdot pr - p_{out}

        - :py:meth:`tespy.components.component.Component.zeta_func`

        - :py:meth:`tespy.components.heat_exchangers.HeatExchangerSimple.darcy_func`
          or :py:meth:`tespy.components.heat_exchangers.HeatExchangerSimple.hw_func`

        **additional equations**

        - :py:meth:`tespy.components.heat_exchangers.SolarCollector.additional_equations`

    Inlets/Outlets

        - in1
        - out1

    Image

        .. image:: _images/SolarCollector.svg
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

    Q : str/float/tespy.tools.data_containers.dc_cp
        Heat transfer, :math:`Q/\text{W}`.

    pr : str/float/tespy.tools.data_containers.dc_cp
        Outlet to inlet pressure ratio, :math:`pr/1`.

    zeta : str/float/tespy.tools.data_containers.dc_cp
        Geometry independent friction coefficient,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    D : str/float/tespy.tools.data_containers.dc_cp
        Diameter of the pipes, :math:`D/\text{m}`.

    L : str/float/tespy.tools.data_containers.dc_cp
        Length of the pipes, :math:`L/\text{m}`.

    ks : str/float/tespy.tools.data_containers.dc_cp
        Pipe's roughness, :math:`ks/\text{m}` for darcy friction,
        :math:`ks/\text{1}` for hazen-williams equation.

    hydro_group : str/tespy.tools.data_containers.dc_gcp
        Parametergroup for pressure drop calculation based on pipes dimensions.
        Choose 'HW' for hazen-williams equation, else darcy friction factor is
        used.

    E : str/float/tespy.tools.data_containers.dc_cp
        irradiance at tilted collector surface area,
        :math:`E/\frac{\text{W}}{\text{m}^2}`.

    eta_opt : str/float/tespy.tools.data_containers.dc_cp
        optical loss at surface cover,
        :math:`\eta_{opt}`.

    lkf_lin : str/float/tespy.tools.data_containers.dc_cp
        Linear thermal loss key figure,
        :math:`\alpha_1/\frac{\text{W}}{\text{K} \cdot \text{m}^2}`.

    lkf_quad : str/float/tespy.tools.data_containers.dc_cp
        Quadratic thermal loss key figure,
        :math:`\alpha_2/\frac{\text{W}}{\text{K}^2 \cdot \text{m}^2}`.

    A : str/float/tespy.tools.data_containers.dc_cp
        Collector surface area :math:`A/\text{m}^2`.

    Tamb : float/tespy.tools.data_containers.dc_simple
        Ambient temperature, provide parameter in network's temperature unit.

    energy_group : tespy.tools.data_containers.dc_gcp
        Parametergroup for energy balance of solarthermal collector.

    Example
    -------
    The solar collector is used to calculate heat transferred to the heating
    system from irradiance on a tilted plane. For instance, it is possible to
    calculate the collector surface area required to transfer a specific amount
    of heat at a given irradiance. The collector parameters are the linear and
    the quadratic loss keyfigure as well as the optical effifiency.

    >>> from tespy.components import Sink, Source, SolarCollector
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> import shutil
    >>> fluids = ['H2O']
    >>> nw = Network(fluids=fluids)
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg', iterinfo=False)
    >>> so = Source('source')
    >>> si = Sink('sink')
    >>> sc = SolarCollector('solar collector')
    >>> sc.component()
    'solar collector'
    >>> sc.set_attr(pr=0.95, Q=1e4, design=['pr', 'Q'], offdesign=['zeta'],
    ...     Tamb=25, A='var', eta_opt=0.92, lkf_lin=1, lkf_quad=0.005, E=8e2)
    >>> inc = Connection(so, 'out1', sc, 'in1')
    >>> outg = Connection(sc, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)

    The outlet temperature should be at 90 °C at a constant mass flow, which
    is determined in the design calculation. In offdesign operation (at a
    different irradiance) using the calculated surface area and mass flow, it
    is possible to predict the outlet temperature. It would instead be
    possible to calulate the change in mass flow required to hold the
    specified outlet temperature, too.

    >>> inc.set_attr(fluid={'H2O': 1}, T=40, p=3, offdesign=['m'])
    >>> outg.set_attr(T=90, design=['T'])
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(sc.A.val, 1)
    14.5
    >>> sc.set_attr(A=sc.A.val, E=5e2, Tamb=20)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(sc.Q.val, 1)
    6083.8
    >>> round(outg.T.val, 1)
    70.5
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    @staticmethod
    def component():
        return 'solar collector'

    @staticmethod
    def attr():
        return {
            'Q': dc_cp(),
            'pr': dc_cp(min_val=1e-4, max_val=1), 'zeta': dc_cp(min_val=0),
            'D': dc_cp(min_val=1e-2, max_val=2, d=1e-4),
            'L': dc_cp(min_val=1e-1, d=1e-3),
            'ks': dc_cp(val=1e-4, min_val=1e-7, max_val=1e-4, d=1e-8),
            'E': dc_cp(min_val=0), 'A': dc_cp(min_val=0),
            'eta_opt': dc_cp(min_val=0, max_val=1),
            'lkf_lin': dc_cp(min_val=0), 'lkf_quad': dc_cp(min_val=0),
            'Tamb': dc_simple(),
            'Q_loss': dc_cp(min_val=0), 'SQ': dc_simple(),
            'hydro_group': dc_gcp(), 'energy_group': dc_gcp()
        }

    def comp_init(self, nw):

        Component.comp_init(self, nw)

        self.Tamb.val_SI = ((self.Tamb.val + nw.T[nw.T_unit][0]) *
                            nw.T[nw.T_unit][1])

        # parameters for hydro group
        self.hydro_group.set_attr(elements=[self.L, self.ks, self.D])

        is_set = True
        for e in self.hydro_group.elements:
            if e.is_set is False:
                is_set = False

        if is_set is True:
            self.hydro_group.set_attr(is_set=True)
        elif self.hydro_group.is_set is True:
            msg = ('All parameters of the component group have to be '
                   'specified! This component group uses the following '
                   'parameters: L, ks, D at ' + self.label + '. '
                   'Group will be set to False.')
            logging.warning(msg)
            self.hydro_group.set_attr(is_set=False)
        else:
            self.hydro_group.set_attr(is_set=False)

        # parameters for energy group
        self.energy_group.set_attr(
            elements=[self.E, self.eta_opt, self.lkf_lin,
                      self.lkf_quad, self.A, self.Tamb])

        is_set = True
        for e in self.energy_group.elements:
            if e.is_set is False:
                is_set = False

        if is_set is True:
            self.energy_group.set_attr(is_set=True)
        elif self.energy_group.is_set is True:
            msg = ('All parameters of the component group have to be '
                   'specified! This component group uses the following '
                   'parameters: E, eta_opt, lkf_lin, lkf_quad, A, Tamb at ' +
                   self.label + '. Group will be set to False.')
            logging.warning(msg)
            self.energy_group.set_attr(is_set=False)
        else:
            self.energy_group.set_attr(is_set=False)

        # number of mandatroy equations for
        # fluid balance: num_fl
        # mass flow: 1
        self.num_eq = self.num_nw_fluids + 1
        for var in [self.Q, self.pr, self.zeta, self.hydro_group,
                    self.energy_group]:
            if var.is_set is True:
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

            - :py:meth:`tespy.components.heat_exchangers.SolarCollector.energy_func`
        """
        ######################################################################
        # equation for specified energy-group paremeters
        if self.energy_group.is_set:
            if (np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0 or
                    self.always_all_equations):
                self.residual[k] = self.energy_func()

    def additional_derivatives(self, increment_filter, k):
        r"""Calculate partial derivatives for given additional equations."""
        ######################################################################
        # derivatives for specified energy-group paremeters
        if self.energy_group.is_set:
            f = self.energy_func
            self.jacobian[k, 0, 0] = (
                self.outl[0].h.val_SI - self.inl[0].h.val_SI)
            if not increment_filter[0, 1]:
                self.jacobian[k, 0, 1] = self.numeric_deriv(f, 'p', 0)
            if not increment_filter[0, 2]:
                self.jacobian[k, 0, 2] = self.numeric_deriv(f, 'h', 0)
            if not increment_filter[1, 1]:
                self.jacobian[k, 1, 1] = self.numeric_deriv(f, 'p', 1)
            if not increment_filter[1, 2]:
                self.jacobian[k, 1, 2] = self.numeric_deriv(f, 'h', 1)
            # custom variables for the energy-group
            for var in self.energy_group.elements:
                if var == self.Tamb:
                    continue
                if var.is_var:
                    self.jacobian[k, 2 + var.var_pos, 0] = (
                        self.numeric_deriv(f, self.vars[var], 2))
            k += 1

    def energy_func(self):
        r"""
        Equation for parabolic trough energy balance.

        Returns
        -------
        res : float
            Residual value of equation.

        Note
        ----
        .. math::

            \begin{split}
            T_m = & \frac{T_{out} + T_{in}}{2}\\
            0 = & \dot{m} \cdot \left( h_{out} - h_{in} \right)\\
            & - A \cdot \left[E \cdot \eta_{opt} - \alpha_1 \cdot
            \left(T_m - T_{amb} \right) - \alpha_2 \cdot
            \left(T_m - T_{amb}\right)^2 \right]
            \end{split}

        Reference: :cite:`Quaschning2013`.
        """
        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()

        T_m = (T_mix_ph(i, T0=self.inl[0].T.val_SI) +
               T_mix_ph(o, T0=self.outl[0].T.val_SI)) / 2

        return (i[0] * (o[2] - i[2]) -
                self.A.val * (
                    self.E.val * self.eta_opt.val -
                    (T_m - self.Tamb.val_SI) * self.lkf_lin.val -
                    self.lkf_quad.val * (T_m - self.Tamb.val_SI) ** 2))

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()

        self.SQ.val = i[0] * (s_mix_ph(o) - s_mix_ph(i))
        self.Q.val = i[0] * (o[2] - i[2])
        self.pr.val = o[1] / i[1]
        self.zeta.val = ((i[1] - o[1]) * np.pi ** 2 /
                         (8 * i[0] ** 2 * (v_mix_ph(i) + v_mix_ph(o)) / 2))
        if self.energy_group.is_set is True:
            self.Q_loss.val = self.E.val * self.A.val - self.Q.val

        self.check_parameter_bounds()
