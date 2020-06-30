# -*- coding: utf-8

"""Module for components of type heat exchanger.

Components in this module:

- :func:`tespy.components.heat_exchangers.condenser`
- :func:`tespy.components.heat_exchangers.desuperheater`
- :func:`tespy.components.heat_exchangers.heat_exchanger`
- :func:`tespy.components.heat_exchangers.heat_exchanger_simple`
- :func:`tespy.components.heat_exchangers.parabolic_trough`
- :func:`tespy.components.heat_exchangers.solar_collector`


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/heat_exchangers.py

SPDX-License-Identifier: MIT
"""

import logging
import warnings

import numpy as np

from tespy.components.components import component
from tespy.tools.data_containers import dc_cc
from tespy.tools.data_containers import dc_cp
from tespy.tools.data_containers import dc_gcp
from tespy.tools.data_containers import dc_simple
from tespy.tools.fluid_properties import T_bp_p
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.fluid_properties import dh_mix_dpQ
from tespy.tools.fluid_properties import h_mix_pQ
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.fluid_properties import s_mix_ph
from tespy.tools.fluid_properties import v_mix_ph
from tespy.tools.fluid_properties import visc_mix_ph
from tespy.tools.global_vars import err
from tespy.tools.helpers import lamb

# %%


class heat_exchanger_simple(component):
    r"""
    A basic heat exchanger representing a heat source or heat sink.

    The component heat_exchanger_simple is the parent class for the components
    solar_collector and pipe (of module piping).

    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_func`
        - :func:`tespy.components.components.component.mass_flow_func`

        **optional equations**

        .. math::

            0 = \dot{m}_{in} \cdot \left(h_{out} - h_{in} \right) -
            \dot{Q}

            0 = p_{in} \cdot pr - p_{out}

        - :func:`tespy.components.components.component.zeta_func`

        - :func:`tespy.components.heat_exchangers.heat_exchanger_simple.darcy_func`
          or :func:`tespy.components.heat_exchangers.heat_exchanger_simple.hw_func`

        **additional equations**

        - :func:`tespy.components.heat_exchangers.heat_exchanger_simple.additional_equations`

    Inlets/Outlets

        - in1
        - out1

    Image

        .. image:: _images/pipe.svg
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

    kA : str/float/tespy.tools.data_containers.dc_cp
        Area independent heat transition coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char : tespy.tools.charactersitics.char_line/tespy.tools.data_containers.dc_cc
        Characteristic line for heat transfer coefficient.

    Tamb : float/tespy.tools.data_containers.dc_simple
        Ambient temperature, provide parameter in network's temperature unit.

    kA_group : tespy.tools.data_containers.dc_gcp
        Parametergroup for heat transfer calculation from ambient temperature
        and area independent heat transfer coefficient kA.

    Example
    -------
    The heat_exchanger_simple can be used as a sink or source of heat. This
    component does not simulate the secondary side of the heat exchanger. It
    is possible to calculate the pressure ratio with the Darcy-Weisbach
    equation or in case of liquid water use the Hazen-Williams equation.
    Also, given ambient temperature and the heat transfer coeffiecient, it is
    possible to predict heat transfer.

    >>> from tespy.components import sink, source, heat_exchanger_simple
    >>> from tespy.connections import connection
    >>> from tespy.networks import network
    >>> import shutil
    >>> fluids = ['N2']
    >>> nw = network(fluids=fluids)
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg', iterinfo=False)
    >>> so1 = source('source 1')
    >>> si1 = sink('sink 1')
    >>> heat_sink = heat_exchanger_simple('heat sink')
    >>> heat_sink.component()
    'heat exchanger simple'
    >>> heat_sink.set_attr(Tamb=10, pr=0.95, design=['pr'],
    ... offdesign=['zeta', 'kA_char'])
    >>> inc = connection(so1, 'out1', heat_sink, 'in1')
    >>> outg = connection(heat_sink, 'out1', si1, 'in1')
    >>> nw.add_conns(inc, outg)

    It is possible to determine the amount of heat transferred when the fluid
    enters the heat sink at a temperature of 200 °C and is cooled down to
    150 °C. Given an ambient temperature of 10 °C this also determines the heat
    transfer coefficient to the ambient. Assuming a characteristic function
    for the heat transfer coefficient we can predict the heat transferred at
    variable flow rates.

    >>> inc.set_attr(fluid={'N2': 1}, m=1, T=200, p=5)
    >>> outg.set_attr(T=150, design=['T'])
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(heat_sink.Q.val, 0)
    -52581.0
    >>> round(heat_sink.kA.val, 0)
    321.0
    >>> inc.set_attr(m=1.25)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(heat_sink.Q.val, 0)
    -56599.0
    >>> round(outg.T.val, 1)
    156.9
    >>> inc.set_attr(m=0.75)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(heat_sink.Q.val, 1)
    -47275.8
    >>> round(outg.T.val, 1)
    140.0
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    @staticmethod
    def component():
        return 'heat exchanger simple'

    @staticmethod
    def attr():
        return {
            'Q': dc_cp(),
            'pr': dc_cp(min_val=1e-4, max_val=1), 'zeta': dc_cp(min_val=0),
            'D': dc_cp(min_val=1e-2, max_val=2, d=1e-4),
            'L': dc_cp(min_val=1e-1, d=1e-3),
            'ks': dc_cp(val=1e-4, min_val=1e-7, max_val=1e-3, d=1e-8),
            'kA': dc_cp(min_val=0, d=1),
            'kA_char': dc_cc(param='m'), 'Tamb': dc_simple(),
            'SQ1': dc_simple(), 'SQ2': dc_simple(), 'Sirr': dc_simple(),
            'hydro_group': dc_gcp(), 'kA_group': dc_gcp(),
            'kA_char_group': dc_gcp()
        }

    @staticmethod
    def inlets():
        return ['in1']

    @staticmethod
    def outlets():
        return ['out1']

    def comp_init(self, nw):

        component.comp_init(self, nw)

        self.Tamb.val_SI = ((self.Tamb.val + nw.T[nw.T_unit][0]) *
                            nw.T[nw.T_unit][1])

        # parameters for hydro group
        self.hydro_group.set_attr(elements=[self.L, self.ks, self.D])

        is_set = self.L.is_set and self.ks.is_set and self.D.is_set

        if is_set:
            self.hydro_group.set_attr(is_set=True)
            if self.hydro_group.method == 'HW':
                method = 'Hazen-Williams equation'
            else:
                method = 'darcy friction factor'
            msg = (
                'Pressure loss calculation from pipe dimensions method is set '
                'to ' + method + '.')
            logging.debug(msg)

        elif self.hydro_group.is_set:
            msg = (
                'All parameters of the component group have to be specified! '
                'This component group uses the following parameters: L, ks, D '
                'at ' + self.label + '. Group will be set to False.')
            logging.warning(msg)
            self.hydro_group.set_attr(is_set=False)
        else:
            self.hydro_group.set_attr(is_set=False)

        # parameters for kA group
        self.kA_group.set_attr(elements=[self.kA, self.Tamb])

        is_set = self.kA.is_set and self.Tamb.is_set

        if is_set:
            self.kA_group.set_attr(is_set=True)
        elif self.kA_group.is_set:
            msg = (
                'All parameters of the component group have to be specified! '
                'This component group uses the following parameters: kA, Tamb '
                'at ' + self.label + '. Group will be set to False.')
            logging.warning(msg)
            self.kA_group.set_attr(is_set=False)
        else:
            self.kA_group.set_attr(is_set=False)

        # parameters for kA_char group
        self.kA_char_group.set_attr(elements=[self.kA_char, self.Tamb])

        is_set = self.kA_char.is_set and self.Tamb.is_set

        if is_set:
            self.kA_char_group.set_attr(is_set=True)
        elif self.kA_char_group.is_set:
            msg = (
                'All parameters of the component group have to be specified! '
                'This component group uses the following parameters: kA_char, '
                'Tamb at ' + self.label + '. Group will be set to False.')
            logging.warning(msg)
            self.kA_char_group.set_attr(is_set=False)
        else:
            self.kA_char_group.set_attr(is_set=False)

        # number of mandatroy equations for
        # fluid balance: num_fl
        # mass flow: 1
        self.num_eq = self.num_nw_fluids + 1
        for var in [self.Q, self.pr, self.zeta, self.hydro_group,
                    self.kA_group, self.kA_char_group]:
            if var.is_set is True:
                self.num_eq += 1

        if self.kA.is_set:
            msg = (
                'The usage of the parameter kA has changed for offdesign '
                'calculation. Specifying kA will keep a constant value for kA '
                'in the calculation. If you want to use the value adaption of '
                'kA by the characteristic line, please use kA_char as '
                'parameter instead (occurred at ' + self.label + '). This '
                'warning will disappear in TESPy version 0.4.0.')
            warnings.warn(msg, FutureWarning, stacklevel=2)

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
        # equations for fluid balance
        self.residual[k:k + self.num_nw_fluids] = self.fluid_func()
        k += self.num_nw_fluids

        ######################################################################
        # equations for mass flow balance
        self.residual[k] = self.mass_flow_func()
        k += 1

        ######################################################################
        # equations for specified heta transfer
        if self.Q.is_set:
            self.residual[k] = self.inl[0].m.val_SI * (
                self.outl[0].h.val_SI - self.inl[0].h.val_SI) - self.Q.val
            k += 1

        ######################################################################
        # equations for specified pressure ratio
        if self.pr.is_set:
            self.residual[k] = (
                self.inl[0].p.val_SI * self.pr.val - self.outl[0].p.val_SI)
            k += 1

        ######################################################################
        # equations for specified zeta
        if self.zeta.is_set:
            if np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0:
                self.residual[k] = self.zeta_func(zeta='zeta')
            k += 1

        ######################################################################
        # equation for specified hydro-group paremeters
        if self.hydro_group.is_set:
            if np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0:
                # hazen williams equation
                if self.hydro_group.method == 'HW':
                    func = self.hw_func
                # darcy friction factor
                else:
                    func = self.darcy_func
                self.residual[k] = func()
            k += 1

        ######################################################################
        # additional equations
        self.additional_equations(k)

    def additional_equations(self, k):
        r"""
        Calculate results of additional equations.

        Equations

            **optional equations**

            - :func:`tespy.components.heat_exchangers.heat_exchanger_simple.kA_func`
            - :func:`tespy.components.heat_exchangers.heat_exchanger_simple.kA_char_func`
        """
        ######################################################################
        # equation for specified kA_group paremeters
        if self.kA_group.is_set:
            if np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0:
                self.residual[k] = self.kA_func()
            k += 1

        ######################################################################
        # equation for specified kA_char_group paremeters
        if self.kA_char_group.is_set:
            if np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0:
                self.residual[k] = self.kA_char_func()
            k += 1

    def derivatives(self, increment_filter):
        r"""Calculate partial derivatives for given equations."""
        ######################################################################
        # derivatives fluid and mass balance are static
        k = self.num_nw_fluids + 1

        ######################################################################
        # derivatives for specified heat transfer
        if self.Q.is_set:
            self.jacobian[k, 0, 0] = (
                self.outl[0].h.val_SI - self.inl[0].h.val_SI)
            self.jacobian[k, 0, 2] = -self.inl[0].m.val_SI
            self.jacobian[k, 1, 2] = self.inl[0].m.val_SI
            # custom variable Q
            if self.Q.is_var:
                self.jacobian[k, 2 + self.Q.var_pos, 0] = -1
            k += 1

        ######################################################################
        # derivatives for specified pressure ratio
        if self.pr.is_set:
            self.jacobian[k, 0, 1] = self.pr.val
            self.jacobian[k, 1, 1] = -1
            # custom variable pr
            if self.pr.is_var:
                self.jacobian[k, 2 + self.pr.var_pos, 0] = (
                    self.inl[0].p.val_SI)
            k += 1

        ######################################################################
        # derivatives for specified zeta
        if self.zeta.is_set:
            f = self.zeta_func
            if not increment_filter[0, 0]:
                self.jacobian[k, 0, 0] = self.numeric_deriv(
                    f, 'm', 0, zeta='zeta')
            if not increment_filter[0, 2]:
                self.jacobian[k, 0, 1] = self.numeric_deriv(
                    f, 'p', 0, zeta='zeta')
            if not increment_filter[0, 2]:
                self.jacobian[k, 0, 2] = self.numeric_deriv(
                    f, 'h', 0, zeta='zeta')
            if not increment_filter[1, 1]:
                self.jacobian[k, 1, 1] = self.numeric_deriv(
                    f, 'p', 1, zeta='zeta')
            if not increment_filter[1, 2]:
                self.jacobian[k, 1, 2] = self.numeric_deriv(
                    f, 'h', 1, zeta='zeta')
            # custom variable zeta
            if self.zeta.is_var:
                self.jacobian[k, 2 + self.zeta.var_pos, 0] = (
                    self.numeric_deriv(f, 'zeta', 2, zeta='zeta'))
            k += 1

        ######################################################################
        # derivatives for specified hydro-group parameters
        if self.hydro_group.is_set:
            # hazen williams equation
            if self.hydro_group.method == 'HW':
                func = self.hw_func
            # darcy friction factor
            else:
                func = self.darcy_func

            if not increment_filter[0, 0]:
                self.jacobian[k, 0, 0] = self.numeric_deriv(func, 'm', 0)
            if not increment_filter[0, 1]:
                self.jacobian[k, 0, 1] = self.numeric_deriv(func, 'p', 0)
            if not increment_filter[0, 2]:
                self.jacobian[k, 0, 2] = self.numeric_deriv(func, 'h', 0)
            if not increment_filter[1, 1]:
                self.jacobian[k, 1, 1] = self.numeric_deriv(func, 'p', 1)
            if not increment_filter[1, 2]:
                self.jacobian[k, 1, 2] = self.numeric_deriv(func, 'h', 1)
            # custom variables of hydro group
            for var in self.hydro_group.elements:
                if var.is_var:
                    self.jacobian[k, 2 + var.var_pos, 0] = (
                        self.numeric_deriv(func, self.vars[var], 2))
            k += 1

        ######################################################################
        # derivatives for additional equations
        self.additional_derivatives(increment_filter, k)

    def additional_derivatives(self, increment_filter, k):
        r"""Calculategit partial derivatives for given additional equations."""
        ######################################################################
        # derivatives for specified kA_group paremeters
        if self.kA_group.is_set:
            f = self.kA_func
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
            if self.kA.is_var:
                self.jacobian[k, 2 + self.kA.var_pos, 0] = (
                    self.numeric_deriv(f, self.vars[self.kA], 2))
            k += 1

        ######################################################################
        # derivatives for specified kA_char_group paremeters

        if self.kA_char_group.is_set:
            f = self.kA_char_func
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

    def darcy_func(self):
        r"""
        Equation for pressure drop calculation from darcy friction factor.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                Re = \frac{4 \cdot |\dot{m}_{in}|}{\pi \cdot D \cdot
                \frac{\eta_{in}+\eta_{out}}{2}}\\

                0 = p_{in} - p_{out} - \frac{8 \cdot |\dot{m}_{in}| \cdot
                \dot{m}_{in} \cdot \frac{v_{in}+v_{out}}{2} \cdot L \cdot
                \lambda\left(Re, ks, D\right)}{\pi^2 \cdot D^5}\\

                \eta: \text{dynamic viscosity}\\
                v: \text{specific volume}\\
                \lambda: \text{darcy friction factor}
        """
        i, o = self.inl[0].to_flow(), self.outl[0].to_flow()

        if abs(i[0]) < 1e-4:
            return i[1] - o[1]

        visc_i = visc_mix_ph(i, T0=self.inl[0].T.val_SI)
        visc_o = visc_mix_ph(o, T0=self.outl[0].T.val_SI)
        v_i = v_mix_ph(i, T0=self.inl[0].T.val_SI)
        v_o = v_mix_ph(o, T0=self.outl[0].T.val_SI)

        re = 4 * abs(i[0]) / (np.pi * self.D.val * (visc_i + visc_o) / 2)

        return ((i[1] - o[1]) - 8 * abs(i[0]) * i[0] * (v_i + v_o) / 2 *
                self.L.val * lamb(re, self.ks.val, self.D.val) /
                (np.pi ** 2 * self.D.val ** 5))

    def hw_func(self):
        r"""
        Equation for pressure drop calculation from Hazen-Williams equation.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                res = \left(p_{in} - p_{out} \right) \cdot \left(-1\right)^i -
                \frac{10.67 \cdot |\dot{m}_{in}| ^ {1.852}
                \cdot L}{ks^{1.852} \cdot D^{4.871}} \cdot g \cdot
                \left(\frac{v_{in} + v_{out}}{2}\right)^{0.852}

                i = \begin{cases}
                0 & \dot{m}_{in} \geq 0\\
                1 & \dot{m}_{in} < 0
                \end{cases}

        Note
        ----
        Gravity g is set to :math:`9.81 \frac{m}{s^2}`
        """
        i, o = self.inl[0].to_flow(), self.outl[0].to_flow()

        if abs(i[0]) < 1e-4:
            return i[1] - o[1]

        v_i = v_mix_ph(i, T0=self.inl[0].T.val_SI)
        v_o = v_mix_ph(o, T0=self.outl[0].T.val_SI)
        flow_dir = np.sign(i[0])

        return ((i[1] - o[1]) * flow_dir -
                (10.67 * abs(i[0]) ** 1.852 * self.L.val /
                 (self.ks.val ** 1.852 * self.D.val ** 4.871)) *
                (9.81 * ((v_i + v_o) / 2) ** 0.852))

    def kA_func(self):
        r"""
        Calculate heat transfer from heat transfer coefficient.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                ttd_u = \begin{cases}
                T_{in} - T_{amb} & \dot{m} \geq 0\\
                T_{out} - T_{amb} & \dot{m} < 0
                \end{cases}

                ttd_l = \begin{cases}
                T_{in} - T_{amb} & \dot{m} < 0\\
                T_{out} - T_{amb} & \dot{m} \geq 0
                \end{cases}

                0 = \dot{m}_{in} \cdot \left( h_{out} - h_{in}\right) +
                kA \cdot \frac{ttd_u - ttd_l}
                {\ln{\frac{ttd_u}{ttd_l}}}

                T_{amb}: \text{ambient temperature}
        """
        i, o = self.inl[0].to_flow(), self.outl[0].to_flow()

        ttd_1 = T_mix_ph(i, T0=self.inl[0].T.val_SI) - self.Tamb.val_SI
        ttd_2 = T_mix_ph(o, T0=self.outl[0].T.val_SI) - self.Tamb.val_SI

        if ttd_1 > ttd_2:
            td_log = (ttd_1 - ttd_2) / np.log(ttd_1 / ttd_2)
        elif ttd_1 < ttd_2:
            td_log = (ttd_2 - ttd_1) / np.log(ttd_2 / ttd_1)
        else:
            td_log = 0

        return i[0] * (o[2] - i[2]) + self.kA.val * td_log

    def kA_char_func(self):
        r"""
        Calculate heat transfer from heat transfer coefficient characteristic.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                ttd_u = \begin{cases}
                T_{in} - T_{amb} & \dot{m} \geq 0\\
                T_{out} - T_{amb} & \dot{m} < 0
                \end{cases}

                ttd_l = \begin{cases}
                T_{in} - T_{amb} & \dot{m} < 0\\
                T_{out} - T_{amb} & \dot{m} \geq 0
                \end{cases}

                0 = \dot{m}_{in} \cdot \left( h_{out} - h_{in}\right) +
                kA_{ref} \cdot f_{kA} \cdot \frac{ttd_u - ttd_l}
                {\ln{\frac{ttd_u}{ttd_l}}}

                f_{kA} = \frac{2}{1 + \frac{1}
                {f\left(\frac{m_1}{m_{1,ref}}\right)}}

                T_{amb}: \text{ambient temperature}

        Note
        ----
        For standard function of f\ :subscript:`1` \ see module
        :func:`tespy.data`.
        """
        i, o = self.inl[0].to_flow(), self.outl[0].to_flow()

        ttd_1 = T_mix_ph(i, T0=self.inl[0].T.val_SI) - self.Tamb.val_SI
        ttd_2 = T_mix_ph(o, T0=self.outl[0].T.val_SI) - self.Tamb.val_SI

        if ttd_1 > ttd_2:
            td_log = (ttd_1 - ttd_2) / np.log(ttd_1 / ttd_2)
        elif ttd_1 < ttd_2:
            td_log = (ttd_2 - ttd_1) / np.log(ttd_2 / ttd_1)
        else:
            td_log = 0

        f = 1
        if not np.isnan(self.inl[0].m.design):
            if self.kA_char.param == 'm':
                f = self.kA_char.func.evaluate(i[0] / self.inl[0].m.design)

        fkA = 2 / (1 + 1 / f)

        return i[0] * (o[2] - i[2]) + self.kA.design * fkA * td_log

    def bus_func(self, bus):
        r"""
        Calculate the value of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus
            TESPy bus object.

        Returns
        -------
        val : float
            Value of energy transfer :math:`\dot{E}`. This value is passed to
            :py:meth:`tespy.components.components.component.calc_bus_value`
            for value manipulation according to the specified characteristic
            line of the bus.

            .. math::

                \dot{E} = \dot{m}_{in} \cdot \left( h_{out} - h_{in} \right)
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
        bus : tespy.connections.bus
            TESPy bus object.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 2, self.num_nw_vars))
        f = self.calc_bus_value
        deriv[0, 0, 0] = self.numeric_deriv(f, 'm', 0, bus=bus)
        deriv[0, 0, 2] = self.numeric_deriv(f, 'h', 0, bus=bus)
        deriv[0, 1, 2] = self.numeric_deriv(f, 'h', 1, bus=bus)
        return deriv

    def initialise_source(self, c, key):
        r"""
        Return a starting value for pressure and enthalpy the outlets.

        Parameters
        ----------
        c : tespy.connections.connection
            Connection to perform initialisation on.

        key : str
            Fluid property to retrieve.

        Returns
        -------
        val : float
            Starting value for pressure/enthalpy in SI units.

            .. math::

                val = \begin{cases}
                \begin{cases}
                1 \cdot 10^5 \; \frac{\text{J}}{\text{kg}} & \dot{Q} < 0\\
                3 \cdot 10^5 \; \frac{\text{J}}{\text{kg}} & \dot{Q} = 0\\
                5 \cdot 10^5 \; \frac{\text{J}}{\text{kg}} & \dot{Q} > 0
                \end{cases} & \text{key = 'h'}\\
                \; \; \; \; 10^5 \text{Pa} & \text{key = 'p'}
                \end{cases}

        """
        if key == 'p':
            return 1e5
        elif key == 'h':
            if self.Q.val < 0 and self.Q.is_set:
                return 1e5
            elif self.Q.val > 0 and self.Q.is_set:
                return 5e5
            else:
                return 3e5

    def initialise_target(self, c, key):
        r"""
        Return a starting value for pressure and enthalpy the inlets.

        Parameters
        ----------
        c : tespy.connections.connection
            Connection to perform initialisation on.

        key : str
            Fluid property to retrieve.

        Returns
        -------
        val : float
            Starting value for pressure/enthalpy in SI units.

            .. math::

                val = \begin{cases}
                1 \cdot 10^5 & \text{key = 'p'}\\
                \begin{cases}
                5 \cdot 10^5 & \dot{Q} < 0\\
                3 \cdot 10^5 & \dot{Q} = 0\\
                1 \cdot 10^5 & \dot{Q} > 0
                \end{cases} & \text{key = 'h'}\\
                \end{cases}
        """
        if key == 'p':
            return 1e5
        elif key == 'h':
            if self.Q.val < 0 and self.Q.is_set:
                return 5e5
            elif self.Q.val > 0 and self.Q.is_set:
                return 1e5
            else:
                return 3e5

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()
        v_i = v_mix_ph(i, T0=self.inl[0].T.val_SI)
        v_o = v_mix_ph(o, T0=self.outl[0].T.val_SI)

        self.SQ1.val = i[0] * (s_mix_ph(o) - s_mix_ph(i))
        self.Q.val = i[0] * (o[2] - i[2])
        self.pr.val = o[1] / i[1]
        self.zeta.val = ((i[1] - o[1]) * np.pi ** 2 /
                         (8 * i[0] ** 2 * (v_i + v_o) / 2))

        if self.Tamb.is_set:
            self.SQ2.val = -i[0] * (o[2] - i[2]) / self.Tamb.val_SI
            self.Sirr.val = self.SQ1.val + self.SQ2.val

            ttd_1 = T_mix_ph(i, T0=self.inl[0].T.val_SI) - self.Tamb.val_SI
            ttd_2 = T_mix_ph(o, T0=self.outl[0].T.val_SI) - self.Tamb.val_SI

            if ttd_1 > ttd_2:
                td_log = (ttd_1 - ttd_2) / np.log(ttd_1 / ttd_2)
            elif ttd_1 < ttd_2:
                td_log = (ttd_2 - ttd_1) / np.log(ttd_2 / ttd_1)
            else:
                td_log = 0

            self.kA.val = abs(i[0] * (o[2] - i[2]) / td_log)

        if self.kA.is_set:
            # get bound errors for kA characteristic line
            if self.kA_char.param == 'm':
                self.kA_char.func.get_bound_errors(i[0] / self.inl[0].m.design,
                                                   self.label)

        self.check_parameter_bounds()

# %%


class parabolic_trough(heat_exchanger_simple):
    r"""
    The parabolic trough calculates heat output from irradiance.

    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_func`
        - :func:`tespy.components.components.component.mass_flow_func`

        **optional equations**

        .. math::

            0 = \dot{m}_{in} \cdot \left(h_{out} - h_{in} \right) -
            \dot{Q}

            0 = p_{in} \cdot pr - p_{out}

        - :func:`tespy.components.components.component.zeta_func`

        - :func:`tespy.components.heat_exchangers.heat_exchanger_simple.darcy_func`
          or :func:`tespy.components.heat_exchangers.heat_exchanger_simple.hw_func`

        **additional equations**

        - :func:`tespy.components.heat_exchangers.parabolic_trough.additional_equations`

    Inlets/Outlets

        - in1
        - out1

    Image

        .. image:: _images/parabolic_trough.svg
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
        Diameter of the absorber tube, :math:`D/\text{m}`.

    L : str/float/tespy.tools.data_containers.dc_cp
        Length of the absorber tube, :math:`L/\text{m}`.

    ks : str/float/tespy.tools.data_containers.dc_cp
        Tube's roughness, :math:`ks/\text{m}` for darcy friction,
        :math:`ks/\text{1}` for hazen-williams equation.

    hydro_group : str/tespy.tools.data_containers.dc_gcp
        Parametergroup for pressure drop calculation based on pipes dimensions.
        Choose 'HW' for hazen-williams equation, else darcy friction factor is
        used.

    E : str/float/tespy.tools.data_containers.dc_cp
        Direct irradiance to tilted collector,
        :math:`E/\frac{\text{W}}{\text{m}^2}`.

    aoi : str/float/tespy.tools.data_containers.dc_cp
        Angle of incidience, :math:`aoi/^\circ`.

    doc : str/float/tespy.tools.data_containers.dc_cp
        Degree of cleanliness (1: full absorption, 0: no absorption),
        :math:`X`.

    eta_opt : str/float/tespy.tools.data_containers.dc_cp
        (constant) optical losses due to surface reflection,
        :math:`\eta_{opt}`.

    c_1 : str/float/tespy.tools.data_containers.dc_cp
        Linear thermal loss key figure,
        :math:`c_1/\frac{\text{W}}{\text{K} \cdot \text{m}^2}`.

    c_2 : str/float/tespy.tools.data_containers.dc_cp
        Quadratic thermal loss key figure,
        :math:`c_2/\frac{\text{W}}{\text{K}^2 \cdot \text{m}^2}`.

    iam_1 : str/float/tespy.tools.data_containers.dc_cp
        Linear incidence angle modifier,
        :math:`iam_1/\frac{1}{^\circ}`.

    iam_2 : str/float/tespy.tools.data_containers.dc_cp
        Quadratic incidence angle modifier,
        :math:`iam_2/\left(\frac{1}{^\circ}\right)^2`.

    A : str/float/tespy.tools.data_containers.dc_cp
        Collector aperture surface area :math:`A/\text{m}^2`.

    Tamb : float/tespy.tools.data_containers.dc_simple
        Ambient temperature, provide parameter in network's temperature unit.

    energy_group : tespy.tools.data_containers.dc_gcp
        Parametergroup for energy balance of solarthermal collector.

    Example
    -------
    A parabolic trough is installed using S800 as thermo-fluid.
    First, the operation conditions from :cite:`Janotte2014` are reproduced.
    Therefore, the direct normal irradiance :math:`\dot{E}_\mathrm{DNI}` is at
    1000 :math:`\frac{\text{W}}{\text{m}^2}` at an angle of incidence
    :math:`aoi` at 20 °. This means, the direct irradiance to the parabolic
    trough :math:`E` is at
    :math:`\dot{E}_{DNI} \cdot cos\left(20^\circ\right)`.

    >>> from tespy.components import sink, source, parabolic_trough
    >>> from tespy.connections import connection
    >>> from tespy.networks import network
    >>> import numpy as np
    >>> import shutil
    >>> fluids = ['INCOMP::S800']
    >>> nw = network(fluids=fluids)
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg', iterinfo=False)
    >>> so = source('source')
    >>> si = sink('sink')
    >>> pt = parabolic_trough('parabolic trough collector')
    >>> pt.component()
    'parabolic trough'
    >>> inc = connection(so, 'out1', pt, 'in1')
    >>> outg = connection(pt, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)

    The pressure ratio is at a constant level of 1. However, it is possible to
    specify the pressure losses from the absorber tube length, roughness and
    diameter, too. The aperture surface :math:`A` is specified to 1
    :math:`\text{m}^2` for simplicity reasons.

    >>> aoi = 20
    >>> E = 1000 * np.cos(aoi / 180 * np.pi)
    >>> pt.set_attr(pr=1, aoi=aoi, doc=1,
    ... Tamb=20, A=1, eta_opt=0.816, c_1=0.0622, c_2=0.00023, E=E,
    ... iam_1=-1.59e-3, iam_2=9.77e-5)
    >>> inc.set_attr(fluid={'S800': 1}, T=220, p=2)
    >>> outg.set_attr(T=260)
    >>> nw.solve('design')
    >>> round(pt.Q.val, 0)
    736.0

    For example, it is possible to calculate the aperture area of the parabolic
    trough given the total heat production, outflow temperature and mass flow.

    >>> pt.set_attr(A='var', Q=5e6, Tamb=25)
    >>> inc.set_attr(T=None)
    >>> outg.set_attr(T=350, m=20)
    >>> nw.solve('design')
    >>> round(inc.T.val)
    229.0
    >>> round(pt.A.val)
    6862.0

    Given this design, it is possible to calculate the outlet temperature as
    well as the heat transfer at different operating points.

    >>> aoi = 30
    >>> E = 800 * np.cos(aoi / 180 * np.pi)
    >>> pt.set_attr(A=pt.A.val, aoi=aoi, Q=None, E=E)
    >>> inc.set_attr(T=150)
    >>> outg.set_attr(T=None)
    >>> nw.solve('design')
    >>> round(outg.T.val)
    244.0
    >>> round(pt.Q.val)
    3603027.0
    """

    @staticmethod
    def component():
        return 'parabolic trough'

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
            'c_1': dc_cp(min_val=0), 'c_2': dc_cp(min_val=0),
            'iam_1': dc_cp(), 'iam_2': dc_cp(),
            'aoi': dc_cp(min_val=-90, max_val=90),
            'doc': dc_cp(min_val=0, max_val=1),
            'Tamb': dc_simple(),
            'Q_loss': dc_cp(min_val=0), 'SQ': dc_simple(),
            'hydro_group': dc_gcp(), 'energy_group': dc_gcp()
        }

    def comp_init(self, nw):

        component.comp_init(self, nw)

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
            msg = (
                'All parameters of the component group have to be specified! '
                'This component group uses the following parameters: L, ks, D '
                'at ' + self.label + '. Group will be set to False.')
            logging.warning(msg)
            self.hydro_group.set_attr(is_set=False)
        else:
            self.hydro_group.set_attr(is_set=False)

        # parameters for energy group
        self.energy_group.set_attr(
            elements=[
                self.E, self.eta_opt, self.aoi, self.doc, self.c_1, self.c_2,
                self.iam_1, self.iam_2, self.A, self.Tamb])

        is_set = True
        for e in self.energy_group.elements:
            if e.is_set is False:
                is_set = False

        if is_set is True:
            self.energy_group.set_attr(is_set=True)
        elif self.energy_group.is_set is True:
            msg = (
                'All parameters of the component group have to be specified! '
                'This component group uses the following parameters: E, '
                'eta_opt, aoi, doc, c_1, c_2, iam_1, iam_2, A, Tamb at ' +
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

            - :func:`tespy.components.heat_exchangers.parabolic_trough.energy_func`
        """
        ######################################################################
        # equation for specified energy-group paremeters
        if self.energy_group.is_set:
            if np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0:
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
        Equation for solar collector energy balance.

        Returns
        -------
        res : float
            Residual value of equation.

        Note
        ----
        .. math::

            \begin{split}
            T_m = & \frac{T_{out} + T_{in}}{2}\\
            iam = & 1 - iam_1 \cdot |aoi| - iam_2 \cdot aoi^2\\
            0 = & \dot{m} \cdot \left( h_{out} - h_{in} \right)\\
            & - A \cdot \left[E \cdot \eta_{opt} \cdot doc^{1.5} \cdot
            iam \right. \\
            & \left. - c_1 \cdot \left(T_m - T_{amb} \right) -
            c_2 \cdot \left(T_m - T_{amb}\right)^2
            \vphantom{ \eta_{opt} \cdot doc^{1.5}} \right]
            \end{split}

        Reference: :cite:`Janotte2014`.
        """
        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()

        T_m = (T_mix_ph(i, T0=self.inl[0].T.val_SI) +
               T_mix_ph(o, T0=self.outl[0].T.val_SI)) / 2

        iam = (
            1 - self.iam_1.val * abs(self.aoi.val) -
            self.iam_2.val * self.aoi.val ** 2)

        return (i[0] * (o[2] - i[2]) -
                self.A.val * (
                    self.E.val * self.eta_opt.val * self.doc.val ** 1.5 * iam -
                    (T_m - self.Tamb.val_SI) * self.c_1.val -
                    self.c_2.val * (T_m - self.Tamb.val_SI) ** 2))

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

# %%


class solar_collector(heat_exchanger_simple):
    r"""
    The solar collector calculates heat output from irradiance.

    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_func`
        - :func:`tespy.components.components.component.mass_flow_func`

        **optional equations**

        .. math::

            0 = \dot{m}_{in} \cdot \left(h_{out} - h_{in} \right) -
            \dot{Q}

            0 = p_{in} \cdot pr - p_{out}

        - :func:`tespy.components.components.component.zeta_func`

        - :func:`tespy.components.heat_exchangers.heat_exchanger_simple.darcy_func`
          or :func:`tespy.components.heat_exchangers.heat_exchanger_simple.hw_func`

        **additional equations**

        - :func:`tespy.components.heat_exchangers.solar_collector.additional_equations`

    Inlets/Outlets

        - in1
        - out1

    Image

        .. image:: _images/solar_collector.svg
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

    >>> from tespy.components import sink, source, solar_collector
    >>> from tespy.connections import connection
    >>> from tespy.networks import network
    >>> import shutil
    >>> fluids = ['H2O']
    >>> nw = network(fluids=fluids)
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg', iterinfo=False)
    >>> so = source('source')
    >>> si = sink('sink')
    >>> sc = solar_collector('solar collector')
    >>> sc.component()
    'solar collector'
    >>> sc.set_attr(pr=0.95, Q=1e4, design=['pr', 'Q'], offdesign=['zeta'],
    ...     Tamb=25, A='var', eta_opt=0.92, lkf_lin=1, lkf_quad=0.005, E=8e2)
    >>> inc = connection(so, 'out1', sc, 'in1')
    >>> outg = connection(sc, 'out1', si, 'in1')
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

        component.comp_init(self, nw)

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

            - :func:`tespy.components.heat_exchangers.solar_collector.energy_func`
        """
        ######################################################################
        # equation for specified energy-group paremeters
        if self.energy_group.is_set:
            if np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0:
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

# %%


class heat_exchanger(component):
    r"""
    Class heat_exchanger is the parent class for condenser and desuperheater.

    The heat exchanger represents counter current heat exchangers. Both, hot
    and cold side of the heat exchanger, are simulated.

    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_func`
        - :func:`tespy.components.heat_exchangers.heat_exchanger.mass_flow_func`

        - :func:`tespy.components.heat_exchangers.heat_exchanger.energy_func`

        **optional equations**

        .. math::

            0 = \dot{m}_{in} \cdot \left(h_{out} - h_{in} \right) - \dot{Q}

        - :func:`tespy.components.heat_exchangers.heat_exchanger.kA_func`
        - :func:`tespy.components.heat_exchangers.condenser.kA_char_func`
        - :func:`tespy.components.heat_exchangers.heat_exchanger.ttd_u_func`
        - :func:`tespy.components.heat_exchangers.heat_exchanger.ttd_l_func`

        .. math::

            0 = p_{1,in} \cdot pr1 - p_{1,out}\\
            0 = p_{2,in} \cdot pr2 - p_{2,out}

        - hot side :func:`tespy.components.components.component.zeta_func`
        - cold side :func:`tespy.components.components.component.zeta_func`

        **additional equations**

        - :func:`tespy.components.heat_exchangers.heat_exchanger.additional_equations`

    Inlets/Outlets

        - in1, in2 (index 1: hot side, index 2: cold side)
        - out1, out2 (index 1: hot side, index 2: cold side)

    Image

        .. image:: _images/heat_exchanger.svg
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

    pr1 : str/float/tespy.tools.data_containers.dc_cp
        Outlet to inlet pressure ratio at hot side, :math:`pr/1`.

    pr2 : str/float/tespy.tools.data_containers.dc_cp
        Outlet to inlet pressure ratio at cold side, :math:`pr/1`.

    zeta1 : str/float/tespy.tools.data_containers.dc_cp
        Geometry independent friction coefficient at hot side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    zeta2 : str/float/tespy.tools.data_containers.dc_cp
        Geometry independent friction coefficient at cold side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    kA : float/tespy.tools.data_containers.dc_cp
        Area independent heat transition coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char : tespy.tools.data_containers.dc_simple
        Area independent heat transition coefficient characteristic.

    kA_char1 : tespy.tools.charactersitics.char_line/tespy.tools.data_containers.dc_cc
        Characteristic line for hot side heat transfer coefficient.

    kA_char2 : tespy.tools.charactersitics.char_line/tespy.tools.data_containers.dc_cc
        Characteristic line for cold side heat transfer coefficient.

    Note
    ----
    The heat exchanger and subclasses (desuperheater, condenser) are
    countercurrent heat exchangers. Equations (kA, ttd_u, ttd_l) do not work
    for directcurrent and crosscurrent or combinations of different types.

    Example
    -------
    A water cooling is installed to transfer heat from hot exhaust air. The
    heat exchanger is designed for a terminal temperature difference of 5 K.
    From this, it is possible to calculate the heat transfer coefficient and
    predict water and air outlet temperature in offdesign operation.

    >>> from tespy.components import sink, source, heat_exchanger
    >>> from tespy.connections import connection
    >>> from tespy.networks import network
    >>> import shutil
    >>> nw = network(fluids=['water', 'air'], T_unit='C', p_unit='bar',
    ... h_unit='kJ / kg', iterinfo=False)
    >>> exhaust_hot = source('Exhaust air outlet')
    >>> exhaust_cold = sink('Exhaust air inlet')
    >>> cw_cold = source('cooling water inlet')
    >>> cw_hot = sink('cooling water outlet')
    >>> he = heat_exchanger('waste heat exchanger')
    >>> he.component()
    'heat exchanger'
    >>> ex_he = connection(exhaust_hot, 'out1', he, 'in1')
    >>> he_ex = connection(he, 'out1', exhaust_cold, 'in1')
    >>> cw_he = connection(cw_cold, 'out1', he, 'in2')
    >>> he_cw = connection(he, 'out2', cw_hot, 'in1')
    >>> nw.add_conns(ex_he, he_ex, cw_he, he_cw)

    The volumetric flow of the air is at 100 l/s. After designing the component
    it is possible to predict the temperature at different flow rates or
    different inlet temperatures of the exhaust air.

    >>> he.set_attr(pr1=0.98, pr2=0.98, ttd_u=5,
    ... design=['pr1', 'pr2', 'ttd_u'], offdesign=['zeta1', 'zeta2', 'kA_char'])
    >>> cw_he.set_attr(fluid={'air': 0, 'water': 1}, T=10, p=3,
    ... offdesign=['m'])
    >>> ex_he.set_attr(fluid={'air': 1, 'water': 0}, v=0.1, T=35)
    >>> he_ex.set_attr(T=17.5, p=1, design=['T'])
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(ex_he.T.val - he_cw.T.val, 0)
    5.0
    >>> ex_he.set_attr(v=0.075)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(he_cw.T.val, 1)
    27.5
    >>> round(he_ex.T.val, 1)
    14.4
    >>> ex_he.set_attr(v=0.1, T=40)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(he_cw.T.val, 1)
    33.9
    >>> round(he_ex.T.val, 1)
    18.8
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    @staticmethod
    def component():
        return 'heat exchanger'

    @staticmethod
    def attr():
        return {
            'Q': dc_cp(max_val=0),
            'kA': dc_cp(min_val=0),
            'td_log': dc_cp(min_val=0),
            'ttd_u': dc_cp(min_val=0), 'ttd_l': dc_cp(min_val=0),
            'pr1': dc_cp(max_val=1), 'pr2': dc_cp(max_val=1),
            'zeta1': dc_cp(min_val=0), 'zeta2': dc_cp(min_val=0),
            'kA_char': dc_simple(),
            'kA_char1': dc_cc(param='m'), 'kA_char2': dc_cc(param='m'),
            'SQ1': dc_simple(), 'SQ2': dc_simple(), 'Sirr': dc_simple()
        }

    @staticmethod
    def inlets():
        return ['in1', 'in2']

    @staticmethod
    def outlets():
        return ['out1', 'out2']

    def comp_init(self, nw):

        component.comp_init(self, nw)

        # number of mandatroy equations for
        # fluid balance: num_fl * 2
        # mass flow: 2
        # energy balance: 1
        self.num_eq = self.num_nw_fluids * 2 + 3
        for var in [self.Q, self.kA, self.kA_char, self.ttd_u, self.ttd_l,
                    self.pr1, self.pr2, self.zeta1, self.zeta2]:
            if var.is_set is True:
                self.num_eq += 1

        if self.kA.is_set:
            msg = (
                'The usage of the parameter kA has changed for offdesign '
                'calculation. Specifying kA will keep a constant value for kA '
                'in the calculation. If you want to use the value adaption of '
                'kA by the characteristic line, please use kA_char as '
                'parameter instead (occurred at ' + self.label + '). This '
                'warning will disappear in TESPy version 0.4.0.')
            warnings.warn(msg, FutureWarning, stacklevel=2)

        self.jacobian = np.zeros((
            self.num_eq,
            self.num_i + self.num_o + self.num_vars,
            self.num_nw_vars))

        self.residual = np.zeros(self.num_eq)
        pos = self.num_nw_fluids * 2
        self.jacobian[0:pos] = self.fluid_deriv()
        self.jacobian[pos:pos + 2] = self.mass_flow_deriv()

    def equations(self):
        r"""Calculate residual vector with results of equations."""
        k = 0
        ######################################################################
        # equations for fluid balance
        self.residual[k:k + self.num_nw_fluids * 2] = self.fluid_func()
        k += self.num_nw_fluids * 2

        ######################################################################
        # equations for mass flow balance
        self.residual[k:k + 2] = self.mass_flow_func()
        k += 2

        ######################################################################
        # equations for energy balance
        self.residual[k] = self.energy_func()
        k += 1

        ######################################################################
        # equations for specified heat transfer
        if self.Q.is_set:
            self.residual[k] = (
                self.inl[0].m.val_SI * (
                    self.outl[0].h.val_SI - self.inl[0].h.val_SI) - self.Q.val)
            k += 1

        ######################################################################
        # equations for specified heat transfer coefficient
        if self.kA.is_set:
            if np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0:
                self.residual[k] = self.kA_func()
            k += 1

        ######################################################################
        # equations for specified heat transfer coefficient characteristic
        if self.kA_char.is_set:
            if np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0:
                self.residual[k] = self.kA_char_func()
            k += 1

        ######################################################################
        # equations for specified upper terminal temperature difference
        if self.ttd_u.is_set:
            self.residual[k] = self.ttd_u_func()
            k += 1

        ######################################################################
        # equations for specified lower terminal temperature difference
        if self.ttd_l.is_set:
            self.residual[k] = self.ttd_l_func()
            k += 1

        ######################################################################
        # equations for specified pressure ratio at hot side
        if self.pr1.is_set:
            self.residual[k] = (
                self.pr1.val * self.inl[0].p.val_SI - self.outl[0].p.val_SI)
            k += 1

        ######################################################################
        # equations for specified pressure ratio at cold side
        if self.pr2.is_set:
            self.residual[k] = (
                self.pr2.val * self.inl[1].p.val_SI - self.outl[1].p.val_SI)
            k += 1

        ######################################################################
        # equations for specified zeta at hot side
        if self.zeta1.is_set:
            if np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0:
                self.residual[k] = self.zeta_func(
                    zeta='zeta1', inconn=0, outconn=0)
            k += 1

        ######################################################################
        # equations for specified zeta at cold side
        if self.zeta2.is_set:
            if np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0:
                self.residual[k] = self.zeta_func(
                    zeta='zeta2', inconn=1, outconn=1)
            k += 1

        ######################################################################
        # additional equations
        self.additional_equations(k)

    def additional_equations(self, k):
        r"""Calculate results of additional equations."""
        return

    def derivatives(self, increment_filter):
        r"""
        Calculate partial derivatives for given equations.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        ######################################################################
        # derivatives fluid and mass balance are static
        k = self.num_nw_fluids * 2 + 2

        ######################################################################
        # derivatives for energy balance equation
        for i in range(2):
            self.jacobian[k, i, 0] = (
                self.outl[i].h.val_SI - self.inl[i].h.val_SI)
            self.jacobian[k, i, 2] = -self.inl[i].m.val_SI

        self.jacobian[k, 2, 2] = self.inl[0].m.val_SI
        self.jacobian[k, 3, 2] = self.inl[1].m.val_SI
        k += 1

        ######################################################################
        # derivatives for specified heat transfer
        if self.Q.is_set:
            self.jacobian[k, 0, 0] = (
                self.outl[0].h.val_SI - self.inl[0].h.val_SI)
            self.jacobian[k, 0, 2] = -self.inl[0].m.val_SI
            self.jacobian[k, 2, 2] = self.inl[0].m.val_SI
            k += 1

        ######################################################################
        # derivatives for specified heat transfer coefficient
        if self.kA.is_set:
            f = self.kA_func
            self.jacobian[k, 0, 0] = (
                self.outl[0].h.val_SI - self.inl[0].h.val_SI)
            for i in range(4):
                if not increment_filter[i, 1]:
                    self.jacobian[k, i, 1] = self.numeric_deriv(f, 'p', i)
                if not increment_filter[i, 2]:
                    self.jacobian[k, i, 2] = self.numeric_deriv(f, 'h', i)
            k += 1

        ######################################################################
        # derivatives for specified heat transfer coefficient
        if self.kA_char.is_set:
            f = self.kA_char_func
            if not increment_filter[0, 0]:
                self.jacobian[k, 0, 0] = self.numeric_deriv(f, 'm', 0)
            if not increment_filter[1, 0]:
                self.jacobian[k, 1, 0] = self.numeric_deriv(f, 'm', 1)
            for i in range(4):
                if not increment_filter[i, 1]:
                    self.jacobian[k, i, 1] = self.numeric_deriv(f, 'p', i)
                if not increment_filter[i, 2]:
                    self.jacobian[k, i, 2] = self.numeric_deriv(f, 'h', i)
            k += 1

        ######################################################################
        # derivatives for specified upper terminal temperature difference
        if self.ttd_u.is_set:
            f = self.ttd_u_func
            for i in [0, 3]:
                if not increment_filter[i, 1]:
                    self.jacobian[k, i, 1] = self.numeric_deriv(f, 'p', i)
                if not increment_filter[i, 2]:
                    self.jacobian[k, i, 2] = self.numeric_deriv(f, 'h', i)
            k += 1

        ######################################################################
        # derivatives for specified lower terminal temperature difference
        if self.ttd_l.is_set:
            f = self.ttd_l_func
            for i in [1, 2]:
                if not increment_filter[i, 1]:
                    self.jacobian[k, i, 1] = self.numeric_deriv(f, 'p', i)
                if not increment_filter[i, 2]:
                    self.jacobian[k, i, 2] = self.numeric_deriv(f, 'h', i)
            k += 1

        ######################################################################
        # derivatives for specified pressure ratio at hot side
        if self.pr1.is_set:
            self.jacobian[k, 0, 1] = self.pr1.val
            self.jacobian[k, 2, 1] = -1
            k += 1

        ######################################################################
        # derivatives for specified pressure ratio at cold side
        if self.pr2.is_set:
            self.jacobian[k, 1, 1] = self.pr2.val
            self.jacobian[k, 3, 1] = -1
            k += 1

        ######################################################################
        # derivatives for specified zeta at hot side
        if self.zeta1.is_set:
            f = self.zeta_func
            if not increment_filter[0, 0]:
                self.jacobian[k, 0, 0] = self.numeric_deriv(
                    f, 'm', 0, zeta='zeta1', inconn=0, outconn=0)
            if not increment_filter[0, 1]:
                self.jacobian[k, 0, 1] = self.numeric_deriv(
                    f, 'p', 0, zeta='zeta1', inconn=0, outconn=0)
            if not increment_filter[0, 2]:
                self.jacobian[k, 0, 2] = self.numeric_deriv(
                    f, 'h', 0, zeta='zeta1', inconn=0, outconn=0)
            if not increment_filter[2, 1]:
                self.jacobian[k, 2, 1] = self.numeric_deriv(
                    f, 'p', 2, zeta='zeta1', inconn=0, outconn=0)
            if not increment_filter[2, 2]:
                self.jacobian[k, 2, 2] = self.numeric_deriv(
                    f, 'h', 2, zeta='zeta1', inconn=0, outconn=0)
            k += 1

        ######################################################################
        # derivatives for specified zeta at cold side
        if self.zeta2.is_set:
            f = self.zeta_func
            if not increment_filter[1, 0]:
                self.jacobian[k, 1, 0] = self.numeric_deriv(
                    f, 'm', 1, zeta='zeta2', inconn=1, outconn=1)
            if not increment_filter[1, 1]:
                self.jacobian[k, 1, 1] = self.numeric_deriv(
                    f, 'p', 1, zeta='zeta2', inconn=1, outconn=1)
            if not increment_filter[1, 2]:
                self.jacobian[k, 1, 2] = self.numeric_deriv(
                    f, 'h', 1, zeta='zeta2', inconn=1, outconn=1)
            if not increment_filter[3, 1]:
                self.jacobian[k, 3, 1] = self.numeric_deriv(
                    f, 'p', 3, zeta='zeta2', inconn=1, outconn=1)
            if not increment_filter[3, 2]:
                self.jacobian[k, 3, 2] = self.numeric_deriv(
                    f, 'h', 3, zeta='zeta2', inconn=1, outconn=1)
            k += 1

        ######################################################################
        # derivatives for additional equations
        self.additional_derivatives(increment_filter, k)

    def additional_derivatives(self, increment_filter, k):
        r"""Calculate partial derivatives for given additional equations."""
        return

    def mass_flow_func(self):
        r"""
        Calculate the residual value for mass flow balance equation.

        Returns
        -------
        residual : list
            Vector with residual value for component's mass flow balance.

            .. math::

                0 = \dot{m}_{in,i} - \dot{m}_{out,i} \;
                \forall i \in inlets/outlets
        """
        residual = []
        for i in range(self.num_i):
            residual += [self.inl[i].m.val_SI - self.outl[i].m.val_SI]
        return residual

    def mass_flow_deriv(self):
        r"""
        Calculate partial derivatives for all mass flow balance equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the mass flow balance
            equations.
        """
        deriv = np.zeros((2, 4 + self.num_vars, self.num_nw_vars))
        for i in range(self.num_i):
            deriv[i, i, 0] = 1
        for j in range(self.num_o):
            deriv[j, j + i + 1, 0] = -1
        return deriv

    def energy_func(self):
        r"""
        Equation for heat exchanger energy balance.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_{1,in} \cdot \left(h_{1,out} - h_{1,in} \right) +
                \dot{m}_{2,in} \cdot \left(h_{2,out} - h_{2,in} \right)
        """
        return (
            self.inl[0].m.val_SI * (
                self.outl[0].h.val_SI - self.inl[0].h.val_SI) +
            self.inl[1].m.val_SI * (
                self.outl[1].h.val_SI - self.inl[1].h.val_SI))

    def kA_func(self):
        r"""
        Calculate heat transfer from heat transfer coefficient.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                res = \dot{m}_{1,in} \cdot \left( h_{1,out} - h_{1,in}\right) +
                kA \cdot \frac{T_{1,out} -
                T_{2,in} - T_{1,in} + T_{2,out}}
                {\ln{\frac{T_{1,out} - T_{2,in}}{T_{1,in} - T_{2,out}}}}

        Note
        ----
        For standard functions f\ :subscript:`1` \ and f\ :subscript:`2` \ see
        module :func:`tespy.data`.

        - Calculate temperatures at inlets and outlets.
        - Perform value manipulation, if temperature levels are not physically
          feasible.
        """
        i1 = self.inl[0].to_flow()
        i2 = self.inl[1].to_flow()
        o1 = self.outl[0].to_flow()
        o2 = self.outl[1].to_flow()

        T_i1 = T_mix_ph(i1, T0=self.inl[0].T.val_SI)
        T_i2 = T_mix_ph(i2, T0=self.inl[1].T.val_SI)
        T_o1 = T_mix_ph(o1, T0=self.outl[0].T.val_SI)
        T_o2 = T_mix_ph(o2, T0=self.outl[1].T.val_SI)

        if T_i1 <= T_o2:
            T_i1 = T_o2 + 0.01
        if T_i1 <= T_o2:
            T_o2 = T_i1 - 0.01
        if T_i1 <= T_o2:
            T_o1 = T_i2 + 0.02
        if T_o1 <= T_i2:
            T_i2 = T_o1 - 0.02

        td_log = ((T_o1 - T_i2 - T_i1 + T_o2) /
                  np.log((T_o1 - T_i2) / (T_i1 - T_o2)))

        return i1[0] * (o1[2] - i1[2]) + self.kA.val * td_log

    def kA_char_func(self):
        r"""
        Calculate heat transfer from heat transfer coefficient characteristic.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                res = \dot{m}_{1,in} \cdot \left( h_{1,out} - h_{1,in}\right) +
                kA_{ref} \cdot f_{kA} \cdot \frac{T_{1,out} -
                T_{2,in} - T_{1,in} + T_{2,out}}
                {\ln{\frac{T_{1,out} - T_{2,in}}{T_{1,in} - T_{2,out}}}}

                f_{kA} = \frac{2}{
                \frac{1}{f_1\left(\frac{m_1}{m_{1,ref}}\right)} +
                \frac{1}{f_2\left(\frac{m_2}{m_{2,ref}}\right)}}

        Note
        ----
        For standard functions f\ :subscript:`1` \ and f\ :subscript:`2` \ see
        module :func:`tespy.data`.

        - Calculate temperatures at inlets and outlets.
        - Perform value manipulation, if temperature levels are not physically
          feasible.
        """
        i1 = self.inl[0].to_flow()
        i2 = self.inl[1].to_flow()
        o1 = self.outl[0].to_flow()
        o2 = self.outl[1].to_flow()

        i1_d = self.inl[0].to_flow_design()
        i2_d = self.inl[1].to_flow_design()

        T_i1 = T_mix_ph(i1, T0=self.inl[0].T.val_SI)
        T_i2 = T_mix_ph(i2, T0=self.inl[1].T.val_SI)
        T_o1 = T_mix_ph(o1, T0=self.outl[0].T.val_SI)
        T_o2 = T_mix_ph(o2, T0=self.outl[1].T.val_SI)

        if T_i1 <= T_o2:
            T_i1 = T_o2 + 0.01
        if T_i1 <= T_o2:
            T_o2 = T_i1 - 0.01
        if T_i1 <= T_o2:
            T_o1 = T_i2 + 0.02
        if T_o1 <= T_i2:
            T_i2 = T_o1 - 0.02

        fkA1 = 1
        if self.kA_char1.param == 'm':
            fkA1 = self.kA_char1.func.evaluate(i1[0] / i1_d[0])

        fkA2 = 1
        if self.kA_char2.param == 'm':
            fkA2 = self.kA_char2.func.evaluate(i2[0] / i2_d[0])

        fkA = 2 / (1 / fkA1 + 1 / fkA2)

        td_log = ((T_o1 - T_i2 - T_i1 + T_o2) /
                  np.log((T_o1 - T_i2) / (T_i1 - T_o2)))
        return i1[0] * (o1[2] - i1[2]) + self.kA.design * fkA * td_log

    def ttd_u_func(self):
        r"""
        Equation for upper terminal temperature difference.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                res = ttd_{u} - T_{1,in} + T_{2,out}
        """
        T_i1 = T_mix_ph(self.inl[0].to_flow(), T0=self.inl[0].T.val_SI)
        T_o2 = T_mix_ph(self.outl[1].to_flow(), T0=self.outl[1].T.val_SI)
        return self.ttd_u.val - T_i1 + T_o2

    def ttd_l_func(self):
        r"""
        Equation for upper terminal temperature difference.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                res = ttd_{l} - T_{1,out} + T_{2,in}
        """
        i2 = self.inl[1].to_flow()
        o1 = self.outl[0].to_flow()
        return (self.ttd_l.val - T_mix_ph(o1, T0=self.outl[0].T.val_SI) +
                T_mix_ph(i2, T0=self.inl[1].T.val_SI))

    def bus_func(self, bus):
        r"""
        Calculate the value of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus
            TESPy bus object.

        Returns
        -------
        val : float
            Value of energy transfer :math:`\dot{E}`. This value is passed to
            :py:meth:`tespy.components.components.component.calc_bus_value`
            for value manipulation according to the specified characteristic
            line of the bus.

            .. math::

                \dot{E} = \dot{m}_{1,in} \cdot \left(
                h_{1,out} - h_{1,in} \right)
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
        bus : tespy.connections.bus
            TESPy bus object.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 4, self.num_nw_vars))
        f = self.calc_bus_value
        deriv[0, 0, 0] = self.numeric_deriv(f, 'm', 0, bus=bus)
        deriv[0, 0, 2] = self.numeric_deriv(f, 'h', 0, bus=bus)
        deriv[0, 2, 2] = self.numeric_deriv(f, 'h', 2, bus=bus)
        return deriv

    def initialise_source(self, c, key):
        r"""
        Return a starting value for pressure and enthalpy at outlet.

        Parameters
        ----------
        c : tespy.connections.connection
            Connection to perform initialisation on.

        key : str
            Fluid property to retrieve.

        Returns
        -------
        val : float
            Starting value for pressure/enthalpy in SI units.

            .. math::

                val = \begin{cases}
                4 \cdot 10^5 & \text{key = 'p'}\\
                h\left(p, 200 \text{K} \right) & \text{key = 'h' at outlet 1}\\
                h\left(p, 250 \text{K} \right) & \text{key = 'h' at outlet 2}
                \end{cases}
        """
        if key == 'p':
            return 50e5
        elif key == 'h':
            flow = c.to_flow()
            if c.source_id == 'out1':
                T = 200 + 273.15
                return h_mix_pT(flow, T)
            else:
                T = 250 + 273.15
                return h_mix_pT(flow, T)

    def initialise_target(self, c, key):
        r"""
        Return a starting value for pressure and enthalpy at inlet.

        Parameters
        ----------
        c : tespy.connections.connection
            Connection to perform initialisation on.

        key : str
            Fluid property to retrieve.

        Returns
        -------
        val : float
            Starting value for pressure/enthalpy in SI units.

            .. math::

                val = \begin{cases}
                4 \cdot 10^5 & \text{key = 'p'}\\
                h\left(p, 300 \text{K} \right) & \text{key = 'h' at inlet 1}\\
                h\left(p, 220 \text{K} \right) & \text{key = 'h' at outlet 2}
                \end{cases}
        """
        if key == 'p':
            return 50e5
        elif key == 'h':
            flow = c.to_flow()
            if c.target_id == 'in1':
                T = 300 + 273.15
                return h_mix_pT(flow, T)
            else:
                T = 220 + 273.15
                return h_mix_pT(flow, T)

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        # connection information
        i1 = self.inl[0].to_flow()
        i2 = self.inl[1].to_flow()
        o1 = self.outl[0].to_flow()
        o2 = self.outl[1].to_flow()

        # temperatures
        if isinstance(self, condenser):
            T_i1 = T_bp_p(i1)
        else:
            T_i1 = T_mix_ph(i1, T0=self.inl[0].T.val_SI)
        T_i2 = T_mix_ph(i2, T0=self.inl[1].T.val_SI)
        T_o1 = T_mix_ph(o1, T0=self.outl[0].T.val_SI)
        T_o2 = T_mix_ph(o2, T0=self.outl[1].T.val_SI)

        # specific volume
        v_i1 = v_mix_ph(i1, T0=T_i1)
        v_i2 = v_mix_ph(i2, T0=T_i2)
        v_o1 = v_mix_ph(o1, T0=T_o1)
        v_o2 = v_mix_ph(o2, T0=T_o2)

        # specific entropy
        s_i1 = s_mix_ph(i1, T0=T_i1)
        s_i2 = s_mix_ph(i2, T0=T_i2)
        s_o1 = s_mix_ph(o1, T0=T_o1)
        s_o2 = s_mix_ph(o2, T0=T_o2)

        # component parameters
        self.ttd_u.val = T_i1 - T_o2
        self.ttd_l.val = T_o1 - T_i2
        self.Q.val = i1[0] * (o1[2] - i1[2])

        self.pr1.val = o1[1] / i1[1]
        self.pr2.val = o2[1] / i2[1]
        self.zeta1.val = ((i1[1] - o1[1]) * np.pi ** 2 /
                          (8 * i1[0] ** 2 * (v_i1 + v_o1) / 2))
        self.zeta2.val = ((i2[1] - o2[1]) * np.pi ** 2 /
                          (8 * i2[0] ** 2 * (v_i2 + v_o2) / 2))

        self.SQ1.val = self.inl[0].m.val_SI * (s_o1 - s_i1)
        self.SQ2.val = self.inl[1].m.val_SI * (s_o2 - s_i2)
        self.Sirr.val = self.SQ1.val + self.SQ2.val

        # kA and logarithmic temperature difference
        if T_i1 <= T_o2 or T_o1 <= T_i2:
            self.td_log.val = np.nan
            self.kA.val = np.nan
        else:
            self.td_log.val = ((T_o1 - T_i2 - T_i1 + T_o2) /
                               np.log((T_o1 - T_i2) / (T_i1 - T_o2)))
            self.kA.val = -(i1[0] * (o1[2] - i1[2]) / self.td_log.val)

        if self.kA_char.is_set:
            # get bound errors for kA hot side characteristics
            if self.kA_char1.param == 'm':
                i1_d = self.inl[0].to_flow_design()
                if not np.isnan(i1_d[0]):
                    if not i1[0] == 0:
                        self.kA_char1.func.get_bound_errors(i1[0] / i1_d[0],
                                                            self.label)

            # get bound errors for kA copld side characteristics
            if self.kA_char2.param == 'm':
                i2_d = self.inl[1].to_flow_design()
                if not np.isnan(i2_d[0]):
                    if not i1[0] == 0:
                        self.kA_char2.func.get_bound_errors(i2[0] / i2_d[0],
                                                            self.label)

        self.check_parameter_bounds()

# %%


class condenser(heat_exchanger):
    r"""
    A condenser cools a fluid until it is in liquid state.

    The condensing fluid is cooled by the cold side fluid. The fluid on the hot
    side of the condenser must be pure. Subcooling is available.

    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_func`
        - :func:`tespy.components.heat_exchangers.heat_exchanger.mass_flow_func`
        - :func:`tespy.components.heat_exchangers.condenser.energy_func`

        **optional equations**

        .. math::

            0 = \dot{m}_{in} \cdot \left(h_{out} - h_{in} \right) - \dot{Q}

        - :func:`tespy.components.heat_exchangers.condenser.kA_func`
        - :func:`tespy.components.heat_exchangers.condenser.kA_char_func`
        - :func:`tespy.components.heat_exchangers.condenser.ttd_u_func`

        .. math::

            0 = p_{1,in} \cdot pr1 - p_{1,out}\\
            0 = p_{2,in} \cdot pr2 - p_{2,out}

        - hot side :func:`tespy.components.components.component.zeta_func`
        - cold side :func:`tespy.components.components.component.zeta_func`

        **additional equations**

        - :func:`tespy.components.heat_exchangers.condenser.additional_equations`

    Inlets/Outlets

        - in1, in2 (index 1: hot side, index 2: cold side)
        - out1, out2 (index 1: hot side, index 2: cold side)

    Image

        .. image:: _images/condenser.svg
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

    pr1 : str/float/tespy.tools.data_containers.dc_cp
        Outlet to inlet pressure ratio at hot side, :math:`pr/1`.

    pr2 : str/float/tespy.tools.data_containers.dc_cp
        Outlet to inlet pressure ratio at cold side, :math:`pr/1`.

    zeta1 : str/float/tespy.tools.data_containers.dc_cp
        Geometry independent friction coefficient at hot side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    zeta2 : str/float/tespy.tools.data_containers.dc_cp
        Geometry independent friction coefficient at cold side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    kA : float/tespy.tools.data_containers.dc_cp
        Area independent heat transition coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char : tespy.tools.data_containers.dc_simple
        Area independent heat transition coefficient characteristic.

    kA_char1 : tespy.tools.charactersitics.char_line/tespy.tools.data_containers.dc_cc
        Characteristic line for hot side heat transfer coefficient.

    kA_char2 : tespy.tools.charactersitics.char_line/tespy.tools.data_containers.dc_cc
        Characteristic line for cold side heat transfer coefficient.

    subcooling : bool
        Enable/disable subcooling, default value: disabled.

    Note
    ----
    The condenser has an additional equation for enthalpy at hot side outlet:
    The fluid leaves the component in saturated liquid state. If subcooling
    is activated, it possible to specify the enthalpy at the outgoing
    connection manually.

    It has different calculation method for given heat transfer coefficient and
    upper terminal temperature dierence: These parameters refer to the
    **condensing** temperature, even if the fluid on the hot side enters the
    component in superheated state.

    Example
    -------
    Air is used to condensate water in a condenser. 1 kg/s waste steam is
    chilled with a terminal temperature difference of 15 K.

    >>> from tespy.components import sink, source, condenser
    >>> from tespy.connections import connection
    >>> from tespy.networks import network
    >>> from tespy.tools.fluid_properties import T_bp_p
    >>> import shutil
    >>> nw = network(fluids=['water', 'air'], T_unit='C', p_unit='bar',
    ... h_unit='kJ / kg', m_range=[0.01, 1000], iterinfo=False)
    >>> amb_in = sink('ambient air inlet')
    >>> amb_out = source('air outlet')
    >>> waste_steam = source('waste steam')
    >>> c = sink('condensate sink')
    >>> cond = condenser('condenser')
    >>> cond.component()
    'condenser'
    >>> amb_he = connection(amb_out, 'out1', cond, 'in2')
    >>> he_amb = connection(cond, 'out2', amb_in, 'in1')
    >>> ws_he = connection(waste_steam, 'out1', cond, 'in1')
    >>> he_c = connection(cond, 'out1', c, 'in1')
    >>> nw.add_conns(amb_he, he_amb, ws_he, he_c)

    The air flow can not be controlled, thus is constant in offdesign
    operation. If the waste steam mass flow or the ambient air temperature
    change, the outlet temperature of the air will change, too.

    >>> cond.set_attr(pr1=0.98, pr2=0.999, ttd_u=15, design=['pr2', 'ttd_u'],
    ... offdesign=['zeta2', 'kA_char'])
    >>> ws_he.set_attr(fluid={'water': 1, 'air': 0}, h=2700, m=1)
    >>> amb_he.set_attr(fluid={'water': 0, 'air': 1}, T=20, offdesign=['v'])
    >>> he_amb.set_attr(p=1, T=40, design=['T'])
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(amb_he.v.val, 2)
    103.17
    >>> round(ws_he.T.val - he_amb.T.val, 1)
    66.9
    >>> round(T_bp_p(ws_he.to_flow()) - 273.15 - he_amb.T.val, 1)
    15.0
    >>> ws_he.set_attr(m=0.7)
    >>> amb_he.set_attr(T=30)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(ws_he.T.val - he_amb.T.val, 1)
    62.5
    >>> round(T_bp_p(ws_he.to_flow()) - 273.15 - he_amb.T.val, 1)
    11.3

    It is possible to activate subcooling. The difference to boiling point
    temperature is specified to 5 K.

    >>> cond.set_attr(subcooling=True)
    >>> he_c.set_attr(Td_bp=-5)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(ws_he.T.val - he_amb.T.val, 1)
    62.5
    >>> round(T_bp_p(ws_he.to_flow()) - 273.15 - he_amb.T.val, 1)
    13.4
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    @staticmethod
    def component():
        return 'condenser'

    @staticmethod
    def attr():
        return {
            'Q': dc_cp(max_val=0), 'kA': dc_cp(min_val=0),
            'td_log': dc_cp(min_val=0),
            'ttd_u': dc_cp(min_val=0), 'ttd_l': dc_cp(min_val=0),
            'pr1': dc_cp(max_val=1), 'pr2': dc_cp(max_val=1),
            'zeta1': dc_cp(min_val=0), 'zeta2': dc_cp(min_val=0),
            'subcooling': dc_simple(val=False),
            'kA_char': dc_simple(),
            'kA_char1': dc_cc(param='m'), 'kA_char2': dc_cc(param='m'),
            'SQ1': dc_simple(), 'SQ2': dc_simple(), 'Sirr': dc_simple(),
        }

    def comp_init(self, nw):

        component.comp_init(self, nw)

        # number of mandatroy equations for
        # fluid balance: num_fl * 2
        # mass flow: 2
        # energy balance: 1
        self.num_eq = self.num_nw_fluids * 2 + 3
        # enthalpy hot side outlet (if not subcooling): 1
        if self.subcooling.val is False:
            self.num_eq += 1
        for var in [self.Q, self.kA, self.kA_char, self.ttd_u, self.ttd_l,
                    self.pr1, self.pr2, self.zeta1, self.zeta2]:
            if var.is_set is True:
                self.num_eq += 1

        if self.kA.is_set:
            msg = (
                'The usage of the parameter kA has changed for offdesign '
                'calculation. Specifying kA will keep a constant value for kA '
                'in the calculation. If you want to use the value adaption of '
                'kA by the characteristic line, please use kA_char as '
                'parameter instead (occurred at ' + self.label + '). This '
                'warning will disappear in TESPy version 0.4.0.')
            warnings.warn(msg, FutureWarning, stacklevel=2)

        self.jacobian = np.zeros((
            self.num_eq,
            self.num_i + self.num_o + self.num_vars,
            self.num_nw_vars))

        self.residual = np.zeros(self.num_eq)
        pos = self.num_nw_fluids * 2
        self.jacobian[0:pos] = self.fluid_deriv()
        self.jacobian[pos:pos + 2] = self.mass_flow_deriv()

    def additional_equations(self, k):
        r"""
        Calculate results of additional equations.

        Equations

            **mandatory equations**

            .. math::

                0 = h_{1,out} - h\left(p, x=0 \right)\\
                x: \text{vapour mass fraction}
        """
        ######################################################################
        # equation for saturated liquid at hot side outlet
        if self.subcooling.val is False:
            o1 = self.outl[0].to_flow()
            self.residual[k] = o1[2] - h_mix_pQ(o1, 0)
            k += 1

    def additional_derivatives(self, increment_filter, k):
        r"""Calculate partial derivatives for given additional equations."""
        ######################################################################
        # derivatives for saturated liquid at hot side outlet equation
        if self.subcooling.val is False:
            o1 = self.outl[0].to_flow()
            self.jacobian[k, 2, 1] = -dh_mix_dpQ(o1, 0)
            self.jacobian[k, 2, 2] = 1
            k += 1

    def kA_func(self):
        r"""
        Calculate heat transfer from heat transfer coefficient.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                res = \dot{m}_{1,in} \cdot \left( h_{1,out} - h_{1,in}\right) +
                kA \cdot \frac{T_{1,out} -
                T_{2,in} - T_s \left(p_{1,in}\right) +
                T_{2,out}}
                {\ln{\frac{T_{1,out} - T_{2,in}}
                {T_s \left(p_{1,in}\right) - T_{2,out}}}}
        """

        i1 = self.inl[0].to_flow()
        i2 = self.inl[1].to_flow()
        o1 = self.outl[0].to_flow()
        o2 = self.outl[1].to_flow()

        T_i1 = T_bp_p(i1)
        T_i2 = T_mix_ph(i2, T0=self.inl[1].T.val_SI)
        T_o1 = T_mix_ph(o1, T0=self.outl[0].T.val_SI)
        T_o2 = T_mix_ph(o2, T0=self.outl[1].T.val_SI)

        if T_i1 <= T_o2 and not self.inl[0].T.val_set:
            T_i1 = T_o2 + 0.5
        if T_i1 <= T_o2 and not self.outl[1].T.val_set:
            T_o2 = T_i1 - 0.5

        if T_o1 <= T_i2 and not self.outl[0].T.val_set:
            T_o1 = T_i2 + 1
        if T_o1 <= T_i2 and not self.inl[1].T.val_set:
            T_i2 = T_o1 - 1

        td_log = ((T_o1 - T_i2 - T_i1 + T_o2) /
                  np.log((T_o1 - T_i2) / (T_i1 - T_o2)))

        return i1[0] * (o1[2] - i1[2]) + self.kA.val * td_log

    def kA_char_func(self):
        r"""
        Calculate heat transfer from heat transfer coefficient characteristic.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                res = \dot{m}_{1,in} \cdot \left( h_{1,out} - h_{1,in}\right) +
                kA_{ref} \cdot f_{kA} \cdot \frac{T_{1,out} -
                T_{2,in} - T_s \left(p_{1,in}\right) +
                T_{2,out}}
                {\ln{\frac{T_{1,out} - T_{2,in}}
                {T_s \left(p_{1,in}\right) - T_{2,out}}}}

                f_{kA} = \frac{2}{
                \frac{1}{f_1\left(\frac{m_1}{m_{1,ref}}\right)} +
                \frac{1}{f_2\left(\frac{m_2}{m_{2,ref}}\right)}}

        Note
        ----
        For standard functions f\ :subscript:`1` \ and f\ :subscript:`2` \ see
        module :func:`tespy.data`.

        - Calculate temperatures at inlets and outlets.
        - Perform value manipulation, if temperature levels are physically
          infeasible.
        """

        i1 = self.inl[0].to_flow()
        i2 = self.inl[1].to_flow()
        o1 = self.outl[0].to_flow()
        o2 = self.outl[1].to_flow()

        i1_d = self.inl[0].to_flow_design()
        i2_d = self.inl[1].to_flow_design()

        T_i1 = T_bp_p(i1)
        T_i2 = T_mix_ph(i2, T0=self.inl[1].T.val_SI)
        T_o1 = T_mix_ph(o1, T0=self.outl[0].T.val_SI)
        T_o2 = T_mix_ph(o2, T0=self.outl[1].T.val_SI)

        if T_i1 <= T_o2 and not self.inl[0].T.val_set:
            T_i1 = T_o2 + 0.5
        if T_i1 <= T_o2 and not self.outl[1].T.val_set:
            T_o2 = T_i1 - 0.5

        if T_o1 <= T_i2 and not self.outl[0].T.val_set:
            T_o1 = T_i2 + 1
        if T_o1 <= T_i2 and not self.inl[1].T.val_set:
            T_i2 = T_o1 - 1

        fkA1 = 1
        if self.kA_char1.param == 'm':
            fkA1 = self.kA_char1.func.evaluate(i1[0] / i1_d[0])

        fkA2 = 1
        if self.kA_char2.param == 'm':
            fkA2 = self.kA_char2.func.evaluate(i2[0] / i2_d[0])

        fkA = 2 / (1 / fkA1 + 1 / fkA2)

        td_log = ((T_o1 - T_i2 - T_i1 + T_o2) /
                  np.log((T_o1 - T_i2) / (T_i1 - T_o2)))

        return i1[0] * (o1[2] - i1[2]) + self.kA.design * fkA * td_log

    def ttd_u_func(self):
        r"""
        Equation for upper terminal temperature difference.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                res = ttd_{u} - T_s \left(p_{1,in}\right) + T_{2,out}

        Note
        ----
        The upper terminal temperature difference ttd_u refers to boiling
        temperature at hot side inlet.
        """
        i1 = self.inl[0].to_flow()
        o2 = self.outl[1].to_flow()
        T_o2 = T_mix_ph(o2, T0=self.outl[1].T.val_SI)
        return self.ttd_u.val - T_bp_p(i1) + T_o2

# %%


class desuperheater(heat_exchanger):
    r"""
    The desuperheater cools a fluid to the saturated gas state.

    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_func`
        - :func:`tespy.components.heat_exchangers.heat_exchanger.mass_flow_func`
        - :func:`tespy.components.heat_exchangers.heat_exchanger.energy_func`

        **optional equations**

        .. math::

            0 = \dot{m}_{in} \cdot \left(h_{out} - h_{in} \right) - \dot{Q}

        - :func:`tespy.components.heat_exchangers.heat_exchanger.kA_func`
        - :func:`tespy.components.heat_exchangers.heat_exchanger.ttd_u_func`

        .. math::

            0 = p_{1,in} \cdot pr1 - p_{1,out}\\
            0 = p_{2,in} \cdot pr2 - p_{2,out}

        - hot side :func:`tespy.components.components.component.zeta_func`
        - cold side :func:`tespy.components.components.component.zeta_func`

        **additional equations**

        - :func:`tespy.components.heat_exchangers.desuperheater.additional_equations`

    Inlets/Outlets

        - in1, in2 (index 1: hot side, index 2: cold side)
        - out1, out2 (index 1: hot side, index 2: cold side)

    Image

        .. image:: _images/heat_exchanger.svg
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

    pr1 : str/float/tespy.tools.data_containers.dc_cp
        Outlet to inlet pressure ratio at hot side, :math:`pr/1`.

    pr2 : str/float/tespy.tools.data_containers.dc_cp
        Outlet to inlet pressure ratio at cold side, :math:`pr/1`.

    zeta1 : str/float/tespy.tools.data_containers.dc_cp
        Geometry independent friction coefficient at hot side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    zeta2 : str/float/tespy.tools.data_containers.dc_cp
        Geometry independent friction coefficient at cold side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    kA : str/float/tespy.tools.data_containers.dc_cp
        Area independent heat transition coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char1 : tespy.tools.charactersitics.char_line/tespy.tools.data_containers.dc_cc
        Characteristic line for hot side heat transfer coefficient.

    kA_char2 : tespy.tools.charactersitics.char_line/tespy.tools.data_containers.dc_cc
        Characteristic line for cold side heat transfer coefficient.

    Note
    ----
    The desuperheater has an additional equation for enthalpy at hot side
    outlet: The fluid leaves the component in saturated gas state.

    Example
    -------
    Overheated enthanol is cooled with water in a heat exchanger until it
    reaches the state of saturated gas.

    >>> from tespy.components import sink, source, desuperheater
    >>> from tespy.connections import connection
    >>> from tespy.networks import network
    >>> from tespy.tools.fluid_properties import T_bp_p
    >>> import shutil
    >>> nw = network(fluids=['water', 'ethanol'], T_unit='C', p_unit='bar',
    ... h_unit='kJ / kg', v_unit='l / s', m_range=[0.001, 10], iterinfo=False)
    >>> et_in = source('ethanol inlet')
    >>> et_out = sink('ethanol outlet')
    >>> cw_in = source('cooling water inlet')
    >>> cw_out = sink('cooling water outlet')
    >>> desu = desuperheater('desuperheater')
    >>> desu.component()
    'desuperheater'
    >>> et_de = connection(et_in, 'out1', desu, 'in1')
    >>> de_et = connection(desu, 'out1', et_out, 'in1')
    >>> cw_de = connection(cw_in, 'out1', desu, 'in2')
    >>> de_cw = connection(desu, 'out2', cw_out, 'in1')
    >>> nw.add_conns(et_de, de_et, cw_de, de_cw)

    The cooling water enters the component at 15 °C. 10 l/s of ethanol is
    cooled from 100 K above boiling point. The water flow rate is at 1 l/s.
    Knowing the component's design parameters it is possible to predict
    behavior at different inlet temperatures or different volumetric flow of
    ethanol. Controlling the ethanol's state at the outlet is only possible,
    if the cooling water flow rate is adjusted accordingly.

    >>> desu.set_attr(pr1=0.99, pr2=0.98, design=['pr1', 'pr2'],
    ... offdesign=['zeta1', 'zeta2', 'kA_char'])
    >>> cw_de.set_attr(fluid={'water': 1, 'ethanol': 0}, T=15, v=1,
    ... design=['v'])
    >>> de_cw.set_attr(p=1)
    >>> et_de.set_attr(fluid={'water': 0, 'ethanol': 1}, Td_bp=100, v=10)
    >>> de_et.set_attr(p=1)
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(de_cw.T.val, 1)
    15.5
    >>> round(de_et.x.val, 1)
    1.0
    >>> et_de.set_attr(v=12)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(cw_de.v.val, 2)
    1.94
    >>> et_de.set_attr(v=7)
    >>> nw.solve('offdesign', design_path='tmp', init_path='tmp')
    >>> round(cw_de.v.val, 2)
    0.41
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    @staticmethod
    def component():
        return 'desuperheater'

    def comp_init(self, nw):

        component.comp_init(self, nw)

        # number of mandatroy equations for
        # fluid balance: num_fl * 2
        # mass flow: 2
        # energy balance: 1
        # enthalpy hot side outlet: 1
        self.num_eq = self.num_nw_fluids * 2 + 4
        for var in [self.Q, self.kA, self.kA_char, self.ttd_u, self.ttd_l,
                    self.pr1, self.pr2, self.zeta1, self.zeta2]:
            if var.is_set is True:
                self.num_eq += 1

        if self.kA.is_set:
            msg = (
                'The usage of the parameter kA has changed for offdesign '
                'calculation. Specifying kA will keep a constant value for kA '
                'in the calculation. If you want to use the value adaption of '
                'kA by the characteristic line, please use kA_char as '
                'parameter instead (occurred at ' + self.label + '). This '
                'warning will disappear in TESPy version 0.4.0.')
            warnings.warn(msg, FutureWarning, stacklevel=2)

        self.jacobian = np.zeros((
            self.num_eq,
            self.num_i + self.num_o + self.num_vars,
            self.num_nw_vars))

        self.residual = np.zeros(self.num_eq)
        pos = self.num_nw_fluids * 2
        self.jacobian[0:pos] = self.fluid_deriv()
        self.jacobian[pos:pos + 2] = self.mass_flow_deriv()

    def additional_equations(self, k):
        r"""
        Calculate results of additional equations.

        Equations

            **mandatory equations**

            .. math::

                0 = h_{1,out} - h\left(p, x=1 \right)\\
                x: \text{vapour mass fraction}
        """
        ######################################################################
        # equation for saturated gas at hot side outlet
        o1 = self.outl[0].to_flow()
        self.residual[k] = o1[2] - h_mix_pQ(o1, 1)

    def additional_derivatives(self, increment_filter, k):
        r"""Calculate partial derivatives for given additional equations."""
        ######################################################################
        # derivatives for saturated gas at hot side outlet equation
        o1 = self.outl[0].to_flow()
        self.jacobian[k, 2, 1] = -dh_mix_dpQ(o1, 1)
        self.jacobian[k, 2, 2] = 1
        k += 1
