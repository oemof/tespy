# -*- coding: utf-8

"""Module of class HeatExchangerSimple.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/heat_exchangers/heat_exchanger_simple.py

SPDX-License-Identifier: MIT
"""

import logging
import warnings

import numpy as np

from tespy.components.component import Component
from tespy.tools.data_containers import dc_cc
from tespy.tools.data_containers import dc_cp
from tespy.tools.data_containers import dc_gcp
from tespy.tools.data_containers import dc_simple
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.fluid_properties import s_mix_ph
from tespy.tools.fluid_properties import v_mix_ph
from tespy.tools.fluid_properties import visc_mix_ph
from tespy.tools.global_vars import err
from tespy.tools.helpers import darcy_friction_factor as dff


class HeatExchangerSimple(Component):
    r"""
    A basic heat exchanger representing a heat source or heat sink.

    The component HeatExchangerSimple is the parent class for the components
    SolarCollector and pipe (of module piping).

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

        - :py:meth:`tespy.components.heat_exchangers.HeatExchangerSimple.additional_equations`

    Inlets/Outlets

        - in1
        - out1

    Image

        .. image:: _images/Pipe.svg
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

    hydro_group : str, tespy.tools.data_containers.dc_gcp
        Parametergroup for pressure drop calculation based on pipes dimensions.
        Choose 'HW' for hazen-williams equation, else darcy friction factor is
        used.

    kA : str, float, tespy.tools.data_containers.dc_cp
        Area independent heat transition coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char : tespy.tools.charactersitics.CharLine, tespy.tools.data_containers.dc_cc
        Characteristic line for heat transfer coefficient.

    Tamb : float, tespy.tools.data_containers.dc_simple
        Ambient temperature, provide parameter in network's temperature unit.

    kA_group : tespy.tools.data_containers.dc_gcp
        Parametergroup for heat transfer calculation from ambient temperature
        and area independent heat transfer coefficient kA.

    Example
    -------
    The HeatExchangerSimple can be used as a sink or source of heat. This
    component does not simulate the secondary side of the heat exchanger. It
    is possible to calculate the pressure ratio with the Darcy-Weisbach
    equation or in case of liquid water use the Hazen-Williams equation.
    Also, given ambient temperature and the heat transfer coeffiecient, it is
    possible to predict heat transfer.

    >>> from tespy.components import Sink, Source, HeatExchangerSimple
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> import shutil
    >>> fluids = ['N2']
    >>> nw = Network(fluids=fluids)
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg', iterinfo=False)
    >>> so1 = Source('source 1')
    >>> si1 = Sink('sink 1')
    >>> heat_sink = HeatExchangerSimple('heat sink')
    >>> heat_sink.component()
    'heat exchanger simple'
    >>> heat_sink.set_attr(Tamb=10, pr=0.95, design=['pr'],
    ... offdesign=['zeta', 'kA_char'])
    >>> inc = Connection(so1, 'out1', heat_sink, 'in1')
    >>> outg = Connection(heat_sink, 'out1', si1, 'in1')
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
            'pr': dc_cp(min_val=1e-4, max_val=1),
            'zeta': dc_cp(min_val=0, max_val=1e15),
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

        Component.comp_init(self, nw)

        self.Tamb.val_SI = ((self.Tamb.val + nw.T[nw.T_unit][0]) *
                            nw.T[nw.T_unit][1])

        # parameters for hydro group
        self.hydro_group.set_attr(elements=[self.L, self.ks, self.D])

        is_set = self.L.is_set and self.ks.is_set and self.D.is_set

        if is_set:
            self.hydro_group.set_attr(is_set=True)
            if self.hydro_group.method == 'HW':
                method = 'Hazen-Williams equation'
                self.ks.set_attr(max_val=200)
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
            if (np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0 or
                    self.always_all_equations):
                self.residual[k] = self.zeta_func(zeta='zeta')
            k += 1

        ######################################################################
        # equation for specified hydro-group paremeters
        if self.hydro_group.is_set:
            if (np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0 or
                    self.always_all_equations):
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

            - :py:meth:`tespy.components.heat_exchangers.HeatExchangerSimple.kA_func`
            - :py:meth:`tespy.components.heat_exchangers.HeatExchangerSimple.kA_char_func`
        """
        ######################################################################
        # equation for specified kA_group paremeters
        if self.kA_group.is_set:
            if (np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0 or
                    self.always_all_equations):
                self.residual[k] = self.kA_func()
            k += 1

        ######################################################################
        # equation for specified kA_char_group paremeters
        if self.kA_char_group.is_set:
            if (np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0 or
                    self.always_all_equations):
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
                self.L.val * dff(re, self.ks.val, self.D.val) /
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

        # For numerical stability: If one temperature difference is negative
        # use mean difference to avoid negative logarithm otherwise use lmtd.
        if (ttd_1 / ttd_2) < 0:
            td_log = (ttd_2 + ttd_1) / 2
        elif ttd_1 > ttd_2:
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
        :py:mod:`tespy.data`.
        """
        i, o = self.inl[0].to_flow(), self.outl[0].to_flow()

        # For numerical stability: If one temperature difference is negative
        # use mean difference to avoid negative logarithm otherwise use lmtd.

        ttd_1 = T_mix_ph(i, T0=self.inl[0].T.val_SI) - self.Tamb.val_SI
        ttd_2 = T_mix_ph(o, T0=self.outl[0].T.val_SI) - self.Tamb.val_SI

        if (ttd_1 / ttd_2) < 0:
            td_log = (ttd_2 + ttd_1) / 2
        elif ttd_1 > ttd_2:
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
            :py:meth:`tespy.components.component.Component.calc_bus_value`
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
                self.kA_char.func.get_bound_errors(
                    i[0] / self.inl[0].m.design, self.label)

        self.check_parameter_bounds()

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
                'isoline_value': self.inl[0].p.val,
                'isoline_value_end': self.outl[0].p.val,
                'starting_point_property': 's',
                'starting_point_value': self.inl[0].s.val,
                'ending_point_property': 's',
                'ending_point_value': self.outl[0].s.val
            }
        }
