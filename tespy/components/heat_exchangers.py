# -*- coding: utf-8

"""This module contains components of type heat exchanger:
        heat_exchanger_simple, heat_exchanger, condenser, desuperheater,
        solar_collector


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/heat_exchangers.py

SPDX-License-Identifier: MIT
"""

import logging

import numpy as np

from tespy.components.components import component

from tespy.tools.data_containers import dc_cc, dc_cp, dc_simple, dc_gcp
from tespy.tools.fluid_properties import (
        h_mix_pT, s_mix_ph, v_mix_ph, visc_mix_ph, T_mix_ph,
        dh_mix_dpQ, h_mix_pQ, T_bp_p, memorise
        )
from tespy.tools.helpers import lamb, single_fluid

# %%


class heat_exchanger_simple(component):
    r"""
    The component heat_exchanger_simple is the parent class for pipe and
    solar_collector.

    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_func`
        - :func:`tespy.components.components.component.mass_flow_func`

        **optional equations**

        - :func:`tespy.components.components.heat_exchanger_simple.Q_func`

        .. math::

            0 = p_{in} \cdot pr - p_{out}

        - :func:`tespy.components.components.component.zeta_func`

        - :func:`tespy.components.components.heat_exchanger_simple.darcy_func`
          or :func:`tespy.components.components.heat_exchanger_simple.hw_func`

        **additional equations**

        - :func:`tespy.components.components.heat_exchanger_simple.additional_equations`

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

    Q : Sring/float/tespy.helpers.dc_cp
        Heat transfer, :math:`Q/\text{W}`.

    pr : Sring/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio, :math:`pr/1`.

    zeta : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    D : str/float/tespy.helpers.dc_cp
        Diameter of the pipes, :math:`D/\text{m}`.

    L : str/float/tespy.helpers.dc_cp
        Length of the pipes, :math:`L/\text{m}`.

    ks : str/float/tespy.helpers.dc_cp
        Pipes roughness, :math:`ks/\text{m}` for darcy friction,
        :math:`ks/\text{1}` for hazen-williams equation.

    hydro_group : Sring/tespy.helpers.dc_gcp
        Parametergroup for pressure drop calculation based on pipes dimensions.
        Choose 'HW' for hazen-williams equation, else darcy friction factor is
        used.

    kA : str/float/tespy.helpers.dc_cp
        Area independent heat transition coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient, provide x and y
        values or use generic values (e. g. calculated from design case).
        Standard method 'HE_COLD', Parameter 'm'.

    Tamb : float/tespy.helpers.dc_cp
        Ambient temperature, provide parameter in network's temperature
        unit.

    Tamb_ref : float/tespy.helpers.dc_cp
         Ambient temperature for reference in offdesign case, provide
         parameter in network's temperature unit.

    kA_group : tespy.helpers.dc_gcp
        Parametergroup for heat transfer calculation from ambient temperature
        and area independent heat transfer coefficient kA.

    Example
    -------
    >>> from tespy import cmp, con, nwk
    >>> import shutil
    >>> fluids = ['H2O']
    >>> nw = nwk.network(fluids=fluids)
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')
    >>> nw.set_printoptions(print_level='none')
    >>> so1 = cmp.source('source 1')
    >>> si1 = cmp.sink('sink 1')
    >>> pi = cmp.pipe('pipe')
    >>> pi.component()
    'pipe'
    >>> pi.set_attr(Tamb=10, pr=0.95, design=['pr'], offdesign=['zeta', 'kA'])
    >>> inc = con.connection(so1, 'out1', pi, 'in1')
    >>> outg = con.connection(pi, 'out1', si1, 'in1')
    >>> nw.add_conns(inc, outg)
    >>> inc.set_attr(fluid={'H2O': 1}, m=1, T=200, p=12)
    >>> outg.set_attr(T=190, design=['T'])
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(pi.Q.val, 1)
    -22252.3
    >>> inc.set_attr(m=1.2)
    >>> pi.set_attr(Tamb=-10)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(pi.kA.val, 1)
    126.5
    >>> round(pi.Q.val, 1)
    -25890.6
    >>> round(outg.T.val, 1)
    189.5
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    def component(self):
        return 'heat exchanger simple'

    def attr(self):
        return {'Q': dc_cp(),
                'pr': dc_cp(min_val=1e-4, max_val=1),
                'zeta': dc_cp(min_val=1e-4),
                'D': dc_cp(min_val=1e-2, max_val=2, d=1e-3),
                'L': dc_cp(min_val=1e-1, d=1e-3),
                'ks': dc_cp(val=1e-4, min_val=1e-7, max_val=1e-4, d=1e-8),
                'kA': dc_cp(min_val=0, d=1),
                'Tamb': dc_cp(),
                'kA_char': dc_cc(method='HE_HOT', param='m'),
                'SQ1': dc_simple(), 'SQ2': dc_simple(), 'Sirr': dc_simple(),
                'hydro_group': dc_gcp(), 'kA_group': dc_gcp()}

    def inlets(self):
        return ['in1']

    def outlets(self):
        return ['out1']

    def comp_init(self, nw):

        component.comp_init(self, nw)

        self.fl_deriv = self.fluid_deriv()
        self.m_deriv = self.mass_flow_deriv()

        self.Tamb.val_SI = ((self.Tamb.val + nw.T[nw.T_unit][0]) *
                            nw.T[nw.T_unit][1])
        self.Tamb.design = ((self.Tamb.design + nw.T[nw.T_unit][0]) *
                            nw.T[nw.T_unit][1])

        # parameters for hydro group
        self.hydro_group.set_attr(elements=[self.L, self.ks, self.D])

        is_set = True
        for e in self.hydro_group.elements:
            if not e.is_set:
                is_set = False

        if is_set:
            self.hydro_group.set_attr(is_set=True)
            if self.hydro_group.method == 'HW':
                method = 'Hazen-Williams equation.'
            else:
                method = 'darcy friction factor.'
            msg = ('Pressure loss calculation from pipe dimensions method is '
                   'set to ' + method + '.')
            logging.debug(msg)

        elif self.hydro_group.is_set:
            msg = ('All parameters of the component group have to be '
                   'specified! This component group uses the following '
                   'parameters: L, ks, D at ' + self.label + '. '
                   'Group will be set to False.')
            logging.info(msg)
            self.hydro_group.set_attr(is_set=False)
        else:
            self.hydro_group.set_attr(is_set=False)

        # parameters for kA group
        self.kA_group.set_attr(elements=[self.kA, self.Tamb])

        is_set = True
        for e in self.kA_group.elements:
            if not e.is_set:
                is_set = False

        if is_set:
            self.kA_group.set_attr(is_set=True)
        elif self.kA_group.is_set:
            msg = ('All parameters of the component group have to be '
                   'specified! This component group uses the following '
                   'parameters: kA, Tamb at ' + self.label + '. '
                   'Group will be set to False.')
            logging.info(msg)
            self.kA_group.set_attr(is_set=False)
        else:
            self.kA_group.set_attr(is_set=False)

    def equations(self):
        r"""
        Calculates vector vec_res with results of equations for this component.

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # equations for fluid balance
        vec_res += self.fluid_func()

        ######################################################################
        # equations for mass flow balance
        vec_res += self.mass_flow_func()

        ######################################################################
        # equations for specified heta transfer
        if self.Q.is_set:
            vec_res += [self.Q_func()]

        ######################################################################
        # equations for specified pressure ratio
        if self.pr.is_set:
            vec_res += [self.inl[0].p.val_SI * self.pr.val -
                        self.outl[0].p.val_SI]

        ######################################################################
        # equations for specified zeta
        if self.zeta.is_set:
            vec_res += [self.zeta_func()]

        ######################################################################
        # equation for specified hydro-group paremeters
        if self.hydro_group.is_set:
            # hazen williams equation
            if self.hydro_group.method == 'HW':
                func = self.hw_func
            # darcy friction factor
            else:
                func = self.darcy_func
            vec_res += [func()]

        ######################################################################
        # additional equations
        vec_res += self.additional_equations()

        return vec_res

    def additional_equations(self):
        r"""
        Calculates vector vec_res with results of additional equations for this
        component.

        Equations

            **optional equations**

            - :func:`tespy.components.components.heat_exchanger_simple.kA_func`

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # equation for specified kA-group paremeters
        if self.kA_group.is_set:
            vec_res += [self.kA_func()]

        return vec_res

    def derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given equations.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        mat_deriv = []

        ######################################################################
        # derivatives for fluid balance
        mat_deriv += self.fl_deriv

        ######################################################################
        # derivatives for mass flow balance
        mat_deriv += self.m_deriv

        ######################################################################
        # derivatives for specified heat transfer
        if self.Q.is_set:
            mat_deriv += self.Q_deriv()

        ######################################################################
        # derivatives for specified pressure ratio
        if self.pr.is_set:
            pr_deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))
            pr_deriv[0, 0, 1] = self.pr.val
            pr_deriv[0, 1, 1] = -1
            # custom variable pr
            if self.pr.is_var:
                pr_deriv[0, 2 + self.pr.var_pos, 0] = self.inl[0].p.val_SI
            mat_deriv += pr_deriv.tolist()

        ######################################################################
        # derivatives for specified zeta
        if self.zeta.is_set:
            zeta_deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))
            zeta_deriv[0, 0, 0] = self.numeric_deriv(self.zeta_func, 'm', 0)
            zeta_deriv[0, 0, 1] = self.numeric_deriv(self.zeta_func, 'p', 0)
            zeta_deriv[0, 0, 2] = self.numeric_deriv(self.zeta_func, 'h', 0)
            zeta_deriv[0, 1, 1] = self.numeric_deriv(self.zeta_func, 'p', 1)
            zeta_deriv[0, 1, 2] = self.numeric_deriv(self.zeta_func, 'h', 1)
            # custom variable zeta
            if self.zeta.is_var:
                zeta_deriv[0, 2 + self.zeta.var_pos, 0] = (
                        self.numeric_deriv(self.zeta_func, 'zeta', 2))
            mat_deriv += zeta_deriv.tolist()

        ######################################################################
        # derivatives for specified hydro-group parameters
        if self.hydro_group.is_set:
            # hazen williams equation
            if self.hydro_group.method == 'HW':
                func = self.hw_func
            # darcy friction factor
            else:
                func = self.darcy_func

            deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))
            deriv[0, 0, 0] = self.numeric_deriv(func, 'm', 0)
            deriv[0, 0, 1] = self.numeric_deriv(func, 'p', 0)
            deriv[0, 0, 2] = self.numeric_deriv(func, 'h', 0)
            deriv[0, 1, 1] = self.numeric_deriv(func, 'p', 1)
            deriv[0, 1, 2] = self.numeric_deriv(func, 'h', 1)
            # custom variables of hydro group
            for var in self.hydro_group.elements:
                if var.is_var:
                    deriv[0, 2 + var.var_pos, 0] = (
                            self.numeric_deriv(func, self.vars[var], 2))
            mat_deriv += deriv.tolist()

        ######################################################################
        # derivatives for additional equations
        mat_deriv += self.additional_derivatives()

        return np.asarray(mat_deriv)

    def additional_derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given additional
        equations.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        mat_deriv = []

        ######################################################################
        # derivatives for specified kA-group paremeters
        if self.kA_group.is_set:
            deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))
            deriv[0, 0, 0] = self.numeric_deriv(self.kA_func, 'm', 0)
            deriv[0, 0, 1] = self.numeric_deriv(self.kA_func, 'p', 0)
            deriv[0, 0, 2] = self.numeric_deriv(self.kA_func, 'h', 0)
            deriv[0, 1, 1] = self.numeric_deriv(self.kA_func, 'p', 1)
            deriv[0, 1, 2] = self.numeric_deriv(self.kA_func, 'h', 1)
            #
            for var in self.kA_group.elements:
                if var.is_var:
                    deriv[0, 2 + var.var_pos, 0] = (
                            self.numeric_deriv(self.kA_func, self.vars[var], 2)
                            )
            mat_deriv += deriv.tolist()

        return mat_deriv

    def Q_func(self):
        r"""
        Equation for heat transfer of the simple heat exchanger.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                res = \dot{m}_{in} \cdot \left(h_{out} - h_{in} \right) -
                \dot{Q}
        """
        return self.inl[0].m.val_SI * (
                self.outl[0].h.val_SI - self.inl[0].h.val_SI) - self.Q.val

    def Q_deriv(self):
        r"""
        Calculates the matrix of partial derivatives for heat transfer equation.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))
        deriv[0, 0, 0] = self.outl[0].h.val_SI - self.inl[0].h.val_SI
        deriv[0, 0, 2] = -self.inl[0].m.val_SI
        deriv[0, 1, 2] = self.inl[0].m.val_SI
        # custom variable Q
        if self.Q.is_var:
            deriv[0, 2 + self.Q.var_pos, 0] = -1

        return deriv.tolist()

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
        Equation for heat transfer calculation from ambient conditions and heat
        transfer coefficient.

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
                kA \cdot f_{kA} \cdot \frac{ttd_u - ttd_l}
                {\ln{\frac{ttd_u}{ttd_l}}}

                f_{kA} = f_1\left(\frac{m_1}{m_{1,ref}}\right)

                T_{amb}: \text{ambient temperature}

            for f\ :subscript:`1` \ see class
            :func:`tespy.components.characteristics.characteristics`
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

        fkA = 1
        if not np.isnan(self.inl[0].m.design):
            if self.kA_char.param == 'm':
                fkA = self.kA_char.func.f_x(i[0] / self.inl[0].m.design)

        return i[0] * (o[2] - i[2]) + self.kA.val * fkA * td_log

    def bus_func(self, bus):
        r"""
        Calculates the residual value of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus
            TESPy bus object.

        Returns
        -------
        val : float
            Residual value of equation.

            .. math::

                val = P \cdot f\left( \frac{P}{P_{ref}}\right)

                P = \dot{m}_{in} \cdot \left( h_{out} - h_{in} \right)
        """
        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()

        val = i[0] * (o[2] - i[2])
        if np.isnan(bus.P_ref):
            expr = 1
        else:
            expr = abs(val / bus.P_ref)
        return val * bus.char.f_x(expr)

    def bus_deriv(self, bus):
        r"""
        Calculates the matrix of partial derivatives of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus
            TESPy bus object.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 2, len(self.inl[0].fluid.val) + 3))
        deriv[0, 0, 0] = self.numeric_deriv(self.bus_func, 'm', 0, bus=bus)
        deriv[0, 0, 2] = self.numeric_deriv(self.bus_func, 'h', 0, bus=bus)
        deriv[0, 1, 2] = self.numeric_deriv(self.bus_func, 'h', 1, bus=bus)
        return deriv

    def initialise_source(self, c, key):
        r"""
        Returns a starting value for pressure and enthalpy at component's
        outlet.

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
        Returns a starting value for pressure and enthalpy at component's
        inlet.

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
        r"""
        Postprocessing parameter calculation.
        """
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


class heat_exchanger(component):
    r"""
    Class heat_exchanger is the parent class for condenser and desuperheater.

    Equations

        **mandatory equations**

        - :func:`tespy.components.components.heat_exchanger.fluid_func`
        - :func:`tespy.components.components.heat_exchanger.mass_flow_func`

        **heat exchanger**
        - :func:`tespy.components.components.heat_exchanger.energy_func`

        **optional equations**

        .. math::

            0 = \dot{m}_{in} \cdot \left(h_{out} - h_{in} \right) - \dot{Q}

        **heat exchanger**

        - :func:`tespy.components.components.heat_exchanger.kA_func`
        - :func:`tespy.components.components.heat_exchanger.ttd_u_func`
        - :func:`tespy.components.components.heat_exchanger.ttd_l_func`

        .. math::

            0 = p_{1,in} \cdot pr1 - p_{1,out}\\
            0 = p_{2,in} \cdot pr2 - p_{2,out}

        - :func:`tespy.components.components.heat_exchanger.zeta_func`
        - :func:`tespy.components.components.heat_exchanger.zeta2_func`

        **additional equations**

        - :func:`tespy.components.components.heat_exchanger.additional_equations`

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

    Q : Sring/float/tespy.helpers.dc_cp
        Heat transfer, :math:`Q/\text{W}`.

    pr1 : Sring/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio at hot side, :math:`pr/1`.

    pr2 : Sring/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio at cold side, :math:`pr/1`.

    zeta1 : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient at hot side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    zeta2 : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient at cold side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    kA : str/float/tespy.helpers.dc_cp
        Area independent heat transition coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char1 : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient at hot side, provide
        x and y values or use generic values (e. g. calculated from design
        case). Standard method 'HE_HOT', Parameter 'm'.

    kA_char2 : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient at cold side,
        provide x and y values or use generic values (e. g. calculated from
        design case). Standard method 'HE_COLD', Parameter 'm'.

    Note
    ---
    The heat exchanger and subclasses (desuperheater, condenser) are
    countercurrent heat exchangers. Equations (kA, ttd_u, ttd_l) do not work
    for directcurrent and crosscurrent or combinations of different types.

    Example
    -------
    >>> from tespy import cmp, con, nwk
    >>> import shutil
    >>> nw = nwk.network(fluids=['water'], T_unit='C', p_unit='bar',
    ...     h_unit='kJ / kg')
    >>> nw.set_printoptions(print_level='none')
    >>> tesin = cmp.sink('TES in')
    >>> tesout = cmp.source('TES out')
    >>> hsin = cmp.sink('HS in')
    >>> hsout = cmp.source('HS out')
    >>> he = cmp.heat_exchanger('heat exchanger')
    >>> he.component()
    'heat exchanger'
    >>> tes_he = con.connection(tesout, 'out1', he, 'in2')
    >>> he_tes = con.connection(he, 'out2', tesin, 'in1')
    >>> hs_he = con.connection(hsout, 'out1', he, 'in1')
    >>> he_hs = con.connection(he, 'out1', hsin, 'in1')
    >>> nw.add_conns(tes_he, he_tes, hs_he, he_hs)
    >>> he.set_attr(pr1=0.98, pr2=0.98, ttd_u=5,
    ...     design=['pr1', 'pr2', 'ttd_u'], offdesign=['zeta1', 'zeta2', 'kA'])
    >>> hs_he.set_attr(Td_bp=-10, p=3, fluid={'water': 1})
    >>> he_hs.set_attr(T=70)
    >>> tes_he.set_attr(p=5, fluid={'water': 1})
    >>> tes_he.set_attr(T=40)
    >>> he.set_attr(Q=-80e3)
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(tes_he.m.val, 2)
    0.24
    >>> round(he_tes.T.val, 1)
    118.5
    >>> he.set_attr(Q=-60e3)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(tes_he.m.val, 2)
    0.18
    >>> round(he_tes.T.val, 1)
    119.4
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    def component(self):
        return 'heat exchanger'

    def attr(self):
        return {'Q': dc_cp(max_val=0),
                'kA': dc_cp(min_val=0),
                'td_log': dc_cp(min_val=0),
                'ttd_u': dc_cp(min_val=0), 'ttd_l': dc_cp(min_val=0),
                'pr1': dc_cp(max_val=1), 'pr2': dc_cp(max_val=1),
                'zeta1': dc_cp(min_val=0), 'zeta2': dc_cp(min_val=0),
                'kA_char1': dc_cc(method='HE_HOT', param='m'),
                'kA_char2': dc_cc(method='HE_COLD', param='m'),
                'SQ1': dc_simple(), 'SQ2': dc_simple(), 'Sirr': dc_simple(),
                'zero_flag': dc_simple()}

    def inlets(self):
        return ['in1', 'in2']

    def outlets(self):
        return ['out1', 'out2']

    def comp_init(self, nw):

        component.comp_init(self, nw)

        self.fl_deriv = self.fluid_deriv()
        self.m_deriv = self.mass_flow_deriv()

    def equations(self):
        r"""
        Calculates vector vec_res with results of equations for this component.

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # equations for fluid balance
        vec_res += self.fluid_func()

        ######################################################################
        # equations for mass flow balance
        vec_res += self.mass_flow_func()

        ######################################################################
        # equations for energy balance
        vec_res += [self.energy_func()]

        ######################################################################
        # equations for specified heat transfer
        if self.Q.is_set:
            vec_res += [self.inl[0].m.val_SI *
                        (self.outl[0].h.val_SI - self.inl[0].h.val_SI) -
                        self.Q.val]

        ######################################################################
        # equations for specified heat transfer coefficient
        if self.kA.is_set:
            vec_res += [self.kA_func()]

        ######################################################################
        # equations for specified upper terminal temperature difference
        if self.ttd_u.is_set:
            vec_res += [self.ttd_u_func()]

        ######################################################################
        # equations for specified lower terminal temperature difference
        if self.ttd_l.is_set:
            vec_res += [self.ttd_l_func()]

        ######################################################################
        # equations for specified pressure ratio at hot side
        if self.pr1.is_set:
            vec_res += [self.pr1.val * self.inl[0].p.val_SI -
                        self.outl[0].p.val_SI]

        ######################################################################
        # equations for specified pressure ratio at cold side
        if self.pr2.is_set:
            vec_res += [self.pr2.val * self.inl[1].p.val_SI -
                        self.outl[1].p.val_SI]

        ######################################################################
        # equations for specified zeta at hot side
        if self.zeta1.is_set:
            vec_res += [self.zeta_func()]

        ######################################################################
        # equations for specified zeta at cold side
        if self.zeta2.is_set:
            vec_res += [self.zeta2_func()]

        ######################################################################
        # additional equations
        vec_res += self.additional_equations()

        return vec_res

    def additional_equations(self):
        r"""
        Calculates vector vec_res with results of additional equations for
        this component.

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        return []

    def derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given equations.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        mat_deriv = []

        ######################################################################
        # derivatives for fluid balance equations
        mat_deriv += self.fl_deriv
        ######################################################################
        # derivatives for mass flow balance equations
        mat_deriv += self.m_deriv

        ######################################################################
        # derivatives for energy balance equation
        mat_deriv += self.energy_deriv()

        ######################################################################
        # derivatives for specified heat transfer
        if self.Q.is_set:
            deriv = np.zeros((1, 4, self.num_fl + 3))
            deriv[0, 0, 0] = self.outl[0].h.val_SI - self.inl[0].h.val_SI
            deriv[0, 0, 2] = -self.inl[0].m.val_SI
            deriv[0, 2, 2] = self.inl[0].m.val_SI
            mat_deriv += deriv.tolist()

        ######################################################################
        # derivatives for specified heat transfer coefficient
        if self.kA.is_set:
            kA_deriv = np.zeros((1, 4, self.num_fl + 3))
            kA_deriv[0, 0, 0] = self.numeric_deriv(self.kA_func, 'm', 0)
            kA_deriv[0, 1, 0] = self.numeric_deriv(self.kA_func, 'm', 1)
            for i in range(4):
                kA_deriv[0, i, 1] = self.numeric_deriv(self.kA_func, 'p', i)
                kA_deriv[0, i, 2] = self.numeric_deriv(self.kA_func, 'h', i)
            mat_deriv += kA_deriv.tolist()

        ######################################################################
        # derivatives for specified upper terminal temperature difference
        if self.ttd_u.is_set:
            mat_deriv += self.ttd_u_deriv()

        ######################################################################
        # derivatives for specified lower terminal temperature difference
        if self.ttd_l.is_set:
            mat_deriv += self.ttd_l_deriv()

        ######################################################################
        # derivatives for specified pressure ratio at hot side
        if self.pr1.is_set:
            pr1_deriv = np.zeros((1, 4, self.num_fl + 3))
            pr1_deriv[0, 0, 1] = self.pr1.val
            pr1_deriv[0, 2, 1] = -1
            mat_deriv += pr1_deriv.tolist()

        ######################################################################
        # derivatives for specified pressure ratio at cold side
        if self.pr2.is_set:
            pr2_deriv = np.zeros((1, 4, self.num_fl + 3))
            pr2_deriv[0, 1, 1] = self.pr2.val
            pr2_deriv[0, 3, 1] = -1
            mat_deriv += pr2_deriv.tolist()

        ######################################################################
        # derivatives for specified zeta at hot side
        if self.zeta1.is_set:
            zeta1_deriv = np.zeros((1, 4, self.num_fl + 3))
            zeta1_deriv[0, 0, 0] = self.numeric_deriv(self.zeta_func, 'm', 0)
            zeta1_deriv[0, 0, 1] = self.numeric_deriv(self.zeta_func, 'p', 0)
            zeta1_deriv[0, 0, 2] = self.numeric_deriv(self.zeta_func, 'h', 0)
            zeta1_deriv[0, 2, 1] = self.numeric_deriv(self.zeta_func, 'p', 2)
            zeta1_deriv[0, 2, 2] = self.numeric_deriv(self.zeta_func, 'h', 2)
            mat_deriv += zeta1_deriv.tolist()

        ######################################################################
        # derivatives for specified zeta at cold side
        if self.zeta2.is_set:
            zeta2_deriv = np.zeros((1, 4, self.num_fl + 3))
            zeta2_deriv[0, 1, 0] = self.numeric_deriv(self.zeta2_func, 'm', 1)
            zeta2_deriv[0, 1, 1] = self.numeric_deriv(self.zeta2_func, 'p', 1)
            zeta2_deriv[0, 1, 2] = self.numeric_deriv(self.zeta2_func, 'h', 1)
            zeta2_deriv[0, 3, 1] = self.numeric_deriv(self.zeta2_func, 'p', 3)
            zeta2_deriv[0, 3, 2] = self.numeric_deriv(self.zeta2_func, 'h', 3)
            mat_deriv += zeta2_deriv.tolist()

        ######################################################################
        # derivatives for additional equations
        mat_deriv += self.additional_derivatives()

        return np.asarray(mat_deriv)

    def additional_derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given additional
        equations.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        return []

    def fluid_func(self):
        r"""
        Calculates the vector of residual values for component's fluid balance
        equations.

        Returns
        -------
        vec_res : list
            Vector of residual values for component's fluid balance.

            .. math::

                0 = fluid_{i,in_{j}} - fluid_{i,out_{j}} \;
                \forall i \in \mathrm{fluid}, \; \forall j \in inlets/outlets
        """
        vec_res = []

        for i in range(self.num_i):
            for fluid, x in self.inl[i].fluid.val.items():
                vec_res += [x - self.outl[i].fluid.val[fluid]]
        return vec_res

    def mass_flow_func(self):
        r"""
        Calculates the residual value for component's mass flow balance
        equation.

        Returns
        -------
        vec_res : list
            Vector with residual value for component's mass flow balance.

            .. math::

                0 = \dot{m}_{in,i} - \dot{m}_{out,i} \;
                \forall i \in inlets/outlets
        """
        vec_res = []
        for i in range(self.num_i):
            vec_res += [self.inl[i].m.val_SI - self.outl[i].m.val_SI]
        return vec_res

    def fluid_deriv(self):
        r"""
        Calculates the partial derivatives for all fluid balance equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((self.num_fl * 2, 4 + self.num_vars, 3 + self.num_fl))
        # hot side
        i = 0
        for fluid in self.fluids:
            deriv[i, 0, i + 3] = 1
            deriv[i, 2, i + 3] = -1
            i += 1
        # cold side
        j = 0
        for fluid in self.fluids:
            deriv[i + j, 1, j + 3] = 1
            deriv[i + j, 3, j + 3] = -1
            j += 1
        return deriv.tolist()

    def mass_flow_deriv(self):
        r"""
        Calculates the partial derivatives for all mass flow balance equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the mass flow balance
            equations.
        """
        deriv = np.zeros((2, 4 + self.num_vars, self.num_fl + 3))
        for i in range(self.num_i):
            deriv[i, i, 0] = 1
        for j in range(self.num_o):
            deriv[j, j + i + 1, 0] = -1
        return deriv.tolist()

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
#        if self.zero_flag.val_set:
#            c = self.zero_flag.val
#            if c[0] > 0 and c[1] < 3:
#                return self.inl[0].m.val_SI
#
#            elif ((c[0] == 0 and c[1] < 3) or
#                  (c[0] > 1 and c[1] > 2 and c[1] < 5)):
#                return self.outl[0].h.val_SI - self.inl[0].h.val_SI
#
#            elif ((c[0] < 2 and c[1] > 2 and c[1] < 5) or
#                  (c[0] == 3 and c[1] == 5)):
#                return self.inl[1].m.val_SI
#            else:
#                return self.outl[1].h.val_SI - self.inl[1].h.val_SI
#
#        else:
        return (self.inl[0].m.val_SI * (self.outl[0].h.val_SI -
                                        self.inl[0].h.val_SI) +
                self.inl[1].m.val_SI * (self.outl[1].h.val_SI -
                                        self.inl[1].h.val_SI))

    def energy_deriv(self):
        r"""
        Calculates the matrix of partial derivatives for energy balance
        equation.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 4, len(self.inl[0].fluid.val) + 3))

#        if self.zero_flag.val_set:
#            c = self.zero_flag.val
#            if c[0] > 0 and c[1] < 3:
#                deriv[0, 0, 0] = 1
#
#            elif ((c[0] == 0 and c[1] < 3) or
#                  (c[0] > 1 and c[1] > 2 and c[1] < 5)):
#                deriv[0, 0, 2] = -1
#                deriv[0, 2, 2] = 1
#
#            elif ((c[0] < 2 and c[1] > 2 and c[1] < 5) or
#                  (c[0] == 3 and c[1] == 5)):
#                deriv[0, 1, 0] = 1
#            else:
#                deriv[0, 1, 2] = -1
#                deriv[0, 3, 2] = 1
#
#        else:
        for k in range(2):
            deriv[0, k, 0] = self.outl[k].h.val_SI - self.inl[k].h.val_SI
            deriv[0, k, 2] = -self.inl[k].m.val_SI

        deriv[0, 2, 2] = self.inl[0].m.val_SI
        deriv[0, 3, 2] = self.inl[1].m.val_SI
        return deriv.tolist()

    def kA_func(self):
        r"""
        Equation for heat transfer from conditions on both sides of heat
        exchanger.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                res = \dot{m}_{1,in} \cdot \left( h_{1,out} - h_{1,in}\right) +
                kA \cdot f_{kA} \cdot \frac{T_{1,out} -
                T_{2,in} - T_{1,in} + T_{2,out}}
                {\ln{\frac{T_{1,out} - T_{2,in}}{T_{1,in} - T_{2,out}}}}

                f_{kA} = f_1\left(\frac{m_1}{m_{1,ref}}\right) \cdot
                f_2\left(\frac{m_2}{m_{2,ref}}\right)

        Note
        ----
        For standard functions f\ :subscript:`1` \ and f\ :subscript:`2` \ see
        class :func:`tespy.components.characteristics.characteristics`.

        - Calculate temperatures at inlets and outlets.
        - Perform value manipulation, if temperature levels are not physically
          feasible.
        """

#        if self.zero_flag.val_set:
#            c = self.zero_flag.val
#            if c[1] == 2 or c[1] == 4 or c[1] == 5:
#                T_i1 = T_mix_ph(self.inl[0].to_flow(), T0=self.inl[0].T.val_SI)
#                T_i2 = T_mix_ph(self.inl[1].to_flow(), T0=self.inl[1].T.val_SI)
#                T_o1 = T_mix_ph(self.outl[0].to_flow(),
#                                T0=self.outl[0].T.val_SI)
#                T_o2 = T_mix_ph(self.outl[1].to_flow(),
#                                T0=self.outl[1].T.val_SI)
#                return T_o1 - T_i2 - T_i1 + T_o2
#
#            elif c[0] < 3 and (c[1] == 1 or c[1] == 3):
#                return self.outl[1].h.val_SI - self.inl[1].h.val_SI
#
#            elif ((c[0] < 2 and c[1] == 0) or
#                  (c[0] == 3 and (c[1] == 1 or c[1] == 3))):
#                return self.inl[1].m.val_SI
#
#            else:
#                return self.outl[0].h.val_SI - self.inl[0].h.val_SI

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

#        if T_i1 <= T_o2 and self.inl[0].T.val_set is False:
        if T_i1 <= T_o2:
            T_i1 = T_o2 + 0.01
#        if T_i1 <= T_o2 and self.outl[1].T.val_set is False:
        if T_i1 <= T_o2:
            T_o2 = T_i1 - 0.01
#        if T_i1 < T_o2 and self.inl[0].T.val_set and self.outl[1].T.val_set:
#            msg = ('Infeasibility at ' + str(self.label) + ': Value for upper '
#                   'temperature difference is ' + str(round(T_i1 - T_o2)) +
#                   '.')
#            logging.error(msg)
#            raise ValueError(msg)

#        if T_i1 <= T_o2 and self.outl[1].T.val_set is False:
        if T_i1 <= T_o2:
            T_o1 = T_i2 + 0.02
#        if T_o1 <= T_i2 and self.inl[1].T.val_set is False:
        if T_o1 <= T_i2:
            T_i2 = T_o1 - 0.02
#        if T_o1 < T_i2 and self.inl[1].T.val_set and self.outl[0].T.val_set:
#            msg = ('Infeasibility at ' + str(self.label) + ': Value for lower '
#                   'temperature difference is ' + str(round(T_o1 - T_i2)) +
#                   '.')
#            logging.error(msg)
#            raise ValueError(msg)

        fkA1 = 1
        if self.kA_char1.param == 'm':
            if not np.isnan(i1_d[0]):
                if not i1[0] == 0:
                    fkA1 = self.kA_char1.func.f_x(i1[0] / i1_d[0])

        fkA2 = 1
        if self.kA_char2.param == 'm':
            if not np.isnan(i2_d[0]):
                if not i2[0] == 0:
                    fkA2 = self.kA_char2.func.f_x(i2[0] / i2_d[0])

        td_log = ((T_o1 - T_i2 - T_i1 + T_o2) /
                  np.log((T_o1 - T_i2) / (T_i1 - T_o2)))
        return i1[0] * (o1[2] - i1[2]) + self.kA.val * fkA1 * fkA2 * td_log

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

    def ttd_u_deriv(self):
        r"""
        Calculates the matrix of partial derivatives for upper temperature
        difference equation.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 4, len(self.inl[0].fluid.val) + 3))
        for i in range(2):
            deriv[0, i * 3, 1] = self.numeric_deriv(self.ttd_u_func,
                                                    'p', i * 3)
            deriv[0, i * 3, 2] = self.numeric_deriv(self.ttd_u_func,
                                                    'h', i * 3)
        return deriv.tolist()

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

    def ttd_l_deriv(self):
        r"""
        Calculates the matrix of partial derivatives for lower temperature
        difference equation.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 4, len(self.inl[0].fluid.val) + 3))
        for i in range(2):
            deriv[0, i + 1, 1] = self.numeric_deriv(self.ttd_l_func,
                                                    'p', i + 1)
            deriv[0, i + 1, 2] = self.numeric_deriv(self.ttd_l_func,
                                                    'h', i + 1)
        return deriv.tolist()

    def bus_func(self, bus):
        r"""
        Calculates the residual value of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus
            TESPy bus object.

        Returns
        -------
        val : float
            Residual value of equation.

            .. math::

                val = P \cdot f\left( \frac{P}{P_{ref}}\right)

                P = \dot{m}_{1,in} \cdot \left( h_{1,out} - h_{1,in} \right)
        """
        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()

        val = i[0] * (o[2] - i[2])
        if np.isnan(bus.P_ref):
            expr = 1
        else:
            expr = abs(val / bus.P_ref)
        return val * bus.char.f_x(expr)

    def bus_deriv(self, bus):
        r"""
        Calculates the matrix of partial derivatives of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus
            TESPy bus object.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 4, len(self.inl[0].fluid.val) + 3))
        deriv[0, 0, 0] = self.numeric_deriv(self.bus_func, 'm', 0, bus=bus)
        deriv[0, 0, 2] = self.numeric_deriv(self.bus_func, 'h', 0, bus=bus)
        deriv[0, 2, 2] = self.numeric_deriv(self.bus_func, 'h', 2, bus=bus)
        return deriv

    def convergence_check(self, nw):
        r"""
        Performs a convergence check.

        Parameters
        ----------
        nw : tespy.networks.network
            The network object using this component.

        Note
        ----
        Manipulate enthalpies/pressure at inlet and outlet if not specified by
        user to match physically feasible constraints, keep fluid composition
        within feasible range and then propagates it towards the outlet.
        """
        i, o = self.inl, self.outl

        if self.ttd_l.is_set or self.ttd_u.is_set:
            fl_i1 = single_fluid(i[0].fluid.val)
            fl_i2 = single_fluid(i[1].fluid.val)
            fl_o1 = single_fluid(o[0].fluid.val)
            fl_o2 = single_fluid(o[1].fluid.val)

        if self.ttd_l.is_set:
            if isinstance(fl_o1, str):
                T_min_o1 = memorise.vrange[fl_o1][2] * 1.1
            else:
                T_min_o1 = nw.T_range_SI[0] * 1.1
            if isinstance(fl_i2, str):
                T_min_i2 = memorise.vrange[fl_i2][2] * 1.1
            else:
                T_min_i2 = nw.T_range_SI[0] * 1.1
            h_min_o1 = h_mix_pT(o[0].to_flow(), T_min_o1)
            h_min_i2 = h_mix_pT(i[1].to_flow(), T_min_i2)
            if not o[0].h.val_set and o[0].h.val_SI < h_min_o1 * 2:
                o[0].h.val_SI = h_min_o1 * 2
            if not i[1].h.val_set and i[1].h.val_SI < h_min_i2:
                i[1].h.val_SI = h_min_i2 * 1.1

        if self.ttd_u.is_set:
            if isinstance(fl_i1, str):
                T_min_i1 = memorise.vrange[fl_i1][2] * 1.1
            else:
                T_min_i1 = nw.T_range_SI[0] * 1.1
            if isinstance(fl_o2, str):
                T_min_o2 = memorise.vrange[fl_o2][2] * 1.1
            else:
                T_min_o2 = nw.T_range_SI[0] * 1.1
            h_min_i1 = h_mix_pT(i[0].to_flow(), T_min_i1)
            h_min_o2 = h_mix_pT(o[1].to_flow(), T_min_o2)
            if not i[0].h.val_set and i[0].h.val_SI < h_min_i1 * 2:
                i[0].h.val_SI = h_min_i1 * 2
            if not o[1].h.val_set and o[1].h.val_SI < h_min_o2:
                o[1].h.val_SI = h_min_o2 * 1.1

    def initialise_source(self, c, key):
        r"""
        Returns a starting value for pressure and enthalpy at component's
        outlet.

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
            flow = [c.m.val0, c.p.val_SI, c.h.val_SI, c.fluid.val]
            if c.s_id == 'out1':
                T = 200 + 273.15
                return h_mix_pT(flow, T)
            else:
                T = 250 + 273.15
                return h_mix_pT(flow, T)

    def initialise_target(self, c, key):
        r"""
        Returns a starting value for pressure and enthalpy at component's
        inlet.

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
            flow = [c.m.val0, c.p.val_SI, c.h.val_SI, c.fluid.val]
            if c.t_id == 'in1':
                T = 300 + 273.15
                return h_mix_pT(flow, T)
            else:
                T = 220 + 273.15
                return h_mix_pT(flow, T)

    def calc_parameters(self):
        r"""
        Postprocessing parameter calculation.
        """
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

        if self.kA.is_set:
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
    Equations

        **mandatory equations**

        - :func:`tespy.components.components.heat_exchanger.fluid_func`
        - :func:`tespy.components.components.heat_exchanger.mass_flow_func`
        - :func:`tespy.components.components.condenser.energy_func`

        **optional equations**

        .. math::

            0 = \dot{m}_{in} \cdot \left(h_{out} - h_{in} \right) - \dot{Q}

        - :func:`tespy.components.components.condenser.kA_func`
        - :func:`tespy.components.components.condenser.ttd_u_func`

        .. math::

            0 = p_{1,in} \cdot pr1 - p_{1,out}\\
            0 = p_{2,in} \cdot pr2 - p_{2,out}

        - :func:`tespy.components.components.heat_exchanger.zeta_func`
        - :func:`tespy.components.components.heat_exchanger.zeta2_func`

        **additional equations**

        - :func:`tespy.components.components.condenser.additional_equations`

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

    Q : Sring/float/tespy.helpers.dc_cp
        Heat transfer, :math:`Q/\text{W}`.

    pr1 : Sring/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio at hot side, :math:`pr/1`.

    pr2 : Sring/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio at cold side, :math:`pr/1`.

    zeta1 : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient at hot side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    zeta2 : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient at cold side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    kA : str/float/tespy.helpers.dc_cp
        Area independent heat transition coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char1 : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient at hot side, provide
        x and y values or use generic values (e. g. calculated from design
        case). Standard method 'COND_HOT', Parameter 'm'.

    kA_char2 : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient at cold side,
        provide x and y values or use generic values (e. g. calculated from
        design case). Standard method 'COND_COLD', Parameter 'm'.

    subcooling : bool
        Enable/disable subcooling, default value: disabled.

    Note
    ----

    - The condenser has an additional equation for enthalpy at hot side outlet.
    - The pressure drop via zeta1 at hot side is not an offdesign parameter.
    - It has different calculation method for given heat transfer coefficient
      and upper terminal temperature difference.

    Example
    -------
    >>> from tespy import cmp, con, nwk
    >>> import shutil
    >>> nw = nwk.network(fluids=['water', 'air'], T_unit='C', p_unit='bar',
    ...     h_unit='kJ / kg', m_range=[0.01, 10])
    >>> nw.set_printoptions(print_level='none')
    >>> amb_in = cmp.sink('ambient in')
    >>> amb_out = cmp.source('ambient out')
    >>> hsin = cmp.sink('HS in')
    >>> hsout = cmp.source('HS out')
    >>> he = cmp.condenser('condenser')
    >>> he.component()
    'condenser'
    >>> amb_he = con.connection(amb_out, 'out1', he, 'in2')
    >>> he_amb = con.connection(he, 'out2', amb_in, 'in1')
    >>> hs_he = con.connection(hsout, 'out1', he, 'in1')
    >>> he_hs = con.connection(he, 'out1', hsin, 'in1')
    >>> nw.add_conns(amb_he, he_amb, hs_he, he_hs)
    >>> he.set_attr(pr1=0.98, pr2=0.999, design=['pr2'],
    ...     offdesign=['zeta2', 'kA'])
    >>> hs_he.set_attr(Td_bp=20, p=1, fluid={'water': 1, 'air': 0})
    >>> amb_he.set_attr(fluid={'water': 0, 'air': 1}, T=20)
    >>> he_amb.set_attr(p=1, T=40, design=['T'])
    >>> he.set_attr(Q=-80e3)
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(hs_he.m.val, 2)
    0.03
    >>> round(amb_he.m.val, 2)
    3.97
    >>> round(he_amb.T.val, 1)
    40.0
    >>> he.set_attr(Q=-60e3)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(amb_he.m.val, 2)
    2.78
    >>> round(he_amb.T.val, 1)
    41.5
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    def component(self):
        return 'condenser'

    def attr(self):
        return {'Q': dc_cp(max_val=0),
                'kA': dc_cp(min_val=0),
                'td_log': dc_cp(min_val=0),
                'ttd_u': dc_cp(min_val=0), 'ttd_l': dc_cp(min_val=0),
                'pr1': dc_cp(max_val=1), 'pr2': dc_cp(max_val=1),
                'zeta1': dc_cp(min_val=0), 'zeta2': dc_cp(min_val=0),
                'subcooling': dc_simple(val=False),
                'kA_char1': dc_cc(method='COND_HOT', param='m'),
                'kA_char2': dc_cc(method='COND_COLD', param='m'),
                'SQ1': dc_simple(), 'SQ2': dc_simple(), 'Sirr': dc_simple(),
                'zero_flag': dc_simple()}

    def additional_equations(self):
        r"""
        Calculates vector vec_res with results of additional equations for this
        component.

        Equations

            **mandatory equations**

            .. math::

                0 = h_{1,out} - h\left(p, x=0 \right)\\
                x: \text{vapour mass fraction}

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # equation for saturated liquid at hot side outlet
        if not self.subcooling.val:
            outl = self.outl
            o1 = outl[0].to_flow()
            vec_res += [o1[2] - h_mix_pQ(o1, 0)]

        return vec_res

    def additional_derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given additional equations.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        mat_deriv = []

        ######################################################################
        # derivatives for saturated liquid at hot side outlet equation
        if not self.subcooling.val:
            o1 = self.outl[0].to_flow()
            x_deriv = np.zeros((1, 4, self.num_fl + 3))
            x_deriv[0, 2, 1] = -dh_mix_dpQ(o1, 0)
            x_deriv[0, 2, 2] = 1
            mat_deriv += x_deriv.tolist()

        return mat_deriv

    def energy_func(self):
        r"""
        Equation for condenser energy balance.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_{1,in} \cdot \left(h_{1,out} - h_{1,in} \right) +
                \dot{m}_{2,in} \cdot \left(h_{2,out} - h_{2,in} \right)
        """
        return (self.inl[0].m.val_SI * (self.outl[0].h.val_SI -
                                        self.inl[0].h.val_SI) +
                self.inl[1].m.val_SI * (self.outl[1].h.val_SI -
                                        self.inl[1].h.val_SI))

    def energy_deriv(self):
        r"""
        Calculates the matrix of partial derivatives for energy balance
        equation.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 4, len(self.inl[0].fluid.val) + 3))
        for k in range(2):
            deriv[0, k, 0] = self.outl[k].h.val_SI - self.inl[k].h.val_SI
            deriv[0, k, 2] = -self.inl[k].m.val_SI

        deriv[0, 2, 2] = self.inl[0].m.val_SI
        deriv[0, 3, 2] = self.inl[1].m.val_SI
        return deriv.tolist()

    def kA_func(self):
        r"""
        Equation for heat transfer from conditions on both sides of condenser.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                res = \dot{m}_{1,in} \cdot \left( h_{1,out} - h_{1,in}\right) +
                kA \cdot f_{kA} \cdot \frac{T_{1,out} -
                T_{2,in} - T_s \left(p_{1,in}\right) +
                T_{2,out}}
                {\ln{\frac{T_{1,out} - T_{2,in}}
                {T_s \left(p_{1,in}\right) - T_{2,out}}}}

                f_{kA} = f_1\left(\frac{m_1}{m_{1,ref}}\right) \cdot
                f_2\left(\frac{m_2}{m_{2,ref}}\right)

        Note
        ----
        For standard functions f\ :subscript:`1` \ and f\ :subscript:`2` \ see
        class :func:`tespy.components.characteristics.characteristics`.

        - Calculate temperatures at inlets and outlets.
        - Perform value manipulation, if temperature levels are physically
          infeasible.
        """
        if self.zero_flag.val_set:
            return self.inl[0].p.val_SI - self.inl[0].p.design

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
            if not np.isnan(i1_d[0]):
                fkA1 = self.kA_char1.func.f_x(i1[0] / i1_d[0])

        fkA2 = 1
        if self.kA_char2.param == 'm':
            if not np.isnan(i2_d[0]):
                fkA2 = self.kA_char2.func.f_x(i2[0] / i2_d[0])

        td_log = ((T_o1 - T_i2 - T_i1 + T_o2) /
                  np.log((T_o1 - T_i2) / (T_i1 - T_o2)))
        return i1[0] * (o1[2] - i1[2]) + self.kA.val * fkA1 * fkA2 * td_log

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
    Equations

        **mandatory equations**

        - :func:`tespy.components.components.heat_exchanger.fluid_func`
        - :func:`tespy.components.components.heat_exchanger.mass_flow_func`
        - :func:`tespy.components.components.heat_exchanger.energy_func`

        **optional equations**

        .. math::

            0 = \dot{m}_{in} \cdot \left(h_{out} - h_{in} \right) - \dot{Q}

        - :func:`tespy.components.components.heat_exchanger.kA_func`
        - :func:`tespy.components.components.heat_exchanger.ttd_u_func`

        .. math::

            0 = p_{1,in} \cdot pr1 - p_{1,out}\\
            0 = p_{2,in} \cdot pr2 - p_{2,out}

        - :func:`tespy.components.components.heat_exchanger.zeta_func`
        - :func:`tespy.components.components.heat_exchanger.zeta2_func`

        **additional equations**

        - :func:`tespy.components.components.desuperheater.additional_equations`

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

    Q : Sring/float/tespy.helpers.dc_cp
        Heat transfer, :math:`Q/\text{W}`.

    pr1 : Sring/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio at hot side, :math:`pr/1`.

    pr2 : Sring/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio at cold side, :math:`pr/1`.

    zeta1 : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient at hot side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    zeta2 : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient at cold side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    kA : str/float/tespy.helpers.dc_cp
        Area independent heat transition coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char1 : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient at hot side, provide
        x and y values or use generic values (e. g. calculated from design
        case). Standard method 'COND_HOT', Parameter 'm'.

    kA_char2 : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient at cold side,
        provide x and y values or use generic values (e. g. calculated from
        design case). Standard method 'COND_COLD', Parameter 'm'.

    Note
    ----
    The desuperheater has an additional equation for enthalpy at hot side
    outlet.

    Example
    -------
    >>> from tespy import cmp, con, nwk
    >>> import shutil
    >>> nw = nwk.network(fluids=['water', 'air'], T_unit='C', p_unit='bar',
    ...     h_unit='kJ / kg')
    >>> nw.set_printoptions(print_level='none')
    >>> amb_in = cmp.sink('ambient in')
    >>> amb_out = cmp.source('ambient out')
    >>> hsin = cmp.sink('HS in')
    >>> hsout = cmp.source('HS out')
    >>> he = cmp.desuperheater('desuperheater')
    >>> he.component()
    'desuperheater'
    >>> amb_he = con.connection(amb_out, 'out1', he, 'in2')
    >>> he_amb = con.connection(he, 'out2', amb_in, 'in1')
    >>> hs_he = con.connection(hsout, 'out1', he, 'in1')
    >>> he_hs = con.connection(he, 'out1', hsin, 'in1')
    >>> nw.add_conns(amb_he, he_amb, hs_he, he_hs)
    >>> he.set_attr(pr1=0.98, pr2=0.999, design=['pr1', 'pr2'],
    ...     offdesign=['zeta1', 'zeta2', 'kA'])
    >>> hs_he.set_attr(T=200, p=1, fluid={'water': 1, 'air': 0})
    >>> amb_he.set_attr(fluid={'water': 0, 'air': 1}, T=20)
    >>> he_amb.set_attr(p=1, T=40, design=['T'])
    >>> he.set_attr(Q=-80e3)
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(hs_he.m.val, 1)
    0.4
    >>> round(amb_he.m.val, 2)
    3.97
    >>> round(he_amb.T.val, 1)
    40.0
    >>> he.set_attr(Q=-60e3)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(hs_he.m.val, 1)
    0.3
    >>> round(amb_he.m.val, 2)
    2.56
    >>> round(he_amb.T.val, 1)
    43.3
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    def component(self):
        return 'desuperheater'

    def additional_equations(self):
        r"""
        Calculates vector vec_res with results of additional equations for this
        component.

        Equations

            **mandatory equations**

            .. math::

                0 = h_{1,out} - h\left(p, x=1 \right)\\
                x: \text{vapour mass fraction}

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # equation for saturated gas at hot side outlet
        o1 = self.outl[0].to_flow()
        vec_res += [o1[2] - h_mix_pQ(o1, 1)]

        return vec_res

    def additional_derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given additional
        equations.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        mat_deriv = []

        ######################################################################
        # derivatives for saturated gas at hot side outlet equation
        o1 = self.outl[0].to_flow()
        deriv = np.zeros((1, 4, self.num_fl + 3))
        deriv[0, 2, 1] = -dh_mix_dpQ(o1, 1)
        deriv[0, 2, 2] = 1
        mat_deriv += deriv.tolist()

        return mat_deriv

# %%


class solar_collector(heat_exchanger_simple):
    r"""
    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_func`
        - :func:`tespy.components.components.component.mass_flow_func`

        **optional equations**

        - :func:`tespy.components.components.heat_exchanger_simple.Q_func`

        .. math::

            0 = p_{in} \cdot pr - p_{out}

        - :func:`tespy.components.components.component.zeta_func`

        - :func:`tespy.components.components.heat_exchanger_simple.darcy_func`
          or :func:`tespy.components.components.heat_exchanger_simple.hw_func`

        **additional equations**

        - :func:`tespy.components.components.solar_collector.additional_equations`

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

    Q : Sring/float/tespy.helpers.dc_cp
        Heat transfer, :math:`Q/\text{W}`.

    pr : Sring/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio, :math:`pr/1`.

    zeta : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    D : str/float/tespy.helpers.dc_cp
        Diameter of the pipes, :math:`D/\text{m}`.

    L : str/float/tespy.helpers.dc_cp
        Length of the pipes, :math:`L/\text{m}`.

    ks : str/float/tespy.helpers.dc_cp
        Pipes roughness, :math:`ks/\text{m}` for darcy friction,
        :math:`ks/\text{1}` for hazen-williams equation.

    hydro_group : Sring/tespy.helpers.dc_gcp
        Parametergroup for pressure drop calculation based on pipes dimensions.
        Choose 'HW' for hazen-williams equation, else darcy friction factor is
        used.

    E : str/float/tespy.helpers.dc_cp
        Radiation at tilted collector surface area,
        :math:`E/\frac{\text{W}}{\text{m}^2}`.

    eta_opt : str/float/tespy.helpers.dc_cp
        optical loss at surface cover,
        :math:`\eta_{opt}`.

    lkf_lin : str/float/tespy.helpers.dc_cp
        Linear loss key figure,
        :math:`\alpha_1/\frac{\text{W}}{\text{K} \cdot \text{m}^2}`.

    lkf_quad : str/float/tespy.helpers.dc_cp
        Quadratic loss key figure,
        :math:`\alpha_2/\frac{\text{W}}{\text{K}^2 \cdot \text{m}^2}`.

    A : str/float/tespy.helpers.dc_cp
        Collector surface area :math:`A/\text{m}^2`.

    Tamb : float/tespy.helpers.dc_cp
        Ambient temperature, provide parameter in network's temperature unit.

    energy_group : tespy.helpers.dc_gcp
        Parametergroup for energy balance of solarthermal collector.

    Example
    -------
    The solar collector is used to calculate heat transferred to the heating
    system from radiation on a tilted plane. For instance, it is possible to
    calculate the collector surface area required to transfer a specific amount
    of heat at a given radiation. The collector parameters are the linear and
    the quadratic loss keyfigure as well as the optical effifiency.

    >>> from tespy.components.basics import sink, source
    >>> from tespy.components.heat_exchangers import solar_collector
    >>> from tespy.connections import connection
    >>> from tespy.networks.networks import network
    >>> import shutil
    >>> fluids = ['H2O']
    >>> nw = network(fluids=fluids)
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')
    >>> nw.set_printoptions(print_level='none')
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
    different radiation) using the calculated surface area and mass flow, it
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

    def component(self):
        return 'solar collector'

    def attr(self):
        return {'Q': dc_cp(),
                'pr': dc_cp(min_val=1e-4, max_val=1),
                'zeta': dc_cp(min_val=1e-4),
                'D': dc_cp(min_val=1e-2, max_val=2, d=1e-3),
                'L': dc_cp(min_val=1e-1, d=1e-3),
                'ks': dc_cp(val=1e-4, min_val=1e-7, max_val=1e-4, d=1e-8),
                'E': dc_cp(min_val=0),
                'eta_opt': dc_cp(min_val=0, max_val=1),
                'lkf_lin': dc_cp(min_val=0),
                'lkf_quad': dc_cp(min_val=0),
                'A': dc_cp(min_val=0),
                'Tamb': dc_cp(),
                'SQ': dc_simple(),
                'hydro_group': dc_gcp(), 'energy_group': dc_gcp()}

    def inlets(self):
        return ['in1']

    def outlets(self):
        return ['out1']

    def comp_init(self, nw):

        component.comp_init(self, nw)

        self.fl_deriv = self.fluid_deriv()
        self.m_deriv = self.mass_flow_deriv()

        self.Tamb.val_SI = ((self.Tamb.val + nw.T[nw.T_unit][0]) *
                            nw.T[nw.T_unit][1])

        # parameters for hydro group
        self.hydro_group.set_attr(elements=[self.L, self.ks, self.D])

        is_set = True
        for e in self.hydro_group.elements:
            if not e.is_set:
                is_set = False

        if is_set:
            self.hydro_group.set_attr(is_set=True)
        elif self.hydro_group.is_set:
            msg = ('All parameters of the component group have to be '
                   'specified! This component group uses the following '
                   'parameters: L, ks, D at ' + self.label + '. '
                   'Group will be set to False.')
            logging.info(msg)
            self.hydro_group.set_attr(is_set=False)
        else:
            self.hydro_group.set_attr(is_set=False)

        # parameters for energy group
        self.energy_group.set_attr(elements=[self.E, self.eta_opt, self.lkf_lin,
                                             self.lkf_quad, self.A, self.Tamb])

        is_set = True
        for e in self.energy_group.elements:
            if not e.is_set:
                is_set = False

        if is_set:
            self.energy_group.set_attr(is_set=True)
        elif self.energy_group.is_set:
            msg = ('All parameters of the component group have to be '
                   'specified! This component group uses the following '
                   'parameters: E, eta_opt, lkf_lin, lkf_quad, A, Tamb at ' +
                   self.label + '. Group will be set to False.')
            logging.info(msg)
            self.energy_group.set_attr(is_set=False)
        else:
            self.energy_group.set_attr(is_set=False)

    def additional_equations(self):
        r"""
        Calculates vector vec_res with results of additional equations for this
        component.

        Equations

            **optional equations**

            - :func:`tespy.components.components.solar_collector.energy_func`

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # equation for specified energy-group paremeters
        if self.energy_group.is_set:
            vec_res += [self.energy_func()]

        return vec_res

    def additional_derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given additional
        equations.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        mat_deriv = []

        ######################################################################
        # derivatives for specified energy-group paremeters
        if self.energy_group.is_set:
            deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))
            deriv[0, 0, 0] = self.outl[0].h.val_SI - self.inl[0].h.val_SI
            deriv[0, 0, 1] = self.numeric_deriv(self.energy_func, 'p', 0)
            deriv[0, 0, 2] = self.numeric_deriv(self.energy_func, 'h', 0)
            deriv[0, 1, 1] = self.numeric_deriv(self.energy_func, 'p', 1)
            deriv[0, 1, 2] = self.numeric_deriv(self.energy_func, 'h', 1)
            # custom variables for the energy-group
            for var in self.energy_group.elements:
                if var.is_var:
                    deriv[0, 2 + var.var_pos, 0] = (
                            self.numeric_deriv(self.energy_func,
                                               self.vars[var], 2))
            mat_deriv += deriv.tolist()

        return mat_deriv

    def energy_func(self):
        r"""
        Equation for solar collector energy balance.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                T_m = \frac{T_{out} + T_{in}}{2}\\

                \begin{split}
                0 = & \dot{m} \cdot \left( h_{out} - h_{in} \right)\\
                & - A \cdot \left[E \cdot \eta_{opt} - \alpha_1 \cdot
                \left(T_m - T_{amb} \right) - \alpha_2 \cdot
                \left(T_m - T_{amb}\right)^2 \right]
                \end{split}
        """

        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()

        T_m = (T_mix_ph(i, T0=self.inl[0].T.val_SI) +
               T_mix_ph(o, T0=self.outl[0].T.val_SI)) / 2

        return (i[0] * (o[2] - i[2]) - self.A.val * (self.E.val *
                self.eta_opt.val -(T_m - self.Tamb.val_SI) * self.lkf_lin.val -
                self.lkf_quad.val * (T_m - self.Tamb.val_SI) ** 2))

    def calc_parameters(self):
        r"""
        Postprocessing parameter calculation.
        """
        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()

        self.SQ.val = i[0] * (s_mix_ph(o) - s_mix_ph(i))
        self.Q.val = i[0] * (o[2] - i[2])
        self.pr.val = o[1] / i[1]
        self.zeta.val = ((i[1] - o[1]) * np.pi ** 2 /
                         (8 * i[0] ** 2 * (v_mix_ph(i) + v_mix_ph(o)) / 2))

        self.check_parameter_bounds()
