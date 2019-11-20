# -*- coding: utf-8

"""This module contains reactors: water_electrolyzer


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/reactors.py

SPDX-License-Identifier: MIT
"""

import CoolProp.CoolProp as CP

import logging

import numpy as np

from tespy.components.components import component

from tespy.tools.data_containers import dc_cc, dc_cp, dc_simple
from tespy.tools.fluid_properties import (
        h_mix_pT, T_mix_ph, dT_mix_dph, dT_mix_pdh, v_mix_ph)
from tespy.tools.global_vars import molar_masses

# %%


class water_electrolyzer(component):
    r"""
    Equations

        **mandatory equations**

        .. math::

            0  = x_{i,in1} - x_{i,out1} \forall i \in \text{fluids}

            \forall i \in \text{network fluids}:

            0 = \begin{cases}
                1 - x_{i,in2} & \text{i=}H_{2}O\\
                x_{i,in2} & \text{else}
            \end{cases}\\

            0 = \begin{cases}
                1 - x_{i,out2} & \text{i=}O_{2}\\
                x_{i,out2} & \text{else}
            \end{cases}\\

            0 = \begin{cases}
                1 - x_{i,out3} & \text{i=}H_{2}\\
                x_{i,out3} & \text{else}
            \end{cases}\\

            O_2 = \frac{M_{O_2}}{M_{O_2} + 2 \cdot M_{H_2}}\\

            0 = \dot{m}_{H_{2}O,in1} - \dot{m}_{H_{2}O,out1}\\
            0 = O_2 \cdot \dot{m}_{H_{2}O,in2} - \dot{m}_{O_2,out2}\\
            0 = \left(1 - O_2\right) \cdot \dot{m}_{H_{2}O,in2} -
            \dot{m}_{H_2,out3}\\

            0 = p_{H_{2}O,in2} - p_{O_2,out2}\\
            0 = p_{H_{2}O,in2} - p_{H_2,out3}

            0 = P - f_{eb}\left( \right)

        For energy balance (f_eb) calculation see
        :func:`tespy.components.components.water_electrolyzer.energy_balance`.

        .. math::

            0 = T_{O_2,out2} - T_{H_2,out3}

        **optional equations**

        .. math::

            0 = P - \dot{m}_{H_2,out3} \cdot e\\

            0 = p_{H_{2}O,in1} \cdot pr - p_{H_{2}O,out1}

        - :func:`tespy.components.components.component.zeta_func`

        .. math::

            0 = \dot{Q} - \dot{m}_{in1} \cdot \left(h_{in1} -
            h_{out1}\right) \\

            0 = P - \dot{m}_{H_2,out3} \cdot \frac{e_0}{\eta}

        - :func:`tespy.components.components.water_electrolyzer.eta_char_func`

    Inlets/Outlets

        - in1 (cooling inlet), in2 (feed water inlet)
        - out1 (cooling outlet), out2 (hydrogen outlet), out3 (oxigen outlet)

    Image

        .. image:: _images/water_electrolyzer.svg
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

    P : float/tespy.helpers.dc_cp
        Power input, :math:`P/\text{W}`.

    Q : float/tespy.helpers.dc_cp
        Heat output of cooling, :math:`Q/\text{W}`

    e : float/tespy.helpers.dc_cp
        Electrolysis specific energy consumption,
        :math:`e/(\text{J}/\text{m}^3)`.

    eta : float/tespy.helpers.dc_cp
        Electrolysis efficiency, :math:`\eta/1`.

    eta_char : str/tespy.helpers.dc_cc
        Electrolysis efficiency characteristic line.

    pr : float/tespy.helpers.dc_cp
        Cooling loop pressure ratio, :math:`pr/1`.

    zeta : float/tespy.helpers.dc_cp
        Geometry independent friction coefficient for cooling loop pressure
        drop, :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    Note
    ----

    Example
    -------
    >>> from tespy import cmp, con, nwk
    >>> import shutil
    >>> fluid_list = ['O2', 'water', 'H2']
    >>> nw = nwk.network(fluids=fluid_list, T_unit='C', p_unit='bar',
    ... h_unit='kJ / kg')
    >>> nw.set_printoptions(print_level='none')

    >>> fw = cmp.source('feed water')
    >>> oxy = cmp.sink('oxygen sink')
    >>> hydro = cmp.sink('hydrogen sink')
    >>> cw = cmp.source('cooling water')
    >>> cw_hot = cmp.sink('cooling water out')

    >>> el = cmp.water_electrolyzer('electrolyzer 1', eta=0.8, design=['eta'],
    ... offdesign=['eta_char'])
    >>> el.component()
    'water electrolyzer'
    >>> comp = cmp.compressor('compressor', eta_s=0.9)

    >>> fw_el = con.connection(fw, 'out1', el, 'in2', m=0.1, p=10, T=15)
    >>> el_o = con.connection(el, 'out2', oxy, 'in1')
    >>> el_cmp = con.connection(el, 'out3', comp, 'in1', T=50)
    >>> cmp_h = con.connection(comp, 'out1', hydro, 'in1', p=50)
    >>> cw_el = con.connection(cw, 'out1', el, 'in1', p=5, T=15,
    ... fluid={'water': 1, 'H2': 0, 'O2': 0})
    >>> el_cw = con.connection(el, 'out1', cw_hot, 'in1', T=45, p=4.9)
    >>> nw.add_conns(fw_el, el_o, el_cmp, cmp_h, cw_el, el_cw)
    >>> nw.solve('design')
    >>> round(el.eta.val, 1)
    0.8
    >>> nw.save('tmp')
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(el.eta.val, 1)
    0.8
    >>> fw_el.set_attr(m=0.05)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(el.eta.val, 2)
    0.82
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    def component(self):
        return 'water electrolyzer'

    def attr(self):
        return {'P': dc_cp(min_val=0),
                'Q': dc_cp(max_val=0),
                'eta': dc_cp(min_val=0, max_val=1),
                'e': dc_cp(),
                'pr_c': dc_cp(max_val=1),
                'zeta': dc_cp(min_val=0),
                'eta_char': dc_cc(method='GENERIC'),
                'S': dc_simple()}

    def inlets(self):
        return ['in1', 'in2']

    def outlets(self):
        return ['out1', 'out2', 'out3']

    def comp_init(self, nw):

        if not self.P.is_set:
            self.set_attr(P='var')
            msg = ('The power output of cogeneration units must be set! '
                   'We are adding the power output of component ' +
                   self.label + ' as custom variable of the system.')
            logging.info(msg)

        component.comp_init(self, nw)

        o2 = [x for x in nw.fluids if x in [a.replace(' ', '')
              for a in CP.get_aliases('O2')]]
        if len(o2) == 0:
            msg = ('Missing oxygen in network fluids, component ' +
                   self.label + ' of type ' + self.component() +
                   ' requires oxygen in network fluids.')
            logging.error(msg)
            raise ValueError(msg)
        else:
            self.o2 = o2[0]

        h2o = [x for x in nw.fluids if x in [a.replace(' ', '')
               for a in CP.get_aliases('H2O')]]
        if len(h2o) == 0:
            msg = ('Missing water in network fluids, component ' +
                   self.label + ' of type ' + self.component() +
                   ' requires water in network fluids.')
            logging.error(msg)
            raise ValueError(msg)
        else:
            self.h2o = h2o[0]

        h2 = [x for x in nw.fluids if x in [a.replace(' ', '')
              for a in CP.get_aliases('H2')]]
        if len(h2) == 0:
            msg = ('Missing hydrogen in network fluids, component ' +
                   self.label + ' of type ' + self.component() +
                   ' requires hydrogen in network fluids.')
            logging.error(msg)
            raise ValueError(msg)
        else:
            self.h2 = h2[0]

        self.e0 = self.calc_e0()

    def calc_e0(self):
        r"""
        Calculates the minimum specific energy required for electrolysis.

        Returns
        -------
        val : float
            Minimum specific energy.

            .. math::
                LHV = -\frac{\sum_i {\Delta H_f^0}_i -
                \sum_j {\Delta H_f^0}_j }
                {M_{fuel}}\\
                \forall i \in \text{reation products},\\
                \forall j \in \text{reation educts},\\
                \Delta H_f^0: \text{molar formation enthalpy}
        """

        hf = {}
        hf['H2O'] = -286
        hf['H2'] = 0
        hf['O2'] = 0
        M = molar_masses['H2']
        e0 = -(2 * hf['H2O'] - 2 * hf['H2'] + hf['O2']) / (2 * M)

        return e0 * 1000

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
        # equations for fluids

        # equations for fluid composition in cooling water
        for fluid, x in self.inl[0].fluid.val.items():
            vec_res += [x - self.outl[0].fluid.val[fluid]]

        # equations to constrain fluids to inlets/outlets
        vec_res += [1 - self.inl[1].fluid.val[self.h2o]]
        vec_res += [1 - self.outl[1].fluid.val[self.o2]]
        vec_res += [1 - self.outl[2].fluid.val[self.h2]]

        # equations to ban fluids off inlets/outlets
        for fluid in self.inl[1].fluid.val.keys():
            if fluid != self.h2o:
                vec_res += [0 - self.inl[1].fluid.val[fluid]]
            if fluid != self.o2:
                vec_res += [0 - self.outl[1].fluid.val[fluid]]
            if fluid != self.h2:
                vec_res += [0 - self.outl[2].fluid.val[fluid]]

        ######################################################################
        # eqations for mass flow balance
        # equation to calculate the ratio of o2 in water
        o2 = molar_masses[self.o2] / (molar_masses[self.o2] +
                                      2 * molar_masses[self.h2])

        # equation for mass flow balance cooling water
        vec_res += [self.inl[0].m.val_SI - self.outl[0].m.val_SI]

        # equations for mass flow balance electrolyzer
        vec_res += [o2 * self.inl[1].m.val_SI - self.outl[1].m.val_SI]
        vec_res += [(1 - o2) * self.inl[1].m.val_SI - self.outl[2].m.val_SI]

        ######################################################################
        # equations for pressure to set o2 and h2 output equal
        vec_res += [self.inl[1].p.val_SI - self.outl[1].p.val_SI]
        vec_res += [self.inl[1].p.val_SI - self.outl[2].p.val_SI]

        ######################################################################
        # equation for energy balance
        vec_res += [self.P.val + self.energy_balance()]

        ######################################################################
        # temperature electrolyzer outlet
        vec_res += [T_mix_ph(self.outl[1].to_flow()) -
                    T_mix_ph(self.outl[2].to_flow())]

        ######################################################################
        # power vs hydrogen production
        if self.e.is_set:
            vec_res += [self.P.val - self.outl[2].m.val_SI * self.e.val]

        ######################################################################
        #pr_c.val = pressure ratio Druckverlust (als Faktor vorgegeben)
        if self.pr_c.is_set:
            vec_res += [self.inl[0].p.val_SI * self.pr_c.val -
                        self.outl[0].p.val_SI]

        if self.zeta.is_set:
            vec_res += [self.zeta_func()]

        # equation for heat transfer

        if self.Q.is_set:
            vec_res += [self.Q.val - self.inl[0].m.val_SI *
                        (self.inl[0].h.val_SI - self.outl[0].h.val_SI)]

        ######################################################################
        # specified efficiency (efficiency definition: e0 / e)
        if self.eta.is_set:
            vec_res += [self.P.val - self.outl[2].m.val_SI *
                        self.e0 / self.eta.val]

        ######################################################################
        # specified characteristic line for efficiency
        if self.eta_char.is_set:
            vec_res += [self.eta_char_func()]

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
        # derivatives for cooling liquid composition
        deriv = np.zeros((self.num_fl, 5 + self.num_vars, self.num_fl + 3))

        j = 0
        for fluid, x in self.inl[0].fluid.val.items():
            deriv[j, 0, 3 + j] = 1
            deriv[j, 2, 3 + j] = -1
            j += 1

        mat_deriv += deriv.tolist()

        # derivatives to constrain fluids to inlets/outlets
        deriv = np.zeros((3, 5 + self.num_vars, self.num_fl + 3))

        i = 0
        for fluid in self.inl[0].fluid.val.keys():
            if fluid == self.h2o:
                deriv[0, 1, 3 + i] = -1
            elif fluid == self.o2:
                deriv[1, 3, 3 + i] = -1
            elif fluid == self.h2:
                deriv[2, 4, 3 + i] = -1
            i += 1

        mat_deriv += deriv.tolist()

        # derivatives to ban fluids off inlets/outlets
        deriv = np.zeros((3 * len(self.inl[1].fluid.val.keys()) - 3,
                          5 + self.num_vars, self.num_fl + 3))

        i = 0
        j = 0
        for fluid in self.inl[1].fluid.val.keys():
            if fluid != self.h2o:
                deriv[j, 1, 3 + i] = -1
                j += 1
            if fluid != self.o2:
                deriv[j, 3, 3 + i] = -1
                j += 1
            if fluid != self.h2:
                deriv[j, 4, 3 + i] = -1
                j += 1
            i += 1

        mat_deriv += deriv.tolist()

        ######################################################################
        # derivatives for mass balance equations

        # deritatives for mass flow balance in the heat exchanger
        deriv = np.zeros((3, 5 + self.num_vars, self.num_fl + 3))

        deriv[0, 0, 0] = 1
        deriv[0, 2, 0] = -1

        # derivatives for mass flow balance for oxygen output
        o2 = molar_masses[self.o2] / (molar_masses[self.o2] +
                                      2 * molar_masses[self.h2])
        deriv[1, 1, 0] = o2
        deriv[1, 3, 0] = -1

        # derivatives for mass flow balance for hydrogen output
        deriv[2, 1, 0] = (1 - o2)
        deriv[2, 4, 0] = -1

        mat_deriv += deriv.tolist()

        ######################################################################
        # derivatives for pressure equations

        # derivatives for pressure oxygen outlet
        deriv = np.zeros((2, 5 + self.num_vars, self.num_fl + 3))

        deriv[0, 1, 1] = 1
        deriv[0, 3, 1] = -1

        # derivatives for pressure hydrogen outlet
        deriv[1, 1, 1] = 1
        deriv[1, 4, 1] = -1

        mat_deriv += deriv.tolist()

        ######################################################################
        # derivatives for energy balance equations

        deriv = np.zeros((1, 5 + self.num_vars, self.num_fl + 3))

        T_ref = 293.15
        p_ref = 1e5

        h_refh2o = h_mix_pT([1, p_ref, 0, self.inl[1].fluid.val], T_ref)
        h_refh2 = h_mix_pT([1, p_ref, 0, self.outl[2].fluid.val], T_ref)
        h_refo2 = h_mix_pT([1, p_ref, 0, self.outl[1].fluid.val], T_ref)

        # derivatives cooling water inlet
        deriv[0, 0, 0] = - (self.outl[0].h.val_SI - self.inl[0].h.val_SI)
        deriv[0, 0, 2] = self.inl[0].m.val_SI

        # derivatives feed water inlet
        deriv[0, 1, 0] = (self.inl[1].h.val_SI - h_refh2o)
        deriv[0, 1, 2] = self.inl[1].m.val_SI

        # derivative cooling water outlet
        deriv[0, 2, 2] = - self.inl[0].m.val_SI

        # derivatives oxygen outlet
        deriv[0, 3, 0] = - (self.outl[1].h.val_SI - h_refo2)
        deriv[0, 3, 2] = - self.outl[1].m.val_SI

        # derivatives hydrogen outlet
        deriv[0, 4, 0] = - self.e0 - (self.outl[2].h.val_SI - h_refh2)
        deriv[0, 4, 2] = - self.outl[2].m.val_SI

        # derivatives for variable P
        if self.P.is_var:
            deriv[0, 5 + self.P.var_pos, 0] = 1

        mat_deriv += deriv.tolist()

        ######################################################################
        # derivatives for temperature at gas outlets

        deriv = np.zeros((1, 5 + self.num_vars, self.num_fl + 3))

        # derivatives for outlet 1
        deriv[0, 3, 1] = dT_mix_dph(self.outl[1].to_flow())
        deriv[0, 3, 2] = dT_mix_pdh(self.outl[1].to_flow())

        # derivatives for outlet 2
        deriv[0, 4, 1] = - dT_mix_dph(self.outl[2].to_flow())
        deriv[0, 4, 2] = - dT_mix_pdh(self.outl[2].to_flow())

        mat_deriv += deriv.tolist()

        ######################################################################
        # derivatives for power vs. hydrogen production

        if self.e.is_set:
            deriv = np.zeros((1, 5 + self.num_vars, self.num_fl + 3))

            deriv[0, 4, 0] = - self.e.val

            # derivatives for variable P
            if self.P.is_var:
                deriv[0, 5 + self.P.var_pos, 0] = 1

            # derivatives for variable e
            if self.e.is_var:
                deriv[0, 5 + self.e.var_pos, 0] = - self.outl[2].m.val_SI

            mat_deriv += deriv.tolist()

        ######################################################################
        # derivatives for pressure ratio
        if self.pr_c.is_set:

            deriv = np.zeros((1, 5 + self.num_vars, self.num_fl + 3))

            deriv[0, 0, 1] = self.pr_c.val
            deriv[0, 2, 1] = - 1

            mat_deriv += deriv.tolist()

        ######################################################################
        #pr_c.val = pressure ratio Druckverlust (als Faktor vorgegeben)
        # derivatives for zeta value
        if self.zeta.is_set:

            deriv = np.zeros((1, 5 + self.num_vars, self.num_fl + 3))
            deriv[0, 0, 0] = self.numeric_deriv(self.zeta_func, 'm', 0)
            deriv[0, 0, 1] = self.numeric_deriv(self.zeta_func, 'p', 0)
            deriv[0, 0, 2] = self.numeric_deriv(self.zeta_func, 'h', 0)
            deriv[0, 2, 1] = self.numeric_deriv(self.zeta_func, 'p', 2)
            deriv[0, 2, 2] = self.numeric_deriv(self.zeta_func, 'h', 2)

            # derivatives for variable zeta
            if self.zeta.is_var:
                deriv[0, 5 + self.zeta.var_pos, 0] = (
                        self.numeric_deriv(self.zeta_func, 'zeta', 5))

            mat_deriv += deriv.tolist()

        ######################################################################
        # derivatives for heat flow
        if self.Q.is_set:

            deriv = np.zeros((1, 5 + self.num_vars, self.num_fl + 3))

            deriv[0, 0, 0] = - (self.inl[0].h.val_SI - self.outl[0].h.val_SI)
            deriv[0, 0, 2] = - self.inl[0].m.val_SI
            deriv[0, 2, 2] = self.inl[0].m.val_SI

            mat_deriv += deriv.tolist()

        ######################################################################
        # specified efficiency (efficiency definition: e0 / e)
        if self.eta.is_set:

            deriv = np.zeros((1, 5 + self.num_vars, self.num_fl + 3))

            deriv[0, 4, 0] = - self.e0 / self.eta.val

            # derivatives for variable P
            if self.P.is_var:
                deriv[0, 5 + self.P.var_pos, 0] = 1

            mat_deriv += deriv.tolist()

        ######################################################################
        # specified characteristic line for efficiency
        if self.eta_char.is_set:

            mat_deriv += self.eta_char_deriv()

        ######################################################################

        return np.asarray(mat_deriv)

    def eta_char_func(self):
        r"""
        Equation for given efficiency characteristic of a water electrolyzer.
        Efficiency is linked to hydrogen production.

        Returns
        -------
        res : ndarray
            Residual value of equation.

            .. math::

                0 = P - \dot{m}_{H_2,out3} \cdot \frac{e_0}{\eta_0 \cdot
                f\left(\frac{\dot{m}_{H_2,out3}}{\dot{m}_{H_2,out3,0}} \right)}
        """
        expr = self.outl[2].m.val_SI / self.outl[2].m.design

        return (self.P.val - self.outl[2].m.val_SI * self.e0 / (
                self.eta.design * self.eta_char.func.f_x(expr)))

    def eta_char_deriv(self):
        r"""
        Calculates the matrix of partial derivatives of the efficiency
        characteristic function.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """

        deriv = np.zeros((1, 5 + self.num_vars, self.num_fl + 3))

        deriv[0, 4, 0] = self.numeric_deriv(self.eta_char_func, 'm', 4)

        # derivatives for variable P
        if self.P.is_var:
            deriv[0, 5 + self.P.var_pos, 0] = 1

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

                val = \begin{cases}
                P \cdot f_{char}\left( \frac{P}{P_{ref}}\right) &
                \text{key = 'P'}\\
                \dot{Q} \cdot f_{char}\left( \frac{\dot{Q}}
                {\dot{Q}_{ref}}\right) & \text{key = 'Q'}\\
                \end{cases}\\
                \dot{Q} = - \dot{m}_{1,in} \cdot
                \left(h_{out,1} - h_{in,1} \right)\\
        """
        ######################################################################
        # equations for power on bus
        if bus.param == 'P':
            P = - self.energy_balance()
            if np.isnan(bus.P_ref):
                expr = 1
            else:
                expr = abs(P / bus.P_ref)
            return P * bus.char.f_x(expr)

        ######################################################################
        # equations for heat on bus

        elif bus.param == 'Q':
            val = - self.inl[0].m.val_SI * (self.outl[0].h.val_SI -
                                            self.inl[0].h.val_SI)
            if np.isnan(bus.P_ref):
                expr = 1
            else:
                expr = abs(val / bus.P_ref)
            return val * bus.char.f_x(expr)

        ######################################################################
        # missing/invalid bus parameter

        else:
            msg = ('The parameter ' + str(bus.param) + ' is not a valid '
                   'parameter for a component of type ' + self.component() +
                   '. Please specify a bus parameter (P/Q) for component ' +
                   self.label + '.')
            logging.error(msg)
            raise ValueError(msg)

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
        deriv = np.zeros((1, 5 + self.num_vars, self.num_fl + 3))

        ######################################################################
        # derivatives for power on bus
        if bus.param == 'P':
            deriv[0, 0, 0] = self.numeric_deriv(self.bus_func, 'm', 0, bus=bus)
            deriv[0, 0, 2] = self.numeric_deriv(self.bus_func, 'h', 0, bus=bus)

            deriv[0, 1, 0] = self.numeric_deriv(self.bus_func, 'm', 1, bus=bus)
            deriv[0, 1, 2] = self.numeric_deriv(self.bus_func, 'h', 1, bus=bus)

            deriv[0, 2, 2] = self.numeric_deriv(self.bus_func, 'h', 2, bus=bus)

            deriv[0, 3, 0] = self.numeric_deriv(self.bus_func, 'm', 3, bus=bus)
            deriv[0, 3, 2] = self.numeric_deriv(self.bus_func, 'h', 3, bus=bus)

            deriv[0, 4, 0] = self.numeric_deriv(self.bus_func, 'm', 4, bus=bus)
            deriv[0, 4, 2] = self.numeric_deriv(self.bus_func, 'h', 4, bus=bus)
            # variable power
            if self.P.is_var:
                deriv[0, 5 + self.P.var_pos, 0] = (
                        self.numeric_deriv(self.bus_func, 'P', 5, bus=bus))

        ######################################################################
        # derivatives for heat on bus
        elif bus.param == 'Q':

            deriv = np.zeros((1, 5 + self.num_vars, self.num_fl + 3))

            deriv[0, 0, 0] = self.numeric_deriv(self.bus_func, 'm', 0, bus=bus)
            deriv[0, 0, 2] = self.numeric_deriv(self.bus_func, 'h', 0, bus=bus)
            deriv[0, 2, 2] = self.numeric_deriv(self.bus_func, 'h', 2, bus=bus)

        ######################################################################
        # missing/invalid bus parameter

        else:
            msg = ('The parameter ' + str(bus.param) + ' is not a valid '
                   'parameter for a component of type ' + self.component() +
                   '. Please specify a bus parameter (P/Q) for component ' +
                   self.label + '.')
            logging.error(msg)
            raise ValueError(msg)

        return deriv

    def energy_balance(self):
        r"""
        Calculates the residual in energy balance of the adiabatic water
        electrolyzer. The residual is the negative to the necessary power
        input.

        Returns
        -------
        res : float
            Residual value.

            .. math::

                \begin{split}
                res = & \dot{m}_{in,2} \cdot \left( h_{in,2} - h_{in,2,ref}
                \right)\\ & - \dot{m}_{out,3} \cdot e_0\\
                & -\dot{m}_{in,1} \cdot \left( h_{out,1} - h_{in,1} \right)\\
                & - \dot{m}_{out,2} \cdot \left( h_{out,2} - h_{out,2,ref}
                \right)\\
                & - \dot{m}_{out,3} \cdot \left( h_{out,3} - h_{out,3,ref}
                \right)\\
                \end{split}

        Note
        ----
        The temperature for the reference state is set to 20 °C, thus
        the feed water must be liquid as proposed in the calculation of
        the minimum specific energy consumption for electrolysis:
        :func:`tespy.components.components.water_electrolyzer.calc_e0`.
        The part of the equation regarding the cooling water is implemented
        with negative sign as the energy for cooling is extracted from the
        reactor.

        - Reference temperature: 293.15 K.
        - Reference pressure: 1 bar.
        """
        T_ref = 293.15
        p_ref = 1e5

        # equations to set a reference point for each h2o, h2 and o2
        h_refh2o = h_mix_pT([1, p_ref, 0, self.inl[1].fluid.val], T_ref)
        h_refh2 = h_mix_pT([1, p_ref, 0, self.outl[2].fluid.val], T_ref)
        h_refo2 = h_mix_pT([1, p_ref, 0, self.outl[1].fluid.val], T_ref)

        val = (self.inl[1].m.val_SI * (self.inl[1].h.val_SI - h_refh2o) -
               self.outl[2].m.val_SI * self.e0 -
               self.inl[0].m.val_SI * (self.outl[0].h.val_SI -
                                       self.inl[0].h.val_SI) -
               self.outl[1].m.val_SI * (self.outl[1].h.val_SI - h_refo2) -
               self.outl[2].m.val_SI * (self.outl[2].h.val_SI - h_refh2))
        return val

    def initialise_fluids(self, nw):
        r"""
        Sets values to pure fluid on water inlet and gas outlets.

        Parameters
        ----------
        nw : tespy.networks.network
            Network using this component object.
        """
        self.outl[1].fluid.val[self.o2] = 1
        self.outl[2].fluid.val[self.h2] = 1
        self.inl[1].fluid.val[self.h2o] = 1

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
                5  \cdot 10^5 & \text{key = 'p'}\\
                h\left(T=323.15, p=5  \cdot 10^5\right) & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 5e5
        elif key == 'h':
            flow = [c.m.val0, 5e5, c.h.val_SI, c.fluid.val]
            T = 50 + 273.15
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
                5  \cdot 10^5 & \text{key = 'p'}\\
                h\left(T=293.15, p=5  \cdot 10^5\right) & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 5e5
        elif key == 'h':
            flow = [c.m.val0, 5e5, c.h.val_SI, c.fluid.val]
            T = 20 + 273.15
            return h_mix_pT(flow, T)

    def calc_parameters(self):
        r"""
        Postprocessing parameter calculation.
        """
        self.Q.val = - self.inl[0].m.val_SI * (self.outl[0].h.val_SI -
                                               self.inl[0].h.val_SI)
        self.pr_c.val = self.outl[0].p.val_SI / self.inl[0].p.val_SI
        self.e.val = self.P.val / self.outl[2].m.val_SI
        self.eta.val = self.e0 / self.e.val

        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()
        self.zeta.val = ((i[1] - o[1]) * np.pi ** 2 /
                         (8 * i[0] ** 2 * (v_mix_ph(i) + v_mix_ph(o)) / 2))

        if self.eta_char.is_set:
            # get bound errors for kA hot side characteristics
            expr = self.outl[2].m.val_SI / self.outl[2].m.design
            self.eta_char.func.get_bound_errors(expr, self.label)

        self.check_parameter_bounds()
