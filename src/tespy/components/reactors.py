# -*- coding: utf-8

"""Module for components of type reactor.

Components in this module:

- :func:`tespy.components.reactors.water_electrolyzer`


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/reactors.py

SPDX-License-Identifier: MIT
"""

import logging

import CoolProp.CoolProp as CP
import numpy as np

from tespy.components.components import component
from tespy.tools.data_containers import dc_cc
from tespy.tools.data_containers import dc_cp
from tespy.tools.data_containers import dc_simple
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.fluid_properties import dT_mix_dph
from tespy.tools.fluid_properties import dT_mix_pdh
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.fluid_properties import v_mix_ph
from tespy.tools.global_vars import err
from tespy.tools.global_vars import molar_masses
from tespy.tools.helpers import TESPyComponentError

# %%


class water_electrolyzer(component):
    r"""
    The water electrolyzer produces hydrogen and oxygen from water and power.

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

        Energy balance equation

        :func:`tespy.components.reactors.water_electrolyzer.energy_balance`.

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

        - :func:`tespy.components.reactors.water_electrolyzer.eta_char_func`

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

    P : float/tespy.helpers.dc_cp
        Power input, :math:`P/\text{W}`.

    Q : float/tespy.tools.data_containers.dc_cp
        Heat output of cooling, :math:`Q/\text{W}`

    e : str/float/tespy.tools.data_containers.dc_cp
        Electrolysis specific energy consumption,
        :math:`e/(\text{J}/\text{m}^3)`.

    eta : str/float/tespy.tools.data_containers.dc_cp
        Electrolysis efficiency, :math:`\eta/1`.

    eta_char : tespy.tools.characteristics.char_line/tespy.tools.data_containers.dc_cc
        Electrolysis efficiency characteristic line.

    pr : float/tespy.tools.data_containers.dc_cp
        Cooling loop pressure ratio, :math:`pr/1`.

    zeta : str/float/tespy.tools.data_containers.dc_cp
        Geometry independent friction coefficient for cooling loop pressure
        drop, :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    Note
    ----
    Other than usual components, the water electrolyzer has the fluid
    composition built into its equations for the feed water inlet and the
    hydrogen and oxygen outlet. Thus, the user must not specify the fluid
    composition at these connections!

    Example
    -------
    Create a water electrolyzer and compress the hydrogen, e.g. for a hydrogen
    storage.

    >>> from tespy.components import (sink, source, compressor,
    ... water_electrolyzer)
    >>> from tespy.connections import connection
    >>> from tespy.networks import network
    >>> from tespy.tools.data_containers import dc_cc
    >>> import shutil
    >>> fluid_list = ['O2', 'H2O', 'H2']
    >>> nw = network(fluids=fluid_list, T_unit='C', p_unit='bar',
    ... v_unit='l / s', iterinfo=False)
    >>> fw = source('feed water')
    >>> oxy = sink('oxygen sink')
    >>> hydro = sink('hydrogen sink')
    >>> cw_cold = source('cooling water source')
    >>> cw_hot = sink('cooling water sink')
    >>> comp = compressor('compressor', eta_s=0.9)
    >>> el = water_electrolyzer('electrolyzer')
    >>> el.component()
    'water electrolyzer'

    The electrolyzer should produce 100 l/s of hydrogen at an operating
    pressure of 10 bars and an outlet temperature of 50 °C. The fluid
    composition needs to be specified for the cooling liquid only. The storage
    pressure is 25 bars. The electrolysis efficiency is at 80 % and the
    compressor isentropic efficiency at 85 %. After designing the plant the
    offdesign electrolysis efficiency is predicted by the characteristic line.
    The default characteristic line can be found here: :py:mod:`tespy.data`.

    >>> fw_el = connection(fw, 'out1', el, 'in2')
    >>> el_o = connection(el, 'out2', oxy, 'in1')
    >>> el_cmp = connection(el, 'out3', comp, 'in1')
    >>> cmp_h = connection(comp, 'out1', hydro, 'in1')
    >>> cw_el = connection(cw_cold, 'out1', el, 'in1')
    >>> el_cw = connection(el, 'out1', cw_hot, 'in1')
    >>> nw.add_conns(fw_el, el_o, el_cmp, cmp_h, cw_el, el_cw)
    >>> fw_el.set_attr(p=10, T=15)
    >>> cw_el.set_attr(p=5, T=15, fluid={'H2O': 1, 'H2': 0, 'O2': 0})
    >>> el_cw.set_attr(T=45)
    >>> cmp_h.set_attr(p=25)
    >>> el_cmp.set_attr(v=100, T=50)
    >>> el.set_attr(eta=0.8, pr_c=0.99, design=['eta', 'pr_c'],
    ... offdesign=['eta_char', 'zeta'])
    >>> comp.set_attr(eta_s=0.85)
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(el.e0 / el.P.val * el_cmp.m.val_SI, 1)
    0.8
    >>> P_design = el.P.val / 1e6
    >>> round(P_design, 1)
    13.2
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(el.eta.val, 1)
    0.8
    >>> el_cmp.set_attr(v=np.nan)
    >>> el.set_attr(P=P_design * 0.66)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(el.eta.val, 2)
    0.88
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    @staticmethod
    def component():
        return 'water electrolyzer'

    @staticmethod
    def attr():
        return {
            'P': dc_cp(min_val=0), 'Q': dc_cp(max_val=0),
            'pr_c': dc_cp(max_val=1), 'zeta': dc_cp(min_val=0),
            'eta': dc_cp(min_val=0, max_val=1), 'eta_char': dc_cc(),
            'e': dc_cp(),
            'S': dc_simple()
        }

    @staticmethod
    def inlets():
        return ['in1', 'in2']

    @staticmethod
    def outlets():
        return ['out1', 'out2', 'out3']

    def comp_init(self, nw):

        if not self.P.is_set:
            self.set_attr(P='var')
            msg = ('The power output of cogeneration units must be set! '
                   'We are adding the power output of component ' +
                   self.label + ' as custom variable of the system.')
            logging.info(msg)

        component.comp_init(self, nw)

        o2 = [x for x in nw.fluids if x in [
            a.replace(' ', '') for a in CP.get_aliases('O2')]]
        if len(o2) == 0:
            msg = ('Missing oxygen in network fluids, component ' +
                   self.label + ' of type ' + self.component() +
                   ' requires oxygen in network fluids.')
            logging.error(msg)
            raise TESPyComponentError(msg)
        else:
            self.o2 = o2[0]

        h2o = [x for x in nw.fluids if x in [
            a.replace(' ', '') for a in CP.get_aliases('H2O')]]
        if len(h2o) == 0:
            msg = ('Missing water in network fluids, component ' +
                   self.label + ' of type ' + self.component() +
                   ' requires water in network fluids.')
            logging.error(msg)
            raise TESPyComponentError(msg)
        else:
            self.h2o = h2o[0]

        h2 = [x for x in nw.fluids if x in [
            a.replace(' ', '') for a in CP.get_aliases('H2')]]
        if len(h2) == 0:
            msg = ('Missing hydrogen in network fluids, component ' +
                   self.label + ' of type ' + self.component() +
                   ' requires hydrogen in network fluids.')
            logging.error(msg)
            raise TESPyComponentError(msg)
        else:
            self.h2 = h2[0]

        self.e0 = self.calc_e0()

        # number of mandatroy equations for
        # cooling loop fluids: num_fl
        # fluid composition at reactor inlets/outlets: 3 * num_fl
        # mass flow balances: 3
        # pressure: 2
        # energy balance: 1
        # reactor outlet temperatures: 1
        self.num_eq = self.num_nw_fluids * 4 + 7
        for var in [self.e, self.eta, self.eta_char, self.Q, self.pr_c,
                    self.zeta]:
            if var.is_set is True:
                self.num_eq += 1

        self.jacobian = np.zeros((
            self.num_eq,
            self.num_i + self.num_o + self.num_vars,
            self.num_nw_vars))

        self.residual = np.zeros(self.num_eq)
        pos = self.num_nw_fluids * 4
        self.jacobian[0:pos] = self.fluid_deriv()
        self.jacobian[pos:pos + 3] = self.mass_flow_deriv()
        self.jacobian[pos + 3:pos + 5] = self.pressure_deriv()

    def calc_e0(self):
        r"""
        Calculate the minimum specific energy required for electrolysis.

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
        r"""Calculate residual vector with results of equations."""
        k = 0
        ######################################################################
        # equations for fluids

        # equations for fluid composition in cooling water
        for fluid, x in self.inl[0].fluid.val.items():
            self.residual[k] = x - self.outl[0].fluid.val[fluid]
            k += 1

        # equations to constrain fluids to inlets/outlets
        self.residual[k] = 1 - self.inl[1].fluid.val[self.h2o]
        k += 1
        self.residual[k] = 1 - self.outl[1].fluid.val[self.o2]
        k += 1
        self.residual[k] = 1 - self.outl[2].fluid.val[self.h2]
        k += 1

        # equations to ban fluids off inlets/outlets
        for fluid in self.inl[1].fluid.val.keys():
            if fluid != self.h2o:
                self.residual[k] = 0 - self.inl[1].fluid.val[fluid]
                k += 1
            if fluid != self.o2:
                self.residual[k] = 0 - self.outl[1].fluid.val[fluid]
                k += 1
            if fluid != self.h2:
                self.residual[k] = 0 - self.outl[2].fluid.val[fluid]
                k += 1

        ######################################################################
        # eqations for mass flow balance
        # equation to calculate the ratio of o2 in water
        o2 = molar_masses[self.o2] / (
            molar_masses[self.o2] + 2 * molar_masses[self.h2])

        # equation for mass flow balance cooling water
        self.residual[k] = self.inl[0].m.val_SI - self.outl[0].m.val_SI
        k += 1

        # equations for mass flow balance electrolyzer
        self.residual[k] = o2 * self.inl[1].m.val_SI - self.outl[1].m.val_SI
        k += 1
        self.residual[k] = (
            (1 - o2) * self.inl[1].m.val_SI - self.outl[2].m.val_SI)
        k += 1

        ######################################################################
        # equations for pressure to set o2 and h2 output equal
        self.residual[k] = self.inl[1].p.val_SI - self.outl[1].p.val_SI
        k += 1
        self.residual[k] = self.inl[1].p.val_SI - self.outl[2].p.val_SI
        k += 1

        ######################################################################
        # equation for energy balance
        if np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0:
            self.residual[k] = self.P.val + self.energy_balance()
        k += 1

        ######################################################################
        # temperature electrolyzer outlet
        if np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0:
            self.residual[k] = (
                T_mix_ph(self.outl[1].to_flow()) -
                T_mix_ph(self.outl[2].to_flow()))
        k += 1

        ######################################################################
        # power vs hydrogen production
        if self.e.is_set:
            self.residual[k] = self.P.val - self.outl[2].m.val_SI * self.e.val
            k += 1

        ######################################################################
        # specified pressure ratio
        if self.pr_c.is_set:
            self.residual[k] = (
                self.inl[0].p.val_SI * self.pr_c.val - self.outl[0].p.val_SI)
            k += 1

        ######################################################################
        # specified zeta value
        if self.zeta.is_set:
            if np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0:
                self.residual[k] = self.zeta_func(zeta='zeta')
            k += 1

        # equation for heat transfer

        if self.Q.is_set:
            self.residual[k] = (
                self.Q.val - self.inl[0].m.val_SI *
                (self.inl[0].h.val_SI - self.outl[0].h.val_SI))
            k += 1

        ######################################################################
        # specified efficiency (efficiency definition: e0 / e)
        if self.eta.is_set:
            self.residual[k] = (
                self.P.val - self.outl[2].m.val_SI *
                self.e0 / self.eta.val)
            k += 1

        ######################################################################
        # specified characteristic line for efficiency
        if self.eta_char.is_set:
            self.residual[k] = self.eta_char_func()
            k += 1

    def derivatives(self, increment_filter):
        r"""Calculate partial derivatives for given equations."""
        ######################################################################
        # derivatives fluid, mass flow and reactor pressure are static
        k = self.num_nw_fluids * 4 + 5
        ######################################################################
        # derivatives for energy balance equations
        T_ref = 293.15
        p_ref = 1e5
        h_refh2o = h_mix_pT([1, p_ref, 0, self.inl[1].fluid.val], T_ref)
        h_refh2 = h_mix_pT([1, p_ref, 0, self.outl[2].fluid.val], T_ref)
        h_refo2 = h_mix_pT([1, p_ref, 0, self.outl[1].fluid.val], T_ref)

        # derivatives cooling water inlet
        self.jacobian[k, 0, 0] = -(
            self.outl[0].h.val_SI - self.inl[0].h.val_SI)
        self.jacobian[k, 0, 2] = self.inl[0].m.val_SI

        # derivatives feed water inlet
        self.jacobian[k, 1, 0] = (self.inl[1].h.val_SI - h_refh2o)
        self.jacobian[k, 1, 2] = self.inl[1].m.val_SI

        # derivative cooling water outlet
        self.jacobian[k, 2, 2] = - self.inl[0].m.val_SI

        # derivatives oxygen outlet
        self.jacobian[k, 3, 0] = - (self.outl[1].h.val_SI - h_refo2)
        self.jacobian[k, 3, 2] = - self.outl[1].m.val_SI

        # derivatives hydrogen outlet
        self.jacobian[k, 4, 0] = - self.e0 - (self.outl[2].h.val_SI - h_refh2)
        self.jacobian[k, 4, 2] = - self.outl[2].m.val_SI

        # derivatives for variable P
        if self.P.is_var:
            self.jacobian[k, 5 + self.P.var_pos, 0] = 1

        k += 1

        ######################################################################
        # derivatives for temperature at gas outlets

        # derivatives for outlet 1
        if not increment_filter[3, 1]:
            self.jacobian[k, 3, 1] = dT_mix_dph(self.outl[1].to_flow())
        if not increment_filter[3, 2]:
            self.jacobian[k, 3, 2] = dT_mix_pdh(self.outl[1].to_flow())

        # derivatives for outlet 2
        if not increment_filter[4, 1]:
            self.jacobian[k, 4, 1] = - dT_mix_dph(self.outl[2].to_flow())
        if not increment_filter[4, 2]:
            self.jacobian[k, 4, 2] = - dT_mix_pdh(self.outl[2].to_flow())

        k += 1

        ######################################################################
        # derivatives for power vs. hydrogen production
        if self.e.is_set:
            self.jacobian[k, 4, 0] = - self.e.val
            # derivatives for variable P
            if self.P.is_var:
                self.jacobian[k, 5 + self.P.var_pos, 0] = 1
            # derivatives for variable e
            if self.e.is_var:
                self.jacobian[k, 5 + self.e.var_pos, 0] = -(
                    self.outl[2].m.val_SI)
            k += 1

        ######################################################################
        # derivatives for pressure ratio
        if self.pr_c.is_set:
            self.jacobian[k, 0, 1] = self.pr_c.val
            self.jacobian[k, 2, 1] = - 1
            k += 1

        ######################################################################
        # derivatives for zeta value
        if self.zeta.is_set:
            f = self.zeta_func
            if not increment_filter[0, 0]:
                self.jacobian[k, 0, 0] = self.numeric_deriv(
                    f, 'm', 0, zeta='zeta')
            if not increment_filter[0, 1]:
                self.jacobian[k, 0, 1] = self.numeric_deriv(
                    f, 'p', 0, zeta='zeta')
            if not increment_filter[0, 2]:
                self.jacobian[k, 0, 2] = self.numeric_deriv(
                    f, 'h', 0, zeta='zeta')
            if not increment_filter[2, 1]:
                self.jacobian[k, 2, 1] = self.numeric_deriv(
                    f, 'p', 2, zeta='zeta')
            if not increment_filter[2, 2]:
                self.jacobian[k, 2, 2] = self.numeric_deriv(
                    f, 'h', 2, zeta='zeta')

            # derivatives for variable zeta
            if self.zeta.is_var:
                self.jacobian[k, 5 + self.zeta.var_pos, 0] = (
                    self.numeric_deriv(
                        f, 'zeta', 5, zeta='zeta'))
            k += 1

        ######################################################################
        # derivatives for heat flow
        if self.Q.is_set:
            self.jacobian[k, 0, 0] = -(
                self.inl[0].h.val_SI - self.outl[0].h.val_SI)
            self.jacobian[k, 0, 2] = - self.inl[0].m.val_SI
            self.jacobian[k, 2, 2] = self.inl[0].m.val_SI
            k += 1

        ######################################################################
        # specified efficiency (efficiency definition: e0 / e)
        if self.eta.is_set:
            self.jacobian[k, 4, 0] = - self.e0 / self.eta.val
            # derivatives for variable P
            if self.P.is_var:
                self.jacobian[k, 5 + self.P.var_pos, 0] = 1
            k += 1

        ######################################################################
        # specified characteristic line for efficiency
        if self.eta_char.is_set:
            self.jacobian[k, 4, 0] = self.numeric_deriv(
                self.eta_char_func, 'm', 4)
            # derivatives for variable P
            if self.P.is_var:
                self.jacobian[k, 5 + self.P.var_pos, 0] = 1
            k += 1

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
                self.eta.design * self.eta_char.func.evaluate(expr)))

    def fluid_deriv(self):
        r"""
        Calculate the partial derivatives for cooling loop fluid balance.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the fluid equations.
        """
        # derivatives for cooling liquid composition
        deriv = np.zeros((
            self.num_nw_fluids * 4,
            5 + self.num_vars,
            self.num_nw_vars))

        k = 0
        for fluid, x in self.inl[0].fluid.val.items():
            deriv[k, 0, 3 + k] = 1
            deriv[k, 2, 3 + k] = -1
            k += 1

        # derivatives to constrain fluids to inlets/outlets
        i = 0
        for fluid in self.nw_fluids:
            if fluid == self.h2o:
                deriv[k, 1, 3 + i] = -1
            elif fluid == self.o2:
                deriv[k + 1, 3, 3 + i] = -1
            elif fluid == self.h2:
                deriv[k + 2, 4, 3 + i] = -1
            i += 1
        k += 3

        # derivatives to ban fluids off inlets/outlets
        i = 0
        for fluid in self.nw_fluids:
            if fluid != self.h2o:
                deriv[k, 1, 3 + i] = -1
                k += 1
            if fluid != self.o2:
                deriv[k, 3, 3 + i] = -1
                k += 1
            if fluid != self.h2:
                deriv[k, 4, 3 + i] = -1
                k += 1
            i += 1

        return deriv

    def mass_flow_deriv(self):
        r"""
        Calculate the partial derivatives for all mass flow balance equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the mass flow equations.
        """
        # deritatives for mass flow balance in the heat exchanger
        deriv = np.zeros((3, 5 + self.num_vars, self.num_nw_vars))
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

        return deriv

    def pressure_deriv(self):
        r"""
        Calculate the partial derivatives for combustion pressure equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the pressure equations.
        """
        deriv = np.zeros((2, 5 + self.num_vars, self.num_nw_vars))
        # derivatives for pressure oxygen outlet
        deriv[0, 1, 1] = 1
        deriv[0, 3, 1] = -1
        # derivatives for pressure hydrogen outlet
        deriv[1, 1, 1] = 1
        deriv[1, 4, 1] = -1

        return deriv

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

                \dot{E} = \begin{cases}
                P & \text{key = 'P'}\\
                - \dot{m}_{1,in} \cdot \left(h_{out,1} - h_{in,1} \right) &
                \text{key = 'Q'}\\
                \end{cases}
        """
        ######################################################################
        # equations for power on bus
        if bus['param'] == 'P':
            val = - self.energy_balance()

        ######################################################################
        # equations for heat on bus

        elif bus['param'] == 'Q':
            val = - self.inl[0].m.val_SI * (
                self.outl[0].h.val_SI - self.inl[0].h.val_SI)

        ######################################################################
        # missing/invalid bus parameter

        else:
            msg = ('The parameter ' + str(bus['param']) + ' is not a valid '
                   'parameter for a component of type ' + self.component() +
                   '. Please specify a bus parameter (P/Q) for component ' +
                   self.label + '.')
            logging.error(msg)
            raise ValueError(msg)

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
        deriv : ndarray
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 5 + self.num_vars, self.num_nw_vars))
        f = self.calc_bus_value
        b = bus.comps.loc[self]

        ######################################################################
        # derivatives for power on bus
        if b['param'] == 'P':
            deriv[0, 0, 0] = self.numeric_deriv(f, 'm', 0, bus=bus)
            deriv[0, 0, 2] = self.numeric_deriv(f, 'h', 0, bus=bus)

            deriv[0, 1, 0] = self.numeric_deriv(f, 'm', 1, bus=bus)
            deriv[0, 1, 2] = self.numeric_deriv(f, 'h', 1, bus=bus)

            deriv[0, 2, 2] = self.numeric_deriv(f, 'h', 2, bus=bus)

            deriv[0, 3, 0] = self.numeric_deriv(f, 'm', 3, bus=bus)
            deriv[0, 3, 2] = self.numeric_deriv(f, 'h', 3, bus=bus)

            deriv[0, 4, 0] = self.numeric_deriv(f, 'm', 4, bus=bus)
            deriv[0, 4, 2] = self.numeric_deriv(f, 'h', 4, bus=bus)
            # variable power
            if self.P.is_var:
                deriv[0, 5 + self.P.var_pos, 0] = (
                        self.numeric_deriv(f, 'P', 5, bus=bus))

        ######################################################################
        # derivatives for heat on bus
        elif b['param'] == 'Q':

            deriv[0, 0, 0] = self.numeric_deriv(f, 'm', 0, bus=bus)
            deriv[0, 0, 2] = self.numeric_deriv(f, 'h', 0, bus=bus)
            deriv[0, 2, 2] = self.numeric_deriv(f, 'h', 2, bus=bus)

        ######################################################################
        # missing/invalid bus parameter

        else:
            msg = ('The parameter ' + str(b['param']) + ' is not a valid '
                   'parameter for a component of type ' + self.component() +
                   '. Please specify a bus parameter (P/Q) for component ' +
                   self.label + '.')
            logging.error(msg)
            raise ValueError(msg)

        return deriv

    def energy_balance(self):
        r"""
        Calculate the residual in energy balance of the water electrolyzer.

        The residual is the negative to the necessary power input.

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
        :func:`tespy.components.reactors.water_electrolyzer.calc_e0`.
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
        Set values to pure fluid on water inlet and gas outlets.

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
                5  \cdot 10^5 & \text{key = 'p'}\\
                h\left(T=323.15, p=5  \cdot 10^5\right) & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 5e5
        elif key == 'h':
            flow = c.to_flow()
            T = 50 + 273.15
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
                5  \cdot 10^5 & \text{key = 'p'}\\
                h\left(T=293.15, p=5  \cdot 10^5\right) & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 5e5
        elif key == 'h':
            flow = c.to_flow()
            T = 20 + 273.15
            return h_mix_pT(flow, T)

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
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
            # get bound errors for efficiency characteristics
            expr = self.outl[2].m.val_SI / self.outl[2].m.design
            self.eta_char.func.get_bound_errors(expr, self.label)

        self.check_parameter_bounds()
