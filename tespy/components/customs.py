# -*- coding: utf-8

"""Module for custom components.

Components in this module:

    - :func:`tespy.components.customs.orc_evaporator`

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/customs.py

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


class orc_evaporator(component):
    r"""
    Class orc_evaporator is the evaporator component in the Organic Rankine Cycle (ORC).
    Generally, the hot side of the geo-fluid from the geothermal wells keeps 2-phase.
    In order to fully use the energy in the geo-fluid, there are 2 inlets at the hot side.

    The ORC evaporator represents counter current evaporators. Both, 2 hot
    and 1 cold side of the evaporator, are simulated.

    Equations

        **mandatory equations**

        - :func:`tespy.components.customs.orc_evaporator.fluid_func`
        - :func:`tespy.components.customs.orc_evaporator.mass_flow_func`

        - :func:`tespy.components.customs.orc_evaporator.energy_func`

        **optional equations**

        .. math::

            0 = \dot{m}_{in3} \cdot \left(h_{out3} - h_{in3} \right) - \dot{Q}

        - :func:`tespy.components.customs.orc_evaporator.kA_func`

        .. math::

            0 = p_{1,in} \cdot pr1 - p_{1,out}\\
            0 = p_{2,in} \cdot pr2 - p_{2,out}\\
            0 = p_{3,in} \cdot pr3 - p_{3,out}

        - :func:`tespy.components.customs.orc_evaporator.zeta_func`
        - :func:`tespy.components.customs.orc_evaporator.zeta2_func`
        - :func:`tespy.components.customs.orc_evaporator.zeta3_func`

        **additional equations**

        - :func:`tespy.components.customs.orc_evaporator.additional_equations`

    Inlets/Outlets

        - in1, in2, in3 (index 1: hot side 1, index 2: hot side 2, index 3: cold side)
        - out1, out2, out3 (index 1: hot side 1, index 2: hot side 2, index 3: cold side)

    Image

        .. image:: _images/orc_evaporator.svg
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

    Q : String/float/tespy.helpers.dc_cp
        Heat transfer, :math:`Q/\text{W}`.

    pr1 : String/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio at hot side 1, :math:`pr/1`.

    pr2 : String/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio at hot side 2, :math:`pr/1`.

    pr3 : String/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio at cold side, :math:`pr/1`.

    zeta1 : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient at hot side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    zeta2 : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient at cold side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    zeta3 : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient at cold side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    kA : str/float/tespy.helpers.dc_cp
        Area independent heat transition coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char1 : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient at hot side 1,
        provide x and y values or use generic values (e. g. calculated
        from design case). Standard method 'HE_HOT', Parameter 'm'.

    kA_char2 : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient at hot side 2,
        provide x and y values or use generic values (e. g. calculated from
        design case). Standard method 'HE_HOT', Parameter 'm'.

    kA_char3 : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient at cold side,
        provide x and y values or use generic values (e. g. calculated from
        design case). Standard method 'HE_COLD', Parameter 'm'.

    Note
    ----
    The ORC evaporator are countercurrent heat exchangers. Equation kA do not work
    for directcurrent and crosscurrent or combinations of different types.

    Example
    -------
    A 2-phase geo-fluid is used as the heat source for evaporating the working fluid.
    The evaporator is designed for calculate the mass flow rate of the working fluid
    with known steam and brine mass flow rate. The state of the brine, steam and is fixed.
    From this, it is possible to calculate the mass flow rate of the working fluid that is
    fully evaporated through the ORC evaporator and its heat transfer coefficient.

    >>> from tespy.connections import connection
    >>> from tespy.networks import network
    >>> from tespy.components import source, sink
    >>> from tespy.components.customs import orc_evaporator
    >>> from tespy.tools import logger
    >>> import logging
    >>> mypath = logger.define_logging(
    ... log_path=True, log_version=True, timed_rotating={'backupCount': 4},
    ... screen_level=logging.WARNING, screen_datefmt = "no_date")
    >>> fluids = ['water', 'Isopentane']
    >>> nw = network(fluids=fluids)
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')
    >>> evaporator = orc_evaporator('orc_evaporator')
    >>> source_wf = source('working fluid source')
    >>> sink_wf = sink('working fluid sink')
    >>> source_s = source('steam source')
    >>> source_b = source('brine source')
    >>> sink_s = sink('steam sink')
    >>> sink_b = sink('brine sink')
    >>> evaporator_wf_in = connection(source_wf, 'out1', evaporator, 'in3')
    >>> evaporator_wf_out = connection(evaporator, 'out3', sink_wf, 'in1')
    >>> evaporator_steam_in = connection(source_s, 'out1', evaporator, 'in1')
    >>> evaporator_sink_s = connection(evaporator, 'out1', sink_s, 'in1')
    >>> evaporator_brine_in = connection(source_b, 'out1', evaporator, 'in2')
    >>> evaporator_sink_b = connection(evaporator, 'out2', sink_b, 'in1')
    >>> nw.add_conns(evaporator_wf_in, evaporator_wf_out)
    >>> nw.add_conns(evaporator_steam_in, evaporator_sink_s)
    >>> nw.add_conns(evaporator_brine_in, evaporator_sink_b)
    >>> evaporator.set_attr(pr1=0.93181818, pr2=0.970588, pr3=1)
    >>> evaporator_wf_in.set_attr(T=111.6, p=10.8, fluid={'water': 0, 'Isopentane': 1})
    >>> evaporator_steam_in.set_attr(T=146.6, p=4.34, m=20.4, state='g', fluid={'water': 1, 'Isopentane': 0})
    >>> evaporator_brine_in.set_attr(T=146.6, p=10.2, m=190.8, fluid={'water': 1, 'Isopentane': 0})
    >>> evaporator_sink_b.set_attr(T=118.6)
    >>> mode = 'design'
    >>> file = 'orc_evaporator'
    >>> nw.solve(mode=mode)
    >>> nw.print_results()
    >>> nw.save(file)
    """

    @staticmethod
    def component():
        return 'orc_evaporator'

    @staticmethod
    def attr():
        return {'Q': dc_cp(max_val=0),
                'kA': dc_cp(min_val=0),
                'td_log': dc_cp(min_val=0),
                'pr1': dc_cp(max_val=1), 'pr2': dc_cp(max_val=1), 'pr3': dc_cp(max_val=1),
                'zeta1': dc_cp(min_val=0), 'zeta2': dc_cp(min_val=0), 'zeta3': dc_cp(min_val=0),
                'subcooling': dc_simple(val=False), 'overheating': dc_simple(val=False),
                'kA_char1': dc_cc(param='m'),
                'kA_char2': dc_cc(param='m'),
                'kA_char3': dc_cc(param='m'),
                'SQ1': dc_simple(), 'SQ2': dc_simple(), 'SQ3': dc_simple(), 'Sirr': dc_simple(),
                'zero_flag': dc_simple()}

    @staticmethod
    def inlets():
        return ['in1', 'in2', 'in3']

    @staticmethod
    def outlets():
        return ['out1', 'out2', 'out3']

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
            vec_res += [self.inl[3].m.val_SI *
                        (self.outl[3].h.val_SI - self.inl[3].h.val_SI) +
                        self.Q.val]

        ######################################################################
        # equations for specified heat transfer coefficient
        if self.kA.is_set:
            vec_res += [self.kA_func()]

        ######################################################################
        # equations for specified pressure ratio at hot side 1
        if self.pr1.is_set:
            vec_res += [self.pr1.val * self.inl[0].p.val_SI -
                        self.outl[0].p.val_SI]

        ######################################################################
        # equations for specified pressure ratio at hot side 2
        if self.pr2.is_set:
            vec_res += [self.pr2.val * self.inl[1].p.val_SI -
                        self.outl[1].p.val_SI]

        ######################################################################
        # equations for specified pressure ratio at cold side
        if self.pr3.is_set:
            vec_res += [self.pr3.val * self.inl[2].p.val_SI -
                        self.outl[2].p.val_SI]

        ######################################################################
        # equations for specified zeta at hot side 1
        if self.zeta1.is_set:
            vec_res += [self.zeta_func()]

        ######################################################################
        # equations for specified zeta at hot side 2
        if self.zeta2.is_set:
            vec_res += [self.zeta2_func()]

        ######################################################################
        # equations for specified zeta at cold side
        if self.zeta3.is_set:
            vec_res += [self.zeta3_func()]

        ######################################################################
        # additional equations
        vec_res += self.additional_equations()

        return vec_res

    def additional_equations(self):
        r"""
        Calculates vector vec_res with results of additional equations for this
        component.

        Equations

            **mandatory equations at outlet 1 of the hot side**

            .. math::

                0 = h_{1,out} - h\left(p, x=0 \right)\\
                x: \text{vapour mass fraction}

            **mandatory equations at outlet of the cold side**

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
        # equation for saturated liquid at hot side 1 outlet
        if not self.subcooling.val:
            outl = self.outl
            o1 = outl[0].to_flow()
            vec_res += [o1[2] - h_mix_pQ(o1, 0)]

        ######################################################################
        # equation for saturated gas at cold side outlet
        if not self.overheating.val:
            outl = self.outl
            o3 = outl[2].to_flow()
            vec_res += [o3[2] - h_mix_pQ(o3, 1)]

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
        # derivatives for fluid balance equations
        mat_deriv += self.fl_deriv
        ######################################################################
        # derivatives for mass flow balance equations
        mat_deriv += self.m_deriv

        ######################################################################
        # derivatives for energy balance equation
        mat_deriv += self.energy_deriv()

        ######################################################################
        # derivatives for specified pressure ratio at hot side 1
        if self.pr1.is_set:
            pr1_deriv = np.zeros((1, 6, self.num_fl + 3))
            pr1_deriv[0, 0, 1] = self.pr1.val
            pr1_deriv[0, 3, 1] = -1
            mat_deriv += pr1_deriv.tolist()

        ######################################################################
        # derivatives for specified pressure ratio at hot side 2
        if self.pr2.is_set:
            pr2_deriv = np.zeros((1, 6, self.num_fl + 3))
            pr2_deriv[0, 1, 1] = self.pr2.val
            pr2_deriv[0, 4, 1] = -1
            mat_deriv += pr2_deriv.tolist()

        ######################################################################
        # derivatives for specified pressure ratio at cold side
        if self.pr3.is_set:
            pr3_deriv = np.zeros((1, 6, self.num_fl + 3))
            pr3_deriv[0, 2, 1] = self.pr3.val
            pr3_deriv[0, 5, 1] = -1
            mat_deriv += pr3_deriv.tolist()

        ######################################################################
        # derivatives for specified zeta at hot side 1
        if self.zeta1.is_set:
            zeta1_deriv = np.zeros((1, 6, self.num_fl + 3))
            zeta1_deriv[0, 0, 0] = self.numeric_deriv(self.zeta_func, 'm', 0)
            zeta1_deriv[0, 0, 1] = self.numeric_deriv(self.zeta_func, 'p', 0)
            zeta1_deriv[0, 0, 2] = self.numeric_deriv(self.zeta_func, 'h', 0)
            zeta1_deriv[0, 3, 1] = self.numeric_deriv(self.zeta_func, 'p', 3)
            zeta1_deriv[0, 3, 2] = self.numeric_deriv(self.zeta_func, 'h', 3)
            mat_deriv += zeta1_deriv.tolist()

        ######################################################################
        # derivatives for specified zeta at hot side 2
        if self.zeta2.is_set:
            zeta2_deriv = np.zeros((1, 6, self.num_fl + 3))
            zeta2_deriv[0, 1, 0] = self.numeric_deriv(self.zeta2_func, 'm', 1)
            zeta2_deriv[0, 1, 1] = self.numeric_deriv(self.zeta2_func, 'p', 1)
            zeta2_deriv[0, 1, 2] = self.numeric_deriv(self.zeta2_func, 'h', 1)
            zeta2_deriv[0, 4, 1] = self.numeric_deriv(self.zeta2_func, 'p', 4)
            zeta2_deriv[0, 4, 2] = self.numeric_deriv(self.zeta2_func, 'h', 4)
            mat_deriv += zeta2_deriv.tolist()

        ######################################################################
        # derivatives for specified zeta at cold side
        if self.zeta3.is_set:
            zeta3_deriv = np.zeros((1, 6, self.num_fl + 3))
            zeta3_deriv[0, 2, 0] = self.numeric_deriv(self.zeta3_func, 'm', 2)
            zeta3_deriv[0, 2, 1] = self.numeric_deriv(self.zeta3_func, 'p', 2)
            zeta3_deriv[0, 2, 2] = self.numeric_deriv(self.zeta3_func, 'h', 2)
            zeta3_deriv[0, 5, 1] = self.numeric_deriv(self.zeta3_func, 'p', 5)
            zeta3_deriv[0, 5, 2] = self.numeric_deriv(self.zeta3_func, 'h', 5)
            mat_deriv += zeta3_deriv.tolist()

        ######################################################################
        # derivatives for additional equations
        mat_deriv += self.additional_derivatives()

        return np.asarray(mat_deriv)

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
        # derivatives for saturated liquid at hot side 1 outlet equation
        if not self.subcooling.val:
            o1 = self.outl[0].to_flow()
            x_deriv = np.zeros((1, 6, self.num_fl + 3))
            x_deriv[0, 3, 1] = -dh_mix_dpQ(o1, 0)
            x_deriv[0, 3, 2] = 1
            mat_deriv += x_deriv.tolist()

        ######################################################################
        # derivatives for saturated gas at cold side outlet 3 equation
        if not self.overheating.val:
            o3 = self.outl[2].to_flow()
            deriv = np.zeros((1, 6, self.num_fl + 3))
            deriv[0, 5, 1] = -dh_mix_dpQ(o3, 1)
            deriv[0, 5, 2] = 1
            mat_deriv += deriv.tolist()

        return mat_deriv

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
        deriv = np.zeros((self.num_fl * 3, 6 + self.num_vars, 3 + self.num_fl))
        # hot side 1
        i = 0
        for fluid in self.fluids:
            deriv[i, 0, i + 3] = 1
            deriv[i, 3, i + 3] = -1
            i += 1
        # hot side 2
        j = 0
        for fluid in self.fluids:
            deriv[i + j, 1, j + 3] = 1
            deriv[i + j, 4, j + 3] = -1
            j += 1
        # cold side
        k = 0
        for fluid in self.fluids:
            deriv[i + j + k, 2, k + 3] = 1
            deriv[i + j + k, 5, k + 3] = -1
            k += 1
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
        deriv = np.zeros((3, 6 + self.num_vars, self.num_fl + 3))
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
                \dot{m}_{2,in} \cdot \left(h_{2,out} - h_{2,in} \right) +
                \dot{m}_{3,in} \cdot \left(h_{3,out} - h_{3,in} \right)
        """
#        if self.zero_flag.is_set:
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
                                        self.inl[1].h.val_SI) +
                self.inl[2].m.val_SI * (self.outl[2].h.val_SI -
                                        self.inl[2].h.val_SI))

    def energy_deriv(self):
        r"""
        Calculates the matrix of partial derivatives for energy balance
        equation.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 6, len(self.inl[0].fluid.val) + 3))

#        if self.zero_flag.is_set:
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
        for k in range(3):
            deriv[0, k, 0] = self.outl[k].h.val_SI - self.inl[k].h.val_SI
            deriv[0, k, 2] = -self.inl[k].m.val_SI

        deriv[0, 3, 2] = self.inl[0].m.val_SI
        deriv[0, 4, 2] = self.inl[1].m.val_SI
        deriv[0, 5, 2] = self.inl[2].m.val_SI
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

                res = \dot{m}_{3,in} \cdot \left( h_{3,out} - h_{3,in}\right) +
                kA \cdot f_{kA} \cdot \frac{T_{1,out} + T_{2,out} - T_{3,in} -
                T_{2,in} - T_{1,in} + T_{3,out}}
                {\ln{\frac{T_{1,out} + T_{2,out} - T_{3,in}}
                {T_{1,in} + T_{2,in} - T_{3,out}}}}

                f_{kA} = f_1\left(\frac{m_1}{m_{1,ref}}\right) \cdot
                f_2\left(\frac{m_2}{m_{2,ref}} \cdot
                f_3\left(\frac{m_3}{m_{3,ref}}\right)

        Note
        ----
        For standard functions f\ :subscript:`1` \, f\ :subscript:`2` \
        and f\ :subscript:`2` \ see module :func:`tespy.data`.

        - Calculate temperatures at inlets and outlets.
        - Perform value manipulation, if temperature levels are not physically
          feasible.
        """

#        if self.zero_flag.is_set:
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
        i3 = self.inl[2].to_flow()
        o1 = self.outl[0].to_flow()
        o2 = self.outl[1].to_flow()
        o3 = self.outl[2].to_flow()

        i1_d = self.inl[0].to_flow_design()
        i2_d = self.inl[1].to_flow_design()
        i3_d = self.inl[2].to_flow_design()

        T_i1 = T_mix_ph(i1, T0=self.inl[0].T.val_SI)
        T_i2 = T_mix_ph(i2, T0=self.inl[1].T.val_SI)
        T_i3 = T_mix_ph(i3, T0=self.inl[2].T.val_SI)
        T_o1 = T_mix_ph(o1, T0=self.outl[0].T.val_SI)
        T_o2 = T_mix_ph(o2, T0=self.outl[1].T.val_SI)
        T_o3 = T_mix_ph(o3, T0=self.outl[2].T.val_SI)

#        if T_i1 <= T_o2 and self.inl[0].T.val_set is False:
        if T_i1 <= T_o3:
            T_i1 = T_o3 + 0.01
#        if T_i1 <= T_o2 and self.outl[1].T.val_set is False:
        if T_i1 <= T_o3:
            T_o3 = T_i1 - 0.01
#        if T_i1 < T_o2 and self.inl[0].T.val_set and self.outl[1].T.val_set:
#            msg = ('Infeasibility at ' + str(self.label) + ': Value for upper '
#                   'temperature difference is ' + str(round(T_i1 - T_o2)) +
#                   '.')
#            logging.error(msg)
#            raise ValueError(msg)

#        if T_i1 <= T_o2 and self.outl[1].T.val_set is False:
        if T_i1 <= T_o3:
            T_o1 = T_i3 + 0.02

#        if T_o1 <= T_i2 and self.inl[1].T.val_set is False:
        if T_o1 <= T_i3:
            T_i3 = T_o1 - 0.02
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
                    fkA1 = self.kA_char1.func.evaluate(i1[0] / i1_d[0])

        fkA2 = 1
        if self.kA_char2.param == 'm':
            if not np.isnan(i2_d[0]):
                if not i2[0] == 0:
                    fkA2 = self.kA_char2.func.evaluate(i2[0] / i2_d[0])

        fkA3 = 1
        if self.kA_char3.param == 'm':
            if not np.isnan(i3_d[0]):
                if not i3[0] == 0:
                    fkA3 = self.kA_char3.func.evaluate(i3[0] / i3_d[0])

        td_log = ((T_o1 + T_o2 - T_i3 - T_i1 - T_i2 + T_o3) /
                  np.log((T_o1 + T_o2 - T_i3) / (T_i1 + T_i2 - T_o3)))
        return i3[0] * (o3[2] - i3[2]) + self.kA.val * fkA1 * fkA2 * fkA3 * td_log

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

                P = \dot{m}_{3,in} \cdot \left( h_{3,out} - h_{3,in} \right)
        """
        i = self.inl[3].to_flow()
        o = self.outl[3].to_flow()

        val = i[0] * (o[2] - i[2])
        if np.isnan(bus.P_ref):
            expr = 1
        else:
            expr = abs(val / bus.P_ref)
        return val * bus.char.evaluate(expr)

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
        deriv = np.zeros((1, 6, len(self.inl[0].fluid.val) + 3))
        deriv[0, 0, 0] = self.numeric_deriv(self.bus_func, 'm', 0, bus=bus)
        deriv[0, 0, 2] = self.numeric_deriv(self.bus_func, 'h', 0, bus=bus)
        deriv[0, 3, 2] = self.numeric_deriv(self.bus_func, 'h', 3, bus=bus)
        return deriv

    def convergence_check(self, nw):
        r"""
        Performs a convergence check.
        """

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
            elif c.s_id == 'out2':
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
            elif c.t_id == 'in2':
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
        i3 = self.inl[2].to_flow()
        o1 = self.outl[0].to_flow()
        o2 = self.outl[1].to_flow()
        o3 = self.outl[2].to_flow()

        # temperatures
        if isinstance(self, orc_evaporator):
            T_i1 = T_bp_p(i1)
        else:
            T_i1 = T_mix_ph(i1, T0=self.inl[0].T.val_SI)
        T_i2 = T_mix_ph(i2, T0=self.inl[1].T.val_SI)
        T_i3 = T_mix_ph(i3, T0=self.inl[2].T.val_SI)
        T_o1 = T_mix_ph(o1, T0=self.outl[0].T.val_SI)
        T_o2 = T_mix_ph(o2, T0=self.outl[1].T.val_SI)
        T_o3 = T_mix_ph(o3, T0=self.outl[2].T.val_SI)

        # specific volume
        v_i1 = v_mix_ph(i1, T0=T_i1)
        v_i2 = v_mix_ph(i2, T0=T_i2)
        v_i3 = v_mix_ph(i3, T0=T_i3)
        v_o1 = v_mix_ph(o1, T0=T_o1)
        v_o2 = v_mix_ph(o2, T0=T_o2)
        v_o3 = v_mix_ph(o3, T0=T_o3)

        # specific entropy
        s_i1 = s_mix_ph(i1, T0=T_i1)
        s_i2 = s_mix_ph(i2, T0=T_i2)
        s_i3 = s_mix_ph(i3, T0=T_i3)
        s_o1 = s_mix_ph(o1, T0=T_o1)
        s_o2 = s_mix_ph(o2, T0=T_o2)
        s_o3 = s_mix_ph(o3, T0=T_o3)

        # component parameters
        self.Q.val = -i3[0] * (o3[2] - i3[2])

        self.pr1.val = o1[1] / i1[1]
        self.pr2.val = o2[1] / i2[1]
        self.pr3.val = o2[2] / i2[2]
        self.zeta1.val = ((i1[1] - o1[1]) * np.pi ** 2 /
                          (8 * i1[0] ** 2 * (v_i1 + v_o1) / 2))
        self.zeta2.val = ((i2[1] - o2[1]) * np.pi ** 2 /
                          (8 * i2[0] ** 2 * (v_i2 + v_o2) / 2))
        self.zeta3.val = ((i3[1] - o3[1]) * np.pi ** 2 /
                          (8 * i3[0] ** 2 * (v_i3 + v_o3) / 2))

        self.SQ1.val = self.inl[0].m.val_SI * (s_o1 - s_i1)
        self.SQ2.val = self.inl[1].m.val_SI * (s_o2 - s_i2)
        self.SQ3.val = self.inl[2].m.val_SI * (s_o3 - s_i3)
        self.Sirr.val = self.SQ1.val + self.SQ2.val + self.SQ3.val

        # kA and logarithmic temperature difference
        if T_i1 <= T_o3 or T_o1 <= T_i3:
            self.td_log.val = np.nan
            self.kA.val = np.nan
        else:
            self.td_log.val = ((T_o1 - T_i3 - T_i1 + T_o3) /
                               np.log((T_o1 - T_i3) / (T_i1 - T_o3)))
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

            # get bound errors for kA copld side characteristics
            if self.kA_char3.param == 'm':
                i3_d = self.inl[2].to_flow_design()
                if not np.isnan(i3_d[0]):
                    if not i1[0] == 0:
                        self.kA_char3.func.get_bound_errors(i3[0] / i3_d[0],
                                                            self.label)

        self.check_parameter_bounds()
