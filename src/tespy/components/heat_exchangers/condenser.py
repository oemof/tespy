# -*- coding: utf-8

"""Module of class Condenser.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/heat_exchangers/condenser.py

SPDX-License-Identifier: MIT
"""

import warnings

import numpy as np

from tespy.components.component import Component
from tespy.components.heat_exchangers.heat_exchanger import HeatExchanger
from tespy.tools.data_containers import dc_cc
from tespy.tools.data_containers import dc_cp
from tespy.tools.data_containers import dc_simple
from tespy.tools.fluid_properties import T_bp_p
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.fluid_properties import dh_mix_dpQ
from tespy.tools.fluid_properties import h_mix_pQ
from tespy.tools.fluid_properties import s_mix_ph
from tespy.tools.fluid_properties import v_mix_ph


class Condenser(HeatExchanger):
    r"""
    A Condenser cools a fluid until it is in liquid state.

    The condensing fluid is cooled by the cold side fluid. The fluid on the hot
    side of the condenser must be pure. Subcooling is available.

    Equations

        **mandatory equations**

        - :py:meth:`tespy.components.component.Component.fluid_func`
        - :py:meth:`tespy.components.heat_exchangers.HeatExchanger.mass_flow_func`
        - :py:meth:`tespy.components.heat_exchangers.HeatExchanger.energy_func`

        **optional equations**

        .. math::

            0 = \dot{m}_{in} \cdot \left(h_{out} - h_{in} \right) - \dot{Q}

        - :py:meth:`tespy.components.heat_exchangers.Condenser.kA_func`
        - :py:meth:`tespy.components.heat_exchangers.Condenser.kA_char_func`
        - :py:meth:`tespy.components.heat_exchangers.Condenser.ttd_u_func`

        .. math::

            0 = p_{1,in} \cdot pr1 - p_{1,out}\\
            0 = p_{2,in} \cdot pr2 - p_{2,out}

        - hot side :py:meth:`tespy.components.component.Component.zeta_func`
        - cold side :py:meth:`tespy.components.component.Component.zeta_func`

        **additional equations**

        - :py:meth:`tespy.components.heat_exchangers.Condenser.additional_equations`

    Inlets/Outlets

        - in1, in2 (index 1: hot side, index 2: cold side)
        - out1, out2 (index 1: hot side, index 2: cold side)

    Image

        .. image:: _images/Condenser.svg
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

    Q : str, float, tespy.tools.data_containers.dc_cp
        Heat transfer, :math:`Q/\text{W}`.

    pr1 : str, float, tespy.tools.data_containers.dc_cp
        Outlet to inlet pressure ratio at hot side, :math:`pr/1`.

    pr2 : str, float, tespy.tools.data_containers.dc_cp
        Outlet to inlet pressure ratio at cold side, :math:`pr/1`.

    zeta1 : str, float, tespy.tools.data_containers.dc_cp
        Geometry independent friction coefficient at hot side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    zeta2 : str, float, tespy.tools.data_containers.dc_cp
        Geometry independent friction coefficient at cold side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    kA : float, tespy.tools.data_containers.dc_cp
        Area independent heat transition coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char : tespy.tools.data_containers.dc_simple
        Area independent heat transition coefficient characteristic.

    kA_char1 : tespy.tools.charactersitics.CharLine, tespy.tools.data_containers.dc_cc
        Characteristic line for hot side heat transfer coefficient.

    kA_char2 : tespy.tools.charactersitics.CharLine, tespy.tools.data_containers.dc_cc
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

    >>> from tespy.components import Sink, Source, Condenser
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> from tespy.tools.fluid_properties import T_bp_p
    >>> import shutil
    >>> nw = Network(fluids=['water', 'air'], T_unit='C', p_unit='bar',
    ... h_unit='kJ / kg', m_range=[0.01, 1000], iterinfo=False)
    >>> amb_in = Sink('ambient air inlet')
    >>> amb_out = Source('air outlet')
    >>> waste_steam = Source('waste steam')
    >>> c = Sink('condensate sink')
    >>> cond = Condenser('condenser')
    >>> cond.component()
    'condenser'
    >>> amb_he = Connection(amb_out, 'out1', cond, 'in2')
    >>> he_amb = Connection(cond, 'out2', amb_in, 'in1')
    >>> ws_he = Connection(waste_steam, 'out1', cond, 'in1')
    >>> he_c = Connection(cond, 'out1', c, 'in1')
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

        Component.comp_init(self, nw)

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
        module :py:mod:`tespy.data`.

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

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        # connection information
        i1 = self.inl[0].to_flow()
        i2 = self.inl[1].to_flow()
        o1 = self.outl[0].to_flow()
        o2 = self.outl[1].to_flow()

        # temperatures
        T_i1 = T_bp_p(i1)
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
