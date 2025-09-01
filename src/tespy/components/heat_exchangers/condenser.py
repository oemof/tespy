# -*- coding: utf-8

"""Module of class Condenser.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/heat_exchangers/condenser.py

SPDX-License-Identifier: MIT
"""

import math

import numpy as np

from tespy.components.component import component_registry
from tespy.components.heat_exchangers.base import HeatExchanger
from tespy.tools import logger
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
from tespy.tools.fluid_properties import h_mix_pQ
from tespy.tools.fluid_properties import single_fluid


@component_registry
class Condenser(HeatExchanger):
    r"""
    A Condenser cools a fluid until it is in liquid state.

    The condensing fluid is cooled by the cold side fluid. The fluid on the hot
    side of the condenser must be pure. Subcooling is available.

    **Mandatory Equations**

    - fluid: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`
    - mass flow: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_func`
    - condensate outlet state, function can be disabled by specifying
      :code:`set_attr(subcooling=True)`
      :py:meth:`tespy.components.heat_exchangers.condenser.Condenser.subcooling_func`

    **Optional Equations**

    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_hot_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.kA_func`
       (utilizes the :code:`td_log` calculation of the Condenser class)
    - :py:meth:`tespy.components.heat_exchangers.condenser.Condenser.kA_char_func`
    - :py:meth:`tespy.components.heat_exchangers.condenser.Condenser.ttd_u_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.ttd_l_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.ttd_min_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.eff_cold_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.eff_hot_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.eff_max_func`
    - hot side :py:meth:`tespy.components.component.Component.pr_func`
    - cold side :py:meth:`tespy.components.component.Component.pr_func`
    - hot side :py:meth:`tespy.components.component.Component.zeta_func`
    - cold side :py:meth:`tespy.components.component.Component.zeta_func`
    - hot side :py:meth:`tespy.components.component.Component.dp_func`
    - cold side :py:meth:`tespy.components.component.Component.dp_func`

    Inlets/Outlets

    - in1, in2 (index 1: hot side, index 2: cold side)
    - out1, out2 (index 1: hot side, index 2: cold side)

    Image

    .. image:: /api/_images/Condenser.svg
       :alt: flowsheet of the condenser
       :align: center
       :class: only-light

    .. image:: /api/_images/Condenser_darkmode.svg
       :alt: flowsheet of the condenser
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

    Q : float, dict
        Heat transfer, :math:`Q/\text{W}`.

    pr1 : float, dict, :code:`"var"`
        Outlet to inlet pressure ratio at hot side, :math:`pr/1`.

    pr2 : float, dict, :code:`"var"`
        Outlet to inlet pressure ratio at cold side, :math:`pr/1`.

    dp1 : float, dict, :code:`"var"`
        Inlet to outlet pressure delta at hot side, unit is the network's
        pressure unit!.

    dp2 : float, dict, :code:`"var"`
        Inlet to outlet pressure delta at cold side, unit is the network's
        pressure unit!.

    zeta1 : float, dict, :code:`"var"`
        Geometry independent friction coefficient at hot side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    zeta2 : float, dict, :code:`"var"`
        Geometry independent friction coefficient at cold side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    ttd_l : float, dict
        Lower terminal temperature difference :math:`ttd_\mathrm{l}/\text{K}`.

    ttd_u : float, dict
        Upper terminal temperature difference (referring to saturation
        temprature of condensing fluid) :math:`ttd_\mathrm{u}/\text{K}`.

    ttd_min : float, dict
        Minumum terminal temperature difference :math:`ttd_\mathrm{min}/\text{K}`.

    eff_cold : float, dict
        Cold side heat exchanger effectiveness :math:`eff_\text{cold}/\text{1}`.

    eff_hot : float, dict
        Hot side heat exchanger effectiveness :math:`eff_\text{hot}/\text{1}`.

    eff_max : float, dict
        Max value of hot and cold side heat exchanger effectiveness values
        :math:`eff_\text{max}/\text{1}`.

    kA : float, dict
        Area independent heat transfer coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char : tespy.tools.data_containers.SimpleDataContainer
        Area independent heat transfer coefficient characteristic.

    kA_char1 : tespy.tools.characteristics.CharLine, dict
        Characteristic line for hot side heat transfer coefficient.

    kA_char2 : tespy.tools.characteristics.CharLine, dict
        Characteristic line for cold side heat transfer coefficient.

    subcooling : boolean
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
    >>> from tespy.tools.fluid_properties import T_sat_p
    >>> import os
    >>> nw = Network(m_range=[0.01, 1000], iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "temperature": "degC", "enthalpy": "kJ/kg"
    ... })
    >>> amb_in = Source('ambient air inlet')
    >>> amb_out = Sink('air outlet')
    >>> waste_steam = Source('waste steam')
    >>> c = Sink('condensate sink')
    >>> cond = Condenser('condenser')
    >>> amb_he = Connection(amb_in, 'out1', cond, 'in2')
    >>> he_amb = Connection(cond, 'out2', amb_out, 'in1')
    >>> ws_he = Connection(waste_steam, 'out1', cond, 'in1')
    >>> he_c = Connection(cond, 'out1', c, 'in1')
    >>> nw.add_conns(amb_he, he_amb, ws_he, he_c)

    The air flow can not be controlled, thus is constant in offdesign
    operation. If the waste steam mass flow or the ambient air temperature
    change, the outlet temperature of the air will change, too.

    >>> cond.set_attr(pr1=0.98, pr2=0.999, ttd_u=15, design=['pr2', 'ttd_u'],
    ... offdesign=['zeta2', 'kA_char'])
    >>> ws_he.set_attr(fluid={'water': 1}, h=2700, m=1)
    >>> amb_he.set_attr(fluid={'air': 1}, T=20, offdesign=['v'])
    >>> he_amb.set_attr(p=1, T=40, design=['T'])
    >>> nw.solve('design')
    >>> nw.save('tmp.json')
    >>> round(amb_he.v.val, 2)
    103.17
    >>> round(ws_he.T.val - he_amb.T.val, 1)
    66.9
    >>> round(ws_he.calc_T_sat() - 273.15 - he_amb.T.val, 1)
    15.0
    >>> ws_he.set_attr(m=0.7)
    >>> amb_he.set_attr(T=30)
    >>> nw.solve('offdesign', design_path='tmp.json')
    >>> round(ws_he.T.val - he_amb.T.val, 1)
    62.5
    >>> round(ws_he.calc_T_sat() - 273.15 - he_amb.T.val, 1)
    11.3

    It is possible to activate subcooling. The difference to boiling point
    temperature is specified to 5 K.

    >>> cond.set_attr(subcooling=True)
    >>> he_c.set_attr(Td_bp=-5)
    >>> nw.solve('offdesign', design_path='tmp.json')
    >>> round(ws_he.T.val - he_amb.T.val, 1)
    62.5
    >>> round(ws_he.calc_T_sat() - 273.15 - he_amb.T.val, 1)
    13.4
    >>> os.remove('tmp.json')
    """

    def get_parameters(self):
        params = super().get_parameters()
        params.update({
            'subcooling': dc_simple(
                _val=False, num_eq_sets=1,
                func=self.subcooling_func,
                dependents=self.subcooling_dependents,
            )
        })
        return params

    def _preprocess(self, row_idx):

        # if subcooling is True, outlet state method must not be calculated
        self.subcooling.is_set = not self.subcooling.val
        super()._preprocess(row_idx)

    def subcooling_func(self):
        r"""
        Equation for hot side outlet state.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0=h_{out,1} -h\left(p_{out,1}, x=0 \right)

        Note
        ----
        This equation is applied in case subcooling is False!
        """
        o = self.outl[0]
        return o.h.val_SI - h_mix_pQ(o.p.val_SI, 0, o.fluid_data)

    def subcooling_dependents(self):
        return [
            self.outl[0].p,
            self.outl[0].h
        ]

    def calculate_td_log(self):

        i1 = self.inl[0]
        i2 = self.inl[1]
        o1 = self.outl[0]
        o2 = self.outl[1]

        T_i1 = i1.calc_T_sat()
        T_i2 = i2.calc_T()
        T_o1 = o1.calc_T()
        T_o2 = o2.calc_T()

        if T_i1 <= T_o2 and not i1.T.is_set:
            T_i1 = T_o2 + 0.5
        if T_i1 <= T_o2 and not o2.T.is_set:
            T_o2 = T_i1 - 0.5
        if T_o1 <= T_i2 and not o1.T.is_set:
            T_o1 = T_i2 + 1
        if T_o1 <= T_i2 and not i2.T.is_set:
            T_i2 = T_o1 - 1

        ttd_u = T_i1 - T_o2
        ttd_l = T_o1 - T_i2
        if round(ttd_u, 6) == round(ttd_l, 6):
            td_log = ttd_l
        else:
            td_log = (ttd_l - ttd_u) / math.log((ttd_l) / (ttd_u))

        return td_log

    def kA_char_func(self):
        r"""
        Calculate heat transfer from heat transfer coefficient characteristic.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_{in,1} \cdot \left( h_{out,1} - h_{in,1}\right) +
                kA_{design} \cdot f_{kA} \cdot \frac{T_{out,1} -
                T_{in,2} - T_{sat} \left(p_{in,1}\right) + T_{out,2}}
                {\ln{\frac{T_{out,1} - T_{in,2}}
                {T_{sat} \left(p_{in,1}\right) - T_{out,2}}}}

                f_{kA} = \frac{2}{\frac{1}{f_1 \left( expr_1\right)} +
                \frac{1}{f_2 \left( expr_2\right)}}

        Note
        ----
        For standard functions f\ :subscript:`1` \ and f\ :subscript:`2` \ see
        module :ref:`tespy.data <tespy_data_label>`.
        """
        return super().kA_char_func()

    def ttd_u_func(self):
        r"""
        Equation for upper terminal temperature difference.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = ttd_{u} - T_{sat} \left(p_{in,1}\right) + T_{out,2}

        Note
        ----
        The upper terminal temperature difference ttd_u refers to boiling
        temperature at hot side inlet.
        """
        i = self.inl[0]
        o = self.outl[1]
        T_i1 = i.calc_T_sat()
        T_o2 = o.calc_T()
        return self.ttd_u.val - T_i1 + T_o2

    def ttd_u_dependents(self):
        return [
            self.inl[0].p,
            self.outl[1].p,
            self.outl[1].h,
        ]

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        self.Q.val_SI = self.inl[0].m.val_SI * (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI
        )
        self.ttd_u.val_SI = self.inl[0].calc_T_sat() - self.outl[1].T.val_SI
        self.ttd_l.val_SI = self.outl[0].T.val_SI - self.inl[1].T.val_SI
        self.ttd_min.val_SI = min(self.ttd_u.val_SI, self.ttd_l.val_SI)

        # pr and zeta
        for i in range(2):
            self.get_attr(f'pr{i + 1}').val_SI = (
                self.outl[i].p.val_SI / self.inl[i].p.val_SI
            )
            self.get_attr(f'zeta{i + 1}').val_SI = self.calc_zeta(
                self.inl[i], self.outl[i]
            )
            self.get_attr(f'dp{i + 1}').val_SI = (
                self.inl[i].p.val_SI - self.outl[i].p.val_SI
            )

        # kA and logarithmic temperature difference
        if self.ttd_u.val_SI < 0 or self.ttd_l.val_SI < 0:
            self.td_log.val_SI = np.nan
        elif round(self.ttd_l.val_SI, 6) == round(self.ttd_u.val_SI, 6):
            self.td_log.val_SI = self.ttd_l.val_SI
        elif round(self.ttd_l.val_SI, 6) == 0 or round(self.ttd_u.val_SI, 6) == 0:
            self.td_log.val_SI = np.nan
        else:
            self.td_log.val_SI = (
                (self.ttd_l.val_SI - self.ttd_u.val_SI)
                / math.log(self.ttd_l.val_SI / self.ttd_u.val_SI)
            )

        self.kA.val_SI = -self.Q.val_SI / self.td_log.val_SI

        # heat exchanger efficiencies
        try:
            self.eff_hot.val_SI = (
                (self.outl[0].h.val_SI - self.inl[0].h.val_SI)
                / self.calc_dh_max_hot()
            )
        except ValueError:
            self.eff_hot.val_SI = np.nan
            msg = (
                f"Cannot calculate {self.label} hot side effectiveness "
                "because cold side inlet temperature is out of bounds for hot "
                "side fluid."
            )
            logger.warning(msg)
        try:
            self.eff_cold.val_SI = (
                (self.outl[1].h.val_SI - self.inl[1].h.val_SI)
                / self.calc_dh_max_cold()
            )
        except ValueError:
            self.eff_cold.val_SI = np.nan
            msg = (
                f"Cannot calculate {self.label} cold side effectiveness "
                "because hot side inlet temperature is out of bounds for cold "
                "side fluid."
            )
            logger.warning(msg)
        self.eff_max.val_SI = max(self.eff_hot.val_SI, self.eff_cold.val_SI)

    def convergence_check(self):
        o = self.outl[0]
        if o.p.is_var:
            fluid = single_fluid(o.fluid_data)
            p_crit = o.fluid.wrapper[fluid]._p_crit
            if o.p.val_SI > p_crit:
                o.p.set_reference_val_SI(p_crit * 0.9)

    def initialise_source(self, c, key):
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
        """
        if c.source_id == 'out1':
            if key == 'h':
                return h_mix_pQ(c.p.val_SI, 0, c.fluid_data, c.mixing_rule)

        return super().initialise_source(c, key)
