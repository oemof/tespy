# -*- coding: utf-8

"""Module of class Condenser.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/heat_exchangers/condenser.py

SPDX-License-Identifier: MIT
"""

import math

from tespy.components.component import component_registry
from tespy.components.heat_exchangers.base import HeatExchanger
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
from tespy.tools.fluid_properties import h_mix_pQ
from tespy.tools.fluid_properties import single_fluid


@component_registry
class Condenser(HeatExchanger):
    r"""
    A Condenser cools a fluid until it is in liquid state.

    The condensing fluid is cooled by the cold side fluid. The fluid on the hot
    side of the condenser must be pure. Subcooling is available.

    .. image:: /api/_images/components/Condenser.svg
       :alt: flowsheet of the condenser
       :align: center
       :class: only-light

    .. image:: /api/_images/components/Condenser_darkmode.svg
       :alt: flowsheet of the condenser
       :align: center
       :class: only-dark

    Ports
    -----

    - Fluid inlets: in1, in2
    - Fluid outlets: out1, out2

    Mandatory Equations
    -------------------

    - mass flow equality constraint(s): :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
    - fluid composition equality constraint(s): :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
    - hot side to cold side heat transfer equation: :py:meth:`energy_balance_func <tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_func>`

    Parameters
    ----------

    char_warnings : bool
        Ignore warnings on default characteristics usage for this component.

    design : list
        List containing design parameters (stated as String).

    design_path : str
        Path to the components design case.

    dp1 : float, dict
        Hot side inlet to outlet absolute pressure change. Quantity:
        :code:`pressure_difference`.
        Equation: :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`.

    dp2 : float, dict
        Cold side inlet to outlet absolute pressure change. Quantity:
        :code:`pressure_difference`.
        Equation: :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`.

    eff_cold : float, dict
        Heat exchanger effectiveness for cold side. Quantity:
        :code:`efficiency`.
        Equation: :py:meth:`eff_cold_func <tespy.components.heat_exchangers.base.HeatExchanger.eff_cold_func>`.

    eff_hot : float, dict
        Heat exchanger effectiveness for hot side. Quantity: :code:`efficiency`.
        Equation: :py:meth:`eff_hot_func <tespy.components.heat_exchangers.base.HeatExchanger.eff_hot_func>`.

    eff_max : float, dict
        Maximum heat exchanger effectiveness. Quantity: :code:`efficiency`.
        Equation: :py:meth:`eff_max_func <tespy.components.heat_exchangers.base.HeatExchanger.eff_max_func>`.

    kA : float, dict
        Deprecated, use :code:`UA` instead. Quantity:
        :code:`heat_transfer_coefficient`.

    kA_char : GroupedComponentCharacteristics
        Deprecated, use :code:`UA_char` instead. Elements: :code:`kA_char1`,
        :code:`kA_char2`.

    kA_char1 : tespy.tools.characteristics.CharLine, dict
        Deprecated, use :code:`UA_char1` instead.

    kA_char2 : tespy.tools.characteristics.CharLine, dict
        Deprecated, use :code:`UA_char2` instead.

    label : str
        The label of the component.

    lmtd : float, dict
        Effective logarithmic mean temperature difference |Q|/UA. Quantity:
        :code:`temperature_difference`.

    local_design : bool
        Treat this component in design mode in an offdesign calculation.

    local_offdesign : bool
        Treat this component in offdesign mode in a design calculation.

    offdesign : list
        List containing offdesign parameters (stated as String).

    pr1 : float, dict
        Hot side outlet to inlet pressure ratio. Quantity: :code:`ratio`.
        Equation: :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`.

    pr2 : float, dict
        Cold side outlet to inlet pressure ratio. Quantity: :code:`ratio`.
        Equation: :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`.

    printout : bool
        Include this component in the network's results printout.

    Q : float, dict
        Heat transfer from hot side. Quantity: :code:`heat`.
        Equation: :py:meth:`energy_balance_hot_func <tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_hot_func>`.

    subcooling : bool
        Allow subcooling in the condenser.
        Equation: :py:meth:`subcooling_func <tespy.components.heat_exchangers.condenser.Condenser.subcooling_func>`.

    td_log : float, dict
        Deprecated, use :code:`lmtd` instead. Quantity:
        :code:`temperature_difference`.

    ttd_l : float, dict
        Terminal temperature difference at hot side outlet to cold side inlet.
        Quantity: :code:`temperature_difference`.
        Equation: :py:meth:`ttd_l_func <tespy.components.heat_exchangers.base.HeatExchanger.ttd_l_func>`.

    ttd_min : float, dict
        Minimum terminal temperature difference. Quantity:
        :code:`temperature_difference`.
        Equation: :py:meth:`ttd_min_func <tespy.components.heat_exchangers.base.HeatExchanger.ttd_min_func>`.

    ttd_u : float, dict
        Terminal temperature difference at hot side inlet to cold side outlet.
        Quantity: :code:`temperature_difference`.
        Equation: :py:meth:`ttd_u_func <tespy.components.heat_exchangers.condenser.Condenser.ttd_u_func>`.

    UA : float, dict
        Heat transfer coefficient considering terminal temperature differences.
        Quantity: :code:`heat_transfer_coefficient`.
        Equation: :py:meth:`UA_func <tespy.components.heat_exchangers.base.HeatExchanger.UA_func>`.

    UA_char : GroupedComponentCharacteristics
        Equation for heat transfer based on UA and modification factor.
        Elements: :code:`UA_char1`, :code:`UA_char2`.
        Equation: :py:meth:`UA_char_func <tespy.components.heat_exchangers.condenser.Condenser.UA_char_func>`.

    UA_char1 : tespy.tools.characteristics.CharLine, dict
        Hot side UA modification lookup table for offdesign.

    UA_char2 : tespy.tools.characteristics.CharLine, dict
        Cold side UA modification lookup table for offdesign.

    zeta1 : float, dict
        Deprecated, use :code:`zeta1_d4` instead.

    zeta1_d4 : float, dict
        Hot side geometry-independent friction coefficient zeta/D^4 for pressure
        loss calculation.
        Equation: :py:meth:`zeta_d4_func <tespy.components.component.Component.zeta_d4_func>`.

    zeta2 : float, dict
        Deprecated, use :code:`zeta2_d4` instead.

    zeta2_d4 : float, dict
        Cold side geometry-independent friction coefficient zeta/D^4 for
        pressure loss calculation.
        Equation: :py:meth:`zeta_d4_func <tespy.components.component.Component.zeta_d4_func>`.

    Notes
    -----

    .. note::

        The condenser has an additional equation for enthalpy at hot side outlet:
        The fluid leaves the component in saturated liquid state. If subcooling
        is activated, it is possible to specify the enthalpy at the outgoing
        connection manually.

        It has a different calculation method for given heat transfer coefficient and
        upper terminal temperature difference: These parameters refer to the
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
    >>> nw = Network(m_range=[0.01, 1000], iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "pressure_difference": "bar",
    ...     "temperature": "degC", "enthalpy": "kJ/kg"
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
    ... offdesign=['zeta2_d4', 'UA_char'])
    >>> ws_he.set_attr(fluid={'water': 1}, h=2700, m=1)
    >>> amb_he.set_attr(fluid={'air': 1}, T=20, offdesign=['v'])
    >>> he_amb.set_attr(p=1, T=40, design=['T'])
    >>> nw.solve('design')
    >>> design_state = nw.save(as_dict=True)
    >>> round(amb_he.v.val, 2)
    103.17
    >>> round(ws_he.T.val - he_amb.T.val, 1)
    66.9
    >>> round(ws_he.calc_T_sat() - 273.15 - he_amb.T.val, 1)
    15.0
    >>> ws_he.set_attr(m=0.7)
    >>> amb_he.set_attr(T=30)
    >>> nw.solve('offdesign', design_path=design_state)
    >>> round(ws_he.T.val - he_amb.T.val, 1)
    62.5
    >>> round(ws_he.calc_T_sat() - 273.15 - he_amb.T.val, 1)
    11.3

    It is possible to activate subcooling. The difference to boiling point
    temperature is specified to 5 K.

    >>> cond.set_attr(subcooling=True)
    >>> he_c.set_attr(td_bubble=5)
    >>> nw.solve('offdesign', design_path=design_state)
    >>> round(ws_he.T.val - he_amb.T.val, 1)
    62.5
    >>> round(ws_he.calc_T_sat() - 273.15 - he_amb.T.val, 1)
    13.4
    """

    def get_parameters(self):
        params = super().get_parameters()
        params.update({
            'subcooling': dc_simple(
                _val=False, dtype="bool", num_eq_sets=1,
                func=self.subcooling_func,
                dependents=self.subcooling_dependents,
                description="allow subcooling in the condenser"
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
        T_i1 = self.inl[0].calc_T_sat()
        T_i2 = self.inl[1].calc_T()
        T_o1 = self.outl[0].calc_T()
        T_o2 = self.outl[1].calc_T()

        ttd_u = T_i1 - T_o2
        ttd_l = T_o1 - T_i2
        min_ttd = min(ttd_u, ttd_l)
        if min_ttd <= 0:
            return min_ttd
        if round(ttd_u, 6) == round(ttd_l, 6):
            return ttd_l
        return (ttd_l - ttd_u) / math.log(ttd_l / ttd_u)

    def UA_char_func(self):
        r"""
        Calculate heat transfer from heat transfer coefficient characteristic.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_{in,1} \cdot \left( h_{out,1} - h_{in,1}\right) +
                UA_{design} \cdot f_{UA} \cdot \frac{T_{out,1} -
                T_{in,2} - T_{sat} \left(p_{in,1}\right) + T_{out,2}}
                {\ln{\frac{T_{out,1} - T_{in,2}}
                {T_{sat} \left(p_{in,1}\right) - T_{out,2}}}}

                f_{UA} = \frac{2}{\frac{1}{f_1 \left( expr_1\right)} +
                \frac{1}{f_2 \left( expr_2\right)}}

        Note
        ----
        For standard functions f\ :subscript:`1` \ and f\ :subscript:`2` \ see
        module :ref:`tespy.data <data_label>`.
        """
        return super().UA_char_func()

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
        return self.ttd_u.val_SI - T_i1 + T_o2

    def ttd_u_dependents(self):
        return [
            self.inl[0].p,
            self.outl[1].p,
            self.outl[1].h,
        ]

    def _calc_ttd_u(self):
        return self.inl[0].T_dew.val_SI - self.outl[1].T.val_SI

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
