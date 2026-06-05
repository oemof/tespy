# -*- coding: utf-8

"""Module of class HeatExchanger.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/heat_exchangers/base.py

SPDX-License-Identifier: MIT
"""
import math
import warnings

import numpy as np

from tespy.components.component import Component
from tespy.components.component import component_registry
from tespy.tools import logger
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentCharacteristics as dc_gcc
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.fluid_properties import s_mix_ph
from tespy.tools.fluid_properties import single_fluid
from tespy.tools.data_containers import ComponentArrayProperties as dc_cap
from tespy.tools.helpers import _get_dependents
from tespy.tools.helpers import _numeric_deriv


@component_registry
class HeatExchanger(Component):
    r"""
    Class for counter flow heat exchanger.

    The component HeatExchanger is the parent class for the components:

    - :py:class:`tespy.components.heat_exchangers.condenser.Condenser`
    - :py:class:`tespy.components.heat_exchangers.desuperheater.Desuperheater`
    - :py:class:`tespy.components.heat_exchangers.movingboundary.MovingBoundaryHeatExchanger`

    .. image:: /api/_images/components/HeatExchanger.svg
       :alt: flowsheet of the heatexchanger
       :align: center
       :class: only-light

    .. image:: /api/_images/components/HeatExchanger_darkmode.svg
       :alt: flowsheet of the heatexchanger
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
        Effective logarithmic mean temperature difference :code:`Q/UA`. Quantity:
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
        Equation: :py:meth:`ttd_u_func <tespy.components.heat_exchangers.base.HeatExchanger.ttd_u_func>`.

    UA : float, dict
        Heat transfer coefficient considering terminal temperature differences.
        Quantity: :code:`heat_transfer_coefficient`.
        Equation: :py:meth:`UA_func <tespy.components.heat_exchangers.base.HeatExchanger.UA_func>`.

    UA_char : GroupedComponentCharacteristics
        Equation for heat transfer based on UA and modification factor.
        Elements: :code:`UA_char1`, :code:`UA_char2`.
        Equation: :py:meth:`UA_char_func <tespy.components.heat_exchangers.base.HeatExchanger.UA_char_func>`.

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

        The HeatExchanger and subclasses (
        :py:class:`tespy.components.heat_exchangers.condenser.Condenser`,
        :py:class:`tespy.components.heat_exchangers.desuperheater.Desuperheater`)
        are countercurrent heat exchangers. Equations (:code:`kA`, :code:`ttd_u`,
        :code:`ttd_l`) do not work for directcurrent and crosscurrent or
        combinations of different types.

    Example
    -------
    A water cooling is installed to transfer heat from hot exhaust air. The
    heat exchanger is designed for a terminal temperature difference of 5 K.
    From this, it is possible to calculate the heat transfer coefficient and
    predict water and air outlet temperature in offdesign operation.

    >>> from tespy.components import Sink, Source, HeatExchanger
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "pressure_difference": "bar",
    ...     "temperature": "degC", "enthalpy": "kJ/kg",
    ...     "heat_transfer_coefficient": "kW/K"
    ... })
    >>> exhaust_hot = Source('Exhaust air outlet')
    >>> exhaust_cold = Sink('Exhaust air inlet')
    >>> cw_cold = Source('cooling water inlet')
    >>> cw_hot = Sink('cooling water outlet')
    >>> he = HeatExchanger('waste heat exchanger')
    >>> ex_he = Connection(exhaust_hot, 'out1', he, 'in1')
    >>> he_ex = Connection(he, 'out1', exhaust_cold, 'in1')
    >>> cw_he = Connection(cw_cold, 'out1', he, 'in2')
    >>> he_cw = Connection(he, 'out2', cw_hot, 'in1')
    >>> nw.add_conns(ex_he, he_ex, cw_he, he_cw)

    The volumetric flow of the air is at 100 l/s. After designing the component
    it is possible to predict the temperature at different flow rates or
    different inlet temperatures of the exhaust air.

    >>> he.set_attr(
    ...     pr1=0.98, pr2=0.98, ttd_u=5,
    ...     design=['pr1', 'pr2', 'ttd_u'],
    ...     offdesign=['zeta1_d4', 'zeta2_d4', 'UA_char']
    ... )
    >>> cw_he.set_attr(
    ...     fluid={'water': 1}, T=10, p=3, offdesign=['m']
    ... )
    >>> ex_he.set_attr(fluid={'air': 1}, v=0.1, T=35)
    >>> he_ex.set_attr(T=17.5, p=1, design=['T'])
    >>> nw.solve('design')
    >>> design_state = nw.save(as_dict=True)
    >>> round(ex_he.T.val - he_cw.T.val, 0)
    5.0
    >>> ex_he.set_attr(v=0.075)
    >>> nw.solve('offdesign', design_path=design_state)
    >>> round(he_cw.T.val, 1)
    27.5
    >>> round(he_ex.T.val, 1)
    14.4
    >>> ex_he.set_attr(v=0.1, T=40)
    >>> nw.solve('offdesign', design_path=design_state)
    >>> round(he_cw.T.val, 1)
    33.9
    >>> round(he_ex.T.val, 1)
    18.8
    """

    _parameter_aliases = {
        'kA': 'UA',
        'kA_char': 'UA_char',
        'kA_char1': 'UA_char1',
        'kA_char2': 'UA_char2',
        'zeta1': 'zeta1_d4',
        'zeta2': 'zeta2_d4',
        'td_log': 'lmtd',
    }

    def get_parameters(self):
        return {
            'Q': dc_cp(
                max_val=0, num_eq_sets=1,
                func=self.energy_balance_hot_func,
                dependents=self.energy_balance_hot_dependents,
                quantity="heat",
                description="heat transfer from hot side",
                calc=self._calc_Q
            ),
            'UA': dc_cp(
                min_val=0, num_eq_sets=1,
                func=self.UA_func,
                dependents=self.UA_dependents,
                deriv=self.UA_deriv,
                quantity="heat_transfer_coefficient",
                description="heat transfer coefficient considering terminal temperature differences",
                calc=self._calc_UA, calc_deps=['Q', 'ttd_u', 'ttd_l']
            ),
            'kA': dc_cp(
                min_val=0, is_result=True,
                quantity="heat_transfer_coefficient",
                description="deprecated, use :code:`UA` instead",
                calc=self._calc_UA, calc_deps=['Q', 'ttd_u', 'ttd_l']
            ),
            'td_log': dc_cp(
                min_val=0, is_result=True, quantity="temperature_difference",
                description="deprecated, use :code:`lmtd` instead",
                calc=self._calc_lmtd, calc_deps=['lmtd']
            ),
            'lmtd': dc_cp(
                min_val=0, is_result=True, quantity="temperature_difference",
                description="effective logarithmic mean temperature difference :code:`Q/UA`",
                calc=self._calc_lmtd, calc_deps=['Q', 'UA']
            ),
            'ttd_u': dc_cp(
                min_val=0, num_eq_sets=1,
                func=self.ttd_u_func,
                dependents=self.ttd_u_dependents,
                quantity="temperature_difference",
                description="terminal temperature difference at hot side inlet to cold side outlet",
                calc=self._calc_ttd_u
            ),
            'ttd_l': dc_cp(
                min_val=0, num_eq_sets=1,
                func=self.ttd_l_func,
                dependents=self.ttd_l_dependents,
                quantity="temperature_difference",
                description="terminal temperature difference at hot side outlet to cold side inlet",
                calc=self._calc_ttd_l
            ),
            'ttd_min': dc_cp(
                min_val=0, num_eq_sets=1,
                func=self.ttd_min_func,
                dependents=self.ttd_min_dependents,
                quantity="temperature_difference",
                description="minimum terminal temperature difference",
                calc=self._calc_ttd_min, calc_deps=['ttd_u', 'ttd_l']
            ),
            'pr1': dc_cp(
                min_val=1e-4, max_val=1, num_eq_sets=1,
                structure_matrix=self.pr_structure_matrix,
                func_params={'pr': 'pr1'},
                quantity="ratio",
                description="hot side outlet to inlet pressure ratio",
                calc=self._calc_pr
            ),
            'pr2': dc_cp(
                min_val=1e-4, max_val=1, num_eq_sets=1,
                structure_matrix=self.pr_structure_matrix,
                func_params={'pr': 'pr2', 'inconn': 1, 'outconn': 1},
                quantity="ratio",
                description="cold side outlet to inlet pressure ratio",
                calc=self._calc_pr, calc_params={'inconn': 1, 'outconn': 1}
            ),
            'dp1': dc_cp(
                min_val=0, max_val=1e15, num_eq_sets=1,
                structure_matrix=self.dp_structure_matrix,
                func_params={'dp': 'dp1', 'inconn': 0, 'outconn': 0},
                quantity="pressure_difference",
                description="hot side inlet to outlet absolute pressure change",
                calc=self._calc_dp
            ),
            'dp2': dc_cp(
                min_val=0, max_val=1e15, num_eq_sets=1,
                structure_matrix=self.dp_structure_matrix,
                func_params={'dp': 'dp2', 'inconn': 1, 'outconn': 1},
                quantity="pressure_difference",
                description="cold side inlet to outlet absolute pressure change",
                calc=self._calc_dp, calc_params={'inconn': 1, 'outconn': 1}
            ),
            'zeta1_d4': dc_cp(
                min_val=0, max_val=1e15, num_eq_sets=1,
                func=self.zeta_d4_func,
                dependents=self.zeta_d4_dependents,
                func_params={'zeta': 'zeta1_d4'},
                description="hot side geometry-independent friction coefficient zeta/D^4 for pressure loss calculation",
                calc=self._calc_zeta_d4
            ),
            'zeta1': dc_cp(
                min_val=0, max_val=1e15, is_result=True,
                description="deprecated, use :code:`zeta1_d4` instead",
                calc=self._calc_zeta_d4
            ),
            'zeta2_d4': dc_cp(
                min_val=0, max_val=1e15, num_eq_sets=1,
                func=self.zeta_d4_func,
                dependents=self.zeta_d4_dependents,
                func_params={'zeta': 'zeta2_d4', 'inconn': 1, 'outconn': 1},
                description="cold side geometry-independent friction coefficient zeta/D^4 for pressure loss calculation",
                calc=self._calc_zeta_d4, calc_params={'inconn': 1, 'outconn': 1}
            ),
            'zeta2': dc_cp(
                min_val=0, max_val=1e15, is_result=True,
                description="deprecated, use :code:`zeta2_d4` instead",
                calc=self._calc_zeta_d4, calc_params={'inconn': 1, 'outconn': 1}
            ),
            'UA_char': dc_gcc(
                elements=['UA_char1', 'UA_char2'],
                num_eq_sets=1,
                func=self.UA_char_func,
                dependents=self.UA_char_dependents,
                description="equation for heat transfer based on UA and modification factor"
            ),
            'UA_char1': dc_cc(
                param='m',
                description="hot side UA modification lookup table for offdesign"
            ),
            'UA_char2': dc_cc(
                param='m',
                char_params={'type': 'rel', 'inconn': 1, 'outconn': 1},
                description="cold side UA modification lookup table for offdesign"
            ),
            'kA_char': dc_gcc(
                elements=['kA_char1', 'kA_char2'],
                description="deprecated, use :code:`UA_char` instead"
            ),
            'kA_char1': dc_cc(
                param='m',
                description="deprecated, use :code:`UA_char1` instead"
            ),
            'kA_char2': dc_cc(
                param='m',
                char_params={'type': 'rel', 'inconn': 1, 'outconn': 1},
                description="deprecated, use :code:`UA_char2` instead"
            ),
            'eff_hot': dc_cp(
                min_val=0, max_val=1, num_eq_sets=1,
                func=self.eff_hot_func,
                dependents=self.eff_hot_dependents,
                quantity="efficiency",
                description="heat exchanger effectiveness for hot side",
                calc=self._calc_eff_hot
            ),
            'eff_cold': dc_cp(
                min_val=0, max_val=1, num_eq_sets=1,
                func=self.eff_cold_func,
                dependents=self.eff_cold_dependents,
                quantity="efficiency",
                description="heat exchanger effectiveness for cold side",
                calc=self._calc_eff_cold
            ),
            'eff_max': dc_cp(
                min_val=0, max_val=1, num_eq_sets=1,
                func=self.eff_max_func,
                dependents=self.eff_max_dependents,
                quantity="efficiency",
                description="maximum heat exchanger effectiveness",
                calc=self._calc_eff_max, calc_deps=['eff_hot', 'eff_cold']
            ),
            'Q_sections': dc_cap(quantity="heat"),
            'T_hot_sections': dc_cap(quantity="temperature"),
            'T_cold_sections': dc_cap(quantity="temperature"),
            'Q_per_section': dc_cap(quantity="heat"),
            'lmtd_per_section': dc_cap(quantity="temperature_difference"),
        }

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        constraints.update({
            'energy_balance_constraints': dc_cmc(**{
                'func': self.energy_balance_func,
                'dependents': self.energy_balance_dependents,
                'num_eq_sets': 1,
                "description": "hot side to cold side heat transfer equation"
            })
        })
        return constraints

    def get_bypass_constraints(self):
        return {
            'mass_flow_constraints': dc_cmc(**{
                'structure_matrix': self.variable_equality_structure_matrix,
                'num_eq_sets': self.num_i,
                'func_params': {'variable': 'm'}
            }),
            'pressure_constraints': dc_cmc(**{
                'structure_matrix': self.variable_equality_structure_matrix,
                'num_eq_sets': self.num_i,
                'func_params': {'variable': 'p'}
            }),
            'enthalpy_constraints': dc_cmc(**{
                'structure_matrix': self.variable_equality_structure_matrix,
                'num_eq_sets': self.num_i,
                'func_params': {'variable': 'h'}
            }),
            'fluid_constraints': dc_cmc(**{
                'structure_matrix': self.variable_equality_structure_matrix,
                'num_eq_sets': self.num_i,
                'func_params': {'variable': 'fluid'}
            })
        }

    @staticmethod
    def inlets():
        return ['in1', 'in2']

    @staticmethod
    def outlets():
        return ['out1', 'out2']

    def energy_balance_func(self):
        r"""
        Equation for heat exchanger energy balance.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_{in,1} \cdot \left(h_{out,1} - h_{in,1} \right) +
                \dot{m}_{in,2} \cdot \left(h_{out,2} - h_{in,2} \right)
        """
        return (
            self.inl[0].m.val_SI
            * (self.outl[0].h.val_SI - self.inl[0].h.val_SI)
            + self.inl[1].m.val_SI
            * (self.outl[1].h.val_SI - self.inl[1].h.val_SI)
        )

    def energy_balance_dependents(self):
        return [
            self.inl[0].m,
            self.inl[1].m,
            self.inl[0].h,
            self.inl[1].h,
            self.outl[0].h,
            self.outl[1].h
        ]

    def energy_balance_hot_func(self):
        r"""
        Equation for hot side heat exchanger energy balance.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 =\dot{m}_{in,1} \cdot \left(h_{out,1}-h_{in,1}\right)-\dot{Q}
        """
        return self._calc_Q() - self.Q.val_SI

    def energy_balance_hot_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].h,
            self.outl[0].h
        ]

    def _calc_Q(self):
        return self.inl[0].m.val_SI * (self.outl[0].h.val_SI - self.inl[0].h.val_SI)

    def _calc_ttd_u(self):
        return self.inl[0].T.val_SI - self.outl[1].T.val_SI

    def _calc_ttd_l(self):
        return self.outl[0].T.val_SI - self.inl[1].T.val_SI

    def _calc_ttd_min(self):
        return min(self.ttd_u.val_SI, self.ttd_l.val_SI)

    def _calc_lmtd_from_ttd(self):
        ttd_u = self.ttd_u.val_SI
        ttd_l = self.ttd_l.val_SI
        if ttd_u < 0 or ttd_l < 0:
            return np.nan
        elif round(ttd_l, 6) == round(ttd_u, 6):
            return ttd_l
        elif round(ttd_l, 6) == 0 or round(ttd_u, 6) == 0:
            return np.nan
        return (ttd_l - ttd_u) / math.log(ttd_l / ttd_u)

    def _calc_UA(self):
        return -self.Q.val_SI / self._calc_lmtd_from_ttd()

    def _calc_lmtd(self):
        if self.UA.val_SI == 0:
            return np.nan
        return abs(self.Q.val_SI) / self.UA.val_SI

    def _calc_eff_hot(self):
        try:
            return (
                (self.outl[0].h.val_SI - self.inl[0].h.val_SI)
                / self.calc_dh_max_hot()
            )
        except ValueError:
            logger.debug(
                f"Cannot calculate {self.label} hot side effectiveness "
                "because cold side inlet temperature is out of bounds for hot "
                "side fluid."
            )
            return np.nan

    def _calc_eff_cold(self):
        try:
            return (
                (self.outl[1].h.val_SI - self.inl[1].h.val_SI)
                / self.calc_dh_max_cold()
            )
        except ValueError:
            logger.debug(
                f"Cannot calculate {self.label} cold side effectiveness "
                "because hot side inlet temperature is out of bounds for cold "
                "side fluid."
            )
            return np.nan

    def _calc_eff_max(self):
        return max(self.eff_hot.val_SI, self.eff_cold.val_SI)

    def calculate_td_log(self):
        """Method to calculate logarithmic temperature difference during
        iteration. It returns the minimal temperature difference value instead
        of the logarithmic temperature difference if the minimal temperature
        difference is negative during iteration to progress in convergence
        """
        T_i1 = self.inl[0].calc_T()
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

    def UA_func(self):
        r"""
        Calculate heat transfer from heat transfer coefficient.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_{in,1} \cdot \left( h_{out,1} - h_{in,1}\right) +
                UA \cdot \frac{T_{out,1} -
                T_{in,2} - T_{in,1} + T_{out,2}}
                {\ln{\frac{T_{out,1} - T_{in,2}}{T_{in,1} - T_{out,2}}}}
        """
        Q = self.inl[0].m.val_SI * (self.outl[0].h.val_SI - self.inl[0].h.val_SI)
        return Q + self.UA.val_SI * self.calculate_td_log()

    def UA_deriv(self, increment_filter, k, dependents=None):
        r"""
        Partial derivatives of heat transfer coefficient function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        dependents = dependents["scalars"][0]
        f = self.UA_func
        i = self.inl[0]
        o = self.outl[0]
        if i.m.is_var:
            self.jacobian[k, i.m.J_col] = o.h.val_SI - i.h.val_SI

        for var in dependents.difference(_get_dependents([i.m])[0]):
            self._partial_derivative(var, k, f, increment_filter)

    def UA_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].p,
            self.outl[0].h,
            self.inl[1].p,
            self.inl[1].h,
            self.outl[1].p,
            self.outl[1].h,
        ]

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
                T_{in,2} - T_{in,1} + T_{out,2}}
                {\ln{\frac{T_{out,1} - T_{in,2}}{T_{in,1} - T_{out,2}}}}

                f_{UA} = \frac{2}{\frac{1}{f_1\left( expr_1\right)} +
                \frac{1}{f_2\left( expr_2\right)}}

        Note
        ----
        For standard functions f\ :subscript:`1` \ and f\ :subscript:`2` \ see
        module :ref:`tespy.data <data_label>`.
        """
        p1 = self.UA_char1.param
        p2 = self.UA_char2.param
        if self.local_offdesign:
            design_value = self._connection_offdesign[self.inl[0].label][p1]
            actual_value = getattr(self.inl[0], p1).val_SI
            f1 = actual_value / design_value

            design_value = self._connection_offdesign[self.inl[1].label][p2]
            actual_value = getattr(self.inl[1], p2).val_SI
            f2 = actual_value / design_value
        else:
            f1 = self.get_char_expr(p1, **self.UA_char1.char_params)
            f2 = self.get_char_expr(p2, **self.UA_char2.char_params)

        fUA1 = self.UA_char1.char_func.evaluate(f1)
        fUA2 = self.UA_char2.char_func.evaluate(f2)
        fUA = 2 / (1 / fUA1 + 1 / fUA2)

        Q = self.inl[0].m.val_SI * (self.outl[0].h.val_SI - self.inl[0].h.val_SI)
        return Q + self.UA.design * fUA * self.calculate_td_log()

    def UA_char_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].p,
            self.outl[0].h,
            self.inl[1].m,
            self.inl[1].p,
            self.inl[1].h,
            self.outl[1].p,
            self.outl[1].h,
        ]

    def ttd_u_func(self):
        r"""
        Equation for upper terminal temperature difference.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = ttd_{u} - T_{in,1} + T_{out,2}
        """
        i = self.inl[0]
        o = self.outl[1]
        T_i1 = i.calc_T()
        T_o2 = o.calc_T()
        return self.ttd_u.val_SI - T_i1 + T_o2

    def ttd_u_dependents(self):
        return [
            self.inl[0].p,
            self.inl[0].h,
            self.outl[1].p,
            self.outl[1].h,
        ]

    def ttd_l_func(self):
        r"""
        Equation for lower terminal temperature difference.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = ttd_{l} - T_{out,1} + T_{in,2}
        """
        i = self.inl[1]
        o = self.outl[0]
        T_i2 = i.calc_T()
        T_o1 = o.calc_T()
        return self.ttd_l.val_SI - T_o1 + T_i2

    def ttd_l_dependents(self):
        return [
            self.inl[1].p,
            self.inl[1].h,
            self.outl[0].p,
            self.outl[0].h,
        ]

    def ttd_min_func(self):
        r"""
        Equation for minimum terminal temperature difference.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                ttd_{l} = T_{out,1} - T_{in,2}
                ttd_{u} = T_{in,1} - T_{out,2}
                0 = \text{min}\left(ttd_{u}, ttd_{l}\right)

        """
        i = self.inl[1]
        o = self.outl[0]
        T_i2 = i.calc_T()
        T_o1 = o.calc_T()

        i = self.inl[0]
        o = self.outl[1]
        T_i1 = i.calc_T()
        T_o2 = o.calc_T()

        ttd_l = T_o1 - T_i2
        ttd_u = T_i1 - T_o2

        return self.ttd_min.val_SI - min(ttd_l, ttd_u)

    def ttd_min_dependents(self):
        return [
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].p,
            self.outl[0].h,
            self.inl[1].p,
            self.inl[1].h,
            self.outl[1].p,
            self.outl[1].h,
        ]

    def calc_dh_max_cold(self):
        r"""Calculate the theoretical maximum enthalpy increase on the cold side

        Returns
        -------
        float
            Maximum cold side enthalpy increase.

            .. math::

                h\left(p_{out,2}, T_{in,1}\right) - h_{in,2}
        """
        o2 = self.outl[1]
        T_in_hot = self.inl[0].calc_T()
        h_at_T_in_hot = h_mix_pT(
            o2.p.val_SI, T_in_hot, o2.fluid_data, o2.mixing_rule
        )
        return h_at_T_in_hot - self.inl[1].h.val_SI

    def eff_cold_func(self):
        r"""Equation for cold side heat exchanger effectiveness.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \text{eff}_\text{cold} \cdot
                \left(h\left(p_{out,2}, T_{in,1} \right) - h_{in,2}\right)
                - \left( h_{out,2} - h_{in,2} \right)
        """
        return (
            self.eff_cold.val_SI * self.calc_dh_max_cold()
            - (self.outl[1].h.val_SI - self.inl[1].h.val_SI)
        )

    def eff_cold_dependents(self):
        return [
            self.inl[0].p,
            self.inl[0].h,
            self.inl[1].h,
            self.outl[1].p,
            self.outl[1].h,
        ]

    def calc_dh_max_hot(self):
        r"""Calculate the theoretical maximum enthalpy decrease on the hot side

        Returns
        -------
        float
            Maximum hot side enthalpy decrease.

            .. math::

                h\left(p_{out,1}, T_{in,2}\right) - h_{in,1}
        """
        o1 = self.outl[0]
        T_in_cold = self.inl[1].calc_T()
        h_at_T_in_cold = h_mix_pT(
            o1.p.val_SI, T_in_cold, o1.fluid_data, o1.mixing_rule
        )
        return h_at_T_in_cold - self.inl[0].h.val_SI

    def eff_hot_func(self):
        r"""Equation for hot side heat exchanger effectiveness.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \text{eff}_\text{hot} \cdot
                \left(h\left(p_{out,1}, T_{in,2}\right) - h_{in,1}\right)
                - \left( h_{out,1} - h_{in,1}\right)
        """
        return (
            self.eff_hot.val_SI * self.calc_dh_max_hot()
            - (self.outl[0].h.val_SI - self.inl[0].h.val_SI)
        )

    def eff_hot_dependents(self):
        return [
            self.inl[0].h,
            self.inl[1].p,
            self.inl[1].h,
            self.outl[0].p,
            self.outl[0].h,
        ]

    def eff_max_func(self):
        r"""Equation for maximum heat exchanger effectiveness.

        .. note::

            This functions works on what is larger: hot side or cold side
            effectiveness. It may cause numerical issues, if applied, when one
            of both sides' effectiveness is already predetermined, e.g. by
            temperature specifications.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \text{eff}_\text{max} - \text{max}
                \left(\text{eff}_\text{hot},\text{eff}_\text{cold}\right)
        """
        eff_hot = (
            (self.outl[0].h.val_SI - self.inl[0].h.val_SI)
            / self.calc_dh_max_hot()
        )
        eff_cold = (
            (self.outl[1].h.val_SI - self.inl[1].h.val_SI)
            / self.calc_dh_max_cold()
        )
        return self.eff_max.val_SI - max(eff_hot, eff_cold)

    def eff_max_dependents(self):
        return [
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].p,
            self.outl[0].h,
            self.inl[1].p,
            self.inl[1].h,
            self.outl[1].p,
            self.outl[1].h,
        ]

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
        if key == 'p':
            fluid = single_fluid(c.fluid_data)
            if fluid is not None:
                return c.fluid.wrapper[fluid]._p_crit / 2
            else:
                return 1e5

        elif key == 'h':
            fluid = single_fluid(c.fluid_data)
            if fluid is not None:
                temp = c.fluid.wrapper[fluid]._T_crit
                if temp is None:
                    temp = c.fluid.wrapper[fluid]._T_max

                dT = temp - c.fluid.wrapper[fluid]._T_min

                if c.source_id == 'out1':
                    temp = temp - dT * 2 / 3
                else:
                    temp = temp - dT / 3
            else:
                if c.source_id == 'out1':
                    temp = 600
                else:
                    temp = 500
            return h_mix_pT(c.p.val_SI, temp, c.fluid_data, c.mixing_rule)

    def initialise_target(self, c, key):
        r"""
        Return a starting value for pressure and enthalpy at inlet.

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
        if key == 'p':
            fluid = single_fluid(c.fluid_data)
            if fluid is not None:
                return c.fluid.wrapper[fluid]._p_crit / 2
            else:
                return 1e5
        elif key == 'h':
            fluid = single_fluid(c.fluid_data)
            if fluid is not None:
                temp = c.fluid.wrapper[fluid]._T_crit
                if temp is None:
                    temp = c.fluid.wrapper[fluid]._T_max

                dT = temp - c.fluid.wrapper[fluid]._T_min

                if c.target_id == 'in1':
                    temp = temp - dT / 4
                else:
                    temp = temp - dT / 3

            else:
                if c.target_id == 'in1':
                    temp = 650
                else:
                    temp = 550

            return h_mix_pT(c.p.val_SI, temp, c.fluid_data, c.mixing_rule)

    def entropy_balance(self):
        r"""
        Calculate entropy balance of a heat exchanger.

        The allocation of the entropy streams due to heat exchanged and due to
        irreversibility is performed by solving for T on both sides of the heat
        exchanger:

        .. math::

            h_\text{out} - h_\text{in} = \int_\text{in}^\text{out} v
            \cdot dp - \int_\text{in}^\text{out} T \cdot ds

        As solving :math:`\int_\text{in}^\text{out} v \cdot dp` for non
        isobaric processes would require perfect process knowledge (the path)
        on how specific volume and pressure change throughout the component, the
        heat transfer is split into three separate virtual processes for
        both sides:

        - in->in*: decrease pressure to
          :math:`p_\text{in*}=p_\text{in}\cdot\sqrt{\frac{p_\text{out}}{p_\text{in}}}`
          without changing enthalpy.
        - in*->out* transfer heat without changing pressure.
          :math:`h_\text{out*}-h_\text{in*}=h_\text{out}-h_\text{in}`
        - out*->out decrease pressure to outlet pressure :math:`p_\text{out}`
          without changing enthalpy.

        Note
        ----
        The entropy balance makes the following parameter available:

        .. math::

            \text{S\_Q1}=\dot{m} \cdot \left(s_\text{out*,1}-s_\text{in*,1}
            \right)\\
            \text{S\_Q2}=\dot{m} \cdot \left(s_\text{out*,2}-s_\text{in*,2}
            \right)\\
            \text{S\_Qirr}=\text{S\_Q2} - \text{S\_Q1}\\
            \text{S\_irr1}=\dot{m} \cdot \left(s_\text{out,1}-s_\text{in,1}
            \right) - \text{S\_Q1}\\
            \text{S\_irr2}=\dot{m} \cdot \left(s_\text{out,2}-s_\text{in,2}
            \right) - \text{S\_Q2}\\
            \text{S\_irr}=\sum \dot{S}_\text{irr}\\
            \text{T\_mQ1}=\frac{\dot{Q}}{\text{S\_Q1}}\\
            \text{T\_mQ2}=\frac{\dot{Q}}{\text{S\_Q2}}
        """
        self.S_irr = 0
        for i in range(2):
            inl = self.inl[i]
            out = self.outl[i]
            p_star = inl.p.val_SI * (
                self.get_attr('pr' + str(i + 1)).val) ** 0.5
            s_i_star = s_mix_ph(
                p_star, inl.h.val_SI, inl.fluid_data, inl.mixing_rule,
                T0=inl.T.val_SI
            )
            s_o_star = s_mix_ph(
                p_star, out.h.val_SI, out.fluid_data, out.mixing_rule,
                T0=out.T.val_SI
            )

            setattr(
                self, 'S_Q' + str(i + 1),
                inl.m.val_SI * (s_o_star - s_i_star)
            )
            S_Q = self.get_attr('S_Q' + str(i + 1))
            setattr(
                self, 'S_irr' + str(i + 1),
                inl.m.val_SI * (out.s.val_SI - inl.s.val_SI) - S_Q
            )
            setattr(
                self, 'T_mQ' + str(i + 1),
                inl.m.val_SI * (out.h.val_SI - inl.h.val_SI) / S_Q
            )

            self.S_irr += self.get_attr('S_irr' + str(i + 1))

        self.S_irr += self.S_Q1 + self.S_Q2

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
            i + 1: {
                'isoline_property': 'p',
                'isoline_value': self.inl[i].p.val,
                'isoline_value_end': self.outl[i].p.val,
                'starting_point_property': 'vol',
                'starting_point_value': self.inl[i].vol.val,
                'ending_point_property': 'vol',
                'ending_point_value': self.outl[i].vol.val
            } for i in range(2)}

    @staticmethod
    def _get_Q_sections(h_at_steps, mass_flow):
        """Calculate the heat exchange of every section given steps of
        enthalpy and mass flow.

        Parameters
        ----------
        h_at_steps : list
            Enthalpy values at sections (inlet, phase change points, outlet)
        mass_flow : float
            Mass flow value

        Returns
        -------
        float
            Heat exchanged between defined steps of enthalpy.
        """
        return np.diff(h_at_steps) * mass_flow

    @staticmethod
    def _assign_to_steps(start, end, steps):
        return start + steps * (end - start)

    def _get_Q_cumsum_steps(self, steps):
        """Assign the sections of the heat exchanger

        Returns
        -------
        list
            List of cumulative sum of heat exchanged defining the heat exchanger
            sections.
        """
        start = self.outl[0].h.val_SI
        end = self.inl[0].h.val_SI

        h_steps_hot = self._assign_to_steps(start, end, steps)
        Q_sections_hot = self._get_Q_sections(h_steps_hot, self.inl[0].m.val_SI)
        return np.insert(np.cumsum(Q_sections_hot), 0, 0.0)

    def _assign_steps(self):
        """Assign the sections of the heat exchanger

        Returns
        -------
        list
            List of cumulative sum of heat exchanged defining the heat exchanger
            sections.
        """
        return np.array([0, 1])

    def _preprocess(self, row_idx):
        self._T_cache_hot = {}
        self._T_cache_cold = {}
        super()._preprocess(row_idx)

    def _get_T_at_steps(self, steps):
        """Calculate the temperature values for the provided sections.

        Parameters
        ----------
        Q_sections : list
            Cumulative heat exchanged from the hot side to the cold side
            defining the sections of the heat exchanger.

        Returns
        -------
        tuple
            Lists of cold side and hot side temperature
        """
        h_steps_hot = self._assign_to_steps(
            self.outl[0].h.val_SI, self.inl[0].h.val_SI, steps
        )
        p_steps_hot = self._assign_to_steps(
            self.outl[0].p.val_SI, self.inl[0].p.val_SI, steps
        )
        h_steps_cold = self._assign_to_steps(
            self.inl[1].h.val_SI, self.outl[1].h.val_SI, steps
        )
        p_steps_cold = self._assign_to_steps(
            self.inl[1].p.val_SI, self.outl[1].p.val_SI, steps
        )

        T_steps_hot = np.empty(len(steps))
        for i, (p, h) in enumerate(zip(p_steps_hot, h_steps_hot)):
            key = (p, h)
            if key not in self._T_cache_hot:
                self._T_cache_hot[key] = T_mix_ph(
                    p, h, self.inl[0].fluid_data, self.inl[0].mixing_rule,
                )
            T_steps_hot[i] = self._T_cache_hot[key]

        T_steps_cold = np.empty(len(steps))
        for i, (p, h) in enumerate(zip(p_steps_cold, h_steps_cold)):
            key = (p, h)
            if key not in self._T_cache_cold:
                self._T_cache_cold[key] = T_mix_ph(
                    p, h, self.inl[1].fluid_data, self.inl[1].mixing_rule,
                )
            T_steps_cold[i] = self._T_cache_cold[key]

        return T_steps_hot, T_steps_cold

    @staticmethod
    def _calc_lmtd_per_section(T_steps_hot, T_steps_cold, postprocess=False):
        """Calculate the logarithmic temperature difference values per section
        of heat exchanged.

        Parameters
        ----------
        T_steps_hot : list
            Temperature hot side at beginning and end of sections.

        T_steps_cold : list
            Temperature cold side at beginning and end of sections.

        Returns
        -------
        list
            Lists of temperature differences per section of heat exchanged.
        """
        td_at_steps = T_steps_hot - T_steps_cold
        if postprocess:
            if (td_at_steps <= 0).any():
                return np.ones(len(td_at_steps) - 1) * np.nan
        td_at_steps[td_at_steps <= 0] = 1e-3
        return np.array([
            (td_at_steps[i + 1] - td_at_steps[i])
            / math.log(td_at_steps[i + 1] / td_at_steps[i])
            # round is required because tiny differences may cause
            # inconsistencies due to rounding errors
            if round(td_at_steps[i + 1], 6) != round(td_at_steps[i], 6)
            else td_at_steps[i + 1]
            for i in range(len(td_at_steps) - 1)
        ])

    def calc_parameters(self):
        self._store_sections()
        super().calc_parameters()

    def _calc_sections_SI(self, postprocess=True):
        """Compute section data in SI units. Used internally during solving
        and as the basis for :py:meth:`calc_sections`.
        """
        steps = self._assign_steps()
        Q_sections = self._get_Q_cumsum_steps(steps)
        T_steps_hot, T_steps_cold = self._get_T_at_steps(steps)
        Q_per_section = np.diff(Q_sections)
        lmtd_per_section = self._calc_lmtd_per_section(
            T_steps_hot, T_steps_cold, postprocess
        )
        return Q_sections, T_steps_hot, T_steps_cold, Q_per_section, lmtd_per_section

    def _store_sections(self):
        """Compute section data and store results as
        :py:class:`ComponentArrayProperties <tespy.tools.data_containers.ComponentArrayProperties>`
        attributes. Each attribute exposes :code:`.val` in network units and
        :code:`.val_SI` in SI units.

        Attributes set
        --------------
        Q_sections, T_hot_sections, T_cold_sections, Q_per_section,
        lmtd_per_section
        """
        Q_si, T_hot_si, T_cold_si, Q_per_si, lmtd_si = self._calc_sections_SI(postprocess=True)
        self.Q_sections.val_SI = Q_si
        self.T_hot_sections.val_SI = T_hot_si
        self.T_cold_sections.val_SI = T_cold_si
        self.Q_per_section.val_SI = Q_per_si
        self.lmtd_per_section.val_SI = lmtd_si

    def calc_sections(self, postprocess=True):
        r"""
        .. deprecated::
            Use the component attributes :code:`Q_sections`, :code:`T_hot_sections`,
            :code:`T_cold_sections`, :code:`Q_per_section`, :code:`lmtd_per_section`
            instead. These are populated automatically after each solve. The return
            value of this method will be removed in a future version.
        """
        warnings.warn(
            f"The return value of {self.__class__.__name__}.calc_sections() is deprecated. "
            "Access section data via the component attributes Q_sections, T_hot_sections, "
            "T_cold_sections, Q_per_section, lmtd_per_section instead. Each attribute "
            "exposes .val (network units) and .val_SI (SI units). The return value will "
            "be removed in a future version.",
            FutureWarning,
            stacklevel=2,
        )
        return self._calc_sections_SI(postprocess=postprocess)
