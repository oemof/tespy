# -*- coding: utf-8

"""Module of class SimpleHeatExchanger.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/heat_exchangers/simple.py

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
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.fluid_properties import s_mix_ph
from tespy.tools.fluid_properties.helpers import darcy_friction_factor as dff
from tespy.tools.helpers import TESPyNetworkError
from tespy.tools.helpers import _numeric_deriv


@component_registry
class SimpleHeatExchanger(Component):
    r"""
    A basic heat exchanger representing a heat source or heat sink.

    The component SimpleHeatExchanger is the parent class for the components:

    - :py:class:`tespy.components.heat_exchangers.solar_collector.SolarCollector`
    - :py:class:`tespy.components.heat_exchangers.parabolic_trough.ParabolicTrough`
    - :py:class:`tespy.components.piping.pipe.Pipe`

    **Mandatory Equations**

    - fluid: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`
    - mass flow: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`

    **Optional Equations**

    - :py:meth:`tespy.components.component.Component.pr_structure_matrix`
    - :py:meth:`tespy.components.component.Component.dp_structure_matrix`
    - :py:meth:`tespy.components.component.Component.zeta_d4_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.energy_balance_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.darcy_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.hazen_williams_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.UA_group_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.UA_char_group_func`

    Inlets/Outlets

    - in1
    - out1

    Image

    .. image:: /api/_images/Pipe.svg
       :alt: flowsheet of the simple heat exchanger
       :align: center
       :class: only-light

    .. image:: /api/_images/Pipe_darkmode.svg
       :alt: flowsheet of the simple heat exchanger
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

    Q : float, dict, :code:`"var"`
        Heat transfer, :math:`Q/\text{W}`.

    pr : float, dict, :code:`"var"`
        Outlet to inlet pressure ratio, :math:`pr/1`.

    zeta_d4 : float, dict, :code:`"var"`
        Geometry independent friction coefficient,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    D : float, dict, :code:`"var"`
        Diameter of the pipes, :math:`D/\text{m}`.

    L : float, dict, :code:`"var"`
        Length of the pipes, :math:`L/\text{m}`.

    ks : float, dict, :code:`"var"`
        Pipe's roughness, :math:`ks/\text{m}`.

    darcy_group : str, dict
        Parametergroup for pressure drop calculation based on pipes dimensions
        using darcy weissbach equation.

    ks_HW : float, dict, :code:`"var"`
        Pipe's roughness, :math:`ks/\text{1}`.

    hw_group : str, dict
        Parametergroup for pressure drop calculation based on pipes dimensions
        using hazen williams equation.

    UA : float, dict, :code:`"var"`
        Area independent heat transfer coefficient,
        :math:`UA/\frac{\text{W}}{\text{K}}`.

    UA_char : tespy.tools.characteristics.CharLine, dict
        Characteristic line for heat transfer coefficient.

    Tamb : float, dict
        Ambient temperature, provide parameter in network's temperature unit.

    lmtd : float
        Effective logarithmic mean temperature difference :math:`lmtd/\text{K}`,
        defined as :math:`|\dot{Q}|/UA`.

    UA_group : str, dict
        Parametergroup for heat transfer calculation from ambient temperature
        and area independent heat transfer coefficient UA.

    Example
    -------
    The SimpleHeatExchanger can be used as a sink or source of heat. This
    component does not simulate the secondary side of the heat exchanger. It
    is possible to calculate the pressure ratio with the Darcy-Weisbach
    equation or in case of liquid water use the Hazen-Williams equation.
    Also, given ambient temperature and the heat transfer coefficient, it is
    possible to predict heat transfer.

    >>> from tespy.components import Sink, Source, SimpleHeatExchanger
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "pressure_difference": "bar",
    ...     "temperature": "degC", "enthalpy": "kJ/kg"
    ... })
    >>> so1 = Source('source 1')
    >>> si1 = Sink('sink 1')
    >>> heat_sink = SimpleHeatExchanger('heat sink')
    >>> heat_sink.set_attr(Tamb=10, pr=0.95, design=['pr'],
    ... offdesign=['zeta_d4', 'UA_char'])
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
    >>> design_state = nw.save(as_dict=True)
    >>> round(heat_sink.Q.val, 0)
    -52581.0
    >>> round(heat_sink.kA.val, 0)
    321.0
    >>> inc.set_attr(m=1.25)
    >>> nw.solve('offdesign', design_path=design_state)
    >>> round(heat_sink.Q.val, 0)
    -56599.0
    >>> round(outg.T.val, 1)
    156.9
    >>> inc.set_attr(m=0.75)
    >>> nw.solve('offdesign', design_path=design_state)
    >>> round(heat_sink.Q.val, 1)
    -47275.8
    >>> round(outg.T.val, 1)
    140.0

    Use of the HeatConnection
    -------------------------

    Single sided heat exchangers can also connect to a
    :code:`HeatConnection`. The component auto-detects whether the connection
    is on the inlet or outlet side - no prior declaration is needed.

    >>> from tespy.connections import HeatConnection
    >>> from tespy.components import HeatSink
    >>> ambient = HeatSink('ambient heat dissipation')

    Create and add the :code:`HeatConnection` as an outlet. We run a new
    design calculation, because the old design case did not include the
    :code:`HeatConnection`. The energy value will be identical to the heat
    transfer of the component.

    >>> h1 = HeatConnection(heat_sink, 'heat', ambient, 'heat', label='h1')
    >>> nw.add_conns(h1)
    >>> nw.solve('design')
    >>> round(h1.E.val) == round(-heat_sink.Q.val)
    True
    """

    _parameter_aliases = {
        'kA': 'UA',
        'kA_char': 'UA_char',
        'kA_group': 'UA_group',
        'kA_char_group': 'UA_char_group',
        'zeta': 'zeta_d4',
    }

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        if self.power_inl + self.power_outl + self.heat_inl + self.heat_outl:
            constraints["energy_connector_balance"] = dc_cmc(**{
                "func": self.energy_connector_balance_func,
                "dependents": self.energy_connector_dependents,
                "num_eq_sets": 1
            })

        return constraints

    def set_attr(self, **kwargs):
        if 'power_connector_location' in kwargs:
            warnings.warn(
                "The parameter 'power_connector_location' is deprecated and has no "
                "effect. Connect the component directly on either the inlet or outlet "
                "side without prior declaration.",
                FutureWarning,
                stacklevel=2,
            )
        super().set_attr(**kwargs)

    def _calc_Q(self):
        return self.inl[0].m.val_SI * (self.outl[0].h.val_SI - self.inl[0].h.val_SI)

    def _calc_UA(self):
        if not self.Tamb.is_set:
            return np.nan
        ttd_1 = self.inl[0].T.val_SI - self.Tamb.val_SI
        ttd_2 = self.outl[0].T.val_SI - self.Tamb.val_SI
        if ttd_1 / ttd_2 < 0:
            return np.nan
        return abs(self.Q.val_SI / self._calculate_td_log())

    def _calc_lmtd(self):
        if self.UA.val_SI == 0:
            return np.nan
        return abs(self.Q.val_SI) / self.UA.val_SI

    def get_parameters(self):
        return {
            'power_connector_location': dc_simple(),
            'Q': dc_cp(
                num_eq_sets=1,
                func=self.energy_balance_func,
                dependents=self.energy_balance_dependents,
                quantity="heat",
                description="heat transfer",
                calc=self._calc_Q
            ),
            'pr': dc_cp(
                min_val=1e-4, max_val=1, num_eq_sets=1,
                structure_matrix=self.pr_structure_matrix,
                func_params={'pr': 'pr'},
                quantity="ratio",
                description="outlet to inlet pressure ratio",
                calc=self._calc_pr
            ),
            'dp': dc_cp(
                min_val=0, max_val=1e15, num_eq_sets=1,
                structure_matrix=self.dp_structure_matrix,
                func_params={'dp': 'dp'},
                quantity="pressure_difference",
                description="inlet to outlet absolute pressure change",
                calc=self._calc_dp
            ),
            'zeta_d4': dc_cp(
                min_val=0, max_val=1e15, num_eq_sets=1,
                func=self.zeta_d4_func,
                dependents=self.zeta_d4_dependents,
                func_params={'zeta': 'zeta_d4'},
                description="geometry-independent friction coefficient zeta/D^4 for pressure loss calculation",
                calc=self._calc_zeta_d4
            ),
            'zeta': dc_cp(
                min_val=0, is_result=True,
                description="deprecated, use :code:`zeta_d4` instead",
                calc=self._calc_zeta_d4
            ),
            'D': dc_cp(
                min_val=1e-2, max_val=2, d=1e-5, quantity="length",
                description="diameter of channel",
                _potential_var=True
            ),
            'L': dc_cp(
                min_val=1e-1, quantity="length",
                description="length of channel",
                _potential_var=True
            ),
            'ks': dc_cp(
                _val=1e-4, min_val=1e-7, max_val=1e-3,
                quantity="length", description="roughness of wall material",
                _potential_var=True
            ),
            'ks_HW': dc_cp(
                _val=10, min_val=1e-1, max_val=1e3,
                description="Hazen-Williams roughness",
                _potential_var=True
            ),
            'UA': dc_cp(
                min_val=0, quantity="heat_transfer_coefficient",
                description="heat transfer coefficient considering ambient temperature",
                _potential_var=True,
                calc=self._calc_UA, calc_deps=['Q']
            ),
            'kA': dc_cp(
                min_val=0, quantity="heat_transfer_coefficient",
                description="deprecated, use :code:`UA` instead",
                _potential_var=True,
                calc=self._calc_UA, calc_deps=['Q']
            ),
            'lmtd': dc_cp(
                min_val=0, is_result=True, quantity="temperature_difference",
                description="effective logarithmic mean temperature difference |Q|/UA",
                calc=self._calc_lmtd, calc_deps=['Q', 'UA']
            ),
            'UA_char': dc_cc(
                param='m',
                description="heat transfer coefficient lookup table for offdesign"
            ),
            'kA_char': dc_cc(
                param='m',
                description="deprecated, use :code:`UA_char` instead"
            ),
            'Tamb': dc_cp(
                quantity="temperature",
                description="ambient temperature"
            ),
            'dissipative': dc_simple(_val=None),
            'darcy_group': dc_gcp(
                elements=['L', 'ks', 'D'], num_eq_sets=1,
                func=self.darcy_func,
                dependents=self.darcy_dependents,
                description="Darcy-Weißbach equation for pressure loss"
            ),
            'hw_group': dc_gcp(
                elements=['L', 'ks_HW', 'D'], num_eq_sets=1,
                func=self.hazen_williams_func,
                dependents=self.hazen_williams_dependents,
                description="Hazen-Williams equation for pressure loss"
            ),
            'UA_group': dc_gcp(
                elements=['UA', 'Tamb'], num_eq_sets=1,
                func=self.UA_group_func,
                dependents=self.UA_group_dependents,
                description="equation for heat transfer based on ambient temperature and heat transfer coefficient"
            ),
            'kA_group': dc_gcp(
                elements=['kA', 'Tamb'],
                description="deprecated, use :code:`UA_group` instead"
            ),
            'UA_char_group': dc_gcp(
                elements=['UA_char', 'Tamb'], num_eq_sets=1,
                func=self.UA_char_group_func,
                dependents=self.UA_char_group_dependents,
                description="heat transfer from design heat transfer coefficient, modifier lookup table and ambient temperature"
            ),
            'kA_char_group': dc_gcp(
                elements=['kA_char', 'Tamb'],
                description="deprecated, use :code:`UA_char_group` instead"
            )
        }

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
        return ['in1']

    @staticmethod
    def outlets():
        return ['out1']

    @staticmethod
    def powerinlets():
        return ['heat']

    @staticmethod
    def poweroutlets():
        return ['heat']

    @staticmethod
    def heatinlets():
        return ['heat']

    @staticmethod
    def heatoutlets():
        return ['heat']

    def _validate_connections(self):
        super()._validate_connections()
        all_energy = self.power_inl + self.power_outl + self.heat_inl + self.heat_outl
        if len(all_energy) > 1:
            msg = (
                f"Component {self.label} has more than one energy connection. "
                "Connect to exactly one side using one connection type."
            )
            raise TESPyNetworkError(msg)
        if self.power_inl or self.power_outl:
            warnings.warn(
                f"Component {self.label} is connected via PowerConnection. "
                "Please use HeatConnection instead. PowerConnection support for "
                "SimpleHeatExchanger will be removed in a future version.",
                FutureWarning,
                stacklevel=2,
            )

    def _get_energy_connector_location(self):
        """Return (connector, side, val_attr) for the active energy connection."""
        if self.heat_inl:
            return self.heat_inl[0], "inlet"
        if self.heat_outl:
            return self.heat_outl[0], "outlet"
        if self.power_inl:
            return self.power_inl[0], "inlet"
        return self.power_outl[0], "outlet"

    def energy_connector_balance_func(self):
        connector, side = self._get_energy_connector_location()
        energy_flow = (
            -connector.E.val_SI if side == "inlet"
            else connector.E.val_SI
        )
        return energy_flow + self.inl[0].m.val_SI * (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI
        )

    def energy_connector_dependents(self):
        connector, _ = self._get_energy_connector_location()
        return [connector.E, self.inl[0].m, self.outl[0].h, self.inl[0].h]

    def energy_balance_func(self):
        r"""
        Equation for pressure drop calculation.

        Returns
        -------
        residual : float
            Residual value of equation:

            .. math::

                0 =\dot{m}_{in}\cdot\left( h_{out}-h_{in}\right) -\dot{Q}
        """
        return self._calc_Q() - self.Q.val_SI

    def energy_balance_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].h,
            self.outl[0].h
        ]

    def darcy_func(self):
        r"""
        Equation for pressure drop calculation from darcy friction factor.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = p_{in} - p_{out} - \frac{8 \cdot |\dot{m}_{in}| \cdot
                \dot{m}_{in} \cdot \frac{v_{in}+v_{out}}{2} \cdot L \cdot
                \lambda\left(Re, ks, D\right)}{\pi^2 \cdot D^5}\\

                Re = \frac{4 \cdot |\dot{m}_{in}|}{\pi \cdot D \cdot
                \frac{\eta_{in}+\eta_{out}}{2}}\\
                \eta: \text{dynamic viscosity}\\
                v: \text{specific volume}\\
                \lambda: \text{darcy friction factor}
        """
        i = self.inl[0]
        o = self.outl[0]

        if abs(i.m.val_SI) < 1e-4:
            return i.p.val_SI - o.p.val_SI

        visc_i = i.calc_viscosity(T0=i.T.val_SI)
        visc_o = o.calc_viscosity(T0=o.T.val_SI)
        v_i = i.calc_vol(T0=i.T.val_SI)
        v_o = o.calc_vol(T0=o.T.val_SI)
        Re = (
            4 * abs(i.m.val_SI)
            / (math.pi * self.D.val_SI * (visc_i + visc_o) / 2)
        )

        return (
            (i.p.val_SI - o.p.val_SI)
            - 8 * abs(i.m.val_SI) * i.m.val_SI * (v_i + v_o)
            / 2 * self.L.val_SI * dff(Re, self.ks.val_SI, self.D.val_SI)
            / (math.pi ** 2 * self.D.val_SI ** 5)
        )

    def darcy_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].p,
            self.outl[0].h,
        ] + [self.get_attr(element) for element in self.darcy_group.elements]

    def hazen_williams_func(self):
        r"""
        Equation for pressure drop calculation from Hazen-Williams equation.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \left(p_{in} - p_{out} \right) \cdot \left(-1\right)^i -
                \frac{10.67 \cdot |\dot{m}_{in}| ^ {1.852}
                \cdot L}{ks^{1.852} \cdot D^{4.871}} \cdot g \cdot
                \left(\frac{v_{in} + v_{out}}{2}\right)^{0.852}

                i = \begin{cases}
                0 & \dot{m}_{in} \geq 0\\
                1 & \dot{m}_{in} < 0
                \end{cases}

        Note
        ----
        Gravity :math:`g` is set to :math:`9.81 \frac{m}{s^2}`
        """
        i = self.inl[0]
        o = self.outl[0]

        if abs(i.m.val_SI) < 1e-4:
            return i.p.val_SI - o.p.val_SI

        v_i = i.calc_vol(T0=i.T.val_SI)
        v_o = o.calc_vol(T0=o.T.val_SI)

        return (
            math.copysign(i.p.val_SI - o.p.val_SI, i.m.val_SI)
            - (
                10.67 * abs(i.m.val_SI) ** 1.852 * self.L.val_SI /
                (self.ks_HW.val_SI ** 1.852 * self.D.val_SI ** 4.871)
            ) * (9.81 * ((v_i + v_o) / 2) ** 0.852)
        )

    def hazen_williams_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].p,
            self.outl[0].h,
        ] + [self.get_attr(element) for element in self.hw_group.elements]

    def _calculate_td_log(self):
        r"""
        Calculation of mean logarithmic temperature difference.

        For numerical stability: If temperature differences have
        different sign use mean difference to avoid negative logarithm.

        Returns
        -------
        deltaT_log : float
            Mean logarithmic temperature difference.

            .. math::

                \Delta T_{log} = \begin{cases}
                \frac{T_{in}-T_{out}}{\ln{\frac{T_{in}-T_{amb}}
                {T_{out}-T_{amb}}}} & T_{in} > T_{out} \\
                \frac{T_{out}-T_{in}}{\ln{\frac{T_{out}-T_{amb}}
                {T_{in}-T_{amb}}}} & T_{in} < T_{out}\\
                0 & T_{in} = T_{out}
                \end{cases}

                T_{amb}: \text{ambient temperature}
        """
        i = self.inl[0]
        o = self.outl[0]

        ttd_1 = i.calc_T() - self.Tamb.val_SI
        ttd_2 = o.calc_T() - self.Tamb.val_SI

        # For numerical stability: If temperature differences have
        # different sign use mean difference to avoid negative logarithm.
        if (ttd_1 / ttd_2) < 0:
            if ttd_1 > 0:
                if o.h.is_var and self.it < 10:
                    h_out = h_mix_pT(
                        o.p.val_SI,
                        self.Tamb.val_SI + 0.0001,
                        o.fluid_data,
                        o.mixing_rule
                    )
                    o.h.set_reference_val_SI(h_out)
                ttd_2 = 0.1
            elif ttd_1 < 0:
                if o.h.is_var and self.it < 10:
                    h_out = h_mix_pT(
                        o.p.val_SI,
                        self.Tamb.val_SI - 0.0001,
                        o.fluid_data,
                        o.mixing_rule
                    )
                    o.h.set_reference_val_SI(h_out)
                ttd_2 = -0.1

        if round(ttd_1, 6) == round(ttd_2, 6):
            td_log = ttd_2
        elif ttd_1 > ttd_2:
            td_log = (ttd_1 - ttd_2) / math.log(ttd_1 / ttd_2)
        else:
            td_log = (ttd_2 - ttd_1) / math.log(ttd_2 / ttd_1)

        return td_log

    def UA_group_func(self):
        r"""
        Calculate heat transfer from heat transfer coefficient.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_{in} \cdot \left( h_{out} - h_{in}\right) +
                UA \cdot \Delta T_{log}
        """
        i = self.inl[0]
        o = self.outl[0]
        Q = i.m.val_SI * (o.h.val_SI - i.h.val_SI)
        ttd_1 = i.calc_T() - self.Tamb.val_SI
        ttd_2 = o.calc_T() - self.Tamb.val_SI
        if ttd_1 * ttd_2 <= 0:
            return Q + self.UA.val_SI * ttd_2
        return Q + self.UA.val_SI * self._calculate_td_log()

    def UA_group_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].p,
            self.outl[0].h,
            self.UA
        ]

    def UA_char_group_func(self):
        r"""
        Calculate heat transfer from heat transfer coefficient characteristic.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_{in} \cdot \left( h_{out} - h_{in}\right) +
                UA_{design} \cdot f_{UA} \cdot \Delta T_{log}
        """
        i = self.inl[0]
        o = self.outl[0]
        p = self.UA_char.param

        expr = self.get_char_expr(p, **self.UA_char.char_params)
        fkA = 2 / (1 + 1 / self.UA_char.char_func.evaluate(expr))

        Q = i.m.val_SI * (o.h.val_SI - i.h.val_SI)
        ttd_1 = i.calc_T() - self.Tamb.val_SI
        ttd_2 = o.calc_T() - self.Tamb.val_SI
        if ttd_1 * ttd_2 <= 0:
            return Q + self.UA.design * fkA * ttd_2
        return Q + self.UA.design * fkA * self._calculate_td_log()

    def UA_char_group_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].p,
            self.outl[0].h,
        ]

    def convergence_check(self):
        if self.UA_group.is_set:
            i = self.inl[0]
            o = self.outl[0]
            T_in = i.calc_T()
            T_out = o.calc_T()
            if T_in > self.Tamb.val_SI:
                if T_out < self.Tamb.val_SI:
                    if o.h.is_var:
                        h_out = h_mix_pT(
                            o.p.val_SI,
                            self.Tamb.val_SI + 0.0001,
                            o.fluid_data,
                            o.mixing_rule
                        )
                        o.h.set_reference_val_SI(h_out)
            elif T_in < self.Tamb.val_SI:
                if T_out > self.Tamb.val_SI:
                    if o.h.is_var:
                        h_out = h_mix_pT(
                            o.p.val_SI,
                            self.Tamb.val_SI - 0.0001,
                            o.fluid_data,
                            o.mixing_rule
                        )
                        o.h.set_reference_val_SI(h_out)

    def initialise_source(self, c, key):
        r"""
        Return a starting value for pressure and enthalpy the outlets.

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
        c : tespy.connections.connection.Connection
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
        super().calc_parameters()
        if "UA" not in self.parameters:
            return
        self.UA.is_result = self.Tamb.is_set
        self.kA.is_result = self.Tamb.is_set
        self.lmtd.is_result = self.Tamb.is_set

    def entropy_balance(self):
        r"""
        Calculate entropy balance of a simple heat exchanger.

        The allocation of the entropy streams due to heat exchanged and due to
        irreversibility is performed by solving for T:

        .. math::

            h_\text{out} - h_\text{in} = \int_\text{out}^\text{in}
            v \cdot dp - \int_\text{out}^\text{in} T \cdot ds

        As solving :math:`\int_\text{out}^\text{in} v \cdot dp` for non
        isobaric processes would require perfect process knowledge (the path)
        on how specific volume and pressure change throught the component, the
        heat transfer is split into three separate virtual processes:

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

            \text{S\_Q}=\dot{m} \cdot \left(s_\text{out*}-s_\text{in*}
            \right)\\
            \text{S\_irr}=\dot{m} \cdot \left(s_\text{out}-s_\text{in}
            \right) - \text{S\_Q}\\
            \text{T\_mQ}=\frac{\dot{Q}}{\text{S\_Q}}
        """
        i = self.inl[0]
        o = self.outl[0]

        p1_star = i.p.val_SI * (o.p.val_SI / i.p.val_SI) ** 0.5
        s1_star = s_mix_ph(
            p1_star, i.h.val_SI, i.fluid_data, i.mixing_rule, T0=i.T.val_SI
        )
        s2_star = s_mix_ph(
            p1_star, o.h.val_SI, o.fluid_data, o.mixing_rule, T0=o.T.val_SI
        )
        self.S_Q = i.m.val_SI * (s2_star - s1_star)
        self.S_irr = i.m.val_SI * (o.s.val_SI - i.s.val_SI) - self.S_Q
        self.T_mQ = (o.h.val_SI - i.h.val_SI) / (s2_star - s1_star)

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
