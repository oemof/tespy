# -*- coding: utf-8

"""Module of class SimpleHeatExchanger.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/heat_exchangers/simple.py

SPDX-License-Identifier: MIT
"""

import math

import numpy as np

from tespy.components.component import Component
from tespy.components.component import component_registry
from tespy.tools import logger
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
from tespy.tools.fluid_properties import s_mix_ph
from tespy.tools.fluid_properties.helpers import darcy_friction_factor as dff
from tespy.tools.helpers import _numeric_deriv
from tespy.tools.helpers import convert_to_SI


@component_registry
class SimpleHeatExchanger(Component):
    r"""
    A basic heat exchanger representing a heat source or heat sink.

    The component SimpleHeatExchanger is the parent class for the components:

    - :py:class:`tespy.components.heat_exchangers.solar_collector.SolarCollector`
    - :py:class:`tespy.components.heat_exchangers.parabolic_trough.ParabolicTrough`
    - :py:class:`tespy.components.piping.pipe.Pipe`

    **Mandatory Equations**

    - :py:meth:`tespy.components.component.Component.fluid_func`
    - :py:meth:`tespy.components.component.Component.mass_flow_func`

    **Optional Equations**

    - :py:meth:`tespy.components.component.Component.pr_func`
    - :py:meth:`tespy.components.component.Component.zeta_func`
    - :py:meth:`tespy.components.component.Component.dp_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.energy_balance_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.darcy_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.hazen_williams_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.kA_group_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.kA_char_group_func`

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

    zeta : float, dict, :code:`"var"`
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

    kA : float, dict, :code:`"var"`
        Area independent heat transfer coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char : tespy.tools.characteristics.CharLine, dict
        Characteristic line for heat transfer coefficient.

    Tamb : float, dict
        Ambient temperature, provide parameter in network's temperature unit.

    kA_group : str, dict
        Parametergroup for heat transfer calculation from ambient temperature
        and area independent heat transfer coefficient kA.

    Example
    -------
    The SimpleHeatExchanger can be used as a sink or source of heat. This
    component does not simulate the secondary side of the heat exchanger. It
    is possible to calculate the pressure ratio with the Darcy-Weisbach
    equation or in case of liquid water use the Hazen-Williams equation.
    Also, given ambient temperature and the heat transfer coeffiecient, it is
    possible to predict heat transfer.

    >>> from tespy.components import Sink, Source, SimpleHeatExchanger
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> import os
    >>> nw = Network()
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg', iterinfo=False)
    >>> so1 = Source('source 1')
    >>> si1 = Sink('sink 1')
    >>> heat_sink = SimpleHeatExchanger('heat sink')
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
    >>> nw.save('tmp.json')
    >>> round(heat_sink.Q.val, 0)
    -52581.0
    >>> round(heat_sink.kA.val, 0)
    321.0
    >>> inc.set_attr(m=1.25)
    >>> nw.solve('offdesign', design_path='tmp.json')
    >>> round(heat_sink.Q.val, 0)
    -56599.0
    >>> round(outg.T.val, 1)
    156.9
    >>> inc.set_attr(m=0.75)
    >>> nw.solve('offdesign', design_path='tmp.json')
    >>> round(heat_sink.Q.val, 1)
    -47275.8
    >>> round(outg.T.val, 1)
    140.0
    >>> os.remove('tmp.json')

    Use of the PowerConnection
    --------------------------

    Single sided heat exchangers can also connect to a
    :code:`PowerConnection`. In this case the heat exchanger represents a heat
    source from the perspective of the :code:`PowerConnection`, meaning, we
    need a :code:`PowerSink` component as target of that connection.
    To utilize the component in this connects, it is required to set the
    :code:`power_connector_location`, in this case it should be the
    :code:`'outlet'`.

    >>> from tespy.connections import PowerConnection
    >>> from tespy.components import PowerSink
    >>> ambient = PowerSink('ambient heat dissipation')
    >>> heat_sink.set_attr(power_connector_location='outlet')

    Then we can create and add the :code:`PowerConnection`. We run a new
    design calculation, because the old design case did not inlcude the
    :code:`PowerConnection`. The energy value will be identical to the heat
    transfer of the pipe.

    >>> h1 = PowerConnection(heat_sink, 'heat', ambient, 'power', label='h1')
    >>> nw.add_conns(h1)
    >>> nw.solve('design')
    >>> round(h1.E.val) == round(-heat_sink.Q.val)
    True
    """

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        if len(self.power_inl + self.power_outl) > 0:
            constraints["energy_connector_balance"] = dc_cmc(**{
                "func": self.energy_connector_balance_func,
                "dependents": self.energy_connector_dependents,
                "num_eq_sets": 1
            })

        return constraints

    def get_parameters(self):
        return {
            'power_connector_location': dc_simple(),
            'Q': dc_cp(
                num_eq_sets=1,
                func=self.energy_balance_func,
                dependents=self.energy_balance_dependents,
            ),
            'pr': dc_cp(
                min_val=1e-4, max_val=1, num_eq_sets=1,
                func=self.pr_func,
                structure_matrix=self.pr_structure_matrix,
                dependents=self.pr_dependents,
                func_params={'pr': 'pr'},
            ),
            'dp': dc_cp(
                min_val=0, max_val=1e15, num_eq_sets=1,
                func=self.dp_func,
                dependents=self.dp_dependents,
                structure_matrix=self.dp_structure_matrix,
                func_params={'dp': 'dp'}
            ),
            'zeta': dc_cp(
                min_val=0, max_val=1e15, num_eq_sets=1,
                func=self.zeta_func,
                dependents=self.zeta_dependents,
                func_params={'zeta': 'zeta'}
            ),
            'D': dc_cp(min_val=1e-2, max_val=2, d=1e-4),
            'L': dc_cp(min_val=1e-1, d=1e-3),
            'ks': dc_cp(_val=1e-4, min_val=1e-7, max_val=1e-3, d=1e-8),
            'ks_HW': dc_cp(_val=10, min_val=1e-1, max_val=1e3, d=1e-2),
            'kA': dc_cp(min_val=0, d=1),
            'kA_char': dc_cc(param='m'), 'Tamb': dc_cp(),
            'dissipative': dc_simple(_val=None),
            'darcy_group': dc_gcp(
                elements=['L', 'ks', 'D'], num_eq_sets=1,
                func=self.darcy_func,
                dependents=self.darcy_dependents
            ),
            'hw_group': dc_gcp(
                elements=['L', 'ks_HW', 'D'], num_eq_sets=1,
                func=self.hazen_williams_func,
                dependents=self.hazen_williams_dependents
            ),
            'kA_group': dc_gcp(
                elements=['kA', 'Tamb'], num_eq_sets=1,
                func=self.kA_group_func,
                dependents=self.kA_group_dependents
            ),
            'kA_char_group': dc_gcp(
                elements=['kA_char', 'Tamb'], num_eq_sets=1,
                func=self.kA_char_group_func,
                dependents=self.kA_char_group_dependents
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

    def powerinlets(self):
        if self.power_connector_location.val == "inlet":
            return ['heat']
        else:
            return []

    def poweroutlets(self):
        if self.power_connector_location.val == "outlet":
            return ['heat']
        else:
            return []

    def _preprocess(self, row_idx):
        self.Tamb.val_SI = convert_to_SI('T', self.Tamb.val, self.inl[0].T.unit)

        if self.dp.is_set:
            self.dp.val_SI = convert_to_SI('p', self.dp.val, self.inl[0].p.unit)

        super()._preprocess(row_idx)

    def _get_power_connector_location(self):
        if self.power_connector_location.val == "inlet":
            return self.power_inl[0]
        else:
            return self.power_outl[0]

    def energy_connector_balance_func(self):
        connector = self._get_power_connector_location()
        if self.power_connector_location.val == "inlet":
            energy_flow = -connector.E.val_SI
        else:
            energy_flow = connector.E.val_SI

        return energy_flow + self.inl[0].m.val_SI * (
                self.outl[0].h.val_SI - self.inl[0].h.val_SI
            )

    def energy_connector_dependents(self):
        connector = self._get_power_connector_location()
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
        return self.inl[0].m.val_SI * (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI
        ) - self.Q.val

    def energy_balance_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].h,
            self.outl[0].h,
            self.Q,
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
        Re = 4 * abs(i.m.val_SI) / (math.pi * self.D.val * (visc_i + visc_o) / 2)

        return (
            (i.p.val_SI - o.p.val_SI)
            - 8 * abs(i.m.val_SI) * i.m.val_SI * (v_i + v_o)
            / 2 * self.L.val * dff(Re, self.ks.val, self.D.val)
            / (math.pi ** 2 * self.D.val ** 5)
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
                10.67 * abs(i.m.val_SI) ** 1.852 * self.L.val /
                (self.ks_HW.val ** 1.852 * self.D.val ** 4.871)
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
        """
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
            td_log = (ttd_2 + ttd_1) / 2
        elif round(ttd_1, 6) == round(ttd_2, 6):
            td_log = ttd_2
        elif ttd_1 > ttd_2:
            td_log = (ttd_1 - ttd_2) / math.log(ttd_1 / ttd_2)
        else:
            td_log = (ttd_2 - ttd_1) / math.log(ttd_2 / ttd_1)

        return td_log

    def kA_group_func(self):
        r"""
        Calculate heat transfer from heat transfer coefficient.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_{in} \cdot \left( h_{out} - h_{in}\right) +
                kA \cdot \Delta T_{log}

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
        td_log = self._calculate_td_log()
        return i.m.val_SI * (o.h.val_SI - i.h.val_SI) + self.kA.val * td_log

    def kA_group_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].p,
            self.outl[0].h,
            self.kA
        ]

    def kA_char_group_func(self):
        r"""
        Calculate heat transfer from heat transfer coefficient characteristic.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_{in} \cdot \left( h_{out} - h_{in}\right) +
                kA_{design} \cdot f_{kA} \cdot \Delta T_{log}

                \Delta T_{log} = \begin{cases}
                \frac{T_{in}-T_{out}}{\ln{\frac{T_{in}-T_{amb}}
                {T_{out}-T_{amb}}}} & T_{in} > T_{out} \\
                \frac{T_{out}-T_{in}}{\ln{\frac{T_{out}-T_{amb}}
                {T_{in}-T_{amb}}}} & T_{in} < T_{out}\\
                0 & T_{in} = T_{out}
                \end{cases}

                f_{kA} = \frac{2}{1 + \frac{1}{f\left( expr\right)}}

                T_{amb}: \text{ambient temperature}

        Note
        ----
        For standard function of f\ :subscript:`kA` \ see module
        :py:mod:`tespy.data`.
        """
        i = self.inl[0]
        o = self.outl[0]
        p = self.kA_char.param

        expr = self.get_char_expr(p, **self.kA_char.char_params)
        fkA = 2 / (1 + 1 / self.kA_char.char_func.evaluate(expr))

        td_log = self._calculate_td_log()

        return (
            i.m.val_SI * (o.h.val_SI - i.h.val_SI)
            + self.kA.design * fkA * td_log
        )

    def kA_char_group_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].p,
            self.outl[0].h,
        ]

    def bus_func(self, bus):
        r"""
        Calculate the value of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus.Bus
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
        return self.inl[0].m.val_SI * (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI)

    def bus_deriv(self, bus):
        r"""
        Calculate partial derivatives of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus.Bus
            TESPy bus object.

        Returns
        -------
        deriv : ndarray
            Matrix of partial derivatives.
        """
        f = self.calc_bus_value
        if self.inl[0].m.is_var:
            if self.inl[0].m.J_col not in bus.jacobian:
                bus.jacobian[self.inl[0].m.J_col] = 0
            bus.jacobian[self.inl[0].m.J_col] -= _numeric_deriv(self.inl[0].m._reference_container, f, bus=bus)

        if self.inl[0].h.is_var:
            if self.inl[0].h.J_col not in bus.jacobian:
                bus.jacobian[self.inl[0].h.J_col] = 0
            bus.jacobian[self.inl[0].h.J_col] -= _numeric_deriv(self.inl[0].h._reference_container, f, bus=bus)

        if self.outl[0].h.is_var:
            if self.outl[0].h.J_col not in bus.jacobian:
                bus.jacobian[self.outl[0].h.J_col] = 0
            bus.jacobian[self.outl[0].h.J_col] -= _numeric_deriv(self.outl[0].h._reference_container, f, bus=bus)

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
        i = self.inl[0]
        o = self.outl[0]

        self.Q.val = i.m.val_SI * (o.h.val_SI - i.h.val_SI)
        self.pr.val = o.p.val_SI / i.p.val_SI
        self.dp.val_SI = i.p.val_SI - o.p.val_SI
        self.dp.val = i.p.val - o.p.val
        self.zeta.val = self.calc_zeta(i, o)

        if self.Tamb.is_set:
            ttd_1 = i.T.val_SI - self.Tamb.val_SI
            ttd_2 = o.T.val_SI - self.Tamb.val_SI
            if ttd_1 / ttd_2 < 0:
                self.kA.val = np.nan
            else:
                self.kA.val = abs(self.Q.val / self._calculate_td_log())

            self.kA.is_result = True
        else:
            self.kA.is_result = False

    def entropy_balance(self):
        r"""
        Calculate entropy balance of a simple heat exchanger.

        The allocation of the entropy streams due to heat exchanged and due to
        irreversibility is performed by solving for T:

        .. math::

            h_\mathrm{out} - h_\mathrm{in} = \int_\mathrm{out}^\mathrm{in}
            v \cdot dp - \int_\mathrm{out}^\mathrm{in} T \cdot ds

        As solving :math:`\int_\mathrm{out}^\mathrm{in} v \cdot dp` for non
        isobaric processes would require perfect process knowledge (the path)
        on how specific volume and pressure change throught the component, the
        heat transfer is splitted into three separate virtual processes:

        - in->in*: decrease pressure to
          :math:`p_\mathrm{in*}=p_\mathrm{in}\cdot\sqrt{\frac{p_\mathrm{out}}{p_\mathrm{in}}}`
          without changing enthalpy.
        - in*->out* transfer heat without changing pressure.
          :math:`h_\mathrm{out*}-h_\mathrm{in*}=h_\mathrm{out}-h_\mathrm{in}`
        - out*->out decrease pressure to outlet pressure :math:`p_\mathrm{out}`
          without changing enthalpy.

        Note
        ----
        The entropy balance makes the follwing parameter available:

        .. math::

            \text{S\_Q}=\dot{m} \cdot \left(s_\mathrm{out*}-s_\mathrm{in*}
            \right)\\
            \text{S\_irr}=\dot{m} \cdot \left(s_\mathrm{out}-s_\mathrm{in}
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

    def exergy_balance(self, T0):
        r"""
        Calculate exergy balance of a simple heat exchanger.

        The exergy of heat is calculated by allocation of thermal and
        mechanical share of exergy in the physical exergy. Depending on the
        temperature levels at the inlet and outlet of the heat exchanger as
        well as the direction of heat transfer (input or output) fuel and
        product exergy are calculated as follows.

        Parameters
        ----------
        T0 : float
            Ambient temperature T0 / K.

        Note
        ----
        If the fluid transfers heat to the ambient, you can specify
        :code:`mysimpleheatexchanger.set_attr(dissipative=False)` if you do
        NOT want the exergy production nan (only applicable in case
        :math:`\dot{Q}<0`).

        .. math ::

            \dot{E}_\mathrm{P} =
            \begin{cases}
            \begin{cases}
            \begin{cases}
            \text{not defined (nan)} & \text{if dissipative}\\
            \dot{E}_\mathrm{in}^\mathrm{T} - \dot{E}_\mathrm{out}^\mathrm{T} &
            \text{else}\\
            \end{cases}
            & T_\mathrm{in}, T_\mathrm{out} \geq T_0\\
            \dot{E}_\mathrm{out}^\mathrm{T}
            & T_\mathrm{in} \geq T_0 > T_\mathrm{out}\\
            \dot{E}_\mathrm{out}^\mathrm{T} - \dot{E}_\mathrm{in}^\mathrm{T}
            & T_0 \geq T_\mathrm{in}, T_\mathrm{out}\\
            \end{cases} & \dot{Q} < 0\\

            \begin{cases}
            \dot{E}_\mathrm{out}^\mathrm{PH} - \dot{E}_\mathrm{in}^\mathrm{PH}
            & T_\mathrm{in}, T_\mathrm{out} \geq T_0\\
            \dot{E}_\mathrm{in}^\mathrm{T} + \dot{E}_\mathrm{out}^\mathrm{T}
            & T_\mathrm{out} > T_0 \geq T_\mathrm{in}\\
            \dot{E}_\mathrm{in}^\mathrm{T} - \dot{E}_\mathrm{out}^\mathrm{T} +
            \dot{E}_\mathrm{out}^\mathrm{M} - \dot{E}_\mathrm{in}^\mathrm{M} +
            & T_0 \geq T_\mathrm{in}, T_\mathrm{out}\\
            \end{cases} & \dot{Q} > 0\\
            \end{cases}

            \dot{E}_\mathrm{F} =
            \begin{cases}
            \begin{cases}
            \dot{E}_\mathrm{in}^\mathrm{PH} - \dot{E}_\mathrm{out}^\mathrm{PH}
            & T_\mathrm{in}, T_\mathrm{out} \geq T_0\\
            \dot{E}_\mathrm{in}^\mathrm{T} + \dot{E}_\mathrm{in}^\mathrm{M} +
            \dot{E}_\mathrm{out}^\mathrm{T} - \dot{E}_\mathrm{out}^\mathrm{M}
            & T_\mathrm{in} \geq T_0 > T_\mathrm{out}\\
            \dot{E}_\mathrm{out}^\mathrm{T} - \dot{E}_\mathrm{in}^\mathrm{T} +
            \dot{E}_\mathrm{in}^\mathrm{M} - \dot{E}_\mathrm{out}^\mathrm{M} +
            & T_0 \geq T_\mathrm{in}, T_\mathrm{out}\\
            \end{cases} & \dot{Q} < 0\\

            \begin{cases}
            \dot{E}_\mathrm{out}^\mathrm{T} - \dot{E}_\mathrm{in}^\mathrm{T}
            & T_\mathrm{in}, T_\mathrm{out} \geq T_0\\
            \dot{E}_\mathrm{in}^\mathrm{T} + \dot{E}_\mathrm{in}^\mathrm{M} -
            \dot{E}_\mathrm{out}^\mathrm{M}
            & T_\mathrm{out} > T_0 \geq T_\mathrm{in}\\
            \dot{E}_\mathrm{in}^\mathrm{T}-\dot{E}_\mathrm{out}^\mathrm{T}
            & T_0 \geq T_\mathrm{in}, T_\mathrm{out}\\
            \end{cases} & \dot{Q} > 0\\
            \end{cases}

            \dot{E}_\mathrm{bus} =
            \begin{cases}
            \begin{cases}
            \dot{E}_\mathrm{P} & \text{other cases}\\
            \dot{E}_\mathrm{in}^\mathrm{T}
            & T_\mathrm{in} \geq T_0 > T_\mathrm{out}\\
            \end{cases} & \dot{Q} < 0\\
            \dot{E}_\mathrm{F} & \dot{Q} > 0\\
            \end{cases}
        """
        if self.dissipative.val is None:
            self.dissipative.val = True
            msg = (
                "In a future version of TESPy, the dissipative property must "
                "explicitly be set to True or False in the context of the "
                f"exergy analysis for component {self.label}."
            )
            logger.warning(msg)
        if self.Q.val < 0:
            if self.inl[0].T.val_SI >= T0 and self.outl[0].T.val_SI >= T0:
                if self.dissipative.val:
                    self.E_P = np.nan
                else:
                    self.E_P = self.inl[0].Ex_therm - self.outl[0].Ex_therm
                self.E_F = self.inl[0].Ex_physical - self.outl[0].Ex_physical
                self.E_bus = {
                    "chemical": 0, "physical": 0, "massless": self.E_P
                }
            elif self.inl[0].T.val_SI >= T0 and self.outl[0].T.val_SI < T0:
                self.E_P = self.outl[0].Ex_therm
                self.E_F = self.inl[0].Ex_therm + self.outl[0].Ex_therm + (
                    self.inl[0].Ex_mech - self.outl[0].Ex_mech)
                self.E_bus = {
                    "chemical": 0, "physical": 0,
                    "massless": self.inl[0].Ex_therm + self.outl[0].Ex_therm
                }
            elif self.inl[0].T.val_SI <= T0 and self.outl[0].T.val_SI <= T0:
                self.E_P = self.outl[0].Ex_therm - self.inl[0].Ex_therm
                self.E_F = self.outl[0].Ex_therm - self.outl[0].Ex_therm + (
                    self.inl[0].Ex_mech - self.outl[0].Ex_mech)
                self.E_bus = {
                    "chemical": 0, "physical": 0, "massless": self.E_P
                }
            else:
                msg = ('Exergy balance of simple heat exchangers, where '
                       'outlet temperature is higher than inlet temperature '
                       'with heat extracted is not implmented.')
                logger.warning(msg)
                self.E_P = np.nan
                self.E_F = np.nan
                self.E_bus = {
                    "chemical": np.nan, "physical": np.nan, "massless": np.nan
                }
        elif self.Q.val > 0:
            if self.inl[0].T.val_SI >= T0 - 1e-6 and self.outl[0].T.val_SI >= T0 - 1e-6:
                self.E_P = self.outl[0].Ex_physical - self.inl[0].Ex_physical
                self.E_F = self.outl[0].Ex_therm - self.inl[0].Ex_therm
                self.E_bus = {
                    "chemical": 0, "physical": 0, "massless": self.E_F
                }
            elif self.inl[0].T.val_SI <= T0 and self.outl[0].T.val_SI > T0:
                self.E_P = self.outl[0].Ex_therm + self.inl[0].Ex_therm
                self.E_F = self.inl[0].Ex_therm + (
                    self.inl[0].Ex_mech - self.outl[0].Ex_mech)
                self.E_bus = {
                    "chemical": 0, "physical": 0,
                    "massless": self.inl[0].Ex_therm
                }
            elif self.inl[0].T.val_SI < T0 and self.outl[0].T.val_SI < T0:
                if self.dissipative.val:
                    self.E_P = np.nan
                else:
                    self.E_P = self.inl[0].Ex_therm - self.outl[0].Ex_therm + (
                        self.outl[0].Ex_mech - self.inl[0].Ex_mech
                    )
                self.E_F = self.inl[0].Ex_therm - self.outl[0].Ex_therm
                self.E_bus = {
                    "chemical": 0, "physical": 0, "massless": self.E_F
                }
            else:
                msg = ('Exergy balance of simple heat exchangers, where '
                       'inlet temperature is higher than outlet temperature '
                       'with heat injected is not implmented.')
                logger.warning(msg)
                self.E_P = np.nan
                self.E_F = np.nan
                self.E_bus = {
                    "chemical": np.nan, "physical": np.nan, "massless": self.E_F
                }
        else:
            # fully dissipative
            self.E_P = np.nan
            self.E_F = self.inl[0].Ex_physical - self.outl[0].Ex_physical
            self.E_bus = {
                "chemical": np.nan, "physical": np.nan, "massless": np.nan
            }

        if np.isnan(self.E_P):
            self.E_D = self.E_F
        else:
            self.E_D = self.E_F - self.E_P
        self.epsilon = self._calc_epsilon()

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
