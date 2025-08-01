# -*- coding: utf-8

"""Module of class WaterElectrolyzer.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/reactors/water_electrolyzer.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import Component
from tespy.components.component import component_registry
from tespy.tools import logger
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.fluid_properties import dT_mix_dph
from tespy.tools.fluid_properties import dT_mix_pdh
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.helpers import _numeric_deriv
from tespy.tools.helpers import convert_to_SI


@component_registry
class WaterElectrolyzer(Component):
    r"""
    The water electrolyzer produces hydrogen and oxygen from water and power.

    **Mandatory Equations**

    - :py:meth:`tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.fluid_func`
    - :py:meth:`tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.mass_flow_func`
    - :py:meth:`tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.reactor_pressure_func`
    - :py:meth:`tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.energy_balance_func`
    - :py:meth:`tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.gas_temperature_func`

    **Optional Equations**

    - cooling loop:

      - :py:meth:`tespy.components.component.Component.zeta_func`
      - :py:meth:`tespy.components.component.Component.pr_func`

    - :py:meth:`tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.eta_func`
    - :py:meth:`tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.eta_char_func`
    - :py:meth:`tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.heat_func`
    - :py:meth:`tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.specific_energy_consumption_func`

    Inlets/Outlets

    - in1 (cooling inlet), in2 (feed water inlet)
    - out1 (cooling outlet), out2 (oxygen outlet), out3 (hydrogen outlet)

    Image

    .. image:: /api/_images/WaterElectrolyzer.svg
       :alt: flowsheet of the water electrolyzer
       :align: center
       :class: only-light

    .. image:: /api/_images/WaterElectrolyzer_darkmode.svg
       :alt: flowsheet of the water electrolyzer
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

    P : float, dict, :code:`"var"`
        Power input, :math:`P/\text{W}`.

    Q : float, dict
        Heat output of cooling, :math:`Q/\text{W}`

    e : float, dict, :code:`"var"`
        Electrolysis specific energy consumption,
        :math:`e/(\text{J}/\text{m}^3)`.

    eta : float, dict
        Electrolysis efficiency (referring to H2 higher heating value),
        :math:`\eta/1`.

    eta_char : tespy.tools.characteristics.CharLine, dict
        Electrolysis efficiency characteristic line.

    pr : float, dict, :code:`"var"`
        Cooling loop pressure ratio, :math:`pr/1`.

    zeta : float, dict, :code:`"var"`
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

    >>> from tespy.components import (Sink, Source, Compressor,
    ... WaterElectrolyzer)
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> import os
    >>> nw = Network(T_unit='C', p_unit='bar', v_unit='l / s', iterinfo=False)
    >>> fw = Source('feed water')
    >>> oxy = Sink('oxygen sink')
    >>> hydro = Sink('hydrogen sink')
    >>> cw_cold = Source('cooling water source')
    >>> cw_hot = Sink('cooling water sink')
    >>> comp = Compressor('compressor', eta_s=0.9)
    >>> el = WaterElectrolyzer('electrolyzer')

    The electrolyzer should produce 100 l/s of hydrogen at an operating
    pressure of 10 bars and an outlet temperature of 50 °C. The fluid
    composition needs to be specified for the cooling liquid only. The storage
    pressure is 25 bars. The electrolysis efficiency is at 80 % and the
    compressor isentropic efficiency at 85 %. After designing the plant the
    offdesign electrolysis efficiency is predicted by the characteristic line.
    The default characteristic line can be found here: :py:mod:`tespy.data`.

    >>> fw_el = Connection(fw, 'out1', el, 'in2')
    >>> el_o = Connection(el, 'out2', oxy, 'in1')
    >>> el_cmp = Connection(el, 'out3', comp, 'in1')
    >>> cmp_h = Connection(comp, 'out1', hydro, 'in1')
    >>> cw_el = Connection(cw_cold, 'out1', el, 'in1')
    >>> el_cw = Connection(el, 'out1', cw_hot, 'in1')
    >>> nw.add_conns(fw_el, el_o, el_cmp, cmp_h, cw_el, el_cw)
    >>> fw_el.set_attr(p=10, T=15)
    >>> cw_el.set_attr(p=5, T=15, fluid={'H2O': 1})
    >>> el_cw.set_attr(T=45)
    >>> cmp_h.set_attr(p=25)
    >>> el_cmp.set_attr(v=100, T=50)
    >>> el.set_attr(eta=0.8, pr=0.99, design=['eta', 'pr'],
    ... offdesign=['eta_char', 'zeta'])
    >>> comp.set_attr(eta_s=0.85)
    >>> nw.solve('design')
    >>> nw.save('tmp.json')
    >>> round(el.e0 / el.P.val * el_cmp.m.val_SI, 1)
    0.8
    >>> P_design = el.P.val
    >>> round(P_design / 1e6, 1)
    13.2
    >>> nw.solve('offdesign', design_path='tmp.json')
    >>> round(el.eta.val, 1)
    0.8
    >>> el_cmp.set_attr(v=None)
    >>> el.set_attr(P=P_design * 0.2)
    >>> nw.solve('offdesign', design_path='tmp.json')
    >>> round(el.eta.val, 2)
    0.84
    >>> os.remove('tmp.json')
    """

    def get_parameters(self):
        return {
            'P': dc_cp(min_val=0),
            'Q': dc_cp(
                max_val=0, num_eq_sets=1,
                func=self.heat_func,
                dependents=self.heat_dependents
            ),
            'pr': dc_cp(
                max_val=1, num_eq_sets=1,
                structure_matrix=self.pr_structure_matrix,
                func=self.pr_func,
                func_params={'pr': 'pr'}
            ),
            'dp': dc_cp(
                min_val=0,
                structure_matrix=self.dp_structure_matrix,
                func=self.dp_func,
                num_eq_sets=1,
                func_params={"inconn": 0, "outconn": 0, "dp": "dp"}
            ),
            'zeta': dc_cp(
                min_val=0,
                num_eq_sets=1,
                dependents=self.zeta_dependents,
                func=self.zeta_func,
                func_params={'zeta': 'zeta'}
            ),
            'eta': dc_cp(
                min_val=0, max_val=1, num_eq_sets=1,
                func=self.eta_func,
                dependents=self.eta_dependents
            ),
            'eta_char': dc_cc(
                deriv=self.eta_char_deriv,
                func=self.eta_char_func,
                dependents=self.eta_char_dependents,
                num_eq_sets=1,
                param='m_out',
                char_params={'type': 'rel', 'outconn': 2}
            ),
            'e': dc_cp(
                min_val=0, num_eq_sets=1,
                func=self.specific_energy_func,
                dependents=self.specific_energy_dependents
            )
        }

    def get_mandatory_constraints(self):
        return {
            'mass_flow_constraints': dc_cmc(**{
                'func': self.reactor_mass_flow_func,
                'deriv': self.reactor_mass_flow_deriv,
                'constant_deriv': True,
                'num_eq_sets': 2
            }),
            'cooling_mass_flow_constraints': dc_cmc(**{
                'func': self.cooling_mass_flow_func,
                'structure_matrix': self.cooling_mass_flow_structure_matrix,
                'num_eq_sets': 1
            }),
            'cooling_fluid_constraints': dc_cmc(**{
                'func': self.cooling_fluid_func,
                'structure_matrix': self.cooling_fluid_structure_matrix,
                'num_eq_sets': 1
            }),
            'energy_balance_constraints': dc_cmc(**{
                'func': self.energy_balance_func,
                'deriv': self.energy_balance_deriv,
                'constant_deriv': False,
                'num_eq_sets': 1
            }),
            'reactor_pressure_constraints': dc_cmc(**{
                'func': self.reactor_pressure_func,
                'structure_matrix': self.reactor_pressure_structure_matrix,
                'constant_deriv': True,
                'num_eq_sets': 2
            }),
            'gas_temperature_constraints': dc_cmc(**{
                'func': self.gas_temperature_func,
                'dependents': self.gas_temperature_dependents,
                'constant_deriv': False,
                'num_eq_sets': 1
            })
        }

    @staticmethod
    def get_bypass_constraints():
        return {}

    @staticmethod
    def inlets():
        return ['in1', 'in2']

    @staticmethod
    def outlets():
        return ['out1', 'out2', 'out3']

    def _add_missing_fluids(self, connections):
        if self.inl[1] in connections:
            return ["H2O"]
        elif self.outl[1] in connections:
            return ["O2"]
        elif self.outl[2] in connections:
            return ["H2"]
        else:
            return super()._add_missing_fluids(connections)

    def get_variables(self):
        if not self.P.is_set:
            self.set_attr(P='var')
            msg = (
                f'The power input of water electrolyzers ({self.label}) must '
                'either be a set value or part of the system variables. Since '
                'it has not been set to a fixed value it will be added to the '
                'system\'s variables.'
            )
            logger.info(msg)
        return super().get_variables()

    def _preprocess(self, num_nw_vars):
        if self.dp.is_set:
            self.dp.val_SI = convert_to_SI('p', self.dp.val, self.inl[0].p.unit)

        self.o2 = "O2"
        self.h2 = "H2"
        self.h2o = "H2O"

        self.e0 = self.calc_e0()

        T_ref = 298.15
        p_ref = 1e5

        # equations to set a reference point for each h2o, h2 and o2
        self.h_refh2o = h_mix_pT(p_ref, T_ref, self.inl[1].fluid_data, self.inl[1].mixing_rule)
        self.h_refo2 = h_mix_pT(p_ref, T_ref, self.outl[1].fluid_data, self.outl[1].mixing_rule)
        self.h_refh2 = h_mix_pT(p_ref, T_ref, self.outl[2].fluid_data, self.outl[2].mixing_rule)


        super()._preprocess(num_nw_vars)

    def calc_e0(self):
        r"""
        Calculate the minimum specific energy required for electrolysis.

        Returns
        -------
        float
            Minimum specific energy.

            .. math::
                e0 = -\frac{\sum_i {\Delta H_f^0}_i -
                \sum_j {\Delta H_f^0}_j }
                {M_{H_2}}\\
                \forall i \in \text{reation products},\\
                \forall j \in \text{reation educts},\\
                \Delta H_f^0: \text{molar formation enthalpy}
        """
        hf = {}
        hf['H2O'] = -286000
        hf['H2'] = 0
        hf['O2'] = 0
        M = self.outl[2].fluid.wrapper[self.h2]._molar_mass
        e0 = -(2 * hf['H2O'] - 2 * hf['H2'] - hf['O2']) / (2 * M)

        return e0

    def gas_temperature_func(self):
        r"""
        Equation for temperature equality of product gases.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = T_\mathrm{out,2} - T_\mathrm{out,3}
        """
        return self.outl[1].calc_T() -  self.outl[2].calc_T()

    def gas_temperature_dependents(self):
        return [
            self.outl[1].p,
            self.outl[1].h,
            self.outl[2].p,
            self.outl[2].h,
        ]

    def eta_func(self):
        r"""
        Equation for electrolysis efficiency.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = P \cdot \eta - \dot{m}_{H_2,out,3} \cdot e_0
        """
        return self.P.val * self.eta.val - self.outl[2].m.val_SI * self.e0

    def eta_dependents(self):
        return [self.outl[2].m, self.P]

    def heat_func(self):
        r"""
        Equation for heat output.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \dot{Q}-\dot{m}_{in,1}\cdot \left(h_{in,1}-h_{out,1}\right)
        """
        return self.Q.val + self.inl[0].m.val_SI * (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI
        )

    def heat_dependents(self):
        return [self.inl[0].m, self.inl[0].h, self.outl[0].h]

    def specific_energy_func(self):
        r"""
        Equation for specific energy consumption.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = P - \dot{m}_{H_2,out3} \cdot e
        """
        return self.P.val - self.outl[2].m.val_SI * self.e.val

    def specific_energy_dependents(self):
        return [self.outl[2].m, self.P, self.e]

    def energy_balance_func(self):
        r"""
        Calculate the residual in energy balance.

        Returns
        -------
        residual : float
            Residual value of energy balance equation.

            .. math::

                \begin{split}
                0=&P + \dot{m}_\mathrm{in,2}\cdot\left(h_\mathrm{in,2}-
                h_\mathrm{in,2,ref}\right)\\
                &-\dot{m}_\mathrm{in,1}\cdot\left( h_\mathrm{out,1} -
                h_\mathrm{in,1} \right)\\
                & -\dot{m}_\mathrm{out,2} \cdot \left( h_\mathrm{out,2} -
                h_\mathrm{out,2,ref} \right)\\
                & +\dot{m}_\mathrm{out,3} \cdot \left( h_\mathrm{out,3} -
                h_\mathrm{out,3,ref} + e_0\right)\\
                \end{split}

            - Reference temperature: 298.15 K.
            - Reference pressure: 1 bar.
        """
        return self.P.val - self.calc_P()

    def energy_balance_deriv(self, increment_filter, k, dependents=None):
        r"""
        Partial derivatives for reactor energy balance.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        # derivatives cooling water inlet
        i = self.inl[0]
        if i.m.is_var:
            self.jacobian[k, i.m.J_col] = i.h.val_SI - self.outl[0].h.val_SI
        if i.h.is_var:
            self.jacobian[k, i.h.J_col] = i.m.val_SI

        # derivatives feed water inlet
        i = self.inl[1]
        if i.m.is_var:
            self.jacobian[k, i.m.J_col] = i.h.val_SI - self.h_refh2o
        if i.h.is_var:
            self.jacobian[k, i.h.J_col] = i.m.val_SI

        # derivative cooling water outlet
        o = self.outl[0]
        if o.h.is_var:
            self.jacobian[k, o.h.J_col] = -self.inl[0].m.val_SI

        # derivatives oxygen outlet
        o = self.outl[1]
        if o.m.is_var:
            self.jacobian[k, o.m.J_col] = -(o.h.val_SI - self.h_refo2)
        if o.h.is_var:
            self.jacobian[k, o.h.J_col] = -o.m.val_SI

        # derivatives hydrogen outlet
        o = self.outl[2]
        if o.m.is_var:
            self.jacobian[k, o.m.J_col] = -(o.h.val_SI - self.h_refh2 + self.e0)
        if o.h.is_var:
            self.jacobian[k, o.h.J_col] = -o.m.val_SI

        # derivatives for variable P
        if self.P.is_var:
            self.jacobian[k, self.P.J_col] = 1

    def energy_balance_dependents(self):
        return [
            [var for c in self.inl + self.outl for var in [c.m, c.h]] + [self.P]
        ]

    def eta_char_func(self):
        r"""
        Equation for given efficiency characteristic of a water electrolyzer.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = P - \dot{m}_{H_2,out,3} \cdot \frac{e_0}{\eta_{design}\cdot
                f\left(expr \right)}
        """
        p = self.eta_char.param
        expr = self.get_char_expr(p, **self.eta_char.char_params)
        if not expr:
            msg = (
                'Please choose a valid parameter, you want to link the '
                f'efficiency to at component {self.label}.'
            )
            logger.error(msg)
            raise ValueError(msg)

        return (
            self.P.val - self.outl[2].m.val_SI * self.e0 /
            (self.eta.design * self.eta_char.char_func.evaluate(expr))
        )

    def eta_char_deriv(self, increment_filter, k, dependents=None):
        r"""
        Partial derivatives electrolysis efficiency characteristic.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        o = self.outl[2]
        if o.m.is_var:
            self._partial_derivative(o.m, k, self.eta_char_func, increment_filter)

        if self.P.is_var:
            self.jacobian[k, self.P.J_col] = 1

    def eta_char_dependents(self):
        return [self.outl[2].m, self.P]

    def reactor_mass_flow_func(self):
        r"""
        Equations for mass conservation.

        Returns
        -------
        residual : list
            Residual values of equation.

            .. math::

                O_2 = \frac{M_{O_2}}{M_{O_2} + 2 \cdot M_{H_2}}\\
                0 =\dot{m}_\mathrm{in,1}-\dot{m}_\mathrm{out,1}\\
                0=O_2\cdot\dot{m}_\mathrm{H_{2}O,in,2}-
                \dot{m}_\mathrm{O_2,out,2}\\
                0 = \left(1 - O_2\right) \cdot \dot{m}_\mathrm{H_{2}O,in,2} -
                \dot{m}_\mathrm{H_2,out,3}
        """
        # calculate the ratio of o2 in water
        M_o2 = self.outl[1].fluid.wrapper[self.o2]._molar_mass
        M_h2 = self.outl[2].fluid.wrapper[self.h2]._molar_mass

        o2 = M_o2 / (M_o2 + 2 * M_h2)
        # equations for mass flow balance electrolyzer
        residual = []
        residual += [o2 * self.inl[1].m.val_SI - self.outl[1].m.val_SI]
        residual += [(1 - o2) * self.inl[1].m.val_SI - self.outl[2].m.val_SI]
        return residual

    def reactor_mass_flow_deriv(self, increment_filter, k, dependents=None):
        r"""
        Calculate the partial derivatives for all mass flow balance equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the mass flow equations.
        """
        M_o2 = self.outl[1].fluid.wrapper[self.o2]._molar_mass
        M_h2 = self.outl[2].fluid.wrapper[self.h2]._molar_mass

        o2 = M_o2 / (M_o2 + 2 * M_h2)

        # oxygen output to water input
        if self.inl[1].m.is_var:
            self.jacobian[k, self.inl[1].m.J_col] = o2
        if self.outl[1].m.is_var:
            self.jacobian[k, self.outl[1].m.J_col] = -1
        k += 1

        # hydrogen output to water input
        if self.inl[1].m.is_var:
            self.jacobian[k, self.inl[1].m.J_col] = (1 - o2)
        if self.outl[2].m.is_var:
            self.jacobian[k, self.outl[2].m.J_col] = -1

    def reactor_mass_flow_dependents(self):
        return [
            [self.inl[1].m, self.outl[1].m],
            [self.inl[1].m, self.outl[2].m]
        ]

    def cooling_mass_flow_func(self):
        return self.inl[0].m.val_SI - self.outl[0].m.val_SI

    def cooling_mass_flow_structure_matrix(self, k):
        self._structure_matrix[k, self.inl[0].m.sm_col] = 1
        self._structure_matrix[k, self.outl[0].m.sm_col] = -1

    def cooling_fluid_func(self):
        return 0

    def cooling_fluid_structure_matrix(self, k):
        self._structure_matrix[k, self.inl[0].fluid.sm_col] = 1
        self._structure_matrix[k, self.outl[0].fluid.sm_col] = -1

    def reactor_pressure_func(self):
        r"""
        Equations for reactor pressure balance.

        Returns
        -------
        residual : list
            Residual values of equation.

            .. math::

                0 = p_\mathrm{in,2} - p_\mathrm{out,2}\\
                0 = p_\mathrm{in,3} - p_\mathrm{out,2}
        """
        return [
            self.inl[1].p.val_SI - self.outl[1].p.val_SI,
            self.inl[1].p.val_SI - self.outl[2].p.val_SI
        ]

    def reactor_pressure_structure_matrix(self, k):
        r"""
        Equations for reactor pressure balance.

        Returns
        -------
        residual : list
            Residual values of equations.

            .. math::

                0 = p_\mathrm{in,2} - p_\mathrm{out,2}\\
                0 = p_\mathrm{in,3} - p_\mathrm{out,2}
        """
        self._structure_matrix[k, self.inl[1].p.sm_col] = 1
        self._structure_matrix[k, self.outl[1].p.sm_col] = -1

        self._structure_matrix[k + 1, self.inl[1].p.sm_col] = 1
        self._structure_matrix[k + 1, self.outl[2].p.sm_col] = -1

    def calc_P(self):
        r"""
        Calculate water electrolyzer power input.

        Returns
        -------
        P : float
            Value of power input.

            .. math::

                \begin{split}
                P = & -\dot{m}_{in,2} \cdot \left( h_{in,2} - h_{in,2,ref}
                \right)\\
                & + \dot{m}_{in,1} \cdot \left( h_{out,1} - h_{in,1} \right)\\
                & + \dot{m}_{out,2} \cdot \left( h_{out,2} - h_{out,2,ref}
                \right)\\
                & - \dot{m}_{out,3} \cdot \left( h_{out,3} - h_{out,3,ref}
                + e_0\right)\\
                \end{split}

        Note
        ----
        The temperature for the reference state is set to 25 °C, thus
        the feed water must be liquid as proposed in the calculation of
        the minimum specific energy consumption for electrolysis:
        :py:meth:`tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.calc_e0`.
        The part of the equation regarding the cooling water is implemented
        with negative sign as the energy for cooling is extracted from the
        reactor.

        - Reference temperature: 298.15 K.
        - Reference pressure: 1 bar.
        """
        val = (
            -self.inl[1].m.val_SI * (self.inl[1].h.val_SI - self.h_refh2o)
            + self.inl[0].m.val_SI * (self.outl[0].h.val_SI - self.inl[0].h.val_SI)
            + self.outl[1].m.val_SI * (self.outl[1].h.val_SI - self.h_refo2)
            + self.outl[2].m.val_SI * (self.outl[2].h.val_SI - self.h_refh2 + self.e0)
        )
        return val

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

                \dot{E} = \begin{cases}
                P & \text{key = 'P'}\\
                - \dot{m}_{in,1} \cdot \left(h_{out,1} - h_{in,1} \right) &
                \text{key = 'Q'}\\
                \end{cases}
        """
        ######################################################################
        # equations for power on bus
        if bus['param'] == 'P':
            val = self.calc_P()

        ######################################################################
        # equations for heat on bus

        elif bus['param'] == 'Q':
            val = -self.inl[0].m.val_SI * (
                self.outl[0].h.val_SI - self.inl[0].h.val_SI
            )

        ######################################################################
        # missing/invalid bus parameter

        else:
            msg = (
                f'The parameter {bus["param"]} is not a valid parameter for a '
                f'component of type {self.__class__.__name__}. Please specify '
                f'a bus parameter (P/Q) for component {self.label}.'
            )
            logger.error(msg)
            raise ValueError(msg)

        return val

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
        b = bus.comps.loc[self]

        ######################################################################
        # derivatives for power on bus
        if b['param'] == 'P':
            for c in self.inl + self.outl:
                if c.m.is_var and c != self.outl[0]:
                    if c.m.J_col not in bus.jacobian:
                        bus.jacobian[c.m.J_col] = 0
                    bus.jacobian[c.m.J_col] -= _numeric_deriv(c.m._reference_container, f, bus=bus)

                if c.h.is_var:
                    if c.h.J_col not in bus.jacobian:
                        bus.jacobian[c.h.J_col] = 0
                    bus.jacobian[c.h.J_col] -= _numeric_deriv(c.h._reference_container, f, bus=bus)

            # variable power
            if self.P.is_var:
                if self.P.J_col not in bus.jacobian:
                    bus.jacobian[self.P.J_col] = 0
                bus.jacobian[self.P.J_col] -= _numeric_deriv(self.P._reference_container, f, bus=bus)

        ######################################################################
        # derivatives for heat on bus
        elif b['param'] == 'Q':

            i = self.inl[0]
            o = self.outl[0]
            if i.m.is_var:
                if i.m.J_col not in bus.jacobian:
                    bus.jacobian[i.m.J_col] = 0
                bus.jacobian[i.m.J_col] -= _numeric_deriv(i.m._reference_container, f, bus=bus)

            if i.h.is_var:
                if i.h.J_col not in bus.jacobian:
                    bus.jacobian[i.h.J_col] = 0
                bus.jacobian[i.h.J_col] -= _numeric_deriv(i.h._reference_container, f, bus=bus)

            if o.h.is_var:
                if o.h.J_col not in bus.jacobian:
                    bus.jacobian[o.h.J_col] = 0
                bus.jacobian[o.h.J_col] -= _numeric_deriv(o.h._reference_container, f, bus=bus)

    def start_fluid_wrapper_branch(self):
        branches = {}
        for outconn in self.outl[1:]:
            branch = {
                "connections": [outconn],
                "components": [self]
            }
            outconn.target.propagate_wrapper_to_target(branch)
            branches.update({outconn.label: branch})

        return branches

    def propagate_wrapper_to_target(self, branch):
        inconn = branch["connections"][-1]
        if inconn == self.inl[0]:
            conn_idx = self.inl.index(inconn)
            outconn = self.outl[conn_idx]

            branch["connections"] += [outconn]
            branch["components"] += [self]

            outconn.target.propagate_wrapper_to_target(branch)
        else:
            branch["components"] += [self]
            return

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
            return 5e5
        elif key == 'h':
            temp = 50 + 273.15
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
            return 5e5
        elif key == 'h':
            temp = 20 + 273.15
            return h_mix_pT(c.p.val_SI, temp, c.fluid_data, c.mixing_rule)

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        self.Q.val = -self.inl[0].m.val_SI * (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI
        )
        self.pr.val = self.outl[0].p.val_SI / self.inl[0].p.val_SI
        self.dp.val_SI = self.inl[0].p.val_SI - self.outl[0].p.val_SI
        self.dp.val = self.inl[0].p.val - self.outl[0].p.val
        self.e.val = self.P.val / self.outl[2].m.val_SI
        self.eta.val = self.e0 / self.e.val

        i = self.inl[0]
        o = self.outl[0]
        self.zeta.val = self.calc_zeta(i, o)

    def exergy_balance(self, T0):
        self.E_P = (
            self.outl[1].Ex_chemical + self.outl[2].Ex_chemical
            - self.inl[1].Ex_chemical + self.outl[0].Ex_physical
            + self.inl[0].Ex_physical
        )
        self.E_F = self.P.val

        self.E_D = self.E_F - self.E_P
        self.epsilon = self._calc_epsilon()
        self.E_bus = self.P.val
