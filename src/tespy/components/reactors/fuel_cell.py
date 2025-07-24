# -*- coding: utf-8

"""Module of class FuelCell.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/reactors/fuel_cell.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import Component
from tespy.components.component import component_registry
from tespy.tools import logger
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.helpers import convert_to_SI


@component_registry
class FuelCell(Component):
    r"""
    The fuel cell produces power by oxidation of hydrogen.

    **Mandatory Equations**

    - :py:meth:`tespy.components.reactors.fuel_cell.FuelCell.cooling_fluid_func`
    - :py:meth:`tespy.components.reactors.fuel_cell.FuelCell.reactor_mass_flow_func`
    - :py:meth:`tespy.components.reactors.fuel_cell.FuelCell.cooling_mass_flow_func`
    - :py:meth:`tespy.components.reactors.fuel_cell.FuelCell.reactor_pressure_func`
    - :py:meth:`tespy.components.reactors.fuel_cell.FuelCell.energy_balance_func`

    **Optional Equations**

    - cooling loop:

      - :py:meth:`tespy.components.component.Component.zeta_func`
      - :py:meth:`tespy.components.component.Component.pr_func`

    - :py:meth:`tespy.components.reactors.fuel_cell.FuelCell.eta_func`
    - :py:meth:`tespy.components.reactors.fuel_cell.FuelCell.heat_func`
    - :py:meth:`tespy.components.reactors.fuel_cell.FuelCell.specific_energy_func`

    Inlets/Outlets

    - in1 (cooling inlet), in2 (oxygen inlet), in3 (hydrogen inlet)
    - out1 (cooling outlet), out2 (water outlet)

    Image

    .. image:: _images/FuelCell.svg
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
        Electrolysis efficiency, :math:`\eta/1`.

    pr : float, dict, :code:`"var"`
        Cooling loop pressure ratio, :math:`pr/1`.

    zeta : float, dict, :code:`"var"`
        Geometry independent friction coefficient for cooling loop pressure
        drop, :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    Note
    ----
    Other than usual components, the fuel cell has the fluid composition
    built into its equations for the feed hydrogen and oxygen inlets as well
    as the water outlet. Thus, the user must not specify the fluid composition
    at these connections!

    Example
    -------
    The example shows a simple adaptation of the fuel cell. It works with water
    as cooling fluid.

    >>> from tespy.components import (Sink, Source, FuelCell)
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> from tespy.tools import ComponentCharacteristics as dc_cc
    >>> nw = Network(T_unit='C', p_unit='bar', v_unit='l / s', iterinfo=False)
    >>> fc = FuelCell('fuel cell')
    >>> oxygen_source = Source('oxygen_source')
    >>> hydrogen_source = Source('hydrogen_source')
    >>> cw_source = Source('cw_source')
    >>> cw_sink = Sink('cw_sink')
    >>> water_sink = Sink('water_sink')
    >>> cw_in = Connection(cw_source, 'out1', fc, 'in1')
    >>> cw_out = Connection(fc, 'out1', cw_sink, 'in1')
    >>> oxygen_in = Connection(oxygen_source, 'out1', fc, 'in2')
    >>> hydrogen_in = Connection(hydrogen_source, 'out1', fc, 'in3')
    >>> water_out = Connection(fc, 'out2', water_sink, 'in1')
    >>> nw.add_conns(cw_in, cw_out, oxygen_in, hydrogen_in, water_out)

    The fuel cell produces 200kW of electrical power with an efficiency of 0.45.
    The thermodynamic parameters of the input oxygen and hydrogen are given,
    the mass flow rates are calculated out of the given power output. The
    temperature of the water at the outlet should be 50 °C. The cooling fluid is
    pure water and is heated up from 25 °C to 40 °C.

    >>> fc.set_attr(eta=0.45, P=-200e03, pr=0.9)
    >>> cw_in.set_attr(T=25, p=1, fluid={'H2O': 1})
    >>> cw_out.set_attr(T=40)
    >>> oxygen_in.set_attr(T=25, p=1)
    >>> hydrogen_in.set_attr(T=25)
    >>> water_out.set_attr(T=50)
    >>> nw.solve('design')
    >>> round(cw_in.m.val, 1)
    10.2
    >>> Q = fc.Q.val / 1e3
    >>> round(Q, 0)
    -642.0
    >>> round(fc.eta.val, 2)
    0.45
    """

    def get_parameters(self):
        return {
            'P': dc_cp(max_val=0),
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
            'e': dc_cp(
                max_val=0, num_eq_sets=1,
                func=self.specific_energy_func,
                dependents=self.specific_energy_dependents
            )
        }

    def get_mandatory_constraints(self):
        return {
            'mass_flow_constraints': dc_cmc(**{
                'func': self.reactor_mass_flow_func,
                'deriv': self.reactor_mass_flow_deriv,
                'dependents': self.reactor_mass_flow_dependents,
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
                'dependents': self.energy_balance_dependents,
                'constant_deriv': False,
                'num_eq_sets': 1
            }),
            'reactor_pressure_constraints': dc_cmc(**{
                'func': self.reactor_pressure_func,
                'structure_matrix': self.reactor_pressure_structure_matrix,
                'constant_deriv': True,
                'num_eq_sets': 2
            })
        }

    @staticmethod
    def get_bypass_constraints():
        return {}

    @staticmethod
    def inlets():
        return ['in1', 'in2', 'in3']

    @staticmethod
    def outlets():
        return ['out1', 'out2']

    def _add_missing_fluids(self, connections):
        if self.inl[1] in connections:
            return ["O2"]
        elif self.inl[2] in connections:
            return ["H2"]
        elif self.outl[1] in connections:
            return ["H2O"]
        else:
            return super()._add_missing_fluids(connections)

    def get_variables(self):
        if not self.P.is_set:
            self.set_attr(P='var')
            msg = (
                f'The power output of fuel cells ({self.label}) must either '
                'be a set value or part of the system variables. Since it has '
                'not been set to a fixed value it will be added to the '
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

        # derivatives determined from calc_P function
        T_ref = 298.15
        p_ref = 1e5
        self.h_refh2o = h_mix_pT(p_ref, T_ref, self.outl[1].fluid_data, self.outl[1].mixing_rule)
        self.h_refo2 = h_mix_pT(p_ref, T_ref, self.inl[1].fluid_data, self.inl[1].mixing_rule)
        self.h_refh2 = h_mix_pT(p_ref, T_ref, self.inl[2].fluid_data, self.inl[2].mixing_rule)

        super()._preprocess(num_nw_vars)

    def calc_e0(self):
        r"""
        Calculate the specific energy output of the fuel cell.

        Returns
        -------
        float
            Specific energy.

            .. math::

                e0 = \frac{\sum_i {\Delta H_f^0}_i -
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
        M = self.inl[2].fluid.wrapper[self.h2]._molar_mass
        e0 = (2 * hf['H2O'] - 2 * hf['H2'] - hf['O2']) / (2 * M)

        return e0

    def eta_func(self):
        r"""
        Equation for efficiency.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = P - \eta \cdot \dot{m}_{H_2,in} \cdot e_0
        """
        return self.P.val - self.eta.val * self.inl[2].m.val_SI * self.e0

    def eta_dependents(self):
        return [self.inl[2].m, self.P]

    def heat_func(self):
        r"""
        Equation for heat output.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \dot{Q}-\dot{m}_{in,1}\cdot \left(h_{out,1}-h_{in,1}\right)
        """
        return self.Q.val + self.inl[0].m.val_SI * (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI
        )

    def heat_dependents(self):
        return [self.inl[0].m, self.inl[0].h, self.outl[0].h]

    def specific_energy_func(self):
        r"""
        Equation for specific energy output.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = P - \dot{m}_{H_2,in} \cdot e
        """
        return self.P.val - self.inl[2].m.val_SI * self.e.val

    def specific_energy_dependents(self):
        return [self.inl[2].m, self.P, self.e]

    def energy_balance_func(self):
        r"""
        Calculate the residual in energy balance.

        Returns
        -------
        residual : float
            Residual value of energy balance equation.

            .. math::

                \begin{split}
                0=&P + \dot{m}_\mathrm{out,2}\cdot\left(h_\mathrm{out,2}-
                h_\mathrm{out,2,ref}\right)\\
                &+\dot{m}_\mathrm{in,1}\cdot\left( h_\mathrm{out,1} -
                h_\mathrm{in,1} \right)\\
                & -\dot{m}_\mathrm{in,2} \cdot \left( h_\mathrm{in,2} -
                h_\mathrm{in,2,ref} \right)\\
                & -\dot{m}_\mathrm{in,3} \cdot \left( h_\mathrm{in,3} -
                h_\mathrm{in,3,ref} - e_0\right)\\
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
            self.jacobian[k, i.m.J_col] = self.outl[0].h.val_SI - i.h.val_SI
        if i.h.is_var:
            self.jacobian[k, i.h.J_col] = -i.m.val_SI
        # derivative cooling water outlet
        o = self.outl[0]
        if o.h.is_var:
            self.jacobian[k, o.h.J_col] = self.inl[0].m.val_SI

        # derivatives water outlet
        o = self.outl[1]
        if o.m.is_var:
            self.jacobian[k, o.m.J_col] = o.h.val_SI - self.h_refh2o
        if o.h.is_var:
            self.jacobian[k, o.h.J_col] = o.m.val_SI

        # derivatives oxygen inlet
        i = self.inl[1]
        if i.m.is_var:
            self.jacobian[k, i.m.J_col] = -(i.h.val_SI - self.h_refo2)
        if i.h.is_var:
            self.jacobian[k, i.h.J_col] = -i.m.val_SI

        # derivatives hydrogen inlet
        i = self.inl[2]
        if i.m.is_var:
            self.jacobian[k, i.m.J_col] = -(i.h.val_SI - self.h_refh2 - self.e0)
        if i.h.is_var:
            self.jacobian[k, i.h.J_col] = -i.m.val_SI

        # derivatives for variable P
        if self.P.is_var:
            self.jacobian[k, self.p.J_col] = 1

    def energy_balance_dependents(self):
        return [
            [var for c in self.inl + self.outl for var in [c.m, c.h]] + [self.P]
        ]

    def reactor_mass_flow_func(self):
        r"""
        Equations for mass conservation.

        Returns
        -------
        residual : list
            Residual values of equation.

            .. math::

                O_2 = \frac{M_{O_2}}{M_{O_2} + 2 \cdot M_{H_2}}\\
                0=O_2\cdot\dot{m}_\mathrm{H_{2}O,out,1}-
                \dot{m}_\mathrm{O_2,in,2}\\
                0 = \left(1 - O_2\right) \cdot \dot{m}_\mathrm{H_{2}O,out,1} -
                \dot{m}_\mathrm{H_2,in,1}
        """
        # calculate the ratio of o2 in water
        M_o2 = self.inl[1].fluid.wrapper[self.o2]._molar_mass
        M_h2 = self.inl[2].fluid.wrapper[self.h2]._molar_mass
        o2 = M_o2 / (M_o2 + 2 * M_h2)
        # equations for mass flow balance of the fuel cell
        residual = []
        residual += [o2 * self.outl[1].m.val_SI - self.inl[1].m.val_SI]
        residual += [(1 - o2) * self.outl[1].m.val_SI - self.inl[2].m.val_SI]
        return residual

    def reactor_mass_flow_deriv(self, increment_filter, k, dependents=None):
        r"""
        Calculate the partial derivatives for all mass flow balance equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the mass flow equations.
        """
        M_o2 = self.inl[1].fluid.wrapper[self.o2]._molar_mass
        M_h2 = self.inl[2].fluid.wrapper[self.h2]._molar_mass
        o2 = M_o2 / (M_o2 + 2 * M_h2)
        # oxygen input to water output
        if self.inl[1].m.is_var:
            self.jacobian[k, self.inl[1].m.J_col] = -1
        if self.outl[1].m.is_var:
            self.jacobian[k, self.outl[1].m.J_col] = o2
        k += 1

        # hydrogen input to water output
        if self.inl[2].m.is_var:
            self.jacobian[k, self.inl[2].m.J_col] = -1
        if self.outl[1].m.is_var:
            self.jacobian[k, self.outl[1].m.J_col] = (1 - o2)

    def reactor_mass_flow_dependents(self):
        return [
            [self.inl[1].m, self.outl[1].m],
            [self.inl[2].m, self.outl[1].m]
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
            self.outl[1].p.val_SI - self.inl[1].p.val_SI,
            self.outl[1].p.val_SI - self.inl[2].p.val_SI
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

        self._structure_matrix[k + 1, self.inl[2].p.sm_col] = 1
        self._structure_matrix[k + 1, self.outl[1].p.sm_col] = -1

    def calc_P(self):
        r"""
        Calculate fuel cell power output.

        Returns
        -------
        P : float
            Value of power output.

            .. math::

                \begin{split}
                P = & +\dot{m}_{in,2} \cdot \left( h_{in,2} - h_{in,2,ref}
                \right)\\
                & + \dot{m}_{in,3} \cdot \left( h_{in,3} - h_{in,3,ref} - e_0
                \right)\\
                & - \dot{m}_{in,1} \cdot \left( h_{out,1} - h_{in,1} \right)\\
                & - \dot{m}_{out,2} \cdot \left( h_{out,2} - h_{out,2,ref}
                \right)\\
                \end{split}

        Note
        ----
        The temperature for the reference state is set to 25 °C, thus
        the produced water must be liquid as proposed in the calculation of
        the minimum specific energy for oxidation:
        :py:meth:`tespy.components.reactors.fuel_cell.FuelCell.calc_e0`.
        The part of the equation regarding the cooling water is implemented
        with negative sign as the energy for cooling is extracted from the
        reactor.
        - Reference temperature: 298.15 K.
        - Reference pressure: 1 bar.
        """
        val = (
            self.inl[2].m.val_SI * (self.inl[2].h.val_SI - self.h_refh2 - self.e0)
            + self.inl[1].m.val_SI * (self.inl[1].h.val_SI - self.h_refo2)
            - self.inl[0].m.val_SI * (self.outl[0].h.val_SI - self.inl[0].h.val_SI)
            - self.outl[1].m.val_SI * (self.outl[1].h.val_SI - self.h_refh2o)
        )
        return val

    def start_fluid_wrapper_branch(self):
        outconn = self.outl[1]
        branch = {
            "connections": [outconn],
            "components": [self]
        }
        outconn.target.propagate_wrapper_to_target(branch)
        return {outconn.label: branch}

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

    def initialise_target(self, c, key):
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

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        self.Q.val = -self.inl[0].m.val_SI * (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI
        )
        self.pr.val = self.outl[0].p.val_SI / self.inl[0].p.val_SI
        self.dp.val_SI = self.inl[0].p.val_SI - self.outl[0].p.val_SI
        self.dp.val = self.inl[0].p.val - self.outl[0].p.val
        self.e.val = self.P.val / self.inl[2].m.val_SI
        self.eta.val = self.e.val / self.e0

        i = self.inl[0]
        o = self.outl[0]
        self.zeta.val = self.calc_zeta(i, o)
