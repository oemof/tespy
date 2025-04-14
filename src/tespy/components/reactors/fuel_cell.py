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
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.document_models import generate_latex_eq
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.helpers import convert_to_SI


@component_registry
class FuelCell(Component):
    r"""
    The fuel cell produces power by oxidation of hydrogen.

    **Mandatory Equations**

    - :py:meth:`tespy.components.reactors.fuel_cell.FuelCell.fluid_func`
    - :py:meth:`tespy.components.reactors.fuel_cell.FuelCell.mass_flow_func`
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
    >>> import shutil
    >>> nw = Network(T_unit='C', p_unit='bar', v_unit='l / s', iterinfo=False)
    >>> fc = FuelCell('fuel cell')
    >>> fc.component()
    'fuel cell'
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
    @staticmethod
    def component():
        return 'fuel cell'

    def get_parameters(self):
        return {
            'P': dc_cp(max_val=0),
            'Q': dc_cp(
                max_val=0, num_eq=1,
                deriv=self.heat_deriv, func=self.heat_func,
                latex=self.heat_func_doc),
            'pr': dc_cp(
                max_val=1, num_eq=1,
                deriv=self.pr_deriv, func=self.pr_func,
                func_params={'pr': 'pr'}, latex=self.pr_func_doc),
            'dp': dc_cp(
                min_val=0, deriv=self.dp_deriv,
                func=self.dp_func,
                num_eq=1, func_params={"inconn": 0, "outconn": 0, "dp": "dp"}
            ),
            'zeta': dc_cp(
                min_val=0, num_eq=1,
                deriv=self.zeta_deriv, func=self.zeta_func,
                func_params={'zeta': 'zeta'}, latex=self.zeta_func_doc),
            'eta': dc_cp(
                min_val=0, max_val=1, num_eq=1, latex=self.eta_func_doc,
                deriv=self.eta_deriv, func=self.eta_func),
            'e': dc_cp(
                max_val=0, num_eq=1,
                deriv=self.specific_energy_deriv,
                func=self.specific_energy_func,
                latex=self.specific_energy_func_doc)
        }

    def get_mandatory_constraints(self):
        num_mass_eq = (
            (self.inl[1].m.is_var or self.outl[1].m.is_var)
            + (self.inl[1].m.is_var or self.outl[2].m.is_var)
        )
        return {
            'mass_flow_constraints': {
                'func': self.mass_flow_func, 'deriv': self.mass_flow_deriv,
                'constant_deriv': True, 'latex': self.mass_flow_func_doc,
                'num_eq': num_mass_eq},
            'energy_balance_constraints': {
                'func': self.energy_balance_func,
                'deriv': self.energy_balance_deriv,
                'constant_deriv': False, 'latex': self.energy_balance_func_doc,
                'num_eq': 1},
            'reactor_pressure_constraints': {
                'func': self.reactor_pressure_func,
                'deriv': self.reactor_pressure_deriv,
                'constant_deriv': True,
                'latex': self.reactor_pressure_func_doc,
                'num_eq': 2},
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

    def preprocess(self, num_nw_vars):

        if not self.P.is_set:
            self.set_attr(P='var')
            msg = ('The power output of a fuel cell must be set! '
                   'We are adding the power output of component ' +
                   self.label + ' as custom variable of the system.')
            logger.info(msg)

        self.o2 = "O2"
        self.h2 = "H2"
        self.h2o = "H2O"
        self.e0 = self.calc_e0()

        super().preprocess(num_nw_vars)

        if self.dp.is_set:
            self.dp.val_SI = convert_to_SI('p', self.dp.val, self.inl[0].p.unit)

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

    def eta_func_doc(self, label):
        r"""
        Equation for efficiency.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        latex = r'0 = P - \eta \cdot \dot{m}_\mathrm{H_2,in,3} \cdot e_0'
        return generate_latex_eq(self, latex, label)

    def eta_deriv(self, increment_filter, k):
        r"""
        Partial derivatives for efficiency function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        if self.inl[2].m.is_var:
            self.jacobian[k, self.inl[2].m.J_col] = -self.eta.val * self.e0
        # derivatives for variable P
        if self.P.is_var:
            self.jacobian[k, self.P.J_col] = 1

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
            self.outl[0].h.val_SI - self.inl[0].h.val_SI)

    def heat_func_doc(self, label):
        r"""
        Equation for heat output.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        latex = (
            r'0=\dot{Q}+\dot{m}_\mathrm{in,1}\cdot\left(h_\mathrm{out,1}-'
            r'h_\mathrm{in,1}\right)')
        return generate_latex_eq(self, latex, label)

    def heat_deriv(self, increment_filter, k):
        r"""
        Partial derivatives for heat output function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        i = self.inl[0]
        o = self.outl[0]
        if i.m.is_var:
            self.jacobian[k, i.m.J_col] = o.h.val_SI - i.h.val_SI
        if i.h.is_var:
            self.jacobian[k, i.h.J_col] = -i.m.val_SI
        if o.h.is_var:
            self.jacobian[k, o.h.J_col] = i.m.val_SI

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

    def specific_energy_func_doc(self, label):
        r"""
        Equation for specific energy output.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        latex = r'0=P - \dot{m}_\mathrm{H_2,in} \cdot e'
        return generate_latex_eq(self, latex, label)

    def specific_energy_deriv(self, increment_filter, k):
        r"""
        Partial derivatives for specific energy function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        if self.inl[2].m.is_var:
            self.jacobian[k, self.inl[2].m.J_col] = -self.e.val
        # derivatives for variable P
        if self.P.is_var:
            self.jacobian[k, self.P.J_col] = 1
        # derivatives for variable e
        if self.e.is_var:
            self.jacobian[k, self.e.J_col] = -self.inl[2].m.val_SI

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

    def energy_balance_func_doc(self, label):
        r"""
        Calculate the residual in energy balance.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        latex = (
            r'\begin{split}' + '\n'
            r'0=&P + \dot{m}_\mathrm{out,2}\cdot\left(h_\mathrm{out,2}-'
            r'h_\mathrm{out,2,ref}\right)\\' + '\n'
            r'&+\dot{m}_\mathrm{in,1}\cdot\left( h_\mathrm{out,1} -'
            r'h_\mathrm{in,1} \right)\\' + '\n'
            r'& - \dot{m}_\mathrm{in,2} \cdot \left( h_\mathrm{in,2} -'
            r'h_\mathrm{in,2,ref} \right)\\' + '\n'
            r'& - \dot{m}_\mathrm{in,3} \cdot \left( h_\mathrm{in,3} -'
            r'h_\mathrm{in,3,ref} - e_0\right)\\' + '\n'
            r'&p_\mathrm{ref}=\unit[1]{bar},'
            r'\;T_\mathrm{ref}=\unit[25]{^\circ C}\\' + '\n'
            r'\end{split}'
        )
        return generate_latex_eq(self, latex, label)

    def energy_balance_deriv(self, increment_filter, k):
        r"""
        Partial derivatives for reactor energy balance.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        # derivatives determined from calc_P function
        T_ref = 298.15
        p_ref = 1e5
        h_refh2o = h_mix_pT(p_ref, T_ref, self.outl[1].fluid_data, self.outl[1].mixing_rule)
        h_refo2 = h_mix_pT(p_ref, T_ref, self.inl[1].fluid_data, self.inl[1].mixing_rule)
        h_refh2 = h_mix_pT(p_ref, T_ref, self.inl[2].fluid_data, self.inl[2].mixing_rule)

        # derivatives cooling water inlet
        i = self.inl[0]
        if i.m.is_var:
            self.jacobian[k, i.m.J_col] = self.outl[0].h.val_SI - i.h.val_SI
        if i.h.is_var:
            self.jacobian[k, i.h.J_col] = -i.m.val_SI

        # derivatives water outlet
        o = self.outl[1]
        if o.m.is_var:
            self.jacobian[k, o.m.J_col] = o.h.val_SI - h_refh2o
        if o.h.is_var:
            self.jacobian[k, o.h.J_col] = o.m.val_SI

        # derivative cooling water outlet
        o = self.outl[0]
        if o.h.is_var:
            self.jacobian[k, o.h.J_col] = self.inl[0].m.val_SI

        # derivatives oxygen inlet
        i = self.inl[1]
        if i.m.is_var:
            self.jacobian[k, i.m.J_col] = -(i.h.val_SI - h_refo2)
        if i.h.is_var:
            self.jacobian[k, i.h.J_col] = -i.m.val_SI

        # derivatives hydrogen inlet
        i = self.inl[2]
        if i.m.is_var:
            self.jacobian[k, i.m.J_col] = -(i.h.val_SI - h_refh2 - self.e0)
        if i.h.is_var:
            self.jacobian[k, i.h.J_col] = -i.m.val_SI

        # derivatives for variable P
        if self.P.is_var:
            self.jacobian[k, self.P.J_col] = 1

    def mass_flow_func(self):
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
        if self.inl[1].m.is_var or self.outl[1].m.is_var:
            residual += [o2 * self.outl[1].m.val_SI - self.inl[1].m.val_SI]
        if self.inl[2].m.is_var or self.outl[1].m.is_var:
            residual += [(1 - o2) * self.outl[1].m.val_SI - self.inl[2].m.val_SI]
        return residual

    def mass_flow_func_doc(self, label):
        r"""
        Equations for mass conservation.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        latex = (
            r'\begin{split}' + '\n'
            r'O_2 = &\frac{M_{O_2}}{M_{O_2} + 2 \cdot M_{H_2}}\\' + '\n'
            r'0=&O_2\cdot\dot{m}_\mathrm{H_{2}O,out,1}-'
            r'\dot{m}_\mathrm{O_2,in,2}\\' + '\n'
            r'0 =&\left(1 - O_2\right) \cdot \dot{m}_\mathrm{H_{2}O,out,2}-'
            r'\dot{m}_\mathrm{H_2,in,3}\\' + '\n'
            r'\end{split}'
        )
        return generate_latex_eq(self, latex, label)

    def mass_flow_deriv(self, increment_filter, k):
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
        # number of equations may vary here
        if self.inl[1].m.is_var or self.outl[1].m.is_var:
            if self.inl[1].m.is_var:
                self.jacobian[k, self.inl[1].m.J_col] = -1
            if self.outl[1].m.is_var:
                self.jacobian[k, self.outl[1].m.J_col] = o2
            k += 1

        # derivatives for mass flow balance for hydrogen input
        if self.outl[1].m.is_var:
            self.jacobian[k, self.outl[1].m.J_col] = (1 - o2)
        if self.inl[2].m.is_var:
            self.jacobian[k, self.inl[2].m.J_col] = -1

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

    def reactor_pressure_func_doc(self, label):
        r"""
        Equations for reactor pressure balance.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        latex = (
            r'\begin{split}' + '\n'
            r'0 = & p_\mathrm{in,2} - p_\mathrm{out,2}\\' + '\n'
            r'0 = & p_\mathrm{in,3} - p_\mathrm{out,2}\\' + '\n'
            r'\end{split}')
        return generate_latex_eq(self, latex, label)

    def reactor_pressure_deriv(self, increment_filter, k):
        r"""
        Calculate the partial derivatives for combustion pressure equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the pressure equations.
        """
        o = self.outl[1]
        for i in self.inl[1:]:
            if i.p.is_var:
                self.jacobian[k, i.p.J_col] = -1
            if o.p.is_var:
                self.jacobian[k, o.p.J_col] = 1
            k += 1

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
        T_ref = 298.15
        p_ref = 1e5

        # equations to set a reference point for each h2o, h2 and o2        # derivatives determined from calc_P function
        h_refh2o = h_mix_pT(p_ref, T_ref, self.outl[1].fluid_data, self.outl[1].mixing_rule)
        h_refo2 = h_mix_pT(p_ref, T_ref, self.inl[1].fluid_data, self.inl[1].mixing_rule)
        h_refh2 = h_mix_pT(p_ref, T_ref, self.inl[2].fluid_data, self.inl[2].mixing_rule)

        val = (
            self.inl[2].m.val_SI * (self.inl[2].h.val_SI - h_refh2 - self.e0)
            + self.inl[1].m.val_SI * (self.inl[1].h.val_SI - h_refo2)
            - self.inl[0].m.val_SI * (self.outl[0].h.val_SI - self.inl[0].h.val_SI)
            - self.outl[1].m.val_SI * (self.outl[1].h.val_SI - h_refh2o)
        )
        return val



    @staticmethod
    def is_branch_source():
        return True

    def start_branch(self):
        outconn = self.outl[1]
        if "H2O" not in outconn.fluid.val:
            outconn.fluid.val["H2O"] = 1
        branch = {
            "connections": [outconn],
            "components": [self, outconn.target],
            "subbranches": {}
        }
        outconn.target.propagate_to_target(branch)
        return {outconn.label: branch}

    def start_fluid_wrapper_branch(self):
        outconn = self.outl[1]
        branch = {
            "connections": [outconn],
            "components": [self]
        }
        outconn.target.propagate_wrapper_to_target(branch)
        return {outconn.label: branch}

    def propagate_to_target(self, branch):
        inconn = branch["connections"][-1]
        if inconn == self.inl[0]:
            conn_idx = self.inl.index(inconn)
            outconn = self.outl[conn_idx]

            branch["connections"] += [outconn]
            branch["components"] += [outconn.target]

            outconn.target.propagate_to_target(branch)
        else:
            if inconn == self.inl[1] and "O2" not in inconn.fluid.val:
                inconn.fluid.val["O2"] = 1
            if inconn == self.inl[2] and "H2" not in inconn.fluid.val:
                inconn.fluid.val["H2"] = 1
            return

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

            .. math::

                val = \begin{cases}
                5  \cdot 10^5 & \text{key = 'p'}\\
                h\left(T=293.15, p=5  \cdot 10^5\right) & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 5e5
        elif key == 'h':
            T = 20 + 273.15
            return h_mix_pT(c.p.val_SI, T, c.fluid_data, c.mixing_rule)

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

            .. math::

                val = \begin{cases}
                5  \cdot 10^5 & \text{key = 'p'}\\
                h\left(T=323.15, p=5  \cdot 10^5\right) & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 5e5
        elif key == 'h':
            T = 50 + 273.15
            return h_mix_pT(c.p.val_SI, T, c.fluid_data, c.mixing_rule)

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
