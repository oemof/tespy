# -*- coding: utf-8

"""Module of class FuelCell.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/reactors/fuel_cell.py

SPDX-License-Identifier: MIT
"""


import CoolProp.CoolProp as CP
import numpy as np

from tespy.components.component import Component
from tespy.tools import logger
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.document_models import generate_latex_eq
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.global_vars import molar_masses
from tespy.tools.helpers import TESPyComponentError


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
    >>> fluid_list = ['H2O', 'O2', 'H2']
    >>> nw = Network(fluids=fluid_list, T_unit='C', p_unit='bar',
    ... v_unit='l / s', iterinfo=False)
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

    The fuel cell shall produce 200kW of electrical power and 200kW of heat
    with an efficiency of 0.45. The thermodynamic parameters of the input
    oxygen and hydrogen are given, the mass flow rates are calculated out of
    the given power output. The cooling fluid is pure water.

    >>> fc.set_attr(eta=0.45, P=-200e03, Q=-200e03, pr=0.9)
    >>> cw_in.set_attr(T=25, p=1, m=1, fluid={'H2O': 1, 'O2': 0, 'H2': 0})
    >>> oxygen_in.set_attr(T=25, p=1)
    >>> hydrogen_in.set_attr(T=25)
    >>> nw.solve('design')
    >>> P_design = fc.P.val / 1e3
    >>> round(P_design, 0)
    -200.0
    >>> round(fc.eta.val, 2)
    0.45
    """
    @staticmethod
    def component():
        return 'fuel cell'

# %% Variables

    def get_variables(self):
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

# %% Mandatory constraints

    def get_mandatory_constraints(self):
        return {
            'mass_flow_constraints': {
                'func': self.mass_flow_func, 'deriv': self.mass_flow_deriv,
                'constant_deriv': True, 'latex': self.mass_flow_func_doc,
                'num_eq': 3},
            'fluid_constraints': {
                'func': self.fluid_func, 'deriv': self.fluid_deriv,
                'constant_deriv': True, 'latex': self.fluid_func_doc,
                'num_eq': self.num_nw_fluids * 4},
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

# %% Inlets and outlets

    def inlets(self):
        return ['in1', 'in2', 'in3']

    def outlets(self):
        return ['out1', 'out2']

# %% Equations and derivatives

    def preprocess(self, nw):

        if not self.P.is_set:
            self.set_attr(P='var')
            msg = ('The power output of a fuel cell must be set! '
                   'We are adding the power output of component ' +
                   self.label + ' as custom variable of the system.')
            logger.info(msg)

        for fluid in ['H2', 'H2O', 'O2']:
            try:
                setattr(
                    self, fluid, [x for x in nw.fluids if x in [
                        a.replace(' ', '') for a in
                        CP.get_aliases(fluid.upper())
                    ]][0])
            except IndexError:
                msg = (
                    'The component ' + self.label + ' (class ' +
                    self.__class__.__name__ + ') requires that the fluid '
                    '[fluid] is in the network\'s list of fluids.')
                aliases = ', '.join(CP.get_aliases(fluid.upper()))
                msg = msg.replace(
                    '[fluid]', fluid.upper() + ' (aliases: ' + aliases + ')')
                logger.error(msg)
                raise TESPyComponentError(msg)

        self.e0 = self.calc_e0()

        super().preprocess(nw)

    def calc_e0(self):
        r"""
        Calculate the specific energy output of the fuel cell.

        Returns
        -------
        val : float
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
        M = molar_masses[self.H2]
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
        # derivative for m_H2,in:
        self.jacobian[k, 2, 0] = -self.eta.val * self.e0
        # derivatives for variable P:
        if self.P.is_var:
            self.jacobian[k, 5 + self.P.var_pos, 0] = 1

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
        self.jacobian[k, 0, 0] = -(
            self.inl[0].h.val_SI - self.outl[0].h.val_SI
        )
        self.jacobian[k, 0, 2] = -self.inl[0].m.val_SI
        self.jacobian[k, 3, 2] = self.inl[0].m.val_SI

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
        self.jacobian[k, 2, 0] = -self.e.val
        # derivatives for variable P
        if self.P.is_var:
            self.jacobian[k, 5 + self.P.var_pos, 0] = 1
        # derivatives for variable e
        if self.e.is_var:
            self.jacobian[k, 5 + self.e.var_pos, 0] = -self.inl[2].m.val_SI

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
        h_refh2o = h_mix_pT([1, p_ref, 0, self.outl[1].fluid.val], T_ref)
        h_refh2 = h_mix_pT([1, p_ref, 0, self.inl[2].fluid.val], T_ref)
        h_refo2 = h_mix_pT([1, p_ref, 0, self.inl[1].fluid.val], T_ref)

        # derivatives cooling water inlet
        self.jacobian[k, 0, 0] = self.outl[0].h.val_SI - self.inl[0].h.val_SI
        self.jacobian[k, 0, 2] = -self.inl[0].m.val_SI

        # derivatives water outlet
        self.jacobian[k, 4, 0] = (self.outl[1].h.val_SI - h_refh2o)
        self.jacobian[k, 4, 2] = self.outl[1].m.val_SI

        # derivative cooling water outlet
        self.jacobian[k, 3, 2] = self.inl[0].m.val_SI

        # derivatives oxygen inlet
        self.jacobian[k, 1, 0] = -(self.inl[1].h.val_SI - h_refo2)
        self.jacobian[k, 1, 2] = -self.inl[1].m.val_SI

        # derivatives hydrogen inlet
        self.jacobian[k, 2, 0] = -(self.inl[2].h.val_SI - h_refh2 - self.e0)
        self.jacobian[k, 2, 2] = -self.inl[2].m.val_SI

        # derivatives for variable P
        if self.P.is_var:
            self.jacobian[k, 5 + self.P.var_pos, 0] = 1

    def fluid_func(self):
        r"""
        Equations for fluid composition.

        Returns
        -------
        residual : list
            Residual values of equation.

            .. math::

                0  = x_\mathrm{i,in,1} - x_\mathrm{i,out,1}
                \forall i \in \text{network fluids}\\
                0 = \begin{cases}
                    1 - x_\mathrm{i,in2} & \text{i=}H_{2}O\\
                    x_\mathrm{i,in2} & \text{else}
                \end{cases} \forall i \in \text{network fluids}\\
                0 = \begin{cases}
                    1 - x_\mathrm{i,out,2} & \text{i=}O_{2}\\
                    x_\mathrm{i,out,2} & \text{else}
                \end{cases} \forall i \in \text{network fluids}\\
                0 = \begin{cases}
                    1 - x_\mathrm{i,out,3} & \text{i=}H_{2}\\
                    x_\mathrm{i,out,3} & \text{else}
                \end{cases} \forall i \in \text{network fluids}
        """
        residual = []
        # equations for fluid composition in cooling loop
        for fluid, x in self.inl[0].fluid.val.items():
            residual += [x - self.outl[0].fluid.val[fluid]]

        # equations to constrain fluids to inlets/outlets
        residual += [1 - self.inl[1].fluid.val[self.O2]]
        residual += [1 - self.inl[2].fluid.val[self.H2]]
        residual += [1 - self.outl[1].fluid.val[self.H2O]]

        # equations to ban other fluids off inlets/outlets
        for fluid in self.inl[1].fluid.val.keys():
            if fluid != self.H2O:
                residual += [0 - self.outl[1].fluid.val[fluid]]
            if fluid != self.O2:
                residual += [0 - self.inl[1].fluid.val[fluid]]
            if fluid != self.H2:
                residual += [0 - self.inl[2].fluid.val[fluid]]

        return residual

    def fluid_func_doc(self, label):
        r"""
        Equations for fluid composition.

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
            r'0 = &x_\mathrm{i,in,1} - x_\mathrm{i,out,1}\\' + '\n'
            r'0 = &\begin{cases}' + '\n'
            r'1 - x_\mathrm{i,out,2} & \text{i=}H_{2}O\\' + '\n'
            r'x_\mathrm{i,out,2} & \text{else}\\' + '\n'
            r'\end{cases}\\' + '\n'
            r'0 =&\begin{cases}' + '\n'
            r'1 - x_\mathrm{i,in,2} & \text{i=}O_{2}\\' + '\n'
            r'x_\mathrm{i,in,2} & \text{else}\\' + '\n'
            r'\end{cases}\\' + '\n'
            r'0 =&\begin{cases}' + '\n'
            r'1 - x_\mathrm{i,in,3} & \text{i=}H_{2}\\' + '\n'
            r'x_\mathrm{i,in,3} & \text{else}\\' + '\n'
            r'\end{cases}\\' + '\n'
            r'&\forall i \in \text{network fluids}' + '\n'
            r'\end{split}')
        return generate_latex_eq(self, latex, label)

    def fluid_deriv(self):
        r"""
        Calculate the partial derivatives for cooling loop fluid balance.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the fluid equations.
        """
        # derivatives for cooling fluid composition
        deriv = np.zeros((
            self.num_nw_fluids * 4,
            5 + self.num_vars,
            self.num_nw_vars))

        k = 0
        for fluid, x in self.inl[0].fluid.val.items():
            deriv[k, 0, 3 + k] = 1
            deriv[k, 3, 3 + k] = -1
            k += 1

        # derivatives to constrain fluids to inlets/outlets
        i = 0
        for fluid in self.nw_fluids:
            if fluid == self.H2O:
                deriv[k, 4, 3 + i] = -1
            elif fluid == self.O2:
                deriv[k + 1, 1, 3 + i] = -1
            elif fluid == self.H2:
                deriv[k + 2, 2, 3 + i] = -1
            i += 1
        k += 3

        # derivatives to ban fluids off inlets/outlets
        i = 0
        for fluid in self.nw_fluids:
            if fluid != self.H2O:
                deriv[k, 4, 3 + i] = -1
                k += 1
            if fluid != self.O2:
                deriv[k, 1, 3 + i] = -1
                k += 1
            if fluid != self.H2:
                deriv[k, 2, 3 + i] = -1
                k += 1
            i += 1

        return deriv

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
        o2 = molar_masses[self.O2] / (
            molar_masses[self.O2] + 2 * molar_masses[self.H2])
        # equation for mass flow balance cooling water
        residual = []
        residual += [self.inl[0].m.val_SI - self.outl[0].m.val_SI]
        # equations for mass flow balance of the fuel cell
        residual += [o2 * self.outl[1].m.val_SI - self.inl[1].m.val_SI]
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

    def mass_flow_deriv(self):
        r"""
        Calculate the partial derivatives for all mass flow balance equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the mass flow equations.
        """
        # derivatives for mass flow balance in the heat exchanger
        deriv = np.zeros((3, 5 + self.num_vars, self.num_nw_vars))
        deriv[0, 0, 0] = 1
        deriv[0, 3, 0] = -1
        # derivatives for mass flow balance for oxygen input
        o2 = molar_masses[self.O2] / (
            molar_masses[self.O2] + 2 * molar_masses[self.H2])
        deriv[1, 4, 0] = o2
        deriv[1, 1, 0] = -1
        # derivatives for mass flow balance for hydrogen input
        deriv[2, 4, 0] = (1 - o2)
        deriv[2, 2, 0] = -1

        return deriv

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
            self.outl[1].p.val_SI - self.inl[2].p.val_SI]

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

    def reactor_pressure_deriv(self):
        r"""
        Calculate the partial derivatives for combustion pressure equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the pressure equations.
        """
        deriv = np.zeros((2, 5 + self.num_vars, self.num_nw_vars))
        # derivatives for pressure oxygen inlet
        deriv[0, 1, 1] = -1
        deriv[0, 4, 1] = 1
        # derivatives for pressure hydrogen inlet
        deriv[1, 2, 1] = -1
        deriv[1, 4, 1] = 1

        return deriv

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

        # equations to set a reference point for each h2o, h2 and o2
        h_refh2o = h_mix_pT([1, p_ref, 0, self.outl[1].fluid.val], T_ref)
        h_refh2 = h_mix_pT([1, p_ref, 0, self.inl[2].fluid.val], T_ref)
        h_refo2 = h_mix_pT([1, p_ref, 0, self.inl[1].fluid.val], T_ref)

        val = (
            self.inl[2].m.val_SI * (
                self.inl[2].h.val_SI - h_refh2 - self.e0
            )
            + self.inl[1].m.val_SI * (self.inl[1].h.val_SI - h_refo2)
            - self.inl[0].m.val_SI * (
                self.outl[0].h.val_SI - self.inl[0].h.val_SI
            )
            - self.outl[1].m.val_SI * (self.outl[1].h.val_SI - h_refh2o)
        )

        return val

    def initialise_fluids(self):
        # Set values to pure fluid on gas inlets and water outlet.
        self.inl[1].fluid.val[self.O2] = 1
        self.inl[2].fluid.val[self.H2] = 1
        self.outl[1].fluid.val[self.H2O] = 1
        self.inl[1].source.propagate_fluid_to_source(
            self.inl[1], self.inl[1].source)
        self.inl[2].source.propagate_fluid_to_source(
            self.inl[2], self.inl[2].source)
        self.outl[1].target.propagate_fluid_to_target(
            self.outl[1], self.outl[1].target)

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
            flow = c.get_flow()
            T = 20 + 273.15
            return h_mix_pT(flow, T)

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
            flow = c.get_flow()
            T = 50 + 273.15
            return h_mix_pT(flow, T)

    def propagate_fluid_to_target(self, inconn, start):
        r"""
        Propagate the fluids towards connection's target in recursion.

        Parameters
        ----------
        inconn : tespy.connections.connection.Connection
            Connection to initialise.

        start : tespy.components.component.Component
            This component is the fluid propagation starting point.
            The starting component is saved to prevent infinite looping.
        """
        if inconn == self.inl[0]:
            outconn = self.outl[0]

            for fluid, x in inconn.fluid.val.items():
                if (outconn.fluid.val_set[fluid] is False and
                        outconn.good_starting_values is False):
                    outconn.fluid.val[fluid] = x

            outconn.target.propagate_fluid_to_target(outconn, start)

    def propagate_fluid_to_source(self, outconn, start):
        r"""
        Propagate the fluids towards connection's source in recursion.

        Parameters
        ----------
        outconn : tespy.connections.connection.Connection
            Connection to initialise.

        start : tespy.components.component.Component
            This component is the fluid propagation starting point.
            The starting component is saved to prevent infinite looping.
        """
        if outconn == self.outl[0]:
            inconn = self.inl[0]

            for fluid, x in outconn.fluid.val.items():
                if (inconn.fluid.val_set[fluid] is False and
                        inconn.good_starting_values is False):
                    inconn.fluid.val[fluid] = x

            inconn.source.propagate_fluid_to_source(inconn, start)

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        self.Q.val = - self.inl[0].m.val_SI * (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI)
        self.pr.val = self.outl[0].p.val_SI / self.inl[0].p.val_SI
        self.e.val = self.P.val / self.inl[2].m.val_SI
        self.eta.val = self.e.val / self.e0

        i = self.inl[0].get_flow()
        o = self.outl[0].get_flow()
        self.zeta.val = (
            (i[1] - o[1]) * np.pi ** 2 / (
                4 * i[0] ** 2 *
                (self.inl[0].vol.val_SI + self.outl[0].vol.val_SI)
            )
        )
