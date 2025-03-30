# -*- coding: utf-8

"""Module of class NodeBase.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/nodes/base.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import Component
from tespy.components.component import component_registry
from tespy.tools.document_models import generate_latex_eq


@component_registry
class NodeBase(Component):
    """Class NodeBase is parent class for all components of submodule nodes."""

    @staticmethod
    def get_bypass_constraints():
        return {}

    def mass_flow_func(self):
        r"""
        Calculate the residual value for mass flow balance equation.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                0 = \sum \dot{m}_{in,i} - \sum \dot{m}_{out,j} \;
                \forall i \in inlets, \forall j \in outlets
        """
        res = 0
        for i in self.inl:
            res += i.m.val_SI
        for o in self.outl:
            res -= o.m.val_SI
        return res

    def mass_flow_func_doc(self, label):
        r"""
        Calculate the residual value for mass flow balance equation.

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
            r'0 =\sum\dot{m}_{\mathrm{in},i}-\sum\dot{m}_{\mathrm{out},j}'
            r'\;\forall i \in \text{inlets}, \forall j \in \text{outlets}')
        return generate_latex_eq(self, latex, label)

    def mass_flow_deriv(self, increment_filter, k):
        r"""
        Calculate partial derivatives for mass flow equation.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        for i in self.inl:
            if i.m.is_var:
                self.jacobian[k, i.m.J_col] = 1
        for o in self.outl:
            if o.m.is_var:
                self.jacobian[k, o.m.J_col] = -1

    def pressure_equality_func(self):
        r"""
        Calculate the residual values of pressure equality equations.

        Returns
        -------
        residual : list
            Vector with residual value for pressure equality equations.

            .. math::

                0 = p_{in,1} - p_{in,i}\forall i \in \text{inlets > 1}\\
                0 = p_{in,1} - p_{out,j}\forall j \in \text{outlets}
        """
        residual = []
        inl = []
        if self.num_i > 1:
            inl = self.inl[1:]
        for c in inl + self.outl:
            residual += [self.inl[0].p.val_SI - c.p.val_SI]
        return residual

    def pressure_equality_func_doc(self, label):
        r"""
        Calculate the residual values of pressure equality equations.

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
            r'0 = p_\mathrm{in,1} - p_{\mathrm{in,}i} '
            r'& \; \forall i \in \text{inlets} \setminus '
            r'\left\lbrace 1\right\rbrace\\' + '\n'
            r'0 = p_\mathrm{in,1} - p_{\mathrm{out,}j} '
            r'& \; \forall j \in \text{outlets}\\' + '\n'
            r'\end{split}'
        )
        return generate_latex_eq(self, latex, label)

    def pressure_equality_deriv(self, increment_filter, k):
        r"""
        Calculate partial derivatives for all pressure equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the fluid equations.
        """
        if self.num_i > 1:
            conns = self.inl[1:] + self.outl
        else:
            conns = self.outl

        for eq, o in enumerate(conns):
            if self.inl[0].p.is_var:
                self.jacobian[k + eq, self.inl[0].p.J_col] = 1
            if o.p.is_var:
                self.jacobian[k + eq, o.p.J_col] = -1

    @staticmethod
    def initialise_source(c, key):
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
                10^5 & \text{key = 'p'}\\
                5 \cdot 10^5 & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 1e5
        elif key == 'h':
            return 5e5

    @staticmethod
    def initialise_target(c, key):
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
                10^5 & \text{key = 'p'}\\
                5 \cdot 10^5 & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 1e5
        elif key == 'h':
            return 5e5

    def propagate_to_target(self, branch):

        for outconn in self.outl:
            subbranch = {
                "connections": [outconn],
                "components": [self, outconn.target],
                "subbranches": {}
            }
            outconn.target.propagate_to_target(subbranch)
            branch["subbranches"][outconn.label] = subbranch
