# -*- coding: utf-8

"""Module of class NodeBase.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/nodes/base.py

SPDX-License-Identifier: MIT
"""

import numpy as np

from tespy.components.component import Component


class NodeBase(Component):
    """Class NodeBase is parent class for all components of submodule nodes."""

    def mandatory_equations(self, doc=False):
        r"""
        Calculate residual vector of mandatory equations.

        Parameters
        ----------
        doc : boolean
            Return equation in LaTeX format instead of value.

        Returns
        -------
        k : int
            Position of last equation in residual value vector (k-th equation).
        """
        k = 0
        ######################################################################
        # eqation for mass flow balance
        self.residual[k] = self.mass_flow_func()
        if doc:
            self.equation_docs[k:k + 1] = self.mass_flow_func(doc=doc)
        k += 1

        ######################################################################
        # eqation for pressure balance
        num_eq = self.num_i + self.num_o - 1
        self.residual[k:k + num_eq] = self.pressure_equality_func()
        if doc:
            self.equation_docs[k:k + num_eq] = (
                self.pressure_equality_func(doc=doc))
        k += num_eq
        return k

    def mandatory_derivatives(self, increment_filter):
        r"""
        Calculate partial derivatives for mandatory equations.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        Returns
        -------
        k : int
            Position of last equation in residual value vector (k-th equation).
        """
        ######################################################################
        # derivatives are static
        return self.num_i + self.num_o

    def mass_flow_func(self, doc=False):
        r"""
        Calculate the residual value for mass flow balance equation.

        Parameters
        ----------
        doc : boolean
            Return equation in LaTeX format instead of value.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                0 = \sum \dot{m}_{in,i} - \sum \dot{m}_{out,j} \;
                \forall i \in inlets, \forall j \in outlets
        """
        if not doc:
            res = 0
            for i in self.inl:
                res += i.m.val_SI
            for o in self.outl:
                res -= o.m.val_SI
            return res
        else:
            latex = (
                r'0 =\sum\dot{m}_{\mathrm{in},i}-\sum\dot{m}_{\mathrm{out},j}'
                r'\;\forall i \in \text{inlets}, \forall j \in \text{outlets}')
            return [self.generate_latex(latex, 'mass_flow_func')]

    def mass_flow_deriv(self):
        r"""
        Calculate partial derivatives for mass flow equation.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((1, self.num_i + self.num_o, self.num_nw_vars))
        for i in range(self.num_i):
            deriv[0, i, 0] = 1
        for j in range(self.num_o):
            deriv[0, j + i + 1, 0] = -1
        return deriv

    def pressure_equality_func(self, doc=False):
        r"""
        Calculate the residual values of pressure equality equations.

        Parameters
        ----------
        doc : boolean
            Return equation in LaTeX format instead of value.

        Returns
        -------
        residual : list
            Vector with residual value for pressure equality equations.

            .. math::

                0 = p_{in,1} - p_{in,i}\forall i \in \text{inlets > 1}\\
                0 = p_{in,1} - p_{out,j}\forall j \in \text{outlets}
        """
        if not doc:
            residual = []
            inl = []
            if self.num_i > 1:
                inl = self.inl[1:]
            for c in inl + self.outl:
                residual += [self.inl[0].p.val_SI - c.p.val_SI]
            return residual
        else:
            latex = (
                r'\begin{split}' + '\n'
                r'0 = p_\mathrm{in,1} - p_{\mathrm{in,}i} '
                r'& \; \forall i \in \text{inlets} \setminus '
                r'\left\lbrace 1\right\rbrace\\' + '\n'
                r'0 = p_\mathrm{in,1} - p_{\mathrm{out,}j} '
                r'& \; \forall j \in \text{outlets}\\' + '\n'
                r'\end{split}'
            )
            return (
                [self.generate_latex(latex, 'pressure_equality_func')] +
                (self.num_i + self.num_o - 2) * [''])

    def pressure_equality_deriv(self):
        r"""
        Calculate partial derivatives for all pressure equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((
            self.num_i + self.num_o - 1,
            self.num_i + self.num_o,
            self.num_nw_vars))

        inl = []
        if self.num_i > 1:
            inl = self.inl[1:]
        for k in range(len(inl + self.outl)):
            deriv[k, 0, 1] = 1
            deriv[k, k + 1, 1] = -1
        return deriv

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
