# -*- coding: utf-8

"""Module of class NodeBase.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/nodes/base.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import Component
from tespy.components.component import component_registry


@component_registry
class NodeBase(Component):
    """
    Class NodeBase is parent class for all components of submodule nodes.

    Mandatory Equations
    -------------------

    - mass flow equality constraint(s): :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
    - fluid composition equality constraint(s): :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`

    Parameters
    ----------

    char_warnings : bool
        Ignore warnings on default characteristics usage for this component.

    design : list
        List containing design parameters (stated as String).

    design_path : str
        Path to the components design case.

    label : str
        The label of the component.

    local_design : bool
        Treat this component in design mode in an offdesign calculation.

    local_offdesign : bool
        Treat this component in offdesign mode in a design calculation.

    offdesign : list
        List containing offdesign parameters (stated as String).

    printout : bool
        Include this component in the network's results printout.
    """

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

    def mass_flow_dependents(self):
        return [c.m for c in self.inl + self.outl]

    def pressure_structure_matrix(self, k):
        r"""
        Set up structure matrix for all pressure constraints representing the
        following equations:

        .. math::

            0 = p_{in,1} - p_{in,i}\forall i \in \text{inlets > 1}\\
            0 = p_{in,1} - p_{out,j}\forall j \in \text{outlets}
        """
        if self.num_i > 1:
            conns = self.inl[1:] + self.outl
        else:
            conns = self.outl

        for eq, conn in enumerate(conns):
            self._structure_matrix[k + eq, self.inl[0].p.sm_col] = 1
            self._structure_matrix[k + eq, conn.p.sm_col] = -1
