# -*- coding: utf-8

"""Module of class Splitter.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/nodes/splitter.py

SPDX-License-Identifier: MIT
"""

import numpy as np

from tespy.components.nodes.base import NodeBase
from tespy.tools.data_containers import DataContainerSimple as dc_simple
from tespy.tools.document_models import generate_latex_eq


class Splitter(NodeBase):
    r"""
    Split up a mass flow in parts of identical enthalpy and fluid composition.

    **Mandatory Equations**

    - :py:meth:`tespy.components.nodes.base.NodeBase.mass_flow_func`
    - :py:meth:`tespy.components.nodes.base.NodeBase.pressure_equality_func`
    - :py:meth:`tespy.components.nodes.splitter.Splitter.fluid_func`
    - :py:meth:`tespy.components.nodes.splitter.Splitter.energy_balance_func`

    Inlets/Outlets

    - in1
    - specify number of outlets with :code:`num_out` (default value: 2)

    Image

    .. image:: _images/Splitter.svg
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

    num_out : float, dict
        Number of outlets for this component, default value: 2.

    Example
    -------
    A splitter is used to split up a single mass flow into a specified number
    of different parts at identical pressure, enthalpy and fluid composition.

    >>> from tespy.components import Sink, Source, Splitter
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> import shutil
    >>> import numpy as np
    >>> fluid_list = ['O2', 'N2']
    >>> nw = Network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ... iterinfo=False)
    >>> so = Source('source')
    >>> si1 = Sink('sink1')
    >>> si2 = Sink('sink2')
    >>> si3 = Sink('sink3')
    >>> s = Splitter('splitter', num_out=3)
    >>> s.component()
    'splitter'
    >>> inc = Connection(so, 'out1', s, 'in1')
    >>> outg1 = Connection(s, 'out1', si1, 'in1')
    >>> outg2 = Connection(s, 'out2', si2, 'in1')
    >>> outg3 = Connection(s, 'out3', si3, 'in1')
    >>> nw.add_conns(inc, outg1, outg2, outg3)

    An Air (simplified) mass flow is split up into three mass flows. The total
    incoming mass flow is 5 kg/s, 3 kg/s and 1 kg/s respectively are leaving
    the splitter into the first two outlets. The residual mass flow will
    drain in the last outlet. Temperature and fluid composition will not
    change.

    >>> inc.set_attr(fluid={'O2': 0.23, 'N2': 0.77}, p=1, T=20, m=5)
    >>> outg1.set_attr(m=3)
    >>> outg2.set_attr(m=1)
    >>> nw.solve('design')
    >>> round(outg3.m.val_SI, 1)
    1.0
    >>> round(inc.T.val, 1)
    20.0
    >>> round(outg3.T.val, 1)
    20.0
    """

    @staticmethod
    def component():
        return 'splitter'

    @staticmethod
    def get_variables():
        return {'num_out': dc_simple()}

    def get_mandatory_constraints(self):
        return {
            'mass_flow_constraints': {
                'func': self.mass_flow_func, 'deriv': self.mass_flow_deriv,
                'constant_deriv': True, 'latex': self.mass_flow_func_doc,
                'num_eq': 1},
            'fluid_constraints': {
                'func': self.fluid_func, 'deriv': self.fluid_deriv,
                'constant_deriv': True, 'latex': self.fluid_func_doc,
                'num_eq': self.num_o * self.num_nw_fluids},
            'energy_balance_constraints': {
                'func': self.energy_balance_func,
                'deriv': self.energy_balance_deriv,
                'constant_deriv': True, 'latex': self.energy_balance_func_doc,
                'num_eq': self.num_o},
            'pressure_constraints': {
                'func': self.pressure_equality_func,
                'deriv': self.pressure_equality_deriv,
                'constant_deriv': True,
                'latex': self.pressure_equality_func_doc,
                'num_eq': self.num_i + self.num_o - 1}
        }

    @staticmethod
    def inlets():
        return ['in1']

    def outlets(self):
        if self.num_out.is_set:
            return ['out' + str(i + 1) for i in range(self.num_out.val)]
        else:
            self.set_attr(num_out=2)
            return self.outlets()

    def fluid_func(self):
        r"""
        Calculate the vector of residual values for fluid balance equations.

        Returns
        -------
        residual : list
            Vector of residual values for component's fluid balance.

            .. math::

                0 = x_{fl,in} - x_{fl,out,j} \;
                \forall fl \in \text{network fluids,} \; \forall j \in
                \text{outlets}
        """
        residual = []
        for o in self.outl:
            for fluid, x in self.inl[0].fluid.val.items():
                residual += [x - o.fluid.val[fluid]]
        return residual

    def fluid_func_doc(self, label):
        r"""
        Calculate the vector of residual values for fluid balance equations.

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
            r'0 = x_{fl\mathrm{,in}} - x_{fl\mathrm{,out,}j}'
            r'\; \forall fl \in \text{network fluids,} \; \forall j \in'
            r'\text{outlets}'
        )
        return generate_latex_eq(self, latex, label)

    def fluid_deriv(self):
        r"""
        Calculate partial derivatives for all fluid balance equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((self.num_nw_fluids * self.num_o, 1 + self.num_o,
                          self.num_nw_vars))
        k = 0
        for o in self.outl:
            i = 0
            for fluid in self.nw_fluids:
                deriv[i + k * self.num_nw_fluids, 0, i + 3] = 1
                deriv[i + k * self.num_nw_fluids, k + 1, i + 3] = -1
                i += 1
            k += 1
        return deriv

    def energy_balance_func(self):
        r"""
        Calculate energy balance.

        Returns
        -------
        residual : list
            Residual value of energy balance.

            .. math::

                0 = h_{in} - h_{out,j} \;
                \forall j \in \mathrm{outlets}\\
        """
        residual = []
        for o in self.outl:
            residual += [self.inl[0].h.val_SI - o.h.val_SI]
        return residual

    def energy_balance_func_doc(self, label):
        r"""
        Calculate energy balance.

        Parameters
        ----------
        label : str
            Label for equation.
        """
        latex = r'0=h_{in}-h_{\mathrm{out,}j}\;\forall j \in\text{outlets}'
        return generate_latex_eq(self, latex, label)

    def energy_balance_deriv(self):
        r"""
        Calculate partial derivatives for energy balance equation.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        deriv = np.zeros((self.num_o, 1 + self.num_o, self.num_nw_vars))
        k = 0
        for o in self.outl:
            deriv[k, 0, 2] = 1
            deriv[k, k + 1, 2] = -1
            k += 1
        return deriv

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
        for outconn in self.outl:
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
        inconn = self.inl[0]
        for fluid, x in outconn.fluid.val.items():
            if (inconn.fluid.val_set[fluid] is False and
                    inconn.good_starting_values is False):
                inconn.fluid.val[fluid] = x

        inconn.source.propagate_fluid_to_source(inconn, start)
