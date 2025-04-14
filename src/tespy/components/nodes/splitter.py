# -*- coding: utf-8

"""Module of class Splitter.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/nodes/splitter.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import component_registry
from tespy.components.nodes.base import NodeBase
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
from tespy.tools.document_models import generate_latex_eq


@component_registry
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

    .. image:: /api/_images/Splitter.svg
       :alt: flowsheet of the splitter
       :align: center
       :class: only-light

    .. image:: /api/_images/Splitter_darkmode.svg
       :alt: flowsheet of the splitter
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
    >>> nw = Network(p_unit='bar', T_unit='C', iterinfo=False)
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
    def get_parameters():
        return {'num_out': dc_simple()}

    def get_mandatory_constraints(self):
        return {
            'mass_flow_constraints': {
                'func': self.mass_flow_func, 'deriv': self.mass_flow_deriv,
                'constant_deriv': True, 'latex': self.mass_flow_func_doc,
                'num_eq': 1},
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

    def propagate_wrapper_to_target(self, branch):
        branch["components"] += [self]
        for outconn in self.outl:
            branch["connections"] += [outconn]
            outconn.target.propagate_wrapper_to_target(branch)

    def preprocess(self, num_nw_vars):
        super().preprocess(num_nw_vars)
        self._propagation_start = False

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

    def energy_balance_deriv(self, increment_filter, k):
        r"""
        Calculate partial derivatives for energy balance equation.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        for eq, o in enumerate(self.outl):
            if self.inl[0].h.is_var:
                self.jacobian[k + eq, self.inl[0].h.J_col] = 1
            if o.h.is_var:
                self.jacobian[k + eq, o.h.J_col] = -1
