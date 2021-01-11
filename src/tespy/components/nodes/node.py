# -*- coding: utf-8

"""Module of class Node.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/nodes/node.py

SPDX-License-Identifier: MIT
"""

import logging

import numpy as np

from tespy.components.component import Component
from tespy.components.nodes.base import NodeBase
from tespy.tools.data_containers import DataContainerSimple as dc_simple
from tespy.tools.fluid_properties import s_mix_pT
from tespy.tools.helpers import num_fluids


class Node(NodeBase):
    r"""
    Class for merge points with multiple inflows and multiple outflows.

    This class is in experimental state, backflow is possible and taken into
    account. As backflow alters the system of equations of this component, the
    usage might be difficult.

    **Mandatory Equations**

    - :py:meth:`tespy.components.nodes.base.NodeBase.mass_flow_func`
    - :py:meth:`tespy.components.nodes.base.NodeBase.pressure_equality_func`
    - :py:meth:`tespy.components.nodes.node.Node.fluid_func`
    - :py:meth:`tespy.components.nodes.node.Node.energy_balance_func`

    Inlets/Outlets

    - specify number of outlets with :code:`num_in` (default value: 2)
    - specify number of outlets with :code:`num_out` (default value: 2)

    Image

    .. image:: _images/Node.svg
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

    num_in : float, tespy.tools.data_containers.DataContainerSimple
        Number of inlets for this component, default value: 2.

    num_out : float, tespy.tools.data_containers.DataContainerSimple
        Number of outlets for this component, default value: 2.

    Note
    ----
    - Node: Fluid composition and enthalpy at all **outgoing** connections
      (mass flow leaves the node) is result of mixture of the properties of
      the incoming connections (mass flow enters node).
      Incoming and outgoing connections can be a result of the calculation and
      are not identical to the inlets and outlets!
    - Splitter: Fluid composition and enthalpy at all outlets is the same as
      the inlet's properties.
    - Separator: Fluid composition is variable for all outlets, temperature at
      all outlets is the same as the inlet's temperature.
    - Merge: Fluid composition and enthalpy at outlet is result of mixture of
      the inlet's properties.

    Example
    -------
    The node can serve as merge or as splitter at the same time: The sum of
    all mass flow going into the node is identical to the mass flow leaving it.
    All incoming fluids are mixed in the node (mass flow weighted fluid mass
    fractions and enthalpy). All outgoing fluids have the composition of the
    mixture at the mixtures enthalpy/temperature.

    >>> from tespy.components import Sink, Source, Node
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> import shutil
    >>> import numpy as np
    >>> fluid_list = ['O2', 'N2']
    >>> nw = Network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ... h_unit='kJ / kg', iterinfo=False)
    >>> so1 = Source('source1')
    >>> so2 = Source('source2')
    >>> si1 = Sink('sink1')
    >>> si2 = Sink('sink2')
    >>> n = Node('node', num_in=2, num_out=2)
    >>> n.component()
    'node'
    >>> inc1 = Connection(so1, 'out1', n, 'in1')
    >>> inc2 = Connection(so2, 'out1', n, 'in2')
    >>> outg1 = Connection(n, 'out1', si1, 'in1')
    >>> outg2 = Connection(n, 'out2', si2, 'in1')
    >>> nw.add_conns(inc1, inc2, outg1, outg2)

    2 kg/s of pure oxygen is mixed with 5 kg/s mixture (50 % nitrogen and
    50 % oxygen). One of the outgoing connections should be at 3 kg/s.

    >>> inc1.set_attr(fluid={'O2': 1, 'N2': 0}, p=1, T=20, m=2)
    >>> inc2.set_attr(fluid={'O2': 0.5, 'N2': 0.5}, T=50, m=5)
    >>> outg1.set_attr(m=3)
    >>> nw.solve('design')
    >>> (round(outg1.fluid.val['O2'], 3), round(outg1.fluid.val['N2'], 3))
    (0.643, 0.357)
    >>> round(outg1.T.val, 1)
    41.8

    Now calculate at what mass flow of the 50/50 mixture the total
    oxygen fraction of the nodes outlets will be at 80 %.

    >>> inc2.set_attr(m=np.nan)
    >>> outg1.set_attr(fluid={'O2': 0.8})
    >>> nw.solve('design')
    >>> round(inc2.m.val_SI, 3)
    1.333
    """

    @staticmethod
    def component():
        return 'node'

    @staticmethod
    def attr():
        return {'num_in': dc_simple(), 'num_out': dc_simple()}

    def inlets(self):
        if self.num_in.is_set:
            return ['in' + str(i + 1) for i in range(self.num_in.val)]
        else:
            self.set_attr(num_in=2)
            return self.inlets()

    def outlets(self):
        if self.num_out.is_set:
            return ['out' + str(i + 1) for i in range(self.num_out.val)]
        else:
            self.set_attr(num_out=2)
            return self.outlets()

    def comp_init(self, nw):

        # number of mandatroy equations for
        # mass flow: 1
        # pressure: number of inlets + number of outlets - 1
        # the number of fluid and enthalpy equations is variable for a node,
        # initial values are:
        # fluid: number of fluids * number of outlets
        # energy: number of outlets
        num_eq = self.num_i + self.num_o * (2 + len(nw.fluids))
        Component.comp_init(self, nw, num_eq=num_eq)
        # constant derivatives
        self.jacobian[0:1] = self.mass_flow_deriv()
        pos = self.num_i + self.num_o
        self.jacobian[1:pos] = self.pressure_equality_deriv()

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
        k = NodeBase.mandatory_equations(self, doc=doc)
        ######################################################################
        # check for incoming/outgoing mass flows in inlets and outlets
        loc = 0
        # total incoming enthalpy
        self.h_inc = 0
        # total incoming mass flow (constant within every iteration)
        self.m_inc = 0

        self.inc = []
        self.outg = []
        for c in self.inl:
            # incoming
            if c.m.val_SI >= 0:
                self.inc += [[c, loc]]
                self.m_inc += c.m.val_SI
                self.h_inc += c.m.val_SI * c.h.val_SI
            # outgoing
            else:
                self.outg += [[c, c.loc]]
            loc += 1

        for c in self.outl:
            # inconming
            if c.m.val_SI < 0:
                self.inc += [[c, loc]]
                self.m_inc -= c.m.val_SI
                self.h_inc -= c.m.val_SI * c.h.val_SI
            # outgoing
            else:
                self.outg += [[c, loc]]
            loc += 1

        num_out = len(self.outg)

        self.jacobian = self.jacobian[0:k]
        self.jacobian = np.append(self.jacobian, np.zeros((
            num_out * (self.num_nw_fluids + 1),
            self.num_i + self.num_o + self.num_vars,
            self.num_nw_vars)), axis=0)

        ######################################################################
        # equations for fluid composition
        num_eq = self.num_nw_fluids * num_out
        self.residual[k:k + num_eq] = self.fluid_func()
        if doc:
            self.equation_docs[k:k + num_eq] = self.fluid_func(doc=doc)
        k += num_eq
        ######################################################################
        # equations for energy balance
        self.residual[k:k + num_out] = self.energy_balance_func()
        if doc:
            self.equation_docs[k:k + num_out] = (
                self.energy_balance_func(doc=doc))
        k += num_out
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
        k = NodeBase.mandatory_derivatives(self, increment_filter)
        num_out = len(self.outg)
        ######################################################################
        # derivatives for fluid balance equations
        self.fluid_deriv(increment_filter, k)
        k += self.num_nw_fluids * num_out
        ######################################################################
        # derivatives for energy balance equations
        self.energy_balance_deriv(increment_filter, k)
        k += num_out
        return k

    def fluid_func(self, doc=False):
        r"""
        Calculate the vector of residual values for fluid balance equations.

        Parameters
        ----------
        doc : boolean
            Return equation in LaTeX format instead of value.

        Returns
        -------
        residual : list
            Vector of residual values for component's fluid balance.

            .. math::

                0 = \sum_i \left(|\dot{m}_{inc,i}| \cdot x_{fl,inc,i}\right)
                - x_{fl,outg,j} \cdot  \sum_i |\dot{m}_{inc,i}|\\
                \forall fl \in \text{fluids}\\
                \forall i \in \text{incomning mass flows}\\
                \forall j \in \text{outgoing mass flows}
        """
        if not doc:
            residual = []
            for fluid in self.nw_fluids:
                m = 0
                for i in self.inc:
                    m += abs(i[0].m.val_SI) * i[0].fluid.val[fluid]
                for o in self.outg:
                    residual += [m - o[0].fluid.val[fluid] * self.m_inc]
            return residual
        else:
            latex = (
                r'\begin{split}' + '\n'
                r'0 =& \sum_i \left(|\dot{m}_{\mathrm{inc,}i}| \cdot '
                r'x_{fl\mathrm{,inc,}i}\right) - x_{fl\mathrm{,outg,}j} \cdot '
                r'\sum_i |\dot{m}_{\mathrm{inc,}i}|\\' + '\n'
                r'\forall & fl \in \text{fluids}\\' + '\n'
                r'\forall & i \in \text{incomning mass flows}\\' + '\n'
                r'\forall & j \in \text{outgoing mass flows}\\' + '\n'
                r'\end{split}'
            )
            return (
                [self.generate_latex(latex, 'fluid_func')] +
                (self.num_nw_fluids - 1) * [''])

    def fluid_deriv(self, increment_filter, k):
        r"""
        Calculate partial derivatives of fluid balance.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        j = 0
        for fluid in self.nw_fluids:
            for o in self.outg:
                self.jacobian[k, o[1], j + 3] = -self.m_inc
                for i in self.inc:
                    self.jacobian[k, i[1], 0] = -i[0].fluid.val[fluid]
                    self.jacobian[k, i[1], j + 3] = -abs(i[0].m.val_SI)
                k += 1
            j += 1

    def energy_balance_func(self, doc=False):
        r"""
        Calculate energy balance.

        Parameters
        ----------
        doc : boolean
            Return equation in LaTeX format instead of value.

        Returns
        -------
        residual : list
            Residual value of energy balance.

            .. math::

                0 = \sum_i \left(|\dot{m}_{inc,i}| \cdot h_{inc,i}\right)
                - h_{outg,j} \cdot  \sum_i |\dot{m}_{inc,i}|\\
                \forall i \in \text{incomning mass flows}\\
                \forall j \in \text{outgoing mass flows}
        """
        if not doc:
            residual = []
            for o in self.outg:
                residual += [self.h_inc - o[0].h.val_SI * self.m_inc]
            return residual
        else:
            latex = (
                r'\begin{split}' + '\n'
                r'0 =& \sum_i \left(|\dot{m}_{\mathrm{inc,}i}| \cdot '
                r'h_{\mathrm{inc,}i}\right) - h_{\mathrm{outg,}j} \cdot '
                r'\sum_i |\dot{m}_{\mathrm{inc,}i}|\\' + '\n'
                r'\forall & i \in \text{incomning mass flows}\\' + '\n'
                r'\forall & j \in \text{outgoing mass flows}\\' + '\n'
                r'\end{split}'
            )
            return [self.generate_latex(latex, 'energy_balance_func')]

    def energy_balance_deriv(self, increment_filter, k):
        r"""
        Calculate partial derivatives of energy balance.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        for o in self.outg:
            self.jacobian[k, o[1], 2] = -self.m_inc
            for i in self.inc:
                self.jacobian[k, i[1], 0] = i[0].h.val_SI - o[0].h.val_SI
                self.jacobian[k, i[1], 2] = abs(i[0].m.val_SI)
            k += 1

    def initialise_fluids(self):
        """Fluid initialisation for fluid mixture at outlet of the node."""
        num_fl = {}
        for o in self.outl:
            num_fl[o] = num_fluids(o.fluid.val)

        for i in self.inl:
            num_fl[i] = num_fluids(i.fluid.val)

        ls = []
        if any(num_fl.values()) and not all(num_fl.values()):
            for conn, num in num_fl.items():
                if num == 1:
                    ls += [conn]

            for c in ls:
                for fluid in self.nw_fluids:
                    for o in self.outl:
                        if not o.fluid.val_set[fluid]:
                            o.fluid.val[fluid] = c.fluid.val[fluid]
                    for i in self.inl:
                        if not i.fluid.val_set[fluid]:
                            i.fluid.val[fluid] = c.fluid.val[fluid]
            for o in self.outl:
                o.target.propagate_fluid_to_target(o, o.target)

    def propagate_fluid_to_target(self, inconn, start):
        r"""
        Fluid propagation stops here.

        Parameters
        ----------
        inconn : tespy.connections.connection.Connection
            Connection to initialise.

        start : tespy.components.component.Component
            This component is the fluid propagation starting point.
            The starting component is saved to prevent infinite looping.
        """
        return

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
        for inconn in self.inl:
            for fluid, x in outconn.fluid.val.items():
                if (not inconn.fluid.val_set[fluid] and
                        not inconn.good_starting_values):
                    inconn.fluid.val[fluid] = x

            inconn.source.propagate_fluid_to_source(inconn, start)

    def entropy_balance(self):
        r"""
        Calculate entropy balance of a node.

        Note
        ----
        A definition of reference points is included for compensation of
        differences in zero point definitions of different fluid compositions.

        - Reference temperature: 298.15 K.
        - Reference pressure: 1 bar.

        .. math::

            \dot{S}_\mathrm{irr}= \sum_{i} \dot{m}_{\mathrm{outg,}i} \cdot
            \left( s_{\mathrm{outg,}i} - s_{\mathrm{outg,ref,}i} \right)
            - \sum_{i} \dot{m}_{\mathrm{inc,}i} \cdot
            \left( s_{\mathrm{inc,}i} - s_{\mathrm{inc,ref,}i} \right)\\
        """
        # check if reference point definition is necessary (it should not)
        T_ref = 298.15
        p_ref = 1e5
        self.S_irr = 0
        for o in self.outg:
            self.S_irr += o[0].m.val_SI * (
                o[0].s.val_SI -
                s_mix_pT([0, p_ref, 0, o[0].fluid.val], T_ref))
        for i in self.inc:
            self.S_irr -= i[0].m.val_SI * (
                i[0].s.val_SI -
                s_mix_pT([0, p_ref, 0, i[0].fluid.val], T_ref))

    def exergy_balance(self, T0):
        r"""
        Calculate exergy balance of a node.

        Parameters
        ----------
        T0 : float
            Ambient temperature T0 / K.

        Note
        ----
        Please note, that the exergy balance accounts for physical exergy only.

        .. math::

            \dot{E}_\mathrm{P} = \sum E_{\mathrm{out,}i}^\mathrm{PH}\\
            \dot{E}_\mathrm{F} = \sum E_{\mathrm{in,}j}^\mathrm{PH}
        """
        self.E_P = 0
        self.E_F = 0
        for i in self.inc:
            self.E_F += i[0].Ex_physical

        for o in self.outg:
            self.E_P += o[0].Ex_physical

        self.E_bus = np.nan
        self.E_D = self.E_F - self.E_P
        self.epsilon = self.E_P / self.E_F

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
            i + 1: {
                'isoline_property': 'p',
                'isoline_value': self.inc[i][0].p.val,
                'isoline_value_end': self.outg[0][0].p.val,
                'starting_point_property': 'v',
                'starting_point_value': self.inc[i][0].vol.val,
                'ending_point_property': 'v',
                'ending_point_value': self.outg[0][0].vol.val
            } for i in range(len(self.inc))}
