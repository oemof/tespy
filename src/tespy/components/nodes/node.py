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
from tespy.tools.data_containers import DataContainerSimple as dc_simple
from tespy.tools.helpers import num_fluids


class Node(Component):
    r"""
    Class Node is the parent class for Splitter, Separator and Merge.

    Equations

        **mandatory equations**

        - :py:meth:`tespy.components.component.Component.mass_flow_func`

        .. math::

            0 = p_{in,1} - p_{out,i}\\
            \forall i \in \mathrm{outlets}

        **additional equations**

        - :py:meth:`tespy.components.nodes.splitter.additional_equations`
        - :py:meth:`tespy.components.nodes.separator.additional_equations`
        - :py:meth:`tespy.components.nodes.merge.additional_equations`

    Inlets/Outlets

        - specify number of outlets with :code:`num_in` (default value: 2)
        - specify number of outlets with :code:`num_out` (default value: 2)

    Image

        .. image:: _images/Node.svg
           :scale: 100 %
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

    design_path: str
        Path to the components design case.

    local_offdesign : boolean
        Treat this component in offdesign mode in a design calculation.

    local_design : boolean
        Treat this component in design mode in an offdesign calculation.

    char_warnings: boolean
        Ignore warnings on default characteristics usage for this component.

    printout: boolean
        Include this component in the network's results printout.

    num_in : float/tespy.tools.data_containers.dc_simple
        Number of inlets for this component, default value: 2.

    num_out : float/tespy.tools.data_containers.dc_simple
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

        Component.comp_init(self, nw)

        # number of mandatroy equations for
        # mass flow: 1
        # pressure: number of inlets + number of outlets - 1
        # the number of fluid and enthalpy equations is variable for a node!

        self.num_eq = 1 + self.num_i + self.num_o - 1

        # number of fluid equations
        self.num_eq += self.num_nw_fluids * self.num_o
        # number of enthalpy equations
        self.num_eq += self.num_o

        self.jacobian = np.zeros((
            self.num_eq,
            self.num_i + self.num_o + self.num_vars,
            self.num_nw_vars))

        self.residual = np.zeros(self.num_eq)
        self.jacobian[0:1] = self.mass_flow_deriv()
        self.jacobian[1:self.num_i + self.num_o] = self.pressure_deriv()

    def equations(self):
        r"""Calculate residual vector with results of equations."""
        k = 0
        ######################################################################
        # eqation for mass flow balance
        self.residual[k] = self.mass_flow_func()
        k += 1

        ######################################################################
        # equations for pressure
        inl = []
        if self.num_i > 1:
            inl = self.inl[1:]
        for c in inl + self.outl:
            self.residual[k] = self.inl[0].p.val_SI - c.p.val_SI
            k += 1

        ######################################################################
        # additional eqations
        self.additional_equations(k)

    def derivatives(self, increment_filter):
        r"""Calculate partial derivatives for given equations."""
        ######################################################################
        # derivatives for mass and pressure are static
        k = self.num_i + self.num_o
        ######################################################################
        # additional derivatives
        self.additional_derivatives(increment_filter, k)

    def additional_equations(self, k):
        r"""
        Calculate the residual vector of the additional equations.

        Equations

            **mandatroy equations**

            .. math::

                0 = \sum_i \left(\dot{m}_{i} \cdot x_{i,j}\right) - x_{o,j}
                \cdot  \sum_i \dot{m}_{i}\\
                \forall j \in \text{fluids}\\
                \forall o \in \text{outgoing mass flows}\\
                \text{i: incoming mass flows}

            .. math::

                0 = \sum_i \dot{m}_{i} \cdot h_{i}
                - \sum_o \dot{m}_{o} \cdot h_{o}
                \forall o \in \text{outgoing mass flows}\\
                \text{i: incoming mass flows}
        """
        ######################################################################
        # check for incoming/outgoing mass flows in inlets and outlets
        loc = 0
        # total incoming enthalpy
        h = 0
        # total incoming mass flow (constant within every iteration)
        self.m_inc = 0

        self.inc = []
        self.outg = []
        for c in self.inl:
            # incoming
            if c.m.val_SI >= 0:
                self.inc += [[c, loc]]
                self.m_inc += c.m.val_SI
                h += c.m.val_SI * c.h.val_SI
            # outgoing
            else:
                self.outg += [[c, loc]]
            loc += 1

        for c in self.outl:
            # inconming
            if c.m.val_SI < 0:
                self.inc += [[c, loc]]
                self.m_inc -= c.m.val_SI
                h -= c.m.val_SI * c.h.val_SI
            # outgoing
            else:
                self.outg += [[c, loc]]
            loc += 1

        self.jacobian = self.jacobian[0:k]
        self.jacobian = np.append(self.jacobian, np.zeros((
            len(self.outg) * (self.num_nw_fluids + 1),
            self.num_i + self.num_o + self.num_vars,
            self.num_nw_vars)), axis=0)

        ######################################################################
        # equations for fluid composition

        for fluid in self.nw_fluids:
            m = 0
            for i in self.inc:
                m += abs(i[0].m.val_SI) * i[0].fluid.val[fluid]
            for o in self.outg:
                self.residual[k] = m - o[0].fluid.val[fluid] * self.m_inc
                k += 1

        ######################################################################
        # equations for energy balance
        for o in self.outg:
            self.residual[k] = h - o[0].h.val_SI * self.m_inc
            k += 1

    def additional_derivatives(self, increment_filter, k):
        r"""Calculate partial derivatives for given additional equations."""
        ######################################################################
        # derivatives for fluid balance equations
        j = 0
        for fluid in self.nw_fluids:
            for o in self.outg:
                self.jacobian[k, o[1], j + 3] = -self.m_inc
                for i in self.inc:
                    self.jacobian[k, i[1], 0] = -i[0].fluid.val[fluid]
                    self.jacobian[k, i[1], j + 3] = -abs(i[0].m.val_SI)
                k += 1
            j += 1

        ######################################################################
        # derivatives for energy balance equations
        for o in self.outg:
            self.jacobian[k, o[1], 2] = -self.m_inc
            for i in self.inc:
                self.jacobian[k, i[1], 0] = i[0].h.val_SI - o[0].h.val_SI
                self.jacobian[k, i[1], 2] = abs(i[0].m.val_SI)
            k += 1

        logging.info(str(self.jacobian))

    def pressure_deriv(self):
        r"""
        Calculate partial derivatives for all pressure equations.

        Returns
        -------
        deriv : list
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

    def initialise_fluids(self):
        r"""
        Fluid initialisation for fluid mixture at outlet of the node.
        """
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

    @staticmethod
    def initialise_source(c, key):
        r"""
        Return a starting value for pressure and enthalpy at outlet.

        Parameters
        ----------
        c : tespy.connections.connection
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
        c : tespy.connections.connection
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
        Fluid propagation stops here.

        Parameters
        ----------
        outconn : tespy.connections.connection.Connection
            Connection to initialise.

        start : tespy.components.component.Component
            This component is the fluid propagation starting point.
            The starting component is saved to prevent infinite looping.
        """
        return

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
                'starting_point_property': 's',
                'starting_point_value': self.inc[i][0].s.val,
                'ending_point_property': 's',
                'ending_point_value': self.outg[0][0].s.val
            } for i in range(len(self.inc))}
