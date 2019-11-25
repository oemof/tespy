# -*- coding: utf-8

"""Module for basic components.

Components in this module:

    - :func:`tespy.components.basics.source`
    - :func:`tespy.components.basics.sink`
    - :func:`tespy.components.basics.subsystem_interface`
    - :func:`tespy.components.basics.cycle_closer`


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/basics.py

SPDX-License-Identifier: MIT
"""

import numpy as np

from tespy.components.components import component

from tespy.tools.data_containers import dc_simple

# %%


class source(component):
    r"""
    A flow originates from a source.

    Equations
        This component is unconstrained.

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    Example
    -------
    Create a source and specify a label.

    >>> from tespy.components.basics import source
    >>> so = source('a labeled source')
    >>> so.component()
    'source'
    >>> so.label
    'a labeled source'
    """

    def component(self):
        return 'source'

    def outlets(self):
        return ['out1']

# %%


class sink(component):
    r"""
    A flow drains in a sink.

    Equations
        This component is unconstrained.

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    Example
    -------
    Create a sink and specify a label.

    >>> from tespy.components.basics import sink
    >>> si = sink('a labeled sink')
    >>> si.component()
    'sink'
    >>> si.label
    'a labeled sink'
    """

    def component(self):
        return 'sink'

    def inlets(self):
        return ['in1']

# %%


class subsystem_interface(component):
    r"""
    Equations

        **mandatory equations**

        .. math:: 0 = fluid_{i,in_{j}} - fluid_{i,out_{j}} \;
            \forall i \in \mathrm{fluid}, \; \forall j \in inlets/outlets

        .. math:: 0 = \dot{m}_{in_{j}} - \dot{m}_{out_{j}} \;
            \forall j \in inlets/outlets

        .. math:: 0 = p_{in_{j}} - p_{out_{j}} \;
            \forall j \in inlets/outlets

        .. math:: 0 = h_{in_{j}} - h_{out_{j}} \;
            \forall j \in inlets/outlets

    Inlets/Outlets

        - Specify number of inlets and outlets with :code:`num_inter`,
          predefined value: 1.

    Image

        .. image:: _images/subsys_interface.svg
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

    num_inter : float/tespy.helpers.dc_simple
        Number of interfaces for subsystem.

    Note
    ----
    This component passes all fluid properties and mass flow from its inlet to
    the outlet.

    Example
    -------
    As connections can only connect a component with a different
    component, the subsystem interface is used to connect subsystems with the
    rest of your network. It is necessary to specify the number of interfaces
    of the subsystem interface, if you want any number other than 1. We will
    not go in depth of subsystem usage in this example. Please refer to
    TODO: PLACELINKHERE for more information on building your own subsystems.

    >>> from tespy.components.basics import sink, source, subsystem_interface
    >>> from tespy.connections import connection
    >>> from tespy.networks.networks import network
    >>> fluids = ['H2O', 'N2']
    >>> nw = network(fluids=fluids)
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')
    >>> nw.set_printoptions(print_level='none')
    >>> so = source('source 1')
    >>> si = sink('sink 1')
    >>> IF = subsystem_interface('subsystem interface')
    >>> IF.component()
    'subsystem interface'
    >>> len(IF.inlets())
    1

    Add the subsystem the interface to a minimal network containing only a
    source and a sink. The interface does not change the fluid properties
    in any way.

    >>> inc = connection(so, 'out1', IF, 'in1')
    >>> outg = connection(IF, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)
    >>> inc.set_attr(fluid={'H2O': 1, 'N2': 0}, T=40, p=3, m=100)
    >>> nw.solve('design')
    >>> inc.m.val_SI == outg.m.val_SI
    True
    >>> inc.h.val_SI == outg.h.val_SI
    True
    """

    def component(self):
        return 'subsystem interface'

    def attr(self):
        return {'num_inter': dc_simple()}

    def inlets(self):
        if self.num_inter.val_set:
            return ['in' + str(i + 1) for i in range(self.num_inter.val)]
        else:
            return ['in1']

    def outlets(self):
        if self.num_inter.val_set:
            return ['out' + str(i + 1) for i in range(self.num_inter.val)]
        else:
            return ['out1']

    def comp_init(self, nw):

        component.comp_init(self, nw)

        # retrieve always constant derivatives
        self.fl_deriv = self.fluid_deriv()
        self.m_deriv = self.inout_deriv(0)
        self.p_deriv = self.inout_deriv(1)
        self.h_deriv = self.inout_deriv(2)

    def equations(self):
        r"""
        Calculates vector vec_res with results of equations for this component.

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # eqations for fluids
        for i in range(self.num_i):
            for fluid, x in self.inl[i].fluid.val.items():
                vec_res += [x - self.outl[i].fluid.val[fluid]]

        ######################################################################
        # equations for mass flow
        for i in range(self.num_i):
            vec_res += [self.inl[i].m.val_SI - self.outl[i].m.val_SI]

        ######################################################################
        # equations for pressure
        for i in range(self.num_i):
            vec_res += [self.inl[i].p.val_SI - self.outl[i].p.val_SI]

        ######################################################################
        # equations for enthalpy
        for i in range(self.num_i):
            vec_res += [self.inl[i].h.val_SI - self.outl[i].h.val_SI]

        ######################################################################

        return vec_res

    def derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given equations.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        ######################################################################
        # derivatives with constant value (all for this component)
        mat_deriv = self.fl_deriv + self.m_deriv + self.p_deriv + self.h_deriv

        return np.asarray(mat_deriv)

    def fluid_deriv(self):
        r"""
        Calculates the partial derivatives for all fluid balance equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((
                self.num_fl * self.num_i, 2 * self.num_i, 3 + self.num_fl))
        for i in range(self.num_i):
            for j in range(self.num_fl):
                deriv[i * self.num_fl + j, i, j + 3] = 1
                deriv[i * self.num_fl + j, self.num_i + i, j + 3] = -1
        return deriv.tolist()

    def inout_deriv(self, pos):
        r"""
        Calculates the partial derivatives for all mass flow, pressure and
        enthalpy equations.

        Parameters
        ----------
        pos : int
            Position of the variable in the matrix of derivatives.
            mass flow: 0, pressure: 1, enthalpy: 2.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((self.num_i, 2 * self.num_i, self.num_fl + 3))
        for i in range(self.num_i):
            deriv[i, i, pos] = 1
        for j in range(self.num_i):
            deriv[j, j + self.num_i, pos] = -1
        return deriv.tolist()


# %%


class cycle_closer(component):
    r"""
    Equations

        **mandatory equations**

        .. math::

            0 = p_{in} - p_{out}

            0 = h_{in} - h_{out}

    Image not available

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    Note
    ----
    This component can be used to close a cycle process. The system of
    equations describing your plant will overdetermined, if you close a cycle
    without this component or a cut the cycle with a sink and a source at
    some point of the cycle. This component can be used instead of cutting
    the cycle.

    Example
    -------
    Create a cycle containing a pump and a pipe. The pump increases pressure
    the pipe cools the liquid and destroys the pressure rise. The heat
    extracted at the pipe must be the same value of the power input at the
    pump (but negative), as there is no other in- or outputs of energy in the
    system.

    >>> from tespy.components.basics import cycle_closer
    >>> from tespy.components.piping import pipe
    >>> from tespy.components.turbomachinery import pump
    >>> from tespy.connections import connection
    >>> from tespy.networks.networks import network
    >>> nw = network(['water'], p_unit='bar', T_unit='C')
    >>> pi = pipe('pipe')
    >>> pu = pump('pump')
    >>> cc = cycle_closer('cycle closing component')
    >>> pu_pi = connection(pu, 'out1', pi, 'in1')
    >>> pi_cc = connection(pi, 'out1', cc, 'in1')
    >>> cc_pu = connection(cc, 'out1', pu, 'in1')
    >>> nw.add_conns(pu_pi, pi_cc, cc_pu)
    >>> pi_cc.set_attr(p=1, T=20, fluid={'water': 1})
    >>> pu_pi.set_attr(p=10)
    >>> pu.set_attr(eta_s=0.8, P=1000)
    >>> nw.set_printoptions(print_level='none')
    >>> nw.solve('design')
    >>> round(pi.Q.val, 1) == -round(pu.P.val, 1)
    True
    """

    def component(self):
        return 'cycle closer'

    def inlets(self):
        return ['in1']

    def outlets(self):
        return ['out1']

    def comp_init(self, nw):

        component.comp_init(self, nw)

        # all derivatives are constants
        self.mat_deriv = np.zeros((2, 2, 3 + self.num_fl))
        # derivatives for pressure
        self.mat_deriv[0, 0, 1] = 1
        self.mat_deriv[0, 1, 1] = -1
        # derivatives for enthalpy
        self.mat_deriv[1, 0, 2] = 1
        self.mat_deriv[1, 1, 2] = -1

    def equations(self):
        r"""
        Calculates vector vec_res with results of equations for this component.

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # equation for pressure
        vec_res += [self.inl[0].p.val_SI - self.outl[0].p.val_SI]

        ######################################################################
        # equation for enthalpy
        vec_res += [self.inl[0].h.val_SI - self.outl[0].h.val_SI]

        ######################################################################

        return vec_res

    def derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given equations.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        ######################################################################
        # derivatives with constant value (all for this component)

        return self.mat_deriv
