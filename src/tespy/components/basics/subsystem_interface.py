# -*- coding: utf-8

"""Module for class SubsystemInterface.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/basics/subsystem_interface.py

SPDX-License-Identifier: MIT
"""

import numpy as np

from tespy.components.component import Component
from tespy.tools.data_containers import DataContainerSimple as dc_simple


class SubsystemInterface(Component):
    r"""
    The subsystem interface does not change fluid properties.

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

        .. image:: _images/SubsystemInterface.svg
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

    num_inter : float, tespy.tools.data_containers.DataContainerSimple
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
    :ref:`this section <tespy_subsystems_label>` for more information on
    building your own subsystems.

    >>> from tespy.components import Sink, Source, SubsystemInterface
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> fluids = ['H2O', 'N2']
    >>> nw = Network(fluids=fluids)
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg', iterinfo=False)
    >>> so = Source('source 1')
    >>> si = Sink('sink 1')
    >>> IF = SubsystemInterface('subsystem interface', num_inter=1)
    >>> IF.component()
    'subsystem interface'
    >>> len(IF.inlets())
    1

    Add the subsystem the interface to a minimal network containing only a
    source and a sink. The interface does not change the fluid properties
    in any way.

    >>> inc = Connection(so, 'out1', IF, 'in1')
    >>> outg = Connection(IF, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)
    >>> inc.set_attr(fluid={'H2O': 1, 'N2': 0}, T=40, p=3, m=100)
    >>> nw.solve('design')
    >>> inc.m.val_SI == outg.m.val_SI
    True
    >>> inc.h.val_SI == outg.h.val_SI
    True
    """

    @staticmethod
    def component():
        return 'subsystem interface'

    @staticmethod
    def attr():
        return {'num_inter': dc_simple()}

    def inlets(self):
        if self.num_inter.is_set:
            return ['in' + str(i + 1) for i in range(self.num_inter.val)]
        else:
            return ['in1']

    def outlets(self):
        if self.num_inter.is_set:
            return ['out' + str(i + 1) for i in range(self.num_inter.val)]
        else:
            return ['out1']

    def comp_init(self, nw):

        Component.comp_init(self, nw)

        # number of mandatroy equations for
        # fluid: num_inter * num_nw_fluids
        # mass flow: num_inter
        # pressure: num_inter
        # enthalpy: num_inter
        self.num_eq = (self.num_nw_fluids + 3) * self.num_i

        self.jacobian = np.zeros((
            self.num_eq,
            2 * self.num_i,
            self.num_nw_vars))

        self.residual = np.ones(self.num_eq)
        stop = self.num_nw_fluids * self.num_i
        self.jacobian[0:stop] = self.fluid_deriv()
        start = stop
        stop = start + self.num_i
        self.jacobian[start:stop] = self.inout_deriv(0)
        start = stop
        stop = start + self.num_i
        self.jacobian[start:stop] = self.inout_deriv(1)
        start = stop
        stop = start + self.num_i
        self.jacobian[start:stop] = self.inout_deriv(2)

    def equations(self):
        r"""Calculate residual vector with results of equations."""
        k = 0
        ######################################################################
        # eqations for fluids
        for i in range(self.num_i):
            for fluid, x in self.inl[i].fluid.val.items():
                self.residual[k] = x - self.outl[i].fluid.val[fluid]
                k += 1

        ######################################################################
        # equations for mass flow
        for i in range(self.num_i):
            self.residual[k] = self.inl[i].m.val_SI - self.outl[i].m.val_SI
            k += 1

        ######################################################################
        # equations for pressure
        for i in range(self.num_i):
            self.residual[k] = self.inl[i].p.val_SI - self.outl[i].p.val_SI
            k += 1

        ######################################################################
        # equations for enthalpy
        for i in range(self.num_i):
            self.residual[k] = self.inl[i].h.val_SI - self.outl[i].h.val_SI
            k += 1

        ######################################################################

    def derivatives(self, vek_z):
        r"""Calculate partial derivatives for given equations."""
        ######################################################################
        # all derivatives are static

    def inout_deriv(self, pos):
        r"""
        Calculate partial derivatives.

        Method applies for all mass flow, pressure and enthalpy equations.

        Parameters
        ----------
        pos : int
            Position of the variable in the matrix of derivatives.

            - mass flow: 0
            - pressure: 1
            - enthalpy: 2

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((self.num_i, 2 * self.num_i, self.num_nw_vars))
        for i in range(self.num_i):
            deriv[i, i, pos] = 1
        for j in range(self.num_i):
            deriv[j, j + self.num_i, pos] = -1
        return deriv
