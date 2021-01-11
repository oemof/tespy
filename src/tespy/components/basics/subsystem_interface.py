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

    **Mandatory Equations**

    - :py:meth:`tespy.components.component.Component.fluid_func`
    - :py:meth:`tespy.components.component.Component.mass_flow_func`
    - Pressure:
      :py:meth:`tespy.components.basics.subsystem_interface.SubsystemInterface.variable_equality_func`
    - Enthalpy:
      :py:meth:`tespy.components.basics.subsystem_interface.SubsystemInterface.variable_equality_func`

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

        # number of mandatroy equations for
        # fluid: num_inter * num_nw_fluids
        # mass flow: num_inter
        # pressure: num_inter
        # enthalpy: num_inter
        Component.comp_init(self, nw, num_eq=(len(nw.fluids) + 3) * self.num_i)
        # all derivatives are constant
        pos = self.num_nw_fluids * self.num_i
        self.jacobian[0:pos] = self.fluid_deriv()
        self.jacobian[pos:pos + self.num_i] = self.mass_flow_deriv()
        pos += self.num_i
        self.jacobian[pos:pos + self.num_i] = self.variable_equality_deriv(1)
        pos += self.num_i
        self.jacobian[pos:pos + self.num_i] = self.variable_equality_deriv(2)

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
        # equations for fluids
        num_eq = self.num_nw_fluids * self.num_i
        self.residual[k:k + num_eq] = self.fluid_func()
        if doc:
            self.equation_docs[k:k + num_eq] = self.fluid_func(doc=doc)
        k += num_eq

        ######################################################################
        # equations for mass flow
        self.residual[k:k + self.num_i] = self.mass_flow_func()
        if doc:
            self.equation_docs[k:k + self.num_i] = self.mass_flow_func(doc=doc)
        k += self.num_i

        ######################################################################
        # equations for pressure and enthalpy:
        for param in ['p', 'h']:
            self.residual[k:k + self.num_i] = (
                self.variable_equality_func(param))
            if doc:
                self.equation_docs[k:k + self.num_i] = (
                    self.variable_equality_func(param, doc=doc))
            k += self.num_i

        return k

    def variable_equality_func(self, param, doc=False):
        r"""
        Calculate the residual value for primary variable equality equation.

        This equation defines equatilty of mass flow, pressure or enthalpy
        at an inlet and its corresponding outlet.

        Returns
        -------
        residual : list
            Vector with residual values.
        """
        if not doc:
            residual = []
            for i in range(self.num_i):
                residual += [
                    self.inl[i].get_attr(param).val_SI -
                    self.outl[i].get_attr(param).val_SI]
            return residual
        else:
            indices = list(range(1, self.num_i + 1))
            if len(indices) > 1:
                indices = ', '.join(str(idx) for idx in indices)
            else:
                indices = str(indices[0])
            latex = (
                r'0=' + param + r'_{\mathrm{in,}i}-' + param +
                r'_{\mathrm{out,}i}\; \forall i \in [' + indices + r']')
            return (
                [self.generate_latex(latex, param + '_equality_func')] +
                (self.num_i - 1) * [''])

    def variable_equality_deriv(self, pos):
        r"""
        Calculate partial derivatives for pressure and enthalpy equations.

        Parameters
        ----------
        pos : int
            Position of the variable in the matrix of derivatives.

            - pressure: 1
            - enthalpy: 2

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((self.num_i, 2 * self.num_i, self.num_nw_vars))
        for i in range(self.num_i):
            deriv[i, i, pos] = 1
        for j in range(self.num_i):
            deriv[j, j + self.num_i, pos] = -1
        return deriv
