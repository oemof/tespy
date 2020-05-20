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
from tespy.tools.data_containers import dc_cp
from tespy.tools.data_containers import dc_simple

# %%


class cycle_closer(component):
    r"""
    Component for closing cycles.

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

    >>> from tespy.components import cycle_closer, pipe, pump
    >>> from tespy.connections import connection
    >>> from tespy.networks import network
    >>> nw = network(['water'], p_unit='bar', T_unit='C', iterinfo=False)
    >>> pi = pipe('pipe')
    >>> pu = pump('pump')
    >>> cc = cycle_closer('cycle closing component')
    >>> cc.component()
    'cycle closer'
    >>> pu_pi = connection(pu, 'out1', pi, 'in1')
    >>> pi_cc = connection(pi, 'out1', cc, 'in1')
    >>> cc_pu = connection(cc, 'out1', pu, 'in1')
    >>> nw.add_conns(pu_pi, pi_cc, cc_pu)
    >>> pi_cc.set_attr(p=1, T=20, fluid={'water': 1})
    >>> pu_pi.set_attr(p=10)
    >>> pu.set_attr(eta_s=0.8, P=1000)
    >>> nw.solve('design')
    >>> round(pi.Q.val, 1) == -round(pu.P.val, 1)
    True
    """

    @staticmethod
    def component():
        return 'cycle closer'

    @staticmethod
    def attr():
        return {
            'mass_deviation': dc_cp(val=0, max_val=1e-3, printout=False),
            'fluid_deviation': dc_cp(val=0, max_val=1e-5, printout=False)
        }

    @staticmethod
    def inlets():
        return ['in1']

    @staticmethod
    def outlets():
        return ['out1']

    def comp_init(self, nw):

        component.comp_init(self, nw)

        # number of mandatroy equations for
        # pressure: 1
        # enthalpy: 1
        self.num_eq = 2

        self.jacobian = np.zeros((
            self.num_eq,
            self.num_i + self.num_o + self.num_vars,
            self.num_nw_vars))

        self.residual = np.ones(self.num_eq)
        # derivatives for pressure
        self.jacobian[0, 0, 1] = 1
        self.jacobian[0, 1, 1] = -1
        # derivatives for enthalpy
        self.jacobian[1, 0, 2] = 1
        self.jacobian[1, 1, 2] = -1

    def equations(self):
        r"""Calculate residual vector with results of equations."""
        k = 0
        ######################################################################
        # equation for pressure
        self.residual[k] = self.inl[0].p.val_SI - self.outl[0].p.val_SI
        k += 1

        ######################################################################
        # equation for enthalpy
        self.residual[k] = self.inl[0].h.val_SI - self.outl[0].h.val_SI
        k += 1

    def derivatives(self, vek_z):
        r"""Calculate partial derivatives for given equations."""
        ######################################################################
        # all derivatives are static

    def calc_parameters(self):

        # calculate deviation in mass flow
        self.mass_deviation.val = np.abs(
            self.inl[0].m.val_SI - self.outl[0].m.val_SI)

        # calculate deviation in fluid composition
        d1 = self.inl[0].fluid.val
        d2 = self.outl[0].fluid.val
        diff = [d1[key] - d2[key] for key in d1.keys()]
        self.fluid_deviation.val = np.linalg.norm(diff)

        self.check_parameter_bounds()

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

    Example
    -------
    Create a sink and specify a label.

    >>> from tespy.components import sink
    >>> si = sink('a labeled sink')
    >>> si.component()
    'sink'
    >>> si.label
    'a labeled sink'
    """

    @staticmethod
    def component():
        return 'sink'

    @staticmethod
    def inlets():
        return ['in1']

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

    Example
    -------
    Create a source and specify a label.

    >>> from tespy.components import source
    >>> so = source('a labeled source')
    >>> so.component()
    'source'
    >>> so.label
    'a labeled source'
    """

    @staticmethod
    def component():
        return 'source'

    @staticmethod
    def outlets():
        return ['out1']

# %%


class subsystem_interface(component):
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

    num_inter : float/tespy.tools.data_containers.dc_simple
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

    >>> from tespy.components import sink, source, subsystem_interface
    >>> from tespy.connections import connection
    >>> from tespy.networks import network
    >>> fluids = ['H2O', 'N2']
    >>> nw = network(fluids=fluids)
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg', iterinfo=False)
    >>> so = source('source 1')
    >>> si = sink('sink 1')
    >>> IF = subsystem_interface('subsystem interface', num_inter=1)
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

        component.comp_init(self, nw)

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
            mass flow: 0, pressure: 1, enthalpy: 2.

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
