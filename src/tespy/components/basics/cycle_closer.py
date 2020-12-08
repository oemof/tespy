# -*- coding: utf-8

"""Module for class CycleCloser


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/basics/cycle_closer.py

SPDX-License-Identifier: MIT
"""

import numpy as np

from tespy.components.component import Component
from tespy.tools.data_containers import ComponentProperties as dc_cp

# %%


class CycleCloser(Component):
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

    >>> from tespy.components import CycleCloser, Pipe, Pump
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> nw = Network(['water'], p_unit='bar', T_unit='C', iterinfo=False)
    >>> pi = Pipe('pipe')
    >>> pu = Pump('pump')
    >>> cc = CycleCloser('cycle closing component')
    >>> cc.component()
    'cycle closer'
    >>> pu_pi = Connection(pu, 'out1', pi, 'in1')
    >>> pi_cc = Connection(pi, 'out1', cc, 'in1')
    >>> cc_pu = Connection(cc, 'out1', pu, 'in1')
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

        Component.comp_init(self, nw)

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

    def propagate_fluid_to_target(self, inconn, start):
        r"""
        Fluid propagation to target stops here.

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
        Fluid propagation to source stops here.

        Parameters
        ----------
        outconn : tespy.connections.connection.Connection
            Connection to initialise.

        start : tespy.components.component.Component
            This component is the fluid propagation starting point.
            The starting component is saved to prevent infinite looping.
        """
        return

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
