# -*- coding: utf-8

"""Module for class CycleCloser


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/basics/cycle_closer.py

SPDX-License-Identifier: MIT
"""

import numpy as np

from tespy.components.component import Component
from tespy.components.component import component_registry
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import ComponentProperties as dc_cp


@component_registry
class CycleCloser(Component):
    r"""
    Component for closing cycles.

    **Mandatory Equations**

    - pressure: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`
    - enthalpy: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`

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
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "temperature": "degC"
    ... })
    >>> pi = Pipe('pipe')
    >>> pu = Pump('pump')
    >>> cc = CycleCloser('cycle closing component')
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
    def get_parameters():
        return {
            'mass_deviation': dc_cp(
                _val=0, max_val=1e-3, is_result=True, quantity="mass_flow"
            ),
            'fluid_deviation': dc_cp(_val=0, max_val=1e-5, is_result=True)
        }

    def get_mandatory_constraints(self):
        return {
            'pressure_equality_constraint': dc_cmc(**{
                'num_eq_sets': 1,
                'structure_matrix': self.variable_equality_structure_matrix,
                'func_params': {'variable': 'p'}
            }),
            'enthalpy_equality_constraint': dc_cmc(**{
                'num_eq_sets': 1,
                'structure_matrix': self.variable_equality_structure_matrix,
                'func_params': {'variable': 'h'}
            })
        }

    @staticmethod
    def inlets():
        return ['in1']

    @staticmethod
    def outlets():
        return ['out1']

    def start_fluid_wrapper_branch(self):
        outconn = self.outl[0]
        branch = {
            "connections": [outconn],
            "components": [self]
        }
        outconn.target.propagate_wrapper_to_target(branch)

        return {outconn.label: branch}

    def propagate_wrapper_to_target(self, branch):
        branch["components"] += [self]
        return

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        # calculate deviation in mass flow
        self.mass_deviation.val_SI = abs(
            self.inl[0].m.val_SI - self.outl[0].m.val_SI
        )

        # calculate deviation in fluid composition
        d1 = self.inl[0].fluid.val
        d2 = self.outl[0].fluid.val
        diff = [d1[key] - d2[key] for key in d1.keys()]
        self.fluid_deviation.val_SI = np.linalg.norm(diff)
