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

    Ports
    -----

    - Fluid inlets: in1
    - Fluid outlets: out1

    Mandatory Equations
    -------------------

    - pressure equality constraint: :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
    - enthalpy equality constraint: :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`

    Parameters
    ----------

    char_warnings : bool
        Ignore warnings on default characteristics usage for this component.

    design : list
        List containing design parameters (stated as String).

    design_path : str
        Path to the components design case.

    fluid_deviation : float, dict
        Norm of absolute deviation of fluid composition between inlet and
        outlet.

    label : str
        The label of the component.

    local_design : bool
        Treat this component in design mode in an offdesign calculation.

    local_offdesign : bool
        Treat this component in offdesign mode in a design calculation.

    mass_deviation : float, dict
        Absolute deviation of mass flow between inlet and outlet. Quantity:
        :code:`mass_flow`.

    offdesign : list
        List containing offdesign parameters (stated as String).

    printout : bool
        Include this component in the network's results printout.

    Notes
    -----

    .. note::

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
    ...     "pressure": "bar", "pressure_difference": "bar",
    ...     "temperature": "degC"
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

    def _calc_mass_deviation(self):
        return abs(self.inl[0].m.val_SI - self.outl[0].m.val_SI)

    def _calc_fluid_deviation(self):
        d1 = self.inl[0].fluid.val
        d2 = self.outl[0].fluid.val
        return np.linalg.norm([d1[k] - d2[k] for k in d1])

    def get_parameters(self):
        return {
            'mass_deviation': dc_cp(
                _val=0, max_val=1e-3, is_result=True, quantity="mass_flow",
                description="absolute deviation of mass flow between inlet and outlet",
                calc=self._calc_mass_deviation
            ),
            'fluid_deviation': dc_cp(
                _val=0, max_val=1e-5, is_result=True,
                description="norm of absolute deviation of fluid composition between inlet and outlet",
                calc=self._calc_fluid_deviation
            )
        }

    def get_mandatory_constraints(self):
        return {
            "pressure_equality_constraint": dc_cmc(**{
                "num_eq_sets": 1,
                "structure_matrix": self.variable_equality_structure_matrix,
                "func_params": {"variable": "p"},
                "description": "pressure equality constraint"
            }),
            "enthalpy_equality_constraint": dc_cmc(**{
                "num_eq_sets": 1,
                "structure_matrix": self.variable_equality_structure_matrix,
                "func_params": {"variable": "h"},
                "description": "enthalpy equality constraint"
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

