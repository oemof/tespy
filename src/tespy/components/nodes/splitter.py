# -*- coding: utf-8

"""Module of class Splitter.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/nodes/splitter.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import component_registry
from tespy.components.nodes.base import NodeBase
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import SimpleDataContainer as dc_simple


@component_registry
class Splitter(NodeBase):
    r"""
    Split up a mass flow in parts of identical enthalpy and fluid composition.

    .. image:: /api/_images/components/Splitter.svg
       :alt: flowsheet of the splitter
       :align: center
       :class: only-light

    .. image:: /api/_images/components/Splitter_darkmode.svg
       :alt: flowsheet of the splitter
       :align: center
       :class: only-dark

    Ports
    -----

    Fluid inlets: in1

    Fluid outlets: out1, out2, ... (variable, count set by :code:`num_out`)

    Mandatory Equations
    -------------------

    - mass balance constraint: :py:meth:`mass_flow_func <tespy.components.nodes.base.NodeBase.mass_flow_func>`
    - equal enthalpy at all outlets constraint: :py:meth:`enthalpy_structure_matrix <tespy.components.nodes.splitter.Splitter.enthalpy_structure_matrix>`
    - pressure equality constraints: :py:meth:`pressure_structure_matrix <tespy.components.nodes.base.NodeBase.pressure_structure_matrix>`
    - fluid equality constraints: :py:meth:`fluid_structure_matrix <tespy.components.nodes.splitter.Splitter.fluid_structure_matrix>`

    Parameters
    ----------

    char_warnings : bool
        Ignore warnings on default characteristics usage for this component.

    design : list
        List containing design parameters (stated as String).

    design_path : str
        Path to the components design case.

    label : str
        The label of the component.

    local_design : bool
        Treat this component in design mode in an offdesign calculation.

    local_offdesign : bool
        Treat this component in offdesign mode in a design calculation.

    num_out : int
        Number of outlets.

    offdesign : list
        List containing offdesign parameters (stated as String).

    printout : bool
        Include this component in the network's results printout.

    Example
    -------
    A splitter is used to split up a single mass flow into a specified number
    of different parts at identical pressure, enthalpy and fluid composition.

    >>> from tespy.components import Sink, Source, Splitter
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "pressure_difference": "bar",
    ... "temperature": "degC"
    ... })
    >>> so = Source('source')
    >>> si1 = Sink('sink1')
    >>> si2 = Sink('sink2')
    >>> si3 = Sink('sink3')
    >>> s = Splitter('splitter', num_out=3)
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
    def get_parameters():
        return {'num_out': dc_simple(dtype="int", description="number of outlets")}

    @classmethod
    def port_schema(cls):
        return {
            "inlets": {"type": "fixed", "ports": ["in1"]},
            "outlets": {"type": "variable", "parameter": "num_out", "pattern": "out{n}", "min": 2},
            "powerinlets": {"type": "fixed", "ports": []},
            "poweroutlets": {"type": "fixed", "ports": []},
            "heatinlets": {"type": "fixed", "ports": []},
            "heatoutlets": {"type": "fixed", "ports": []},
        }

    def get_mandatory_constraints(self):
        return {
            'mass_flow_constraints': dc_cmc(**{
                'func': self.mass_flow_func,
                'dependents': self.mass_flow_dependents,
                'num_eq_sets': 1,
                'description': 'mass balance constraint'
            }),
            'energy_balance_constraints': dc_cmc(**{
                'structure_matrix': self.enthalpy_structure_matrix,
                'num_eq_sets': self.num_o,
                'description': 'equal enthalpy at all outlets constraint'
            }),
            'pressure_constraints': dc_cmc(**{
                'structure_matrix': self.pressure_structure_matrix,
                'num_eq_sets': self.num_o,
                'description': 'pressure equality constraints'
            }),
            'fluid_constraints': dc_cmc(**{
                'structure_matrix': self.fluid_structure_matrix,
                'num_eq_sets': self.num_o,
                'description': 'fluid equality constraints'
            })
        }

    @staticmethod
    def inlets():
        return ['in1']

    def outlets(self):
        if self.num_out.is_set:
            return [f'out{i + 1}' for i in range(self.num_out.val)]
        else:
            self.set_attr(num_out=2)
            return self.outlets()

    def propagate_wrapper_to_target(self, branch):
        branch["components"] += [self]
        for outconn in self.outl:
            branch["connections"] += [outconn]
            outconn.target.propagate_wrapper_to_target(branch)

    def enthalpy_structure_matrix(self, k):
        r"""
        Calculate partial derivatives for energy balance equation.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        for eq, conn in enumerate(self.outl):
            self._structure_matrix[k + eq, self.inl[0].h.sm_col] = 1
            self._structure_matrix[k + eq, conn.h.sm_col] = -1

    def fluid_structure_matrix(self, k):
        r"""
        Calculate partial derivatives for all pressure equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the fluid equations.
        """
        for eq, conn in enumerate(self.outl):
            self._structure_matrix[k + eq, self.inl[0].fluid.sm_col] = 1
            self._structure_matrix[k + eq, conn.fluid.sm_col] = -1
