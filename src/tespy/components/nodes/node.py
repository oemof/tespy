# -*- coding: utf-8

"""Module of class Node.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/nodes/node.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import component_registry
from tespy.components.nodes.merge import Merge
from tespy.components.nodes.splitter import Splitter
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import SimpleDataContainer as dc_simple


@component_registry
class Node(Splitter, Merge):
    r"""
    Class for combined merge and splitting points with multiple inflows and
    outflows.

    .. image:: /api/_images/components/Node.svg
       :alt: flowsheet of the node
       :align: center
       :class: only-light

    .. image:: /api/_images/components/Node_darkmode.svg
       :alt: flowsheet of the node
       :align: center
       :class: only-dark

    Ports
    -----

    - Fluid inlets: in1, in2, ... (variable, count set by :code:`num_in`)
    - Fluid outlets: out1, out2, ... (variable, count set by :code:`num_out`)

    Mandatory Equations
    -------------------

    - mass balance constraint: :py:meth:`mass_flow_func <tespy.components.nodes.base.NodeBase.mass_flow_func>`
    - pressure equality constraints: :py:meth:`pressure_structure_matrix <tespy.components.nodes.base.NodeBase.pressure_structure_matrix>`
    - equal enthalpy at all outlets constraint(s): :py:meth:`enthalpy_structure_matrix <tespy.components.nodes.node.Node.enthalpy_structure_matrix>`
    - equal fluid at all outlets constraint(s): :py:meth:`fluid_structure_matrix <tespy.components.nodes.node.Node.fluid_structure_matrix>`
    - fluid mass fraction constraints: :py:meth:`fluid_func <tespy.components.nodes.merge.Merge.fluid_func>`
    - energy balance constraint: :py:meth:`energy_balance_func <tespy.components.nodes.merge.Merge.energy_balance_func>`

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

    num_in : int
        Number of inlets.

    num_out : int
        Number of outlets.

    offdesign : list
        List containing offdesign parameters (stated as String).

    printout : bool
        Include this component in the network's results printout.

    Example
    -------
    The :code:`Node` serves as a splitter and a merger at the same time. For
    example, you can use it to represent a preheating vessel, which produces
    saturated liquid water by mixing pumped condensate and a steam stream.
    The saturated liquid is then split into a specified amount of outlet
    streams. In this example it will be just two.

    >>> from tespy.components import Source, Sink, Node
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "pressure_difference": "bar",
    ...     "temperature": "degC"
    ... })
    >>> so1 = Source("source1")
    >>> so2 = Source("source2")
    >>> si1 = Sink("sink1")
    >>> si2 = Sink("sink2")
    >>> node = Node("node")
    >>> node.set_attr(num_in=2, num_out=2)
    >>> c1 = Connection(so1, "out1", node, "in1", label="pumped condensate")
    >>> c2 = Connection(so2, "out1", node, "in2", label="extraction steam")
    >>> c3 = Connection(node, "out1", si1, "in1", label="outlet 1")
    >>> c4 = Connection(node, "out2", si2, "in1", label="outlet 2")
    >>> nw.add_conns(c1, c2, c3, c4)

    We can parametrize the system, for example, to preheat 50 kg/s of water.
    The system will then identify, what amount of extraction steam is required
    to preheat the water to the saturated liquid state. Apart from the inlet
    states we have to add one mass flow specification to define the split
    ratio between the two outlets.

    .. note::

        The enthalpy of the fluid will be identical at all exits. If you want
        to have separation of phases, you have to use a `DropletSeparator`
        downstream of this component.

    >>> c1.set_attr(fluid={"water": 1}, m=50, p=3, T=50)
    >>> c2.set_attr(fluid={"water": 1}, T=200)
    >>> c3.set_attr(x=0)
    >>> c4.set_attr(m=1)
    >>> nw.solve("design")
    """

    @staticmethod
    def get_parameters():
        return {
            'num_out': dc_simple(dtype="int", description="number of outlets"),
            'num_in': dc_simple(dtype="int", description="number of inlets")
        }

    @classmethod
    def port_schema(cls):
        return {
            "inlets": {"type": "variable", "parameter": "num_in", "pattern": "in{n}", "min": 2},
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
            'pressure_constraints': dc_cmc(**{
                'structure_matrix': self.pressure_structure_matrix,
                'num_eq_sets': self.num_i + self.num_o - 1,
                'description': 'pressure equality constraints'
            }),
            'outlet_enthalpy_constraints': dc_cmc(**{
                'structure_matrix': self.enthalpy_structure_matrix,
                'num_eq_sets': self.num_o - 1,
                'description': 'equal enthalpy at all outlets constraint(s)'
            }),
            'outlet_fluid_constraints': dc_cmc(**{
                'structure_matrix': self.fluid_structure_matrix,
                'num_eq_sets': self.num_o - 1,
                'description': 'equal fluid at all outlets constraint(s)'
            }),
            'fluid_constraints': dc_cmc(**{
                'num_eq_sets': 1,
                'func': self.fluid_func,
                'dependents': self.fluid_dependents,
                'description': 'fluid mass fraction constraints'
            }),
            'energy_balance_constraints': dc_cmc(**{
                'num_eq_sets': 1,
                'func': self.energy_balance_func,
                'dependents': self.energy_balance_dependents,
                'description': 'energy balance constraint'
            })
        }

    def inlets(self):
        if self.num_in.is_set:
            return [f'in{i + 1}' for i in range(self.num_in.val)]
        else:
            self.set_attr(num_in=2)
            return self.inlets()

    def outlets(self):
        if self.num_out.is_set:
            return [f'out{i + 1}' for i in range(self.num_out.val)]
        else:
            self.set_attr(num_out=2)
            return self.outlets()

    def enthalpy_structure_matrix(self, k):
        r"""
        Set up the structure matrix for the outlet enthalpy constraints:

        All outlet enthalpy values must be identical.
        """
        for eq, conn in enumerate(self.outl[1:]):
            self._structure_matrix[k + eq, self.outl[0].h.sm_col] = 1
            self._structure_matrix[k + eq, conn.h.sm_col] = -1

    def fluid_structure_matrix(self, k):
        r"""
        Set up the structure matrix for the outlet fluid constraints:

        All outlet fluid compositions must be identical.
        """
        for eq, conn in enumerate(self.outl[1:]):
            self._structure_matrix[k + eq, self.outl[0].fluid.sm_col] = 1
            self._structure_matrix[k + eq, conn.fluid.sm_col] = -1
