# -*- coding: utf-8
"""Module for tespy network class.

The network is the container for every TESPy simulation. The network class
automatically creates the system of equations describing topology and
parametrization of a specific model and solves it.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/networks/networks.py

SPDX-License-Identifier: MIT
"""
import json
import math
import os
from time import time

import numpy as np
import pandas as pd
from numpy.linalg import norm
from tabulate import tabulate

from tespy import connections as con
from tespy.tools import fluid_properties as fp
from tespy.tools import helpers as hlp
from tespy.tools import logger
from tespy.tools.data_containers import ComponentCharacteristicMaps as dc_cm
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import FluidComposition as dc_flu
from tespy.tools.data_containers import FluidProperties as dc_fp
from tespy.tools.data_containers import GroupedComponentCharacteristics as dc_gcc
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.global_vars import ERR
from tespy.tools.global_vars import fluid_property_data as fpd

# Only require cupy if Cuda shall be used
try:
    import cupy as cu
except ModuleNotFoundError:
    cu = None


class Network:
    r"""
    Class component is the base class of all TESPy components.

    Parameters
    ----------
    h_range : list
        List with minimum and maximum values for enthalpy value range.

    h_unit : str
        Specify the unit for enthalpy: 'J / kg', 'kJ / kg', 'MJ / kg'.

    iterinfo : boolean
        Print convergence progress to console.

    m_range : list
        List with minimum and maximum values for mass flow value range.

    m_unit : str
        Specify the unit for mass flow: 'kg / s', 't / h'.

    p_range : list
        List with minimum and maximum values for pressure value range.

    p_unit : str
        Specify the unit for pressure: 'Pa', 'psi', 'bar', 'MPa'.

    s_unit : str
        Specify the unit for specific entropy: 'J / kgK', 'kJ / kgK',
        'MJ / kgK'.

    T_unit : str
        Specify the unit for temperature: 'K', 'C', 'F', 'R'.

    v_unit : str
        Specify the unit for volumetric flow: 'm3 / s', 'm3 / h', 'l / s',
        'l / h'.

    vol_unit : str
        Specify the unit for specific volume: 'm3 / kg', 'l / kg'.

    x_unit : str
        Specify the unit for steam mass fraction: '-', '%'.

    Note
    ----
    Unit specification is optional: If not specified the SI unit (first
    element in above lists) will be applied!

    Range specification is optional, too. The value range is used to stabilize
    the newton algorithm. For more information see the "getting started"
    section in the online-documentation.

    Example
    -------
    Basic example for a setting up a tespy.networks.network.Network object.
    Specifying the fluids is mandatory! Unit systems, fluid property range and
    iterinfo are optional.

    Standard value for iterinfo is :code:`True`. This will print out
    convergence progress to the console. You can stop the printouts by setting
    this property to :code:`False`.

    >>> from tespy.networks import Network
    >>> mynetwork = Network(p_unit='bar', T_unit='C')
    >>> mynetwork.set_attr(p_range=[1, 10])
    >>> type(mynetwork)
    <class 'tespy.networks.network.Network'>
    >>> mynetwork.set_attr(iterinfo=False)
    >>> mynetwork.iterinfo
    False
    >>> mynetwork.set_attr(iterinfo=True)
    >>> mynetwork.iterinfo
    True

    A simple network consisting of a source, a pipe and a sink. This example
    shows how the printout parameter can be used. We specify
    :code:`printout=False` for both connections, the pipe as well as the heat
    bus. Therefore the :code:`.print_results()` method should not print any
    results.

    >>> from tespy.networks import Network
    >>> from tespy.components import Source, Sink, Pipe
    >>> from tespy.connections import Connection, Bus
    >>> nw = Network(T_unit='C', p_unit='bar', v_unit='m3 / s')
    >>> so = Source('source')
    >>> si = Sink('sink')
    >>> p = Pipe('pipe', Q=0, pr=0.95, printout=False)
    >>> a = Connection(so, 'out1', p, 'in1')
    >>> b = Connection(p, 'out1', si, 'in1')
    >>> nw.add_conns(a, b)
    >>> a.set_attr(fluid={'CH4': 1}, T=30, p=10, m=10, printout=False)
    >>> b.set_attr(printout=False)
    >>> b = Bus('heat bus')
    >>> b.add_comps({'comp': p})
    >>> nw.add_busses(b)
    >>> b.set_attr(printout=False)
    >>> nw.set_attr(iterinfo=False)
    >>> nw.solve('design')
    >>> nw.print_results()
    """

    def __init__(self, **kwargs):
        self.set_defaults()
        self.set_attr(**kwargs)

    def _serialize(self):
        return {
            "m_unit": self.m_unit,
            "m_range": list(self.m_range),
            "p_unit": self.p_unit,
            "p_range": list(self.p_range),
            "h_unit": self.h_unit,
            "h_range": list(self.h_range),
            "T_unit": self.T_unit,
            "x_unit": self.x_unit,
            "v_unit": self.v_unit,
            "s_unit": self.s_unit,
        }

    def set_defaults(self):
        """Set default network properties."""
        # connection dataframe

        dtypes={
            "object": object,
            "source": object,
            "source_id": str,
            "target": object,
            "target_id": str
        }
        self.conns = pd.DataFrame(
            columns=list(dtypes.keys())
        ).astype(dtypes)
        self.all_fluids = set()
        # component dataframe
        dtypes = {
            "comp_type": str,
            "object": object,
        }
        self.comps = pd.DataFrame(
            columns=list(dtypes.keys())
        ).astype(dtypes)
        # user defined function dictionary for fast access
        self.user_defined_eq = {}
        # bus dictionary
        self.busses = {}
        # results and specification dictionary
        self.results = {}
        self.specifications = {}

        self.specifications['lookup'] = {
            'properties': 'prop_specifications',
            'chars': 'char_specifications',
            'variables': 'var_specifications',
            'groups': 'group_specifications'
        }

        # in case of a design calculation after an offdesign calculation
        self.redesign = False

        self.checked = False
        self.design_path = None
        self.iterinfo = True

        msg = 'Default unit specifications:\n'
        for prop, data in fpd.items():
            # standard unit set
            self.__dict__.update({prop + '_unit': data['SI_unit']})
            msg += data['text'] + ': ' + data['SI_unit'] + '\n'

        # don't need the last newline
        logger.debug(msg[:-1])

        # generic value range
        self.m_range_SI = [-1e12, 1e12]
        self.p_range_SI = [2e2, 300e5]
        self.h_range_SI = [1e3, 7e6]

        for prop in ['m', 'p', 'h']:
            limits = self.get_attr(prop + '_range_SI')
            msg = (
                f"Default {fpd[prop]['text']} limits\n"
                f"min: {limits[0]} {self.get_attr(prop + '_unit')}\n"
                f"max: {limits[1]} {self.get_attr(prop + '_unit')}"
            )
            logger.debug(msg)

    def set_attr(self, **kwargs):
        r"""
        Set, resets or unsets attributes of a network.

        Parameters
        ----------
        h_range : list
            List with minimum and maximum values for enthalpy value range.

        h_unit : str
            Specify the unit for enthalpy: 'J / kg', 'kJ / kg', 'MJ / kg'.

        iterinfo : boolean
            Print convergence progress to console.

        m_range : list
            List with minimum and maximum values for mass flow value range.

        m_unit : str
            Specify the unit for mass flow: 'kg / s', 't / h'.

        p_range : list
            List with minimum and maximum values for pressure value range.

        p_unit : str
            Specify the unit for pressure: 'Pa', 'psi', 'bar', 'MPa'.

        s_unit : str
            Specify the unit for specific entropy: 'J / kgK', 'kJ / kgK',
            'MJ / kgK'.

        T_unit : str
            Specify the unit for temperature: 'K', 'C', 'F', 'R'.

        v_unit : str
            Specify the unit for volumetric flow: 'm3 / s', 'm3 / h', 'l / s',
            'l / h'.

        vol_unit : str
            Specify the unit for specific volume: 'm3 / kg', 'l / kg'.
        """
        # unit sets
        for prop in fpd.keys():
            unit = f'{prop}_unit'
            if unit in kwargs:
                if kwargs[unit] in fpd[prop]['units']:
                    self.__dict__.update({unit: kwargs[unit]})
                    msg = f'Setting {fpd[prop]["text"]} unit: {kwargs[unit]}.'
                    logger.debug(msg)
                else:
                    keys = ', '.join(fpd[prop]['units'].keys())
                    msg = f'Allowed units for {fpd[prop]["text"]} are: {keys}'
                    logger.error(msg)
                    raise ValueError(msg)

        for prop in ['m', 'p', 'h']:
            if f'{prop}_range' in kwargs:
                if isinstance(kwargs[f'{prop}_range'], list):
                    self.__dict__.update({
                        f'{prop}_range_SI': [hlp.convert_to_SI(
                            prop, value,
                            self.get_attr(f'{prop}_unit')
                        ) for value in kwargs[f'{prop}_range']]
                    })
                else:
                    msg = f'Specify the range as list: [{prop}_min, {prop}_max]'
                    logger.error(msg)
                    raise TypeError(msg)

                limits = self.get_attr(f'{prop}_range_SI')
                msg = (
                    f'Setting {fpd[prop]["text"]} limits\n'
                    f'min: {limits[0]} {fpd[prop]["SI_unit"]}\n'
                    f'max: {limits[1]} {fpd[prop]["SI_unit"]}'
                )
                logger.debug(msg)

        # update non SI value ranges
        for prop in ['m', 'p', 'h']:
            SI_range = self.get_attr(f'{prop}_range_SI')
            self.__dict__.update({
                f'{prop}_range': [hlp.convert_from_SI(
                    prop, SI_value,
                    self.get_attr(f'{prop}_unit')
                ) for SI_value in SI_range]
            })

        self.iterinfo = kwargs.get('iterinfo', self.iterinfo)

        if not isinstance(self.iterinfo, bool):
            msg = ('Network parameter iterinfo must be True or False!')
            logger.error(msg)
            raise TypeError(msg)

    def get_attr(self, key):
        r"""
        Get the value of a network attribute.

        Parameters
        ----------
        key : str
            The attribute you want to retrieve.

        Returns
        -------
        out :
            Specified attribute.
        """
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            msg = f"Network has no attribute '{key}'."
            logger.error(msg)
            raise KeyError(msg)

    def add_subsys(self, *args):
        r"""
        Add one or more subsystems to the network.

        Parameters
        ----------
        c : tespy.components.subsystem.Subsystem
            The subsystem to be added to the network, subsystem objects si
            :code:`network.add_subsys(s1, s2, s3, ...)`.
        """
        for subsys in args:
            for c in subsys.conns.values():
                self.add_conns(c)

    def get_conn(self, label):
        r"""
        Get Connection via label.

        Parameters
        ----------
        label : str
            Label of the Connection object.

        Returns
        -------
        c : tespy.connections.connection.Connection
            Connection object with specified label, None if no Connection of
            the network has this label.
        """
        try:
            return self.conns.loc[label, 'object']
        except KeyError:
            logger.warning(f"Connection with label {label} not found.")
            return None

    def get_comp(self, label):
        r"""
        Get Component via label.

        Parameters
        ----------
        label : str
            Label of the Component object.

        Returns
        -------
        c : tespy.components.component.Component
            Component object with specified label, None if no Component of
            the network has this label.
        """
        try:
            return self.comps.loc[label, 'object']
        except KeyError:
            logger.warning(f"Component with label {label} not found.")
            return None

    def add_conns(self, *args):
        r"""
        Add one or more connections to the network.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            The connection to be added to the network, connections objects ci
            :code:`add_conns(c1, c2, c3, ...)`.
        """
        for c in args:
            if not isinstance(c, con.Connection):
                msg = (
                    'Must provide tespy.connections.connection.Connection '
                    'objects as parameters.'
                )
                logger.error(msg)
                raise TypeError(msg)

            elif c.label in self.conns.index:
                msg = (
                    'There is already a connection with the label '
                    f'{c.label}. The connection labels must be unique!'
                )
                logger.error(msg)
                raise ValueError(msg)

            c.good_starting_values = False

            self.conns.loc[c.label] = [
                c, c.source, c.source_id, c.target, c.target_id
            ]

            msg = f'Added connection {c.label} to network.'
            logger.debug(msg)
            # set status "checked" to false, if connection is added to network.
            self.checked = False
        self._add_comps(*args)

    def del_conns(self, *args):
        """
        Remove one or more connections from the network.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            The connection to be removed from the network, connections objects
            ci :code:`del_conns(c1, c2, c3, ...)`.
        """
        comps = list({cp for c in args for cp in [c.source, c.target]})
        for c in args:
            self.conns.drop(c.label, inplace=True)
            if "Connection" in self.results:
                self.results["Connection"].drop(
                    c.label, inplace=True, errors="ignore"
                )
            msg = f'Deleted connection {c.label} from network.'
            logger.debug(msg)

        self._del_comps(comps)

        # set status "checked" to false, if connection is deleted from network.
        self.checked = False

    def check_conns(self):
        r"""Check connections for multiple usage of inlets or outlets."""
        dub = self.conns.loc[self.conns.duplicated(["source", "source_id"])]
        for c in dub['object']:
            targets = []
            mask = (
                (self.conns["source"].values == c.source)
                & (self.conns["source_id"].values == c.source_id)
            )
            for conns in self.conns.loc[mask, "object"]:
                targets += [f"\"{conns.target.label}\" ({conns.target_id})"]
            targets = ", ".join(targets)

            msg = (
                f"The source \"{c.source.label}\" ({c.source_id}) is attached "
                f"to more than one component on the target side: {targets}. "
                "Please check your network configuration."
            )
            logger.error(msg)
            raise hlp.TESPyNetworkError(msg)

        dub = self.conns.loc[
            self.conns.duplicated(['target', 'target_id'])
        ]
        for c in dub['object']:
            sources = []
            mask = (
                (self.conns["target"].values == c.target)
                & (self.conns["target_id"].values == c.target_id)
            )
            for conns in self.conns.loc[mask, "object"]:
                sources += [f"\"{conns.source.label}\" ({conns.source_id})"]
            sources = ", ".join(sources)
            msg = (
                f"The target \"{c.target.label}\" ({c.target_id}) is attached "
                f"to more than one component on the source side: {sources}. "
                "Please check your network configuration."
            )
            logger.error(msg)
            raise hlp.TESPyNetworkError(msg)

    def _add_comps(self, *args):
        r"""
        Add to network's component DataFrame from added connections.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            The connections, which have been added to the network. The
            components are extracted from these information.
        """
        # get unique components in new connections
        comps = list({cp for c in args for cp in [c.source, c.target]})
        # add to the dataframe of components
        for comp in comps:
            if comp.label in self.comps.index:
                if self.comps.loc[comp.label, 'object'] == comp:
                    continue
                else:
                    comp_type = comp.__class__.__name__
                    other_obj = self.comps.loc[comp.label, "object"]
                    other_comp_type = other_obj.__class__.__name__
                    msg = (
                        f"The component with the label {comp.label} of type "
                        f"{comp_type} cannot be added to the network as a "
                        f"different component of type {other_comp_type} with "
                        "the same label has already been added. All "
                        "components must have unique values!"
                    )
                    raise hlp.TESPyNetworkError(msg)

            comp_type = comp.__class__.__name__
            self.comps.loc[comp.label, 'comp_type'] = comp_type
            self.comps.loc[comp.label, 'object'] = comp

    def _del_comps(self, comps):
        r"""
        Delete from network's component DataFrame from deleted connections.

        For every component it is checked, if it is still part of other
        connections, which have not been deleted. The component is only
        removed if it cannot be found int the remaining connections.

        Parameters
        ----------
        comps : list
            List of components to potentially be deleted.
        """
        for comp in comps:
            if (
                comp not in self.conns["source"].values and
                comp not in self.conns["target"].values
            ):
                self.comps.drop(comp.label, inplace=True)
                self.results[comp.__class__.__name__].drop(
                    comp.label, inplace=True, errors="ignore"
                )
                msg = f"Deleted component {comp.label} from network."
                logger.debug(msg)

    def add_ude(self, *args):
        r"""
        Add a user defined function to the network.

        Parameters
        ----------
        c : tespy.tools.helpers.UserDefinedEquation
            The objects to be added to the network, UserDefinedEquation objects
            ci :code:`add_ude(c1, c2, c3, ...)`.
        """
        for c in args:
            if not isinstance(c, hlp.UserDefinedEquation):
                msg = (
                    'Must provide tespy.tools.helpers.UserDefinedEquation '
                    'objects as parameters.'
                )
                logger.error(msg)
                raise TypeError(msg)

            elif c.label in self.user_defined_eq:
                msg = (
                    'There is already a UserDefinedEquation with the label '
                    f'{c.label} . The UserDefinedEquation labels must be '
                    'unique within a network'
                )
                logger.error(msg)
                raise ValueError(msg)

            self.user_defined_eq[c.label] = c
            msg = f"Added UserDefinedEquation {c.label} to network."
            logger.debug(msg)

    def del_ude(self, *args):
        """
        Remove a user defined function from the network.

        Parameters
        ----------
        c : tespy.tools.helpers.UserDefinedEquation
            The objects to be deleted from the network,
            UserDefinedEquation objects ci :code:`del_ude(c1, c2, c3, ...)`.
        """
        for c in args:
            del self.user_defined_eq[c.label]
            msg = f"Deleted UserDefinedEquation {c.label} from network."
            logger.debug(msg)

    def add_busses(self, *args):
        r"""
        Add one or more busses to the network.

        Parameters
        ----------
        b : tespy.connections.bus.Bus
            The bus to be added to the network, bus objects bi
            :code:`add_busses(b1, b2, b3, ...)`.
        """
        for b in args:
            if self.check_busses(b):
                self.busses[b.label] = b
                msg = f"Added bus {b.label} to network."
                logger.debug(msg)

                self.results[b.label] = pd.DataFrame(
                    columns=[
                        'component value', 'bus value', 'efficiency',
                        'design value'
                    ],
                    dtype='float64'
                )

    def del_busses(self, *args):
        r"""
        Remove one or more busses from the network.

        Parameters
        ----------
        b : tespy.connections.bus.Bus
            The bus to be removed from the network, bus objects bi
            :code:`add_busses(b1, b2, b3, ...)`.
        """
        for b in args:
            if b in self.busses.values():
                del self.busses[b.label]
                msg = f"Deleted bus {b.label} from network."
                logger.debug(msg)

                del self.results[b.label]

    def _convergence_check(self):
        """Check convergence status of a simulation."""
        msg = 'Calculation did not converge!'
        assert (not self.lin_dep) and self.converged, msg

    def check_busses(self, b):
        r"""
        Checksthe busses to be added for type, duplicates and identical labels.

        Parameters
        ----------
        b : tespy.connections.bus.Bus
            The bus to be checked.
        """
        if isinstance(b, con.Bus):
            if len(self.busses) > 0:
                if b in self.busses.values():
                    msg = f"The network contains the bus {b.label} already."
                    logger.error(msg)
                    raise hlp.TESPyNetworkError(msg)
                elif b.label in self.busses:
                    msg = f"The network already has a bus labeld {b.label}."
                    logger.error(msg)
                    raise hlp.TESPyNetworkError(msg)
                else:
                    return True
            else:
                return True
        else:
            msg = 'Only objects of type bus are allowed in *args.'
            logger.error(msg)
            raise TypeError(msg)

    def check_network(self):
        r"""Check if components are connected properly within the network."""
        if len(self.conns) == 0:
            msg = (
                'No connections have been added to the network, please make '
                'sure to add your connections with the .add_conns() method.'
            )
            logger.error(msg)
            raise hlp.TESPyNetworkError(msg)

        self.check_conns()
        self.init_components()
        self.check_components()
        self.create_fluid_wrapper_branches()

        # network checked
        self.checked = True
        msg = 'Networkcheck successful.'
        logger.info(msg)

    def _create_structure_matrix(self):
        self._structure_matrix = {}
        self._rhs = {}
        self._variable_lookup = {}

        for i, conn in enumerate(self.conns["object"]):
            conn.m.sm_col = 4 * i
            conn.p.sm_col = 4 * i + 1
            conn.h.sm_col = 4 * i + 2
            conn.fluid.sm_col = 4 * i + 3
            for prop in ["m", "p", "h", "fluid"]:
                container = conn.get_attr(prop)
                self._variable_lookup[container.sm_col] = {
                    "connection": conn, "property": prop
                }
                container._reference_container = container
                container._factor = 1
                container._offset = 0

        # not sure if useful. In principle component variables should also be
        # part of the structure matrix in order to presolve them. There should
        # however not be many, that can create direct dependencies to just a
        # single other variable. Currently only 1-1 linear dependent variables
        # are used in the presolving stage
        sm_col_counter = conn.fluid.sm_col + 1

        sum_eq = 0
        for conn in self.conns["object"]:
            conn._preprocess(sum_eq)
            self._structure_matrix.update(conn._structure_matrix)
            self._rhs.update(conn._rhs)
            sum_eq += conn.num_eq

        for cp in self.comps["object"]:
            cp._preprocess(sum_eq)

            # see comment on sm_col above
            for var in cp.vars:
                var.sm_col = sm_col_counter
                sm_col_counter += 1

            self._structure_matrix.update(cp._structure_matrix)
            self._rhs.update(cp._rhs)
            sum_eq += cp.num_eq

        self._variable_dependencies = self._find_linear_dependent_variables(
            self._structure_matrix, self._rhs
        )

    def _find_linear_dependent_variables(self, sparse_matrix, rhs):
        edges_with_factors = []
        rhs_offsets = {}
        eq_idx = {}  # The equation indices keep track of which equations to eliminate

        num_rows = 1 + max([k[0] for k in sparse_matrix.keys()])
        num_cols = 1 + max([k[1] for k in sparse_matrix.keys()])

        # Convert sparse matrix to dense form
        dense_matrix = np.zeros((num_rows, num_cols))
        for idx, value in sparse_matrix.items():
            dense_matrix[idx] = value

        # Extract edges and offsets from rows with two non-zero entries
        for row_idx in range(num_rows):
            non_zero_indices = [col_idx for col_idx, value in enumerate(dense_matrix[row_idx]) if value != 0]
            non_zero_values = [dense_matrix[row_idx, col_idx] for col_idx in non_zero_indices]

            if len(non_zero_indices) == 2:
                col1, col2 = non_zero_indices
                val1, val2 = non_zero_values
                factor = -val1 / val2
                offset = rhs[row_idx] / val2
                edges_with_factors.append((col1, col2, factor))
                rhs_offsets[(col1, col2)] = offset
                eq_idx[(col1, col2)] = row_idx

        # Build adjacency list for the graph
        adjacency_list = {}
        for col1, col2, factor in edges_with_factors:
            if col1 not in adjacency_list:
                adjacency_list[col1] = []
            if col2 not in adjacency_list:
                adjacency_list[col2] = []

            # Check for duplicate edges
            if col2 in [neighbor for neighbor, _ in adjacency_list[col1]]:
                msg = f"There already is a linear link between the variables {col1}, {col2}"
                raise hlp.TESPyNetworkError(msg)

            # Add edge with factor
            adjacency_list[col1].append((col2, factor))
            adjacency_list[col2].append((col1, 1 / factor))  # Add reverse edge with reciprocal factor

        print(adjacency_list)
        # Detect cycles (to check for circular dependencies)
        visited = set()
        parent_map = {}

        def dfs_cycle(node, parent):
            visited.add(node)
            for neighbor, _ in adjacency_list.get(node, []):
                if neighbor not in visited:
                    parent_map[neighbor] = node
                    if dfs_cycle(neighbor, node):
                        return True
                elif neighbor != parent:  # A back edge is found
                    return True
            return False

        for node in adjacency_list:
            if node not in visited:
                if dfs_cycle(node, None):
                    variable_names = self._map_column_indices_to_variable_names(list(parent_map.keys()))
                    msg = f"A circular dependency between the following variables has been detected: {variable_names}"
                    raise hlp.TESPyNetworkError(msg)

        # Find connected components and compute factors/offsets
        visited.clear()
        variables_factors_offsets = []

        def dfs_component(node, reference, current_factor, current_offset):
            """DFS to calculate factors and offsets relative to the reference variable."""
            stack = [(node, current_factor, current_offset)]
            factors = {node: current_factor}
            offsets = {node: current_offset}
            equation_indices = {}

            while stack:
                curr_node, curr_factor, curr_offset = stack.pop()
                visited.add(curr_node)

                for neighbor, edge_factor in adjacency_list.get(curr_node, []):
                    if neighbor not in factors:  # Process unvisited neighbor
                        # Calculate edge offset
                        edge_offset = rhs_offsets.get((curr_node, neighbor), 0.0) or -rhs_offsets.get((neighbor, curr_node), 0.0)
                        # Determine which equation to use
                        if (neighbor, curr_node) in rhs_offsets:
                            equation_indices[(neighbor, curr_node)] = eq_idx[(neighbor, curr_node)]
                        else:
                            equation_indices[(curr_node, neighbor)] = eq_idx[(curr_node, neighbor)]

                        # Compute new factor and offset
                        new_factor = curr_factor * edge_factor
                        new_offset = curr_offset * edge_factor + edge_offset

                        # Store and continue traversal
                        factors[neighbor] = new_factor
                        offsets[neighbor] = new_offset
                        stack.append((neighbor, new_factor, new_offset))

            return factors, offsets, equation_indices

        # Process each connected component
        for node in adjacency_list:
            if node not in visited:
                reference = node
                factors, offsets, equation_indices = dfs_component(reference, reference, 1.0, 0.0)

                variables_factors_offsets.append({
                    'variables': list(factors.keys()),
                    'reference': reference,
                    'factors': factors,
                    'offsets': offsets,
                    'equation_indices': equation_indices
                })

        return variables_factors_offsets

    def _map_column_indices_to_variable_names(self, cycle):
        return cycle

    def create_fluid_wrapper_branches(self):

        self.fluid_wrapper_branches = {}
        mask = self.comps["comp_type"].isin(
            ["Source", "CycleCloser", "WaterElectrolyzer", "FuelCell"]
        )
        start_components = self.comps["object"].loc[mask]

        for start in start_components:
            self.fluid_wrapper_branches.update(start.start_fluid_wrapper_branch())

        merged = self.fluid_wrapper_branches.copy()
        for branch_name, branch_data in self.fluid_wrapper_branches.items():
            if branch_name not in merged:
                continue
            for ob_name, ob_data in self.fluid_wrapper_branches.items():
                if ob_name != branch_name:
                    common_connections = list(
                        set(branch_data["connections"])
                        & set(ob_data["connections"])
                    )
                    if len(common_connections) > 0 and ob_name in merged:
                        merged[branch_name]["connections"] = list(
                            set(branch_data["connections"] + ob_data["connections"])
                        )
                        merged[branch_name]["components"] = list(
                            set(branch_data["components"] + ob_data["components"])
                        )
                        del merged[ob_name]
                        break

        self.fluid_wrapper_branches = merged

    def init_components(self):
        r"""Set up necessary component information."""
        for comp in self.comps["object"]:
            # get incoming and outgoing connections of a component
            sources = self.conns[self.conns['source'] == comp]
            sources = sources['source_id'].sort_values().index.tolist()
            targets = self.conns[self.conns['target'] == comp]
            targets = targets['target_id'].sort_values().index.tolist()
            # save the incoming and outgoing as well as the number of
            # connections as component attribute
            comp.inl = self.conns.loc[targets, 'object'].tolist()
            comp.outl = self.conns.loc[sources, 'object'].tolist()
            comp.num_i = len(comp.inlets())
            comp.num_o = len(comp.outlets())

            # set up restults and specification dataframes
            comp_type = comp.__class__.__name__
            if comp_type not in self.results:
                cols = [
                    col for col, data in comp.parameters.items()
                    if isinstance(data, dc_cp)
                ]
                self.results[comp_type] = pd.DataFrame(
                    columns=cols, dtype='float64'
                )
            if comp_type not in self.specifications:

                cols, groups, chars = [], [], []
                for col, data in comp.parameters.items():
                    if isinstance(data, dc_cp):
                        cols += [col]
                    elif isinstance(data, dc_gcp) or isinstance(data, dc_gcc):
                        groups += [col]
                    elif isinstance(data, dc_cc) or isinstance(data, dc_cm):
                        chars += [col]
                self.specifications[comp_type] = {
                    'groups': pd.DataFrame(columns=groups, dtype='bool'),
                    'chars': pd.DataFrame(columns=chars, dtype='object'),
                    'variables': pd.DataFrame(columns=cols, dtype='bool'),
                    'properties': pd.DataFrame(columns=cols, dtype='bool')
                }

    def check_components(self):
        # count number of incoming and outgoing connections and compare to
        # expected values
        for comp in self.comps['object']:
            counts = (self.conns[['source', 'target']] == comp).sum()

            if counts["source"] != comp.num_o:
                msg = (
                    f"The component {comp.label} is missing "
                    f"{comp.num_o - counts['source']} outgoing connections. "
                    "Make sure all outlets are connected and all connections "
                    "have been added to the network."
                )
                logger.error(msg)
                # raise an error in case network check is unsuccesful
                raise hlp.TESPyNetworkError(msg)
            elif counts["target"] != comp.num_i:
                msg = (
                    f"The component {comp.label} is missing "
                    f"{comp.num_i - counts['target']} incoming connections. "
                    "Make sure all inlets are connected and all connections "
                    "have been added to the network."
                )
                logger.error(msg)
                # raise an error in case network check is unsuccesful
                raise hlp.TESPyNetworkError(msg)

    def _map_linear_dependent_variables(self):

        for linear_dependents in self._variable_dependencies:
            is_fluid_list = [
                self._variable_lookup[var]["property"] == "fluid"
                for var in linear_dependents["variables"]
            ]
            # this is suspicious, because it was found, that fluid variables
            # are linear dependent to non-fluid variables
            if not all(is_fluid_list) and any(is_fluid_list):
                msg = "We found linear dependency between fluid and non-fluid variables."
                raise hlp.TESPyNetworkError(msg)

            # no matter what object we use, they will all be non fluid or fluid
            reference = linear_dependents["reference"]
            reference_property = self._variable_lookup[reference]["property"]
            if reference_property == "fluid":
                reference_container = dc_flu()
            else:
                reference_container = dc_fp(_is_var=True)

            self.linear_dependent_variables[reference] = reference_container

            for dependent in linear_dependents["variables"]:
                dependent_object = (
                    self._variable_lookup[dependent]["connection"].get_attr(
                        self._variable_lookup[dependent]["property"]
                        )
                    )
                dependent_object._reference_container = reference_container
                dependent_object._factor = linear_dependents["factors"][dependent]
                dependent_object._offset = linear_dependents["offsets"][dependent]

    def initialise(self):
        r"""
        Initilialise the network depending on calclation mode.

        Design

        - Generic fluid composition and fluid property initialisation.
        - Starting values from initialisation path if provided.

        Offdesign

        - Check offdesign path specification.
        - Set component and connection design point properties.
        - Switch from design/offdesign parameter specification.
        """
        # keep track of the number of bus, component and connection equations
        # as well as number of component variables
        self.num_bus_eq = 0
        self.num_comp_eq = 0
        self.num_conn_eq = 0
        self.num_vars = 0
        self.num_comp_vars = 0
        self.num_conn_vars = 0
        self.variables_dict = {}
        self.linear_dependent_variables = {}

        # in multiprocessing copies are made of all connections
        # the mass flow branches and fluid branches hold references to
        # connections from the original run (where network.checked is False)
        # The assignment of variable spaces etc. is however made on the
        # copies of the connections which do not correspond to the mass flow
        # branches and fluid branches anymore. So the topology simplification
        # does not actually apply to the copied network, therefore the
        # branches have to be recreated for this case. We can detect that by
        # checking whether a network holds a massflow branch with some
        # connections and compare that with the connection object actually
        # present in the network
        self.propagate_fluid_wrappers()
        # self.impose_unit_system()

        if self.mode == 'offdesign':
            self.redesign = True
            if self.design_path is None:
                # must provide design_path
                msg = "Please provide a design_path for offdesign mode."
                logger.error(msg)
                raise hlp.TESPyNetworkError(msg)

            # load design case
            if self.new_design:
                self.init_offdesign_params()

            self.init_offdesign()

        else:
            # reset any preceding offdesign calculation
            self.init_design()

        # self.simplify_problem() will
        # create the structure matrix + rhs of the problem
        # assign linear dependents to groups
        # presolve the fluid vectors
        # determine, which connection equations need to be applied
        # determine, which component equations need to be applied
        # set up jacobian and residual vector structures for connections
        # components, user defined equations and buses
        self.init_set_properties()
        self._create_structure_matrix()
        self._map_linear_dependent_variables()
        self._presolve_fluid_vectors()
        self._prepare_problem()

        # generic fluid property initialisation
        self.init_properties()

        msg = 'Network initialised.'
        logger.info(msg)

    def propagate_fluid_wrappers(self):

        for branch_data in self.fluid_wrapper_branches.values():
            all_connections = [c for c in branch_data["connections"]]

            any_fluids_set = []
            engines = {}
            back_ends = {}
            for c in all_connections:
                for f in c.fluid.is_set:
                    any_fluids_set += [f]
                    if f in c.fluid.engine:
                        engines[f] = c.fluid.engine[f]
                    if f in c.fluid.back_end:
                        back_ends[f] = c.fluid.back_end[f]

            mixing_rules = [
                c.mixing_rule for c in all_connections
                if c.mixing_rule is not None
            ]
            mixing_rule = set(mixing_rules)
            if len(mixing_rule) > 1:
                msg = "You have provided more than one mixing rule."
                raise hlp.TESPyNetworkError(msg)
            elif len(mixing_rule) == 0:
                mixing_rule = set(["ideal-cond"])

            if not any_fluids_set:
                msg = "You are missing fluid specifications."
            any_fluids = [f for c in all_connections for f in c.fluid.val]
            any_fluids0 = [f for c in all_connections for f in c.fluid.val]

            potential_fluids = set(any_fluids_set + any_fluids + any_fluids0)
            num_potential_fluids = len(potential_fluids)
            if num_potential_fluids == 0:
                msg = (
                    "The follwing connections of your network are missing any "
                    "kind of fluid composition information:"
                    f"{', '.join([c.label for c in all_connections])}."
                )
                raise hlp.TESPyNetworkError(msg)

            for c in all_connections:
                c.mixing_rule = list(mixing_rule)[0]
                c._potential_fluids = potential_fluids
                if num_potential_fluids == 1:
                    f = list(potential_fluids)[0]
                    c.fluid.val[f] = 1

                else:
                    for f in potential_fluids:
                        if (f not in c.fluid.is_set and f not in c.fluid.val and f not in c.fluid.val0):
                            c.fluid.val[f] = 1 / len(potential_fluids)
                        elif f not in c.fluid.is_set and f not in c.fluid.val and f in c.fluid.val0:
                            c.fluid.val[f] = c.fluid.val0[f]

                for f, engine in engines.items():
                    c.fluid.engine[f] = engine
                for f, back_end in back_ends.items():
                    c.fluid.back_end[f] = back_end

                c._create_fluid_wrapper()

    def _presolve_fluid_vectors(self):

        # right now, this ignores potential factors and offsets between the
        # fluids.
        # On top of that, branches of constant fluid composition are not
        # caught, if they only consist of a single connection!
        for linear_dependents in self._variable_dependencies:
            reference = linear_dependents["reference"]

            if self._variable_lookup[reference]["property"] != "fluid":
                continue

            all_connections = [
                self._variable_lookup[var]["connection"]
                for var in linear_dependents["variables"]
            ]

            reference_container = self.linear_dependent_variables[reference]
            reference_conn = all_connections[0]

            fluid_specs = [f for c in all_connections for f in c.fluid.is_set]
            if len(fluid_specs) == 0:

                if len(reference_conn._potential_fluids) > 1:
                    reference_container.is_var = {f for f in reference_conn.fluid.val}
                else:
                    reference_container.val[list(reference_conn._potential_fluids)[0]] = 1

            elif len(fluid_specs) != len(set(fluid_specs)):
                msg = (
                    "The mass fraction of a single fluid cannot be specified "
                    "twice within a branch."
                )
                raise hlp.TESPyNetworkError(msg)
            else:
                fixed_fractions = {
                    f: c.fluid.val[f]
                    for c in all_connections
                    for f in fluid_specs
                    if f in c.fluid.is_set
                }
                mass_fraction_sum = sum(fixed_fractions.values())
                if mass_fraction_sum > 1 + ERR:
                    msg = "Total mass fractions within a branch cannot exceed 1"
                    raise ValueError(msg)
                elif mass_fraction_sum < 1 - ERR:
                    # set the fluids with specified mass fraction
                    # remaining fluids are variable, create wrappers for them
                    all_fluids = reference_container.val.keys()
                    num_remaining_fluids = len(all_fluids) - len(fixed_fractions)
                    if num_remaining_fluids == 1:
                        missing_fluid = list(
                            reference_container.val.keys() - fixed_fractions.keys()
                        )[0]
                        fixed_fractions[missing_fluid] = 1 - mass_fraction_sum
                        variable = set()
                    else:
                        missing_fluids = (
                            reference_container.val.keys() - fixed_fractions.keys()
                        )
                        variable = {f for f in missing_fluids}

                else:
                    # fluid mass fraction is 100 %, all other fluids are 0 %
                    all_fluids = reference_container.val.keys()
                    remaining_fluids = (
                        reference_container.val.keys() - fixed_fractions.keys()
                    )
                    for f in remaining_fluids:
                        fixed_fractions[f] = 0

                    variable = set()

                reference_container.val.update(fixed_fractions)
                reference_container.is_set = {f for f in fixed_fractions}
                reference_container.is_var = variable
                # this seems to be a problem in some cases, e.g. the basic
                # gas turbine tutorial
                # num_var = len(variable)
                # for f in variable:
                #     reference_container.val[f] = (1 - mass_fraction_sum) / num_var

            [c.build_fluid_data() for c in all_connections]
            for fluid in reference_container.is_var:
                reference_container._J_col[fluid] = self.num_conn_vars
                self.variables_dict[self.num_conn_vars] = {
                    "obj": reference_container,
                    "variable": "fluid",
                    "fluid": fluid
                }
                self.num_conn_vars += 1

    def init_set_properties(self):
        """Specification of SI values for user set values."""
        self.all_fluids = []
        # fluid property values
        for c in self.conns['object']:
            self.all_fluids += c.fluid.val.keys()

            if not self.init_previous:
                c.good_starting_values = False

            for key in c.property_data:
                # read unit specifications
                prop = key.split("_ref")[0]
                if "fluid" in key:
                    continue
                elif key == 'Td_bp':
                    c.get_attr(key).unit = self.get_attr('T_unit')
                else:
                    c.get_attr(key).unit = self.get_attr(f"{prop}_unit")
                # set SI value
                if c.get_attr(key).is_set:
                    # this could be externalized to either
                    # the connections class or actually the data containers
                    if "ref" in key:
                        if prop == 'T':
                            c.get_attr(key).ref.delta_SI = hlp.convert_to_SI(
                                'Td_bp', c.get_attr(key).ref.delta,
                                c.get_attr(prop).unit
                            )
                        else:
                            c.get_attr(key).ref.delta_SI = hlp.convert_to_SI(
                                prop, c.get_attr(key).ref.delta,
                                c.get_attr(prop).unit
                            )
                    else:
                        c.get_attr(key)._val_SI = hlp.convert_to_SI(
                            key, c.get_attr(key)._val, c.get_attr(key).unit
                        )
                if key in ["m", "p", "h"]:
                    c.get_attr(key)._is_var = not c.get_attr(key).is_set

        # set up results dataframe for connections
        # this should be done based on the connections
        properties = list(fpd.keys())
        self.all_fluids = set(self.all_fluids)
        cols = (
            [col for prop in properties for col in [prop, f"{prop}_unit"]]
            + list(self.all_fluids) + ['phase']
        )
        self.results['Connection'] = pd.DataFrame(columns=cols, dtype='float64')
        # include column for fluid balance in specs dataframe
        self.specifications['Connection'] = pd.DataFrame(
            columns=cols + ['balance'], dtype='bool'
        )
        cols = ["m_ref", "p_ref", "h_ref", "T_ref", "v_ref"]
        self.specifications['Ref'] = pd.DataFrame(columns=cols, dtype='bool')

        msg = (
            "Updated fluid property SI values and fluid mass fraction for user "
            "specified connection parameters."
        )
        logger.debug(msg)

    def _assign_variable_space(self, c):
        for key in ["m", "p", "h"]:
            variable = c.get_attr(key)
            if (
                variable.is_var()
                and variable._reference_container not in self._conn_variables
            ):
                container = variable._reference_container
                container._J_col = self.num_conn_vars
                self.variables_dict[self.num_conn_vars] = {
                    "obj": container, "variable": key
                }
                print(container._J_col, container, container.J_col())
                self._conn_variables += [container]
                self.num_conn_vars += 1

    def init_design(self):
        r"""
        Initialise a design calculation.

        Offdesign parameters are unset, design parameters are set. If
        :code:`local_offdesign` is :code:`True` for connections or components,
        the design point information are read from the .csv-files in the
        respective :code:`design_path`. In this case, the design values are
        unset, the offdesign values set.
        """
        # connections
        self._conn_variables = []
        _local_designs = {}
        for c in self.conns['object']:
            # read design point information of connections with
            # local_offdesign activated from their respective design path
            if c.local_offdesign:
                if c.design_path is None:
                    msg = (
                        "The parameter local_offdesign is True for the "
                        f"connection {c.label}, an individual design_path must "
                        "be specified in this case!"
                    )
                    logger.error(msg)
                    raise hlp.TESPyNetworkError(msg)

                # unset design parameters
                for var in c.design:
                    c.get_attr(var).is_set = False
                # set offdesign parameters
                for var in c.offdesign:
                    c.get_attr(var).is_set = True

                # read design point information
                msg = (
                    "Reading individual design point information for "
                    f"connection {c.label} from {c.design_path}/connections.csv."
                )
                logger.debug(msg)
                if c.design_path not in _local_designs:
                    _local_designs[c.design_path] = self.init_read_connections(
                        c.design_path
                    )
                df = _local_designs[c.design_path]
                # write data to connections
                self.init_conn_design_params(c, df)

            else:
                # unset all design values
                c.m.design = np.nan
                c.p.design = np.nan
                c.h.design = np.nan
                c.fluid.design = {}

                c.new_design = True

                # switch connections to design mode
                if self.redesign:
                    for var in c.design:
                        c.get_attr(var).is_set = True

                    for var in c.offdesign:
                        c.get_attr(var).is_set = False

        # unset design values for busses, count bus equations and
        # reindex bus dictionary
        for b in self.busses.values():
            b.comps['P_ref'] = np.nan

        series = pd.Series(dtype='float64')
        _local_design_paths = {}
        for cp in self.comps['object']:
            c = cp.__class__.__name__
            # read design point information of components with
            # local_offdesign activated from their respective design path
            if cp.local_offdesign:
                if cp.design_path is not None:
                    # read design point information
                    path = os.path.join(cp.design_path, "components", f"{c}.csv")
                    msg = (
                        f"Reading design point information for component "
                        f"{cp.label} of type {c} from path {path}."
                    )
                    logger.debug(msg)
                    if path not in _local_design_paths:
                        _local_design_paths[path] = pd.read_csv(
                            path, sep=';', decimal='.', index_col=0
                        )
                    data = _local_design_paths[path].loc[cp.label]
                    # write data
                    self.init_comp_design_params(cp, data)

                # unset design parameters
                for var in cp.design:
                    cp.get_attr(var).is_set = False

                # set offdesign parameters
                switched = False
                msg = 'Set component attributes '

                for var in cp.offdesign:
                    # set variables provided in .offdesign attribute
                    data = cp.get_attr(var)
                    data.is_set = True

                    # take nominal values from design point
                    if isinstance(data, dc_cp):
                        cp.get_attr(var).val = cp.get_attr(var).design
                        switched = True
                        msg += var + ", "

                if switched:
                    msg = f"{msg[:-2]} to design value at component {cp.label}."
                    logger.debug(msg)

                cp.new_design = False

            else:
                # switch connections to design mode
                if self.redesign:
                    for var in cp.design:
                        cp.get_attr(var).is_set = True

                    for var in cp.offdesign:
                        cp.get_attr(var).is_set = False

                cp.set_parameters(self.mode, series)

    def _prepare_problem(self):
        self._presolve_linear_dependents()

        # iteratively check presolvable fluid properties
        # and distribute presolved variables to all linear dependents
        # until the number of variables does not change anymore
        number_variables = sum([
            conn.get_attr(key).is_var()
            for conn in self.conns['object']
            for key in ["m", "p", "h"]
        ])
        while True:
            for c in self.conns['object']:
                if not c.fluid.is_var:
                    c._presolve()
            self._presolve_linear_dependents()
            reduced_variables = [
                conn.get_attr(key).is_var()
                for conn in self.conns['object']
                for key in ["m", "p", "h"]
            ]
            reduced_variables = sum(reduced_variables)
            if reduced_variables == number_variables:
                break

            number_variables = reduced_variables

        # set up the actual list of equations for connections, components,
        # buses and user defined equations
        presolved_equations = [
            indices
            for dependents in self._variable_dependencies
            for indices in dependents["equation_indices"].values()
        ]

        for c in self.conns['object']:
            c._prepare_for_solver(presolved_equations)
            self._assign_variable_space(c)

        for cp in self.comps['object']:
            c = cp.__class__.__name__
            # component initialisation
            cp._prepare_for_solver(presolved_equations)

            for spec in self.specifications[c].keys():
                if len(cp.get_attr(self.specifications['lookup'][spec])) > 0:
                    self.specifications[c][spec].loc[cp.label] = (
                        cp.get_attr(self.specifications['lookup'][spec])
                    )

            # count number of component equations and variables
            i = self.num_conn_vars + self.num_comp_vars
            for container, name in cp.vars.items():
                self.variables_dict[i] = {"obj": container, "variable": name}
                i += 1
            self.num_comp_vars += cp.num_vars
            self.num_comp_eq += cp.num_eq

        for b in self.busses.values():
            self.busses[b.label] = b
            self.num_bus_eq += b.P.is_set * 1

    def _presolve_linear_dependents(self):
        for linear_dependents in self._variable_dependencies:
            reference = linear_dependents["reference"]

            if self._variable_lookup[reference]["property"] == "fluid":
                continue

            all_containers = [
                self._variable_lookup[var]["connection"].get_attr(
                    self._variable_lookup[var]["property"]
                ) for var in linear_dependents["variables"]
            ]

            number_specifications = sum([not c._is_var for c in all_containers])
            if number_specifications > 1:
                variables_properties = [
                    f"({self._variable_lookup[var]['connection'].label}: "
                    f"{self._variable_lookup[var]['property']})"
                    for var in linear_dependents["variables"]
                ]
                msg = (
                    "You specified more than one variable of the linear "
                    f"dependent variables: {', '.join(variables_properties)}."
                )
                raise hlp.TESPyNetworkError(msg)
            elif number_specifications == 1:
                reference_container = self.linear_dependent_variables[reference]
                reference_container._is_var = False
                specification = [c for c in all_containers if not c._is_var][0]
                reference_container._val_SI = specification.reference_val_SI()

    def init_offdesign_params(self):
        r"""
        Read design point information from specified :code:`design_path`.

        If a :code:`design_path` has been specified individually for components
        or connections, the data will be read from the specified individual
        path instead.

        Note
        ----
        The methods
        :py:meth:`tespy.networks.network.Network.init_comp_design_params`
        (components) and the
        :py:meth:`tespy.networks.network.Network.init_conn_design_params`
        (connections) handle the parameter specification.
        """
        # components without any parameters
        components_with_parameters = [
            cp.label for cp in self.comps["object"] if len(cp.parameters) > 0
        ]
        # fetch all components, reindex with label
        df_comps = self.comps.loc[components_with_parameters].copy()
        # iter through unique types of components (class names)
        for c in df_comps['comp_type'].unique():
            path = os.path.join(self.design_path, "components", f"{c}.csv")
            msg = (
                f"Reading design point information for components of type {c} "
                f"from path {path}."
            )
            logger.debug(msg)
            df = pd.read_csv(path, sep=';', decimal='.', index_col=0)
            df.index = df.index.astype(str)

            # iter through all components of this type and set data
            _individual_design_paths = {}
            for c_label in df.index:
                comp = self.comps.loc[c_label, 'object']
                # read data of components with individual design_path
                if comp.design_path is not None:
                    path_c = os.path.join(
                        comp.design_path, "components", f"{c}.csv"
                    )
                    msg = (
                        f"Reading design point information for component "
                        f"{comp.label} of type {c} from path {path_c}."
                    )
                    logger.debug(msg)
                    if path_c not in _individual_design_paths:
                        _individual_design_paths[path_c] = pd.read_csv(
                            path_c, sep=';', decimal='.', index_col=0
                        )
                        _individual_design_paths[path_c].index = (
                            _individual_design_paths[path_c].index.astype(str)
                        )
                    data = _individual_design_paths[path_c].loc[comp.label]

                else:
                    data = df.loc[comp.label]

                # write data to components
                self.init_comp_design_params(comp, data)

        msg = 'Done reading design point information for components.'
        logger.debug(msg)

        if len(self.busses) > 0:
            path = os.path.join(self.design_path, "busses.json")
            with open(path, "r", encoding="utf-8") as f:
                bus_data = json.load(f)

            for b in bus_data:
                for comp, value in bus_data[b].items():
                    comp = self.get_comp(comp)
                    self.busses[b].comps.loc[comp, "P_ref"] = float(value)

        # read connection design point information
        df = self.init_read_connections(self.design_path)

        # iter through connections
        for c in self.conns['object']:

            # read data of connections with individual design_path
            if c.design_path is not None:
                msg = (
                    "Reading connection design point information for "
                    f"{c.label} from {c.design_path}/connections.csv."
                )
                logger.debug(msg)
                df_c = self.init_read_connections(c.design_path)

                # write data
                self.init_conn_design_params(c, df_c)

            else:
                # write data
                self.init_conn_design_params(c, df)

        msg = 'Done reading design point information for connections.'
        logger.debug(msg)

    def init_comp_design_params(self, component, data):
        r"""
        Write design point information to components.

        Parameters
        ----------
        component : tespy.components.component.Component
            Write design point information to this component.

        data : pandas.core.series.Series, pandas.core.frame.DataFrame
            Design point information.
        """
        # write component design data
        component.set_parameters(self.mode, data)

    def init_conn_design_params(self, c, df):
        r"""
        Write design point information to connections.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            Write design point information to this connection.

        df : pandas.core.frame.DataFrame
            Dataframe containing design point information.
        """
        # match connection (source, source_id, target, target_id) on
        # connection objects of design file
        df.index = df.index.astype(str)
        if c.label not in df.index:
            # no matches in the connections of the network and the design files
            msg = (
                f"Could not find connection '{c.label}' in design case. Please "
                "make sure no connections have been modified or components "
                "have been relabeled for your offdesign calculation."
            )
            logger.exception(msg)
            raise hlp.TESPyNetworkError(msg)

        conn = df.loc[c.label]
        for var in fpd.keys():
            c.get_attr(var).design = hlp.convert_to_SI(
                var, float(conn[var]), conn[f"{var}_unit"]
            )
        if c.m.design != 0.0:
            c.vol.design = c.v.design / c.m.design
        else:
            c.vol.design = math.inf
        for fluid in c.fluid.val:
            c.fluid.design[fluid] = float(conn[fluid])

    def init_conn_params_from_path(self, c, df):
        r"""
        Write parameter information from init_path to connections.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            Write init path information to this connection.

        df : pandas.core.frame.DataFrame
            Dataframe containing init path information.
        """
        # match connection (source, source_id, target, target_id) on
        # connection objects of design file
        df.index = df.index.astype(str)
        if c.label not in df.index:
            # no matches in the connections of the network and the design files
            msg = f"Could not find connection {c.label} in init path file."
            logger.debug(msg)
            return

        conn = df.loc[c.label]

        for prop in ['m', 'p', 'h']:
            data = c.get_attr(prop)
            data.val0 = float(conn[prop])
            data.unit = conn[prop + '_unit']

        for fluid in c.fluid.is_var:
            c.fluid.val[fluid] = float(conn[fluid])
            c.fluid.val0[fluid] = float(c.fluid.val[fluid])

        c.good_starting_values = True

    def init_offdesign(self):
        r"""
        Switch components and connections from design to offdesign mode.

        Note
        ----
        **components**

        All parameters stated in the component's attribute :code:`cp.design`
        will be unset and all parameters stated in the component's attribute
        :code:`cp.offdesign` will be set instead.

        Additionally, all component parameters specified as variables are
        unset and the values from design point are set.

        **connections**

        All parameters given in the connection's attribute :code:`c.design`
        will be unset and all parameters stated in the connections's attribute
        :code:`cp.offdesign` will be set instead. This does also affect
        referenced values!
        """
        self._conn_variables = []
        for c in self.conns['object']:
            if not c.local_design:
                # switch connections to offdesign mode
                for var in c.design:
                    param = c.get_attr(var)
                    param.is_set = False
                    if f"{var}_ref" in c.property_data:
                        c.get_attr(f"{var}_ref").is_set = False

                for var in c.offdesign:
                    param = c.get_attr(var)
                    param.is_set = True
                    param.val_SI = param.design

                c.new_design = False

        msg = 'Switched connections from design to offdesign.'
        logger.debug(msg)

        for cp in self.comps['object']:
            if not cp.local_design:
                # unset variables provided in .design attribute
                for var in cp.design:
                    cp.get_attr(var).is_set = False

                switched = False
                msg = 'Set component attributes '

                for var in cp.offdesign:
                    # set variables provided in .offdesign attribute
                    data = cp.get_attr(var)
                    data.is_set = True

                    # take nominal values from design point
                    if isinstance(data, dc_cp):
                        cp.get_attr(var).val = cp.get_attr(var).design
                        switched = True
                        msg += var + ', '

                if switched:
                    msg = f"{msg[:-2]} to design value at component {cp.label}."
                    logger.debug(msg)

        msg = 'Switched components from design to offdesign.'
        logger.debug(msg)

    def init_properties(self):
        """
        Initialise the fluid properties on every connection of the network.

        - Set generic starting values for mass flow, enthalpy and pressure if
          not user specified, read from :code:`init_path` or available from
          previous calculation.
        - For generic starting values precalculate enthalpy value at points of
          given temperature, vapor mass fraction, temperature difference to
          boiling point or fluid state.
        """
        if self.init_path is not None:
            df = self.init_read_connections(self.init_path)
        # improved starting values for referenced connections,
        # specified vapour content values, temperature values as well as
        # subccooling/overheating and state specification
        for c in self.conns['object']:
            c.build_fluid_data()
            if self.init_path is not None:
                self.init_conn_params_from_path(c, df)

            if sum(c.fluid.val.values()) == 0:
                msg = (
                    'The starting value for the fluid composition of the '
                    f'connection {c.label} is empty. This might lead to issues '
                    'in the initialisation and solving process as fluid '
                    'property functions can not be called. Make sure you '
                    'specified a fluid composition in all parts of the network.'
                )
                logger.warning(msg)

            for key in ['m', 'p', 'h']:
                if c.get_attr(key).is_var():
                    if not c.good_starting_values:
                        self.init_val0(c, key)
                    c.get_attr(key)._val_SI = hlp.convert_to_SI(
                        key, c.get_attr(key)._val0, c.get_attr(key).unit
                    )
                    c.get_attr(key)._reference_container._val_SI = c.get_attr(key).reference_val_SI()

            self.init_count_connections_parameters(c)

        for c in self.conns['object']:
            if not c.good_starting_values:
                self.init_precalc_properties(c)

            # starting values for specified subcooling/overheating
            # and state specification. These should be recalculated even with
            # good starting values, for example, when one exchanges enthalpy
            # with boiling point temperature difference.
            if (c.Td_bp.is_set or c.state.is_set) and c.h.is_var:
                if ((c.Td_bp.val_SI > 0 and c.Td_bp.is_set) or
                        (c.state.val == 'g' and c.state.is_set)):
                    h = fp.h_mix_pQ(c.p.val_SI, 1, c.fluid_data)
                    if c.h.val_SI < h:
                        c.h.val_SI = h * 1.001
                elif ((c.Td_bp.val_SI < 0 and c.Td_bp.is_set) or
                      (c.state.val == 'l' and c.state.is_set)):
                    h = fp.h_mix_pQ(c.p.val_SI, 0, c.fluid_data)
                    if c.h.val_SI > h:
                        c.h.val_SI = h * 0.999

        msg = 'Generic fluid property specification complete.'
        logger.debug(msg)

    def init_count_connections_parameters(self, c):
        """
        Count the number of parameters set on a connection.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            Connection count parameters of.
        """
        # variables 0 to 9: fluid properties
        local_vars = list(fpd.keys())
        row = [c.get_attr(var).is_set for var in local_vars]
        # write information to specifaction dataframe
        self.specifications['Connection'].loc[c.label, local_vars] = row

        row = [c.get_attr(var).is_set for var in self.specifications['Ref'].columns]
        # write refrenced value information to specifaction dataframe
        self.specifications['Ref'].loc[c.label] = row

        # variables 9 to last but one: fluid mass fractions
        fluids = list(self.all_fluids)
        row = [True if f in c.fluid.is_set else False for f in fluids]
        self.specifications['Connection'].loc[c.label, fluids] = row

        # last one: fluid balance specification
        self.specifications['Connection'].loc[
            c.label, 'balance'] = c.fluid_balance.is_set

        # get number of equations
        self.num_conn_eq += c.num_eq

    def init_precalc_properties(self, c):
        """
        Precalculate enthalpy values for connections.

        Precalculation is performed only if temperature, vapor mass fraction,
        temperature difference to boiling point or phase is specified.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            Connection to precalculate values for.
        """
        # starting values for specified vapour content or temperature
        if c.h.is_var():
            if c.x.is_set:
                try:
                    c.h._val_SI = fp.h_mix_pQ(c.p.val_SI(), c.x.val_SI(), c.fluid_data, c.mixing_rule)
                except ValueError:
                    pass

            if c.T.is_set:
                try:
                    c.h._val_SI = fp.h_mix_pT(c.p.val_SI(), c.T.val_SI(), c.fluid_data, c.mixing_rule)
                except ValueError:
                    pass

    def init_val0(self, c, key):
        r"""
        Set starting values for fluid properties.

        The component classes provide generic starting values for their inlets
        and outlets.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            Connection to initialise.
        """
        if np.isnan(c.get_attr(key)._val0):
            # starting value for mass flow is random between 1 and 2 kg/s
            # (should be generated based on some hash maybe?)
            if key == 'm':
                c.get_attr(key)._val0 = float(np.random.random() + 1)

            # generic starting values for pressure and enthalpy
            else:
                # retrieve starting values from component information
                val_s = c.source.initialise_source(c, key)
                val_t = c.target.initialise_target(c, key)

                if val_s == 0 and val_t == 0:
                    if key == 'p':
                        c.get_attr(key)._val0 = 1e5
                    elif key == 'h':
                        c.get_attr(key)._val0 = 1e6

                elif val_s == 0:
                    c.get_attr(key)._val0 = val_t
                elif val_t == 0:
                    c.get_attr(key)._val0 = val_s
                else:
                    c.get_attr(key)._val0 = (val_s + val_t) / 2

                # change value according to specified unit system
                c.get_attr(key)._val0 = hlp.convert_from_SI(
                    key, c.get_attr(key)._val0, self.get_attr(key + '_unit')
                )

    @staticmethod
    def init_read_connections(base_path):
        r"""
        Read connection information from base_path.

        Parameters
        ----------
        base_path : str
            Path to network information.
        """
        path = os.path.join(base_path, 'connections.csv')
        msg = (
            f"Reading design point information for connections from {path}."
        )
        logger.debug(msg)
        df = pd.read_csv(path, index_col=0, delimiter=';', decimal='.')
        return df

    def solve(self, mode, init_path=None, design_path=None,
              max_iter=50, min_iter=4, init_only=False, init_previous=True,
              use_cuda=False, print_results=True, prepare_fast_lane=False):
        r"""
        Solve the network.

        - Check network consistency.
        - Initialise calculation and preprocessing.
        - Perform actual calculation.
        - Postprocessing.

        It is possible to check programatically, if a network was solved
        successfully with the `.converged` property.

        Parameters
        ----------
        mode : str
            Choose from 'design' and 'offdesign'.

        init_path : str
            Path to the folder, where your network was saved to, e.g.
            saving to :code:`nw.save('myplant/tests')` would require loading
            from :code:`init_path='myplant/tests'`.

        design_path : str
            Path to the folder, where your network's design case was saved to,
            e.g. saving to :code:`nw.save('myplant/tests')` would require
            loading from :code:`design_path='myplant/tests'`.

        max_iter : int
            Maximum number of iterations before calculation stops, default: 50.

        min_iter : int
            Minimum number of iterations before calculation stops, default: 4.

        init_only : boolean
            Perform initialisation only, default: :code:`False`.

        init_previous : boolean
            Initialise the calculation with values from the previous
            calculation, default: :code:`True`.

        use_cuda : boolean
            Use cuda instead of numpy for matrix inversion, default:
            :code:`False`.

        Note
        ----
        For more information on the solution process have a look at the online
        documentation at tespy.readthedocs.io in the section "TESPy modules".
        """
        ## to own function
        self.new_design = False
        if self.design_path == design_path and design_path is not None:
            for c in self.conns['object']:
                if c.new_design:
                    self.new_design = True
                    break
            if not self.new_design:
                for cp in self.comps['object']:
                    if cp.new_design:
                        self.new_design = True
                        break

        else:
            self.new_design = True

        self.converged = False
        self.init_path = init_path
        self.design_path = design_path
        self.max_iter = max_iter
        self.min_iter = min_iter
        self.init_previous = init_previous
        self.iter = 0
        self.use_cuda = use_cuda

        if self.use_cuda and cu is None:
            msg = (
                'Specifying use_cuda=True requires cupy to be installed on '
                'your machine. Numpy will be used instead.'
            )
            logger.warning(msg)
            self.use_cuda = False

        if mode not in ['offdesign', 'design']:
            msg = 'Mode must be "design" or "offdesign".'
            logger.error(msg)
            raise ValueError(msg)
        else:
            self.mode = mode

        if not self.checked:
            self.check_network()

        msg = (
            "Solver properties:\n"
            f" - mode: {self.mode}\n"
            f" - init_path: {self.init_path}\n"
            f" - design_path: {self.design_path}\n"
            f" - min_iter: {self.min_iter}\n"
            f" - max_iter: {self.max_iter}\n"
            f" - init_path: {self.init_path}"
        )
        logger.debug(msg)

        msg = (
            "Network information:\n"
            f" - Number of components: {len(self.comps)}\n"
            f" - Number of connections: {len(self.conns)}\n"
            f" - Number of busses: {len(self.busses)}"
        )
        logger.debug(msg)

        self.initialise()

        for data in self.variables_dict.values():
            container = data["obj"]
            print(container._J_col, container._reference_container, container.J_col())

        if init_only:
            return

        msg = 'Starting solver.'
        logger.info(msg)

        self.solve_determination()

        self.solve_loop(print_results=print_results)

        if self.lin_dep:
            msg = (
                'Singularity in jacobian matrix, calculation aborted! Make '
                'sure your network does not have any linear dependencies in '
                'the parametrisation. Other reasons might be\n-> given '
                'temperature with given pressure in two phase region, try '
                'setting enthalpy instead or provide accurate starting value '
                'for pressure.\n-> given logarithmic temperature differences '
                'or kA-values for heat exchangers, \n-> support better '
                'starting values.\n-> bad starting value for fuel mass flow '
                'of combustion chamber, provide small (near to zero, but not '
                'zero) starting value.'
            )
            logger.error(msg)
            return

        if not self.progress:
            msg = (
                'The solver does not seem to make any progress, aborting '
                'calculation. Residual value is '
                '{:.2e}'.format(norm(self.residual)) + '. This frequently '
                'happens, if the solver pushes the fluid properties out of '
                'their feasible range.'
            )
            logger.warning(msg)
            return

        self.postprocessing()

        msg = 'Calculation complete.'
        logger.info(msg)
        return

    def solve_loop(self, print_results=True):
        r"""Loop of the newton algorithm."""
        # parameter definitions
        self.residual_history = np.array([])
        self.residual = np.zeros([self.num_vars])
        self.increment = np.ones([self.num_vars])
        self.jacobian = np.zeros((self.num_vars, self.num_vars))

        self.start_time = time()
        self.progress = True

        if self.iterinfo:
            self.iterinfo_head(print_results)

        for self.iter in range(self.max_iter):
            self.increment_filter = np.absolute(self.increment) < ERR ** 2
            self.solve_control()
            self.residual_history = np.append(
                self.residual_history, norm(self.residual)
            )

            if self.iterinfo:
                self.iterinfo_body(print_results)

            if (
                    (self.iter >= self.min_iter - 1
                     and (self.residual_history[-2:] < ERR ** 0.5).all())
                    or self.lin_dep
                ):
                self.converged = not self.lin_dep
                break

            if self.iter > 40:
                if (
                    all(
                        self.residual_history[(self.iter - 3):] >= self.residual_history[-3] * 0.95
                    ) and self.residual_history[-1] >= self.residual_history[-2] * 0.95
                ):
                    self.progress = False
                    break

        self.end_time = time()

        if self.iterinfo:
            self.iterinfo_tail(print_results)

        if self.iter == self.max_iter - 1:
            msg = (
                f"Reached maximum iteration count ({self.max_iter})), "
                "calculation stopped. Residual value is "
                "{:.2e}".format(norm(self.residual))
            )
            logger.warning(msg)

        return

    def solve_determination(self):
        r"""Check, if the number of supplied parameters is sufficient."""
        # number of user defined functions
        self.num_ude_eq = len(self.user_defined_eq)

        for func in self.user_defined_eq.values():
            # remap connection objects
            func.conns = [
                self.conns.loc[c.label, 'object'] for c in func.conns
            ]
            # remap jacobian
            func.jacobian = {}

        # total number of variables
        self.num_vars = (
            self.num_conn_vars + self.num_comp_vars
        )

        msg = f'Number of connection equations: {self.num_conn_eq}.'
        logger.debug(msg)
        msg = f'Number of bus equations: {self.num_bus_eq}.'
        logger.debug(msg)
        msg = f'Number of component equations: {self.num_comp_eq}.'
        logger.debug(msg)
        msg = f'Number of user defined equations: {self.num_ude_eq}.'
        logger.debug(msg)

        msg = f'Total number of variables: {self.num_vars}.'
        logger.debug(msg)
        msg = f'Number of component variables: {self.num_comp_vars}.'
        logger.debug(msg)
        msg = f"Number of connection variables: {self.num_conn_vars}."
        logger.debug(msg)

        n = (
            self.num_comp_eq + self.num_conn_eq +
            self.num_bus_eq + self.num_ude_eq
        )
        if n > self.num_vars:
            msg = (
                f"You have provided too many parameters: {self.num_vars} "
                f"required, {n} supplied. Aborting calculation!"
            )
            logger.error(msg)
            raise hlp.TESPyNetworkError(msg)
        elif n < self.num_vars:
            msg = (
                f"You have not provided enough parameters: {self.num_vars} "
                f"required, {n} supplied. Aborting calculation!"
            )
            logger.error(msg)
            raise hlp.TESPyNetworkError(msg)

    def iterinfo_head(self, print_results=True):
        """Print head of convergence progress."""
        # Start with defining the format here
        self.iterinfo_fmt = ' {iter:5s} | {residual:10s} | {progress:10s} '
        self.iterinfo_fmt += '| {massflow:10s} | {pressure:10s} | {enthalpy:10s} '
        self.iterinfo_fmt += '| {fluid:10s} | {component:10s} '
        # Use the format to create the first logging entry
        msg = self.iterinfo_fmt.format(
            iter='iter',
            residual='residual',
            progress='progress',
            massflow='massflow',
            pressure='pressure',
            enthalpy='enthalpy',
            fluid='fluid',
            component='component'
        )
        logger.progress(0, msg)
        msg2 = '-' * 7 + '+------------' * 7

        logger.progress(0, msg2)
        if print_results:
            print('\n' + msg + '\n' + msg2)
        return

    def iterinfo_body(self, print_results=True):
        """Print convergence progress."""
        m = [k for k, v in self.variables_dict.items() if v["variable"] == "m"]
        p = [k for k, v in self.variables_dict.items() if v["variable"] == "h"]
        p = [k for k, v in self.variables_dict.items() if v["variable"] == "p"]
        h = [k for k, v in self.variables_dict.items() if v["variable"] == "h"]
        fl = [k for k, v in self.variables_dict.items() if v["variable"] == "fluid"]
        cp = [k for k in self.variables_dict if k not in m + p + h + fl]

        iter_str = str(self.iter + 1)
        residual_norm = norm(self.residual)
        residual = 'NaN'
        progress = 'NaN'
        massflow = 'NaN'
        pressure = 'NaN'
        enthalpy = 'NaN'
        fluid = 'NaN'
        component = 'NaN'

        progress_val = -1

        if not np.isnan(residual_norm):
            residual = '{:.2e}'.format(residual_norm)

            if not self.lin_dep:
                massflow = '{:.2e}'.format(norm(self.increment[m]))
                pressure = '{:.2e}'.format(norm(self.increment[p]))
                enthalpy = '{:.2e}'.format(norm(self.increment[h]))
                fluid = '{:.2e}'.format(norm(self.increment[fl]))
                component  = '{:.2e}'.format(norm(self.increment[cp]))

            # This should not be hardcoded here.
            if residual_norm > np.finfo(float).eps * 100:
                progress_min = math.log(ERR)
                progress_max = math.log(ERR ** 0.5) * -1
                progress_val = math.log(max(residual_norm, ERR)) * -1
                # Scale to 0-1
                progres_scaled = (
                    (progress_val - progress_min)
                    / (progress_max - progress_min)
                )
                progress_val = max(0, min(1, progres_scaled))
                # Scale to 100%
                progress_val = int(progress_val * 100)
            else:
                progress_val = 100

            progress = '{:d} %'.format(progress_val)

        msg = self.iterinfo_fmt.format(
            iter=iter_str,
            residual=residual,
            progress=progress,
            massflow=massflow,
            pressure=pressure,
            enthalpy=enthalpy,
            fluid=fluid,
            component=component
        )
        logger.progress(progress_val, msg)
        if print_results:
            print(msg)
        return

    def iterinfo_tail(self, print_results=True):
        """Print tail of convergence progress."""
        num_iter = self.iter + 1
        clc_time = self.end_time - self.start_time
        num_ips = num_iter / clc_time if clc_time > 1e-10 else np.inf
        msg = '-' * 7 + '+------------' * 7
        logger.progress(100, msg)
        msg = (
            "Total iterations: {0:d}, Calculation time: {1:.2f} s, "
            "Iterations per second: {2:.2f}"
        ).format(num_iter, clc_time, num_ips)
        logger.debug(msg)
        if print_results:
            print(msg)
        return

    def matrix_inversion(self):
        """Invert matrix of derivatives and caluclate increment."""
        self.lin_dep = True
        try:
            # Let the matrix inversion be computed by the GPU if use_cuda in
            # global_vars.py is true.
            if self.use_cuda:
                self.increment = cu.asnumpy(cu.dot(
                    cu.linalg.inv(cu.asarray(self.jacobian)),
                    -cu.asarray(self.residual)
                ))
            else:
                self.increment = np.linalg.inv(
                    self.jacobian
                ).dot(-self.residual)
            self.lin_dep = False
        except np.linalg.LinAlgError:
            self.increment = self.residual * 0

    def update_variables(self):
        # cast dtype to float from numpy float64
        # this is necessary to keep the doctests running and note make them
        # look ugly all over the place
        # I have yet to come up with a better idea, or vectorize all operations
        # which requires major changes in tespy
        increment = [float(val) for val in self.increment]
        # the J_cols here point to actual variables, no need to call to
        # get_J_col yet
        for _, data in self.variables_dict.items():
            if data["variable"] in ["m", "h"]:
                container = data["obj"]
                container._val_SI += increment[container.J_col()]
            elif data["variable"] == "p":
                container = data["obj"]
                relax = max(
                    1, -2 * increment[container.J_col()] / container.val_SI()
                )
                container._val_SI += increment[container.J_col()] / relax
            elif data["variable"] == "fluid":
                container = data["obj"].fluid
                container.val[data["fluid"]] += increment[
                    container.J_col()[data["fluid"]]
                ]

                if container.val[data["fluid"]] < ERR :
                    container.val[data["fluid"]] = 0
                elif container.val[data["fluid"]] > 1 - ERR :
                    container.val[data["fluid"]] = 1
            else:
                # add increment
                data["obj"].val += increment[data["obj"].J_col()]

                # keep value within specified value range
                if data["obj"].val < data["obj"].min_val:
                    data["obj"].val = data["obj"].min_val
                elif data["obj"].val > data["obj"].max_val:
                    data["obj"].val = data["obj"].max_val

    def check_variable_bounds(self):

        for c in self.conns['object']:
            # check the fluid properties for physical ranges
            if len(c.fluid.is_var) > 0:
                total_mass_fractions = sum(c.fluid.val.values())
                for fluid in c.fluid.is_var:
                    c.fluid.val[fluid] /= total_mass_fractions

            c.build_fluid_data()
            self.check_connection_properties(c)

        # second property check for first three iterations without an init_file
        if self.iter < 3:
            for cp in self.comps['object']:
                cp.convergence_check()

            for c in self.conns['object']:
                self.check_connection_properties(c)

    def solve_control(self):
        r"""
        Control iteration step of the newton algorithm.

        - Calculate the residual value for each equation
        - Calculate the jacobian matrix
        - Calculate new values for variables
        - Restrict fluid properties to value ranges
        - Check component parameters for consistency
        """
        self.solve_components()
        self.solve_busses()
        self.solve_connections()
        self.solve_user_defined_eq()
        self.matrix_inversion()

        # check for linear dependency
        if self.lin_dep:
            return

        self.update_variables()
        # self.check_variable_bounds()

    def check_connection_properties(self, c):
        r"""
        Check for invalid fluid property values.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            Connection to check fluid properties.
        """
        fl = fp.single_fluid(c.fluid_data)

        # pure fluid
        if fl is not None:
            # pressure
            if c.p.is_var():
                c.check_pressure_bounds(fl)

            # enthalpy
            if c.h.is_var():
                c.check_enthalpy_bounds(fl)

                # two-phase related
                if (c.Td_bp.is_set or c.state.is_set) and self.iter < 3:
                    c.check_two_phase_bounds(fl)

        # mixture
        elif self.iter < 4 and not c.good_starting_values:
            # pressure
            if c.p.is_var:
                if c.p.val_SI <= self.p_range_SI[0]:
                    c.p.val_SI = self.p_range_SI[0]
                    logger.debug(c._property_range_message('p'))

                elif c.p.val_SI >= self.p_range_SI[1]:
                    c.p.val_SI = self.p_range_SI[1]
                    logger.debug(c._property_range_message('p'))

            # enthalpy
            if c.h.is_var:
                if c.h.val_SI < self.h_range_SI[0]:
                    c.h.val_SI = self.h_range_SI[0]
                    logger.debug(c._property_range_message('h'))

                elif c.h.val_SI > self.h_range_SI[1]:
                    c.h.val_SI = self.h_range_SI[1]
                    logger.debug(c._property_range_message('h'))

                # temperature
                if c.T.is_set:
                    c.check_temperature_bounds()

        # mass flow
        if c.m.val_SI <= self.m_range_SI[0] and c.m.is_var:
            c.m.val_SI = self.m_range_SI[0]
            logger.debug(c._property_range_message('m'))

        elif c.m.val_SI >= self.m_range_SI[1] and c.m.is_var:
            c.m.val_SI = self.m_range_SI[1]
            logger.debug(c._property_range_message('m'))

    def solve_components(self):
        r"""
        Calculate the residual and derivatives of component equations.
        """
        # fetch component equation residuals and component partial derivatives
        sum_eq = 0
        for cp in self.comps['object']:
            print(cp.label)
            cp.solve(self.increment_filter)
            self.residual[sum_eq:sum_eq + cp.num_eq] = cp.residual

            if len(cp.jacobian) > 0:
                rows = [k[0] + sum_eq for k in cp.jacobian]
                columns = [k[1] for k in cp.jacobian]
                data = list(cp.jacobian.values())
                self.jacobian[rows, columns] = data
                sum_eq += cp.num_eq

            cp.it += 1

    def solve_connections(self):
        r"""
        Calculate the residual and derivatives of connection equations.
        """
        sum_eq = self.num_comp_eq
        for c in self.conns['object']:
            c.solve(self.increment_filter)
            self.residual[sum_eq:sum_eq + c.num_eq] = c.residual

            if len(c.jacobian) > 0:
                rows = [k[0] + sum_eq for k in c.jacobian]
                columns = [k[1] for k in c.jacobian]
                data = list(c.jacobian.values())
                self.jacobian[rows, columns] = data
                sum_eq += c.num_eq

            c.it += 1

    def solve_user_defined_eq(self):
        """
        Calculate the residual and jacobian of user defined equations.
        """
        sum_eq = self.num_comp_eq + self.num_conn_eq + self.num_bus_eq
        for ude in self.user_defined_eq.values():
            ude.solve()
            self.residual[sum_eq] = ude.residual

            if len(ude.jacobian) > 0:
                columns = [k for k in ude.jacobian]
                data = list(ude.jacobian.values())
                self.jacobian[sum_eq, columns] = data
                sum_eq += 1

    def solve_busses(self):
        r"""
        Calculate the equations and the partial derivatives for the busses.
        """
        sum_eq = self.num_comp_eq + self.num_conn_eq
        for bus in self.busses.values():
            if bus.P.is_set:

                bus.solve()
                self.residual[sum_eq] = bus.residual

                if len(bus.jacobian) > 0:
                    columns = [k for k in bus.jacobian]
                    data = list(bus.jacobian.values())
                    self.jacobian[sum_eq, columns] = data

                bus.clear_jacobian()
                sum_eq += 1

    def postprocessing(self):
        r"""Calculate connection, bus and component parameters."""
        self.process_connections()
        self.process_components()
        self.process_busses()

        msg = 'Postprocessing complete.'
        logger.info(msg)

    def process_connections(self):
        """Process the Connection results."""
        for c in self.conns['object']:
            c.good_starting_values = True
            c.calc_results()

            self.results['Connection'].loc[c.label] = (
                [
                    _ for key in fpd.keys()
                    for _ in [c.get_attr(key).val, c.get_attr(key).unit]
                ] + [
                    c.fluid.val[fluid] if fluid in c.fluid.val else np.nan
                    for fluid in self.all_fluids
                ] + [
                    c.phase.val
                ]

            )

    def process_components(self):
        """Process the component results."""
        # components
        for cp in self.comps['object']:
            cp.calc_parameters()
            cp.check_parameter_bounds()

            key = cp.__class__.__name__
            for param in self.results[key].columns:
                p = cp.get_attr(param)
                if (p.func is not None or (p.func is None and p.is_set) or
                        p.is_result):
                    self.results[key].loc[cp.label, param] = p.val
                else:
                    self.results[key].loc[cp.label, param] = np.nan

    def process_busses(self):
        """Process the bus results."""
        # busses
        for b in self.busses.values():
            for cp in b.comps.index:
                # get components bus func value
                bus_val = cp.calc_bus_value(b)
                eff = cp.calc_bus_efficiency(b)
                cmp_val = cp.bus_func(b.comps.loc[cp])

                b.comps.loc[cp, 'char'].get_domain_errors(
                    cp.calc_bus_expr(b), cp.label)

                # save as reference value
                if self.mode == 'design':
                    if b.comps.loc[cp, 'base'] == 'component':
                        design_value = cmp_val
                    else:
                        design_value = bus_val

                    b.comps.loc[cp, 'P_ref'] = design_value

                else:
                    design_value = b.comps.loc[cp, 'P_ref']

                result = [cmp_val, bus_val, eff, design_value]
                self.results[b.label].loc[cp.label] = result

            b.P.val = float(self.results[b.label]['bus value'].sum())

    def print_results(self, colored=True, colors=None, print_results=True):
        r"""Print the calculations results to prompt."""
        # Define colors for highlighting values in result table
        if colors is None:
            colors = {}
        result = ""
        coloring = {
            'end': '\033[0m',
            'set': '\033[94m',
            'err': '\033[31m',
            'var': '\033[32m'
        }
        coloring.update(colors)

        if not hasattr(self, 'results'):
            msg = (
                'It is not possible to print the results of a network, that '
                'has never been solved successfully. Results DataFrames are '
                'only available after a full simulation run is performed.')
            raise hlp.TESPyNetworkError(msg)

        for cp in self.comps['comp_type'].unique():
            df = self.results[cp].copy()

            # are there any parameters to print?
            if df.size > 0:
                cols = df.columns
                if len(cols) > 0:
                    for col in cols:
                        df[col] = df.apply(
                            self.print_components, axis=1,
                            args=(col, colored, coloring))

                    df.dropna(how='all', inplace=True)

                    if len(df) > 0:
                        # printout with tabulate
                        result += f"\n##### RESULTS ({cp}) #####\n"
                        result += (
                            tabulate(
                                df, headers='keys', tablefmt='psql',
                                floatfmt='.2e'
                            )
                        )

        # connection properties
        df = self.results['Connection'].loc[:, ['m', 'p', 'h', 'T', 'x', 'phase']].copy()
        df = df.astype(str)
        for c in df.index:
            if not self.get_conn(c).printout:
                df.drop([c], axis=0, inplace=True)

            elif colored:
                conn = self.get_conn(c)
                for col in df.columns:
                    if conn.get_attr(col).is_set:
                        df.loc[c, col] = (
                            coloring['set'] + str(conn.get_attr(col).val) +
                            coloring['end'])

        if len(df) > 0:
            result += ('\n##### RESULTS (Connection) #####\n')
            result += (
                tabulate(df, headers='keys', tablefmt='psql', floatfmt='.3e')
            )

        for b in self.busses.values():
            if b.printout:
                df = self.results[b.label].loc[
                    :, ['component value', 'bus value', 'efficiency']
                ].copy()
                df.loc['total'] = df.sum()
                df.loc['total', 'efficiency'] = np.nan
                if colored:
                    df["bus value"] = df["bus value"].astype(str)
                    if b.P.is_set:
                        df.loc['total', 'bus value'] = (
                            coloring['set'] + str(df.loc['total', 'bus value']) +
                            coloring['end']
                        )
                result += f"\n##### RESULTS (Bus: {b.label}) #####\n"
                result += (
                    tabulate(
                        df, headers='keys', tablefmt='psql',
                        floatfmt='.3e'
                    )
                )
        if len(str(result)) > 0:
            logger.result(result)
            if print_results:
                print(result)
        return

    def print_components(self, c, *args):
        """
        Get the print values for the component data.

        Parameters
        ----------
        c : pandas.core.series.Series
            Series containing the component data.

        param : str
            Component parameter to print.

        colored : bool
            Color the printout.

        coloring : dict
            Coloring information for colored printout.

        Returns
        ----------
        value : str
            String representation of the value to print.
        """
        param, colored, coloring = args
        comp = self.get_comp(c.name)
        if comp.printout:
            # select parameter from results DataFrame
            val = c[param]
            if not colored:
                return str(val)
            # else part
            if (val < comp.get_attr(param).min_val - ERR or
                    val > comp.get_attr(param).max_val + ERR ):
                return f"{coloring['err']} {val} {coloring['end']}"
            if comp.get_attr(args[0]).is_var:
                return f"{coloring['var']} {val} {coloring['end']}"
            if comp.get_attr(args[0]).is_set:
                return f"{coloring['set']} {val} {coloring['end']}"
            return str(val)
        else:
            return np.nan

    def export(self, path):
        """Export the network structure and parametrization."""
        path, path_comps = self._create_export_paths(path)
        self.export_network(path)
        self.export_connections(path)
        self.export_components(path_comps)
        self.export_busses(path)

    def save(self, path):
        r"""
        Save the results to results files.

        Parameters
        ----------
        filename : str
            Path for the results.

        Note
        ----
        Results will be saved to path. The results contain:

        - network.json (network information)
        - connections.csv (connection information)
        - folder components containing .csv files for busses and
          characteristics as well as .csv files for all types of components
          within your network.
        """
        path, path_comps = self._create_export_paths(path)

        # save relevant design point information
        self.save_connections(path)
        self.save_components(path_comps)
        self.save_busses(path)

    def _create_export_paths(self, path):
        logger.debug('Saving network to path %s.', path)
        # creat path, if non existent
        os.makedirs(path, exist_ok=True)

        # create path for component folder if non existent
        path_comps = os.path.join(path, "components")
        os.makedirs(path_comps, exist_ok=True)

        return path, path_comps

    def export_network(self, fn):
        r"""
        Save basic network configuration.

        Parameters
        ----------
        fn : str
            Path/filename for the network configuration file.
        """
        with open(os.path.join(fn, 'network.json'), 'w') as f:
            json.dump(self._serialize(), f, indent=4)

        logger.debug('Network information saved to %s.', fn)

    def save_connections(self, fn):
        r"""
        Save the connection properties.

        Parameters
        ----------
        fn : str
            Path/filename for the file.
        """
        self.results["Connection"].to_csv(
            os.path.join(fn, "connections.csv"), sep=';', decimal='.', index=True, na_rep='nan'
        )
        logger.debug('Connection information saved to %s.', fn)

    def save_components(self, path):
        r"""
        Save the component properties.

        Parameters
        ----------
        path : str
            Path/filename for the file.
        """
        for c in self.comps['comp_type'].unique():
            fn = os.path.join(path, f"{c}.csv")
            self.results[c].to_csv(fn, sep=';', decimal='.', index=True, na_rep='nan')
            logger.debug('Component information (%s) saved to %s.', c, fn)

    def save_busses(self, fn):
        r"""
        Save the bus properties.

        Parameters
        ----------
        fn : str
            Path/filename for the file.
        """
        if len(self.busses) > 0:
            bus_data = {}
            for label, bus in self.busses.items():
                bus_data[label] = self.results[label]["design value"].to_dict()
            fn = os.path.join(fn, "busses.json")
            with open(fn, "w", encoding="utf-8") as f:
                json.dump(bus_data, f, indent=4)
            logger.debug('Bus information saved to %s.', fn)

    def export_connections(self, fn):
        connections = {}
        for c in self.conns["object"]:
            connections.update(c._serialize())

        fn = os.path.join(fn, "connections.json")
        with open(fn, "w", encoding="utf-8") as f:
            json.dump(connections, f, indent=4)
        logger.debug('Connection information exported to %s.', fn)

    def export_components(self, fn):
        for c in self.comps["comp_type"].unique():
            components = {}
            for cp in self.comps.loc[self.comps["comp_type"] == c, "object"]:
                components.update(cp._serialize())

            fname = os.path.join(fn, f"{c}.json")
            with open(fname, "w", encoding="utf-8") as f:
                json.dump(components, f, indent=4)
            logger.debug('Component information exported to %s.', fname)

    def export_busses(self, fn):
        if len(self.busses) > 0:
            busses = {}
            for bus in self.busses.values():
                busses.update(bus._serialize())
            fn = os.path.join(fn, 'busses.json')
            with open(fn, "w", encoding="utf-8") as f:
                json.dump(busses, f, indent=4)
            logger.debug('Bus information exported to %s.', fn)
