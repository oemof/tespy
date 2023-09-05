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
import ast
import json
import os
from collections import OrderedDict
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
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
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

    def __init__(self, fluids=None, **kwargs):
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
            "fluids": list(self.all_fluids),
        }

    def set_defaults(self):
        """Set default network properties."""
        # connection dataframe
        self.conns = pd.DataFrame(
            columns=['object', 'source', 'source_id', 'target', 'target_id'],
            dtype='object'
        )
        self.all_fluids = set()
        # component dataframe
        self.comps = pd.DataFrame(dtype='object')
        # user defined function dictionary for fast access
        self.user_defined_eq = {}
        # bus dictionary
        self.busses = OrderedDict()
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
        self.m_range_SI = np.array([-1e12, 1e12])
        self.p_range_SI = np.array([2e2, 300e5])
        self.h_range_SI = np.array([1e3, 7e6])

        for prop in ['m', 'p', 'h']:
            limits = self.get_attr(prop + '_range_SI')
            msg = (
                'Default ' + fpd[prop]['text'] + ' limits\n'
                'min: ' + str(limits[0]) + ' ' +
                self.get_attr(prop + '_unit') + '\n'
                'max: ' + str(limits[1]) + ' ' + self.get_attr(prop + '_unit'))
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
            unit = prop + '_unit'
            if unit in kwargs:
                if kwargs[unit] in fpd[prop]['units']:
                    self.__dict__.update({unit: kwargs[unit]})
                    msg = (
                        'Setting ' + fpd[prop]['text'] +
                        ' unit: ' + kwargs[unit] + '.')
                    logger.debug(msg)
                else:
                    keys = ', '.join(fpd[prop]['units'].keys())
                    msg = (
                        'Allowed units for ' +
                        fpd[prop]['text'] + ' are: ' + keys)
                    logger.error(msg)
                    raise ValueError(msg)

        for prop in ['m', 'p', 'h']:
            if prop + '_range' in kwargs:
                if isinstance(kwargs[prop + '_range'], list):
                    self.__dict__.update(
                        {prop + '_range_SI': hlp.convert_to_SI(
                            prop, np.array(kwargs[prop + '_range']),
                            self.get_attr(prop + '_unit'))})
                else:
                    msg = (
                        'Specify the value range as list: [' + prop +
                        '_min, ' + prop + '_max]')
                    logger.error(msg)
                    raise TypeError(msg)

                limits = self.get_attr(prop + '_range_SI')
                msg = (
                    'Setting ' + fpd[prop]['text'] +
                    ' limits\nmin: ' + str(limits[0]) + ' ' +
                    fpd[prop]['SI_unit'] + '\n'
                    'max: ' + str(limits[1]) + ' ' +
                    fpd[prop]['SI_unit'])
                logger.debug(msg)

        # update non SI value ranges
        for prop in ['m', 'p', 'h']:
            self.__dict__.update({
                prop + '_range': hlp.convert_from_SI(
                    prop, self.get_attr(prop + '_range_SI'),
                    self.get_attr(prop + '_unit')
                )
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
            msg = 'Network has no attribute \"' + str(key) + '\".'
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
            logger.warning('Connection with label %s not found.', label)
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
            logger.warning('Component with label %s not found.', label)
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
                msg = ('Must provide tespy.connections.connection.Connection '
                       'objects as parameters.')
                logger.error(msg)
                raise TypeError(msg)

            elif c.label in self.conns.index:
                msg = (
                    'There is already a connection with the label ' +
                    c.label + '. The connection labels must be unique!')
                logger.error(msg)
                raise ValueError(msg)

            c.good_starting_values = False

            self.conns.loc[c.label] = [
                c, c.source, c.source_id, c.target, c.target_id
            ]

            msg = 'Added connection ' + c.label + ' to network.'
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
            msg = ('Deleted connection ' + c.label + ' from network.')
            logger.debug(msg)

        self._del_comps(comps)

        # set status "checked" to false, if connection is deleted from network.
        self.checked = False

    def check_conns(self):
        r"""Check connections for multiple usage of inlets or outlets."""
        dub = self.conns.loc[self.conns.duplicated(["source", "source_id"])]
        for c in dub['object']:
            targets = []
            for conns in self.conns.loc[
                    (self.conns["source"].values == c.source) &
                    (self.conns["source_id"].values == c.source_id),
                    "object"]:
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
            for conns in self.conns.loc[
                    (self.conns["target"].values == c.target) &
                    (self.conns["target_id"].values == c.target_id),
                    "object"]:
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
                msg = f"Deleted component {comp.label} from network."
                logger.debug(msg)

    def add_ude(self, *args):
        r"""
        Add a user defined function to the network.

        Parameters
        ----------
        c : tespy.tools.helpers.UserDefinedEquation
            The objects to be added to the network, UserDefinedEquation objects
            ci :code:`del_conns(c1, c2, c3, ...)`.
        """
        for c in args:
            if not isinstance(c, hlp.UserDefinedEquation):
                msg = ('Must provide tespy.connections.connection.Connection '
                       'objects as parameters.')
                logger.error(msg)
                raise TypeError(msg)

            elif c.label in self.user_defined_eq:
                msg = (
                    'There is already a UserDefinedEquation with the label ' +
                    c.label + '. The UserDefinedEquation labels must be '
                    'unique within a network')
                logger.error(msg)
                raise ValueError(msg)

            self.user_defined_eq[c.label] = c
            msg = 'Added UserDefinedEquation ' + c.label + ' to network.'
            logger.debug(msg)

    def del_ude(self, *args):
        """
        Remove a user defined function from the network.

        Parameters
        ----------
        c : tespy.tools.helpers.UserDefinedEquation
            The objects to be added deleted from the network,
            UserDefinedEquation objects ci :code:`del_conns(c1, c2, c3, ...)`.
        """
        for c in args:
            del self.user_defined_eq[c.label]
            msg = 'Deleted UserDefinedEquation ' + c.label + ' from network.'
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
                    dtype='float64')

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
        self.create_branches()

        # network checked
        self.checked = True
        msg = 'Networkcheck successful.'
        logger.info(msg)

    def create_branches(self):
        self.branches = {}
        mask = self.comps["object"].apply(lambda c: c.is_branch_source())
        start_components = self.comps["object"].loc[mask]
        if len(start_components) == 0:
            msg = (
                "You cannot build a system without at least one CycleCloser or "
                "a Source and Sink."
            )
            raise hlp.TESPyNetworkError(msg)

        for start in start_components:
            self.branches.update(start.start_branch())

        self.fluid_wrapper_branches = {}
        mask = self.comps["comp_type"].isin(["Source", "CycleCloser"])
        start_components = self.comps["object"].loc[mask]

        for start in start_components:
            self.fluid_wrapper_branches.update(start.start_fluid_wrapper_branch())

        merged_fluid_wrapper_branches = self.fluid_wrapper_branches.copy()
        for branch_name, branch_data in self.fluid_wrapper_branches.items():
            if branch_name not in merged_fluid_wrapper_branches:
                continue
            merges = [cp for cp in branch_data["components"] if cp.component() == "merge" or cp.component() == "combustion chamber"]
            if any(merges):
                for ob_name, ob_data in self.fluid_wrapper_branches.copy().items():
                    if ob_name != branch_name:
                        for merge in merges:
                            if merge in ob_data["components"]:
                                merged_fluid_wrapper_branches[branch_name]["connections"] = list(set(branch_data["connections"] + ob_data["connections"]))
                                merged_fluid_wrapper_branches[branch_name]["components"] = list(set(branch_data["components"] + ob_data["components"]))
                                del merged_fluid_wrapper_branches[ob_name]
                                break

        self.fluid_wrapper_branches = merged_fluid_wrapper_branches

        self.massflow_branches = hlp.get_all_subdictionaries(self.branches)

        self.fluid_branches = {}
        for branch_name, branch_data in self.branches.items():
            subbranches = hlp.get_all_subdictionaries(branch_data["subbranches"])
            main = {k: v for k, v in branch_data.items() if k != "subbranches"}
            self.fluid_branches[branch_name] = [main] + subbranches

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
                    columns=cols, dtype='float64')
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

        self.propagate_fluid_wrappers()
        self.presolve_massflow_topology()
        self.presolve_fluid_topology()

        self.init_set_properties()

        if self.mode == 'offdesign':
            self.redesign = True
            if self.design_path is None:
                # must provide design_path
                msg = ('Please provide "design_path" for every offdesign '
                       'calculation.')
                logger.error(msg)
                raise hlp.TESPyNetworkError(msg)

            # load design case
            if self.new_design:
                self.init_offdesign_params()

            self.init_offdesign()

        else:
            # reset any preceding offdesign calculation
            self.init_design()
            # generic fluid initialisation
            # for offdesign cases good starting values should be available
            self.init_fluids()

        # generic fluid property initialisation
        self.init_properties()

        msg = 'Network initialised.'
        logger.info(msg)

    def propagate_fluid_wrappers(self):

        for branch_data in self.fluid_wrapper_branches.values():
            all_connections = [c for c in branch_data["connections"]]
            any_fluids_set = [f for c in all_connections for f in c.fluid.is_set]
            fluid_set_wrappers = {f: w for c in all_connections for f, w in c.fluid.wrapper.items() if f in c.fluid.is_set}
            mixing_rules = [c.mixing_rule for c in branch_data["connections"] if c.mixing_rule is not None]
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
            for c in all_connections:
                c.mixing_rule = list(mixing_rule)[0]
                c._potential_fluids = potential_fluids
                for f in potential_fluids:
                    if (f not in c.fluid.is_set and f not in c.fluid.val and f not in c.fluid.val0):
                        c.fluid.val[f] = 0
                    elif f not in c.fluid.is_set and f in c.fluid.val0:
                            c.fluid.val[f] = c.fluid.val0[f]
                    if f not in c.fluid.wrapper and f in fluid_set_wrappers:
                        c.fluid.wrapper[f] = fluid_set_wrappers[f]
                    elif f not in c.fluid.wrapper:
                        c._create_fluid_wrapper(f, fp.CoolPropWrapper, "HEOS")

    def presolve_massflow_topology(self):

        # mass flow is a single variable in each sub branch
        # fluid composition is a single variable in each main branch
        for branch in self.massflow_branches:

            num_massflow_specs = 0
            for c in branch["connections"]:
                # number of specifications cannot exceed 1
                num_massflow_specs += c.m.is_set

                if c.m.is_set:
                    main_conn = c

                # self reference is not allowed
                if c.m_ref.is_set:
                    if c.m_ref.val.obj in branch["connections"]:
                        msg = (
                            "You cannot reference a mass flow in the same "
                            f"linear branch. The connection {c.label} "
                            "references the connection "
                            f"{c.m_ref.val.obj.label}."
                        )
                        raise hlp.TESPyNetworkError(msg)

            if num_massflow_specs == 1:
                # set every mass flow in branch to the specified value
                for c in branch["connections"]:
                    # map all connection's mass flow data containers to first
                    # branch element
                    c._m_tmp = c.m
                    c.m = main_conn.m

                msg = (
                    "Removing "
                    f"{len(branch['connections']) - num_massflow_specs} "
                    "mass flow variables from system variables."
                )
                logger.debug(msg)
            elif num_massflow_specs > 1:
                msg = (
                    "You cannot specify two or more values for mass flow in "
                    "the same linear branch (starting at "
                    f"{branch['components'][0].label} and ending at "
                    f"{branch['components'][-1].label})."
                )
                raise hlp.TESPyNetworkError(msg)

            else:
                main_conn = branch["connections"][0]
                for c in branch["connections"][1:]:
                    # map all connection's mass flow data containers to first
                    # branch element
                    c._m_tmp = c.m
                    c.m = main_conn.m

    def presolve_fluid_topology(self):

        for branch_name, branch in self.fluid_branches.items():
            all_connections = [
                c for b in branch for c in b["connections"]
            ]
            main_conn = all_connections[0]
            fluid_specs = [f for c in all_connections for f in c.fluid.is_set]
            if len(fluid_specs) == 0:
                main_conn._fluid_tmp = dc_flu()
                main_conn._fluid_tmp.val = main_conn.fluid.val.copy()
                main_conn._fluid_tmp.is_set = main_conn.fluid.is_set.copy()
                main_conn._fluid_tmp.is_var = main_conn.fluid.is_var.copy()
                main_conn._fluid_tmp.wrapper = main_conn.fluid.wrapper.copy()

                for c in all_connections[1:]:
                    c._fluid_tmp = c.fluid
                    c.fluid = main_conn.fluid

                if len(main_conn._potential_fluids) > 1:
                    main_conn.fluid.is_var = {f for f in main_conn.fluid.val}
                else:
                    main_conn.fluid.val[list(main_conn._potential_fluids)[0]] = 1

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
                    all_fluids = main_conn.fluid.val.keys()
                    num_remaining_fluids = len(all_fluids) - len(fixed_fractions)
                    if num_remaining_fluids == 1:
                        missing_fluid = list(
                            main_conn.fluid.val.keys() - fixed_fractions.keys()
                        )[0]
                        fixed_fractions[missing_fluid] = 1 - mass_fraction_sum
                        variable = set()
                    else:
                        missing_fluids = (
                            main_conn.fluid.val.keys() - fixed_fractions.keys()
                        )
                        variable = {f for f in missing_fluids}

                else:
                    # fluid mass fraction is 100 %, all other fluids are 0 %
                    all_fluids = main_conn.fluid.val.keys()
                    remaining_fluids = (
                        main_conn.fluid.val.keys() - fixed_fractions.keys()
                    )
                    for f in remaining_fluids:
                        fixed_fractions[f] = 0

                    variable = set()

                main_conn._fluid_tmp = dc_flu()
                main_conn._fluid_tmp.val = main_conn.fluid.val.copy()
                main_conn._fluid_tmp.is_set = main_conn.fluid.is_set.copy()
                main_conn._fluid_tmp.is_var = main_conn.fluid.is_var.copy()
                main_conn._fluid_tmp.wrapper = main_conn.fluid.wrapper.copy()

                for c in all_connections[1:]:
                    c._fluid_tmp = c.fluid
                    c.fluid = main_conn.fluid

                main_conn.fluid.val.update(fixed_fractions)
                main_conn.fluid.is_set = {f: False for f in variable}
                main_conn.fluid.is_set.update({f: True for f in fixed_fractions})
                main_conn.fluid.is_var = variable

            [c.build_fluid_data() for c in all_connections]
            for fluid in main_conn.fluid.is_var:
                main_conn.fluid.J_col[fluid] = self.num_conn_vars
                self.variables_dict[self.num_conn_vars] = {
                    "obj": main_conn, "variable": "fluid", "fluid": fluid
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

            for key in ['m', 'p', 'h', 'T', 'x', 'v', 'Td_bp', 'vol', 's']:
                # read unit specifications
                if key == 'Td_bp':
                    c.get_attr(key).unit = self.get_attr('T_unit')
                else:
                    c.get_attr(key).unit = self.get_attr(key + '_unit')
                # set SI value
                if c.get_attr(key).is_set:
                    c.get_attr(key).val_SI = hlp.convert_to_SI(
                        key, c.get_attr(key).val, c.get_attr(key).unit)
                if c.get_attr(key).ref_set:
                    if key == 'T':
                        c.get_attr(key).ref.delta_SI = hlp.convert_to_SI(
                            'Td_bp', c.get_attr(key).ref.delta,
                            c.get_attr(key).unit)
                    else:
                        c.get_attr(key).ref.delta_SI = hlp.convert_to_SI(
                            key, c.get_attr(key).ref.delta,
                            c.get_attr(key).unit)

        if len(self.all_fluids) == 0:
            msg = (
                'Network has no fluids, please specify a list with fluids on '
                'network creation.'
            )
            logger.error(msg)
            raise hlp.TESPyNetworkError(msg)

        # set up results dataframe for connections
        # this should be done based on the connections
        properties = list(fpd.keys())
        self.all_fluids = set(self.all_fluids)
        cols = (
            [col for prop in properties for col in [prop, f"{prop}_unit"]]
            + list(self.all_fluids)
        )
        self.results['Connection'] = pd.DataFrame(columns=cols, dtype='float64')
        # include column for fluid balance in specs dataframe
        self.specifications['Connection'] = pd.DataFrame(columns=cols + ['balance'], dtype='bool')
        self.specifications['Ref'] = pd.DataFrame(columns=cols, dtype='bool')

        msg = (
            'Updated fluid property SI values and fluid mass fraction for '
            'user specified connection parameters.')
        logger.debug(msg)

    def _assign_variable_space(self, c):
        for key in ["m", "p", "h"]:
            variable = c.get_attr(key)
            if variable.is_var and variable not in self._conn_variables:
                variable.J_col = self.num_conn_vars
                self.variables_dict[self.num_conn_vars] = {
                    "obj": c, "variable": key
                }
                self._conn_variables += [variable]
                self.num_conn_vars += 1

        c.preprocess()

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
        for c in self.conns['object']:
            # read design point information of connections with
            # local_offdesign activated from their respective design path
            if c.local_offdesign:
                if c.design_path is None:
                    msg = (
                        'The parameter local_offdesign is True for the '
                        'connection ' + c.label + ', an individual '
                        'design_path must be specified in this case!')
                    logger.error(msg)
                    raise hlp.TESPyNetworkError(msg)

                # unset design parameters
                for var in c.design:
                    c.get_attr(var).is_set = False
                # set offdesign parameters
                for var in c.offdesign:
                    c.get_attr(var).is_set = True

                # read design point information
                df = self.init_read_connections(c.design_path)
                msg = (
                    'Reading individual design point information for '
                    'connection ' + c.label + ' from path ' + c.design_path +
                    'connections.')
                logger.debug(msg)

                # write data to connections
                self.init_conn_design_params(c, df)

            else:
                # unset all design values
                c.m.design = np.nan
                c.p.design = np.nan
                c.h.design = np.nan
                c.fluid.design = OrderedDict()

                c.new_design = True

                # switch connections to design mode
                if self.redesign:
                    for var in c.design:
                        c.get_attr(var).is_set = True

                    for var in c.offdesign:
                        c.get_attr(var).is_set = False

            self._assign_variable_space(c)

        # unset design values for busses, count bus equations and
        # reindex bus dictionary
        for b in self.busses.values():
            self.busses[b.label] = b
            self.num_bus_eq += b.P.is_set * 1
            for cp in b.comps.index:
                b.comps.loc[cp, 'P_ref'] = np.nan

        series = pd.Series(dtype='float64')
        for cp in self.comps['object']:
            # read design point information of components with
            # local_offdesign activated from their respective design path
            if cp.local_offdesign:
                if cp.design_path is not None:
                    # get type of component (class name)
                    c = cp.__class__.__name__
                    # read design point information
                    path = hlp.modify_path_os(
                        cp.design_path + '/components/' + c + '.csv')
                    df = pd.read_csv(
                        path, sep=';', decimal='.', converters={
                            'busses': ast.literal_eval,
                            'bus_P_ref': ast.literal_eval})
                    df.set_index('label', inplace=True)
                    # write data
                    self.init_comp_design_params(cp, df.loc[cp.label])

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
                        msg += var + ', '

                if switched:
                    msg = (msg[:-2] + ' to design value at component ' +
                           cp.label + '.')
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

            # component initialisation
            cp.preprocess(self.num_conn_vars + self.num_comp_vars)

            ct = cp.__class__.__name__
            for spec in self.specifications[ct].keys():
                if len(cp.get_attr(self.specifications['lookup'][spec])) > 0:
                    self.specifications[ct][spec].loc[cp.label] = (
                        cp.get_attr(self.specifications['lookup'][spec]))

            # count number of component equations and variables
            i = self.num_conn_vars + self.num_comp_vars
            for container, name in cp.vars.items():
                self.variables_dict[i] = {"obj": container, "variable": name}
                i += 1
            self.num_comp_vars += cp.num_vars
            self.num_comp_eq += cp.num_eq

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
            path = hlp.modify_path_os(f"{self.design_path}/components/{c}.csv")
            msg = (
                f"Reading design point information for components of type {c}"
                f"from path {path}."
            )
            logger.debug(msg)

            # read data
            df = pd.read_csv(path, sep=';', decimal='.', index_col=0)
            # iter through all components of this type and set data
            for c_label in df.index:
                comp = df_comps.loc[c_label, 'object']
                # read data of components with individual design_path
                if comp.design_path is not None:
                    path_c = hlp.modify_path_os(
                        comp.design_path + '/components/' + c + '.csv')
                    df_c = pd.read_csv(
                        path_c, sep=';', decimal='.', converters={
                             'busses': ast.literal_eval,
                             'bus_P_ref': ast.literal_eval})
                    df_c.set_index('label', inplace=True)
                    data = df_c.loc[comp.label]

                else:
                    data = df.loc[comp.label]

                # write data to components
                self.init_comp_design_params(comp, data)

        msg = 'Done reading design point information for components.'
        logger.debug(msg)

        # read connection design point information
        df = self.init_read_connections(self.design_path)
        msg = (
            'Reading design point information for connections from path ' +
            self.design_path + '/connections.csv.')
        logger.debug(msg)

        # iter through connections
        for c in self.conns['object']:

            # read data of connections with individual design_path
            if c.design_path is not None:
                df_c = self.init_read_connections(c.design_path)
                msg = (
                    'Reading individual design point information for '
                    'connection ' + c.label + ' from path ' + c.design_path +
                    '/connections.csv.')
                logger.debug(msg)

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
        # write design values to busses
        i = 0
        for b in data.busses:
            bus = self.busses[b].comps
            bus.loc[component, 'P_ref'] = data['bus_P_ref'][i]
            i += 1

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
        if c.label not in df.index:
            # no matches in the connections of the network and the design files
            msg = (
                'Could not find connection %s in design case. '
                'Please, make sure no connections have been modified or '
                'components have been relabeled for your offdesign '
                'calculation.'
            )
            logger.exception(msg, c.label)
            raise hlp.TESPyNetworkError(msg)

        conn = df.loc[c.label]
        for var in fpd.keys():
            c.get_attr(var).design = hlp.convert_to_SI(
                var, conn[var], conn[f"{var}_unit"]
            )
        c.vol.design = c.v.design / c.m.design
        for fluid in c.fluid.val:
            c.fluid.design[fluid] = conn[fluid]

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
                    param.is_var = True

                for var in c.offdesign:
                    param = c.get_attr(var)
                    param.is_set = True
                    param.is_var = False
                    param.val_SI = param.design

                c.new_design = False

            self._assign_variable_space(c)

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
                    msg = (msg[:-2] + ' to design value at component ' +
                           cp.label + '.')
                    logger.debug(msg)

            # start component initialisation
            cp.preprocess(self.num_conn_vars + self.num_comp_vars)
            ct = cp.__class__.__name__
            for spec in self.specifications[ct].keys():
                if len(cp.get_attr(self.specifications['lookup'][spec])) > 0:
                    self.specifications[ct][spec].loc[cp.label] = (
                        cp.get_attr(self.specifications['lookup'][spec]))

            i = self.num_conn_vars + self.num_comp_vars
            for container, name in cp.vars.items():
                self.variables_dict[i] = {"obj": container, "variable": name}
                i += 1
            cp.new_design = False
            self.num_comp_vars += cp.num_vars
            self.num_comp_eq += cp.num_eq

        msg = 'Switched components from design to offdesign.'
        logger.debug(msg)

        # count bus equations and reindex bus dictionary
        for b in self.busses.values():
            self.busses[b.label] = b
            self.num_bus_eq += b.P.is_set * 1

    def init_fluids(self):
        r"""
        Initialise the fluid vector on every connection of the network.

        - Create fluid vector for every component as dict,
          index: nw.fluids,
          values: 0 if not set by user.
        - Create fluid_set vector with same logic,
          index: nw.fluids,
          values: False if not set by user.
        - If there are any combustion chambers in the network, calculate fluid
          vector starting from there.
        - Propagate fluid vector in direction of sources and targets.
        """
        # stop fluid propagation for single fluid networks
        if len(self.all_fluids) == 1:
            return

        # fluid propagation from set values
        for c in self.conns['object']:
            if any(c.fluid.is_set.values()):
                c.target.propagate_fluid_to_target(c, c, entry_point=True)
                c.source.propagate_fluid_to_source(c, c, entry_point=True)

        # To save resources:
        # find empty fluid data and propagate from connections with data that
        # are interfaced directly to those connections by a component
            # if any(c.fluid.val0.values()):
            #     c.target.propagate_fluid_to_target(c, c, entry_point=True)
            #     c.source.propagate_fluid_to_source(c, c, entry_point=True)

        # fluid starting value generation based on components
        for cp in self.comps['object']:
            cp.initialise_fluids()

        msg = 'Fluid initialisation done.'
        logger.debug(msg)

    def init_properties(self):
        """
        Initialise the fluid properties on every connection of the network.

        - Set generic starting values for mass flow, enthalpy and pressure if
          not user specified, read from :code:`ìnit_path` or available from
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
                conn = df.loc[
                    df['source'].isin([c.source.label]) &
                    df['target'].isin([c.target.label]) &
                    df['source_id'].isin([c.source_id]) &
                    df['target_id'].isin([c.target_id])]
                try:
                    conn_id = conn.index[0]
                    # overwrite SI-values with values from init_file,
                    # except user specified values
                    for prop in ['m', 'p', 'h']:
                        data = c.get_attr(prop)
                        data.val0 = df.loc[conn_id, prop]
                        data.unit = df.loc[conn_id, prop + '_unit']

                    for fluid in self.fluids:
                        if not c.fluid.is_set[fluid]:
                            c.fluid.val[fluid] = df.loc[conn_id, fluid]
                        c.fluid.val0[fluid] = c.fluid.val[fluid]

                    c.good_starting_values = True

                except IndexError:
                    msg = (
                        'Could not find connection ' + c.label + ' in '
                        'connections.csv of init_path ' + self.init_path + '.')
                    logger.debug(msg)

            if sum(c.fluid.val.values()) == 0:
                msg = (
                    'The starting value for the fluid composition of the '
                    'connection ' + c.label + ' is empty. This might lead to '
                    'issues in the initialisation and solving process as '
                    'fluid property functions can not be called. Make sure '
                    'you specified a fluid composition in all parts of the '
                    'network.')
                logger.warning(msg)

            for key in ['m', 'p', 'h']:
                if not c.good_starting_values:
                    self.init_val0(c, key)
                if c.get_attr(key).is_var:
                    c.get_attr(key).val_SI = hlp.convert_to_SI(
                        key, c.get_attr(key).val0, c.get_attr(key).unit)

            self.init_count_connections_parameters(c)

        for c in self.conns['object']:
            if not c.good_starting_values:
                for key in ['m', 'p', 'h', 'T']:
                    if (c.get_attr(key).ref_set and
                            not c.get_attr(key).val_set):
                        c.get_attr(key).val_SI = (
                                c.get_attr(key).ref.obj.get_attr(key).val_SI *
                                c.get_attr(key).ref.factor +
                                c.get_attr(key).ref.delta_SI)

                self.init_precalc_properties(c)

            # starting values for specified subcooling/overheating
            # and state specification. These should be recalculated even with
            # good starting values, for example, when one exchanges enthalpy
            # with boiling point temperature difference.
            if ((c.Td_bp.is_set or c.state.is_set) and
                    not c.h.is_set):
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
        local_vars = self.specifications['Connection'].columns[:9]
        row = [c.get_attr(var).is_set for var in fpd.keys()]
        # write information to specifaction dataframe
        self.specifications['Connection'].loc[c.label, local_vars] = row

        # row = [c.get_attr(var).ref_set for var in local_vars]
        # write refrenced value information to specifaction dataframe
        # self.specifications['Ref'].loc[c.label, local_vars] = row

        # variables 9 to last but one: fluid mass fractions
        fluids = self.specifications['Connection'].columns[9:-1]
        row = [True if f in c.fluid.is_set else False for f in fluids]
        self.specifications['Connection'].loc[c.label, fluids] = row

        # last one: fluid balance specification
        self.specifications['Connection'].loc[
            c.label, 'balance'] = c.fluid.balance

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
        if c.x.is_set and not c.h.is_set:
            try:
                c.h.val_SI = fp.h_mix_pQ(c.p.val_SI, c.x.val_SI, c.fluid_data, c.mixing_rule)
            except ValueError:
                pass

        if c.T.is_set and not c.h.is_set:
            try:
                c.h.val_SI = fp.h_mix_pT(c.p.val_SI, c.T.val_SI, c.fluid_data, c.mixing_rule)
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
        if np.isnan(c.get_attr(key).val0):
            # starting value for mass flow is random between 1 and 2 kg/s
            # (should be generated based on some hash maybe?)
            if key == 'm':
                c.get_attr(key).val0 = np.random.random() + 1

            # generic starting values for pressure and enthalpy
            else:
                # retrieve starting values from component information
                val_s = c.source.initialise_source(c, key)
                val_t = c.target.initialise_target(c, key)

                if val_s == 0 and val_t == 0:
                    if key == 'p':
                        c.get_attr(key).val0 = 1e5
                    elif key == 'h':
                        c.get_attr(key).val0 = 1e6

                elif val_s == 0:
                    c.get_attr(key).val0 = val_t
                elif val_t == 0:
                    c.get_attr(key).val0 = val_s
                else:
                    c.get_attr(key).val0 = (val_s + val_t) / 2

                # change value according to specified unit system
                c.get_attr(key).val0 = hlp.convert_from_SI(
                    key, c.get_attr(key).val0, self.get_attr(key + '_unit'))

    @staticmethod
    def init_read_connections(base_path):
        r"""
        Read connection information from base_path.

        Parameters
        ----------
        base_path : str
            Path to network information.
        """
        path = hlp.modify_path_os(base_path + '/connections.csv')
        df = pd.read_csv(path, index_col=0, delimiter=';', decimal='.')
        return df

    def solve(self, mode, init_path=None, design_path=None,
              max_iter=50, min_iter=4, init_only=False, init_previous=True,
              use_cuda=False, always_all_equations=True, print_results=True):
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

        always_all_equations : boolean
            Calculate all equations in every iteration. Disabling this flag,
            will increase calculation speed, especially for mixtures, default:
            :code:`True`.

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
        self.always_all_equations = always_all_equations

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

        self.postprocessing()
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
                     and self.residual_history[-1] < ERR ** 0.5)
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
        component = '' if self.num_comp_vars == 0 else 'component'
        msg = self.iterinfo_fmt.format(
            iter='iter',
            residual='residual',
            progress='progress',
            massflow='massflow',
            pressure='pressure',
            enthalpy='enthalpy',
            fluid='fluid',
            component=component
        )
        logger.progress(0, msg)
        msg2 = '-' * 7 + '+------------' * 6 + "+"
        if self.num_comp_vars > 0:
            msg2 += '-------------'

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

        if not np.isnan(residual_norm):
            residual = '{:.2e}'.format(residual_norm)

        if not self.lin_dep and not np.isnan(residual_norm):
            massflow = '{:.2e}'.format(norm(self.increment[m]))
            pressure = '{:.2e}'.format(norm(self.increment[p]))
            enthalpy = '{:.2e}'.format(norm(self.increment[h]))
            fluid = '{:.2e}'.format(norm(self.increment[fl]))
            component  = '{:.2e}'.format(norm(self.increment[cp]))

        progress_val = -1
        if not np.isnan(residual_norm):
            # This should not be hardcoded here.
            if residual_norm > np.finfo(float).eps * 100:
                progress_min = np.log(ERR)
                progress_max = np.log(ERR ** 0.5) * -1
                progress_val = np.log(max(residual_norm, ERR)) * -1
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
        num_ips = num_iter / clc_time if clc_time > 1e-10 else np.Inf
        msg = '-' * 7 + '+------------' * 7
        logger.progress(100, msg)
        msg = 'Total iterations: {0:d}, Calculation time: {1:.2f} s, Iterations per second: {2:.2f}'.format(num_iter, clc_time, num_ips)
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
                    -cu.asarray(self.residual)))
            else:
                self.increment = np.linalg.inv(
                    self.jacobian).dot(-self.residual)
            self.lin_dep = False
        except np.linalg.linalg.LinAlgError:
            self.increment = self.residual * 0

    def update_variables(self):
        # add the increment
        for data in self.variables_dict.values():
            if data["variable"] in ["m", "h"]:
                container = data["obj"].get_attr(data["variable"])
                container.val_SI += self.increment[container.J_col]
            elif data["variable"] == "p":
                container = data["obj"].p
                increment = self.increment[container.J_col]
                relax = max(1, -2 * increment / container.val_SI)
                container.val_SI += increment / relax
            elif data["variable"] == "fluid":
                container = data["obj"].fluid
                container.val[data["fluid"]] += self.increment[
                    container.J_col[data["fluid"]]
                ]

                if container.val[data["fluid"]] < ERR :
                    container.val[data["fluid"]] = 0
                elif container.val[data["fluid"]] > 1 - ERR :
                    container.val[data["fluid"]] = 1
            else:
                # add increment
                data["obj"].val += self.increment[data["obj"].J_col]

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
        self.check_variable_bounds()

    def check_connection_properties(self, c):
        r"""
        Check for invalid fluid property values.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            Connection to check fluid properties.
        """
        fl = hlp.single_fluid(c.fluid.val)

        # pure fluid
        if fl is not None:
            # pressure
            if c.p.is_var:
                c.check_pressure_bounds(fl)

            # enthalpy
            if c.h.is_var:
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
        if c.m.val_SI <= self.m_range_SI[0] and not c.m.is_set:
            c.m.val_SI = self.m_range_SI[0]
            logger.debug(self._property_range_message(c, 'm'))

        elif c.m.val_SI >= self.m_range_SI[1] and not c.m.is_set:
            c.m.val_SI = self.m_range_SI[1]
            logger.debug(self._property_range_message(c, 'm'))

    def solve_components(self):
        r"""
        Calculate the residual and derivatives of component equations.

        - Iterate through components in network to get residuals and
          derivatives.
        - Place residual values in residual value vector of the network.
        - Place partial derivatives in Jacobian matrix of the network.
        """
        # fetch component equation residuals and component partial derivatives
        sum_eq = 0
        for cp in self.comps['object']:
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

        - Iterate through connections in network to get residuals and
          derivatives.
        - Place residual values in residual value vector of the network.
        - Place partial derivatives in Jacobian matrix of the network.
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

        - Iterate through user defined functions and calculate residual value
          and corresponding jacobian.
        - Place residual values in residual value vector of the network.
        - Place partial derivatives regarding connection parameters in Jacobian
          matrix of the network.
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

        - Iterate through busses in network to get residuals and derivatives.
        - Place residual values in residual value vector of the network.
        - Place partial derivatives in jacobian matrix of the network.
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
            if hasattr(c, "_m_tmp"):
                value = c.m.val_SI
                unit = c.m.unit
                c.m = c._m_tmp
                c.m.val_SI = value
                c.m.unit = unit
                del c._m_tmp
            if hasattr(c, "_fluid_tmp"):
                val = c.fluid.val
                c.fluid = c._fluid_tmp
                c.fluid.val = val
                del c._fluid_tmp

            c.good_starting_values = True
            c.calc_results()

            self.results['Connection'].loc[c.label] = (
                [
                    _ for key in fpd.keys()
                    for _ in [c.get_attr(key).val, c.get_attr(key).unit]
                ] + [
                    c.fluid.val[fluid] if fluid in c.fluid.val else np.nan
                    for fluid in self.all_fluids
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

            b.P.val = self.results[b.label]['bus value'].sum()

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
                        result += ('\n##### RESULTS (' + cp + ') #####\n')
                        result += (
                            tabulate(
                                df, headers='keys', tablefmt='psql',
                                floatfmt='.2e'
                            )
                        )

        # connection properties
        df = self.results['Connection'].loc[:, ['m', 'p', 'h', 'T']]
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
                    :, ['component value', 'bus value', 'efficiency']]
                df.loc['total'] = df.sum()
                df.loc['total', 'efficiency'] = np.nan
                if colored and b.P.is_set:
                    df.loc['total', 'bus value'] = (
                        coloring['set'] + str(df.loc['total', 'bus value']) +
                        coloring['end'])
                result += ('\n##### RESULTS (Bus: ' + b.label + ') #####\n')
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

    def save(self, path, **kwargs):
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
        if path[-1] != '/' and path[-1] != '\\':
            path += '/'
        path = hlp.modify_path_os(path)

        logger.debug('Saving network to path %s.', path)
        # creat path, if non existent
        if not os.path.exists(path):
            os.makedirs(path)

        # create path for component folder if non existent
        path_comps = hlp.modify_path_os(path + 'components/')
        if not os.path.exists(path_comps):
            os.makedirs(path_comps)

        # save all network information
        self.save_network(path + 'network.json')
        self.save_connections(path + 'connections.csv')
        self.save_components(path_comps)
        # self.save_busses(path_comps + 'bus.csv')
        # self.save_characteristics(path_comps)

    def save_network(self, fn):
        r"""
        Save basic network configuration.

        Parameters
        ----------
        fn : str
            Path/filename for the network configuration file.
        """
        with open(fn, 'w') as f:
            f.write(json.dumps(self._serialize(), indent=4))

        logger.debug('Network information saved to %s.', fn)

    def save_connections(self, fn):
        r"""
        Save the connection properties.

        - Uses connections object id as row identifier and saves

            - connections source and target as well as
            - properties with references and
            - fluid vector (including user specification if structure is True).

        - Connections source and target are identified by its labels.

        Parameters
        ----------
        fn : str
            Path/filename for the file.
        """
        self.results["Connection"].to_csv(
            fn, sep=';', decimal='.', index=True, na_rep='nan'
        )
        logger.debug('Connection information saved to %s.', fn)

    def save_components(self, path):
        r"""
        Save the component properties.

        - Uses components labels as row identifier.
        - Writes:

            - component's incomming and outgoing connections (object id) and
            - component's parametrisation.

        Parameters
        ----------
        path : str
            Path/filename for the file.
        """
        for c in self.comps['comp_type'].unique():
            if self.results[c].size > 0:
                fn = path + c + '.csv'
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
            df = pd.DataFrame(
                {'id': self.busses.values()}, index=self.busses.values(),
                dtype='object')
            df['label'] = df.apply(Network.get_props, axis=1, args=('label',))
            df['P'] = df.apply(Network.get_props, axis=1, args=('P', 'val'))
            df['P_set'] = df.apply(Network.get_props, axis=1,
                                   args=('P', 'is_set'))
            df.drop('id', axis=1, inplace=True)

            df.set_index('label', inplace=True)
            df.to_csv(fn, sep=';', decimal='.', index=True, na_rep='nan')
            logger.debug('Bus information saved to %s.', fn)

    def save_characteristics(self, path):
        r"""
        Save the characteristics.

        Parameters
        ----------
        fn : str
            Path/filename for the file.
        """
        # characteristic lines in components
        char_lines = []
        char_maps = []
        for c in self.comps['object']:
            for _col, data in c.variables.items():
                if isinstance(data, dc_cc):
                    char_lines += [data.char_func]
                elif isinstance(data, dc_cm):
                    char_maps += [data.char_func]

        # characteristic lines in busses
        for bus in self.busses.values():
            for c in bus.comps.index:
                ch = bus.comps.loc[c, 'char']
                if ch not in char_lines:
                    char_lines += [ch]

        # characteristic line export
        if len(char_lines) > 0:
            # get id and data
            df = pd.DataFrame(
                {'id': char_lines}, index=char_lines, dtype='object')
            df['id'] = df.apply(Network.get_id, axis=1)
            df['type'] = df.apply(Network.get_class_base, axis=1)

            cols = ['x', 'y', 'extrapolate']
            for val in cols:
                df[val] = df.apply(Network.get_props, axis=1, args=(val,))

            # write to char.csv
            fn = path + 'char_line.csv'
            df.to_csv(fn, sep=';', decimal='.', index=False, na_rep='nan')
            logger.debug('Characteristic line information saved to %s.', fn)

        if len(char_maps) > 0:
            # get id and data
            df = pd.DataFrame(
                {'id': char_maps}, index=char_maps, dtype='object')
            df['id'] = df.apply(Network.get_id, axis=1)
            df['type'] = df.apply(Network.get_class_base, axis=1)

            cols = ['x', 'y', 'z']
            for val in cols:
                df[val] = df.apply(Network.get_props, axis=1, args=(val,))

            # write to char_map.csv
            fn = path + 'char_map.csv'
            df.to_csv(fn, sep=';', decimal='.', index=False, na_rep='nan')
            logger.debug('Characteristic map information saved to %s.', fn)

    @staticmethod
    def get_id(c):
        """Return the id of the python object."""
        return str(c.name)[str(c.name).find(' at ') + 4:-1]

    @staticmethod
    def get_class_base(c):
        """Return the class name."""
        return c.name.__class__.__name__

    @staticmethod
    def get_props(c, *args):
        """Return properties."""
        if hasattr(c.name, args[0]):
            if (not isinstance(c.name.get_attr(args[0]), int) and
                    not isinstance(c.name.get_attr(args[0]), str) and
                    not isinstance(c.name.get_attr(args[0]), float) and
                    not isinstance(c.name.get_attr(args[0]), list) and
                    not isinstance(c.name.get_attr(args[0]), np.ndarray) and
                    not isinstance(c.name.get_attr(args[0]), con.Connection)):
                if len(args) == 1:
                    return c.name.get_attr(args[0])
                elif args[0] == 'fluid' and args[1] != 'balance':
                    return c.name.fluid.get_attr(args[1])[args[2]]
                elif args[1] == 'ref':
                    obj = c.name.get_attr(args[0]).get_attr(args[1])
                    if obj is not None:
                        return obj.get_attr(args[2])
                    else:
                        return np.nan
                else:
                    return c.name.get_attr(args[0]).get_attr(args[1])
            elif isinstance(c.name.get_attr(args[0]), np.ndarray):
                if len(c.name.get_attr(args[0]).shape) > 1:
                    return tuple(c.name.get_attr(args[0]).tolist())
                else:
                    return c.name.get_attr(args[0]).tolist()
            else:
                return c.name.get_attr(args[0])

    @staticmethod
    def get_busses(c, *args):
        """Return the list of busses a component is integrated in."""
        busses = []
        for bus in args[0]:
            if c.name in bus.comps.index:
                busses += [bus.label]
        return busses

    @staticmethod
    def get_bus_data(c, *args):
        """Return bus information of a component."""
        items = []
        if args[1] == 'char':
            for bus in args[0]:
                if c.name in bus.comps.index:
                    val = bus.comps.loc[c.name, args[1]]
                    items += [str(val)[str(val).find(' at ') + 4:-1]]

        else:
            for bus in args[0]:
                if c.name in bus.comps.index:
                    items += [bus.comps.loc[c.name, args[1]]]

        return items
