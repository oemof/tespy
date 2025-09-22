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
import importlib
import json
import math
import os
import warnings
from time import time

import numpy as np
import pandas as pd
from numpy.linalg import norm
from tabulate import tabulate

from tespy.components.component import component_registry
from tespy.connections import Bus
from tespy.connections import Connection
from tespy.connections.connection import ConnectionBase
from tespy.connections.connection import connection_registry
from tespy.tools import helpers as hlp
from tespy.tools import logger
from tespy.tools.characteristics import CharLine
from tespy.tools.characteristics import CharMap
from tespy.tools.data_containers import ComponentCharacteristicMaps as dc_cm
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import DataContainer as dc
from tespy.tools.data_containers import FluidProperties as dc_prop
from tespy.tools.data_containers import ScalarVariable as dc_scavar
from tespy.tools.data_containers import VectorVariable as dc_vecvar
from tespy.tools.global_vars import ERR
from tespy.tools.global_vars import fluid_property_data as fpd
from tespy.tools.units import SI_UNITS
from tespy.tools.units import Units

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
    iterinfo : boolean
        Print convergence progress to console.

    h_range : list
        List with minimum and maximum values for enthalpy value range.

    m_range : list
        List with minimum and maximum values for mass flow value range.

    p_range : list
        List with minimum and maximum values for pressure value range.

    Note
    ----
    Units are specified via the :code:`Network.units.set_defaults` interface.
    The specification is optional and will use SI units by default.

    Range specification is optional, too. The value range is used to stabilize
    the newton algorithm. For more information see the "getting started"
    section in the online-documentation.

    Example
    -------
    Basic example for a setting up a :code:`tespy.networks.network.Network`
    object.

    Standard value for iterinfo is :code:`True`. This will print out
    convergence progress to the console. You can stop the printouts by setting
    this property to :code:`False`.

    >>> from tespy.networks import Network
    >>> mynetwork = Network()
    >>> mynetwork.units.set_defaults(**{
    ...     "pressure": "bar", "temperature": "degC"
    ... })
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
    >>> from tespy.components import Source, Sink, Pipe, PowerSink
    >>> from tespy.connections import Connection, PowerConnection
    >>> nw = Network()
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "temperature": "degC"
    ... })
    >>> so = Source('source')
    >>> si = Sink('sink')
    >>> p = Pipe('pipe', Q=0, pr=0.95, printout=False, power_connector_location="outlet")
    >>> h = PowerSink('heat to ambient')
    >>> a = Connection(so, 'out1', p, 'in1')
    >>> b = Connection(p, 'out1', si, 'in1')
    >>> nw.add_conns(a, b)
    >>> a.set_attr(fluid={'CH4': 1}, T=30, p=10, m=10, printout=False)
    >>> b.set_attr(printout=False)
    >>> e = PowerConnection(p, 'heat', h, 'power', printout=False)
    >>> nw.add_conns(e)
    >>> nw.set_attr(iterinfo=False)
    >>> nw.solve('design')
    >>> nw.print_results()
    """

    def __init__(self, **kwargs):
        self.set_defaults()
        self.set_attr(**kwargs)

    def _serialize(self):
        return {
            "m_range": list(self.m_range.magnitude),
            "p_range": list(self.p_range.magnitude),
            "h_range": list(self.h_range.magnitude),
            "units": self.units._serialize()
        }

    def set_defaults(self):
        """Set default network properties."""
        # connection dataframe

        dtypes={
            "object": object,
            "source": object,
            "source_id": str,
            "target": object,
            "target_id": str,
            "conn_type": str
        }
        self.conns = pd.DataFrame(columns=list(dtypes.keys())).astype(dtypes)
        self.all_fluids = set()
        # component dataframe
        dtypes = {
            "comp_type": str,
            "object": object,
        }
        self.comps = pd.DataFrame(columns=list(dtypes.keys())).astype(dtypes)
        # user defined function dictionary for fast access
        self.user_defined_eq = {}
        # bus dictionary
        self.busses = {}
        self.subsystems = {}
        # results and specification dictionary
        self.results = {}

        # in case of a design calculation after an offdesign calculation
        self.redesign = False

        self.checked = False
        self.design_path = None
        self.iterinfo = True
        self.units = Units()

        msg = 'Default unit specifications:\n'
        for prop, unit in self.units.default.items():
            # standard unit set
            msg += f"{prop}: {unit}" + "\n"

        # don't need the last newline
        logger.debug(msg[:-1])

        # generic value range
        self.m_range_SI = [-1e12, 1e12]
        self.p_range_SI = [2e2, 300e5]
        self.h_range_SI = [1e3, 7e6]

        property_names = {"m": "mass_flow", "p": "pressure", "h": "enthalpy"}
        for prop, name in property_names.items():
            limits = self.get_attr(f"{prop}_range_SI")
            msg = (
                f"Default {name} limits\n"
                f"min: {limits[0]} {self.units._quantities[name]}\n"
                f"max: {limits[1]} {self.units._quantities[name]}"
            )
            logger.debug(msg)

            unit = self.units.default[name]
            key = f"{prop}_range"
            self.__dict__.update({
                key: self.units.ureg.Quantity(
                    np.array(self.get_attr(f"{key}_SI")), unit
                )
            })

    def set_attr(self, **kwargs):
        r"""
        Set, resets or unsets attributes of a network.

        Parameters
        ----------
        iterinfo : boolean
            Print convergence progress to console.

        h_range : list
            List with minimum and maximum values for enthalpy value range.

        m_range : list
            List with minimum and maximum values for mass flow value range.

        p_range : list
            List with minimum and maximum values for pressure value range.
        """
        self.units = kwargs.get('units', self.units)
        unit_replace = {
            "C": "degC",
            "J / kgK": "J / (kg * K)",
            "kJ / kgK": "kJ / (kg * K)",
            "MJ / kgK": "MJ / (kg * K)",
        }
        # unit sets
        msg = None
        for prop in fpd.keys():
            unit = f'{prop}_unit'
            if unit in kwargs:
                if msg is None:
                    msg = (
                        "The API for specification of units in a Network "
                        "changed. The old variant will be removed in the next "
                        "major release. Please use the "
                        "'Network.units.set_defaults' method instead."
                    )
                # for backwards compatibility: Update in the default units
                self.units.set_defaults(**{
                    fpd[prop]["text"].replace(" ", "_"):
                    unit_replace.get(kwargs[unit], kwargs[unit])
                })

        if msg:
            warnings.warn(msg, FutureWarning)

        for prop in ['m', 'p', 'h']:
            key = f"{prop}_range"
            if key in kwargs:
                if isinstance(kwargs[key], list):
                    quantity = fpd[prop]["text"].replace(" ", "_")
                    unit = self.units.default[quantity]
                    self.__dict__.update({
                        key: self.units.ureg.Quantity(
                            np.array(kwargs[key]),
                            unit
                        )
                    })
                    self.__dict__.update({
                        f"{key}_SI":
                        self.get_attr(key).to(SI_UNITS[quantity]).magnitude
                    })
                else:
                    msg = f'Specify the range as list: [{prop}_min, {prop}_max]'
                    logger.error(msg)
                    raise TypeError(msg)

                limits = self.get_attr(f'{key}_SI')
                msg = (
                    f'Setting {fpd[prop]["text"]} limits\n'
                    f'min: {limits[0]} {fpd[prop]["SI_unit"]}\n'
                    f'max: {limits[1]} {fpd[prop]["SI_unit"]}'
                )
                logger.debug(msg)

        self.iterinfo = kwargs.get('iterinfo', self.iterinfo)

        if not isinstance(self.iterinfo, bool):
            msg = 'Network parameter iterinfo must be True or False!'
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

    def add_subsystems(self, *args):
        r"""
        Add one or more subsystems to the network.

        Parameters
        ----------
        c : tespy.components.subsystem.Subsystem
            The subsystem to be added to the network, subsystem objects si
            :code:`network.add_subsystems(s1, s2, s3, ...)`.
        """
        for subsystem in args:
            if subsystem.label in self.subsystems:
                msg = (
                    'There is already a subsystem with the label '
                    f'{subsystem.label}. The labels must be unique!'
                )
                logger.error(msg)
                raise ValueError(msg)

            self.subsystems[subsystem.label] = subsystem

            for c in subsystem.conns.values():
                self.add_conns(c)

    def del_subsystems(self, *args):
        r"""
        Delete one or more subsystems from the network.

        Parameters
        ----------
        c : tespy.components.subsystem.Subsystem
            The subsystem to be deleted from the network, subsystem objects si
            :code:`network.del_subsystems(s1, s2, s3, ...)`.
        """
        for subsystem in args:
            if subsystem.label in self.subsystems:
                for c in subsystem.conns.values():
                    self.del_conns(c)

            del self.subsystems[subsystem.label: str]

    def get_subsystem(self, label):
        r"""
        Get Subsystem via label.

        Parameters
        ----------
        label : str
            Label of the Subsystem object.

        Returns
        -------
        tespy.components.subsystem.Subsystem
            Subsystem objectt with specified label, None if no Subsystem of
            the network has this label.
        """
        try:
            return self.subsystems[label]
        except KeyError:
            logger.warning(f"Subsystem with label {label} not found.")
            return None

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
            if not isinstance(c, ConnectionBase):
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

            conn_type = c.__class__.__name__
            self.conns.loc[c.label] = [
                c, c.source, c.source_id, c.target, c.target_id, conn_type
            ]

            msg = f'Added connection {c.label} to network.'
            logger.debug(msg)
            # set status "checked" to false, if connection is added to network.
            self.checked = False

        self.conns = self.conns.sort_index()

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
            if c.__class__.__name__ in self.results:
                self.results[c.__class__.__name__].drop(
                    c.label, inplace=True, errors="ignore"
                )
            msg = f'Deleted connection {c.label} from network.'
            logger.debug(msg)

        self._del_comps(comps)

        # set status "checked" to false, if connection is deleted from network.
        self.checked = False

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
        self.comps = self.comps.sort_index()

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

    def assert_convergence(self):
        """Check convergence status of a simulation."""
        msg = 'Calculation did not converge!'
        assert self.converged, msg

    @property
    def converged(self):
        if hasattr(self, "status"):
            return self.status == 0 or self.status == 1
        else:
            msg = (
                "The converged attribute can only be accessed after the first "
                "call of the solve method"
            )
            raise AttributeError(msg)

    def check_busses(self, b):
        r"""
        Checksthe busses to be added for type, duplicates and identical labels.

        Parameters
        ----------
        b : tespy.connections.bus.Bus
            The bus to be checked.
        """
        if isinstance(b, Bus):
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

    def check_topology(self):
        r"""Check if components are connected properly within the network."""
        if len(self.conns) == 0:
            msg = (
                'No connections have been added to the network, please make '
                'sure to add your connections with the .add_conns() method.'
            )
            logger.error(msg)
            raise hlp.TESPyNetworkError(msg)

        self._check_connections()
        self._init_components()
        self._check_components()
        self._create_fluid_wrapper_branches()

        # network checked
        self.checked = True
        msg = 'Networkcheck successful.'
        logger.info(msg)

    def _check_connections(self):
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

    def _init_connection_result_datastructure(self):

        for conn_type in self.conns["conn_type"].unique():
            if conn_type in self.results:
                del self.results[conn_type]

        for conn in self.conns["object"]:
            conn_type = conn.__class__.__name__
            # this will move somewhere else!
            # set up results dataframe for connections
            # this should be done based on the connections
            if conn_type not in self.results:
                cols = conn._get_result_cols(set(self.all_fluids))
                self.results[conn_type] = pd.DataFrame(columns=cols, dtype='float64')

    def _init_components(self):
        r"""Set up necessary component information."""
        for comp in self.comps["object"]:
            # get incoming and outgoing connections of a component
            source_mask = self.conns["source"] == comp
            required_connectors_mask = self.conns["source_id"].isin(comp.outlets())
            sources = self.conns[source_mask & required_connectors_mask]
            sources = sources["source_id"].sort_values().index.tolist()
            target_mask = self.conns["target"] == comp
            required_connectors_mask = self.conns["target_id"].isin(comp.inlets())
            targets = self.conns[target_mask & required_connectors_mask]
            targets = targets["target_id"].sort_values().index.tolist()
            # save the incoming and outgoing as well as the number of
            # connections as component attribute
            comp.inl = self.conns.loc[targets, "object"].tolist()
            comp.outl = self.conns.loc[sources, "object"].tolist()
            comp.num_i = len(comp.inlets())
            comp.num_o = len(comp.outlets())

            required_connectors_mask = self.conns["source_id"].isin(comp.poweroutlets())
            sources = self.conns[source_mask & required_connectors_mask]
            sources = sources["source_id"].sort_values().index.tolist()

            required_connectors_mask = self.conns["target_id"].isin(comp.powerinlets())
            targets = self.conns[target_mask & required_connectors_mask]
            targets = targets["target_id"].sort_values().index.tolist()

            comp.power_inl = self.conns.loc[targets, "object"].tolist()
            comp.power_outl = self.conns.loc[sources, "object"].tolist()
            comp.num_power_i = len(comp.powerinlets())
            comp.num_power_o = len(comp.poweroutlets())

            # set up restults and specification dataframes
            comp_type = comp.__class__.__name__
            if comp_type not in self.results:
                cols = [
                    c for col, data in comp.parameters.items()
                    if isinstance(data, dc_cp)
                    for c in [col, f"{col}_unit"]
                ]
                self.results[comp_type] = pd.DataFrame(
                    columns=cols, dtype='float64'
                )

    def _check_components(self):
        # count number of incoming and outgoing connections and compare to
        # expected values
        for comp in self.comps['object']:
            if len(comp.outl) != comp.num_o:
                msg = (
                    f"The component {comp.label} is missing "
                    f"{comp.num_o - len(comp.outl)} outgoing connections. "
                    "Make sure all outlets are connected and all connections "
                    "have been added to the network."
                )
                logger.error(msg)
                # raise an error in case network check is unsuccesful
                raise hlp.TESPyNetworkError(msg)
            elif len(comp.inl) != comp.num_i:
                msg = (
                    f"The component {comp.label} is missing "
                    f"{comp.num_i - len(comp.inl)} incoming connections. "
                    "Make sure all inlets are connected and all connections "
                    "have been added to the network."
                )
                logger.error(msg)
                # raise an error in case network check is unsuccesful
                raise hlp.TESPyNetworkError(msg)

            # this rule only applies, in case there are any power connections
            if len(comp.power_inl) + len(comp.power_outl) > 0:
                if len(comp.power_outl) != comp.num_power_o:
                    msg = (
                        f"The component {comp.label} is missing "
                        f"{comp.num_power_o - len(comp.power_outl)} outgoing "
                        "power connections. Make sure all outlets are "
                        "connected and all connections have been added to the "
                        "network."
                    )
                    logger.error(msg)
                    # raise an error in case network check is unsuccesful
                    raise hlp.TESPyNetworkError(msg)
                elif len(comp.power_inl) != comp.num_power_i:
                    msg = (
                        f"The component {comp.label} is missing "
                        f"{comp.num_power_i - len(comp.power_inl)} incoming "
                        "power connections. Make sure all outlets are "
                        "connected and all connections have been added to the "
                        "network."
                    )
                    logger.error(msg)
                    # raise an error in case network check is unsuccesful
                    raise hlp.TESPyNetworkError(msg)

    def _prepare_problem(self):
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
        self.variable_counter = 0
        self.variables_dict = {}

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
        for k, v in self.fluid_wrapper_branches.items():
            if self.conns.loc[v["connections"][0].label, "object"] != v["connections"][0]:
                self._create_fluid_wrapper_branches()
            continue

        self._propagate_fluid_wrappers()
        self._init_connection_result_datastructure()
        self._prepare_solve_mode()
        # this method will distribute units and set SI values from given values
        # and units
        self._init_set_properties()
        self._create_structure_matrix()

        self._presolve()
        self._prepare_for_solver()

        # generic fluid property initialisation
        self._init_properties()

        msg = 'Network initialised.'
        logger.info(msg)

    def _propagate_fluid_wrappers(self):

        connections_in_wrapper_branches = []
        self.all_fluids = []
        for branch_data in self.fluid_wrapper_branches.values():
            all_connections = [c for c in branch_data["connections"]]

            any_fluids_set = []
            engines = {}
            back_ends = {}
            any_fluids = []
            any_fluids0 = []
            mixing_rules = []
            for c in all_connections:
                for f in c.fluid.is_set:
                    any_fluids_set += [f]
                    if f in c.fluid.engine:
                        engines[f] = c.fluid.engine[f]
                    if f in c.fluid.back_end:
                        back_ends[f] = c.fluid.back_end[f]
                any_fluids += list(c.fluid.val.keys())
                any_fluids0 += list(c.fluid.val0.keys())
                if c.mixing_rule is not None:
                    mixing_rules += [c.mixing_rule]

            mixing_rule = list(set(mixing_rules))
            if len(mixing_rule) > 1:
                msg = (
                    "You have provided more than one mixing rule in the "
                    "branches including the following connections: "
                    f"{', '.join([c.label for c in all_connections])}"
                )
                raise hlp.TESPyNetworkError(msg)
            elif len(mixing_rule) == 0:
                mixing_rule = "ideal-cond"
            else:
                mixing_rule = mixing_rules[0]

            if not any_fluids_set:
                msg = "You are missing fluid specifications."

            all_components = [c for c in branch_data["components"]]
            for cp in all_components:
                any_fluids += cp._add_missing_fluids(all_connections)

            potential_fluids = set(any_fluids_set + any_fluids + any_fluids0)
            self.all_fluids += list(potential_fluids)
            num_potential_fluids = len(potential_fluids)
            if num_potential_fluids == 0:
                msg = (
                    "The follwing connections of your network are missing any "
                    "kind of fluid composition information:"
                    f"{', '.join([c.label for c in all_connections])}."
                )
                raise hlp.TESPyNetworkError(msg)

            for c in all_connections:
                c.mixing_rule = mixing_rule
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

            connections_in_wrapper_branches += all_connections

        mask = self.conns["conn_type"] == "Connection"
        missing_wrappers = (
            set(self.conns.loc[mask, "object"].tolist())
            - set(connections_in_wrapper_branches)
        )
        if len(missing_wrappers) > 0:
            msg = (
                f"The fluid information propagation for the connections "
                f"{', '.join([c.label for c in missing_wrappers])} failed. "
                "The reason for this is likely, that these connections do not "
                "have any Sources or a CycleCloser attached to them."
            )
            logger.error(msg)
            raise hlp.TESPyNetworkError(msg)

        self.all_fluids = set(self.all_fluids)

    def _prepare_solve_mode(self):
        if self.mode == 'offdesign':
            self.redesign = True
            if self.design_path is None:
                # must provide design_path
                msg = "Please provide a design_path for offdesign mode."
                logger.error(msg)
                raise hlp.TESPyNetworkError(msg)

            # load design case
            if self.new_design:
                self._init_offdesign_params()

            self._init_offdesign()

        else:
            # reset any preceding offdesign calculation
            self._init_design()

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
                self._variable_lookup[var]["object"]
                for var in linear_dependents["variables"]
            ]
            reference_container = self._reference_container_lookup[reference]
            reference_conn = all_connections[0]

            fluid_specs = [f for c in all_connections for f in c.fluid.is_set]
            fluid0 = {f: value for c in all_connections for f, value in c.fluid.val0.items()}
            if len(fluid_specs) == 0:

                if len(reference_conn._potential_fluids) > 1:
                    reference_container.is_var = {
                        f for f in reference_conn._potential_fluids
                    }
                    reference_container.val = {
                        f: 1 / len(reference_container.is_var)
                        for f in reference_container.is_var
                    }
                    # load up specification of starting values if any are
                    # available
                    reference_container.val.update(fluid0)
                else:
                    reference_container.val[
                        list(reference_conn._potential_fluids)[0]
                    ] = 1

            elif len(fluid_specs) != len(set(fluid_specs)):
                msg = (
                    "The mass fraction of a single fluid has been been "
                    "specified more than once in the following linear branch "
                    "of connections: "
                    f"{', '.join([c.label for c in all_connections])}."
                )
                raise hlp.TESPyNetworkError(msg)
            else:
                fixed_fractions = {
                    f: c.fluid._val[f]
                    for c in all_connections
                    for f in fluid_specs
                    if f in c.fluid._is_set
                }
                mass_fraction_sum = sum(fixed_fractions.values())
                if mass_fraction_sum > 1 + ERR:
                    msg = (
                        "The mass fraction of fluids within a linear branch "
                        "of connections cannot exceed 1: "
                        f"{', '.join([c.label for c in all_connections])}."
                    )
                    raise ValueError(msg)
                elif mass_fraction_sum < 1 - ERR:
                    # set the fluids with specified mass fraction
                    # remaining fluids are variable, create wrappers for them
                    all_fluids = reference_conn._potential_fluids
                    num_remaining_fluids = len(all_fluids) - len(fixed_fractions)
                    if num_remaining_fluids == 1:
                        missing_fluid = list(
                            set(all_fluids) - set(fixed_fractions.keys())
                        )[0]
                        fixed_fractions[missing_fluid] = 1 - mass_fraction_sum
                        variable = set()
                    else:
                        missing_fluids = (
                            set(all_fluids) - set(fixed_fractions.keys())
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
                num_var = len(variable)
                for f in variable:
                    reference_container.val[f] = (1 - mass_fraction_sum) / num_var
                    if f in fluid0:
                        reference_container.val[f] = fluid0[f]

            for fluid in reference_container.is_var:
                reference_container._J_col[fluid] = self.variable_counter
                self.variables_dict[self.variable_counter] = {
                    "obj": reference_container,
                    "variable": "fluid",
                    "fluid": fluid,
                    "_represents": linear_dependents["variables"]
                }
                self.variable_counter += 1

    def _create_structure_matrix(self):
        self._structure_matrix = {}
        self._rhs = {}
        self._variable_lookup = {}
        self._object_to_variable_lookup = {}
        self._equation_set_lookup = {}
        self._presolved_equations = []
        self._reference_container_lookup = {}
        self._equation_lookup = {}
        self._incidence_matrix = {}

        num_vars = self._prepare_variables()

        self._reassign_ude_objects()

        sum_eq = 0
        sum_eq = self._preprocess_network_parts(self.conns["object"], sum_eq)
        sum_eq = self._preprocess_network_parts(self.comps["object"], sum_eq)
        sum_eq = self._preprocess_network_parts(self.user_defined_eq.values(), sum_eq)

        _linear_dependencies = self._find_linear_dependent_variables(
            self._structure_matrix, self._rhs
        )
        _linear_dependent_variables = [
            var for linear_dependents in _linear_dependencies
            for var in linear_dependents["variables"]
        ]
        _missing_variables = [
            {
                "variables": [var],
                "reference": var,
                "factors": {var: 1.0},
                "offsets": {var: 0.0},
                "equation_indices": {}
            }
            for var in set(range(num_vars)) - set(_linear_dependent_variables)
        ]
        self._variable_dependencies = _missing_variables + _linear_dependencies

        for linear_dependents in self._variable_dependencies:
            reference_variable = self._variable_lookup[
                linear_dependents["reference"]
            ]
            reference = reference_variable["object"].get_attr(
                reference_variable["property"]
            )
            d = reference.d
            if hasattr(reference, "min_val"):
                min_val = reference.min_val
                max_val = reference.max_val
            else:
                min_val = None
                max_val = None

            if reference_variable["property"] != "fluid":
                DataContainer = dc_scavar
                reference_container = DataContainer(
                    _is_var=True,
                    _d=d,
                    min_val=min_val,
                    max_val=max_val,
                )
            else:
                DataContainer = dc_vecvar
                reference_container = DataContainer(
                    _d=d
                )

            self._reference_container_lookup[
                linear_dependents['reference']
            ] = reference_container

            for variable in linear_dependents["variables"]:
                variable_data = self._variable_lookup[variable]
                if variable_data["property"] != reference_variable["property"]:
                    msg = (
                        "There is a direct linear dependency between two "
                        "variables of different properties. This is unexpected "
                        "and might not work as intended: "
                        f"{reference_variable['object'].label}: "
                        f"{reference_variable['property']} and "
                        f"{variable_data['object'].label}: "
                        f"{variable_data['property']}."
                    )
                    logger.warning(msg)

                container = variable_data["object"].get_attr(variable_data["property"])
                container._reference_container = reference_container
                container._factor = linear_dependents["factors"][variable]
                container._offset = linear_dependents["offsets"][variable]

        # impose set values in the reference containers
        for conn in self.conns["object"]:
            for prop, container in conn.get_variables().items():
                if conn.get_attr(prop).is_set:
                    conn.get_attr(prop).set_reference_val_SI(conn.get_attr(prop)._val_SI)

        # collect all presolved equations
        self._presolved_equations = [
            indices
            for dependents in self._variable_dependencies
            for indices in dependents["equation_indices"].values()
        ]

    def _prepare_variables(self):
        num_vars = 0
        for conn in self.conns["object"]:
            for prop, container in conn.get_variables().items():
                # flag potential variables
                container._potential_var = not container.is_set
                container.sm_col = num_vars
                num_vars += 1

                self._variable_lookup[container.sm_col] = {
                    "object": conn, "property": prop
                }
                if conn not in self._object_to_variable_lookup:
                    self._object_to_variable_lookup[conn] = {}
                self._object_to_variable_lookup[conn].update(
                    {prop: container.sm_col}
                )

            if hasattr(conn, "fluid"):
                # fluid is handled separately
                container = conn.fluid
                container.sm_col = num_vars
                num_vars += 1
                self._variable_lookup[container.sm_col] = {
                    "object": conn, "property": "fluid"
                }
                if conn not in self._object_to_variable_lookup:
                    self._object_to_variable_lookup[conn] = {}
                self._object_to_variable_lookup[conn].update(
                    {"fluid": container.sm_col}
                )

        for comp in self.comps["object"]:
            for prop, container in comp.get_variables().items():
                container.sm_col = num_vars
                num_vars += 1

                self._variable_lookup[container.sm_col] = {
                    "object": comp, "property": prop
                }
                if comp not in self._object_to_variable_lookup:
                    self._object_to_variable_lookup[comp] = {}
                self._object_to_variable_lookup[comp].update(
                    {prop: container.sm_col}
                )
        return num_vars

    def _reassign_ude_objects(self):
        for ude in self.user_defined_eq.values():
            ude.conns = [self.get_conn(c.label) for c in ude.conns]
            ude.comps = [self.get_comp(c.label) for c in ude.comps]

    def _preprocess_network_parts(self, parts, eq_counter):

        for obj in parts:
            obj._preprocess(eq_counter)
            self._structure_matrix.update(obj._structure_matrix)
            self._rhs.update(obj._rhs)
            eq_map = {
                eq_num: (obj.label, eq_name)
                for eq_num, eq_name in obj._equation_set_lookup.items()
            }
            self._equation_set_lookup.update(eq_map)
            eq_counter += obj.num_eq

        return eq_counter

    def _find_linear_dependent_variables(self, sparse_matrix, rhs):
        if len(sparse_matrix) == 0:
            return []

        adjacency_list, eq_idx, edges_with_factors, rhs_offsets = (
            self._build_graph(sparse_matrix, rhs)
        )
        # Detect cycles (to check for circular dependencies)
        cycle = self._find_cycles_in_graph(
            {k: [x[0] for x in v] for k, v in adjacency_list.items()}
        )
        if cycle is not None:
            self._raise_error_if_cycle(cycle, edges_with_factors, eq_idx)

        # Find connected components and compute factors/offsets
        visited = set()
        variables_factors_offsets = []

        def dfs_component(node, current_factor, current_offset):
            """DFS to calculate factors and offsets relative to the reference
            variable."""
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
                        idx_f = neighbor, curr_node
                        idx_b = curr_node, neighbor
                        edge_offset = (
                            rhs_offsets.get(idx_b, 0.0)
                            or -rhs_offsets.get(idx_f, 0.0)
                        )
                        # Determine which equation to use
                        if idx_f in rhs_offsets:
                            equation_indices[idx_f] = eq_idx[idx_f]
                        else:
                            equation_indices[idx_b] = eq_idx[idx_b]

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
                factors, offsets, equation_indices = dfs_component(
                    reference, 1.0, 0.0
                )

                variables_factors_offsets.append({
                    'variables': list(factors.keys()),
                    'reference': reference,
                    'factors': factors,
                    'offsets': offsets,
                    'equation_indices': equation_indices
                })

        return variables_factors_offsets

    def _build_graph(self, sparse_matrix, rhs):
        edges_with_factors = []
        rhs_offsets = {}
        eq_idx = {}
        # The equation indices keep track of which equations to eliminate
        # Extract edges and offsets from rows with two non-zero entries
        rows = {k[0] for k in sparse_matrix}
        # sorting needs to be applied to always have same orientation on edges
        # otherwise duplicate edges are not found if one is just in reverse
        rows_with_cols = {
            row: sorted([k[1] for k in sparse_matrix if k[0] == row])
            for row in rows
        }
        for row, cols in rows_with_cols.items():
            if len(cols) == 2:
                non_zero_values = (
                    sparse_matrix[(row, cols[0])], sparse_matrix[(row, cols[1])]
                )
                col1, col2 = cols
                val1, val2 = non_zero_values
                factor = -val1 / val2
                offset = rhs[row] / val2
                edges_with_factors.append((col1, col2, factor))
                rhs_offsets[(col1, col2)] = offset
                if (col1, col2) in eq_idx:
                    variables = self._get_variables_before_presolve_by_number([col1, col2])
                    equations = self._get_equation_sets_by_eq_set_number(
                        [eq_idx[(col1, col2)], row]
                    )
                    msg = (
                        "The variables "
                        f"{', '.join([str(v) for v in variables])} are "
                        "directly linked with two equations "
                        f"{', '.join([str(e) for e in equations])}. This "
                        "overdetermines the problem."
                    )
                    raise hlp.TESPyNetworkError(msg)

                eq_idx[(col1, col2)] = row

        # Build adjacency list for the graph
        adjacency_list = {}
        for col1, col2, factor in edges_with_factors:
            if col1 not in adjacency_list:
                adjacency_list[col1] = []
            if col2 not in adjacency_list:
                adjacency_list[col2] = []

            # Add edge with factor and reverse edge with reciprocal value
            adjacency_list[col1].append((col2, factor))
            adjacency_list[col2].append((col1, 1 / factor))

        return adjacency_list, eq_idx, edges_with_factors, rhs_offsets

    def _find_cycles_in_graph(self, graph):
        visited = set()
        parent = {}

        def dfs(node, prev):
            visited.add(node)
            for neighbor in graph.get(node, []):
                if neighbor not in visited:
                    parent[neighbor] = node
                    result = dfs(neighbor, node)
                    if result:
                        return result
                elif neighbor != prev:
                    # Cycle found, reconstruct it
                    cycle = [neighbor, node]
                    while cycle[-1] != neighbor:
                        cycle.append(parent[cycle[-1]])
                    cycle.reverse()
                    return cycle
            return None

        for node in graph:
            if node not in visited:
                parent[node] = None
                cycle = dfs(node, None)
                if cycle:
                    return set(cycle)

        return None

    def _raise_error_if_cycle(self, cycle, edges_with_factors, eq_idx):
        edge_list = [
            e[:2] for e in edges_with_factors
            if e[0] in cycle or e[1] in cycle
        ]
        cycling_eqs = [v for k, v in eq_idx.items() if k in edge_list]
        variable_names = self._get_variables_before_presolve_by_number(cycle)
        equations = self._get_equation_sets_by_eq_set_number(cycling_eqs)
        msg = (
            "A circular dependency between the variables "
            f"{', '.join([str(v) for v in variable_names])} "
            "caused by the equations "
            f"{', '.join([str(e) for e in equations])} has been "
            "detected. This overdetermines the problem."
        )
        raise hlp.TESPyNetworkError(msg)

    def _create_fluid_wrapper_branches(self):

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

    def _presolve(self):
        # handle the fluid vector variables
        self._presolve_fluid_vectors()
        # set up the actual list of equations for connections, components,

        for c in self.conns['object']:
            self._presolved_equations += c._presolve()

        self._presolve_linear_dependents()

        # iteratively check presolvable fluid properties
        # and distribute presolved variables to all linear dependents
        # until the number of variables does not change anymore
        number_variables = sum([
            variable.is_var
            for conn in self.conns['object']
            for variable in conn.get_variables().values()
        ])
        while True:
            for c in self.conns['object']:
                self._presolved_equations += c._presolve()
            self._presolve_linear_dependents()
            reduced_variables = [
                variable.is_var
                for conn in self.conns['object']
                for variable in conn.get_variables().values()
            ]
            reduced_variables = sum(reduced_variables)
            if reduced_variables == number_variables:
                break

            number_variables = reduced_variables

    def _prepare_for_solver(self):
        for variable in self._variable_dependencies:
            reference = self._variable_lookup[variable["reference"]]
            represents = variable["variables"]
            if reference["property"] != "fluid":
                self._assign_variable_space(reference, represents)

        eq_counter = 0

        _eq_counter = self._prepare_network_parts(self.comps["object"], eq_counter)
        self.num_comp_eq = _eq_counter - eq_counter
        eq_counter = _eq_counter

        _eq_counter = self._prepare_network_parts(self.conns["object"], eq_counter)
        self.num_conn_eq = _eq_counter - eq_counter
        eq_counter = _eq_counter

        _eq_counter = self._prepare_network_parts(self.user_defined_eq.values(), eq_counter)
        self.num_ude_eq = _eq_counter - eq_counter
        eq_counter = _eq_counter

        for b in self.busses.values():
            self.busses[b.label] = b
            self.num_bus_eq += b.P.is_set * 1

    def _prepare_network_parts(self, parts, eq_counter):
        for obj in parts:
            eq_counter = obj._prepare_for_solver(self._presolved_equations, eq_counter)
            eq_map = {
                eq_num: (obj.label, eq_name)
                for eq_num, eq_name in obj._equation_lookup.items()
            }
            self._equation_lookup.update(eq_map)

            dependents_map = {
                eq_num: [dependent.J_col for dependent in dependents]
                for eq_num, dependents in obj._equation_scalar_dependents_lookup.items()
            }
            self._incidence_matrix.update(dependents_map)

            dependents_map = {
                eq_num: [
                    dependent.J_col[key]
                    for dependent, keys in dependents.items()
                    for key in keys
                ]
                for eq_num, dependents in obj._equation_vector_dependents_lookup.items()
            }
            for eq_num, dependents in dependents_map.items():
                if eq_num in self._incidence_matrix:
                    self._incidence_matrix[eq_num] += dependents

                else:
                    self._incidence_matrix[eq_num] = dependents

        return eq_counter

    def _presolve_linear_dependents(self):
        for linear_dependents in self._variable_dependencies:
            reference = linear_dependents["reference"]

            if self._variable_lookup[reference]["property"] == "fluid":
                continue

            all_containers = [
                self._variable_lookup[var]["object"].get_attr(
                    self._variable_lookup[var]["property"]
                ) for var in linear_dependents["variables"]
            ]

            number_specifications = sum(
                [not c._potential_var for c in all_containers]
            )
            if number_specifications > 1:
                variables_properties = [
                    f"({self._variable_lookup[var]['object'].label}: "
                    f"{self._variable_lookup[var]['property']})"
                    for var in linear_dependents["variables"]
                ]
                msg = (
                    "You specified more than one variable of the linear "
                    f"dependent variables: {', '.join(variables_properties)}."
                )
                raise hlp.TESPyNetworkError(msg)
            elif number_specifications == 1:
                reference_data = self._variable_lookup[reference]
                reference_container = reference_data["object"].get_attr(
                    reference_data["property"]
                )._reference_container
                reference_container.is_var = False

    def _assign_variable_space(self, reference, represents):
        container = reference["object"].get_attr(reference["property"])._reference_container
        if container.is_var:
            container.J_col = self.variable_counter
            self.variables_dict[self.variable_counter] = {
                "obj": container,
                "variable": reference["property"],
                "fluid": None,
                "_represents": represents
            }
            self.variable_counter += 1

    def _init_set_properties(self):
        """Specification of SI values for user set values."""
        # fluid property values
        for c in self.conns['object']:

            if not self.init_previous:
                c.good_starting_values = False

            for key in c.property_data:
                if "fluid" in key:
                    continue

                param = c.get_attr(key)
                if param.is_set:
                    if "ref" in key:
                        unit = self.units.default[param.quantity]
                        param.ref.delta_SI = self.units.ureg.Quantity(
                            param.ref.delta,
                            unit
                        ).to(SI_UNITS[param.quantity]).magnitude
                    else:
                        param.set_SI_from_val(self.units)
        msg = (
            "Updated fluid property SI values and fluid mass fraction for user "
            "specified connection parameters."
        )
        logger.debug(msg)

        for cp in self.comps["object"]:
            for param, value in cp.parameters.items():
                if isinstance(value, dc_prop) and (value.is_set or value.is_var):
                    if np.isnan(value._val):
                        value.val = (value.min_val + value.max_val) / 2
                    value.set_SI_from_val(self.units)


    def _init_design(self):
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
                path = c.design_path
                if path is None:
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

                if path not in _local_designs:
                    _local_designs[path] = self._load_network_state(path)

                df = _local_designs[c.design_path][c.__class__.__name__]
                # write data to connections
                self._init_conn_design_params(c, df)

            else:
                c._reset_design(self.redesign)
                # unset all design values
        # unset design values for busses, count bus equations and
        # reindex bus dictionary
        for b in self.busses.values():
            b.comps['P_ref'] = np.nan

        series = pd.Series(dtype='float64')
        for cp in self.comps['object']:
            c = cp.__class__.__name__
            # read design point information of components with
            # local_offdesign activated from their respective design path
            if cp.local_offdesign:
                path = cp.design_path
                if path is None:
                    msg = (
                        "The parameter local_offdesign is True for the "
                        f"component {cp.label}, an individual design_path must "
                        "be specified in this case!"
                    )
                    logger.error(msg)
                    raise hlp.TESPyNetworkError(msg)
                if path not in _local_designs:
                    _local_designs[path] = self._load_network_state(path)

                data = _local_designs[path][c]
                # write data
                self._init_comp_design_params(cp, data)

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

                cp.new_design = True

            else:
                # switch connections to design mode
                if self.redesign:
                    for var in cp.design:
                        cp.get_attr(var).is_set = True

                    for var in cp.offdesign:
                        cp.get_attr(var).is_set = False

                cp._set_design_parameters(self.mode, series)

    def _init_offdesign(self):
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
                    param.set_val_from_SI(self.units)

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
                        data.val_SI = data.design
                        data.set_val_from_SI(self.units)
                        switched = True
                        msg += var + ', '

                if switched:
                    msg = f"{msg[:-2]} to design value at component {cp.label}."
                    logger.debug(msg)

                cp.new_design = False

        msg = 'Switched components from design to offdesign.'
        logger.debug(msg)

    def _init_offdesign_params(self):
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
        # components with offdesign parameters
        components_with_parameters = [
            cp.label for cp in self.comps["object"] if len(cp.parameters) > 0
        ]
        # fetch all components, reindex with label
        df_comps = self.comps.loc[components_with_parameters].copy()
        # iter through unique types of components (class names)
        dfs = self._load_network_state(self.design_path)
        # iter through all components of this type and set data
        ind_designs = {}
        for label, row in df_comps.iterrows():
            df = dfs[row["comp_type"]]
            comp = row["object"]
            path = comp.design_path
            # read data of components with individual design_path
            if path is not None:
                if path not in ind_designs:
                    ind_designs[path] = self._load_network_state(path)
                data = ind_designs[path][row["comp_type"]]

            else:
                data = df

            # write data to components
            self._init_comp_design_params(comp, data)

        msg = 'Done reading design point information for components.'
        logger.debug(msg)

        if len(self.busses) > 0:
            for b, bus in self.busses.items():
                # the bus design data are stored in dfs[b][0] (column is not named)
                if len(bus.comps) > 0:
                    bus.comps.loc[self.get_comp(dfs[b].index), "P_ref"] = dfs[b][0].values

        # iter through connections
        for c in self.conns['object']:
            conn_type = c.__class__.__name__
            df = dfs[conn_type]
            # read data of connections with individual design_path
            path = c.design_path
            if path is not None:
                if path not in ind_designs:
                    ind_designs[path] = self._load_network_state(path)
                data = ind_designs[path][conn_type]
            else:
                data = df

            self._init_conn_design_params(c, data)

        msg = 'Done reading design point information for connections.'
        logger.debug(msg)

    def _init_comp_design_params(self, c, df):
        r"""
        Write design point information to components.

        Parameters
        ----------
        component : tespy.components.component.Component
            Write design point information to this component.

        data : pandas.core.series.Series, pandas.core.frame.DataFrame
            Design point information.
        """
        if c.label not in df.index:
            # no matches in the connections of the network and the design files
            msg = (
                f"Could not find component '{c.label}' in design case file. "
                "This is is critical only to components, which need to load "
                "design values from this case."
            )
            logger.debug(msg)
            return
        # write component design data
        data = df.loc[c.label]
        c._set_design_parameters(self.mode, data)

    def _init_conn_design_params(self, c, df):
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
                f"Could not find connection '{c.label}' in design case. Please "
                "make sure no connections have been modified or components "
                "have been relabeled for your offdesign calculation."
            )
            logger.exception(msg)
            raise hlp.TESPyNetworkError(msg)

        data = df.loc[c.label]
        c._set_design_params(data, self.units)

    def _init_conn_params_from_path(self, c, df):
        r"""
        Write parameter information from init_path to connections.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            Write init path information to this connection.

        df : pandas.core.frame.DataFrame
            Dataframe containing init path information.
        """
        if c.label not in df.index:
            # no matches in the connections of the network and the design files
            msg = f"Could not find connection {c.label} in init path file."
            logger.debug(msg)
            return

        data = df.loc[c.label]
        c._set_starting_values(data, self.units)
        c.good_starting_values = True

    def _init_properties(self):
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
            dfs = self._load_network_state(self.init_path)
        # improved starting values for referenced connections,
        # specified vapour content values, temperature values as well as
        # subccooling/overheating and state specification
        for c in self.conns['object']:
            if self.init_path is not None:
                df = dfs[c.__class__.__name__]
                self._init_conn_params_from_path(c, df)

            if type(c) == Connection:
                # the below part does not work for PowerConnection right now
                if sum(c.fluid.val.values()) == 0:
                    msg = (
                        'The starting value for the fluid composition of the '
                        f'connection {c.label} is empty. This might lead to issues '
                        'in the initialisation and solving process as fluid '
                        'property functions can not be called. Make sure you '
                        'specified a fluid composition in all parts of the network.'
                    )
                    logger.warning(msg)

            for key, variable in c.get_variables().items():
                # for connections variables can be presolved and not be var anymore
                if variable.is_var:
                    if not c.good_starting_values:
                        self.init_val0(c, key)

                    variable.set_SI_from_val0(self.units)
                    # variable.set_SI_from_val0()
                    variable.set_reference_val_SI(variable._val_SI)

        for cp in self.comps["object"]:
            for key, variable in cp.get_variables().items():
                # for components every variable should be an actual variable
                # if variable.is_var:
                if np.isnan(variable.val):
                    variable.val = 1.0
                variable.set_SI_from_val(self.units)
                variable.set_reference_val_SI(variable._val_SI)

        for c in self.conns['object']:
            c._precalc_guess_values()

        msg = 'Generic fluid property specification complete.'
        logger.debug(msg)

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
                seed = abs(hash(c.label)) % (2**32)
                rng = np.random.default_rng(seed=seed)
                value = float(rng.random() + 1)

            # generic starting values for pressure and enthalpy
            elif key in ['p', 'h']:
                # retrieve starting values from component information
                val_s = c.source.initialise_source(c, key)
                val_t = c.target.initialise_target(c, key)

                if val_s == 0 and val_t == 0:
                    if key == 'p':
                        value = 1e5
                    elif key == 'h':
                        value = 1e6

                elif val_s == 0:
                    value = val_t
                elif val_t == 0:
                    value = val_s
                else:
                    value = (val_s + val_t) / 2

            elif key == 'E':
                value = 0.0

            # these values are SI, so they are set to the respective variable
            c.get_attr(key).set_reference_val_SI(value)
            c.get_attr(key).set_val0_from_SI(self.units)

    @staticmethod
    def _load_network_state(json_path):
        r"""
        Read network state from given file.

        Parameters
        ----------
        json_path : str
            Path to network information.
        """
        with open(json_path, "r") as f:
            data = json.load(f)

        dfs = {}
        if "Connection" in data["Connection"]:
            for key, value in data["Connection"].items():
                with pd.option_context("future.no_silent_downcasting", True):
                    dfs[key] = pd.DataFrame.from_dict(value, orient="index").fillna(np.nan)
                dfs[key].index = dfs[key].index.astype(str)
        # TODO: depricate
        # this is for compatibility of older savestates
        else:
            key = "Connection"
            value = data["Connection"]
            with pd.option_context("future.no_silent_downcasting", True):
                dfs[key] = pd.DataFrame.from_dict(value, orient="index").fillna(np.nan)
            dfs[key].index = dfs[key].index.astype(str)

        for key, value in data["Component"].items():
            with pd.option_context("future.no_silent_downcasting", True):
                dfs[key] = pd.DataFrame.from_dict(value, orient="index").fillna(np.nan)
            dfs[key].index = dfs[key].index.astype(str)
        for key, value in data["Bus"].items():
            with pd.option_context("future.no_silent_downcasting", True):
                dfs[key] = pd.DataFrame.from_dict(value, orient="index").fillna(np.nan)
            dfs[key].index = dfs[key].index.astype(str)

        return dfs

    def get_linear_dependent_variables(self) -> list:
        """Get a list with sublists containing linear dependent variables

        Returns
        -------
        list
            List of lists of linear dependent variables
        """
        variable_list = []
        for dependents in self._variable_dependencies:
            variables = [
                self._variable_lookup[v] for v in dependents["variables"]
            ]
            variable_list += [
                [(v["object"].label, v["property"]) for v in variables]
            ]
        return variable_list

    def _get_equation_sets_by_eq_set_number(self, number_list) -> list:
        return [self._equation_set_lookup[num] for num in number_list]

    def _get_variables_before_presolve_by_number(self, number_list) -> list:
        return [
            (v["object"].label, v["property"])
            for k, v in self._variable_lookup.items()
            if k in number_list
        ]

    def get_presolved_equations(self) -> list:
        """Get the list of equations, that has been presolved with their
        respective parent object

        Returns
        -------
        list
            list of presolved equations
        """
        return [
            v for k, v in self._equation_set_lookup.items()
            if k in self._presolved_equations
        ]

    def get_variables_before_presolve(self) -> list:
        """Get the list of variables before presolving.

        Returns
        -------
        list
            list of original variables
        """
        return [
            (v["object"].label, v["property"])
            for v in self._variable_lookup.values()
        ]

    def get_presolved_variables(self) -> list:
        """Get the list of presolved variables with their respective parent
        object and property.

        Returns
        -------
        list
            list of presolved variables
        """
        represented_variables = []
        for v in self.variables_dict.values():
            represented_variables += v["_represents"]
        if len(self.variables_dict) == 0 and len(self._presolved_equations) == 0:
            return []
        return [
            (v["object"].label, v["property"])
            for key, v in self._variable_lookup.items()
            if key not in represented_variables
        ]

    def get_variables(self) -> dict:
        """Get all variables of the presolved problem with their respective
        represented original variables.

        Returns
        -------
        dict
            variable number and property with the list of represented variables
        """
        return {
            (key, data["variable"]):
            [
                (
                    self._variable_lookup[v]["object"].label,
                    self._variable_lookup[v]["property"]
                ) for v in data["_represents"]
            ]
            for key, data in self.variables_dict.items()
        }

    def _get_variables_by_number(self, number_list) -> dict:
        """Get all variables of the presolved problem by variable numbers.

        Returns
        -------
        dict
            variable number and property with the list of represented variables
        """
        return {
            (key, data["variable"]):
            [
                (
                    self._variable_lookup[v]["object"].label,
                    self._variable_lookup[v]["property"]
                ) for v in data["_represents"]
            ]
            for key, data in self.variables_dict.items()
            if key in number_list
        }

    def get_equations(self) -> dict:
        """Get the actual equations after presolving the problem

        Returns
        -------
        dict
            Lookup with equation number as index and tuple of label and
            parameter defining the equation. In case one parameter defines
            multiple equations, the same equation is repeated.
        """
        return self._equation_lookup

    def get_equations_with_dependents(self) -> dict:
        """Get the equations together with the variables they depend on.

        Returns
        -------
        dict
            Lookup with equation (component, (parameter_label, number)) and
            the variables it depends on as a list
            (variable number, variable type)
        """
        dependencies = {}
        for eq_idx, dependents in self._incidence_matrix.items():
            dependencies.update({
                self._equation_lookup[eq_idx]:
                list(self._get_variables_by_number(dependents).keys())
            })
        return dependencies


    def _get_equations_by_number(self, number_list) -> dict:
        """Get the actual equations after presolving the problem by equation
        number

        Returns
        -------
        dict
            Lookup with equation number as index and tuple of label and
            parameter defining the equation. In case one parameter defines
            multiple equations, the same equation is repeated.
        """
        return {
            k: v for k, v in self._equation_lookup.items() if k in number_list
        }

    def get_linear_dependents_by_object(self, obj, prop) -> list:
        """Get the list of linear dependent variables for a specified variable

        Parameters
        ----------
        obj : object
            Parent object holding a variable
        prop : str
            Name of the variable (e.g. 'm' or 'h')

        Returns
        -------
        list
            list of linear dependent variables

        Raises
        ------
        KeyError
            In case the object does not have any variables
        KeyError
            In case the specified property is not a variable
        """
        if obj not in self._object_to_variable_lookup:
            msg = f"The object {obj.label} does not have any variables."
            raise KeyError(msg)

        if prop not in self._object_to_variable_lookup[obj]:
            msg = f"The object {obj.label} does not have a variable {prop}."
            raise KeyError(msg)

        variable_idx = self._object_to_variable_lookup[obj][prop]
        return self._get_linear_dependents_by_variable_index(variable_idx)

    def _get_linear_dependents_by_variable_index(self, idx) -> list:
        """Get the list of linear dependent variables for a specified variable

        Parameters
        ----------
        idx : object
            Index of the variable

        Returns
        -------
        list
            list of linear dependent variables
        """
        for dependents in self._variable_dependencies:
            if idx in dependents["variables"]:
                break
        variables = [self._variable_lookup[v] for v in dependents["variables"]]
        variable_list = [(v["object"].label, v["property"]) for v in variables]
        return variable_list

    def solve(self, mode, init_path=None, design_path=None,
              max_iter=50, min_iter=4, init_only=False, init_previous=True,
              use_cuda=False, print_results=True, robust_relax=False):
        r"""
        Solve the network.

        - Check network consistency.
        - Initialise calculation and preprocessing.
        - Perform actual calculation.
        - Postprocessing.

        It is possible to check programatically, if a network was solved
        successfully with the `.converged` attribute.

        Parameters
        ----------
        mode : str
            Choose from 'design' and 'offdesign'.

        init_path : str
            Path to the folder, where your network was saved to, e.g.
            saving to :code:`nw.save('myplant/test.json')` would require loading
            from :code:`init_path='myplant/test.json'`.

        design_path : str
            Path to the folder, where your network's design case was saved to,
            e.g. saving to :code:`nw.save('myplant/test.json')` would require
            loading from :code:`design_path='myplant/test.json'`.

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
        self.status = 99
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

        self.init_path = init_path
        self.design_path = design_path
        self.max_iter = max_iter
        self.min_iter = min_iter
        self.init_previous = init_previous
        self.iter = 0
        self.use_cuda = use_cuda
        self.robust_relax = robust_relax

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
            self.check_topology()

        msg = (
            "Solver properties:\n"
            f" - mode: {self.mode}\n"
            f" - init_path: {self.init_path}\n"
            f" - design_path: {self.design_path}\n"
            f" - min_iter: {self.min_iter}\n"
            f" - max_iter: {self.max_iter}"
        )
        logger.debug(msg)

        msg = (
            "Network information:\n"
            f" - Number of components: {len(self.comps)}\n"
            f" - Number of connections: {len(self.conns)}\n"
            f" - Number of busses: {len(self.busses)}"
        )
        logger.debug(msg)

        self._prepare_problem()

        if init_only:
            return

        msg = 'Starting solver.'
        logger.info(msg)

        self.solve_determination()

        try:
            self.solve_loop(print_results=print_results)
        except ValueError as e:
            self.status = 99
            msg = f"Simulation crashed due to an unexpected error:\n{e}"
            logger.exception(msg)
            self.unload_variables()
            return

        self.unload_variables()

        if self.status == 3:
            logger.error(self.singularity_msg)
            return

        if self.status == 2:
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
        self.residual = np.zeros([self.variable_counter])
        self.increment = np.ones([self.variable_counter])
        self.jacobian = np.zeros((self.variable_counter, self.variable_counter))

        self.start_time = time()

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

            if self.lin_dep:
                self.status = 3
                break

            elif self.iter > 40:
                if (
                    all(
                        self.residual_history[(self.iter - 3):] >= self.residual_history[-3] * 0.95
                    ) and self.residual_history[-1] >= self.residual_history[-2] * 0.95
                ):
                    self.status = 2
                    break

            if (
                    self.iter >= self.min_iter - 1
                    and (self.residual_history[-2:] < ERR ** 0.5).all()
                    # the increment should also be small
                    and (abs(self.increment) < ERR ** 0.5).all()
                ):
                self.status = 0
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
            self.status = 2

    def solve_determination(self):
        r"""Check, if the number of supplied parameters is sufficient."""
        msg = f'Number of connection equations: {self.num_conn_eq}.'
        logger.debug(msg)
        msg = f'Number of bus equations: {self.num_bus_eq}.'
        logger.debug(msg)
        msg = f'Number of component equations: {self.num_comp_eq}.'
        logger.debug(msg)
        msg = f'Number of user defined equations: {self.num_ude_eq}.'
        logger.debug(msg)

        msg = f'Total number of variables: {self.variable_counter}.'
        logger.debug(msg)

        n = (
            self.num_comp_eq + self.num_conn_eq +
            self.num_bus_eq + self.num_ude_eq
        )
        if n > self.variable_counter:
            msg = (
                f"You have provided too many parameters: {self.variable_counter} "
                f"required, {n} supplied. Aborting calculation!"
            )
            logger.error(msg)
            self.status = 12
            raise hlp.TESPyNetworkError(msg)
        elif n < self.variable_counter:
            msg = (
                f"You have not provided enough parameters: {self.variable_counter} "
                f"required, {n} supplied. Aborting calculation!"
            )
            logger.error(msg)
            self.status = 11
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

    def _invert_jacobian(self):
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

        n = self.variable_counter
        self._incidence_matrix_dense = np.zeros((n, n))
        for row, cols in self._incidence_matrix.items():
            self._incidence_matrix_dense[row, cols] = 1

        if self.lin_dep:
            if self.iter == 0 and np.linalg.det(self._incidence_matrix_dense) == 0.0:
                self.singularity_msg = (
                    "Detected singularity in Jacobian matrix. This singularity "
                    "is most likely caused by the parametrization of your "
                    "problem and NOT a numerical issue. Double check your "
                    "setup.\n"
                )
                self._find_linear_dependencies(self.jacobian)
                return

            # here we can analyse the same things as above but on the
            # Jacobian to give hints what equations/variables are causing
            # the issue.

            # On top it is possible to compare the contents
            # of the Jacobian with the contents of the incidence matrix
            # to maybe pinpoint in more detail, what causes the trouble
            expected_entries = self._incidence_matrix_dense.astype(bool)
            actual_entries = self.jacobian.astype(bool)
            rows, cols = np.where(expected_entries != actual_entries)

            missing_entries = []
            for row, col in zip(rows, cols):
                equation = self._get_equations_by_number([row])
                variable = self._get_variables_by_number([col])
                missing_entries += [f"{equation}: {variable}"]

            _nl = "\n"
            self.singularity_msg = (
                "Found singularity in Jacobian matrix, calculation "
                "aborted! The setup of you problem seems to be solvable. It "
                "failed due to partial derivatives in the Jacobian being "
                "zero, which were expected not to be zero, or the other way "
                "around. The reason for this usually lies in starting value "
                "selection or bad convergence. The following equations (key "
                "of outer dict) may have an unexpected zero/non-zero in the "
                "partial derivative towards the variable (value of outer "
                f"dict) and be the root of evil: {_nl.join(missing_entries)}"
            )
            self._find_linear_dependencies(self.jacobian)
            return

    def _find_linear_dependencies(self, matrix):
        _nl = "\n"
        all_zero_cols = self._check_all_zero_columns(matrix)
        all_zero_rows = self._check_all_zero_rows(matrix)
        if len(all_zero_cols) + len(all_zero_rows) == 0:
            equations = self._cauchy_schwarz_inequality(matrix)
            equations = self._get_equations_by_number(equations)
            self.singularity_msg += (
                "The following equations form a linear dependency in "
                "the : "
                f"{', '.join([str(e) for e in equations.values()])}{_nl}"
            )
        else:
            if len(all_zero_cols) > 0:
                variables = self._get_variables_by_number(all_zero_cols)
                self.singularity_msg += (
                    "The following variables of your problem are not "
                    "in connection with any equation: "
                    f"{', '.join([str(v) for v in variables])}{_nl}"
                )
            if len(all_zero_rows) > 0:
                equations = self._get_equations_by_number(all_zero_rows)
                self.singularity_msg += (
                    "The following equations of your problem do not "
                    "depend on any variable: "
                    f"{', '.join([str(e) for e in equations.values()])}{_nl}"
                )

    def _check_all_zero_columns(self, matrix):
        return np.where((matrix == 0).all(axis=0))[0]

    def _check_all_zero_rows(self, matrix):
        return np.where((matrix == 0).all(axis=1))[0]

    def _cauchy_schwarz_inequality(self, matrix):
        n = matrix.shape[0]
        dependent_equations = []
        for i in range(n):
            for j in range(n):
                if i != j:
                    inner_product = np.inner(
                        matrix[i,:],
                        matrix[j,:]
                    )
                    norm_i = np.linalg.norm(matrix[i,:])
                    norm_j = np.linalg.norm(matrix[j,:])

                    if np.abs(inner_product - norm_j * norm_i) < 1e-5:
                        dependent_equations += [i]
        return list(set(dependent_equations))

    def _update_variables(self):
        # cast dtype to float from numpy float64
        # this is necessary to keep the doctests running and note make them
        # look ugly all over the place
        # I have yet to come up with a better idea, or vectorize all operations
        # which requires major changes in tespy
        increment = [float(val) for val in self.increment]
        # the J_cols here point to actual variables, no need to call to
        # get_J_col yet
        relax = 1
        if self.robust_relax:
            if self.iter < 3:
                relax = 0.25
            elif self.iter < 5:
                relax = 0.5
            elif self.iter < 8:
                relax = 0.75

        for _, data in self.variables_dict.items():
            if data["variable"] in ["m", "h", "E"]:
                container = data["obj"]
                container._val_SI += increment[container.J_col] * relax
            elif data["variable"] == "p":
                container = data["obj"]
                relax = max(
                    1, -2 * increment[container.J_col] / container.val_SI
                )
                container._val_SI += increment[container.J_col] / relax
            elif data["variable"] == "fluid":
                container = data["obj"]
                container.val[data["fluid"]] += increment[
                    container.J_col[data["fluid"]]
                ]

                if container.val[data["fluid"]] < ERR :
                    container.val[data["fluid"]] = 0
                elif container.val[data["fluid"]] > 1 - ERR :
                    container.val[data["fluid"]] = 1
            else:
                # add increment
                data["obj"]._val_SI += increment[data["obj"].J_col] * relax

                # keep value within specified value range
                if data["obj"].val_SI < data["obj"].min_val:
                    data["obj"].val_SI = data["obj"].min_val
                elif data["obj"].val_SI > data["obj"].max_val:
                    data["obj"].val_SI = data["obj"].max_val

    def check_variable_bounds(self):

        # this could be in a different place, its kind of in between
        # network and connection
        for idx, data in self.variables_dict.items():
            if type(data["obj"]) == dc_vecvar:
                total_mass_fractions = sum(data["obj"].val.values())
                for fluid in data["obj"].is_var:
                    data["obj"]._val[fluid] /= total_mass_fractions

        if norm(self.increment) > 1e-1:
            for c in self.conns['object']:
                # check the fluid properties for physical ranges
                c._adjust_to_property_limits(self)

        # second check based on component heuristics
        # - for first three iterations
        # - only if the increment is sufficiently large
        # - only in design case
        if (
                self.iter < 3
                and norm(self.increment) > 1e3
                and self.mode == "design"
            ):
            for cp in self.comps['object']:
                cp.convergence_check()

            for c in self.conns['object']:
                c._adjust_to_property_limits(self)

    def solve_control(self):
        r"""
        Control iteration step of the newton algorithm.

        - Calculate the residual value for each equation
        - Calculate the jacobian matrix
        - Calculate new values for variables
        - Restrict fluid properties to value ranges
        - Check component parameters for consistency
        """
        self.solve_equations()
        self.solve_busses()
        self._invert_jacobian()

        # check for linear dependency
        if self.lin_dep:
            return

        self._update_variables()
        self.check_variable_bounds()

    def solve_equations(self):
        r"""
        Calculate the residual and derivatives of all equations.
        """
        to_solve = (
            self.comps["object"].tolist()
            + self.conns["object"].tolist()
            + list(self.user_defined_eq.values())
        )
        for obj in to_solve:
            hlp.solve(obj, self.increment_filter)
            if len(obj.jacobian) > 0:
                rows = list(obj.residual.keys())
                data = list(obj.residual.values())
                self.residual[rows] = data

                rows = [k[0] for k in obj.jacobian]
                columns = [k[1] for k in obj.jacobian]
                data = list(obj.jacobian.values())
                self.jacobian[rows, columns] = data

            obj.it += 1

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
        _converged = self.process_connections()
        _converged = self.process_components() and _converged
        self.process_busses()

        if self.status == 0 and not _converged:
            self.status = 1

        msg = 'Postprocessing complete.'
        logger.info(msg)

    def unload_variables(self):
        for dependents in self._variable_dependencies:
            for variable_num in dependents["variables"]:
                variable_dict = self._variable_lookup[variable_num]
                variable = variable_dict["object"].get_attr(variable_dict["property"])
                if variable_dict["property"] != "fluid":
                    variable.val_SI = variable.val_SI
                else:
                    variable.val = variable.val
                variable._reference_container = None

    def process_connections(self):
        """Process the Connection results."""
        _converged = True
        for c in self.conns['object']:
            c.good_starting_values = True
            _converged = c.calc_results(self.units) and _converged
            self.results[c.__class__.__name__].loc[c.label] = c.collect_results(self.all_fluids)
        return _converged

    def process_components(self):
        """Process the component results."""
        # components
        _converged = True
        for cp in self.comps['object']:
            cp.calc_parameters()
            _converged = _converged and cp.check_parameter_bounds()
            # this thing could be somewhere else
            for key, value in cp.parameters.items():
                if isinstance(value, dc_prop):
                    result = value._get_val_from_SI(self.units)
                    if (
                        value.is_set
                        and not value.is_var
                        and not np.isclose(result.magnitude, value.val, 1e-3, 1e-3)
                        and not cp.bypass
                    ):
                        _converged = False
                        msg = (
                            "The simulation converged but the calculated "
                            f"result {result} for the fixed input parameter "
                            f"{key} is not equal to the originally specified "
                            f"value: {value.val}. Usually, this can happen, "
                            "when a method internally manipulates the "
                            "associated equation during iteration in order to "
                            "allow progress in situations, when the equation "
                            "is otherwise not well defined for the current"
                            "values of the variables, e.g. in case a negative "
                            "root would need to be evaluated.  Often, this "
                            "can happen during the first iterations and then "
                            "will resolve itself as convergence progresses. "
                            "In this case it did not, meaning convergence was "
                            "not actually achieved."
                        )
                        logger.warning(msg)
                        self.status = 2
                    else:
                        if not value.is_set or value.is_var:
                            value.set_val_from_SI(self.units)

        if self.status == 2:
            return False

        for cp in self.comps['object']:
            key = cp.__class__.__name__
            result = cp.collect_results()
            if len(result) == 0:
                continue
            self.results[cp.__class__.__name__].loc[cp.label] = result

        return _converged

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

    def print_results(self, colored=True, colors=None, print_results=True, subsystem=None):
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
                'only available after a full simulation run is performed.'
            )
            raise hlp.TESPyNetworkError(msg)

        result += self._print_components(colored, coloring, subsystem)
        result += self._print_connections(colored, coloring, subsystem)
        result += self._print_buses(colored, coloring, subsystem)

        if len(str(result)) > 0:
            logger.result(result)
            if print_results:
                print(result)
        return

    def _print_components(self, colored, coloring, subsystem) -> str:
        result = ""
        for cp in self.comps['comp_type'].unique():
            df = self.results[cp].copy()
            for c in df.index:
                if not self.get_comp(c).printout:
                    df = df.drop(c)
            # are there any parameters to print?
            if df.size > 0:
                if subsystem is not None:
                    component_labels = [
                        c.label for c in subsystem.comps.values()
                        if c.label in df.index
                    ]
                    df = df.loc[component_labels]

                c = self.comps.loc[self.comps["comp_type"] == cp, "object"]
                cols = c.iloc[0]._get_result_attributes()
                if len(cols) > 0:
                    for col in cols:
                        df[col] = df.apply(
                            self._color_component_prints, axis=1,
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
        return result

    def _print_connections(self, colored, coloring, subsystem) -> str:
        result = ""

        # connection properties
        for c_type in self.conns["conn_type"].unique():
            cols = connection_registry.items[c_type]._print_attributes()
            df = self.results[c_type].copy().loc[:, cols]

            if subsystem is not None:
                connection_labels = [c.label for c in subsystem.conns.values()]
                df = df.loc[connection_labels]

            df = df.astype(str)
            for c in df.index:
                if not self.get_conn(c).printout:
                    df.drop([c], axis=0, inplace=True)

                elif colored:
                    conn = self.get_conn(c)
                    for col in df.columns:
                        if conn.get_attr(col).is_set:
                            value = conn.get_attr(col).val
                            df.loc[c, col] = (
                                f"{coloring['set']}{value}{coloring['end']}"
                            )

            if len(df) > 0:
                result += (f'\n##### RESULTS ({c_type}) #####\n')
                result += (
                    tabulate(df, headers='keys', tablefmt='psql', floatfmt='.3e')
                )
        return result

    def _print_buses(self, colored, coloring, subsystem) -> str:
        result = ""
        # bus printout only if not subsystem is passed
        if subsystem is None:
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
                            value = df.loc['total', 'bus value']
                            df.loc['total', 'bus value'] = (
                                f"{coloring['set']}{value}{coloring['end']}"
                            )
                    result += f"\n##### RESULTS (Bus: {b.label}) #####\n"
                    result += (
                        tabulate(
                            df, headers='keys', tablefmt='psql',
                            floatfmt='.3e'
                        )
                    )
        return result

    def _color_component_prints(self, c, *args):
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
            param_obj = comp.get_attr(param)
            val = param_obj.val
            val_SI = param_obj.val_SI
            if not colored:
                return str(val)
            # else part
            if (val_SI < comp.get_attr(param).min_val - ERR or
                    val_SI > comp.get_attr(param).max_val + ERR ):
                return f"{coloring['err']}{val}{coloring['end']}"
            if comp.get_attr(args[0]).is_var:
                return f"{coloring['var']}{val}{coloring['end']}"
            if comp.get_attr(args[0]).is_set:
                return f"{coloring['set']}{val}{coloring['end']}"
            return str(val)
        else:
            return np.nan

    @classmethod
    def from_json(cls, json_file_path):
        r"""
        Load a network from a base path.

        Parameters
        ----------
        path : str
            The path to the network data.

        Returns
        -------
        nw : tespy.networks.network.Network
            TESPy networks object.

        Note
        ----
        If you export the network structure of an existing TESPy network, it
        will be saved to the path you specified. The structure of the saved
        data in that path is the structure you need to provide in the path for
        loading the network.

        The structure of the path must be as follows:

        - Folder: path (e.g. 'mynetwork')
        - Component.json
        - Connection.json
        - Bus.json
        - Network.json

        Example
        -------
        Create a network and export it. This is followed by loading the network
        from the exported json file. All network information stored will be
        passed to a new network object. Components, connections and busses will
        be accessible by label. The following example setup is simple gas
        turbine setup with compressor, combustion chamber and turbine. The fuel
        is fed from a pipeline and throttled to the required pressure while
        keeping the temperature at a constant value.

        >>> from tespy.components import (
        ...     Sink, Source, CombustionChamber, Compressor, Turbine,
        ...     SimpleHeatExchanger, PowerBus, PowerSink, Generator
        ... )
        >>> from tespy.connections import Connection, Ref, PowerConnection
        >>> from tespy.networks import Network
        >>> import os
        >>> nw = Network(iterinfo=False)
        >>> nw.units.set_defaults(**{
        ...     "pressure": "bar", "temperature": "degC", "enthalpy": "kJ/kg",
        ...     "power": "MW"
        ... })
        >>> air = Source('air')
        >>> f = Source('fuel')
        >>> compressor = Compressor('compressor')
        >>> combustion = CombustionChamber('combustion')
        >>> turbine = Turbine('turbine')
        >>> preheater = SimpleHeatExchanger('fuel preheater')
        >>> si = Sink('sink')
        >>> shaft = PowerBus('shaft', num_in=1, num_out=2)
        >>> generator = Generator('generator')
        >>> grid = PowerSink('grid')
        >>> c1 = Connection(air, 'out1', compressor, 'in1', label='c01')
        >>> c2 = Connection(compressor, 'out1', combustion, 'in1', label='c02')
        >>> c11 = Connection(f, 'out1', preheater, 'in1', label='c11')
        >>> c12 = Connection(preheater, 'out1', combustion, 'in2', label='c12')
        >>> c3 = Connection(combustion, 'out1', turbine, 'in1', label='c03')
        >>> c4 = Connection(turbine, 'out1', si, 'in1', label='c04')
        >>> nw.add_conns(c1, c2, c11, c12, c3, c4)
        >>> e1 = PowerConnection(turbine, 'power', shaft, 'power_in1', label='e1')
        >>> e2 = PowerConnection(shaft, 'power_out1', compressor, 'power', label='e2')
        >>> e3 = PowerConnection(shaft, 'power_out2', generator, 'power_in', label='e3')
        >>> e4 = PowerConnection(generator, 'power_out', grid, 'power', label='e4')
        >>> nw.add_conns(e1, e2, e3, e4)

        Specify component and connection properties. The intlet pressure at the
        compressor and the outlet pressure after the turbine are identical. For
        the compressor, the pressure ratio and isentropic efficiency are design
        parameters. A compressor map (efficiency vs. mass flow and pressure
        rise vs. mass flow) is selected for the compressor. Fuel is Methane.

        >>> compressor.set_attr(
        ...     pr=10, eta_s=0.88, design=['eta_s', 'pr'],
        ...     offdesign=['char_map_eta_s', 'char_map_pr']
        ... )
        >>> turbine.set_attr(
        ...     eta_s=0.9, design=['eta_s'],
        ...     offdesign=['eta_s_char', 'cone']
        ... )
        >>> combustion.set_attr(lamb=2)
        >>> c1.set_attr(
        ...     fluid={'N2': 0.7556, 'O2': 0.2315, 'Ar': 0.0129}, T=25, p=1
        ... )
        >>> c11.set_attr(fluid={'CH4': 0.96, 'CO2': 0.04}, T=25, p=40)
        >>> c12.set_attr(T=25)
        >>> c4.set_attr(p=Ref(c1, 1, 0))
        >>> generator.set_attr(eta=1)

        For a stable start, we specify the fresh air mass flow.

        >>> c1.set_attr(m=3)
        >>> nw.solve('design')
        >>> nw.assert_convergence()

        The total power output is set to 1 MW, electrical or mechanical
        efficiencies are not considered in this example. The documentation
        example in class :py:class:`tespy.connections.bus.Bus` provides more
        information on efficiencies of generators, for instance.

        >>> combustion.set_attr(lamb=None)
        >>> c3.set_attr(T=1100)
        >>> c1.set_attr(m=None)
        >>> e4.set_attr(E=1)
        >>> nw.solve('design')
        >>> nw.assert_convergence()
        >>> nw.save('design_state.json')
        >>> _ = nw.export('exported_nwk.json')
        >>> mass_flow = round(nw.get_conn('c01').m.val_SI, 1)
        >>> compressor.set_attr(igva='var')
        >>> nw.solve('offdesign', design_path='design_state.json')
        >>> round(turbine.eta_s.val, 1)
        0.9
        >>> e4.set_attr(E=0.75)
        >>> nw.solve('offdesign', design_path='design_state.json')
        >>> nw.assert_convergence()
        >>> eta_s_t = round(turbine.eta_s.val, 3)
        >>> igva = round(compressor.igva.val, 3)
        >>> eta_s_t
        0.898
        >>> igva
        20.138

        The designed network is exported to the path 'exported_nwk'. Now import
        the network and recalculate. Check if the results match with the
        previous calculation in design and offdesign case.

        >>> imported_nwk = Network.from_json('exported_nwk.json')
        >>> imported_nwk.set_attr(iterinfo=False)
        >>> imported_nwk.solve('design')
        >>> imported_nwk.lin_dep
        False
        >>> round(imported_nwk.get_conn('c01').m.val_SI, 1) == mass_flow
        True
        >>> round(imported_nwk.get_comp('turbine').eta_s.val, 3)
        0.9
        >>> imported_nwk.get_comp('compressor').set_attr(igva='var')
        >>> imported_nwk.solve('offdesign', design_path='design_state.json')
        >>> round(imported_nwk.get_comp('turbine').eta_s.val, 3)
        0.9
        >>> imported_nwk.get_conn('e4').set_attr(E=0.75)
        >>> imported_nwk.solve('offdesign', design_path='design_state.json')
        >>> round(imported_nwk.get_comp('turbine').eta_s.val, 3) == eta_s_t
        True
        >>> round(imported_nwk.get_comp('compressor').igva.val, 3) == igva
        True
        >>> os.remove('exported_nwk.json')
        >>> os.remove('design_state.json')
        """
        msg = f'Reading network data from base path {json_file_path}.'
        logger.info(msg)

        with open(json_file_path, "r") as f:
            network_data = json.load(f)
        # create network
        # get method to ensure compatibility with old style export
        units = Units.from_json(network_data["Network"].get("units", {}))
        network_data["Network"]["units"] = units
        nw = cls(**network_data["Network"])

        # load components
        comps = {}

        module_name = "tespy.components"
        _ = importlib.import_module(module_name)

        for component, data in network_data["Component"].items():
            if component not in component_registry.items:
                msg = (
                    f"A class {component} is not available through the "
                    "tespy.components.component.component_registry decorator. "
                    "If you are using a custom component make sure to "
                    "decorate the class."
                )
                logger.error(msg)
                raise hlp.TESPyNetworkError(msg)

            target_class = component_registry.items[component]
            comps.update(_construct_components(target_class, data, nw))

        msg = 'Created network components.'
        logger.info(msg)

        conns = {}
        # load connections
        if "Connection" not in network_data["Connection"]:
            # v0.8 compatibility
            target_class = connection_registry.items["Connection"]
            conns.update(_construct_connections(
                target_class, network_data["Connection"], comps)
            )
        else:
            for connection, data in network_data["Connection"].items():
                if connection not in connection_registry.items:
                    msg = (
                        f"A class {connection} is not available through the "
                        "tespy.connections.connection.connection_registry "
                        "decorator. If you are using a custom connection make "
                        "sure to decorate the class."
                    )
                    logger.error(msg)
                    raise hlp.TESPyNetworkError(msg)

                target_class = connection_registry.items[connection]
                conns.update(_construct_connections(target_class, data, comps))

        # add connections to network
        for c in conns.values():
            nw.add_conns(c)

        msg = 'Created connections.'
        logger.info(msg)

        # load busses
        data = network_data.get("Bus", {})
        if len(data) > 0:
            busses = _construct_busses(data, comps)
            # add busses to network
            for b in busses.values():
                nw.add_busses(b)

            msg = 'Created busses.'
            logger.info(msg)

        else:
            msg = 'No bus data found!'
            logger.debug(msg)

        msg = 'Created network.'
        logger.info(msg)

        nw.check_topology()

        return nw

    def export(self, json_file_path=None):
        """Export the parametrization and structure of the Network instance

        Parameters
        ----------
        json_file_path : str, optional
            Path for exporting to filesystem. If path is None, the data are
            only returned and not written to the filesystem, by default None.

        Returns
        -------
        dict
            Parametrization and structure of the Network instance.
        """
        export = {}
        export["Network"] = self._export_network()
        export["Connection"] = self._export_connections()
        export["Component"] = self._export_components()
        export["Bus"] = self._export_busses()

        if json_file_path:
            with open(json_file_path, "w") as f:
                json.dump(export, f, indent=2)

            logger.debug(f'Model information saved to {json_file_path}.')

        return export

    def to_exerpy(self, Tamb, pamb, exerpy_mappings):
        """Export the network to exerpy

        Parameters
        ----------
        Tamb : float
            Ambient temperature.
        pamb : float
            Ambient pressure.
        exerpy_mappings : dict
            Mappings for tespy components to exerpy components

        Returns
        -------
        dict
            exerpy compatible input dictionary
        """
        component_results = self._save_components()
        component_json = {}
        for comp_type in self.comps["comp_type"].unique():
            if comp_type not in exerpy_mappings.keys():
                msg = f"Component class {comp_type} not available in exerpy."
                logger.warning(msg)
                continue

            key = exerpy_mappings[comp_type]
            if key not in component_json:
                component_json[key] = {}

            result = component_results[comp_type].dropna(axis=1)

            for c in self.comps.loc[self.comps["comp_type"] == comp_type, "object"]:
                parameters = {}
                if c.label in result.index and not result.loc[c.label].dropna().empty:
                    parameters = result.loc[c.label].dropna().to_dict()
                component_json[key][c.label] = {
                    "name": c.label,
                    "type": comp_type,
                    "parameters": parameters
                }

        connection_json = {}
        for c in self.conns["object"]:
            connection_json.update(c._to_exerpy(pamb, Tamb))

        return {
            "components": component_json,
            "connections": connection_json,
            "ambient_conditions": {
                "Tamb": Tamb,
                "Tamb_unit": "K",
                "pamb": pamb,
                "pamb_unit": "Pa"
            }
        }

    def save(self, json_file_path):
        r"""
        Dump the results to a json style output.

        Parameters
        ----------
        json_file_path : str
            Filename to dump results into.

        Note
        ----
        Results will be saved to specified file path
        """
        dump = {}

        # save relevant state information only
        dump["Connection"] = self._save_connections()
        dump["Component"] = self._save_components()
        dump["Bus"] = self._save_busses()

        dump = hlp._nested_dict_of_dataframes_to_dict(dump)

        with open(json_file_path, "w") as f:
            json.dump(dump, f)

    def save_csv(self, folder_path):
        """Export the results in multiple csv files in a folder structure

        - Connection.csv
        - Component/
          - Compressor.csv
          - ....
        - Bus/
          - power input bus.csv
          - ...

        Parameters
        ----------
        folder_path : str
            Path to dump results to
        """
        dump = {}
        # save relevant state information only
        dump["Connection"] = self._save_connections()
        dump["Component"] = self._save_components()
        dump["Bus"] = self._save_busses()
        hlp._nested_dict_of_dataframes_to_filetree(dump, folder_path)

    def _save_connections(self):
        """Save the connection properties.

        Returns
        -------
        pandas.DataFrame
            pandas.Dataframe of the connection results
        """
        dump = {}
        for c in self.conns["conn_type"].unique():
            dump[c] = self.results[c].replace(np.nan, None)
        return dump

    def _save_components(self):
        r"""
        Save the component properties.

        Returns
        -------
        dump : dict
            Dump of the component information.
        """
        dump = {}
        for c in self.comps['comp_type'].unique():
            dump[c] = self.results[c].replace(np.nan, None)
        return dump

    def _save_busses(self):
        r"""
        Save the bus properties.

        Returns
        -------
        dump : dict
            Dump of the component information.
        """
        dump = {}
        for label in self.busses:
            dump[label] = self.results[label]["design value"].replace(np.nan, None)
        return dump

    def _export_network(self):
        r"""Export network information

        Returns
        -------
        dict
            Serialization of network object.
        """
        return self._serialize()

    def _export_connections(self):
        """Export connection information

        Returns
        -------
        dict
            Serialization of connection objects.
        """
        connections = {}
        for c in self.conns["object"]:
            conn_type = c.__class__.__name__
            if conn_type not in connections:
                connections[conn_type] = {}
            connections[conn_type].update(c._serialize())
        return connections

    def _export_components(self):
        """Export component information

        Returns
        -------
        dict
            Dict of dicts with per class serialization of component objects.
        """
        components = {}
        for c in self.comps["comp_type"].unique():
            components[c] = {}
            for cp in self.comps.loc[self.comps["comp_type"] == c, "object"]:
                components[c].update(cp._serialize())

        return components

    def _export_busses(self):
        """Export bus information

        Returns
        -------
        dict
            Serialization of bus objects.
        """
        busses = {}
        for bus in self.busses.values():
            busses.update(bus._serialize())

        return busses


def v07_to_v08_save(path):
    """Transform the v0.7 network save to a dictionary compatible to v0.8.

    Parameters
    ----------
    path : str
        Path to the save structure

    Returns
    -------
    dict
        Dictionary of the v0.8 network save
    """
    data = {
        "Component": {},
        "Bus": {},
        "Connection": {},
    }
    component_files_path = os.path.join(path, "components")
    for file in os.listdir(component_files_path):
        df = pd.read_csv(
            os.path.join(component_files_path, file),
            sep=';', decimal='.', index_col=0
        )
        data["Component"][file.removesuffix(".csv")] = df.to_dict(orient="index")
    df = pd.read_csv(
        os.path.join(path, f"connections.csv"),
        sep=';', decimal='.', index_col=0
    )
    data["Connection"] = df.to_dict(orient="index")

    with open(os.path.join(path, "busses.json"), "r") as f:
        data["Bus"] = json.load(f)

    return data


def v07_to_v08_export(path):
    """Transform the v0.7 network export to a dictionary compatible to v0.8.

    Parameters
    ----------
    path : str
        Path to the export structure

    Returns
    -------
    dict
        Dictionary of the v0.8 network export
    """
    data = {
        "Component": {},
        "Bus": {},
        "Connection": {},
        "Network": {}
    }
    component_files_path = os.path.join(path, "components")
    for file in os.listdir(component_files_path):
        with open(os.path.join(component_files_path, file), "r") as f:
            data["Component"][file.removesuffix(".json")] = json.load(f)
    files = {"busses": "Bus", "connections": "Connection", "network": "Network"}
    for file, key in files.items():
        with open(os.path.join(path, f"{file}.json"), "r") as f:
            data[key] = json.load(f)

    return data


def _construct_components(target_class, data, nw):
    r"""
    Create TESPy component from class name and set parameters.

    Parameters
    ----------
    component : str
        Name of the component class to be constructed.

    data : dict
        Dictionary with component information.

    Returns
    -------
    dict
        Dictionary of all components of the specified type.
    """
    instances = {}
    for cp, cp_data in data.items():
        instances[cp] = target_class(cp)
        for param, param_data in cp_data.items():
            container = instances[cp].get_attr(param)
            if isinstance(container, dc):
                if "char_func" in param_data:
                    if isinstance(container, dc_cc):
                        param_data["char_func"] = CharLine(**param_data["char_func"])
                    elif isinstance(container, dc_cm):
                        param_data["char_func"] = CharMap(**param_data["char_func"])

                if "val" in param_data:
                    if "unit" in param_data and param_data["unit"] is not None:
                        param_data["val"] = nw.units.ureg.Quantity(
                            param_data["val"], param_data["unit"]
                        )
                    if "val0" in param_data:
                        param_data["val0"] = param_data["val"]
                container.set_attr(**param_data)
            else:
                instances[cp].set_attr(**{param: param_data})

    return instances


def _construct_connections(target_class, data, comps):
    r"""
    Create TESPy connection from data in the .json-file and its parameters.

    Parameters
    ----------
    data : dict
        Dictionary with connection data.

    comps : dict
        Dictionary of constructed components.

    Returns
    -------
    dict
        Dictionary of TESPy connection objects.
    """
    conns = {}

    module_name = "tespy.tools.fluid_properties.wrappers"
    _ = importlib.import_module(module_name)

    for label, conn_data in data.items():
        conns[label] = target_class(
            comps[conn_data["source"]], conn_data["source_id"],
            comps[conn_data["target"]], conn_data["target_id"],
            label=label
        )

    for label, conn_data in data.items():
        conns[label]._deserialize(conn_data, conns)

    return conns


def _construct_busses(data, comps):
    r"""
    Create busses of the network.

    Parameters
    ----------
    data : dict
        Bus information from .json file.

    comps : dict
        TESPy components dictionary.

    Returns
    -------
    dict
        Dict with TESPy bus objects.
    """
    busses = {}

    for label, bus_data in data.items():
        busses[label] = Bus(label)
        busses[label].P.set_attr(**bus_data["P"])

        components = [_ for _ in bus_data if _ != "P"]
        for cp in components:
            char = CharLine(**bus_data[cp]["char"])
            component_data = {
                "comp": comps[cp], "param": bus_data[cp]["param"],
                "base": bus_data[cp]["base"], "char": char
            }
            busses[label].add_comps(component_data)

    return busses
