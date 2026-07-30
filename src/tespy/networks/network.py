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
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from tabulate import tabulate

from tespy.components.component import component_registry

from tespy.connections.connection import ConnectionBase
from tespy.connections.connection import connection_registry
from tespy.solver import Problem
from tespy.tools import helpers as hlp
from tespy.tools import logger
from tespy.tools.characteristics import CharLine
from tespy.tools.characteristics import CharMap
from tespy.tools.data_containers import ComponentArrayProperties as dc_cap
from tespy.tools.data_containers import ComponentCharacteristicMaps as dc_cm
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import DataContainer as dc
from tespy.tools.data_containers import FluidProperties as dc_prop
from tespy.tools.global_vars import ERR
from tespy.tools.units import SI_UNITS
from tespy.tools.units import Units

# Only require cupy if Cuda shall be used
try:
    import cupy as cu
except ModuleNotFoundError:
    cu = None


class Network:
    r"""
    The container for TESPy simulations.

    The network collects components, connections and user defined equations,
    builds the system of equations describing the plant topology and
    parametrization and solves it.

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
    ...     "pressure": "bar", "pressure_difference": "bar",
    ...     "temperature": "degC"
    ... })
    >>> mynetwork.p_range = [1, 10]
    >>> type(mynetwork)
    <class 'tespy.networks.network.Network'>
    >>> mynetwork.iterinfo = False
    >>> mynetwork.iterinfo
    False
    >>> mynetwork.iterinfo = True
    >>> mynetwork.iterinfo
    True

    A simple network consisting of a source, a pipe and a sink. This example
    shows how the printout parameter can be used. We specify
    :code:`printout=False` for both connections, the pipe as well as the power
    connection. Therefore the :code:`.print_results()` method should not print
    any results.

    >>> from tespy.networks import Network
    >>> from tespy.components import Source, Sink, Pipe, HeatSink
    >>> from tespy.connections import Connection, HeatConnection
    >>> nw = Network()
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "pressure_difference": "bar",
    ...     "temperature": "degC"
    ... })
    >>> so = Source('source')
    >>> si = Sink('sink')
    >>> p = Pipe('pipe', Q=0, pr=0.95, printout=False)
    >>> h = HeatSink('heat to ambient')
    >>> a = Connection(so, 'out1', p, 'in1')
    >>> b = Connection(p, 'out1', si, 'in1')
    >>> nw.add_conns(a, b)
    >>> a.set_attr(fluid={'CH4': 1}, T=30, p=10, m=10, printout=False)
    >>> b.set_attr(printout=False)
    >>> e = HeatConnection(p, 'heat', h, 'heat', printout=False)
    >>> nw.add_conns(e)
    >>> nw.iterinfo = False
    >>> nw.solve('design')
    >>> nw.print_results()
    """

    def __init__(self, iterinfo=True, units=None, m_range=None, p_range=None, h_range=None, **kwargs):
        self._set_defaults()
        self.iterinfo = iterinfo

        if units is not None:
            self.units = units

        self.set_attr(**kwargs)

        # because the units can still be specified via the deprecated API of
        # set_attr, ranges need to be updated after set_attr!
        if m_range is not None:
            self.m_range = m_range
        if p_range is not None:
            self.p_range = p_range
        if h_range is not None:
            self.h_range = h_range

    def _serialize(self):
        # the ranges are interpreted in the default units at load time, so
        # they have to be exported in exactly those units independent of
        # which defaults were active when the range was assigned
        return {
            "m_range": list(
                self.m_range.m_as(self.units.default["mass_flow"])
            ),
            "p_range": list(
                self.p_range.m_as(self.units.default["pressure"])
            ),
            "h_range": list(
                self.h_range.m_as(self.units.default["enthalpy"])
            ),
            "units": self.units._serialize()
        }

    def _set_defaults(self):
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
        self.user_defined_var = {}
        self.subsystems = {}
        # results and specification dictionary
        self.results = {}

        # in case of a design calculation after an offdesign calculation
        self.redesign = False

        self.checked = False
        self.design_path = None
        self.iterinfo = True
        self._problem = None
        self._presolve_pending = False
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

        self.m_range = self.m_range_SI
        self.p_range = self.p_range_SI
        self.h_range = self.h_range_SI

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
        if kwargs:
            msg = (
                "The set_attr method of Network is deprecated and will be "
                "removed in version 0.12. Please explicitly call "
                "the respective set methods for specification of value "
                "ranges, units or iterinfo."
            )
            warnings.warn(msg, FutureWarning, stacklevel=2)
        self.units = kwargs.get('units', self.units)
        for key in kwargs:
            if "_unit" in key:
                msg = (
                    f"Passing '{key}' to Network.set_attr is no longer "
                    "supported. Use Network.units.set_defaults() instead."
                )
                raise TypeError(msg)

        for prop in ['m', 'p', 'h']:
            key = f"{prop}_range"
            if key in kwargs:
                msg = (
                    "Setting variable ranges through the Network.set_attr "
                    f"is deprecated and will stop working in version 0.12. "
                    f"Please use Network.{key} = [min, max] instead."
                )
                warnings.warn(msg, FutureWarning)
                logger.warning(msg)
                if key == "m_range":
                    self.m_range = kwargs[key]
                elif key == "p_range":
                    self.p_range = kwargs[key]
                else:
                    self.h_range = kwargs[key]

        self.iterinfo = kwargs.get('iterinfo', self.iterinfo)
        if "iterinfo" in kwargs:
            msg = (
                "Setting iterinfo through the Network.set_attr is deprecated "
                "and will stop working in version 0.12. Please directly "
                "specify Network.iterinfo=True/False instead."
            )
            warnings.warn(msg, FutureWarning)
            logger.warning(msg)

    def _set_iterinfo(self, value):
        if not isinstance(value, bool):
            msg = 'Network parameter iterinfo must be True or False!'
            logger.error(msg)
            raise TypeError(msg)
        else:
            self._iterinfo = value

    def _get_iterinfo(self):
        return self._iterinfo

    def _set_units(self, value):
        if not isinstance(value, Units):
            msg = (
                "The units must be an instance of class "
                "tespy.tools.units.Units."
            )
            logger.error(msg)
            raise TypeError(msg)
        else:
            self._units = value

    def _get_units(self):
        return self._units

    def _set_m_range(self, value):
        self._check_range_dtype(value, "mass flow")
        quantity = "mass_flow"
        unit = self.units.default[quantity]
        self._m_range = self.units.ureg.Quantity(np.array(value), unit)
        self.m_range_SI = self.m_range.m_as(SI_UNITS[quantity])

    def _get_m_range(self):
        return self._m_range

    def _set_p_range(self, value):
        self._check_range_dtype(value, "pressure")
        quantity = "pressure"
        unit = self.units.default[quantity]
        self._p_range = self.units.ureg.Quantity(np.array(value), unit)
        self.p_range_SI = self.p_range.m_as(SI_UNITS[quantity])

    def _get_p_range(self):
        return self._p_range

    def _set_h_range(self, value):
        self._check_range_dtype(value, "enthalpy")
        quantity = "enthalpy"
        unit = self.units.default[quantity]
        self._h_range = self.units.ureg.Quantity(np.array(value), unit)
        self.h_range_SI = self.h_range.m_as(SI_UNITS[quantity])

    def _get_h_range(self):
        return self._h_range

    @staticmethod
    def _check_range_dtype(value, property):
        if isinstance(value, list) or isinstance(value, np.ndarray):
            return
        else:
            msg = (
                f"Specify the range for {property} as list: [min, max]."
            )
            logger.error(msg)
            raise TypeError(msg)

    iterinfo = property(_get_iterinfo, _set_iterinfo)
    units = property(_get_units, _set_units)
    m_range = property(_get_m_range, _set_m_range)
    p_range = property(_get_p_range, _set_p_range)
    h_range = property(_get_h_range, _set_h_range)

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
        msg = (
            "The Network.get_attr method is deprecated and will be removed "
            "in version 0.12."
        )
        warnings.warn(msg, FutureWarning)
        logger.warning(msg)

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
                del self.subsystems[subsystem.label]

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
            warnings.warn(
                f"Connection with label {label} not found. Returning None is "
                "deprecated and will raise a KeyError in version 0.12.",
                FutureWarning,
                stacklevel=2,
            )
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
            warnings.warn(
                f"Component with label {label} not found. Returning None is "
                "deprecated and will raise a KeyError in version 0.12.",
                FutureWarning,
                stacklevel=2,
            )
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
                comp_type = comp.__class__.__name__
                if comp_type in self.results:
                    self.results[comp_type].drop(
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

    def get_ude(self, label):
        r"""
        Get UserDefinedEquation via label.

        Parameters
        ----------
        label : str
            Label of the UserDefinedEquation object.

        Returns
        -------
        c : tespy.tools.helpers.UserDefinedEquation
            UserDefinedEquation object with specified label, None if no
            UserDefinedEquation of the network has this label.
        """
        try:
            return self.user_defined_eq[label]
        except KeyError:
            warnings.warn(
                f"UserDefinedEquation with label {label} not found. Returning "
                "None is deprecated and will raise a KeyError in version 0.12.",
                FutureWarning,
                stacklevel=2,
            )
            return None

    def add_udv(self, *args):
        r"""
        Add a user defined variable to the network.

        Parameters
        ----------
        c : tespy.tools.helpers.UserDefinedVariable
            The objects to be added to the network, UserDefinedVariable
            objects ci :code:`add_udv(c1, c2, c3, ...)`.
        """
        for c in args:
            if not isinstance(c, hlp.UserDefinedVariable):
                msg = (
                    'Must provide tespy.tools.helpers.UserDefinedVariable '
                    'objects as parameters.'
                )
                logger.error(msg)
                raise TypeError(msg)

            elif c.label in self.user_defined_var:
                msg = (
                    'There is already a UserDefinedVariable with the label '
                    f'{c.label} . The UserDefinedVariable labels must be '
                    'unique within a network'
                )
                logger.error(msg)
                raise ValueError(msg)

            self.user_defined_var[c.label] = c
            msg = f"Added UserDefinedVariable {c.label} to network."
            logger.debug(msg)

    def del_udv(self, *args):
        """
        Remove a user defined variable from the network.

        Parameters
        ----------
        c : tespy.tools.helpers.UserDefinedVariable
            The objects to be deleted from the network,
            UserDefinedVariable objects ci :code:`del_udv(c1, c2, c3, ...)`.
        """
        for c in args:
            del self.user_defined_var[c.label]
            msg = f"Deleted UserDefinedVariable {c.label} from network."
            logger.debug(msg)

    def get_udv(self, label):
        r"""
        Get UserDefinedVariable via label.

        Parameters
        ----------
        label : str
            Label of the UserDefinedVariable object.

        Returns
        -------
        c : tespy.tools.helpers.UserDefinedVariable
            UserDefinedVariable object with specified label.
        """
        return self.user_defined_var[label]

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

    @property
    def problem(self):
        """Solver :code:`Problem` instance of the most recent solve call.

        Holds the variable space, equation lookups and the state of the
        newton algorithm (residual vector, jacobian, increment).
        """
        if self._problem is None:
            msg = (
                "The network does not hold a problem instance yet. It is "
                "created with the first call of the solve method."
            )
            raise hlp.TESPyNetworkError(msg)
        return self._problem

    _PROBLEM_ATTRIBUTES = frozenset({
        "variables_dict", "variable_counter", "residual", "jacobian",
        "increment", "increment_filter", "lin_dep", "residual_history",
        "singularity_msg", "iter", "max_iter", "min_iter", "use_cuda",
        "robust_relax", "oscillation_damping",
        "num_comp_eq", "num_conn_eq", "num_ude_eq",
        "get_sorted_residual_index",
    })

    def __getattr__(self, name):
        if (
                name in Network._PROBLEM_ATTRIBUTES
                and self.__dict__.get("_problem") is not None
            ):
            msg = (
                f"Accessing Network.{name} is deprecated and will stop "
                f"working in version 0.12, use Network.problem.{name} "
                "instead."
            )
            if name == "residual_history":
                msg += (
                    " Note that its values hold the scaled maximum "
                    "residual norm of every iteration since version "
                    "0.11, previous versions stored the euclidean norm "
                    "of the unscaled residual vector."
                )
            warnings.warn(msg, FutureWarning, stacklevel=2)
            return getattr(self._problem, name)
        raise AttributeError(
            f"'{type(self).__name__}' object has no attribute '{name}'"
        )

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
            source_mask = self.conns["source"] == comp
            target_mask = self.conns["target"] == comp

            comp.inl, comp.outl = self._resolve_comp_conn_domain(
                source_mask, target_mask, comp.inlets(), comp.outlets()
            )
            comp.power_inl, comp.power_outl = self._resolve_comp_conn_domain(
                source_mask, target_mask,
                comp.powerinlets(), comp.poweroutlets(), "PowerConnection"
            )
            comp.num_power_i = len(comp.powerinlets())
            comp.num_power_o = len(comp.poweroutlets())

            comp.heat_inl, comp.heat_outl = self._resolve_comp_conn_domain(
                source_mask, target_mask,
                comp.heatinlets(), comp.heatoutlets(), "HeatConnection"
            )
            comp.num_heat_i = len(comp.heatinlets())
            comp.num_heat_o = len(comp.heatoutlets())

            # set up results and specification dataframes
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

    def _resolve_comp_conn_domain(
        self, source_mask, target_mask, inlet_ids, outlet_ids, conn_type=None
    ):
        """Return :code:`(inl, outl)` connection lists for one domain.

        Parameters
        ----------
        source_mask, target_mask : boolean Series
            Rows in :code:`self.conns` where the component is source / target.
        inlet_ids, outlet_ids : list[str]
            Port IDs returned by the component's :code:`*inlets()` / :code:`*outlets()`.
        conn_type : str, optional
            If given, further restrict to rows whose :code:`conn_type` column
            matches this class name (e.g. :code:`"PowerConnection"`).
        """
        if conn_type is not None:
            type_mask = self.conns["conn_type"] == conn_type
            src = self.conns[source_mask & self.conns["source_id"].isin(outlet_ids) & type_mask]
            tgt = self.conns[target_mask & self.conns["target_id"].isin(inlet_ids) & type_mask]
        else:
            src = self.conns[source_mask & self.conns["source_id"].isin(outlet_ids)]
            tgt = self.conns[target_mask & self.conns["target_id"].isin(inlet_ids)]

        return (
            self.conns.loc[tgt["target_id"].sort_values().index, "object"].tolist(),
            self.conns.loc[src["source_id"].sort_values().index, "object"].tolist(),
        )

    def _check_components(self):
        for comp in self.comps['object']:
            comp._validate_connections()

    def _prepare_problem(self):
        r"""
        Initialise the network depending on calculation mode.

        Design

        - Generic fluid composition and fluid property initialisation.
        - Starting values from initialisation path if provided.

        Offdesign

        - Check offdesign path specification.
        - Set component and connection design point properties.
        - Switch from design/offdesign parameter specification.
        """
        self._problem = Problem(self)

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
        for branch_data in self.fluid_wrapper_branches.values():
            first = branch_data["connections"][0]
            if self.conns.loc[first.label, "object"] != first:
                self._create_fluid_wrapper_branches()
                break

        self._propagate_fluid_wrappers()
        self._init_connection_result_datastructure()

        self._prepare_solve_mode()
        # this method will distribute units and set SI values from given values
        # and units
        self._transform_user_input_to_SI()

        self._problem.build()

        # generic fluid property initialisation
        self._set_starting_values()

        msg = 'Network initialised.'
        logger.info(msg)

    def _propagate_fluid_wrappers(self):

        connections_in_wrapper_branches = []
        self.all_fluids = []
        for branch_name, branch_data in self.fluid_wrapper_branches.items():
            all_connections = [c for c in branch_data["connections"]]

            any_fluids_set = []
            engines = {}
            back_ends = {}
            wrapper_kwargs = {}
            any_fluids = []

            all_components = [c for c in branch_data["components"]]
            for cp in all_components:
                any_fluids += cp._add_missing_fluids(all_connections)

            any_fluids0 = []
            mixing_rules = []
            for c in all_connections:
                any_fluids_set += list(c.fluid.is_set)
                any_fluids += list(c.fluid.val.keys())
                any_fluids0 += list(c.fluid.val0.keys())
                if c.mixing_rule is not None:
                    mixing_rules += [c.mixing_rule]

            for c in all_connections:
                for f in set(any_fluids):
                    if f in c.fluid.engine:
                        if f in engines and engines[f] != c.fluid.engine[f]:
                            raise ValueError("")
                        engines[f] = c.fluid.engine[f]

                    if f in c.fluid.back_end:
                        if f in back_ends and back_ends[f] != c.fluid.back_end[f]:
                            raise ValueError("")
                        back_ends[f] = c.fluid.back_end[f]

                    if f in c.fluid.wrapper_kwargs:
                        if f in wrapper_kwargs and wrapper_kwargs[f] != c.fluid.wrapper_kwargs[f]:
                            raise ValueError("")
                        wrapper_kwargs[f] = c.fluid.wrapper_kwargs[f]

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

            potential_fluids = set(any_fluids_set + any_fluids + any_fluids0)
            self.all_fluids += list(potential_fluids)
            num_potential_fluids = len(potential_fluids)
            if num_potential_fluids == 0:
                msg = (
                    "The following connections of your network are missing any "
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
                for f, w_kwargs in wrapper_kwargs.items():
                    c.fluid.wrapper_kwargs[f] = w_kwargs

                c._create_fluid_wrapper()

            msg = (
                f"Fluid domain {branch_name} ({len(all_connections)} "
                "connections): potential fluids "
                f"[{', '.join(sorted(potential_fluids))}], mixing rule "
                f"{mixing_rule}."
            )
            logger.debug(msg)

            connections_in_wrapper_branches += all_connections

        missing_wrappers = (
            {c for c in self.conns["object"] if c._has_fluid_vector}
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
                self._load_offdesign_state()

            self._prepare_offdesign()

        else:
            # reset any preceding offdesign calculation
            self._prepare_design()

    def _create_fluid_wrapper_branches(self):

        self.fluid_wrapper_branches = {}
        mask = self.comps["object"].apply(lambda x: x._is_wrapper_branch_source)
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

    def _transform_user_input_to_SI(self):
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
                        ).m_as(SI_UNITS[param.quantity])
                    else:
                        param.set_SI_from_val(self.units)

        for cp in self.comps["object"]:
            for param, value in cp.parameters.items():
                if isinstance(value, dc_prop) and value.is_set:
                    value.set_SI_from_val(self.units)

    def _prepare_design(self):
        r"""
        Initialise a design calculation.

        Offdesign parameters are unset, design parameters are set. If
        :code:`local_offdesign` is :code:`True` for connections or components,
        the design point information are read from the .csv-files in the
        respective :code:`design_path`. In this case, the design values are
        unset, the offdesign values set.
        """
        # connections
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

                entries = self._load_network_state(path)[c.__class__.__name__]
                # write data to connections
                self._write_design_state_to_connection(c, entries)

            else:
                c._reset_design(self.redesign)
                # unset all design values

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

                local_design = self._load_network_state(path)
                data = local_design[c]
                # resolve design label (may differ from cp.label)
                label = self._find_isolated_comp_label(cp, data)
                # write data
                self._write_design_state_to_component(cp, data, label)

                # store adjacent connection design values from the component's
                # own design_path for use in offdesign equations
                cp._local_connection_design_state = {}
                for adj_conn in cp.all_connections:
                    conn_type = adj_conn.__class__.__name__
                    if conn_type in local_design:
                        conn_entries = local_design[conn_type]
                        matched_row = self._find_conn_in_isolated_design(
                            adj_conn, cp, label, conn_entries
                        )
                        if matched_row is not None:
                            cp._local_connection_design_state[adj_conn.label] = (
                                adj_conn._get_design_state_SI(matched_row, self.units)
                            )
                        else:
                            msg = (
                                "Could not retrieve connection design point "
                                "data in local_offdesign of component "
                                f"{cp.label} for the connections adjacent to "
                                "the component."
                            )
                            raise KeyError(msg)

                # unset design parameters
                for var in cp.design:
                    cp.get_attr(var).is_set = False

                # set offdesign parameters
                set_values = []
                activated = []
                for var in cp.offdesign:
                    # set variables provided in .offdesign attribute
                    data = cp.get_attr(var)
                    data.is_set = True

                    # take nominal values from design point
                    if isinstance(data, dc_cp):
                        cp.get_attr(var).val = cp.get_attr(var).design
                        set_values.append(var)
                    else:
                        activated.append(var)

                if cp.design or cp.offdesign:
                    parts = [f"unset [{', '.join(cp.design)}]"]
                    if set_values:
                        parts.append(
                            f"set [{', '.join(set_values)}] to design point "
                            "values"
                        )
                    if activated:
                        parts.append(f"activated [{', '.join(activated)}]")
                    msg = (
                        f"Switched component {cp.label} to local offdesign: "
                        f"{', '.join(parts)}."
                    )
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

    def _prepare_offdesign(self):
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

                if c.design or c.offdesign:
                    msg = (
                        f"Switched connection {c.label} to offdesign: unset "
                        f"[{', '.join(c.design)}], set "
                        f"[{', '.join(c.offdesign)}] to design point values."
                    )
                    logger.debug(msg)

                c.new_design = False

        for cp in self.comps['object']:
            if not cp.local_design:
                # unset variables provided in .design attribute
                for var in cp.design:
                    cp.get_attr(var).is_set = False

                set_values = []
                activated = []
                for var in cp.offdesign:
                    # set variables provided in .offdesign attribute
                    data = cp.get_attr(var)
                    data.is_set = True

                    # take nominal values from design point
                    if isinstance(data, dc_cp):
                        data.val_SI = data.design
                        data.set_val_from_SI(self.units)
                        set_values.append(var)
                    else:
                        activated.append(var)

                if cp.design or cp.offdesign:
                    parts = [f"unset [{', '.join(cp.design)}]"]
                    if set_values:
                        parts.append(
                            f"set [{', '.join(set_values)}] to design point "
                            "values"
                        )
                    if activated:
                        parts.append(f"activated [{', '.join(activated)}]")
                    msg = (
                        f"Switched component {cp.label} to offdesign: "
                        f"{', '.join(parts)}."
                    )
                    logger.debug(msg)

                cp.new_design = False

    def _load_offdesign_state(self):
        r"""
        Read design point information from specified :code:`design_path`.

        If a :code:`design_path` has been specified individually for components
        or connections, the data will be read from the specified individual
        path instead.
        """
        # components with offdesign parameters
        components_with_parameters = [
            cp.label for cp in self.comps["object"] if len(cp.parameters) > 0
        ]
        # fetch all components, reindex with label
        df_comps = self.comps.loc[components_with_parameters].copy()
        # iter through unique types of components (class names)
        state = self._load_network_state(self.design_path)
        # iter through all components of this type and set data
        for _, row in df_comps.iterrows():
            entries = state[row["comp_type"]]
            comp = row["object"]
            path = comp.design_path
            # in offdesign mode any individually specified design_path is used
            # to load this component's design reference, regardless of
            # local_offdesign
            if path is not None:
                _individual_design = self._load_network_state(path)
                data = _individual_design[row["comp_type"]]
                label = self._find_isolated_comp_label(comp, data)
                self._write_design_state_to_component(comp, data, label)
                # write adjacent connections design state from individual
                # design_path to the component
                comp._local_connection_design_state = {}
                for adj_conn in comp.all_connections:
                    conn_type = adj_conn.__class__.__name__
                    if conn_type in _individual_design:
                        conn_entries = _individual_design[conn_type]
                        matched_row = self._find_conn_in_isolated_design(
                            adj_conn, comp, label, conn_entries
                        )
                        if matched_row is not None:
                            comp._local_connection_design_state[adj_conn.label] = (
                                adj_conn._get_design_state_SI(matched_row, self.units)
                            )
                        else:
                            msg = (
                                "Could not retrieve connection design point "
                                f"data for component {comp.label}, connection "
                                f"{adj_conn.label}."
                            )
                            raise KeyError(msg)
            else:
                # write data to components
                self._write_design_state_to_component(comp, entries, comp.label)

        # iter through connections
        for c in self.conns['object']:
            conn_type = c.__class__.__name__
            entries = state[conn_type]
            # read data of connections with individual design_path
            path = c.design_path
            if path is not None:
                entries = self._load_network_state(path)[conn_type]

            self._write_design_state_to_connection(c, entries)

    def _find_isolated_comp_label(self, comp, comp_entries):
        """
        Resolve which label in *comp_entries* corresponds to *comp* for
        isolated design loading.

        - Exact match -> return :code:`comp.label`
        - Single-type fallback: label not found but exactly one entry ->
          return that entry's label (the isolated design contains exactly one
          component of that type, so it is unambiguous)
        - Ambiguous (multiple entries, no exact match) -> raise error
        """
        if comp.label in comp_entries:
            return comp.label
        elif len(comp_entries) == 1:
            return next(iter(comp_entries))
        msg = (
            f"Could not unambiguously resolve the label for component "
            f"'{comp.label}' in the isolated design file: multiple entries "
            f"exist ({', '.join(comp_entries)}) and none match exactly."
        )
        raise hlp.TESPyNetworkError(msg)

    def _find_conn_in_isolated_design(self, adj_conn, comp, comp_label, conn_entries):
        """
        Find the entry in *conn_entries* that corresponds to *adj_conn* when
        loading an isolated design file.

        Matching strategy (in order):

        1. Direct label match (:code:`adj_conn.label` in :code:`conn_entries`).
        2. Port-based topology match using the :code:`source` / :code:`target` /
           :code:`source_id` / :code:`target_id` fields stored by
           :py:meth:`tespy.connections.connection.Connection.collect_results`.

        Parameters
        ----------
        adj_conn : tespy.connections.connection.BaseConnection
            BaseConnection type object
        comp : tespy.components.component.Component
            Component type object
        comp_label : str
            Label of the component to look for inside the connection entries.
        conn_entries : dict
            Mapping of connection labels to their data dicts.

        Returns
        -------
        dict or None
            Data dict for the matched connection, or None if not found.
        """
        # --- direct label match ---
        if adj_conn.label in conn_entries:
            return conn_entries[adj_conn.label]

        # --- port-based topology match ---
        if comp_label is None or not conn_entries:
            return None
        any_row = next(iter(conn_entries.values()))
        if 'source' not in any_row or 'target' not in any_row:
            return None

        if adj_conn in comp.all_inlets:
            matches = [
                row for row in conn_entries.values()
                if row.get('target') == comp_label
                and row.get('target_id') == adj_conn.target_id
            ]
        else:
            matches = [
                row for row in conn_entries.values()
                if row.get('source') == comp_label
                and row.get('source_id') == adj_conn.source_id
            ]

        if len(matches) == 1:
            return matches[0]
        return None

    def _write_design_state_to_component(self, c, entries, label):
        r"""
        Write design point information to components.

        Parameters
        ----------
        c : tespy.components.component.Component
            Write design point information to this component.

        entries : dict
            Mapping of component labels to their design point data dicts.

        label : str
            Label of the component inside the data. It can differ under the
            condition of an individual design_path specified for that
            component.
        """
        if label not in entries:
            # no matches in the connections of the network and the design files
            msg = (
                f"Could not find component '{label}' in design case file. "
                "This is is critical only to components, which need to load "
                "design values from this case."
            )
            logger.debug(msg)
            return
        # write component design data
        c._set_design_parameters(self.mode, entries[label])

    def _write_design_state_to_connection(self, c, entries):
        r"""
        Write design point information to connections.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            Write design point information to this connection.

        entries : dict
            Mapping of connection labels to their design point data dicts.
        """
        if c.label not in entries:
            # no matches in the connections of the network and the design files
            msg = (
                f"Could not find connection '{c.label}' in design case. "
                "Please make sure no connections have been modified or "
                "components have been relabeled for your offdesign "
                "calculation."
            )
            logger.error(msg)
            raise hlp.TESPyNetworkError(msg)

        c._set_design_params(entries[c.label], self.units)

    def _write_starting_values_to_connection(self, c, entries):
        r"""
        Write parameter information from init_path to connections.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            Write init path information to this connection.

        entries : dict
            Mapping of connection labels to their state data dicts.
        """
        if c.label not in entries:
            # no matches in the connections of the network and the design files
            msg = f"Could not find connection {c.label} in init path file."
            logger.debug(msg)
            return False

        c._set_starting_values(entries[c.label], self.units)
        c.good_starting_values = True
        return True

    def _set_starting_values(self):
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
            state = self._load_network_state(self.init_path)
        # improved starting values for referenced connections,
        # specified vapour content values, temperature values as well as
        # subccooling/overheating and state specification
        num_init_path = 0
        num_previous = 0
        num_generic = 0
        for c in self.conns['object']:
            if self.init_path is not None and self._write_starting_values_to_connection(
                    c, state[c.__class__.__name__]
                ):
                num_init_path += 1
            elif c.good_starting_values:
                num_previous += 1
            else:
                num_generic += 1

        # known values - specifications, presolved values and user provided
        # starting values - propagate through the approximate affine
        # relations of the components. Pressures and mass flows propagate
        # first, so the component anchors and the temperature and quality
        # based enthalpy precalculation compute from propagated pressures.
        # Their results then act as additional sources of the enthalpy
        # propagation, which fills the gaps between them.
        seeded = []
        for c in self.conns["object"]:
            seeded += c._seed_starting_values(self.units)

        covered = self._problem.propagate_starting_values(seeded, ("m", "p"))
        num_propagated = len(covered)
        covered |= set(seeded)

        # the reconciled temperatures carry known levels across heat
        # exchangers into coupled circuits: two phase positions receive
        # their saturation pressure, declared single phase positions their
        # enthalpy at the reconciled temperature
        seeded_set = set(seeded)
        temperature_field = self._problem.starting_temperature_field()
        h_sources = []
        for c in self.conns["object"]:
            h_sources += c._apply_temperature_field(
                temperature_field, covered, seeded_set
            )

        for c in self.conns["object"]:
            h_sources += c._guess_starting_values(self.units, covered)

        assigned = self._problem.propagate_starting_values(
            seeded + h_sources, ("h",)
        )
        num_propagated += len(assigned)
        covered |= assigned

        seeded = set(seeded)
        for c in self.conns["object"]:
            c._finalize_starting_values(self.units, covered, seeded, self)

        # with the property field complete, the flow determining equations
        # are linear in the mass and energy flow variables and solve their
        # levels and branch ratios exactly with respect to the field
        num_propagated += self._problem.presolve_flow_variables(seeded)

        # only after the flow presolve: an exactly zero enthalpy difference
        # self-guards its energy balance rows there, while a small one would
        # imply arbitrarily large mass flows
        for cp in self.comps["object"]:
            num_propagated += cp._separate_flat_enthalpy_starts(seeded)

        # here reference values can be updated, e.g. a reference temperature
        # if the starting value of the reference connection is not yet updated
        # then the calculation of the reference can cause issues, therefore:
        # first update all of the starting values and only then to
        # precalculation of reference values
        for c in self.conns["object"]:
            c._precalc_guess_values_for_references()

        for cp in self.comps["object"]:
            for key, variable in cp.get_variables().items():
                # for components every variable should be an actual variable
                # if variable.is_var:
                if np.isnan(variable.val):
                    variable.val = (variable.min_val + variable.max_val) / 2
                variable.set_SI_from_val(self.units)
                variable.set_reference_val_SI(variable._val_SI)

        for udv in self.user_defined_var.values():
            for key, variable in udv.get_variables().items():
                # user defined variables are SI only, the local container
                # always holds a value
                variable.set_reference_val_SI(variable._val_SI)

        msg = (
            f"Starting values: {num_init_path} connections from init_path, "
            f"{num_previous} from a previous solution, {num_generic} generic, "
            f"{num_propagated} variables assigned by propagation."
        )
        logger.debug(msg)


    @staticmethod
    def _load_network_state(json_path: str | bytes | bytearray | Path | dict):
        r"""
        Read network state from given file or in-memory dict.

        Parameters
        ----------
        json_path : str | bytes | bytearray | Path | dict
            Path to a saved network state file, a JSON string, or a state
            dict as returned by :meth:`Network.save` with no arguments.
        """
        if isinstance(json_path, dict):
            data = json_path
        else:
            data = None
            if not isinstance(json_path, Path):
                try:
                    data = json.loads(json_path)
                except json.JSONDecodeError as e:
                    msg = (
                        "The provided json_path could not be decoded. If this is not "
                        "a valid json string, please provide a valid file path instead of "
                        "%s"
                    )
                    logger.debug(msg, str(json_path))
            if data is None:
                with open(json_path, "r") as f:
                    data = json.load(f)

        def _row(d):
            return {col: np.nan if val is None else val for col, val in d.items()}

        state = {}
        # TODO: Let this somehow run through connection-registry and not hardcoded names
        if any(k in data["Connection"] for k in ("Connection", "PowerConnection", "HeatConnection")):
            for key, value in data["Connection"].items():
                state[key] = {str(k): _row(v) for k, v in value.items()}
        # TODO: deprecate
        # this is for compatibility of older savestates
        else:
            state["Connection"] = {str(k): _row(v) for k, v in data["Connection"].items()}

        for key, value in data["Component"].items():
            state[key] = {str(k): _row(v) for k, v in value.items()}

        return state

    def get_structure(self) -> dict:
        """Get a serializable description of the mathematical structure.

        The result joins with the class level schema of
        :py:mod:`tespy.tools.schema` through class names, parameter names and
        quantities and with the network serialization through object labels
        and port identifiers. The problem has to be prepared, e.g. by solving
        with :code:`init_only=True`.

        Returns
        -------
        dict
            Dictionary with the keys :code:`variables`, :code:`equations`,
            :code:`connections` and :code:`components`. Variables carry their
            state (:code:`specified`, :code:`presolved` or :code:`variable`),
            their affine relation to the reference variable and the solver
            column. Equations carry their mathematical kind (:code:`affine`,
            :code:`linear` or :code:`nonlinear`), their origin
            (:code:`topology` or :code:`specification`), their state
            (:code:`consumed` or :code:`active`) and the structural variables
            they relate.
        """
        warnings.warn(
            "The structure API is not yet stable and may change without "
            "notice in future releases.",
            FutureWarning,
            stacklevel=2,
        )
        if self._problem is None or self._problem.structure_graph is None:
            msg = (
                "The mathematical structure is only available after "
                "preprocessing, e.g. by calling "
                "nw.solve(mode, init_only=True) first."
            )
            raise hlp.TESPyNetworkError(msg)

        variables = self._problem.structure_variables()
        equations = self._problem.structure_equations()

        connections = [
            {
                "label": c.label,
                "source": c.source.label,
                "source_id": c.source_id,
                "target": c.target.label,
                "target_id": c.target_id,
                "class_name": type(c).__name__,
            }
            for c in self.conns["object"]
        ]
        components = [
            {"label": c.label, "class_name": type(c).__name__}
            for c in self.comps["object"]
        ]

        return {
            "variables": variables,
            "equations": equations,
            "connections": connections,
            "components": components,
        }

    def get_linear_dependent_variables(self) -> list:
        """Get a list with sublists containing linear dependent variables

        Returns
        -------
        list
            List of lists of linear dependent variables
        """
        return self.problem.get_linear_dependent_variables()

    def get_presolved_equations(self) -> list:
        """Get the list of equations, that has been presolved with their
        respective parent object

        Returns
        -------
        list
            list of presolved equations
        """
        return self.problem.get_presolved_equations()

    def print_presolved_equations(self):
        """Print a formatted table of presolved equations."""
        rows = self.get_presolved_equations()
        print(f"Presolved equations ({len(rows)} total):")
        if rows:
            print(tabulate(rows, headers=["Object", "Equation"], tablefmt="simple"))

    def get_variables_before_presolve(self) -> list:
        """Get the list of variables before presolving.

        Returns
        -------
        list
            list of original variables
        """
        return self.problem.get_variables_before_presolve()

    def print_variables_before_presolve(self):
        """Print a formatted table of all variables before presolving."""
        rows = self.get_variables_before_presolve()
        print(f"Variables before presolving ({len(rows)} total):")
        if rows:
            print(tabulate(rows, headers=["Object", "Property"], tablefmt="simple"))

    def get_presolved_variables(self) -> list:
        """Get the list of presolved variables with their respective parent
        object and property.

        Returns
        -------
        list
            list of presolved variables
        """
        return self.problem.get_presolved_variables()

    def print_presolved_variables(self):
        """Print a formatted table of presolved variables."""
        rows = self.get_presolved_variables()
        print(f"Presolved variables ({len(rows)} total):")
        if rows:
            print(tabulate(rows, headers=["Object", "Property"], tablefmt="simple"))

    def get_variables(self) -> dict:
        """Get all variables of the presolved problem with their respective
        represented original variables.

        Returns
        -------
        dict
            variable number and property with the list of represented variables
        """
        return self.problem.get_variables()

    def print_variables(self):
        """Print a formatted table of variables after presolving."""
        variables = self.get_variables()
        print(f"Variables after presolving ({len(variables)} total):")
        rows = [
            (
                var_idx,
                var_type,
                ", ".join(f"{lbl} ({prop})" for lbl, prop in represents),
            )
            for (var_idx, var_type), represents in variables.items()
        ]
        if rows:
            print(tabulate(rows, headers=["#", "Property", "Represents"], tablefmt="simple"))

    def get_variable_values(self, block=None) -> dict:
        """Get the current values of all variables of the presolved problem.

        Every variable is listed with all of the original variables it
        represents and their individual SI values, which can differ through
        the affine relation between linearly dependent variables.

        Parameters
        ----------
        block : int
            Restrict the result to the variables of the block with the
            given id of the block decomposition, default: :code:`None`
            (all variables).

        Returns
        -------
        dict
            Variable number and property with the list of represented
            variables as tuples of object label, property and SI value.
        """
        return self.problem.get_variable_values(block=block)

    def print_variable_values(self, block=None):
        """Print a formatted table of the current variable values,
        optionally restricted to the variables of a single block."""
        variables = self.get_variable_values(block=block)
        if block is not None:
            print(
                f"Current variable values of block {block} "
                f"({len(variables)} total):"
            )
        else:
            print(f"Current variable values ({len(variables)} total):")
        rows = [
            (var_idx, var_type, f"{lbl} ({prop})", value)
            for (var_idx, var_type), represents in variables.items()
            for lbl, prop, value in represents
        ]
        if rows:
            print(tabulate(
                rows,
                headers=["#", "Property", "Represents", "SI value"],
                tablefmt="simple",
                floatfmt=".6e",
            ))

    def set_variable_value(self, label, prop, value):
        """Set the SI value of a variable of the presolved problem.

        The value can be set through any of the linearly dependent variables
        a solver variable represents, it propagates to the underlying
        reference container through the affine relation. For example setting
        the enthalpy of a connection whose enthalpy is linked to another
        connection updates both.

        Parameters
        ----------
        label : str
            Label of the connection or component holding the variable.

        prop : str
            Name of the variable (e.g. :code:`'m'` or :code:`'h'`).

        value : float
            SI value to impose.

        Raises
        ------
        KeyError
            In case no object with the given label exists or the object does
            not have a variable of the given name.
        tespy.tools.helpers.TESPyNetworkError
            In case the property is not part of the variable space, e.g.
            because it is specified or has been presolved.
        """
        if label in self.conns.index:
            obj = self.conns.loc[label, "object"]
        elif label in self.comps.index:
            obj = self.comps.loc[label, "object"]
        elif label in self.user_defined_var:
            obj = self.user_defined_var[label]
        else:
            msg = (
                f"There is no connection, component or user defined variable "
                f"with label {label}."
            )
            raise KeyError(msg)

        self.problem.set_variable_value(obj, prop, value)

    def get_equations(self) -> dict:
        """Get the actual equations after presolving the problem

        Returns
        -------
        dict
            Lookup with equation number as index and tuple of label and
            parameter defining the equation. In case one parameter defines
            multiple equations, the same equation is repeated.
        """
        return self.problem.get_equations()

    def print_equations(self):
        """Print a formatted table of equations after presolving."""
        equations = self.get_equations()
        print(f"Equations after presolving ({len(equations)} total):")
        rows = [
            (eq_num, label, self._problem.format_eq_name(eq_name))
            for eq_num, (label, eq_name) in sorted(equations.items())
        ]
        if rows:
            print(tabulate(rows, headers=["Eq#", "Object", "Equation"], tablefmt="simple"))

    def get_equations_with_dependents(self) -> dict:
        """Get the equations together with the variables they depend on.

        Returns
        -------
        dict
            Lookup with equation (component, (parameter_label, number)) and
            the variables it depends on as a list
            (variable number, variable type)
        """
        return self.problem.get_equations_with_dependents()

    def print_equations_with_dependents(self):
        """Print a formatted table of equations and the variables they depend on."""
        incidence = self.problem.get_incidence()
        equations = self.problem.get_equations()
        variable_order = [key for key, _ in self.problem.get_variables()]
        print(f"Equations with dependent variables ({len(incidence)} total):")
        rows = []
        for eq_idx, dependents in sorted(incidence.items()):
            label, eq_name = equations[eq_idx]
            dep_set = set(dependents)
            dep_str = ", ".join(
                self.problem.format_var_label(v_idx)
                for v_idx in variable_order if v_idx in dep_set
            )
            rows.append((eq_idx, label, self.problem.format_eq_name(eq_name), dep_str))
        if rows:
            print(tabulate(
                rows,
                headers=["Eq#", "Object", "Equation", "Dependent variables"],
                tablefmt="simple",
            ))

    def print_incidence_matrix(self, block_order=False):
        """Print the incidence matrix with equation rows and variable columns.

        Parameters
        ----------
        block_order : boolean
            Sort the equations and variables in the solve order of the block
            decomposition, so the block lower triangular form of the problem
            becomes visible, default: :code:`False`.
        """
        incidence = self.problem.get_incidence()
        equations = self.problem.get_equations()
        if block_order:
            decomposition = self.problem.decompose()
            eq_indices = [
                eq for block in decomposition.blocks for eq in block.equations
            ]
            all_var_indices = [
                col for block in decomposition.blocks
                for col in block.variables
            ]
            block_of_eq = {
                eq: block.id
                for block in decomposition.blocks for eq in block.equations
            }
            block_of_var = {
                col: block.id
                for block in decomposition.blocks for col in block.variables
            }
        else:
            eq_indices = sorted(incidence.keys())
            all_var_indices = sorted({
                v_idx
                for deps in incidence.values()
                for v_idx in deps
            })

        col_labels = [
            self.problem.format_var_label(v_idx) for v_idx in all_var_indices
        ]

        rows = []
        previous_block = None
        for eq_idx in eq_indices:
            label, eq_name = equations[eq_idx]
            row_label = f"{label}.{self.problem.format_eq_name(eq_name)}"
            dep_set = set(incidence[eq_idx])
            if block_order:
                # entries within the equation's own block are marked with X,
                # dependencies on variables of other blocks with x
                block_id = block_of_eq[eq_idx]
                entries = [
                    (
                        "X" if block_of_var.get(v) == block_id else "x"
                    ) if v in dep_set else "-"
                    for v in all_var_indices
                ]
                marker = block_id if block_id != previous_block else ""
                previous_block = block_id
                rows.append([marker, row_label] + entries)
            else:
                rows.append(
                    [row_label]
                    + ["x" if v in dep_set else "-" for v in all_var_indices]
                )

        print("Incidence matrix:")
        headers = ["Block", ""] if block_order else [""]
        print(tabulate(rows, headers=headers + col_labels, tablefmt="simple"))

    def get_blocks(self) -> list:
        """Get the block lower triangular decomposition of the problem.

        Returns
        -------
        list
            One dictionary per block in solve order with the block id, its
            kind, the ids of the blocks it depends on, the equations as
            tuples of object label and equation name, the short variable
            labels and the solve status (:code:`None` before solving).
        """
        decomposition = self.problem.decompose()
        return [
            {
                "id": block.id,
                "kind": block.kind,
                "prerequisites": decomposition.precedence.get(block.id, []),
                "equations": list(block.equation_labels),
                "variables": list(block.variable_labels),
                "status": block.status,
            }
            for block in decomposition.blocks
        ]

    def print_blocks(self):
        """Print a formatted table of the block decomposition in solve order."""
        blocks = self.get_blocks()
        print(f"Block lower triangular decomposition ({len(blocks)} blocks):")
        rows = [
            (
                block["id"],
                block["kind"],
                ", ".join(str(dep) for dep in block["prerequisites"]),
                ", ".join(
                    f"{lbl}.{self._problem.format_eq_name(name)}"
                    for lbl, name in block["equations"]
                ),
                ", ".join(block["variables"]),
            )
            for block in blocks
        ]
        if rows:
            print(tabulate(
                rows,
                headers=["#", "Kind", "Needs", "Equations", "Variables"],
                tablefmt="simple",
            ))

    def get_structural_analysis(self) -> list:
        """Get the over- and under-determined parts of the problem.

        Returns
        -------
        list
            One dictionary per defective part of the maximum matching of
            the incidence with its kind (:code:`"overdetermined"` or
            :code:`"underdetermined"`), the equations as tuples of object
            label and equation name and the short variable labels. An empty
            list means the problem is structurally sound.
        """
        return self.problem.get_structural_analysis()

    def print_structural_analysis(self):
        """Print the over- and under-determined parts of the problem.

        In an over-determined part more equations compete for the involved
        variables than the variables can satisfy - one of the underlying
        specifications must be removed. In an under-determined part the
        involved equations cannot determine all of the variables - a
        specification is missing. Prints that the problem is structurally
        sound in case neither exists.
        """
        analysis = self.get_structural_analysis()
        if not analysis:
            print("The problem is structurally sound.")
            return
        for part in analysis:
            eq_str = ", ".join(
                f"{lbl}.{self._problem.format_eq_name(name)}"
                for lbl, name in part["equations"]
            )
            var_str = ", ".join(part["variables"])
            if part["kind"] == "overdetermined":
                print(
                    "Over-determined part - the following equations compete "
                    "for the involved variables:"
                )
                print(f"  Equations: {eq_str}")
                print(f"  Involved variables: {var_str}")
            else:
                print(
                    "Under-determined part - the following variables cannot "
                    "be determined by the involved equations:"
                )
                print(f"  Variables: {var_str}")
                print(f"  Involved equations: {eq_str}")

    def get_block_jacobian(self, block) -> dict:
        """Get the linear system of a failed block at its state of failure.

        The jacobian and the residual vector restricted to the block's
        equations and variables, recorded as evaluated in the failing
        iteration of the block's last solution attempt. Only available for
        blocks whose solution failed.

        Parameters
        ----------
        block : int
            Id of the block of the block decomposition.

        Returns
        -------
        dict
            Formatted equation and variable labels together with the
            jacobian matrix and the residual vector of the block.
        """
        return self.problem.get_block_jacobian(block)

    def print_block_jacobian(self, block):
        """Print the jacobian and residual of a failed block at its state
        of failure.

        Entries the incidence matrix does not declare are left blank, a
        derivative that evaluated to zero although its dependency is
        declared is marked as :code:`missing`.
        """
        data = self.get_block_jacobian(block)
        print(f"Jacobian of block {block} at the state of failure:")
        rows = [(
            "(variable values)",
            *(f"{value:.6e}" for value in data["values"]),
            ""
        )]
        for i, equation in enumerate(data["equations"]):
            cells = []
            for j in range(len(data["variables"])):
                value = data["jacobian"][i, j]
                if value == 0.0 and data["expected"][i, j]:
                    cells.append("missing")
                elif value == 0.0:
                    cells.append("")
                else:
                    cells.append(f"{value:.6e}")
            rows.append((equation, *cells, f"{data['residual'][i]:.6e}"))
        print(tabulate(
            rows,
            headers=["Equation"] + data["variables"] + ["Residual"],
            tablefmt="simple",
        ))

    def get_block_states(self, block, at="current") -> list:
        """Get the states of the connections a block touches.

        The properties of every connection involved in the block's
        equations and variables. What is reported comes from the connection
        classes, e.g. mass flow, pressure, enthalpy, temperature and phase
        for fluid connections, the energy flow for power connections.

        Parameters
        ----------
        block : int
            Id of the block of the block decomposition.

        at : str
            :code:`"current"` (default) evaluates the states from the
            current variable values - while the solution process is paused
            at the block this is its entry state, including any
            modification made in between. :code:`"failure"` returns the
            states as recorded in the failing iteration of the block's
            last solution attempt, only available for blocks whose
            solution failed.

        Returns
        -------
        list
            One dict per connection with its label, the reported properties
            and the list of properties that are variables of the block.
        """
        return self.problem.get_block_states(block, at=at)

    def print_block_states(self, block, at="current"):
        """Print the states of the connections a block touches, either at
        the current variable values or as recorded in the failing
        iteration (:code:`at="failure"`)."""
        states = self.get_block_states(block, at=at)
        props = []
        for state in states:
            for prop in state:
                if prop in ("label", "block_variables") or prop in props:
                    continue
                props.append(prop)
        rows = []
        for state in states:
            cells = [state["label"]]
            for prop in props:
                value = state.get(prop)
                if value is None or (
                        isinstance(value, float) and np.isnan(value)):
                    cells.append("")
                elif isinstance(value, str):
                    cells.append(value)
                else:
                    mark = " *" if prop in state["block_variables"] else ""
                    cells.append(f"{value:.6e}{mark}")
            rows.append(cells)
        if at == "failure":
            reference = "at the state of failure"
        else:
            reference = "at the current variable values"
        print(
            f"States of the connections of block {block} {reference} "
            "(*: variable of the block):"
        )
        print(tabulate(
            rows, headers=["Connection"] + props, tablefmt="simple"
        ))

    def print_residuals(self):
        """Print a formatted table of equation residuals, sorted by magnitude."""
        if self._problem is None or self._problem.residual is None:
            print("Residuals are not available before the first solve call.")
            return
        equations = self.problem.get_equations()
        scaled = self.problem.scaled_residual()
        rows = []
        for eq_idx in self.problem.get_sorted_residual_index():
            label, eq_name = equations[eq_idx]
            rows.append((
                eq_idx, label, self.problem.format_eq_name(eq_name),
                scaled[eq_idx], self.problem.residual[eq_idx]
            ))
        print(
            f"Residuals per equation ({len(rows)} total, sorted by scaled "
            "magnitude):"
        )
        if rows:
            print(tabulate(
                rows,
                headers=["Eq#", "Object", "Equation", "Scaled", "Residual"],
                tablefmt="simple",
                floatfmt=".3e",
            ))

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
        return self.problem.get_linear_dependents_by_object(obj, prop)

    def solve(self, mode, init_path=None, design_path=None,
              max_iter=50, min_iter=2, init_only=False, init_previous=True,
              use_cuda=False, print_results=True, robust_relax=False, skip_postprocess=False,
              oscillation_damping=False, block_solve=True,
              pause_on_block_failure=False):
        r"""
        Solve the network.

        - Check network consistency.
        - Initialise calculation and preprocessing.
        - Perform actual calculation.
        - Postprocessing.

        It is possible to check programmatically, if a network was solved
        successfully with the `.converged` attribute.

        Parameters
        ----------
        mode : str
            Choose from 'design' and 'offdesign'.

        init_path : str | Path | dict
            Path to a previously saved network state (e.g.
            :code:`nw.save('myplant/test.json')`), or the dict returned by
            :code:`nw.save(as_dict=True)`.

        design_path : str | Path | dict
            Path to the saved design-case state (e.g.
            :code:`nw.save('myplant/test.json')`), or the dict returned by
            :code:`nw.save(as_dict=True)`.

        max_iter : int
            Maximum number of iterations before calculation stops, default: 50.

        min_iter : int
            Minimum number of iterations of the simultaneous solution before
            convergence can be accepted, default: 2. Convergence is only
            accepted in an iteration in which the value bounds and
            convergence heuristics did not modify any variable, so the
            parameter is a hard floor on top of that, not the primary guard.
            Block-wise solving accepts every block individually on the same
            criterion.

        init_only : boolean
            Perform initialisation only, default: :code:`False`.

        init_previous : boolean
            Initialise the calculation with values from the previous
            calculation, default: :code:`True`.

        use_cuda : boolean
            Use cuda instead of numpy for matrix inversion, default:
            :code:`False`.

        robust_relax : boolean
            Apply a ramped relaxation factor that starts near zero and grows to
            1 over the first quarter of :code:`max_iter` iterations. Helps
            avoid divergence from poor starting values, at the cost of slower
            early convergence. Default: :code:`False`.

        oscillation_damping : boolean
            Detect Newton oscillations caused by non-smooth residuals (e.g.
            phase-transition kinks in sectioned heat exchangers) and dampen
            them automatically. When a residual component changes sign between
            two consecutive iterations - indicating an overshoot - the
            increments for all variables that equation depends on are halved
            before being applied. This converts the oscillating Newton step
            into a bisection-like contraction and restores monotone convergence
            without requiring an external bracketing loop. Default:
            :code:`False`.

        block_solve : boolean
            Decompose the equation system into its block lower triangular
            form and solve the blocks in precedence order instead of solving
            the full system simultaneously. Scalar blocks are solved with a
            bracketing fallback on oscillation. Experimental, default:
            :code:`False`.

        pause_on_block_failure : boolean
            Pause the block-wise solution process at the first block that
            does not converge instead of escalating to the coupled solution
            stages, default: :code:`False`. The variable space stays loaded
            (:code:`status` is 20), the variables of the failed block can
            be inspected with
            :py:meth:`~tespy.networks.network.Network.print_variable_values`
            and modified with
            :py:meth:`~tespy.networks.network.Network.set_variable_value`,
            and
            :py:meth:`~tespy.networks.network.Network.solve_continue`
            retries the block and continues the solution process.

        Note
        ----
        For more information on the solution process have a look at the online
        documentation at tespy.readthedocs.io in the section "TESPy modules".
        """
        self.presolve(
            mode, init_path=init_path, design_path=design_path,
            init_previous=init_previous, check=not init_only
        )

        if init_only:
            return

        self.solve_continue(
            max_iter=max_iter, min_iter=min_iter, use_cuda=use_cuda,
            print_results=print_results, robust_relax=robust_relax,
            skip_postprocess=skip_postprocess,
            oscillation_damping=oscillation_damping, block_solve=block_solve,
            pause_on_block_failure=pause_on_block_failure
        )

    def presolve(self, mode, init_path=None, design_path=None,
                 init_previous=True, check=True):
        r"""
        Prepare and presolve the problem without running the solver.

        The network is checked, the problem is built and presolved and the
        starting values are assigned. The variable space stays loaded
        afterwards, so it can be inspected through :code:`Network.problem`
        and the print methods, and variable values can be modified before
        continuing the solution process with
        :py:meth:`~tespy.networks.network.Network.solve_continue`.

        Parameters
        ----------
        mode : str
            Choose from 'design' and 'offdesign'.

        init_path : str | Path | dict
            Path to a previously saved network state (e.g.
            :code:`nw.save('myplant/test.json')`), or the dict returned by
            :code:`nw.save(as_dict=True)`.

        design_path : str | Path | dict
            Path to the saved design-case state (e.g.
            :code:`nw.save('myplant/test.json')`), or the dict returned by
            :code:`nw.save(as_dict=True)`.

        init_previous : boolean
            Initialise the calculation with values from the previous
            calculation, default: :code:`True`.

        check : boolean
            Check whether the number of parameters matches the number of
            variables, default: :code:`True`. The check is skipped when
            preparing with :code:`solve(mode, init_only=True)`, so ill
            determined problems can be inspected.
        """
        self.status = 99
        if self._problem is not None and self._problem.paused:
            self._problem._release_pause()
            self._problem.unload_variables()
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
        self.init_previous = init_previous

        if mode not in ['offdesign', 'design']:
            msg = 'Mode must be "design" or "offdesign".'
            logger.error(msg)
            raise ValueError(msg)
        else:
            self.mode = mode

        if not self.checked:
            self.check_topology()

        msg = (
            "Network information:\n"
            f" - Number of components: {len(self.comps)}\n"
            f" - Number of connections: {len(self.conns)}\n"
        )
        logger.debug(msg)

        self._prepare_problem()
        self._presolve_pending = True

        if not check:
            return

        try:
            self._problem.check_determination()
        except hlp.TESPyNetworkError:
            self.status = self._problem.status
            raise

    def solve_continue(self, max_iter=50, min_iter=2, use_cuda=False,
                       print_results=True, robust_relax=False,
                       skip_postprocess=False, oscillation_damping=False,
                       block_solve=True, pause_on_block_failure=None):
        r"""
        Run the solver on the previously prepared problem.

        Continues after :py:meth:`~tespy.networks.network.Network.presolve`
        or :code:`solve(mode, init_only=True)`, taking into account any
        modification of variable values made in between. The solver
        parameters correspond to the ones of
        :py:meth:`~tespy.networks.network.Network.solve`.

        When the solution process is paused at a failed block, the block is
        retried with the current variable values and the solution process
        continues from there. :code:`pause_on_block_failure` defaults to
        keeping the value of the initiating solve call, passing
        :code:`False` explicitly hands a still failing block over to the
        standard escalation stages instead.
        """
        if self._problem is None or not (
                self._presolve_pending or self._problem.paused):
            msg = (
                "There is no prepared problem to continue from. Call "
                "Network.presolve or Network.solve first."
            )
            raise hlp.TESPyNetworkError(msg)

        self.skip_postprocess = skip_postprocess
        if self.skip_postprocess:
            msg = (
                "Postprocessing will be skipped, violations of "
                "physical/operational are not reported or logged!"
            )
            logger.debug(msg)

        if use_cuda and cu is None:
            msg = (
                'Specifying use_cuda=True requires cupy to be installed on '
                'your machine. Numpy will be used instead.'
            )
            logger.warning(msg)
            use_cuda = False

        msg = 'Starting solver.'
        logger.info(msg)

        self._presolve_pending = False
        try:
            self._problem.solve_loop(
                max_iter=max_iter, min_iter=min_iter, use_cuda=use_cuda,
                robust_relax=robust_relax,
                oscillation_damping=oscillation_damping,
                iterinfo=self.iterinfo, print_results=print_results,
                block_solve=block_solve,
                pause_on_block_failure=pause_on_block_failure
            )
        except ValueError as e:
            self.status = 99
            msg = f"Simulation crashed due to an unexpected error:\n{e}"
            logger.exception(msg)
            self._problem._release_pause()
            self._problem.unload_variables()
            return

        self.status = self._problem.status
        if self._problem.paused:
            # the variable space stays loaded for inspection and
            # modification of the failed block's variables
            return

        self._problem.unload_variables()

        if self.status == 3:
            logger.error(self._problem.singularity_msg)
            return

        elif self.status == 2:
            logger.warning(self._problem.no_progress_message())
            return

        self._postprocess()

        msg = 'Calculation complete.'
        logger.info(msg)
        return

    def _postprocess(self):
        r"""Calculate connection and component parameters."""
        _converged = self._postprocess_connections()
        _converged = self._postprocess_components() and _converged

        if self.status == 0 and not _converged:
            self.status = 1

        msg = 'Postprocessing complete.'
        logger.info(msg)

    def _postprocess_connections(self):
        """Process the Connection results."""
        _converged = True
        buckets = {}
        for c in self.conns['object']:
            c.good_starting_values = True
            _converged = c.calc_results(self.units, self.skip_postprocess) and _converged
            if self.skip_postprocess:
                continue
            conn_type = c.__class__.__name__
            if conn_type not in buckets:
                buckets[conn_type] = ([], [])
            buckets[conn_type][0].append(c.label)
            buckets[conn_type][1].append(c.collect_results(self.all_fluids))
        for conn_type, (labels, rows) in buckets.items():
            cols = self.results[conn_type].columns
            self.results[conn_type] = pd.DataFrame(rows, index=labels, columns=cols)
        return _converged

    def _postprocess_components(self):
        """Process the component results."""
        # components
        _converged = True
        if self.skip_postprocess:
            return _converged

        for cp in self.comps['object']:
            cp.calc_parameters()
            _converged = _converged and cp.check_parameter_bounds()
            # this thing could be somewhere else
            for key, value in cp.parameters.items():
                if isinstance(value, dc_cap):
                    value.set_val_from_SI(self.units)
                elif isinstance(value, dc_prop):
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

        buckets = {}
        for cp in self.comps['object']:
            result = cp.collect_results()
            if len(result) == 0:
                continue
            key = cp.__class__.__name__
            if key not in buckets:
                buckets[key] = ([], [])
            buckets[key][0].append(cp.label)
            buckets[key][1].append(result)
        for key, (labels, rows) in buckets.items():
            cols = self.results[key].columns
            self.results[key] = pd.DataFrame(rows, index=labels, columns=cols)

        return _converged

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
                cols = [
                    col for col in c.iloc[0]._get_result_attributes()
                    if not col.endswith("_unit")
                ]
                if len(cols) > 0:
                    df = df[cols].dropna(axis=1, how="all")
                    for col in df.columns:
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
            if val_SI < param_obj.min_val - ERR or val_SI > param_obj.max_val + ERR:
                return f"{coloring['err']}{val}{coloring['end']}"
            if param_obj.is_var:
                return f"{coloring['var']}{val}{coloring['end']}"
            if param_obj.is_set:
                return f"{coloring['set']}{val}{coloring['end']}"
            return str(val)
        else:
            return np.nan


    @classmethod
    def from_dict(cls, network_data):
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

        msg = 'Created network.'
        logger.info(msg)

        nw.check_topology()

        return nw

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
        - Network.json

        Example
        -------
        Create a network and export it. This is followed by loading the network
        from the exported json file. All network information stored will be
        passed to a new network object. Components and connections will be
        accessible by label. The following example setup is simple gas
        turbine setup with compressor, combustion chamber and turbine. The fuel
        is fed from a pipeline and throttled to the required pressure while
        keeping the temperature at a constant value.

        >>> from tespy.components import (
        ...     Sink, Source, CombustionChamber, TurboCompressor, Turbine,
        ...     SimpleHeatExchanger, PowerBus, PowerSink, Generator
        ... )
        >>> from tespy.connections import Connection, Ref, PowerConnection
        >>> from tespy.networks import Network
        >>> import os
        >>> nw = Network()
        >>> nw.iterinfo = False
        >>> nw.units.set_defaults(**{
        ...     "pressure": "bar", "pressure_difference": "bar",
        ...     "temperature": "degC", "enthalpy": "kJ/kg",
        ...     "power": "MW"
        ... })
        >>> air = Source('air')
        >>> f = Source('fuel')
        >>> compressor = TurboCompressor('compressor')
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
        efficiencies are not considered in this example. See
        :py:class:`tespy.components.power.motor.Motor` and
        :py:class:`tespy.components.power.generator.Generator` for modelling
        conversion efficiencies between mechanical and electrical power.

        >>> combustion.set_attr(lamb=None)
        >>> c3.set_attr(T=1100)
        >>> c1.set_attr(m=None)
        >>> e4.set_attr(E=1)
        >>> nw.solve('design')
        >>> nw.assert_convergence()
        >>> design_state = nw.save(as_dict=True)
        >>> _ = nw.export('exported_nwk.json')
        >>> mass_flow = round(nw.get_conn('c01').m.val_SI, 1)
        >>> compressor.set_attr(igva='var')
        >>> nw.solve('offdesign', design_path=design_state)
        >>> round(turbine.eta_s.val, 1)
        0.9
        >>> e4.set_attr(E=0.75)
        >>> nw.solve('offdesign', design_path=design_state)
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
        >>> imported_nwk.iterinfo = False
        >>> imported_nwk.solve('design')
        >>> imported_nwk.problem.lin_dep
        False
        >>> round(imported_nwk.get_conn('c01').m.val_SI, 1) == mass_flow
        True
        >>> round(imported_nwk.get_comp('turbine').eta_s.val, 3)
        0.9
        >>> imported_nwk.get_comp('compressor').set_attr(igva='var')
        >>> imported_nwk.solve('offdesign', design_path=design_state)
        >>> round(imported_nwk.get_comp('turbine').eta_s.val, 3)
        0.9
        >>> imported_nwk.get_conn('e4').set_attr(E=0.75)
        >>> imported_nwk.solve('offdesign', design_path=design_state)
        >>> round(imported_nwk.get_comp('turbine').eta_s.val, 3) == eta_s_t
        True
        >>> round(imported_nwk.get_comp('compressor').igva.val, 3) == igva
        True
        >>> os.remove('exported_nwk.json')
        """
        msg = f'Reading network data from base path {json_file_path}.'
        logger.info(msg)

        with open(json_file_path, "r") as f:
            network_data = json.load(f)

        return cls.from_dict(network_data)

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
        # a plain numeric specification means "in the network's default
        # units". The unit is only attached to the containers when the
        # network is transformed to SI for solving: exporting before that
        # would serialize such values with the global default units,
        # misinterpreting them on import
        containers = [
            c.get_attr(key)
            for c in self.conns["object"] for key in c.property_data
        ] + [
            cp.get_attr(key)
            for cp in self.comps["object"] for key in cp.parameters
        ]
        for container in containers:
            if not isinstance(container, dc_prop) or not container.is_set:
                continue
            if container._val_is_quantity or container.quantity is None:
                continue
            try:
                container._assign_default_unit_to_val(self.units)
            except KeyError:
                continue

        export = {}
        export["Network"] = self._export_network()
        export["Connection"] = self._export_connections()
        export["Component"] = self._export_components()

        if json_file_path:
            os.makedirs(os.path.dirname(os.path.abspath(json_file_path)), exist_ok=True)
            with open(json_file_path, "w") as f:
                json.dump(export, f, indent=2)

            logger.debug(f'Model information saved to {json_file_path}.')

        return export

    def save(self, json_file_path: str | Path | None = None, as_dict: bool = False) -> None | dict | str:
        r"""
        Dump the results to a json style output.

        Parameters
        ----------
        json_file_path : str | Path | None
            Filename to dump results into. If :code:`None`, the state is returned
            in-memory (as dict when :code:`as_dict=True`, otherwise as JSON string).
        as_dict : bool
            If :code:`True` and :code:`json_file_path` is :code:`None`, return the state as
            a dict that can be passed directly as :code:`design_path` or
            :code:`init_path` in a subsequent :meth:`solve` call. Default
            :code:`False`; the :code:`False` behaviour (returning a JSON string) is
            deprecated and will be removed in version 0.12.

        Returns
        -------
        None
            If a file path is provided, results are saved to file.
        dict
            If :code:`json_file_path` is :code:`None` and :code:`as_dict=True`.
        str
            If :code:`json_file_path` is :code:`None` and :code:`as_dict=False`
            (deprecated).
        """
        dump = {}

        # save relevant state information only
        dump["Connection"] = self._save_connections()
        dump["Component"] = self._save_components()

        dump = hlp._nested_dict_of_dataframes_to_dict(dump)

        if json_file_path is None:
            if as_dict:
                return dump
            msg = (
                "Calling Network.save() without a file path returns a JSON "
                "string, which is deprecated and will be removed in a future "
                "release. Use Network.save(as_dict=True) to get a dict that "
                "can be passed directly as design_path or init_path."
            )
            warnings.warn(msg, FutureWarning)
            return json.dumps(dump, indent=2)

        os.makedirs(os.path.dirname(os.path.abspath(json_file_path)), exist_ok=True)
        with open(json_file_path, "w") as f:
            json.dump(dump, f)

    def save_csv(self, folder_path):
        """Export the results in multiple csv files in a folder structure

        - Connection.csv
        - Component/
          - Compressor.csv
          - ....

        Parameters
        ----------
        folder_path : str
            Path to dump results to
        """
        dump = {}
        # save relevant state information only
        dump["Connection"] = self._save_connections()
        dump["Component"] = self._save_components()
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
            try:
                container = instances[cp].get_attr(param)
            except KeyError:
                msg = (
                    f"The parameter {param} of component {cp} is not "
                    "available anymore, ignoring the data from the file. "
                    "The file may originate from an older version of tespy."
                )
                logger.warning(msg)
                continue
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
