# -*- coding: utf-8

"""Module for tespy network class.

The network is the container for every TESPy simulation. The network class
automatically creates the system of equations describing topology and
parametrisation of a specific model and solves it.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/networks/networks.py

SPDX-License-Identifier: MIT
"""

# reading .csv
import ast
# ordered dicts for fluid composition vector
from collections import Counter, OrderedDict
# calculation of molar masses and gas constants
from CoolProp.CoolProp import PropsSI as CPPSI
# logging messages
import logging
# numpy functions
import numpy as np
from numpy.linalg import inv
from numpy.linalg import norm
# checking/creating folders
import os
# DataFrames for connections and components
import pandas as pd
# printing results
from tabulate import tabulate

from tespy import connections as con

from tespy.components.basics import (sink, source, subsystem_interface,
                                     cycle_closer)
from tespy.components.combustion import combustion_chamber, combustion_engine
from tespy.components.heat_exchangers import heat_exchanger
from tespy.components.nodes import drum, merge, splitter
from tespy.components.reactors import water_electrolyzer

from tespy.tools import data_containers as dc
from tespy.tools import fluid_properties as fp
from tespy.tools import helpers as hlp

from tespy.tools.global_vars import molar_masses, gas_constants, err
# calculation time
from time import time


class network:
    r"""
    Class component is the base class of all TESPy components.

    Parameters
    ----------
    fluids : list
        A list of all fluids within the network container.

    m_unit : str
        Specify the unit for mass flow: 'kg / s', 't / h'.

    v_unit : str
        Specify the unit for volumetric flow: 'm3 / s', 'm3 / h', 'l / s',
        'l / h'.

    p_unit : str
        Specify the unit for pressure: 'Pa', 'psi', 'bar', 'MPa'.

    h_unit : str
        Specify the unit for mass flow: 'J / kg', 'kJ / kg', 'MJ / kg'.

    T_unit : str
        Specify the unit for mass flow: 'K', 'C', 'F'.

    p_range : list
        List with minimum and maximum values for pressure value range.

    h_range : list
        List with minimum and maximum values for enthalpy value range.

    T_range : list
        List with minimum and maximum values for temperature value range.

    iterinfo : boolean
        Print convergence progress to console.

    Note
    ----
    Unit specification is optional: If not specified the SI unit (first
    element in above lists) will be applied!

    Range specification is optional, too. The value range is used to stabilise
    the newton algorith. For more information see the "getting started" section
    in the online-documentation.

    Example
    -------
    Basic example for a setting up a tespy.networks.network object. Specifying
    the fluids is mandatory! Unit systems, fluid property range and iterinfo
    are optional.

    Standard value for iterinfo is True. This will print out convergence
    progress to the console. You can suop the printouts by setting this
    property to false.

    >>> from tespy.networks import network
    >>> fluid_list = ['water', 'air', 'R134a']
    >>> mynetwork = network(fluids=fluid_list, p_unit='bar', T_unit='C')
    >>> mynetwork.set_attr(p_range=[1, 10])
    >>> type(mynetwork)
    <class 'tespy.networks.networks.network'>
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

    >>> from tespy.networks import network
    >>> from tespy.components import source, sink, pipe
    >>> from tespy.connections import connection, bus
    >>> nw = network(['CH4'], T_unit='C', p_unit='bar', v_unit='m3 / s')
    >>> so = source('source')
    >>> si = sink('sink')
    >>> p = pipe('pipe', Q=0, pr=0.95, printout=False)
    >>> a = connection(so, 'out1', p, 'in1')
    >>> b = connection(p, 'out1', si, 'in1')
    >>> nw.add_conns(a, b)
    >>> a.set_attr(fluid={'CH4': 1}, T=30, p=10, m=10, printout=False)
    >>> b.set_attr(printout=False)
    >>> b = bus('heat bus')
    >>> b.add_comps({'c': p})
    >>> nw.add_busses(b)
    >>> b.set_attr(printout=False)
    >>> nw.set_attr(iterinfo=False)
    >>> nw.solve('design')
    >>> nw.print_results()
    """

    def __init__(self, fluids, **kwargs):

        # initialisation of basic properties
        self.checked = False
        # in case of a design calculation after an offdesign calculation
        self.redesign = False
        # connection dataframe
        self.conns = pd.DataFrame(columns=['s', 's_id', 't', 't_id'])
        # list for busses
        self.busses = OrderedDict()
        # default design_path value
        self.design_path = None

        # fluid list and constants
        if isinstance(fluids, list):
            self.fluids = sorted(fluids)
        else:
            msg = ('Please provide a list containing the network\'s fluids on '
                   'creation.')
            logging.error(msg)
            raise TypeError(msg)

        msg = 'Network fluids are: '
        for f in self.fluids:
            msg += f + ', '
            if 'INCOMP::' in f:
                # molar mass and gas constant not available for incompressibles
                molar_masses[f] = 1
                gas_constants[f] = 1

            elif 'TESPy::' not in f:
                # calculating molar masses/gas constants for network's fluids
                # tespy_fluid molar mass/gas constant are added on lut creation
                molar_masses[f] = CPPSI('M', f)
                gas_constants[f] = CPPSI('GAS_CONSTANT', f)

        msg = msg[:-2] + '.'
        logging.debug(msg)

        # initialise fluid property memorisation function for this network
        fp.memorise.add_fluids(self.fluids)

        # available unit systems
        # mass flow
        self.m = {
            'kg / s': 1,
            't / h': 3.6
        }
        # pressure
        self.p = {
            'Pa': 1,
            'psi': 6.8948e3,
            'bar': 1e5,
            'MPa': 1e6
        }
        # enthalpy
        self.h = {
            'J / kg': 1,
            'kJ / kg': 1e3,
            'MJ / kg': 1e6
        }
        # temperature
        self.T = {
            'C': [273.15, 1],
            'F': [459.67, 5 / 9],
            'K': [0, 1]
        }
        # volumetric flow
        self.v = {
            'm3 / s': 1,
            'l / s': 1e-3,
            'm3 / h': 1 / 3600,
            'l / h': 1 / 3.6
        }
        # SI unit specification
        self.SI_units = {
              'm': 'kg / s',
              'p': 'Pa',
              'h': 'J / kg',
              'T': 'K',
              'v': 'm3 / s'
              }

        # iterinfo
        self.iterinfo = True

        # standard unit set
        self.m_unit = self.SI_units['m']
        self.p_unit = self.SI_units['p']
        self.h_unit = self.SI_units['h']
        self.T_unit = self.SI_units['T']
        self.v_unit = self.SI_units['v']

        msg = ('Default unit specifications: '
               'mass flow: ' + self.m_unit + ', ' +
               'pressure: ' + self.p_unit + ', ' +
               'enthalpy: ' + self.h_unit + ', ' +
               'temperature: ' + self.T_unit + ', ' +
               'volumetric flow: ' + self.v_unit + '.')
        logging.debug(msg)

        # generic value range
        self.m_range_SI = np.array([-1e12, 1e12])
        self.p_range_SI = np.array([2e2, 300e5])
        self.h_range_SI = np.array([1e3, 7e6])
        self.T_range_SI = np.array([273.16, 1773.15])

        msg = ('Default mass flow limits, '
               'min: ' + str(self.m_range_SI[0]) + ' ' + self.m_unit +
               ', max: ' + str(self.m_range_SI[1]) + ' ' + self.m_unit + ', ')
        logging.debug(msg)

        msg = ('Default pressure limits, '
               'min: ' + str(self.p_range_SI[0]) + ' ' + self.p_unit +
               ', max: ' + str(self.p_range_SI[1]) + ' ' + self.p_unit + ', ')
        logging.debug(msg)

        msg = ('Default enthalpy limits, '
               'min: ' + str(self.h_range_SI[0]) + ' ' + self.h_unit +
               ', max: ' + str(self.h_range_SI[1]) + ' ' + self.h_unit + ', ')
        logging.debug(msg)

        msg = ('Default temperature limits, '
               'min: ' + str(self.T_range_SI[0]) + ' ' + self.T_unit +
               ', max: ' + str(self.T_range_SI[1]) + ' ' + self.T_unit + ', ')
        logging.debug(msg)

        self.set_attr(**kwargs)

    def set_attr(self, **kwargs):
        r"""
        Set, resets or unsets attributes of a network.

        Parameters
        ----------
        m_unit : str
            Specify the unit for mass flow: 'kg / s', 't / h'.

        v_unit : str
            Specify the unit for volumetric flow: 'm3 / s', 'm3 / h', 'l / s',
            'l / h'.

        p_unit : str
            Specify the unit for pressure: 'Pa', 'psi', 'bar', 'MPa'.

        h_unit : str
            Specify the unit for enthalpy: 'J / kg', 'kJ / kg', 'MJ / kg'.

        T_unit : str
            Specify the unit for temperature: 'K', 'C', 'F'.

        m_range : list
            List with minimum and maximum values for mass flow value range.

        p_range : list
            List with minimum and maximum values for pressure value range.

        h_range : list
            List with minimum and maximum values for enthalpy value range.

        T_range : list
            List with minimum and maximum values for temperature value range.

        iterinfo : boolean
            Print convergence progress to console.
        """
        # unit sets
        if 'm_unit' in kwargs.keys():
            if kwargs['m_unit'] not in self.m.keys():
                msg = ('Allowed units for mass flow are: ' +
                       str(self.m.keys()))
                logging.error(msg)
                raise ValueError(msg)
            else:
                self.m_unit = kwargs['m_unit']
                msg = 'Setting mass flow unit: ' + self.m_unit
                logging.debug(msg)

        if 'p_unit' in kwargs.keys():
            if kwargs['p_unit'] not in self.p.keys():
                msg = ('Allowed units for pressure are: ' + str(self.p.keys()))
                logging.error(msg)
                raise ValueError(msg)
            else:
                self.p_unit = kwargs['p_unit']
                msg = 'Setting pressure unit: ' + self.p_unit
                logging.debug(msg)

        if 'h_unit' in kwargs.keys():
            if kwargs['h_unit'] not in self.h.keys():
                msg = ('Allowed units for pressure are: ' + str(self.h.keys()))
                logging.error(msg)
                raise ValueError(msg)
            else:
                self.h_unit = kwargs['h_unit']
                msg = 'Setting enthalpy unit: ' + self.h_unit
                logging.debug(msg)

        if 'T_unit' in kwargs.keys():
            if kwargs['T_unit'] not in self.T.keys():
                msg = ('Allowed units for pressure are: ' + str(self.T.keys()))
                logging.error(msg)
                raise ValueError(msg)
            else:
                self.T_unit = kwargs['T_unit']
                msg = 'Setting temperature unit: ' + self.T_unit
                logging.debug(msg)

        if 'v_unit' in kwargs.keys():
            if kwargs['v_unit'] not in self.v.keys():
                msg = ('Allowed units for pressure are: ' + str(self.v.keys()))
                logging.error(msg)
                raise ValueError(msg)
            else:
                self.v_unit = kwargs['v_unit']
                msg = 'Setting volumetric flow unit: ' + self.v_unit
                logging.debug(msg)

        # value ranges
        if 'm_range' in kwargs.keys():
            if not isinstance(kwargs['m_range'], list):
                msg = ('Specify the value range as list: [m_min, m_max]')
                logging.error(msg)
                raise TypeError(msg)
            else:
                self.m_range_SI = (np.array(kwargs['m_range']) *
                                   self.m[self.m_unit])

            msg = ('Setting mass flow limits, min: ' +
                   str(self.m_range_SI[0]) + ' ' + self.SI_units['m'] +
                   ', max: ' + str(self.m_range_SI[1]) + ' ' +
                   self.SI_units['m'] + '.')
            logging.debug(msg)

        if 'p_range' in kwargs.keys():
            if not isinstance(kwargs['p_range'], list):
                msg = ('Specify the value range as list: [p_min, p_max]')
                logging.error(msg)
                raise TypeError(msg)
            else:
                self.p_range_SI = (np.array(kwargs['p_range']) *
                                   self.p[self.p_unit])
            msg = ('Setting pressure limits, min: ' +
                   str(self.p_range_SI[0]) + ' ' + self.SI_units['p'] +
                   ', max: ' + str(self.p_range_SI[1]) + ' ' +
                   self.SI_units['p'] + '.')
            logging.debug(msg)

        if 'h_range' in kwargs.keys():
            if not isinstance(kwargs['h_range'], list):
                msg = ('Specify the value range as list: [h_min, h_max]')
                logging.error(msg)
                raise TypeError(msg)
            else:
                self.h_range_SI = (np.array(kwargs['h_range']) *
                                   self.h[self.h_unit])
            msg = ('Setting enthalpy limits, min: ' +
                   str(self.h_range_SI[0]) + ' ' + self.SI_units['h'] +
                   ', max: ' + str(self.h_range_SI[1]) + ' ' +
                   self.SI_units['h'] + '.')
            logging.debug(msg)

        if 'T_range' in kwargs.keys():
            if not isinstance(kwargs['T_range'], list):
                msg = ('Specify the value range as list: [T_min, T_max]')
                logging.error(msg)
                raise TypeError(msg)
            else:
                self.T_range_SI = ((np.array(kwargs['T_range']) +
                                    self.T[self.T_unit][0]) *
                                   self.T[self.T_unit][1])
            msg = ('Setting temperature limits, min: ' +
                   str(self.T_range_SI[0]) + ' ' + self.SI_units['T'] +
                   ', max: ' + str(self.T_range_SI[1]) + ' ' +
                   self.SI_units['T'] + '.')
            logging.debug(msg)

        # update non SI value ranges
        self.m_range = self.m_range_SI / self.m[self.m_unit]
        self.p_range = self.p_range_SI / self.p[self.p_unit]
        self.h_range = self.h_range_SI / self.h[self.h_unit]
        self.T_range = (self.T_range_SI / self.T[self.T_unit][1] -
                        self.T[self.T_unit][0])

        for f in self.fluids:
            if 'TESPy::' in f:
                fp.memorise.vrange[f][0] = self.p_range_SI[0]
                fp.memorise.vrange[f][1] = self.p_range_SI[1]
                fp.memorise.vrange[f][2] = self.T_range_SI[0]
                fp.memorise.vrange[f][3] = self.T_range_SI[1]

        self.iterinfo = kwargs.get('iterinfo', self.iterinfo)

        if not isinstance(self.iterinfo, bool):
            msg = ('Network parameter iterinfo must be True or False!')
            logging.error(msg)
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
            logging.error(msg)
            raise KeyError(msg)

    def add_subsys(self, *args):
        r"""
        Add one or more subsystems to the network.

        Parameters
        ----------
        c : tespy.components.subsystems.subsystem
            The subsystem to be added to the network, subsystem objects si
            :code:`network.add_subsys(s1, s2, s3, ...)`.
        """
        for subsys in args:
            for c in subsys.conns.values():
                self.add_conns(c)

    def add_nwks(self, *args):
        """
        Add connections from a different network.

        :param args: network objects si :code:`add_subsys(s1, s2, s3, ...)`
        :type args: tespy.networks.network
        :returns: no return value
        """
        for nw in args:
            for c in nw.conns.index:
                self.add_conns(c)

    def add_conns(self, *args):
        r"""
        Add one or more connections to the network.

        Parameters
        ----------
        c : tespy.connections.connection
            The connection to be added to the network, connections objects ci
            :code:`add_conns(c1, c2, c3, ...)`.
        """
        for c in args:
            if not isinstance(c, con.connection):
                msg = ('Must provide tespy.connections.connection objects as '
                       'parameters.')
                logging.error(msg)
                raise TypeError(msg)

            c.good_starting_values = False

            self.conns.loc[c] = [c.s, c.s_id, c.t, c.t_id]
            msg = ('Added connection ' + c.s.label + ' (' + c.s_id + ') -> ' +
                   c.t.label + ' (' + c.t_id + ') to network.')
            logging.debug(msg)
            # set status "checked" to false, if conneciton is added to network.
            self.checked = False

    def del_conns(self, *args):
        """
        Remove one or more connections from the network.

        Parameters
        ----------
        c : tespy.connections.connection
            The connection to be removed from the network, connections objects
            ci :code:`del_conns(c1, c2, c3, ...)`.
        """
        for c in args:
            self.conns = self.conns.drop(c)
            msg = ('Deleted connection ' + c.s.label + ' (' + c.s_id +
                   ') -> ' + c.t.label + ' (' + c.t_id + ') from network.')
            logging.debug(msg)
        # set status "checked" to false, if conneciton is deleted from network.
        self.checked = False

    def check_conns(self):
        r"""Check connections for multiple usage of inlets or outlets."""
        dub = self.conns.loc[self.conns.duplicated(['s', 's_id']) == True]
        for c in dub.index:
            targets = ''
            for conns in self.conns[(self.conns.s == c.s) &
                                    (self.conns.s_id == c.s_id)].index:
                targets += conns.t.label + ' (' + conns.t_id + '); '

            msg = ('The source ' + c.s.label + ' (' + c.s_id + ') is attached '
                   'to more than one target: ' + targets[:-2] + '. '
                   'Please check your network.')
            logging.error(msg)
            raise hlp.TESPyNetworkError(msg)

        dub = self.conns.loc[self.conns.duplicated(['t', 't_id']) == True]
        for c in dub.index:
            sources = ''
            for conns in self.conns[(self.conns.t == c.t) &
                                    (self.conns.t_id == c.t_id)].index:
                sources += conns.s.label + ' (' + conns.s_id + '); '

            msg = ('The target ' + c.t.label + ' (' + c.t_id + ') is attached '
                   'to more than one source: ' + sources[:-2] + '. '
                   'Please check your network.')
            logging.error(msg)
            raise hlp.TESPyNetworkError(msg)

    def add_busses(self, *args):
        r"""
        Add one or more busses to the network.

        Parameters
        ----------
        b : tespy.connections.bus
            The bus to be added to the network, bus objects bi
            :code:`add_busses(b1, b2, b3, ...)`.
        """
        for b in args:
            if self.check_busses(b):
                self.busses[b.label] = b
                msg = 'Added bus ' + b.label + ' to network.'
                logging.debug(msg)

    def del_busses(self, *args):
        r"""
        Remove one or more busses from the network.

        Parameters
        ----------
        b : tespy.connections.bus
            The bus to be removed from the network, bus objects bi
            :code:`add_busses(b1, b2, b3, ...)`.
        """
        for b in args:
            if b in self.busses.values():
                del self.busses[b.label]
                msg = 'Deleted bus ' + b.label + ' from network.'
                logging.debug(msg)

    def check_busses(self, b):
        r"""
        Checksthe busses to be added for type, duplicates and identical labels.

        Parameters
        ----------
        b : tespy.connections.bus
            The bus to be checked.
        """
        if isinstance(b, con.bus):
            if len(self.busses) > 0:
                if b in self.busses.values():
                    msg = ('Network contains the bus ' + b.label + ' (' +
                           str(b) + ') already.')
                    logging.error(msg)
                    raise hlp.TESPyNetworkError(msg)
                elif b.label in self.busses.keys():
                    msg = ('Network already has a bus with the name ' +
                           b.label + '.')
                    logging.error(msg)
                    raise hlp.TESPyNetworkError(msg)
                else:
                    return True
            else:
                return True
        else:
            msg = 'Only objects of type bus are allowed in *args.'
            logging.error(msg)
            raise TypeError(msg)

    def check_network(self):
        r"""Check if components are connected properly within the network."""
        self.check_conns()
        # get unique components in connections dataframe
        comps = pd.unique(self.conns[['s', 't']].values.ravel())
        # build the dataframe for components
        self.init_components(comps)
        # count number of incoming and outgoing connections and compare to
        # expected values
        for comp in self.comps.index:
            num_o = (self.conns[['s', 't']] == comp).sum().s
            num_i = (self.conns[['s', 't']] == comp).sum().t
            if num_o != comp.num_o:
                msg = (comp.label + ' is missing ' + str(comp.num_o - num_o) +
                       ' outgoing connections. Make sure all outlets are '
                       ' connected and all connections have been added to the '
                       'network.')
                logging.error(msg)
                # raise an error in case network check is unsuccesful
                raise hlp.TESPyNetworkError(msg)
            elif num_i != comp.num_i:
                msg = (comp.label + ' is missing ' + str(comp.num_i - num_i) +
                       ' incoming connections. Make sure all inlets are '
                       ' connected and all connections have been added to the '
                       'network.')
                logging.error(msg)
                # raise an error in case network check is unsuccesful
                raise hlp.TESPyNetworkError(msg)

        # network checked
        self.checked = True
        msg = 'Networkcheck successful.'
        logging.info(msg)

    def init_components(self, comps):
        r"""
        Set up a dataframe for the network's components.

        Additionally, check, if all components have unique labels.

        Note
        ----
        The dataframe for the components is derived from the network's
        connections. Thus it does not hold any additional information, the
        dataframe is used to simplify the code, only.
        """
        self.comps = pd.DataFrame(index=comps, columns=['i', 'o'])

        labels = []
        for comp in self.comps.index:
            # get for incoming and outgoing connections of a component
            s = self.conns[self.conns.s == comp]
            s = s.s_id.sort_values().index
            t = self.conns[self.conns.t == comp]
            t = t.t_id.sort_values().index
            self.comps.loc[comp] = [t, s]
            # save the incoming and outgoing as well as the number of
            # connections as component attribute
            comp.inl = t.tolist()
            comp.outl = s.tolist()
            comp.num_i = len(comp.inlets())
            comp.num_o = len(comp.outlets())
            labels += [comp.label]

            # save the connection locations to the components
            comp.conn_loc = []
            for c in comp.inl + comp.outl:
                comp.conn_loc += [self.conns.index.get_loc(c)]

        # check for duplicates in the component labels
        if len(labels) != len(list(set(labels))):
            duplicates = [item for item, count in
                          Counter(labels).items() if count > 1]
            msg = ('All Components must have unique labels, duplicates are: '
                   + str(duplicates) + '.')
            logging.error(msg)
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
        if len(self.conns) == 0:
            msg = ('No connections have been added to the network, please '
                   'make sure to add your connections with the '
                   '.add_conns() method.')
            logging.error(msg)
            raise hlp.TESPyNetworkError(msg)

        if len(self.fluids) == 0:
            msg = ('Network has no fluids, please specify a list with fluids '
                   'on network creation.')
            logging.error(msg)
            raise hlp.TESPyNetworkError(msg)

        if self.mode == 'offdesign':
            self.redesign = True
            if self.design_path is None:
                # must provide design_path
                msg = ('Please provide \'design_path\' for every offdesign '
                       'calculation.')
                logging.error(msg)
                raise hlp.TESPyNetworkError(msg)
            else:
                # load design case
                if self.new_design is True:
                    self.init_offdesign_params()

                self.init_offdesign()
        else:
            # load design case
            self.init_design()

        # generic fluid initialisation
        self.init_fluids()
        # generic fluid property initialisation
        self.init_properties()

        # starting values from init path
        if self.init_path is not None:
            self.init_csv()

        msg = 'Network initialised.'
        logging.info(msg)

    def init_design(self):
        r"""
        Initialise a design calculation.

        Offdesign parameters are unset, design parameters are set. If
        :code:`local_offdesign` is True for connections or components, the
        design point information are read from the .csv-files in the respective
        :code:`design_path`. In this case, the design values are unset, the
        offdesign values set.
        """
        # connections
        for c in self.conns.index:
            # read design point information of connections with
            # local_offdesign activated from their respective design path
            if c.local_offdesign is True:
                if c.design_path is None:
                    msg = ('The parameter local_offdesign is True for the '
                           'connection ' + c.s.label + '(' + c.s_id + ') -> ' +
                           c.t.label + '(' + c.t_id + '), an individual '
                           'design_path must be specified in this case!')
                    logging.error(msg)
                    raise hlp.TESPyNetworkError(msg)

                # unset design parameters
                for var in c.design:
                    c.get_attr(var).set_attr(val_set=False)
                # set offdesign parameters
                for var in c.offdesign:
                    c.get_attr(var).set_attr(val_set=True)

                # read design point information
                path = hlp.modify_path_os(c.design_path + '/conn.csv')
                msg = ('Reading individual design point information for '
                       'connection ' + c.s.label + '(' + c.s_id + ') -> ' +
                       c.t.label + '(' + c.t_id + ') from path ' +
                       path + '.')
                logging.debug(msg)
                df = pd.read_csv(path, index_col=0, delimiter=';', decimal='.')

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
                if self.redesign is True:
                    for var in c.design:
                        c.get_attr(var).set_attr(val_set=True)

                    for var in c.offdesign:
                        c.get_attr(var).set_attr(val_set=False)

        # unset design values for busses
        for b in self.busses.values():
            for cp in b.comps.index:
                b.comps.loc[cp].P_ref = np.nan

        series = pd.Series()
        for cp in self.comps.index:
            # read design point information of components with
            # local_offdesign activated from their respective design path
            if cp.local_offdesign is True:
                if cp.design_path is not None:
                    # get type of component (class name)
                    c = cp.__class__.__name__
                    # read design point information
                    path = hlp.modify_path_os(cp.design_path + '/comps/' + c +
                                              '.csv')
                    df = pd.read_csv(path, sep=';', decimal='.',
                                     converters={
                                             'busses': ast.literal_eval,
                                             'bus_P_ref': ast.literal_eval
                                             })
                    df.set_index('label', inplace=True)
                    # write data
                    self.init_comp_design_params(cp, df.loc[cp.label])

                # unset design parameters
                for var in cp.design:
                    cp.get_attr(var).set_attr(is_set=False)

                # set offdesign parameters
                switched = False
                msg = 'Set component attributes '

                for var in cp.offdesign:
                    # set variables provided in .offdesign attribute
                    data = cp.get_attr(var)
                    data.set_attr(is_set=True)

                    # take nominal values from design point
                    if isinstance(data, dc.dc_cp):
                        cp.get_attr(var).val = cp.get_attr(var).design
                        switched = True
                        msg += var + ', '

                if switched:
                    msg = (msg[:-2] + ' to design value at component ' +
                           cp.label + '.')
                    logging.debug(msg)

                cp.new_design = False

            else:
                # switch connections to design mode
                if self.redesign is True:
                    for var in cp.design:
                        cp.get_attr(var).set_attr(is_set=True)

                    for var in cp.offdesign:
                        cp.get_attr(var).set_attr(is_set=False)

                cp.set_parameters(self.mode, series)

            # component initialisation
            cp.comp_init(self)

    def init_offdesign_params(self):
        r"""
        Read design point information from specified :code:`design_path`.

        If a :code:`design_path` has been specified individually for components
        or connections, the data will be read from the specified individual
        path instead.

        Note
        ----
        The methods
        :func:`tespy.networks.network.init_comp_design_params` (components) and
        the :func:`tespy.networks.network.init_conn_design_params`
        (connections) handle the parameter specification.
        """
        # components without any parameters
        not_required = ['source', 'sink', 'node', 'merge', 'splitter',
                        'separator', 'drum', 'subsystem_interface']
        # fetch all components, reindex with label
        cp_sort = self.comps.copy()
        # get class name
        cp_sort['cp'] = cp_sort.apply(network.get_class_base, axis=1)
        cp_sort['label'] = cp_sort.apply(network.get_props, axis=1,
                                         args=('label',))
        cp_sort['comp'] = cp_sort.index
        cp_sort.set_index('label', inplace=True)

        # iter through unique types of components (class names)
        for c in cp_sort.cp.unique():
            if c not in not_required:
                path = hlp.modify_path_os(self.design_path +
                                          '/comps/' + c + '.csv')
                msg = ('Reading design point information for components of '
                       'type ' + c + ' from path ' + path + '.')
                logging.debug(msg)

                # read data
                df = pd.read_csv(path, sep=';', decimal='.',
                                 converters={'busses': ast.literal_eval,
                                             'bus_P_ref': ast.literal_eval})
                df.set_index('label', inplace=True)
                # iter through all components of this type and set data
                for c_label in df.index:
                    comp = cp_sort.loc[c_label].comp
                    # read data of components with individual design_path
                    if comp.design_path is not None:
                        path_c = hlp.modify_path_os(comp.design_path +
                                                    '/comps/' + c + '.csv')
                        df_c = pd.read_csv(path_c, sep=';', decimal='.',
                                           converters={
                                                 'busses': ast.literal_eval,
                                                 'bus_P_ref': ast.literal_eval
                                                 })
                        df_c.set_index('label', inplace=True)
                        data = df_c.loc[comp.label]

                    else:
                        data = df.loc[comp.label]

                    # write data to components
                    self.init_comp_design_params(comp, data)

        msg = 'Read design point information for components.'
        logging.debug(msg)

        # read connection design point information
        path = hlp.modify_path_os(self.design_path + '/conn.csv')
        df = pd.read_csv(path, index_col=0, delimiter=';', decimal='.')
        msg = ('Reading design point information for connections from path ' +
               path + '.')
        logging.debug(msg)

        # iter through connections
        for c in self.conns.index:

            # read data of connections with individual design_path
            if c.design_path is not None:
                path_c = hlp.modify_path_os(c.design_path + '/conn.csv')
                msg = ('Reading individual design point information for '
                       'connection ' + c.s.label + '(' + c.s_id + ') -> ' +
                       c.t.label + '(' + c.t_id + ') from path ' +
                       path_c + '.')
                logging.debug(msg)
                df_c = pd.read_csv(path_c, index_col=0,
                                   delimiter=';', decimal='.')

                # write data
                self.init_conn_design_params(c, df_c)

            else:
                # write data
                self.init_conn_design_params(c, df)

        msg = 'Read design point information for connections.'
        logging.debug(msg)

    def init_comp_design_params(self, component, data):
        r"""
        Write design point information to components.

        Parameters
        ----------
        component : tespy.components.components.component
            Write design point information to this component.

        data : pandas.core.series.Series or pandas.core.frame.DataFrame
            Design point information.
        """
        # write component design data
        component.set_parameters(self.mode, data)
        # write design values to busses
        i = 0
        for b in data.busses:
            bus = self.busses[b].comps
            bus.loc[component].P_ref = data.bus_P_ref[i]
            i += 1

    def init_conn_design_params(self, c, df):
        r"""
        Write design point information to connections.

        Parameters
        ----------
        c : tespy.connections.connection
            Write design point information to this connection.

        df : pandas.core.frame.DataFrame
            Dataframe containing design point information.
        """
        # match connection (source, source_id, target, target_id) on
        # connection objects of design file
        conn = (df.loc[df['s'].isin([c.s.label]) & df['t'].isin([c.t.label]) &
                       df['s_id'].isin([c.s_id]) & df['t_id'].isin([c.t_id])])

        if len(conn.index) > 0:
            # read connection information
            conn_id = conn.index[0]
            c.m.design = df.loc[conn_id].m * self.m[df.loc[conn_id].m_unit]
            c.p.design = df.loc[conn_id].p * self.p[df.loc[conn_id].p_unit]
            c.h.design = df.loc[conn_id].h * self.h[df.loc[conn_id].h_unit]
            c.v.design = df.loc[conn_id].v * self.v[df.loc[conn_id].v_unit]
            c.x.design = df.loc[conn_id].x
            c.T.design = ((df.loc[conn_id]['T'] +
                           self.T[df.loc[conn_id].T_unit][0]) *
                          self.T[df.loc[conn_id].T_unit][1])
            c.Td_bp.design = (df.loc[conn_id].Td_bp *
                              self.T[df.loc[conn_id].T_unit][1])
            for fluid in self.fluids:
                c.fluid.design[fluid] = df.loc[conn_id][fluid]
        else:
            # no matches in the connections of the network and the design files
            msg = ('Could not find connection ' + c.s.label + '(' + c.s_id +
                   ') -> ' + c.t.label + '(' + c.t_id + ') in design case. '
                   'Please, make sure no connections have been modified '
                   'or components have been relabeled for your offdesign '
                   'calculation.')
            logging.error(msg)
            raise hlp.TESPyNetworkError(msg)

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
        for cp in self.comps.index:
            if cp.local_design is False:
                # unset variables provided in .design attribute
                for var in cp.design:
                    cp.get_attr(var).set_attr(is_set=False)

                switched = False
                msg = 'Set component attributes '

                for var in cp.offdesign:
                    # set variables provided in .offdesign attribute
                    data = cp.get_attr(var)
                    data.set_attr(is_set=True)

                    # take nominal values from design point
                    if isinstance(data, dc.dc_cp):
                        cp.get_attr(var).val = cp.get_attr(var).design
                        switched = True
                        msg += var + ', '

                if switched:
                    msg = (msg[:-2] + ' to design value at component ' +
                           cp.label + '.')
                    logging.debug(msg)

            # start component initialisation
            cp.comp_init(self)
            cp.new_design = False

        msg = 'Switched components from design to offdesign.'
        logging.debug(msg)

        for c in self.conns.index:
            if c.local_design is False:
                # switch connections to offdesign mode
                for var in c.design:
                    c.get_attr(var).set_attr(val_set=False)
                    c.get_attr(var).set_attr(ref_set=False)

                for var in c.offdesign:
                    c.get_attr(var).set_attr(val_set=True)
                    c.get_attr(var).val_SI = c.get_attr(var).design

                c.new_design = False

        msg = 'Switched connections from design to offdesign.'
        logging.debug(msg)

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
        # iterate over connectons, create ordered dicts
        for c in self.conns.index:
            c.conn_loc = self.conns.index.get_loc(c)
            tmp = c.fluid.val.copy()
            tmp0 = c.fluid.val0.copy()
            tmp_set = c.fluid.val_set.copy()
            c.fluid.val = OrderedDict()
            c.fluid.val0 = OrderedDict()
            c.fluid.val_set = OrderedDict()

            # if the number if fluids is one
            if len(self.fluids) == 1:
                c.fluid.val[self.fluids[0]] = 1
                c.fluid.val0[self.fluids[0]] = 1

                if self.fluids[0] in tmp_set.keys():
                    c.fluid.val_set[self.fluids[0]] = tmp_set[self.fluids[0]]
                else:
                    c.fluid.val_set[self.fluids[0]] = False

                # jump to next connection
                continue

            for fluid in self.fluids:

                if fluid in tmp.keys() and fluid in tmp_set.keys():
                    # if fluid in keys and is_set
                    c.fluid.val[fluid] = tmp[fluid]
                    c.fluid.val0[fluid] = tmp[fluid]
                    c.fluid.val_set[fluid] = tmp_set[fluid]

                # if there is a starting value
                elif fluid in tmp0.keys():
                    if fluid not in tmp_set.keys():
                        c.fluid.val[fluid] = tmp0[fluid]
                        c.fluid.val0[fluid] = tmp0[fluid]
                        c.fluid.val_set[fluid] = False

                # if fluid not in keys
                else:
                    c.fluid.val[fluid] = 0
                    c.fluid.val0[fluid] = 0
                    c.fluid.val_set[fluid] = False

        # stop fluid propagation for single fluid networks and
        # for offdesign cases, as good starting values are available
        if self.mode == 'offdesign' or len(self.fluids) == 1:
            msg = 'Fluid initialisation done.'
            logging.debug(msg)
            return

        # fluid propagation from set values
        for c in self.conns.index:
            if any(c.fluid.val_set.values()):
                self.init_target(c, c.t)
                self.init_source(c, c.s)

        # fluid propagation for components
        for cp in self.comps.index:
            # combustion chamber
            if isinstance(cp, combustion_chamber):
                cp.initialise_fluids(self)
                for c in self.comps.loc[cp].o:
                    self.init_target(c, c.t)

            # combustion chamber
            elif isinstance(cp, water_electrolyzer):
                cp.initialise_fluids(self)
                for c in self.comps.loc[cp].o:
                    self.init_target(c, c.t)

            # other components (node, merge)
            else:
                cp.initialise_fluids(self)
                for c in self.comps.loc[cp].o:
                    self.init_target(c, c.t)

        msg = 'Fluid initialisation done.'
        logging.debug(msg)

    def init_target(self, c, start):
        r"""
        Propagate the fluids towards connection's target in recursion.

        The propagation stops, if the component is a

        - sink,
        - merge,
        - combustion chamber or
        - cycle_closer.

        Parameters
        ----------
        c : tespy.connections.connection
            Connection to initialise.

        start : tespy.connections.connection
            This connection is the fluid propagation starting point.
            The starting connection is saved to prevent infinite looping.
        """
        if ((len(c.t.inlets()) == 1 and len(c.t.outlets()) == 1 and
                not isinstance(c.t, cycle_closer)) or
                isinstance(c.t, heat_exchanger) or
                isinstance(c.t, subsystem_interface)):

            outc = pd.DataFrame()
            outc['s'] = self.conns.s == c.t
            outc['s_id'] = self.conns.s_id == c.t_id.replace('in', 'out')
            conn, cid = outc['s'] == True, outc['s_id'] == True
            outconn = outc.index[conn & cid][0]

            for fluid, x in c.fluid.val.items():
                if (outconn.fluid.val_set[fluid] is False and
                        outconn.good_starting_values is False):
                    outconn.fluid.val[fluid] = x

            self.init_target(outconn, start)

        if isinstance(c.t, splitter):
            for outconn in self.comps.loc[c.t].o:
                for fluid, x in c.fluid.val.items():
                    if (outconn.fluid.val_set[fluid] is False and
                            outconn.good_starting_values is False):
                        outconn.fluid.val[fluid] = x

                self.init_target(outconn, start)

        if isinstance(c.t, water_electrolyzer):
            if c == self.comps.loc[c.t].i[0]:
                outconn = self.comps.loc[c.t].o[0]

                for fluid, x in c.fluid.val.items():
                    if (outconn.fluid.val_set[fluid] is False and
                            outconn.good_starting_values is False):
                        outconn.fluid.val[fluid] = x

        if isinstance(c.t, combustion_engine):
            for outconn in self.comps.loc[c.t].o[:2]:
                for fluid, x in c.fluid.val.items():
                    if (outconn.fluid.val_set[fluid] is False and
                            outconn.good_starting_values is False):
                        outconn.fluid.val[fluid] = x

                self.init_target(outconn, start)

        if isinstance(c.t, drum) and c.t != start:
            start = c.t
            for outconn in self.comps.loc[c.t].o:
                for fluid, x in c.fluid.val.items():
                    if (outconn.fluid.val_set[fluid] is False and
                            outconn.good_starting_values is False):
                        outconn.fluid.val[fluid] = x

                self.init_target(outconn, start)

    def init_source(self, c, start):
        r"""
        Propagate the fluids towards connection's source in recursion.

        The propagation stops, if the component is a

        - source,
        - combustion chamber or
        - cycle_closer.

        Parameters
        ----------
        c : tespy.connections.connection
            Connection to initialise.

        start : tespy.connections.connection
            This connection is the fluid propagation starting point.
            The starting connection is saved to prevent infinite looping.
        """
        if ((len(c.s.inlets()) == 1 and len(c.s.outlets()) == 1 and
                not isinstance(c.s, cycle_closer)) or
                isinstance(c.s, heat_exchanger) or
                isinstance(c.s, subsystem_interface)):

            inc = pd.DataFrame()
            inc['t'] = self.conns.t == c.s
            inc['t_id'] = self.conns.t_id == c.s_id.replace('out', 'in')
            conn, cid = inc['t'] == True, inc['t_id'] == True
            inconn = inc.index[conn & cid][0]

            for fluid, x in c.fluid.val.items():
                if (inconn.fluid.val_set[fluid] is False and
                        inconn.good_starting_values is False):
                    inconn.fluid.val[fluid] = x

            self.init_source(inconn, start)

        if isinstance(c.s, splitter):
            for inconn in self.comps.loc[c.s].i:
                for fluid, x in c.fluid.val.items():
                    if (inconn.fluid.val_set[fluid] is False and
                            inconn.good_starting_values is False):
                        inconn.fluid.val[fluid] = x

                self.init_source(inconn, start)

        if isinstance(c.s, merge):
            for inconn in self.comps.loc[c.s].i:
                for fluid, x in c.fluid.val.items():
                    if (inconn.fluid.val_set[fluid] is False and
                            inconn.good_starting_values is False):
                        inconn.fluid.val[fluid] = x

                self.init_source(inconn, start)

        if isinstance(c.s, combustion_engine):
            for inconn in self.comps.loc[c.s].i[:2]:
                for fluid, x in c.fluid.val.items():
                    if (inconn.fluid.val_set[fluid] is False and
                            inconn.good_starting_values is False):
                        inconn.fluid.val[fluid] = x

                self.init_source(inconn, start)

        if isinstance(c.s, drum) and c.s != start:
            start = c.s
            for inconn in self.comps.loc[c.s].i:
                for fluid, x in c.fluid.val.items():
                    if (inconn.fluid.val_set[fluid] is False and
                            inconn.good_starting_values is False):
                        inconn.fluid.val[fluid] = x

                self.init_source(inconn, start)

    def init_properties(self):
        r"""
        Initialise the fluid properties on every connection of the network.

        - Sets standard values for :code:`m0, p0, h0` if not user specified
        - Sets :code:`var = var0` if var_set is False.
        - Initialises reference objects.
        - Sets initial values for enthalpy at given vapour mass fraction or
          temperature.
        """
        # fluid properties
        for c in self.conns.index:
            if self.init_previous is False:
                c.good_starting_values = False
            for key in ['m', 'p', 'h', 'T', 'x', 'v', 'Td_bp']:
                if c.get_attr(key).unit_set is False and key != 'x':
                    if key == 'Td_bp':
                        c.get_attr(key).unit = self.get_attr('T_unit')
                    else:
                        c.get_attr(key).unit = self.get_attr(key + '_unit')
                if (key not in ['T', 'x', 'v', 'Td_bp'] and
                        c.get_attr(key).val_set is False):
                    if c.good_starting_values is False:
                        self.init_val0(c, key)
                    c.get_attr(key).val_SI = (
                            c.get_attr(key).val0 * self.get_attr(key)[
                                    c.get_attr(key).unit])
                elif key not in c.offdesign:
                    if (key not in ['T', 'x', 'v', 'Td_bp'] and
                            c.get_attr(key).val_set is True):
                        c.get_attr(key).val_SI = (
                                c.get_attr(key).val * self.get_attr(key)[
                                        c.get_attr(key).unit])
                    elif (key == 'T' and c.T.val_set is True):
                        c.T.val_SI = ((c.T.val + self.T[c.T.unit][0]) *
                                      self.T[c.T.unit][1])
                    elif key == 'Td_bp' and c.Td_bp.val_set is True:
                        c.Td_bp.val_SI = c.Td_bp.val * self.T[c.T.unit][1]
                    elif key == 'x' and c.x.val_set is True:
                        c.x.val_SI = c.x.val
                    elif key == 'v' and c.v.val_set is True:
                        c.v.val_SI = c.v.val * self.v[c.v.unit]

        msg = ('Retrieved generic starting values and specified SI-values of '
               'connection parameters.')
        logging.debug(msg)

        # improved starting values for referenced connections,
        # specified vapour content values, temperature values as well as
        # subccooling/overheating and state specification
        for c in self.conns.index:
            # starting values for fluid properties with referenced objects
            for key in ['m', 'p', 'h', 'T']:
                if (c.get_attr(key).ref_set and
                        c.get_attr(key).val_set is False):
                    c.get_attr(key).val_SI = (
                            c.get_attr(key).ref.obj.get_attr(key).val_SI *
                            c.get_attr(key).ref.f + c.get_attr(key).ref.d)

            # starting values for specified vapour content or temperature
            if c.x.val_set and not c.h.val_set:
                c.h.val_SI = fp.h_mix_pQ(c.to_flow(), c.x.val_SI)

            if c.T.val_set and not c.h.val_set:
                try:
                    c.h.val_SI = fp.h_mix_pT(c.to_flow(), c.T.val_SI)
                except ValueError:
                    pass

            # starting values for specified subcooling/overheating
            # and state specification
            if ((c.Td_bp.val_set is True or c.state.is_set is True) and
                    c.h.val_set is False):
                if ((c.Td_bp.val_SI > 0 and c.Td_bp.val_set is True) or
                        (c.state.val == 'g' and c.state.is_set is True)):
                    h = fp.h_mix_pQ(c.to_flow(), 1)
                    if c.h.val_SI < h:
                        c.h.val_SI = h * 1.001
                elif ((c.Td_bp.val_SI < 0 and c.Td_bp.val_set is True) or
                      (c.state.val == 'l' and c.state.is_set is True)):
                    h = fp.h_mix_pQ(c.to_flow(), 0)
                    if c.h.val_SI > h:
                        c.h.val_SI = h * 0.999

        msg = 'Generic fluid property specification complete.'
        logging.debug(msg)

    def init_val0(self, c, key):
        r"""
        Set starting values for fluid properties.

        The component classes provide generic starting values for their inlets
        and outlets.

        Parameters
        ----------
        c : tespy.connections.connection
            Connection to initialise.
        """
        if np.isnan(c.get_attr(key).val0):
            # starting value for mass flow is 1 kg/s
            if key == 'm':
                c.get_attr(key).val0 = 1

            # generic starting values for pressure and enthalpy
            else:
                # retrieve starting values from component information
                val_s = c.s.initialise_source(c, key)
                val_t = c.t.initialise_target(c, key)

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
                c.get_attr(key).val0 = (
                        c.get_attr(key).val0 / self.get_attr(key)[
                                self.get_attr(key + '_unit')])

    def init_csv(self):
        r"""
        Read starting values from specified init_path.

        Note
        ----
        This method loads fluid property and fluid components starting values
        using the :code:`init_file` as input file.
        """
        # match connection (source, source_id, target, target_id) on
        # connection objects of design file

        path = hlp.modify_path_os(self.init_path + '/conn.csv')

        df = pd.read_csv(path, index_col=0, delimiter=';', decimal='.')
        for c in self.conns.index:
            conn = (df.loc[df['s'].isin([c.s.label]) &
                           df['t'].isin([c.t.label]) &
                           df['s_id'].isin([c.s_id]) &
                           df['t_id'].isin([c.t_id])])
            if len(conn.index) > 0:
                conn_id = conn.index[0]
                # overwrite SI-values with values from init_file,
                # except user specified values
                if c.m.val_set is False:
                    c.m.val_SI = (df.loc[conn_id].m *
                                  self.m[df.loc[conn_id].m_unit])
                if c.p.val_set is False:
                    c.p.val_SI = (df.loc[conn_id].p *
                                  self.p[df.loc[conn_id].p_unit])
                if c.h.val_set is False:
                    c.h.val_SI = (df.loc[conn_id].h *
                                  self.h[df.loc[conn_id].h_unit])
                for fluid in self.fluids:
                    if c.fluid.val_set[fluid] is False:
                        c.fluid.val[fluid] = df.loc[conn_id][fluid]

                # overwrite starting values
                c.m.val0 = c.m.val_SI / self.m[c.m.unit]
                c.p.val0 = c.p.val_SI / self.p[c.p.unit]
                c.h.val0 = c.h.val_SI / self.h[c.h.unit]
                c.fluid.val0 = c.fluid.val.copy()
                c.good_starting_values = True
            else:
                msg = ('Could not find connection ' + c.s.label + ' (' +
                       c.s_id + ') -> ' + c.t.label + ' (' + c.t_id +
                       ') in .csv-file.')
                logging.debug(msg)

        msg = 'Specified starting values from init_path.'
        logging.debug(msg)

    def solve(self, mode, init_path=None, design_path=None,
              max_iter=50, init_only=False, init_previous=True):
        r"""
        Solve the network.

        - Check network consistency.
        - Initialise calculation and preprocessing.
        - Perform actual calculation.
        - Postprocessing.

        Parameters
        ----------
        mode : str
            Choose from 'design' and 'offdesign'.

        init_path : str
            Path to the folder, where your network was saved to, e. g.
            saving to :code:`nw.save('myplant/tests')` would require loading
            from :code:`init_path='myplant/tests'`.

        design_path : str
            Path to the folder, where your network's design case was saved to,
            e. g. saving to :code:`nw.save('myplant/tests')` would require
            loading from :code:`design_path='myplant/tests'`.

        max_iter : int
            Maximum number of iterations before calculation stops, default: 50.

        init_only : boolean
            Perform initialisation only, default: :code:`False`.

        init_previous : boolean
            Initialise the calculation with values from the previous
            calculation, default: :code:`True`.

        Note
        ----
        For more information on the solution process have a look at the online
        documentation at tespy.readthedocs.io in the section "using TESPy".
        """
        if self.design_path == design_path and design_path is not None:
            for c in self.conns.index:
                if c.new_design is True:
                    self.new_design = True
                    break

            for cp in self.comps.index:
                if cp.new_design is True:
                    self.new_design = True
                    break

            self.new_design = False
        else:
            self.new_design = True

        self.init_path = init_path
        self.design_path = design_path
        self.max_iter = max_iter
        self.init_previous = init_previous

        if mode != 'offdesign' and mode != 'design':
            msg = 'Mode must be \'design\' or \'offdesign\'.'
            logging.error(msg)
            raise ValueError(msg)
        else:
            self.mode = mode

        msg = ('Solver properties: '
               'mode=' + self.mode +
               ', init_path=' + str(self.init_path) +
               ', design_path=' + str(self.design_path) +
               ', max_iter=' + str(max_iter) +
               ', init_only=' + str(init_only))
        logging.debug(msg)

        if not self.checked:
            self.check_network()

        msg = ('Network properties: '
               'number of components=' + str(len(self.comps.index)) +
               ', number of connections=' + str(len(self.conns.index)) +
               ', number of busses=' + str(len(self.busses)))
        logging.debug(msg)

        self.initialise()

        if init_only:
            return

        self.res = np.array([])

        msg = 'Starting solver.'
        logging.info(msg)

        self.iter = 0
        # number of variables per connection
        self.num_conn_vars = len(self.fluids) + 3

        # check for network determination
        self.solve_determination()

        self.solve_loop()

        if not self.progress:
            msg = ('The solver does not seem to make any progress, aborting '
                   'calculation. Residual value is '
                   '{:.2e}'.format(norm(self.vec_res)) + '. This frequently '
                   'happens, if the solver pushes the fluid properties out of '
                   'their feasible range.')
            logging.warning(msg)

        if self.lin_dep:
            msg = ('Singularity in jacobian matrix, calculation aborted! Make '
                   'sure your network does not have any linear dependencies '
                   'in the parametrisation. Other reasons might be\n'
                   '-> given Temperature with given pressure in two phase '
                   'region, try setting enthalpy instead or '
                   'provide accurate starting value for pressure.\n'
                   '-> given logarithmic temperature differences '
                   'or kA-values for heat exchangers, \n'
                   '-> support better starting values.\n'
                   '-> bad starting value for fuel mass flow of '
                   'combustion chamber, provide small (near to zero, '
                   'but not zero) starting value.')
            logging.error(msg)
            return

        self.postprocessing()
        fp.memorise.del_memory(self.fluids)

        if not self.progress:
            return

        msg = 'Calculation complete.'
        logging.info(msg)

    def solve_loop(self):
        r"""Loop of the newton algorithm."""
        self.start_time = time()
        self.progress = True

        if self.iterinfo:
            self.print_iterinfo('start')

        for self.iter in range(self.max_iter):

            self.vec_z_filter = np.absolute(self.vec_z) < err ** 2
            self.solve_control()
            self.res = np.append(self.res, norm(self.vec_res))

            if self.iterinfo:
                self.print_iterinfo('solving')

            if ((self.iter > 3 and self.res[-1] < err ** 0.5) or
                    self.lin_dep):
                break

            if self.iter > 20:
                if (all(self.res[(self.iter - 3):] >= self.res[-2] * 0.95) and
                        self.res[-1] >= self.res[-2] * 0.95):
                    self.progress = False
                    break

        self.end_time = time()

        self.print_iterinfo('end')

        if self.iter == self.max_iter - 1:
            msg = ('Reached maximum iteration count (' + str(self.max_iter) +
                   '), calculation stopped. Residual value is '
                   '{:.2e}'.format(norm(self.vec_res)))
            logging.warning(msg)

    def solve_determination(self):
        r"""Check, if the number of supplied parameters is sufficient."""
        # number of equations from components
        # component variables
        self.num_comp_vars = 0
        # component equations
        self.num_comp_eq = 0
        for cp in self.comps.index:
            self.num_comp_vars += cp.num_vars
            self.num_comp_eq += cp.num_eq

        msg = 'Number of component equations: ' + str(self.num_comp_eq) + '.'
        logging.debug(msg)

        # number of equations from connections
        n = 0
        for c in self.conns.index:
            n += [c.m.val_set, c.p.val_set, c.h.val_set, c.T.val_set,
                  c.x.val_set, c.v.val_set, c.Td_bp.val_set].count(True)
            n += [c.m.ref_set, c.p.ref_set, c.h.ref_set,
                  c.T.ref_set].count(True)
            n += list(c.fluid.val_set.values()).count(True)
            n += [c.fluid.balance].count(True)

        msg = 'Number of connection equations: ' + str(n) + '.'
        logging.debug(msg)
        self.num_conn_eq = n

        # number of equations from busses
        self.num_bus_eq = 0
        for b in self.busses.values():
            self.num_bus_eq += [b.P.is_set].count(True)

        msg = 'Number of bus equations: ' + str(self.num_bus_eq) + '.'
        logging.debug(msg)

        # total number of variables
        self.num_vars = (self.num_conn_vars * len(self.conns.index) +
                         self.num_comp_vars)

        self.vec_res = np.zeros([self.num_vars])
        self.vec_z = np.ones([self.num_vars])
        self.mat_deriv = np.zeros((self.num_vars, self.num_vars))

        msg = 'Total number of variables: ' + str(self.num_vars) + '.'
        logging.debug(msg)
        msg = 'Number of component variables: ' + str(self.num_comp_vars) + '.'
        logging.debug(msg)
        msg = ('Number of connection variables: ' +
               str(self.num_conn_vars * len(self.conns.index)) + '.')
        logging.debug(msg)

        n = self.num_comp_eq + self.num_conn_eq + self.num_bus_eq
        if n > self.num_vars:
            msg = ('You have provided too many parameters: ' +
                   str(self.num_vars) + ' required, ' + str(n) +
                   ' supplied. Aborting calculation!')
            logging.error(msg)
            raise hlp.TESPyNetworkError(msg)
        elif n < self.num_vars:
            msg = ('You have not provided enough parameters: '
                   + str(self.num_vars) + ' required, ' + str(n) +
                   ' supplied. Aborting calculation!')
            logging.error(msg)
            raise hlp.TESPyNetworkError(msg)

    def print_iterinfo(self, position):
        """Print out the convergence progress."""
        if position == 'start':
            if self.num_comp_vars == 0:
                # iterinfo printout without any custom variables
                msg = ('iter\t| residual | massflow | pressure | enthalpy | '
                       'fluid\n')
                msg += ('--------+----------+----------+----------+----------'
                        '+---------')

            else:
                # iterinfo printout with custom variables in network
                msg = ('iter\t| residual | massflow | pressure | enthalpy | '
                       'fluid    | custom\n')
                msg += ('--------+----------+----------+----------+----------'
                        '+----------+---------')

            print(msg)

        elif position == 'solving':
            vec = self.vec_z[0:-(self.num_comp_vars + 1)]
            msg = (str(self.iter + 1))
            if not self.lin_dep and not np.isnan(norm(self.vec_res)):
                msg += '\t| ' + '{:.2e}'.format(norm(self.vec_res))
                msg += ' | ' + '{:.2e}'.format(
                        norm(vec[0::self.num_conn_vars]))
                msg += ' | ' + '{:.2e}'.format(
                        norm(vec[1::self.num_conn_vars]))
                msg += ' | ' + '{:.2e}'.format(
                        norm(vec[2::self.num_conn_vars]))
                ls = []
                for f in range(len(self.fluids)):
                    ls += vec[3 + f::self.num_conn_vars].tolist()

                msg += ' | ' + '{:.2e}'.format(norm(ls))
                if self.num_comp_vars > 0:
                    msg += ' | ' + '{:.2e}'.format(norm(
                            self.vec_z[-self.num_comp_vars:]))

            else:
                if np.isnan(norm(self.vec_res)):
                    msg += '\t|      nan'.format(norm(self.vec_res))
                else:
                    msg += '\t| ' + '{:.2e}'.format(norm(self.vec_res))
                msg += ' |      nan'
                msg += ' |      nan'
                msg += ' |      nan'
                msg += ' |      nan'
                if self.num_comp_vars > 0:
                    msg += ' |      nan'

            print(msg)

        elif position == 'end':
            if self.iterinfo:
                if self.num_comp_vars == 0:
                    msg = ('--------+----------+----------+----------'
                           '+----------+---------')
                else:
                    msg = ('--------+----------+----------+----------'
                           '+----------+----------+---------')
                print(msg)

            msg = ('Total iterations: ' + str(self.iter + 1) + ', '
                   'Calculation time: ' +
                   str(round(self.end_time - self.start_time, 1)) + ' s, '
                   'Iterations per second: ' +
                   str(round((self.iter + 1) /
                             (self.end_time - self.start_time), 2)))
            logging.debug(msg)
            if self.iterinfo:
                print(msg)

        else:
            pass

    def matrix_inversion(self):
        """Invert matrix of derivatives and caluclate increment."""
        self.lin_dep = True
        try:
            self.vec_z = inv(self.mat_deriv).dot(-self.vec_res)
            self.lin_dep = False
        except np.linalg.linalg.LinAlgError:
            self.vec_z = self.vec_res * 0
            pass

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
        self.matrix_inversion()

        # check for linear dependency
        if self.lin_dep:
            return

        # add the increment
        i = 0
        for c in self.conns.index:
            # mass flow, pressure and enthalpy
            if c.m.val_set is False:
                c.m.val_SI += self.vec_z[i * (self.num_conn_vars)]
            if c.p.val_set is False:
                # this prevents negative pressures
                relax = max(1, -self.vec_z[i * (self.num_conn_vars) + 1] /
                            (0.5 * c.p.val_SI))
                c.p.val_SI += self.vec_z[i * (self.num_conn_vars) + 1] / relax
            if c.h.val_set is False:
                c.h.val_SI += self.vec_z[i * (self.num_conn_vars) + 2]

            # fluid vector (only if number of fluids is greater than 1)
            if len(self.fluids) > 1:
                j = 0
                for fluid in self.fluids:
                    # add increment
                    if c.fluid.val_set[fluid] is False:
                        c.fluid.val[fluid] += (
                                self.vec_z[i * (self.num_conn_vars) + 3 + j])

                    # keep mass fractions within [0, 1]
                    if c.fluid.val[fluid] < err:
                        c.fluid.val[fluid] = 0
                    if c.fluid.val[fluid] > 1 - err:
                        c.fluid.val[fluid] = 1

                    j += 1

            # check the fluid properties for physical ranges
            self.solve_check_props(c)
            i += 1

        # increment for the custom variables
        if self.num_comp_vars > 0:
            c_vars = 0
            for cp in self.comps.index:
                for var in cp.vars.keys():
                    pos = var.var_pos

                    # add increment
                    var.val += self.vec_z[
                        self.num_conn_vars * len(self.conns) + c_vars + pos]

                    # keep value within specified value range
                    if var.val < var.min_val:
                        var.val = var.min_val
                    if var.val > var.max_val:
                        var.val = var.max_val

                c_vars += cp.num_vars

        # second property check for first three iterations without an init_file
        if self.iter < 3:
            for cp in self.comps.index:
                cp.convergence_check(self)

            for c in self.conns.index:
                self.solve_check_props(c)

    def property_range_message(self, c, prop):
        r"""
        Return debugging message for fluid property range adjustments.

        Parameters
        ----------
        c : tespy.connections.connection
            Connection to check fluid properties.

        prop : str
            Fluid property.

        Returns
        -------
        msg : str
            Debugging message.
        """
        if prop == 'p':
            msg = 'Pressure '
        elif prop == 'h':
            msg = 'Enthalpy '
        elif prop == 'm':
            msg = 'Mass flow '
        else:
            msg = 'Unspecified '
        msg += ('out of fluid property range at connection ' +
                c.s.label + ' (' + c.s_id + ') -> ' + c.t.label + ' (' +
                c.t_id + ') adjusting value to ' +
                str(c.get_attr(prop).val_SI) + ' ' + self.SI_units[prop] + '.')
        return msg

    def solve_check_props(self, c):
        r"""
        Check for invalid fluid property values.

        Parameters
        ----------
        c : tespy.connections.connection
            Connection to check fluid properties.
        """
        fl = hlp.single_fluid(c.fluid.val)

        if isinstance(fl, str):
            # pressure
            if c.p.val_SI < fp.memorise.vrange[fl][0] and not c.p.val_set:
                c.p.val_SI = fp.memorise.vrange[fl][0] * 1.01
                logging.debug(self.property_range_message(c, 'p'))
            if c.p.val_SI > fp.memorise.vrange[fl][1] and not c.p.val_set:
                c.p.val_SI = fp.memorise.vrange[fl][1] * 0.99
                logging.debug(self.property_range_message(c, 'p'))

            # enthalpy
            f = 1.01
            try:
                hmin = fp.h_pT(c.p.val_SI, fp.memorise.vrange[fl][2] * f, fl)
            except ValueError:
                f = 1.1
                hmin = fp.h_pT(c.p.val_SI, fp.memorise.vrange[fl][2] * f, fl)

            hmax = fp.h_pT(c.p.val_SI, fp.memorise.vrange[fl][3] * 0.99, fl)
            if c.h.val_SI < hmin and not c.h.val_set:
                if hmin < 0:
                    c.h.val_SI = hmin / 1.05
                else:
                    c.h.val_SI = hmin * 1.05
                logging.debug(self.property_range_message(c, 'h'))
            if c.h.val_SI > hmax and not c.h.val_set:
                c.h.val_SI = hmax * 0.95
                logging.debug(self.property_range_message(c, 'h'))

            if ((c.Td_bp.val_set is True or c.state.is_set is True) and
                    c.h.val_set is False and self.iter < 3):
                if (c.Td_bp.val_SI > 0 or
                        (c.state.val == 'g' and c.state.is_set is True)):
                    h = fp.h_mix_pQ(c.to_flow(), 1)
                    if c.h.val_SI < h:
                        c.h.val_SI = h * 1.02
                        logging.debug(self.property_range_message(c, 'h'))
                elif (c.Td_bp.val_SI < 0 or
                      (c.state.val == 'l' and c.state.is_set is True)):
                    h = fp.h_mix_pQ(c.to_flow(), 0)
                    if c.h.val_SI > h:
                        c.h.val_SI = h * 0.98
                        logging.debug(self.property_range_message(c, 'h'))

        elif self.iter < 4 and c.good_starting_values is False:
            # pressure
            if c.p.val_SI <= self.p_range_SI[0] and not c.p.val_set:
                c.p.val_SI = self.p_range_SI[0]
                logging.debug(self.property_range_message(c, 'p'))
            if c.p.val_SI >= self.p_range_SI[1] and not c.p.val_set:
                c.p.val_SI = self.p_range_SI[1]
                logging.debug(self.property_range_message(c, 'p'))

            # enthalpy
            if c.h.val_SI < self.h_range_SI[0] and not c.h.val_set:
                c.h.val_SI = self.h_range_SI[0]
                logging.debug(self.property_range_message(c, 'h'))
            if c.h.val_SI > self.h_range_SI[1] and not c.h.val_set:
                c.h.val_SI = self.h_range_SI[1]
                logging.debug(self.property_range_message(c, 'h'))

            # temperature
            if c.T.val_set and not c.h.val_set:
                self.solve_check_temperature(c)

        # mass flow
        if c.m.val_SI <= self.m_range_SI[0] and not c.m.val_set:
            c.m.val_SI = self.m_range_SI[0]
            logging.debug(self.property_range_message(c, 'm'))
        if c.m.val_SI >= self.m_range_SI[1] and not c.m.val_set:
            c.m.val_SI = self.m_range_SI[1]
            logging.debug(self.property_range_message(c, 'm'))

    def solve_check_temperature(self, c):
        r"""
        Check if temperature is within user specified limits.

        Parameters
        ----------
        c : tespy.connections.connection
            Connection to check fluid properties.
        """
        hmin = fp.h_mix_pT(c.to_flow(), self.T_range_SI[0])
        hmax = fp.h_mix_pT(c.to_flow(), self.T_range_SI[1])

        if c.h.val_SI < hmin:
            if c.h.val_SI < 0:
                c.h.val_SI = hmin * 0.9
            else:
                c.h.val_SI = hmin * 1.1
            logging.debug(self.property_range_message(c, 'h'))

        if c.h.val_SI > hmax:
            c.h.val_SI = hmax * 0.95
            logging.debug(self.property_range_message(c, 'h'))

    def solve_components(self):
        r"""
        Calculate the residual and derivatives of component equations.

        - Iterate through components in network to get residuals and
          derivatives.
        - Place residual values in residual value vector of the network.
        - Place partial derivatives in jacobian matrix of the network.
        """
        # fetch component equation residuals and component partial derivatives
        vec_res = []
        sum_eq = 0
        c_var = 0
        for cp in self.comps.index:

            indexes = []
            for c in cp.conn_loc:
                start = c * self.num_conn_vars
                end = (c + 1) * self.num_conn_vars
                indexes += [np.arange(start, end)]

            cp.equations()
            cp.derivatives(self.vec_z_filter[np.array(indexes)])

            self.vec_res[sum_eq:sum_eq + cp.num_eq] = cp.vec_res
            deriv = cp.mat_deriv

            if deriv is not None:
                i = 0
                # place derivatives in jacobian matrix
                for loc in cp.conn_loc:
                    coll_s = loc * self.num_conn_vars
                    coll_e = (loc + 1) * self.num_conn_vars
                    self.mat_deriv[
                        sum_eq:sum_eq + cp.num_eq, coll_s:coll_e] = deriv[:, i]
                    i += 1

                # derivatives for custom variables
                for j in range(cp.num_vars):
                    coll = self.num_vars - self.num_comp_vars + c_var
                    self.mat_deriv[sum_eq:sum_eq + cp.num_eq, coll] = (
                        deriv[:, i + j, :1].transpose()[0])
                    c_var += 1

                sum_eq += cp.num_eq
            cp.it += 1

    def solve_connections(self):
        r"""
        Calculate the residual and derivatives of connection equations.

        - Iterate through connections in network to get residuals and
          derivatives.
        - Place residual values in residual value vector of the network.
        - Place partial derivatives in jacobian matrix of the network.

        Note
        ----
        **Equations**

        **mass flow, pressure and enthalpy**

        .. math::
            val = 0

        **temperatures**

        .. math::
            val = T_{j} - T \left( p_{j}, h_{j}, fluid_{j} \right)

        **volumetric flow**

        .. math::
            val = \dot{V}_{j} - v \left( p_{j}, h_{j} \right) \cdot \dot{m}_j

        **superheating or subcooling** *Works with pure fluids only!*

        .. math::
            val = T_{j} - td_{bp} - T_{bp}\left( p_{j}, fluid_{j} \right)

            \text{td: temperature difference, bp: boiling point}

        **vapour mass fraction** *Works with pure fluids only!*

        .. math::
            val = h_{j} - h \left( p_{j}, x_{j}, fluid_{j} \right)

        **Referenced values**

        **mass flow, pressure and enthalpy**

        .. math::
            val = x_{j} - x_{j,ref} \cdot a + b

        **temperatures**

        .. math::
            val = T \left( p_{j}, h_{j}, fluid_{j} \right) -
            T \left( p_{j}, h_{j}, fluid_{j} \right) \cdot a + b

        **Derivatives**

        **mass flow, pressure and enthalpy**

        .. math::

            J\left(\frac{\partial f_{i}}{\partial m_{j}}\right) = 1\\
            \text{for equation i, connection j}\\
            \text{pressure and enthalpy analogously}

        **temperatures**

        .. math::

            J\left(\frac{\partial f_{i}}{\partial p_{j}}\right) =
            -\frac{\partial T_{j}}{\partial p_{j}}\\
            J(\left(\frac{\partial f_{i}}{\partial h_{j}}\right) =
            -\frac{\partial T_{j}}{\partial h_{j}}\\
            J\left(\frac{\partial f_{i}}{\partial fluid_{j,k}}\right) =
            - \frac{\partial T_{j}}{\partial fluid_{j,k}}

            \forall k \in \text{fluid components}\\
            \text{for equation i, connection j}

        **volumetric flow**

        .. math::

            J\left(\frac{\partial f_{i}}{\partial m_{j}}\right) =
            -v \left( p_{j}, h_{j} \right)\\
            J\left(\frac{\partial f_{i}}{\partial p_{j}}\right) =
            -\frac{\partial v_{j}}{\partial p_{j}} \cdot \dot{m}_j\\
            J(\left(\frac{\partial f_{i}}{\partial h_{j}}\right) =
            -\frac{\partial v_{j}}{\partial h_{j}} \cdot \dot{m}_j\\

            \forall k \in \text{fluid components}\\
            \text{for equation i, connection j}

        **superheating or subcooling** *Works with pure fluids only!*

        .. math::

            J\left(\frac{\partial f_{i}}{\partial p_{j}}\right) =
            \frac{\partial T \left( p_{j}, h_{j}, fluid_{j} \right)}
            {\partial p_{j}} -
            \frac{\partial T_{bp} \left( p_{j}, fluid_{j} \right)}
            {\partial p_{j}} \\
            J\left(\frac{\partial f_{i}}{\partial h_{j}}\right) =
            \frac{\partial T \left( p_{j}, h_{j}, fluid_{j} \right)}
            {\partial h_{j}}\\

            \text{for equation i, connection j}\\
            \text{td: temperature difference, bp: boiling point}

        **vapour mass fraction** *Works with pure fluids only!*

        .. math::

            J\left(\frac{\partial f_{i}}{\partial p_{j}}\right) =
            -\frac{\partial h \left( p_{j}, x_{j}, fluid_{j} \right)}
            {\partial p_{j}}\\
            J\left(\frac{\partial f_{i}}{\partial h_{j}}\right) = 1\\
            \text{for equation i, connection j, x: vapour mass fraction}

        **Referenced values**

        **mass flow, pressure and enthalpy**

        .. math::
            J\left(\frac{\partial f_{i}}{\partial m_{j}}\right) = 1\\
            J\left(\frac{\partial f_{i}}{\partial m_{j,ref}}\right) = - a\\
            \text{for equation i, connection j}\\
            \text{pressure and enthalpy analogously}

        **temperatures**

        .. math::
            J\left(\frac{\partial f_{i}}{\partial p_{j}}\right) =
            \frac{dT_{j}}{dp_{j}}\\
            J\left(\frac{\partial f_{i}}{\partial h_{j}}\right) =
            \frac{dT_{j}}{dh_{j}}\\
            J\left(\frac{\partial f_{i}}{\partial fluid_{j,k}}\right) =
            \frac{dT_{j}}{dfluid_{j,k}}
            \; , \forall k \in \text{fluid components}\\
            J\left(\frac{\partial f_{i}}{\partial p_{j,ref}}\right) =
            \frac{dT_{j,ref}}{dp_{j,ref}} \cdot a \\
            J\left(\frac{\partial f_{i}}{\partial h_{j,ref}}\right) =
            \frac{dT_{j,ref}}{dh_{j,ref}} \cdot a \\
            J\left(\frac{\partial f_{i}}{\partial fluid_{j,k,ref}}\right) =
            \frac{dT_{j}}{dfluid_{j,k,ref}} \cdot a
            \; , \forall k \in \text{fluid components}\\
            \text{for equation i, connection j}
        """
        primary_vars = {'m': 0, 'p': 1, 'h': 2}
        k = self.num_comp_eq
        for c in self.conns.index:
            flow = c.to_flow()
            col = c.conn_loc * self.num_conn_vars

            # referenced mass flow, pressure or enthalpy
            for var, pos in primary_vars.items():
                if c.get_attr(var).ref_set is True:
                    ref = c.get_attr(var).ref
                    ref_col = ref.obj.conn_loc * self.num_conn_vars
                    self.vec_res[k] = (
                        c.get_attr(var).val_SI - (
                            ref.obj.get_attr(var).val_SI * ref.f + ref.d))
                    self.mat_deriv[k, col + pos] = 1
                    self.mat_deriv[k, ref_col + pos] = -c.get_attr(var).ref.f
                    k += 1

            # temperature
            if c.T.val_set is True:
                # if (np.absolute(self.vec_res[k]) > err ** 2 or
                #         self.iter % 2 == 0):
                self.vec_res[k] = c.T.val_SI - fp.T_mix_ph(
                    flow, T0=c.T.val_SI)

                # if not self.vec_z_filter[col + 1]:
                self.mat_deriv[k, col + 1] = (
                    -fp.dT_mix_dph(flow, T0=c.T.val_SI))
                # if not self.vec_z_filter[col + 2]:
                self.mat_deriv[k, col + 2] = (
                    -fp.dT_mix_pdh(flow, T0=c.T.val_SI))
                if len(self.fluids) != 1:
                    col_s = c.conn_loc * self.num_conn_vars + 3
                    col_e = (c.conn_loc + 1) * self.num_conn_vars
                    if not all(self.vec_z_filter[col_s:col_e]):
                        self.mat_deriv[k, col_s:col_e] = -fp.dT_mix_ph_dfluid(
                            flow, T0=c.T.val_SI)
                k += 1

            # referenced temperature
            if c.T.ref_set is True:
                ref = c.T.ref
                flow_ref = ref.obj.to_flow()
                ref_col = ref.obj.conn_loc * self.num_conn_vars
                # if (np.absolute(self.vec_res[k]) > err ** 2 or
                #         self.iter % 2 == 0):
                self.vec_res[k] = fp.T_mix_ph(flow, T0=c.T.val_SI) - (
                    fp.T_mix_ph(flow_ref, T0=ref.obj.T.val_SI) *
                    ref.f + ref.d)

                # if not self.vec_z_filter[col + 1]:
                self.mat_deriv[k, col + 1] = (
                    fp.dT_mix_dph(flow, T0=c.T.val_SI))
                # if not self.vec_z_filter[col + 1]:
                self.mat_deriv[k, col + 2] = (
                    fp.dT_mix_pdh(flow, T0=c.T.val_SI))

                # if not self.vec_z_filter[ref_col + 1]:
                self.mat_deriv[k, ref_col + 1] = -(
                    fp.dT_mix_dph(flow_ref, T0=ref.obj.T.val_SI) * ref.f)
                # if not self.vec_z_filter[ref_col + 2]:
                self.mat_deriv[k, ref_col + 2] = -(
                    fp.dT_mix_pdh(flow_ref, T0=ref.obj.T.val_SI) * ref.f)

                # dT / dFluid
                if len(self.fluids) != 1:
                    col_s = c.conn_loc * self.num_conn_vars + 3
                    col_e = (c.conn_loc + 1) * self.num_conn_vars
                    ref_col_s = ref.obj.conn_loc * self.num_conn_vars + 3
                    ref_col_e = (ref.obj.conn_loc + 1) * self.num_conn_vars
                    if not all(self.vec_z_filter[col_s:col_e]):
                        self.mat_deriv[k, col_s:col_e] = (
                            fp.dT_mix_ph_dfluid(flow, T0=c.T.val_SI))
                    if not all(self.vec_z_filter[ref_col_s:ref_col_e]):
                        self.mat_deriv[k, ref_col_s:ref_col_e] = -(
                            fp.dT_mix_ph_dfluid(flow_ref, T0=ref.obj.T.val_SI))
                k += 1

            # saturated steam fraction
            if c.x.val_set is True:
                if (np.absolute(self.vec_res[k]) > err ** 2 or
                        self.iter % 2 == 0):
                    self.vec_res[k] = c.h.val_SI - (
                        fp.h_mix_pQ(flow, c.x.val_SI))
                if not self.vec_z_filter[col + 1]:
                    self.mat_deriv[k, col + 1] = -(
                        fp.dh_mix_dpQ(flow, c.x.val_SI))
                self.mat_deriv[k, col + 2] = 1
                k += 1

            # volumetric flow
            if c.v.val_set is True:
                if (np.absolute(self.vec_res[k]) > err ** 2 or
                        self.iter % 2 == 0):
                    self.vec_res[k] = (
                        c.v.val_SI - fp.v_mix_ph(flow, T0=c.T.val_SI) *
                        c.m.val_SI)
                self.mat_deriv[k, col] = -fp.v_mix_ph(flow, T0=c.T.val_SI)
                self.mat_deriv[k, col + 1] = -(
                    fp.dv_mix_dph(flow, T0=c.T.val_SI) * c.m.val_SI)
                self.mat_deriv[k, col + 2] = -(
                    fp.dv_mix_pdh(flow, T0=c.T.val_SI) * c.m.val_SI)
                k += 1

            # temperature difference to boiling point
            if c.Td_bp.val_set is True:
                if (np.absolute(self.vec_res[k]) > err ** 2 or
                        self.iter % 2 == 0):
                    self.vec_res[k] = (
                        fp.T_mix_ph(flow, T0=c.T.val_SI) - c.Td_bp.val_SI -
                        fp.T_bp_p(flow))
                if not self.vec_z_filter[col + 1]:
                    self.mat_deriv[k, col + 1] = (
                        fp.dT_mix_dph(flow, T0=c.T.val_SI) - fp.dT_bp_dp(flow))
                if not self.vec_z_filter[col + 2]:
                    self.mat_deriv[k, col + 2] = fp.dT_mix_pdh(
                        flow, T0=c.T.val_SI)
                k += 1

            # fluid composition balance
            if c.fluid.balance:
                j = 0
                res = 1
                for f in self.fluids:
                    res -= c.fluid.val[f]
                    self.mat_deriv[k, c.conn_loc + 3 + j] = -1
                    j += 1

                self.vec_res[k] = res
                k += 1

        # equations and derivatives for specified primary variables are static
        if self.iter == 0:
            for c in self.conns.index:
                flow = c.to_flow()
                col = c.conn_loc * self.num_conn_vars

                # specified mass flow, pressure and enthalpy
                for var, pos in primary_vars.items():
                    if c.get_attr(var).val_set is True:
                        self.vec_res[k] = 0
                        self.mat_deriv[
                            k, col + pos] = 1
                        k += 1

                j = 0
                # specified fluid mass fraction
                for f in self.fluids:
                    if c.fluid.val_set[f]:
                        self.mat_deriv[k, col + 3 + j] = 1
                        k += 1
                    j += 1

    def solve_busses(self):
        r"""
        Calculate the equations and the partial derivatives for the busses.

        - Iterate through busses in network to get residuals and derivatives.
        - Place residual values in residual value vector of the network.
        - Place partial derivatives in jacobian matrix of the network.
        """
        row = self.num_comp_eq + self.num_conn_eq
        for b in self.busses.values():
            if b.P.is_set is True:
                P_res = 0
                for cp in b.comps.index:
                    i = self.comps.loc[cp].i.tolist()
                    o = self.comps.loc[cp].o.tolist()

                    bus = b.comps.loc[cp]

                    P_res += cp.bus_func(bus)
                    deriv = -cp.bus_deriv(bus)

                    j = 0
                    for c in i + o:
                        loc = self.conns.index.get_loc(c)
                        # start collumn index
                        coll_s = loc * self.num_conn_vars
                        # end collumn index
                        coll_e = (loc + 1) * self.num_conn_vars
                        self.mat_deriv[row, coll_s:coll_e] = deriv[:, j]
                        j += 1

                self.vec_res[row] = b.P.val - P_res

                row += 1

    def postprocessing(self):
        r"""Calculate bus, component parameters and connection parameters."""
        # components
        self.comps.apply(network.process_components, axis=1)

        # busses
        for b in self.busses.values():
            b.P.val = 0
            for cp in b.comps.index:
                # get components bus func value
                val = cp.bus_func(b.comps.loc[cp])
                # save as reference value
                if self.mode == 'design':
                    b.comps.loc[cp].P_ref = (
                        cp.bus_func(b.comps.loc[cp]) /
                        abs(b.comps.loc[cp].char.evaluate(1)))
                b.P.val += val

        # connections
        for c in self.conns.index:
            c.good_starting_values = True
            c.T.val_SI = fp.T_mix_ph(c.to_flow(), T0=c.T.val_SI)
            c.v.val_SI = fp.v_mix_ph(c.to_flow(), T0=c.T.val_SI) * c.m.val_SI
            c.T.val = (c.T.val_SI / self.T[c.T.unit][1] - self.T[c.T.unit][0])
            c.m.val = c.m.val_SI / self.m[c.m.unit]
            c.p.val = c.p.val_SI / self.p[c.p.unit]
            c.h.val = c.h.val_SI / self.h[c.h.unit]
            c.v.val = c.v.val_SI / self.v[c.v.unit]
            fluid = hlp.single_fluid(c.fluid.val)
            if isinstance(fluid, str) and not c.x.val_set:
                c.x.val_SI = fp.Q_ph(c.p.val_SI, c.h.val_SI, fluid)
                c.x.val = c.x.val_SI
            c.T.val0 = c.T.val
            c.m.val0 = c.m.val
            c.p.val0 = c.p.val
            c.h.val0 = c.h.val
            c.fluid.val0 = c.fluid.val.copy()

        msg = 'Postprocessing complete.'
        logging.info(msg)

    def process_components(cols):
        cols.name.calc_parameters()

# %% printing and plotting

    def print_results(self):
        r"""Print the calculations results to prompt."""
        cp_sort = self.comps.copy()
        # sort components by component type alphabetically
        cp_sort['cp'] = cp_sort.apply(network.get_class_base, axis=1)
        cp_sort['label'] = cp_sort.apply(network.get_props, axis=1,
                                         args=('label',))
        cp_sort.drop('i', axis=1, inplace=True)
        cp_sort.drop('o', axis=1, inplace=True)

        pd.options.mode.chained_assignment = None
        for c in cp_sort.cp.unique():
            df = cp_sort[cp_sort['cp'] == c]

            # gather parameters to print for components of type c
            cols = []
            for col, val in df.index[0].attr().items():
                if isinstance(val, dc.dc_cp):
                    if val.get_attr('printout'):
                        cols += [col]

            # are there any parameters to print?
            if len(cols) > 0:
                for col in cols:
                    df[col] = df.apply(network.print_components, axis=1,
                                       args=(col,))

                df.set_index('label', inplace=True)
                df.drop('cp', axis=1, inplace=True)
                df = df.dropna(how='all')

                if len(df) > 0:
                    # printout with tabulate
                    print('##### RESULTS (' + c + ') #####')
                    print(tabulate(df, headers='keys', tablefmt='psql',
                                   floatfmt='.2e'))

        # connection properties
        df = pd.DataFrame(columns=['m / (' + self.m_unit + ')',
                                   'p / (' + self.p_unit + ')',
                                   'h / (' + self.h_unit + ')',
                                   'T / (' + self.T_unit + ')'])
        for c in self.conns.index:
            if c.printout is True:
                row = (c.s.label + ':' + c.s_id + ' -> ' +
                       c.t.label + ':' + c.t_id)
                df.loc[row] = (
                        [c.m.val_SI / self.m[self.m_unit],
                         c.p.val_SI / self.p[self.p_unit],
                         c.h.val_SI / self.h[self.h_unit],
                         c.T.val_SI / self.T[self.T_unit][1] -
                         self.T[self.T_unit][0]]
                        )
        if len(df) > 0:
            print('##### RESULTS (connections) #####')
            print(tabulate(df, headers='keys', tablefmt='psql',
                           floatfmt='.3e'))

        for b in self.busses.values():
            df = pd.DataFrame(columns=['component', 'value'])
            if b.printout is True:
                df['cp'] = b.comps.index
                df['ref'] = b.comps['P_ref'].values
                df['component'] = df['cp'].apply(lambda x: x.label)
                df['value'] = df['cp'].apply(
                    lambda x: x.bus_func(b.comps.loc[x]))
                df.loc['total'] = df.sum()
                df.loc['total', 'component'] = 'total'
                df.set_index('component', inplace=True)
                df.drop('cp', axis=1, inplace=True)
                print('##### RESULTS (' + b.label + ') #####')
                print(tabulate(df, headers='keys', tablefmt='psql',
                               floatfmt='.3e'))

    def print_components(c, *args):
        if c.name.printout is True:
            return c.name.get_attr(args[0]).val
        else:
            return np.nan

# %% saving

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

        - netw.csv (network information)
        - conn.csv (connection information)
        - folder comps containing .csv files for busses and characteristics
          as well as .csv files for all types of components within your
          network.
        """
        if path[-1] != '/' and path[-1] != '\\':
            path += '/'
        path = hlp.modify_path_os(path)

        logging.debug('Saving network to path ' + path + '.')
        # creat path, if non existent
        if not os.path.exists(path):
            os.makedirs(path)

        # create path for component folder if non existent
        path_comps = hlp.modify_path_os(path + 'comps/')
        if not os.path.exists(path_comps):
            os.makedirs(path_comps)

        # save all network information
        self.save_network(path + 'netw.csv')
        self.save_connections(path + 'conn.csv')
        self.save_components(path_comps)
        self.save_busses(path_comps + 'bus.csv')
        self.save_characteristics(path_comps)

    def save_network(self, fn):
        r"""
        Save basic network configuration.

        Parameters
        ----------
        fn : str
            Path/filename for the network configuration file.
        """
        data = {}
        data['m_unit'] = self.m_unit
        data['p_unit'] = self.p_unit
        data['p_min'] = self.p_range[0]
        data['p_max'] = self.p_range[1]
        data['h_unit'] = self.h_unit
        data['h_min'] = self.h_range[0]
        data['h_max'] = self.h_range[1]
        data['T_unit'] = self.T_unit
        data['T_min'] = self.T_range[0]
        data['T_max'] = self.T_range[1]
        data['fluids'] = [self.fluids]

        df = pd.DataFrame(data=data)

        df.to_csv(fn, sep=';', decimal='.', index=False, na_rep='nan')
        logging.debug('Network information saved to ' + fn + '.')

    def save_connections(self, fn):
        r"""
        Save the connection properties.

        - Uses connections object id as row identifier and saves
            * connections source and target as well as
            * properties with references and
            * fluid vector (including user specification if structure is True).
        - Connections source and target are identified by its labels.

        Parameters
        ----------
        fn : str
            Path/filename for the file.
        """
        f = network.get_props
        df = pd.DataFrame()
        # connection id
        df['id'] = self.conns.apply(network.get_id, axis=1)

        # general connection parameters
        # source
        df['s'] = self.conns.apply(f, axis=1, args=('s', 'label'))
        df['s_id'] = self.conns.apply(f, axis=1, args=('s_id',))
        # target
        df['t'] = self.conns.apply(f, axis=1, args=('t', 'label'))
        df['t_id'] = self.conns.apply(f, axis=1, args=('t_id',))

        # design and offdesign properties
        cols = ['design', 'offdesign', 'design_path', 'local_design',
                'local_offdesign']
        for key in cols:
            df[key] = self.conns.apply(f, axis=1, args=(key,))

        # fluid properties
        cols = ['m', 'p', 'h', 'T', 'x', 'v', 'Td_bp']
        for key in cols:
            # values and units
            df[key] = self.conns.apply(f, axis=1, args=(key, 'val'))
            df[key + '_unit'] = self.conns.apply(f, axis=1, args=(key, 'unit'))

            # connection parametrisation
            df[key + '_unit_set'] = self.conns.apply(f, axis=1,
                                                     args=(key, 'unit_set'))
            df[key + '0'] = self.conns.apply(f, axis=1, args=(key, 'val0'))
            df[key + '_set'] = self.conns.apply(f, axis=1,
                                                args=(key, 'val_set'))
            df[key + '_ref'] = self.conns.apply(f, axis=1,
                                                args=(key, 'ref', 'obj',)
                                                ).astype(str)
            df[key + '_ref'] = df[key + '_ref'].str.extract(r' at (.*?)>',
                                                            expand=False)
            df[key + '_ref_f'] = self.conns.apply(f, axis=1,
                                                  args=(key, 'ref', 'f',))
            df[key + '_ref_d'] = self.conns.apply(f, axis=1,
                                                  args=(key, 'ref', 'd',))
            df[key + '_ref_set'] = self.conns.apply(f, axis=1,
                                                    args=(key, 'ref_set',))

        # state property
        key = 'state'
        df[key] = self.conns.apply(f, axis=1, args=(key, 'val'))
        df[key + '_set'] = self.conns.apply(f, axis=1, args=(key, 'is_set'))

        # fluid composition
        for val in self.fluids:
            # fluid mass fraction
            df[val] = self.conns.apply(f, axis=1, args=('fluid', 'val', val))

            # fluid mass fraction parametrisation
            df[val + '0'] = self.conns.apply(f, axis=1,
                                             args=('fluid', 'val0', val))
            df[val + '_set'] = self.conns.apply(f, axis=1,
                                                args=('fluid', 'val_set', val))

        # fluid balance
        df['balance'] = self.conns.apply(f, axis=1, args=('fluid', 'balance'))

        df.to_csv(fn, sep=';', decimal='.', index=False, na_rep='nan')
        logging.debug('Connection information saved to ' + fn + '.')

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
        busses = self.busses.values()
        # create / overwrite csv file
        cp_sort = self.comps.copy()
        # component type
        cp_sort['cp'] = cp_sort.apply(network.get_class_base, axis=1)

        # busses
        cp_sort['busses'] = cp_sort.apply(network.get_busses, axis=1,
                                          args=(busses,))

        for var in ['param', 'P_ref', 'char']:
            cp_sort['bus_' + var] = cp_sort.apply(network.get_bus_data, axis=1,
                                                  args=(busses, var))

        pd.options.mode.chained_assignment = None
        f = network.get_props
        for c in cp_sort.cp.unique():
            df = cp_sort[cp_sort['cp'] == c]

            # basic information
            cols = ['label', 'design', 'offdesign', 'design_path',
                    'local_design', 'local_offdesign']
            for col in cols:
                df[col] = df.apply(f, axis=1, args=(col,))

            # attributes
            for col, data in df.index[0].attr().items():
                # component characteristics container
                if isinstance(data, dc.dc_cc):
                    df[col] = df.apply(f, axis=1,
                                       args=(col, 'func')).astype(str)
                    df[col] = df[col].str.extract(r' at (.*?)>', expand=False)
                    df[col + '_set'] = df.apply(f, axis=1,
                                                args=(col, 'is_set'))
                    df[col + '_param'] = df.apply(f, axis=1,
                                                  args=(col, 'param'))

                # component characteristic map container
                elif isinstance(data, dc.dc_cm):
                    df[col] = df.apply(f, axis=1,
                                       args=(col, 'func')).astype(str)
                    df[col] = df[col].str.extract(r' at (.*?)>',
                                                  expand=False)
                    df[col + '_set'] = df.apply(f, axis=1,
                                                args=(col, 'is_set'))
                    df[col + '_param'] = df.apply(f, axis=1,
                                                  args=(col, 'param'))

                # component property container
                elif isinstance(data, dc.dc_cp):
                    df[col] = df.apply(f, axis=1, args=(col, 'val'))
                    df[col + '_set'] = df.apply(f, axis=1,
                                                args=(col, 'is_set'))
                    df[col + '_var'] = df.apply(f, axis=1,
                                                args=(col, 'is_var'))

                # component property container
                elif isinstance(data, dc.dc_simple):
                    df[col] = df.apply(f, axis=1, args=(col, 'val'))
                    df[col + '_set'] = df.apply(f, axis=1,
                                                args=(col, 'is_set'))

                # component property container
                elif isinstance(data, dc.dc_gcp):
                    df[col] = df.apply(f, axis=1, args=(col, 'method'))

            df.set_index('label', inplace=True)
            df.drop('i', axis=1, inplace=True)
            df.drop('o', axis=1, inplace=True)
            fn = path + c + '.csv'
            df.to_csv(fn, sep=';', decimal='.', index=True, na_rep='nan')
            logging.debug('Component information (' + c + ') saved to ' +
                          fn + '.')

    def save_busses(self, fn):
        r"""
        Save the bus properties.

        Parameters
        ----------
        fn : str
            Path/filename for the file.
        """
        if len(self.busses) > 0:
            df = pd.DataFrame({'id': self.busses.values()},
                              index=self.busses.values())
            df['label'] = df.apply(network.get_props, axis=1, args=('label',))
            df['P'] = df.apply(network.get_props, axis=1, args=('P', 'val'))
            df['P_set'] = df.apply(network.get_props, axis=1,
                                   args=('P', 'is_set'))
            df.drop('id', axis=1, inplace=True)

            df.set_index('label', inplace=True)
            df.to_csv(fn, sep=';', decimal='.', index=True, na_rep='nan')
            logging.debug('Bus information saved to ' + fn + '.')

    def save_characteristics(self, path):
        r"""
        Save the characteristics.

        Parameters
        ----------
        fn : str
            Path/filename for the file.
        """
        # components
        cp_sort = self.comps
        cp_sort['cp'] = cp_sort.apply(network.get_class_base, axis=1)

        # characteristic lines in components
        chars = []
        for c in cp_sort.cp.unique():
            df = cp_sort[cp_sort['cp'] == c]

            for col, data in df.index[0].attr().items():
                if isinstance(data, dc.dc_cc):
                    chars += df.apply(network.get_props, axis=1,
                                      args=(col, 'func')).tolist()

        # characteristic lines in busses
        for bus in self.busses.values():
            for c in bus.comps.index:
                ch = bus.comps.loc[c].char
                if ch not in chars:
                    chars += [ch]

        if len(chars) > 0:
            # get id and data
            df = pd.DataFrame({'id': chars}, index=chars)
            df['id'] = df.apply(network.get_id, axis=1)
            df['type'] = df.apply(network.get_class_base, axis=1)

            cols = ['x', 'y', 'extrapolate']
            for val in cols:
                df[val] = df.apply(network.get_props, axis=1, args=(val,))

            # write to char.csv
            fn = path + 'char_line.csv'
            df.to_csv(fn, sep=';', decimal='.', index=False, na_rep='nan')
            logging.debug('Characteristic line information saved to ' + fn +
                          '.')

        # characteristic maps in components
        chars = []
        for c in cp_sort.cp.unique():
            df = cp_sort[cp_sort['cp'] == c]

            for col, data in df.index[0].attr().items():
                if isinstance(data, dc.dc_cm):
                    chars += df.apply(network.get_props, axis=1,
                                      args=(col, 'func')).tolist()

        if len(chars) > 0:
            # get id and data
            df = pd.DataFrame({'id': chars}, index=chars)
            df['id'] = df.apply(network.get_id, axis=1)
            df['type'] = df.apply(network.get_class_base, axis=1)

            cols = ['x', 'y', 'z1', 'z2']
            for val in cols:
                df[val] = df.apply(network.get_props, axis=1, args=(val,))

            # write to char_map.csv
            fn = path + 'char_map.csv'
            df.to_csv(fn, sep=';', decimal='.', index=False, na_rep='nan')
            logging.debug('Characteristic map information saved to ' + fn +
                          '.')

    def get_id(c):
        """Return the id of the python object."""
        return str(c.name)[str(c.name).find(' at ') + 4:-1]

    def get_class_base(c):
        """Return the class name."""
        return c.name.__class__.__name__

    def get_props(c, *args):
        """Return properties."""
        if hasattr(c.name, args[0]):
            if (not isinstance(c.name.get_attr(args[0]), int) and
                    not isinstance(c.name.get_attr(args[0]), str) and
                    not isinstance(c.name.get_attr(args[0]), float) and
                    not isinstance(c.name.get_attr(args[0]), list) and
                    not isinstance(c.name.get_attr(args[0]), np.ndarray) and
                    not isinstance(c.name.get_attr(args[0]), con.connection)):
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
        else:
            return ''

    def get_busses(c, *args):
        """Return the list of busses a component is integrated in."""
        busses = []
        for bus in args[0]:
            if c.name in bus.comps.index:
                busses += [bus.label]
        return busses

    def get_bus_data(c, *args):
        """Return bus information of a component."""
        items = []
        if args[1] == 'char':
            for bus in args[0]:
                if c.name in bus.comps.index:
                    val = bus.comps.loc[c.name][args[1]]
                    items += [str(val)[str(val).find(' at ') + 4:-1]]

        else:
            for bus in args[0]:
                if c.name in bus.comps.index:
                    items += [bus.comps.loc[c.name][args[1]]]

        return items
