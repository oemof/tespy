"""
.. module:: networks
    :synopsis: contains logic of the network

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>
"""

import math

import pandas as pd
from tabulate import tabulate

import numpy as np
from numpy.linalg import inv
from numpy.linalg import norm

from tespy.components import components as cmp
from tespy.components import characteristics as char
from tespy import connections as con
from tespy import helpers as hlp

import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as colors

import collections

import time
import os
from CoolProp.CoolProp import PropsSI as CPPSI


class network:
    r"""
    Class component is the base class of all TESPy components.

    Parameters
    ----------
    fluids : list
        A list of all fluids within the network container.

    m_unit : String
        Specify the unit for mass flow: 'kg / s', 't / h'.

    v_unit : String
        Specify the unit for volumetric flow: 'm3 / s', 'm3 / h', 'l / s', 'l / h'.

    p_unit : String
        Specify the unit for pressure: 'Pa', 'psi', 'bar', 'MPa'.

    h_unit : String
        Specify the unit for mass flow: 'J / kg', 'kJ / kg', 'MJ / kg'.

    T_unit : String
        Specify the unit for mass flow: 'K', 'C', 'F'.

    p_range : list
        List with minimum and maximum values for pressure value range.

    h_range : list
        List with minimum and maximum values for enthalpy value range.

    T_range : list
        List with minimum and maximum values for temperature value range.

    Note
    ----
    Unit specification is optional: If not specified the SI unit (first element in above lists) will be applied!

    Range specification is optional, too. The value range is used to stabilise the newton algorith.
    For more information see the "getting started" section in the online-documentation.

    Printoptions can be specified with the :func:`tespy.networks.network.set_printoptions`-method, see example.

    Example
    -------
    Basic example for a setting up a tespy.networks.network object. Specifying
    the fluids is mandatory! Unit systems, fluid property range and printlevel
    are optional.

    Standard value for printoptions print_level is info. You can modify this with the
    :func:`tespy.networks.network.set_printoptions`-method by specifying a print_level, or specifying the printout manually.

    >>> from tespy import nwk
    >>> fluid_list = ['water', 'air', 'R134a']
    >>> mynetwork = nwk.network(fluids=fluid_list, p_unit='bar', T_unit='C')
    >>> mynetwork.set_attr(p_range=[1, 10])
    >>> type(mynetwork)
    <class 'tespy.networks.network'>
    >>> mynetwork.set_printoptions(print_level='warn')
    >>> (mynetwork.comperr, mynetwork.iterinfo, mynetwork.nwkwarn)
    (True, False, True)
    >>> mynetwork.set_printoptions(print_level='err')
    >>> (mynetwork.comperr, mynetwork.iterinfo, mynetwork.nwkwarn)
    (True, False, False)
    >>> mynetwork.set_printoptions(print_level='none')
    >>> mynetwork.set_printoptions(iterinfo=True, compwarn=True)
    >>> (mynetwork.comperr, mynetwork.iterinfo, mynetwork.nwkwarn, mynetwork.compwarn)
    (False, True, False, True)
    """

    def __init__(self, fluids, **kwargs):

        # initialisation of basic properties
        self.checked = False
        # connection dataframe
        self.conns = pd.DataFrame(columns=['s', 's_id', 't', 't_id'])
        # list for busses
        self.busses = []

        # fluid list and constants
        self.fluids = sorted(fluids)
        for f in self.fluids:
            try:
                hlp.molar_masses[f] = CPPSI('M', f)
            except:
                hlp.molar_masses[f] = 1

            try:
                hlp.gas_constants[f] = CPPSI('GAS_CONSTANT', f)
            except:
                hlp.gas_constants[f] = np.nan

        # initialise fluid property memorisation function for this network
        hlp.memorise.add_fluids(self.fluids)

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

        # processing printoptions
        self.print_level = 'info'
        self.set_printoptions()

        # standard unit set
        self.m_unit = self.SI_units['m']
        self.p_unit = self.SI_units['p']
        self.h_unit = self.SI_units['h']
        self.T_unit = self.SI_units['T']
        self.v_unit = self.SI_units['v']

        # standard value range
        self.p_range_SI = np.array([2e2, 300e5])
        self.h_range_SI = np.array([1e3, 7e6])
        self.T_range_SI = np.array([273.16, 1773.15])

        # add attributes from kwargs
        for key in kwargs:
            if key in self.attr():
                self.__dict__.update({key: kwargs[key]})

        self.set_attr(**kwargs)

    def set_attr(self, **kwargs):
        r"""
        Sets, resets or unsets attributes of a network for provided keyword arguments.

        Parameters
        ----------
        m_unit : String
            Specify the unit for mass flow: 'kg / s', 't / h'.

        v_unit : String
            Specify the unit for volumetric flow: 'm3 / s', 'm3 / h', 'l / s', 'l / h'.

        p_unit : String
            Specify the unit for pressure: 'Pa', 'psi', 'bar', 'MPa'.

        h_unit : String
            Specify the unit for mass flow: 'J / kg', 'kJ / kg', 'MJ / kg'.

        T_unit : String
            Specify the unit for mass flow: 'K', 'C', 'F'.

        p_range : list
            List with minimum and maximum values for pressure value range.

        h_range : list
            List with minimum and maximum values for enthalpy value range.

        T_range : list
            List with minimum and maximum values for temperature value range.

        Note
        ----
        Use the :func:`tespy.networks.network.set_printoptions` method for adjusting printouts.
        """
        # add attributes from kwargs
        for key in kwargs:
            if key in self.attr():
                self.__dict__.update({key: kwargs[key]})

        # unit sets
        if self.m_unit not in self.m.keys():
            msg = ('Allowed units for mass flow are: ' +
                   str(self.m.keys()))
            raise hlp.MyNetworkError(msg)

        if self.p_unit not in self.p.keys():
            msg = ('Allowed units for pressure are: ' +
                   str(self.p.keys()))
            raise hlp.MyNetworkError(msg)

        if self.h_unit not in self.h.keys():
            msg = ('Allowed units for enthalpy are: ' +
                   str(self.h.keys()))
            raise hlp.MyNetworkError(msg)

        if self.T_unit not in self.T.keys():
            msg = ('Allowed units for temperature are: ' +
                   str(self.T.keys()))
            raise hlp.MyNetworkError(msg)

        if self.v_unit not in self.v.keys():
            msg = ('Allowed units for volumetric flow are: ' +
                   str(self.v.keys()))
            raise hlp.MyNetworkError(msg)

        # value ranges
        if 'p_range' in kwargs.keys():
            if not isinstance(self.p_range, list):
                msg = ('Specify the value range as list: [p_min, p_max]')
                raise TypeError(msg)
            else:
                self.p_range_SI = np.array(self.p_range) * self.p[self.p_unit]
        else:
            self.p_range = self.p_range_SI / self.p[self.p_unit]

        if 'h_range' in kwargs.keys():
            if not isinstance(self.h_range, list):
                msg = ('Specify the value range as list: [h_min, h_max]')
                raise TypeError(msg)
            else:
                self.h_range_SI = np.array(self.h_range) * self.h[self.h_unit]
        else:
            self.h_range = self.h_range_SI / self.h[self.h_unit]

        if 'T_range' in kwargs.keys():
            if not isinstance(self.T_range, list):
                msg = ('Specify the value range as list: [T_min, T_max]')
                raise TypeError(msg)
            else:
                self.T_range_SI = ((np.array(self.T_range) +
                                    self.T[self.T_unit][0]) *
                                   self.T[self.T_unit][1])
        else:
            self.T_range = (self.T_range_SI / self.T[self.T_unit][1] -
                            self.T[self.T_unit][0])

        for f in self.fluids:
            if 'TESPy::' in f:
                hlp.memorise.vrange[f][0] = self.p_range_SI[0]
                hlp.memorise.vrange[f][1] = self.p_range_SI[1]
                hlp.memorise.vrange[f][2] = self.T_range_SI[0]
                hlp.memorise.vrange[f][3] = self.T_range_SI[1]

            if 'INCOMP::' in f:
                hlp.memorise.vrange[f][0] = self.p_range_SI[0]
                hlp.memorise.vrange[f][1] = self.p_range_SI[1]

    def get_attr(self, key):
        r"""
        Get the value of a networks attribute.

        Parameters
        ----------
        key : String
            The attribute you want to retrieve.

        Returns
        -------
        out :
            Specified attribute.
        """
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            if self.nwkwarn:
                print('No attribute \"' + str(key) + '\" available!')
            return None

    def attr(self):
        return ['m_unit', 'p_unit', 'h_unit', 'T_unit', 'v_unit',
                'p_range', 'h_range', 'T_range']

    def set_printoptions(self, **kwargs):
        r"""
        Specification of printouts for tespy.networks.network object.

        Parameters
        ----------
        print_level : String
            Select the print level:

            - 'info': all printouts.
            - 'warn': errors and warnings
            - 'err': errors only
            - 'none': no printouts

        compinfo : bool
            Print information for components.

        compwarn : bool
            Print warnings for components.

        comperr : bool
            Print errors for components.

        nwkinfo : bool
            Print information for network.

        nwkwarn : bool
            Print warnings for network.

        nwkerr : bool
            Print errors for network.

        iterinfo : bool
            Printouts of iteration information in solving process.

        Note
        ----
        Specifying specific print leves manually overwrites the print_level information.
        See the example in :func:`tespy.networks.network`.
        """
        self.print_level = kwargs.get('print_level', self.print_level)

        if self.print_level == 'info':
            self.compinfo = True
            self.nwkinfo = True
            self.iterinfo = True
            self.compwarn = True
            self.nwkwarn = True
            self.comperr = True
            self.nwkerr = True

        elif self.print_level == 'warn':
            self.compinfo = False
            self.nwkinfo = False
            self.iterinfo = False
            self.compwarn = True
            self.nwkwarn = True
            self.comperr = True
            self.nwkerr = True

        elif self.print_level == 'err':
            self.compinfo = False
            self.nwkinfo = False
            self.iterinfo = False
            self.compwarn = False
            self.nwkwarn = False
            self.comperr = True
            self.nwkerr = True

        elif self.print_level == 'none':
            self.compinfo = False
            self.nwkinfo = False
            self.iterinfo = False
            self.compwarn = False
            self.nwkwarn = False
            self.comperr = False
            self.nwkerr = False
        else:
            msg = ('Available print leves are: \'info\', \'warn\', \'err\' and'
                   '\'none\'.')
            raise ValueError(msg)

        self.compinfo = kwargs.get('compinfo', self.compinfo)
        self.nwkinfo = kwargs.get('nwkinfo', self.nwkinfo)
        self.iterinfo = kwargs.get('iterinfo', self.iterinfo)
        self.compwarn = kwargs.get('compwarn', self.compwarn)
        self.nwkwarn = kwargs.get('nwkwarn', self.nwkwarn)
        self.comperr = kwargs.get('comperr', self.comperr)
        self.nwkerr = kwargs.get('nwkerr', self.nwkerr)

    def add_subsys(self, *args):
        r"""
        Adds one or more subsystem to the network.

        Parameters
        ----------
        c : tespy.components.subsystems.subsystem
            The subsystem to be added to the network, subsystem objects si :code:`network.add_subsys(s1, s2, s3, ...)`.
        """
        for subsys in args:
            for c in subsys.conns:
                self.add_conns(c)

    def add_conns(self, *args):
        r"""
        Adds one or more connections to the network.

        Parameters
        ----------
        c : tespy.connections.connection
            The connection to be added to the network, connections objects ci :code:`add_conns(c1, c2, c3, ...)`.
        """
        for c in args:
            if not isinstance(c, con.connection):
                raise TypeError('Must provide tespy.connections.connection objects'
                                ' as parameters.')

            self.conns.loc[c] = [c.s, c.s_id, c.t, c.t_id]
            # set status "checked" to false, if conneciton is added to network.
            self.checked = False

    def del_conns(self, c):
        r"""
        Removes one or more connections from the network.

        Parameters
        ----------
        c : tespy.connections.connection
            The connection to be removed from the network, connections objects ci :code:`del_conns(c1, c2, c3, ...)`.
        """
        self.conns.drop(self.conns.index(c))
        # set status "checked" to false, if conneciton is deleted from network.
        self.checked = False

    def check_conns(self):
        r"""
        Checks the networks connections for multiple usage of inlets or outlets of components.
        """
        dub = self.conns.loc[self.conns.duplicated(['s', 's_id']) == True].index
        for c in dub:
            msg = ('The source ' + str(c.s.label) + ' (' + str(c.s_id) + ') is '
                   'attached to more than one connection. Please check your network.')
            raise hlp.MyNetworkError(msg)

        dub = self.conns.loc[self.conns.duplicated(['t', 't_id']) == True].index
        for c in dub:
            msg = ('The target ' + str(c.t.label) + ' (' + str(c.t_id) + ') is '
                   'attached to more than one connection. Please check your network.')
            raise hlp.MyNetworkError(msg)

    def add_busses(self, *args):
        r"""
        Adds one or more busses to the network.

        Parameters
        ----------
        b : tespy.connections.bus
            The bus to be added to the network, bus objects bi :code:`add_busses(b1, b2, b3, ...)`.
        """
        for b in args:
            if self.check_busses(b):
                self.busses += [b]

    def del_busses(self, b):
        r"""
        Removes one or more busses from the network.

        Parameters
        ----------
        b : tespy.connections.bus
            The bus to be removed from the network, bus objects bi :code:`add_busses(b1, b2, b3, ...)`.
        """
        if b in self.busses:
            self.busses.remove(b)

    def check_busses(self, b):
        r"""
        Checks the busses to be added for type, duplicates and identical labels.

        Parameters
        ----------
        b : tespy.connections.bus
            The bus to be checked.
        """
        if isinstance(b, con.bus):
            if b not in self.busses:
                if b.label not in [x.label for x in self.busses]:
                    return True
                else:
                    msg = ('Network already has a bus with the name ' + b.label + '.')
                    raise hlp.MyNetworkError(msg)
            else:
                msg = 'Network contains this bus (' + str(b) + ') already.'
                raise hlp.MyNetworkError(msg)
        else:
            msg = 'Only objects of type bus are allowed in *args.'
            raise TypeError(msg)

        return False

    def check_network(self):
        r"""
        Checks the network for consistency, have all components the correct amount of incoming and outgoing connections?
        """
        self.check_conns()
        # get unique components in connections dataframe
        comps = pd.unique(self.conns[['s', 't']].values.ravel())
        self.init_components(comps)  # build the dataframe for components
        # count number of incoming and outgoing connections and compare to expected values
        for comp in self.comps.index:
            num_o = (self.conns[['s', 't']] == comp).sum().s
            num_i = (self.conns[['s', 't']] == comp).sum().t
            if num_o != comp.num_o:
                msg = (comp.label + ' is missing ' + str(comp.num_o - num_o) +
                       ' outgoing connections. Make sure all outlets are '
                       ' connected and all connections have been added to the '
                       'network.')
            elif num_i != comp.num_i:
                msg = (comp.label + ' is missing ' + str(comp.num_i - num_i) +
                       ' incoming connections. Make sure all inlets are '
                       ' connected and all connections have been added to the '
                       'network.')
            else:
                continue

            # raise an error in case network check is unsuccesful
            raise hlp.MyNetworkError(msg)

        # network checked
        self.checked = True
        msg = 'Networkcheck successfull.'
        if self.nwkinfo:
            print(msg)

    def initialise(self):
        r"""
        !!!!!!!!!!!!!!!!!!!!!
        Initilialises the network in the following steps:

        - Component initialisation for design case: Call comp_init method.
        - Component initialisation for offdesign case: Switch components to offdesign mode.
        - Fluid propagation on all connections.
        - Fluid property initialisation.
        - Initialisiation from .csv-files
        """
        self.errors = []
        if self.nwkinfo:
            msg = ('Have you adjusted the value ranges for pressure, enthalpy'
                   ' and temperature according to the specified unit system?')
            print(msg)

        if len(self.fluids) == 0:
            msg = ('Network has no fluids, please specify a list with fluids on network creation.')
            raise hlp.MyNetworkError(msg)

        if self.mode == 'offdesign':
            # characteristics for offdesign
            self.init_offdesign()
        else:
            # component initialisation for design case
            for cp in self.comps.index:
                cp.comp_init(self)

        self.init_fluids()  # start standard fluid initialisation
        self.init_properties()  # start standard property initialisation

        if self.mode == 'offdesign' and self.design_file is None:
            msg = ('Please provide \'design_file\' for every offdesign calculation.')
            raise hlp.MyNetworkError(msg)  # must provide design_file
        else:
            self.init_csv()  # initialisation from csv

        if self.nwkinfo:
            print('Network initialised.')

    def init_components(self, comps):
        r"""
        Sets up a dataframe for the network's components and checks, if all components have unique labels.

        Note
        ----
        The dataframe for the components is derived from the network's connections.
        Thus it holds no additional information, the dataframe is used to simplify the code.

        dataframe :code:`network.comps`:

        ======================== ============================ =======
         index                    i                            o
        ======================== ============================ =======
         type: component object   type: list                   see i
         value: object id         values: connection objects
        ======================== ============================ =======
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
            # save the incoming and outgoing as well as the number of connections as component attribute
            comp.inl = t.tolist()
            comp.outl = s.tolist()
            comp.num_i = len(comp.inlets())
            comp.num_o = len(comp.outlets())
            labels += [comp.label]

        # check for duplicates in the component labels
        if len(labels) != len(list(set(labels))):
            duplicates = [item for item, count in collections.Counter(labels).items() if count > 1]
            msg = ('All Components must have unique labels, duplicates are: ' + str(duplicates))
            raise hlp.MyNetworkError(msg)

    def init_fluids(self):
        r"""
        Initialises the fluid vector on every connection of the network.

        - Create fluid vector for every component as dict,
          index: nw.fluids,
          values: 0 if not set by user.
        - Create fluid_set vector with same logic,
          index: nw.fluids,
          values: False if not set by user.
        - If there are any combustion chambers in the network, calculate fluid vector starting from there.
        - Propagate fluid vector in direction of sources and targets.
        """
        # iterate over connectons, create ordered dicts
        for c in self.conns.index:
            tmp = c.fluid.val.copy()
            tmp0 = c.fluid.val0.copy()
            tmp_set = c.fluid.val_set.copy()
            c.fluid.val = collections.OrderedDict()
            c.fluid.val0 = collections.OrderedDict()
            c.fluid.val_set = collections.OrderedDict()

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
                    if fluid in tmp_set.keys():
                        if not tmp_set[fluid]:
                            c.fluid.val[fluid] = tmp0[fluid]
                            c.fluid.val0[fluid] = tmp0[fluid]
                            c.fluid.val_set[fluid] = False
                    else:
                        c.fluid.val[fluid] = tmp0[fluid]
                        c.fluid.val0[fluid] = tmp0[fluid]
                        c.fluid.val_set[fluid] = False

                # if fluid not in keys
                else:
                    c.fluid.val[fluid] = 0
                    c.fluid.val0[fluid] = 0
                    c.fluid.val_set[fluid] = False

        # fluid propagation complete for single fluid networks
        if len(self.fluids) == 1:
            return

        # fluid propagation for combustion chambers
        for cp in self.comps.index:
            if isinstance(cp, cmp.combustion_chamber):
                cp.initialise_fluids(self)
                for c in self.comps.loc[cp].o:
                    self.init_target(c, c.t)

        # fluid propagation from set values
        for c in self.conns.index:
            if any(c.fluid.val_set.values()):
                self.init_target(c, c.t)
                self.init_source(c, c.s)

        # fluid propagation starting from all connections
        for c in self.conns.index:
            c.s.initialise_fluids(self)
            c.t.initialise_fluids(self)

    def init_target(self, c, start):
        r"""
        Propagates the fluids towards connection's target with recursive function calls.
        If the target is a sink, a merge or a combustion chamber, the propagation stops.

        Parameters
        ----------
        c : tespy.connections.connection
            Connection to initialise.

        start : tespy.connections.connection
            This connection is the fluid propagation starting point.
            The starting connection is saved to prevent infinite looping.
        """
        if (len(c.t.inlets()) == 1 and len(c.t.outlets()) == 1 or
                isinstance(c.t, cmp.heat_exchanger) or
                isinstance(c.t, cmp.subsys_interface)):

            outc = pd.DataFrame()
            outc['s'] = self.conns.s == c.t
            outc['s_id'] = self.conns.s_id == c.t_id.replace('in', 'out')
            conn, cid = outc['s'] == True, outc['s_id'] == True
            outc = outc.index[conn & cid][0]

            for fluid, x in c.fluid.val.items():
                if not outc.fluid.val_set[fluid]:
                    outc.fluid.val[fluid] = x

            self.init_target(outc, start)

        if isinstance(c.t, cmp.splitter):
            for outconn in self.comps.loc[c.t].o:
                for fluid, x in c.fluid.val.items():
                    if not outconn.fluid.val_set[fluid]:
                        outconn.fluid.val[fluid] = x

                self.init_target(outconn, start)

        if isinstance(c.t, cmp.cogeneration_unit):
            for outconn in self.comps.loc[c.t].o[:2]:
                for fluid, x in c.fluid.val.items():
                    if not outconn.fluid.val_set[fluid]:
                        outconn.fluid.val[fluid] = x

                self.init_target(outconn, start)

        if isinstance(c.t, cmp.drum) and c.t != start:
            start = c.t
            for outconn in self.comps.loc[c.t].o:
                for fluid, x in c.fluid.val.items():
                    if not outconn.fluid.val_set[fluid]:
                        outconn.fluid.val[fluid] = x

                self.init_target(outconn, start)

    def init_source(self, c, start):
        r"""
        Propagates the fluids towards connection's source with recursive function calls.
        If the source is a source or a combustion chamber, the propagation stops.

        Parameters
        ----------
        c : tespy.connections.connection
            Connection to initialise.

        start : tespy.connections.connection
            This connection is the fluid propagation starting point.
            The starting connection is saved to prevent infinite looping.
        """
        if (len(c.s.inlets()) == 1 and len(c.s.outlets()) == 1 or
                isinstance(c.s, cmp.heat_exchanger) or
                isinstance(c.s, cmp.subsys_interface)):

            inc = pd.DataFrame()
            inc['t'] = self.conns.t == c.s
            inc['t_id'] = self.conns.t_id == c.s_id.replace('out', 'in')
            conn, cid = inc['t'] == True, inc['t_id'] == True
            inc = inc.index[conn & cid][0]

            for fluid, x in c.fluid.val.items():
                if not inc.fluid.val_set[fluid]:
                    inc.fluid.val[fluid] = x

            self.init_source(inc, start)

        if isinstance(c.s, cmp.splitter):
            for inconn in self.comps.loc[c.s].i:
                for fluid, x in c.fluid.val.items():
                    if not inconn.fluid.val_set[fluid]:
                        inconn.fluid.val[fluid] = x

                self.init_source(inconn, start)

        if isinstance(c.s, cmp.merge):
            for inconn in self.comps.loc[c.s].i:
                for fluid, x in c.fluid.val.items():
                    if not inconn.fluid.val_set[fluid]:
                        inconn.fluid.val[fluid] = x

                self.init_source(inconn, start)

        if isinstance(c.s, cmp.cogeneration_unit):
            for inconn in self.comps.loc[c.s].i[:2]:
                for fluid, x in c.fluid.val.items():
                    if not inconn.fluid.val_set[fluid]:
                        inconn.fluid.val[fluid] = x

                self.init_source(inconn, start)

        if isinstance(c.s, cmp.drum) and c.s != start:
            start = c.s
            for inconn in self.comps.loc[c.s].i:
                for fluid, x in c.fluid.val.items():
                    if not inconn.fluid.val_set[fluid]:
                        inconn.fluid.val[fluid] = x

                self.init_source(inconn, start)

    def init_properties(self):
        r"""
        Initialises the fluid properties on every connection of the network.

        - Sets standard values for :code:`m0, p0, h0` if not user specified
        - Sets :code:`var = var0` if var_set is False
        - Initialises reference objects
        - Sets initial values for enthalpy at given vapour mass fraction or temperature
        """
        # fluid properties
        for c in self.conns.index:
            for key in ['m', 'p', 'h', 'T', 'x', 'v']:
                if not c.get_attr(key).unit_set and key != 'x':
                    c.get_attr(key).unit = self.get_attr(key + '_unit')
                if key not in ['T', 'x', 'v'] and not c.get_attr(key).val_set:
                    self.init_val0(c, key)
                    c.get_attr(key).val_SI = c.get_attr(key).val0 * self.get_attr(key)[c.get_attr(key).unit]
                elif key not in ['T', 'x', 'v'] and c.get_attr(key).val_set:
                    c.get_attr(key).val_SI = c.get_attr(key).val * self.get_attr(key)[c.get_attr(key).unit]
                elif key == 'T' and c.T.val_set:
                    c.T.val_SI = (c.T.val + self.T[c.T.unit][0]) * self.T[c.T.unit][1]
                elif key == 'x' and c.x.val_set:
                    c.x.val_SI = c.x.val
                elif key == 'v' and c.v.val_set:
                    c.v.val_SI = c.v.val * self.v[c.v.unit]
                else:
                    continue

        # fluid properties with referenced objects
        for c in self.conns.index:
            for key in ['m', 'p', 'h', 'T']:
                if c.get_attr(key).ref_set and not c.get_attr(key).val_set:
                    c.get_attr(key).val_SI = (
                            c.get_attr(key).ref.obj.get_attr(key).val_SI *
                            c.get_attr(key).ref.f + c.get_attr(key).ref.d)

        for c in self.conns.index:
            if c.x.val_set and not c.h.val_set:
                c.h.val_SI = hlp.h_mix_pQ(c.to_flow(), c.x.val_SI)

            if c.T.val_set and not c.h.val_set:
                c.h.val_SI = hlp.h_mix_pT(c.to_flow(), c.T.val_SI)

    def init_val0(self, c, key):
        r"""
        Set starting values for fluid properties. The components classes provide generic starting
        values for its inlets and outlets.

        Parameters
        ----------
        c : tespy.connections.connection
            Connection to initialise.
        """
        if key == 'x' or key == 'T':
            return

        # starting value for mass flow
        if math.isnan(c.get_attr(key).val0) and key == 'm':
            c.get_attr(key).val0 = 1
            return

        # generic starting values for pressure and enthalpy
        if math.isnan(c.get_attr(key).val0):
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

            # change value to specified unit system
            c.get_attr(key).val0 = c.get_attr(key).val0 / self.get_attr(key)[self.get_attr(key + '_unit')]

    def init_csv(self):
        r"""
        Initialise network from .csv file.

        Note
        ----
        This method is used for preprocessing in offdesign-mode using the :code:`design_file` as input file
        and for fluid property and fluid components starting values using the :code:`init_file` as input file.
        """
        if self.mode == 'offdesign':
            self.init_design_file()

        if self.init_file is not None:
            self.init_init_file()

    def init_design_file(self):
        r"""
        Design file reader for preprocessing of offdesign-calculation.
        """
        df = pd.read_csv(self.design_file, index_col=0, delimiter=';', decimal='.')
        for c in self.conns.index:
            # save actual values to temporary variables
            c.m_tmp = c.m.val_SI
            c.p_tmp = c.p.val_SI
            c.h_tmp = c.h.val_SI
            c.fluid_tmp = c.fluid.val.copy()
            # match connection (source, source_id, target, target_id) on
            # connection objects of design file
            conn = (df.loc[df['s'].isin([c.s.label]) & df['t'].isin([c.t.label]) &
                           df['s_id'].isin([c.s_id]) & df['t_id'].isin([c.t_id])])
            if len(conn.index) > 0:
                conn_id = conn.index[0]
                c.m.val_SI = df.loc[conn_id].m * self.m[df.loc[conn_id].m_unit]
                c.p.val_SI = df.loc[conn_id].p * self.p[df.loc[conn_id].p_unit]
                c.h.val_SI = df.loc[conn_id].h * self.h[df.loc[conn_id].h_unit]
                for fluid in self.fluids:
                    c.fluid.val[fluid] = df.loc[conn_id][fluid]

            else:
                msg = 'raise some error'
                ValueError(msg)

        # calculate reference values
        self.processing('pre')

        # reset to actual values
        for c in self.conns.index:
            c.m.val_SI = c.m_tmp
            c.p.val_SI = c.p_tmp
            c.h.val_SI = c.h_tmp
            c.fluid.val = c.fluid_tmp
            del c.m_tmp
            del c.p_tmp
            del c.h_tmp
            del c.fluid_tmp

    def init_init_file(self):
        r"""
        Init file reader for starting value generation of calculation.
        """
        # match connection (source, source_id, target, target_id) on
        # connection objects of design file
        df = pd.read_csv(self.init_file, index_col=0, delimiter=';', decimal='.')
        for c in self.conns.index:
            conn = (df.loc[df['s'].isin([c.s.label]) & df['t'].isin([c.t.label]) &
                           df['s_id'].isin([c.s_id]) & df['t_id'].isin([c.t_id])])
            if len(conn.index) > 0:
                conn_id = conn.index[0]
                # overwrite SI-values with values from init_file, except user specified values
                if not c.m.val_set:
                    c.m.val_SI = df.loc[conn_id].m * self.m[df.loc[conn_id].m_unit]
                if not c.p.val_set:
                    c.p.val_SI = df.loc[conn_id].p * self.p[df.loc[conn_id].p_unit]
                if not c.h.val_set:
                    c.h.val_SI = df.loc[conn_id].h * self.h[df.loc[conn_id].h_unit]
                for fluid in self.fluids:
                    if not c.fluid.val_set[fluid]:
                        c.fluid.val[fluid] = df.loc[conn_id][fluid]

                # overwrite starting values
                c.m.val0 = c.m.val_SI / self.m[c.m.unit]
                c.p.val0 = c.p.val_SI / self.p[c.p.unit]
                c.h.val0 = c.h.val_SI / self.h[c.h.unit]
                c.fluid.val0 = c.fluid.val.copy()

    def init_offdesign(self):
        r"""
        Auto switches components and connections from design to offdesign mode.

        Note
        ----
        **components**

        If :code:`cp.mode == 'auto'` all parameters stated in the components
        attribute :code:`cp.design` will be unset and all parameters stated in
        the components attribute :code:`cp.offdesign` will be set instead.

        The auto-switch can be deactivated by using
        :code:`your_component.set_attr(mode='man')`

        **connections**

        All parameters given in the connections attribute :code:`c.design`
        will be unset.
        """
        for cp in self.comps.index:
            if cp.mode == 'auto':
                for var in cp.design:
                    if cp.get_attr(var).is_set:
                        cp.get_attr(var).set_attr(is_set=False)
                for var in cp.offdesign:
                    if not cp.get_attr(var).is_set:
                        cp.get_attr(var).set_attr(is_set=True)
            cp.comp_init(self)

        for c in self.conns.index:
            for var in c.design:
                if c.get_attr(var).val_set:
                    c.get_attr(var).set_attr(val_set=False)
                if c.get_attr(var).ref_set:
                    c.get_attr(var).set_attr(ref_set=False)

            for var in c.offdesign:
                c.get_attr(var).set_attr(val_set=True)

    def solve(self, mode, init_file=None, design_file=None, max_iter=50, init_only=False):
        r"""
        Solves the network. Tasks:

        - Check network consistency.
        - Initialise calculation and preprocessing.
        - Perform actual calculation.
        - Postprocessing.

        Parameters
        ----------
        mode : String
            Choose from 'design' and 'offdesign'.

        init_file : String
            Specify path to init_file, default: :code:`None`.

        design_file : String
            Specify path to design_file, default: :code:`None`. Required for offdesign calculation!

        max_iter : int
            Maximum number of iterations before calculation stops, default: 50.

        init_only : bool
            Perform initialisation only? default: :code:`False`.

        Note
        ----
        For more information on the solution process have a look at the online documentation
        at tespy.readthedocs.io in the section "using TESPy".
        """
        self.init_file = init_file
        self.design_file = design_file
        self.max_iter = max_iter

        if mode != 'offdesign' and mode != 'design':
            msg = 'Mode must be \'design\' or \'offdesign\'.'
            raise hlp.MyNetworkError(msg)
        else:
            self.mode = mode

        if not self.checked:
            self.check_network()

        self.initialise()

        if init_only:
            return

        self.res = np.array([])

        if self.nwkinfo:
            print('Solving network.')

        self.iter = 0
        # number of variables
        self.num_vars = len(self.fluids) + 3
        self.solve_determination()

        start_time = time.time()
        errmsg = self.solve_loop()
        end_time = time.time()

        if self.iterinfo:
            if self.num_c_vars == 0:
                print('--------+----------+----------+----------+----------+'
                      '---------')
            else:
                print('--------+----------+----------+----------+----------+'
                      '----------+---------')

            msg = ('Total iterations: ' + str(self.iter + 1) + ', '
                   'Calculation time: ' +
                   str(round(end_time - start_time, 1)) + ' s, '
                   'Iterations per second: ' +
                   str(round((self.iter + 1) / (end_time - start_time), 2)))
            print(msg)

        if self.nwkwarn and errmsg is not None:
            print(errmsg)

        if self.lin_dep:
            if self.nwkerr:
                msg = ('##### ERROR #####\n'
                       'singularity in jacobian matrix, frequent reasons are\n'
                       '-> given Temperature with given pressure in two phase '
                       'region, try setting enthalpy instead or '
                       'provide accurate starting value for pressure.\n'
                       '-> given logarithmic temperature differences '
                       'or kA-values for heat exchangers, \n'
                       '-> support better starting values.\n'
                       '-> bad starting value for fuel mass flow of '
                       'combustion chamber, provide small (near to zero, '
                       'but not zero) starting value.')
                print(msg)

        self.processing('post')

        for c in self.conns.index:
            c.T.val_SI = hlp.T_mix_ph(c.to_flow())
            c.v.val_SI = hlp.v_mix_ph(c.to_flow()) * c.m.val_SI
            c.T.val = (c.T.val_SI /
                       self.T[c.T.unit][1] - self.T[c.T.unit][0])
            c.m.val = c.m.val_SI / self.m[c.m.unit]
            c.p.val = c.p.val_SI / self.p[c.p.unit]
            c.h.val = c.h.val_SI / self.h[c.h.unit]
            c.v.val = c.v.val_SI / self.v[c.v.unit]
            c.T.val0 = c.T.val
            c.m.val0 = c.m.val
            c.p.val0 = c.p.val
            c.h.val0 = c.h.val
            c.fluid.val0 = c.fluid.val.copy()

        if self.lin_dep:
            return

        if self.nwkinfo:
            print('Calculation complete.')

    def solve_loop(self):
        r"""
        Loop of the newton algorithm
        """
        if self.iterinfo:
            if self.num_c_vars == 0:
                # iterinfo printout without any custom variables
                msg = ('iter\t| residual | massflow | pressure | enthalpy |'
                       ' fluid\n')
                msg += ('--------+----------+----------+----------+----------+'
                        '---------')

            else:
                # iterinfo printout with custom variables in network
                msg = ('iter\t| residual | massflow | pressure | enthalpy |'
                       ' fluid    | custom\n')
                msg += ('--------+----------+----------+----------+----------+'
                        '----------+---------')

            print(msg)

        self.relax = 1
        for self.iter in range(self.max_iter):

            self.solve_control()
            self.res = np.append(self.res, norm(self.vec_res))

            if self.iterinfo:
                vec = self.vec_z[0:-(self.num_c_vars + 1)]
                msg = (str(self.iter + 1))
                if not self.lin_dep and not math.isnan(norm(self.vec_res)):
                    msg += '\t| ' + '{:.2e}'.format(norm(self.vec_res))
                    msg += ' | ' + '{:.2e}'.format(norm(vec[0::self.num_vars]))
                    msg += ' | ' + '{:.2e}'.format(norm(vec[1::self.num_vars]))
                    msg += ' | ' + '{:.2e}'.format(norm(vec[2::self.num_vars]))
                    ls = []
                    for f in range(len(self.fluids)):
                        ls += vec[3 + f::self.num_vars].tolist()

                    msg += ' | ' + '{:.2e}'.format(norm(ls))
                    if self.num_c_vars > 0:
                        msg += ' | ' + '{:.2e}'.format(norm(
                                self.vec_z[-self.num_c_vars:]))

                else:
                    if math.isnan(norm(self.vec_res)):
                        msg += '\t|      nan'.format(norm(self.vec_res))
                    else:
                        msg += '\t| ' + '{:.2e}'.format(norm(self.vec_res))
                    msg += ' |      nan'
                    msg += ' |      nan'
                    msg += ' |      nan'
                    msg += ' |      nan'
                    if self.num_c_vars > 0:
                        msg += ' |      nan'
                print(msg)

            msg = None

            if ((self.iter > 2 and self.res[-1] < hlp.err ** (1 / 2)) or
                    self.lin_dep):
                return msg

            if self.iter > 20:
                if (all(self.res[(self.iter - 3):] >= self.res[-2]) and
                        self.res[-1] >= self.res[-2]):
                    msg = ('##### WARNING #####\n'
                           'Convergence is making no progress, calculation '
                           'stopped, residual value is '
                           '{:.2e}'.format(norm(self.vec_res)))
                    return msg

        if self.iter == self.max_iter - 1:
            msg = ('##### WARNING #####\n'
                   'Reached maximum iteration count, calculation '
                   'stopped, residual value is '
                   '{:.2e}'.format(norm(self.vec_res)))
            return msg

    def matrix_inversion(self):

        self.lin_dep = True
        try:
            self.vec_z = inv(self.mat_deriv).dot(-self.vec_res)
            self.lin_dep = False
        except np.linalg.linalg.LinAlgError:
            self.vec_z = np.asarray(self.vec_res) * 0
            pass

    def solve_control(self):
        r"""
        Step of the newton algorithm

        - Calculate the residual value for each equation
        - Calculate the jacobian matrix
        - Calculate new values for variables
        - Restrict fluid properties to value ranges
        - Check component parameters for consistency
        """
        num_cols = len(self.conns) * self.num_vars
        self.vec_res = np.zeros([num_cols + self.num_c_vars])
        self.mat_deriv = np.zeros((num_cols + self.num_c_vars,
                                   num_cols + self.num_c_vars,))

        self.solve_connections()
        self.solve_components()
        self.solve_busses()
        self.matrix_inversion()

        # check for linear dependency
        if self.lin_dep:
            return

        # add the increment
        i = 0
        for c in self.conns.index:
            # mass flow, pressure and enthalpy
            if not c.m.val_set:
                c.m.val_SI += self.vec_z[i * (self.num_vars)]
            if not c.p.val_set:
                # this prevents negative pressures
                relax = max(1, -self.vec_z[i * (self.num_vars) + 1] /
                            (0.5 * c.p.val_SI))
                c.p.val_SI += self.vec_z[i * (self.num_vars) + 1] / relax
            if not c.h.val_set:
                c.h.val_SI += self.vec_z[i * (self.num_vars) + 2]

            # fluid vector (only if number of fluids is greater than 1)
            if len(self.fluids) > 1:
                j = 0
                for fluid in self.fluids:
                    # add increment
                    if not c.fluid.val_set[fluid]:
                        c.fluid.val[fluid] += (
                                self.vec_z[i * (self.num_vars) + 3 + j])

                    # keep mass fractions within [0, 1]
                    if c.fluid.val[fluid] < hlp.err:
                        c.fluid.val[fluid] = 0
                    if c.fluid.val[fluid] > 1 - hlp.err:
                        c.fluid.val[fluid] = 1

                    j += 1

            # check the fluid properties for physical ranges
            self.solve_check_props(c)
            i += 1

        # increment for the custom variables
        if self.num_c_vars > 0:
            c_vars = 0
            for cp in self.comps.index:
                for var in cp.vars.keys():
                    pos = var.var_pos

                    # add increment
                    var.val += self.vec_z[self.num_vars * len(self.conns) + c_vars + pos]

                    # keep value within specified value range
                    if var.val < var.min_val:
                        var.val = var.min_val
                    if var.val > var.max_val:
                        var.val = var.max_val

                c_vars += cp.num_c_vars

        # second property check for first three iterations without an init_file
        if self.iter < 3 and self.init_file is None:
            for cp in self.comps.index:
                cp.convergence_check(self)

            for c in self.conns.index:
                self.solve_check_props(c)

    def solve_check_props(self, c):
        r"""
        Checks for invalid fluid properties of pressure, temperature and enthalpy in solution progress and adjusts
        values if necessary.

        Parameters
        ----------
        c : tespy.connections.connection
            Connection to check fluid properties.
        """
        fl = hlp.single_fluid(c.fluid.val)

        if isinstance(fl, str):
            # pressure
            if c.p.val_SI < hlp.memorise.vrange[fl][0] and not c.p.val_set:
                c.p.val_SI = hlp.memorise.vrange[fl][0] * 1.1
            if c.p.val_SI > hlp.memorise.vrange[fl][1] and not c.p.val_set:
                c.p.val_SI = hlp.memorise.vrange[fl][1] * 0.9

            # enthalpy
            hmin = hlp.h_pT(c.p.val_SI, hlp.memorise.vrange[fl][2] * 1.05, fl)
            hmax = hlp.h_pT(c.p.val_SI, hlp.memorise.vrange[fl][3] * 0.99, fl)
            if c.h.val_SI < hmin and not c.h.val_set:
                if c.h.val_SI < 0:
                    c.h.val_SI = hmin / 3
                else:
                    c.h.val_SI = hmin * 3
            if c.h.val_SI > hmax and not c.h.val_set:
                c.h.val_SI = hmax * 0.9

        if self.iter < 3 and self.init_file is None:
            # pressure
            if c.p.val_SI <= self.p_range_SI[0] and not c.p.val_set:
                c.p.val_SI = self.p_range_SI[0]
            if c.p.val_SI >= self.p_range_SI[1] and not c.p.val_set:
                c.p.val_SI = self.p_range_SI[1]

            # enthalpy
            if c.h.val_SI < self.h_range_SI[0] and not c.h.val_set:
                c.h.val_SI = self.h_range_SI[0]
            if c.h.val_SI > self.h_range_SI[1] and not c.h.val_set:
                c.h.val_SI = self.h_range_SI[1]

            # temperature
            if c.T.val_set and not c.h.val_set and not c.p.val_set:
                self.solve_check_temperature(c)

    def solve_check_temperature(self, c):
        r"""
        Checks if temperature is within user specified limits and adjusts
        enthalpy values if necessary.

        Parameters
        ----------
        c : tespy.connections.connection
            Connection to check fluid properties.
        """

        hmin = hlp.h_mix_pT(c.to_flow(), self.T_range_SI[0])
        hmax = hlp.h_mix_pT(c.to_flow(), self.T_range_SI[1])

        if c.h.val_SI < hmin:
            if c.h.val_SI < 0:
                c.h.val_SI = hmin * 0.9
            else:
                c.h.val_SI = hmin * 1.1

        if c.h.val_SI > hmax:
            c.h.val_SI = hmax * 0.95

    def solve_components(self):
        r"""
        Calculates the equations and the partial derivatives of the network's
        components.

        - Iterate through components in network to get residuals and derivatives.
        - Place residual values in residual value vector of the network.
        - Place partial derivatives in jacobian matrix of the network.
        """
        num_cols = len(self.conns) * self.num_vars

        # fetch component equation residuals and component partial derivatives
        data = network.solve_comp(args=([self.comps], ))

        sum_eq = 0
        vec_res = []

        # append residual values to residual value vector
        vec_res += [it for ls in data[0].tolist() for it in ls]
        k = 0
        c_var = 0
        for cp in self.comps.index:

            if (not isinstance(cp, cmp.source) and
                    not isinstance(cp, cmp.sink)):

                i = 0
                num_eq = len(data[1].iloc[k][0])
                inlets = self.comps.loc[cp].i.tolist()
                outlets = self.comps.loc[cp].o.tolist()

                # place derivatives in jacobian matrix
                for c in inlets + outlets:
                    loc = self.conns.index.get_loc(c)
                    self.mat_deriv[sum_eq:sum_eq + num_eq, loc * self.num_vars:(loc + 1) * self.num_vars] = data[1].iloc[k][0][:, i]
                    i += 1

                # derivatives for custom variables
                for j in range(cp.num_c_vars):
                    self.mat_deriv[sum_eq:sum_eq + num_eq, num_cols + c_var] = (data[1].iloc[k][0][:, i + j, :1].transpose()[0])
                    c_var += 1

                sum_eq += num_eq
            k += 1

        self.vec_res[0:self.num_comp_eq] = vec_res

# deprecated
#    def solve_single_component(self, cp):
#        r"""
#        Calculates the equations and the partial derivatives of a single component.
#
#        - Iterate through components in network to get residuals and derivatives.
#        - Place residual values in residual value vector of the network.
#        - Place partial derivatives in jacobian matrix of the network.
#        """
#
#        vec_res = cp.equations()
#        mat_deriv = cp.derivatives(self)
#
#        sum_eq, num_eq, num_cols = self.zero_flag[cp]
#
#        inlets = self.comps.loc[cp].i.tolist()
#        outlets = self.comps.loc[cp].o.tolist()
#
#        i = 0
#        for c in inlets + outlets:
#            loc = self.conns.index.get_loc(c)
#
#            self.mat_deriv[
#                    sum_eq:sum_eq + num_eq,
#                    loc * (self.num_vars):(loc + 1) * self.num_vars
#                    ] = mat_deriv[:, i]
#            i += 1
#
#        c_var = 0
#        for j in range(cp.num_c_vars):
#            self.mat_deriv[
#                    sum_eq:sum_eq + num_eq, num_cols + c_var
#                    ] = mat_deriv[:, i + j, :1].transpose()[0]
#            c_var += 1
#
#        self.vec_res[sum_eq:sum_eq + num_eq] = vec_res

    def solve_comp(args):
        data = args[0]
        return [
                data[0].apply(network.solve_comp_eq, axis=1),
                data[0].apply(network.solve_comp_deriv, axis=1)
        ]

    def solve_comp_eq(cp):
        return cp.name.equations()

    def solve_comp_deriv(cp):
        return [cp.name.derivatives()]

    def solve_connections(self):
        r"""
        Calculates the residual values and the partial derivatives for the network's
        connections equations.

        - Iterate through connections in network to get residuals and derivatives.
        - Place residual values in residual value vector of the network.
        - Place partial derivatives in jacobian matrix of the network.
        """

        # fetch component equation residuals and component partial derivatives
        data = network.solve_conn(args=(self, [self.conns], ))

        row = self.num_comp_eq
        var = {0: 'm', 1: 'p', 2: 'h', 3: 'T', 4: 'x', 5: 'v',
               6: 'm', 7: 'p', 8: 'h', 9: 'T'}
        vec_res = []

        # write data in residual vector and jacobian matrix
        vec_res += [it for ls in data[0].tolist()
                    for it in ls if it is not None]
        k = 0
        for c in self.conns.index:

            # variable counter
            i = 0
            loc = self.conns.index.get_loc(c)
            for it in data[1].iloc[k]:

                if it is not None:

                    # fluid properties
                    self.mat_deriv[row:row + 1, loc * self.num_vars:(loc + 1) * self.num_vars] = it[0, 0]

                    # referenced fluid properties
                    if it[0].shape[0] == 2:
                        c_ref = c.get_attr(var[i]).get_attr('ref')
                        loc_ref = self.conns.index.get_loc(c_ref.obj)
                        self.mat_deriv[row:row + 1, loc_ref * self.num_vars:(loc_ref + 1) * self.num_vars] = it[0, 1]

                    row += 1
                i += 1
            k += 1

        self.vec_res[self.num_comp_eq:row] = vec_res

        # fluid vector
        for c in self.conns.index:

            col = self.conns.index.get_loc(c) * (self.num_vars)
            j = 0
            # specified fluid mass fraction
            for f in self.fluids:
                if c.fluid.val_set[f]:
                    self.mat_deriv[row, col + 3 + j] = 1
                    row += 1
                j += 1

            # specified fluid mass balance
            if c.fluid.balance:
                j = 0
                res = 1
                for f in self.fluids:
                    res -= c.fluid.val[f]
                    self.mat_deriv[row, col + 3 + j] = -1
                    j += 1

                self.vec_res[row] = res
                row += 1

    def solve_conn(args):
        nw, data = args

        return [data[0].apply(network.solve_conn_eq, axis=1, args=(nw,)),
                data[0].apply(network.solve_conn_deriv, axis=1, args=(nw,))]

    def solve_conn_eq(c, nw):
        return [nw.solve_prop_eq(c.name, 'm'),
                nw.solve_prop_eq(c.name, 'p'),
                nw.solve_prop_eq(c.name, 'h'),
                nw.solve_prop_eq(c.name, 'T'),
                nw.solve_prop_eq(c.name, 'x'),
                nw.solve_prop_eq(c.name, 'v'),
                nw.solve_prop_ref_eq(c.name, 'm'),
                nw.solve_prop_ref_eq(c.name, 'p'),
                nw.solve_prop_ref_eq(c.name, 'h'),
                nw.solve_prop_ref_eq(c.name, 'T')]

    def solve_conn_deriv(c, nw):
        return [nw.solve_prop_deriv(c.name, 'm'),
                nw.solve_prop_deriv(c.name, 'p'),
                nw.solve_prop_deriv(c.name, 'h'),
                nw.solve_prop_deriv(c.name, 'T'),
                nw.solve_prop_deriv(c.name, 'x'),
                nw.solve_prop_deriv(c.name, 'v'),
                nw.solve_prop_ref_deriv(c.name, 'm'),
                nw.solve_prop_ref_deriv(c.name, 'p'),
                nw.solve_prop_ref_deriv(c.name, 'h'),
                nw.solve_prop_ref_deriv(c.name, 'T')]

    def solve_busses(self):
        r"""
        Calculates the equations and the partial derivatives for the network's
        busses.

        - Iterate through busses in network to get residuals and derivatives.
        - Place residual values in residual value vector of the network.
        - Place partial derivatives in jacobian matrix of the network.
        """
        row = self.num_comp_eq + self.num_conn_eq
        for b in self.busses:
            if b.P.val_set:
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
                        self.mat_deriv[row, loc * (self.num_vars):
                                       (loc + 1) * self.num_vars] = (
                            deriv[:, j])
                        j += 1

                self.vec_res[row] = b.P.val - P_res

                row += 1

    def solve_prop_eq(self, c, var):
        r"""
        Calculate residuals for given mass flow, volumetric flow, pressure,
        enthalpy, temperature, volumetric flow or vapour mass fraction.

        Parameters
        ----------
        c : tespy.connections.connection
            Connection to calculate the residual value for.

        var : String
            Variable to calculate the residual value for.

        Returns
        -------
        val : float
            Residual value of the corresponding equation (see note).

        Note
        ----
        **mass flow, pressure and enthalpy**

        .. math::
            val = 0

        **temperatures**

        .. math::
            val = T_{j} - T \left( p_{j}, h_{j}, fluid_{j} \right)

        **volumetric flow**

        .. math::
            val = \dot{V}_{j} - v \left( p_{j}, h_{j} \right) \cdot \dot{m}_j

        **vapour mass fraction** *Works with pure fluids only!*

        .. math::
            val = h_{j} - h \left( p_{j}, x_{j}, fluid_{j} \right)
        """
        if var in ['m', 'p', 'h']:
            if c.get_attr(var).get_attr('val_set'):
                return 0
            else:
                return None

        elif var == 'T':
            if c.T.val_set:
                flow = c.to_flow()
                return c.T.val_SI - hlp.T_mix_ph(flow)
            else:
                return None

        elif var == 'v':
            if c.v.val_set:
                flow = c.to_flow()
                return c.v.val_SI - hlp.v_mix_ph(flow) * c.m.val_SI
            else:
                return None

        else:
            if c.x.val_set:
                flow = c.to_flow()
                return c.h.val_SI - hlp.h_mix_pQ(flow, c.x.val_SI)
            else:
                return None

    def solve_prop_ref_eq(self, c, var):
        r"""
        Calculate residuals for referenced mass flow, pressure, enthalpy or temperature.

        Parameters
        ----------
        c : tespy.connections.connection
            Connection to calculate the residual value for.

        var : String
            Variable to calculate the residual value for.

        Returns
        -------
        val : float
            Residual value of the corresponding equation (see note).

        Note
        ----
        **mass flow, pressure and enthalpy**

        .. math::
            val = x_{j} - x_{j,ref} \cdot a + b

        **temperatures**

        .. math::
            val = T \left( p_{j}, h_{j}, fluid_{j} \right) -
            T \left( p_{j}, h_{j}, fluid_{j} \right) \cdot a + b
        """

        if var in ['m', 'p', 'h']:

            if c.get_attr(var).get_attr('ref_set'):
                c_ref = c.get_attr(var).get_attr('ref')
                return (c.get_attr(var).val_SI -
                        (c_ref.obj.get_attr(var).val_SI * c_ref.f + c_ref.d))

            else:
                return None

        else:

            if c.T.ref_set:
                flow = c.to_flow()
                flow_ref = c.T.ref.obj.to_flow()
                return hlp.T_mix_ph(flow) - (hlp.T_mix_ph(flow_ref) *
                                             c.T.ref.f + c.T.ref.d)

            else:
                return None

    def solve_prop_deriv(self, c, var):
        r"""
        Calculate derivatives for given mass flow, pressure, enthalpy,
        temperature, volumetric flow or vapour mass fraction.

        Parameters
        ----------
        c : tespy.connections.connection
            Connection to calculate the residual value for.

        var : String
            Variable to calculate the residual value for.

        Returns
        -------
        deriv : ndarray
            Array of partial derivatives (see note).

        Note
        ----
        **mass flow, pressure and enthalpy**

        .. math::
            J\left(\frac{\partial f_{i}}{\partial m_{j}}\right) = 1\\
            \text{for equation i, connection j}\\
            \text{pressure and enthalpy analogously}

        **temperatures**

        .. math::
            J\left(\frac{\partial f_{i}}{\partial p_{j}}\right) =
            -\frac{dT_{j}}{dp_{j}}\\
            J(\left(\frac{\partial f_{i}}{\partial h_{j}}\right) =
            -\frac{dT_{j}}{dh_{j}}\\
            J\left(\frac{\partial f_{i}}{\partial fluid_{j,k}}\right) =
            - \frac{dT_{j}}{dfluid_{j,k}}
            \; , \forall k \in \text{fluid components}\\
            \text{for equation i, connection j}

        **volumetric flow**

        .. math::
            J\left(\frac{\partial f_{i}}{\partial m_{j}}\right) =
            -v \left( p_{j}, h_{j} \right)\\
            J\left(\frac{\partial f_{i}}{\partial p_{j}}\right) =
            -\frac{dv_{j}}{dp_{j}} \cdot \dot{m}_j\\
            J(\left(\frac{\partial f_{i}}{\partial h_{j}}\right) =
            -\frac{dv_{j}}{dh_{j}} \cdot \dot{m}_j\\

            \; , \forall k \in \text{fluid components}\\
            \text{for equation i, connection j}

        **vapour mass fraction** *Works with pure fluids only!*

        .. math::
            J\left(\frac{\partial f_{i}}{\partial p_{j}}\right) =
            -\frac{\partial h \left( p_{j}, x_{j}, fluid_{j} \right)}
            {\partial p_{j}}\\
            J\left(\frac{\partial f_{i}}{\partial h_{j}}\right) = 1\\
            \text{for equation i, connection j, x: vapour mass fraction}
        """

        if var in ['m', 'p', 'h']:
            if c.get_attr(var).get_attr('val_set'):
                pos = {'m': 0, 'p': 1, 'h': 2}
                deriv = np.zeros((1, 1, self.num_vars))
                deriv[0, 0, pos[var]] = 1
                return deriv
            else:
                return None

        elif var == 'T':
            if c.T.val_set:
                flow = c.to_flow()
                deriv = np.zeros((1, 1, self.num_vars))
                # dT / dp
                deriv[0, 0, 1] = -hlp.dT_mix_dph(flow)
                # dT / dh
                deriv[0, 0, 2] = -hlp.dT_mix_pdh(flow)
                # dT / dFluid
                if len(self.fluids) != 1:
                    deriv[0, 0, 3:] = -hlp.dT_mix_ph_dfluid(flow)
                return deriv
            else:
                return None

        elif var == 'v':
            if c.v.val_set:
                flow = c.to_flow()
                deriv = np.zeros((1, 1, self.num_vars))
                # dv / dm
                deriv[0, 0, 0] = -hlp.v_mix_ph(flow)
                # dv / dp
                deriv[0, 0, 1] = -hlp.dv_mix_dph(flow) * c.m.val_SI
                # dv / dh
                deriv[0, 0, 2] = -hlp.dv_mix_pdh(flow) * c.m.val_SI
                return deriv
            else:
                return None

        else:
            if c.x.val_set:
                flow = c.to_flow()
                deriv = np.zeros((1, 1, self.num_vars))
                # dx / dp
                deriv[0, 0, 1] = -hlp.dh_mix_dpQ(flow, c.x.val_SI)
                # dx / dh
                deriv[0, 0, 2] = 1
                return deriv
            else:
                return None

    def solve_prop_ref_deriv(self, c, var):
        r"""
        Calculate derivatives for referenced mass flow, pressure, enthalpy,
        or temperature.

        Parameters
        ----------
        c : tespy.connections.connection
            Connection to calculate the residual value for.

        var : String
            Variable to calculate the residual value for.

        Returns
        -------
        deriv : ndarray
            Array of partial derivatives (see note).

        Note
        ----
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

        if var in ['m', 'p', 'h']:

            if c.get_attr(var).get_attr('ref_set'):
                pos = {'m': 0, 'p': 1, 'h': 2}
                deriv = np.zeros((1, 2, self.num_vars))
                deriv[0, 0, pos[var]] = 1
                deriv[0, 1, pos[var]] = -c.get_attr(var).ref.f
                return deriv

            else:
                return None

        else:

            if c.T.ref_set:
                flow = c.to_flow()
                flow_ref = c.T.ref.obj.to_flow()
                deriv = np.zeros((1, 2, self.num_vars))
                # dT / dp
                deriv[0, 0, 1] = hlp.dT_mix_dph(flow)
                deriv[0, 1, 1] = -hlp.dT_mix_dph(flow_ref) * c.T.ref.f
                # dT / dh
                deriv[0, 0, 2] = hlp.dT_mix_pdh(flow)
                deriv[0, 1, 2] = -hlp.dT_mix_pdh(flow_ref) * c.T.ref.f
                # dT / dFluid
                if len(self.fluids) != 1:
                    deriv[0, 0, 3:] = hlp.dT_mix_ph_dfluid(flow)
                    deriv[0, 1, 3:] = -hlp.dT_mix_ph_dfluid(flow_ref)
                return deriv

            else:
                return None

    def solve_determination(self):
        r"""
        Checks, if the number of supplied parameters is sufficient for network determination.
        """
        self.num_c_vars = 0
        n = 0
        for cp in self.comps.index:
            self.num_c_vars += cp.num_c_vars
            n += len(cp.equations())

        # number of equations from components
        self.num_comp_eq = n

        n = 0
        for c in self.conns.index:
            n += [c.m.val_set, c.p.val_set, c.h.val_set,
                  c.T.val_set, c.x.val_set, c.v.val_set].count(True)
            n += [c.m.ref_set, c.p.ref_set, c.h.ref_set,
                  c.T.ref_set].count(True)
            n += list(c.fluid.val_set.values()).count(True)
            n += [c.fluid.balance].count(True)

        # number of equations from connections
        self.num_conn_eq = n

        n = 0
        for b in self.busses:
            n += [b.P.val_set].count(True)

        # number of equations from busses
        self.num_bus_eq = n

        n = self.num_comp_eq + self.num_conn_eq + self.num_bus_eq
        if n > self.num_vars * len(self.conns.index) + self.num_c_vars:
            msg = ('You have provided too many parameters: ' +
                   str(self.num_vars * len(self.conns.index) +
                       self.num_c_vars) + ' required, ' + str(n) +
                   ' supplied.')
            raise hlp.MyNetworkError(msg)
        elif n < self.num_vars * len(self.conns.index) + self.num_c_vars:
            msg = ('You have not provided enough parameters: ' +
                   str(self.num_vars * len(self.conns.index) +
                       self.num_c_vars) + ' required, ' + str(n) +
                   ' supplied.')
            raise hlp.MyNetworkError(msg)

# %% pre and post processing

    def processing(self, mode):
        r"""
        Preprocessing or postprocessing for components: calculation of
        components attributes.

        Parameters
        ----------
        mode : String
            Mode selection for pre- or postprocessing: :code:`mode='pre'` or :code:`mode='post'`.
        """
        modes = ['post', 'pre']
        if mode not in modes:
            msg = ('Processing mode must be \'pre\' or \'post\'.')
            raise hlp.MyNetworkError(msg)

        self.comps.apply(network.process_components, axis=1, args=(self, mode,))

        if self.nwkinfo:
            if mode == 'pre':
                print('Preprocessing done.')
            else:
                print('Postprocessing.')

        if mode == 'post':
            # clear fluid property memory
            hlp.memorise.del_memory(self.fluids)
            self.process_busses()
            if self.nwkinfo:
                print('Done.')

    def process_busses(self):
        r"""
        Processing the networks busses
        """
        for b in self.busses:
            b.P.val = 0
            for cp in b.comps.index:

                bus = b.comps.loc[cp]
                val = cp.bus_func(bus)
                b.P.val += val
                if self.mode == 'design':
                    bus.P_ref = val

    def process_components(cols, nw, mode):
        cols.name.calc_parameters(nw, mode)

# %% printing and plotting

    def print_results(self):
        r"""
        Prints the calculations results for components and connections to prompt.
        """
        cp_sort = self.comps.copy()
        # sort components by component type alphabetically
        cp_sort['cp'] = cp_sort.apply(network.get_class_base, axis=1)
        cp_sort['label'] = cp_sort.apply(network.get_props, axis=1, args=('label',))
        cp_sort.drop('i', axis=1, inplace=True)
        cp_sort.drop('o', axis=1, inplace=True)

        pd.options.mode.chained_assignment = None
        for c in cp_sort.cp.unique():
            df = cp_sort[cp_sort['cp'] == c]

            # gather printouts
            cols = []
            for col, val in df.index[0].attr().items():
                if isinstance(val, hlp.dc_cp):
                    if val.get_attr('printout'):
                        cols += [col]

            # any printouts?
            if len(cols) > 0:
                print('##### RESULTS (' + c + ') #####')
                for col in cols:
                    df[col] = df.apply(network.print_components, axis=1, args=(col,))

                df.set_index('label', inplace=True)
                df.drop('cp', axis=1, inplace=True)

                # printout with tabulate
                print(tabulate(df, headers='keys', tablefmt='psql', floatfmt='.2e'))

        # connection properties
        df = pd.DataFrame(columns=['m / (' + self.m_unit + ')',
                                   'p / (' + self.p_unit + ')',
                                   'h / (' + self.h_unit + ')',
                                   'T / (' + self.T_unit + ')'])
        for c in self.conns.index:
            df.loc[c.s.label + ' -> ' + c.t.label] = (
                    [c.m.val_SI / self.m[self.m_unit],
                     c.p.val_SI / self.p[self.p_unit],
                     c.h.val_SI / self.h[self.h_unit],
                     c.T.val_SI / self.T[self.T_unit][1] -
                     self.T[self.T_unit][0]]
                    )
        print(tabulate(df, headers='keys', tablefmt='psql', floatfmt='.3e'))

    def print_components(c, *args):
        return c.name.get_attr(args[0]).val

# %% saving

    def save(self, filename, structure=False):
        r"""
        Saves the results to results file. If structure is True, the network structure is exported.

        Parameters
        ----------
        filename : String
            Path for the results.

        Note
        ----
        File results will be saved to ./filename/results.csv. If you provide :code:`save(structure=True)`,
        all network information will be saved to path ./filename/.
        """
        path = './' + filename + '/'

        # creat path, if non existent
        if not os.path.exists(path):
            os.makedirs(path)

        # save connection properties
        self.save_connections(path + 'results.csv')

        if structure:
            # create path for component folder if non existent
            if not os.path.exists(path + 'comps/'):
                os.makedirs(path + 'comps/')

            # save all network information
            self.save_network(path + 'netw.csv')
            self.save_connections(path + 'conn.csv', structure=True)
            self.save_components(path + 'comps/')
            self.save_busses(path + 'comps/bus.csv')
            self.save_characteristics(path + 'comps/char.csv')

    def save_network(self, fn):
        r"""
        Saves basic network configuration.

        Parameters
        ----------
        fn : String
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

    def save_connections(self, fn, structure=False):
        r"""
        Saves connections to fn, saves network structure data if structure is True.

        - Uses connections object id as row identifier and saves
            * connections source and target as well as
            * properties with references and
            * fluid vector (including user specification if structure is True).
        - Connections source and target are identified by its labels.

        Parameters
        ----------
        fn : String
            Path/filename for the file.
        """
        df = pd.DataFrame()
        # connection id
        df['id'] = self.conns.apply(network.get_id, axis=1)
        # source
        df['s'] = self.conns.apply(network.get_props, axis=1, args=('s', 'label'))
        df['s_id'] = self.conns.apply(network.get_props, axis=1, args=('s_id',))
        # target
        df['t'] = self.conns.apply(network.get_props, axis=1, args=('t', 'label'))
        df['t_id'] = self.conns.apply(network.get_props, axis=1, args=('t_id',))

        if structure:
            # design and offdesign parameters
            df['design'] = self.conns.apply(network.get_props, axis=1, args=('design',))
            df['offdesign'] = self.conns.apply(network.get_props, axis=1, args=('offdesign',))

        cols = ['m', 'p', 'h', 'T', 'x']
        for key in cols:
            # values and units
            df[key] = self.conns.apply(network.get_props, axis=1, args=(key, 'val'))
            df[key + '_unit'] = self.conns.apply(network.get_props, axis=1, args=(key, 'unit'))

            if structure:
                # connection parametrisation
                df[key + '_unit_set'] = self.conns.apply(network.get_props, axis=1, args=(key, 'unit_set'))
                df[key + '0'] = self.conns.apply(network.get_props, axis=1, args=(key, 'val0'))
                df[key + '_set'] = self.conns.apply(network.get_props, axis=1, args=(key, 'val_set'))
                df[key + '_ref'] = self.conns.apply(network.get_props, axis=1, args=(key, 'ref', 'obj',)).astype(str)
                df[key + '_ref'] = df[key + '_ref'].str.extract(r' at (.*?)>', expand=False)
                df[key + '_ref_f'] = self.conns.apply(network.get_props, axis=1, args=(key, 'ref', 'f',))
                df[key + '_ref_d'] = self.conns.apply(network.get_props, axis=1, args=(key, 'ref', 'd',))
                df[key + '_ref_set'] = self.conns.apply(network.get_props, axis=1, args=(key, 'ref_set',))

        for val in self.fluids:
            # fluid mass fraction
            df[val] = self.conns.apply(network.get_props, axis=1, args=('fluid', 'val', val))

            if structure:
                # fluid mass fraction parametrisation
                df[val + '0'] = self.conns.apply(network.get_props, axis=1, args=('fluid', 'val0', val))
                df[val + '_set'] = self.conns.apply(network.get_props, axis=1, args=('fluid', 'val_set', val))

        if structure:
            # fluid balance parametrisation
            df['balance'] = self.conns.apply(network.get_props, axis=1, args=('fluid', 'balance'))

        df.to_csv(fn, sep=';', decimal='.', index=False, na_rep='nan')

    def save_components(self, path):
        r"""
        Saves the components to filename/comps/*.csv

        - Uses components labels as row identifier.
        - Writes:
            * component's incomming and outgoing connections (object id) and
            * component's parametrisation.

        Parameters
        ----------
        path : String
            Path/filename for the file.
        """
        # create / overwrite csv file
        cp_sort = self.comps.copy()
        # component type
        cp_sort['cp'] = cp_sort.apply(network.get_class_base, axis=1)
        # busses
        cp_sort['busses'] = cp_sort.apply(network.get_busses, axis=1, args=(self.busses,))
        cp_sort['bus_param'] = cp_sort.apply(network.get_bus_data, axis=1, args=(self.busses, 'param'))
        cp_sort['bus_P_ref'] = cp_sort.apply(network.get_bus_data, axis=1, args=(self.busses, 'P_ref'))
        cp_sort['bus_char'] = cp_sort.apply(network.get_bus_data, axis=1, args=(self.busses, 'char'))

        pd.options.mode.chained_assignment = None
        for c in cp_sort.cp.unique():
            df = cp_sort[cp_sort['cp'] == c]

            # basic information
            cols = ['label', 'mode', 'design', 'offdesign']
            for col in cols:
                df[col] = df.apply(network.get_props, axis=1, args=(col,))

            # attributes
            for col, dc in df.index[0].attr().items():
                # component characteristics container
                if isinstance(dc, hlp.dc_cc):
                    df[col] = df.apply(network.get_props, axis=1, args=(col, 'func')).astype(str)
                    df[col] = df[col].str.extract(r' at (.*?)>', expand=False)
                    df[col + '_set'] = df.apply(network.get_props, axis=1, args=(col, 'is_set'))
                    df[col + '_method'] = df.apply(network.get_props, axis=1, args=(col, 'method'))
                    df[col + '_param'] = df.apply(network.get_props, axis=1, args=(col, 'param'))

                # component property container
                elif isinstance(dc, hlp.dc_cp):
                    df[col] = df.apply(network.get_props, axis=1, args=(col, 'val'))
                    df[col + '_set'] = df.apply(network.get_props, axis=1, args=(col, 'is_set'))
                    df[col + '_var'] = df.apply(network.get_props, axis=1, args=(col, 'is_var'))

                else:
                    continue

            df.set_index('label', inplace=True)
            df.drop('i', axis=1, inplace=True)
            df.drop('o', axis=1, inplace=True)
            df.to_csv(path + c + '.csv', sep=';', decimal='.', index=True, na_rep='nan')

    def save_busses(self, fn):
        r"""
        Saves the busses parametrisation to filename/comps/bus.csv

        Parameters
        ----------
        fn : String
            Path/filename for the file.
        """
        df = pd.DataFrame({'id': self.busses}, index=self.busses)
        df['id'] = df.apply(network.get_id, axis=1)

        df['label'] = df.apply(network.get_props, axis=1, args=('label',))

        df['P'] = df.apply(network.get_props, axis=1, args=('P', 'val'))
        df['P_set'] = df.apply(network.get_props, axis=1, args=('P', 'val_set'))

        df.to_csv(fn, sep=';', decimal='.', index=False, na_rep='nan')

    def save_characteristics(self, fn):
        r"""
        Saves the busses parametrisation to filename/comps/char.csv

        Parameters
        ----------
        fn : String
            Path/filename for the file.
        """
        cp_sort = self.comps
        cp_sort['cp'] = cp_sort.apply(network.get_class_base, axis=1)

        chars = []
        for c in cp_sort.cp.unique():
            df = cp_sort[cp_sort['cp'] == c]

            for col, dc in df.index[0].attr().items():
                if isinstance(dc, hlp.dc_cc):
                    chars += df.apply(network.get_props, axis=1, args=(col, 'func')).tolist()
                else:
                    continue

        df = pd.DataFrame({'id': self.busses}, index=self.busses)
        for bus in df.index:
            for c in bus.comps.index:
                ch = bus.comps.loc[c].char
                if ch not in chars:
                    chars += [ch]

        df = pd.DataFrame({'id': chars}, index=chars)
        df['id'] = df.apply(network.get_id, axis=1)

        cols = ['x', 'y']
        for val in cols:
            df[val] = df.apply(network.get_props, axis=1, args=(val,))

        df.to_csv(fn, sep=';', decimal='.', index=False, na_rep='nan')

    def get_id(c):
        return str(c.name)[str(c.name).find(' at ') + 4:-1]

    def get_class_base(c):
        return c.name.__class__.__name__

    def get_props(c, *args):
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
                return c.name.get_attr(args[0]).tolist()
            else:
                return c.name.get_attr(args[0])
        else:
            return ''

    def get_busses(c, *args):
        busses = []
        for bus in args[0]:
            if c.name in bus.comps.index:
                busses += [str(bus)[str(bus).find(' at ') + 4:-1]]
        return busses

    def get_bus_data(c, *args):
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
