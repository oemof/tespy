"""
.. module:: networks
    :synopsis: contains logic of the network

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>
"""

import math

import csv

import pandas as pd
from multiprocessing import cpu_count, Pool, freeze_support

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

import sys
import time
import os
from CoolProp.CoolProp import PropsSI as CPPSI


class network:
    """r

    The network class aggregates information on components, connections and
    busses and performs calculation and processing.

    :param fluids: networks fluids
    :type fluids: list
    :returns: no return value
    :raises: - :code:`MyNetworkError`, if the unit system for mass flow
               pressure, enthalpy or temperature is not available
             - :code:`TypeError`, if the ranges for pressure,
               enthalpy or temperature are not stated as list

    **allowed keywords** in kwargs (also see network.attr()):

    - m_unit (*str*)
    - p_unit (*str*), p_range (*list*)
    - h_unit (*str*), h_range (*list*)
    - T_unit (*str*), T_range (*list*)

    **example**

    .. code-block:: python

        from tespy import nwk

        fluid_list = ['Air', 'water']
        nw = nwk.network(fluid_list, p_unit='bar')

    **improvements**

    - add container for units

    """

    def __init__(self, fluids, **kwargs):

        self.checked = False
        self.errors = []
        self.conns = pd.DataFrame(columns=['s', 's_id', 't', 't_id'])

        self.fluids = sorted(fluids)

        # initialise helpers
        for f in self.fluids:
            try:
                hlp.molar_masses[f] = CPPSI('M', f)
            except:
                hlp.molar_masses[f] = 1

            try:
                hlp.gas_constants[f] = CPPSI('GAS_CONSTANT', f)
            except:
                hlp.gas_constants[f] = np.nan

        # initialise memorisation function
        hlp.memorise(self.fluids)

        self.convergence = np.array([0, 0, 0], dtype=object)
        self.busses = []

    # unit systems, calculation is alsways performed with SI-units
        self.m = {
            'kg / s': 1,
            't / h': 3.6
        }
        self.p = {
            'Pa': 1,
            'psi': 6.8948e3,
            'bar': 1e5,
            'MPa': 1e6
        }
        self.h = {
            'J / kg': 1,
            'kJ / kg': 1e3,
            'MJ / kg': 1e6
        }
        self.T = {
            'C': [273.15, 1],
            'F': [459.67, 5 / 9],
            'K': [0, 1]
        }
        self.SI_units = {
              'm': 'kg / s',
              'p': 'Pa',
              'h': 'J / kg',
              'T': 'K'
              }

        # standard unit set
        self.m_unit = self.SI_units['m']
        self.p_unit = self.SI_units['p']
        self.h_unit = self.SI_units['h']
        self.T_unit = self.SI_units['T']

        # standard value range
        self.p_range = [2e3, 300e5]
        self.h_range = [1e3, 7e6]
        self.T_range = [273.16, 1773.15]

        self.p_range_SI = self.p_range[:]
        self.h_range_SI = self.h_range[:]
        self.T_range_SI = self.T_range[:]

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

        # value ranges
        if not isinstance(self.p_range, list):
            msg = ('Specify the value range as list: [p_min, p_max]')
            raise TypeError(msg)
        else:
            self.p_range_SI = np.array(self.p_range) * self.p[self.p_unit]

        if not isinstance(self.h_range, list):
            msg = ('Specify the value range as list: [h_min, h_max]')
            raise TypeError(msg)
        else:
            self.h_range_SI = np.array(self.h_range) * self.h[self.h_unit]

        if not isinstance(self.T_range, list):
            msg = ('Specify the value range as list: [T_min, T_max]')
            raise TypeError(msg)
        else:
            self.T_range_SI = ((np.array(self.T_range) +
                                self.T[self.T_unit][0]) *
                               self.T[self.T_unit][1])

    def __getstate__(self):
        """
        required to pass Pool object within solving loop
        """
        self_dict = self.__dict__.copy()
        del self_dict['pool']
        return self_dict

    def set_attr(self, **kwargs):
        """
        allows adjustments of unit system and fluid property ranges
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

        # value ranges
        if not isinstance(self.p_range, list):
            msg = ('Specify the value range as list: [p_min, p_max]')
            raise TypeError(msg)
        else:
            self.p_range_SI = np.array(self.p_range) * self.p[self.p_unit]

        if not isinstance(self.h_range, list):
            msg = ('Specify the value range as list: [h_min, h_max]')
            raise TypeError(msg)
        else:
            self.h_range_SI = np.array(self.h_range) * self.h[self.h_unit]

        if not isinstance(self.T_range, list):
            msg = ('Specify the value range as list: [T_min, T_max]')
            raise TypeError(msg)
        else:
            self.T_range_SI = ((np.array(self.T_range) +
                                self.T[self.T_unit][0]) *
                               self.T[self.T_unit][1])

    def get_attr(self, key):
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            print('No attribute \"', key, '\" available!')
            return None

    def attr(self):
        return ['m_unit', 'p_unit', 'h_unit', 'T_unit',
                'p_range', 'h_range', 'T_range']

    def add_subsys(self, *args):
        """
        adds connections to the network, calls check_conns method

        :param args: subsystem objects si :code:`add_subsys(s1, s2, s3, ...)`
        :type args: tespy.components.subsystem
        :returns: no return value
        """
        for subsys in args:
            for c in subsys.conns:
                self.add_conns(c)

    def add_conns(self, *args):
        """
        add connections to the network, calls check_conns method

        :param args: connections objects ci :code:`add_conn(c1, c2, c3, ...)`
        :type args: tespy.connection
        :returns: no return value
        """
        for c in args:
            self.check_conns(c)
            self.checked = False

    def del_conns(self, c):
        """
        delets connections from a network

        :param c: connections object to delete
        :type c: tespy.connection
        :returns: no return value
        :raises: :code:`KeyError` if connections object c is not in the network
        """
        self.conns.drop(self.conns.index(c))
        self.checked = False

    def check_conns(self, c):
        """
        checks the networks connections for multiple usage of inlets or outlets
        of components

        :param c: connections object to check
        :type c: tespy.connections.connection
        :returns: no return value
        :raises:
            - :code:`TypeError`, if c is not a connections object
            - :code:`hlp.MyNetworkError`, if components inlet or outlet is
              already connected to another connections object
        """
        if not isinstance(c, con.connection):
            raise TypeError('Must provide tespy.connections.connection objects'
                            ' as parameters.')

        self.conns.loc[c] = [c.s, c.s_id, c.t, c.t_id]

        if self.conns.duplicated(['s', 's_id'])[c]:
            self.conns = self.conns[self.conns.index != c]
            raise hlp.MyNetworkError('Could not add connection to network, '
                                     'source is already in use.')
        if self.conns.duplicated(['t', 't_id'])[c]:
            self.conns = self.conns[self.conns.index != c]
            raise hlp.MyNetworkError('Could not add connection to network, '
                                     'target is already in use.')

    def add_busses(self, *args):
        """
        adds busses to the network, if check_busses returns :code:`True`

        :param args: bus objects bi :code:`add_conn(b1, b2, b3, ...)`
        :type args: tespy.connections.bus
        :returns: no return value
        """
        for b in args:
            if self.check_busses(b):
                self.busses += [b]

    def del_busses(self, b):
        """
        delets busses from a network

        :param b: bus object to delete
        :type b: tespy.connections.bus
        :returns: no return value
        :raises: :code:`KeyError` if bus object b is not in the network
        """
        if b in self.busses:
            self.busses.remove(b)

    def check_busses(self, b):
        """
        checks the networks connections for multiple usage of inlets or outlets
        of components

        :param c: busses object to check
        :type c: tespy.connections.bus
        :returns: bool
        :raises:
            - :code:`TypeError`, if b is not a busses object
            - :code:`hlp.MyNetworkError`, if bus is already in the network
        """
        if isinstance(b, con.bus):
            if b not in self.busses:
                if b.label not in [x.label for x in self.busses]:
                    return True
                else:
                    msg = ('Network already has a bus with the name ' +
                           b.label + '.')
                    raise hlp.MyNetworkError(msg)
            else:
                msg = 'Network contains this bus (' + str(b) + ') already.'
                raise hlp.MyNetworkError(msg)
        else:
            msg = 'Only objects of type bus are allowed in *args.'
            raise TypeError(msg)

        return False

    def check_network(self):
        """
        checks the network consistency: are all components connected?

        - iterates through components of the network
        - substract the number of connections in the network going in
          and out of the component from number of connections the component
          requires.

        :returns: no return value
        :raises: :code:`hlp.MyNetworkError`, if number of connections in the
                 network does not match number of connections required
        """
        for comp in pd.unique(self.conns[['s', 't']].values.ravel()):
            freq = 0
            freq += (self.conns[['s', 't']] == comp).sum().s
            freq += (self.conns[['s', 't']] == comp).sum().t

            if comp.outlets() is not None:
                freq -= len(comp.outlets())
            if comp.inlets() is not None:
                freq -= len(comp.inlets())
            if freq != 0:
                msg = (str(comp) + ' (' + str(comp.label) + ') is missing ' +
                       str(-freq) + ' connections. Make sure all '
                       'inlets and outlets are connected and '
                       'all connections have been added to the '
                       'network.')
                raise hlp.MyNetworkError(msg)

        self.checked = True
        print('Networkcheck successfull.')

    def initialise(self):
        """
        initilialises the network

        - component initlialisation
        - fluid propagation on all connections
        - initilialise fluid properties
        - initialisiation from .csv-files
        - switch components to offdesign mode for offedesign calculation

        :returns: no return value
        """
        msg = ('Have you adjusted the value ranges for pressure, enthalpy'
               ' and temperature according to the specified unit system?')
        print(msg)

        if len(self.fluids) == 0:
            msg = ('Network has no fluids, please specify a list with fluids '
                   'on network creation.')
            raise hlp.MyNetworkError(msg)
        self.init_components()  # build the dataframe for components

        if self.mode == 'offdesign':
            self.init_offdesign()  # characteristics for offdesign

        self.init_fluids()  # start standard fluid initialisation
        self.init_properties()  # start standard property initialisation

        if self.mode == 'offdesign' and self.design_file is None:
            msg = ('Please provide \'design_file\' for every offdesign '
                   'calculation.')
            raise hlp.MyNetworkError(msg)  # must provide design_file
        else:
            self.init_csv()  # initialisation from csv

    def init_components(self):
        """
        writes the networks components into dataframe

        .. note::

            This data is deriven from the network, thus it holds no additional
            information. Instead it is used to simplify the code only.

        dataframe :code:`network.comps`:

        ======================== ============================ =======
         index                    i                            o
        ======================== ============================ =======
         type: component object   type: list                   see i
         value: object id         values: connection objects
        ======================== ============================ =======

        :returns: no return value
        """
        comps = pd.unique(self.conns[['s', 't']].values.ravel())
        self.comps = pd.DataFrame(index=comps, columns=['i', 'o'])

        labels = []
        for comp in self.comps.index:
            comp.comp_init(self)
            s = self.conns[self.conns.s == comp]
            t = self.conns[self.conns.t == comp]
            self.comps.loc[comp] = [t.t_id.sort_values().index,
                                    s.s_id.sort_values().index]
            labels += [comp.label]

        if len(labels) != len(list(set(labels))):
            duplicates = [item for item, count in
                          collections.Counter(labels).items() if count > 1]
            msg = ('All Components must have unique labels, duplicates are: ' +
                   str(duplicates))
            raise hlp.MyNetworkError(msg)

    def init_fluids(self):
        """
        initialises the fluid vector on every connection of the network

        - create fluid vector for every component as dict,
          index: nw.fluids,
          values: 0 if not set by user
        - create fluid_set vector with same logic,
          index: nw.fluids,
          values: False if not set by user
        - calculate fluid vector starting from combustions chambers
        - propagate fluid vector in direction of sources and targets for
          other components

        :returns: no return value
        """

        # iterate over connectons, create ordered dicts
        for c in self.conns.index:
            tmp = c.fluid.val.copy()
            tmp0 = c.fluid.val0.copy()
            tmp_set = c.fluid.val_set.copy()
            c.fluid.val = collections.OrderedDict()
            c.fluid.val0 = collections.OrderedDict()
            c.fluid.val_set = collections.OrderedDict()

            # set values for ordered dict
            for fluid in self.fluids:

                if len(self.fluids) == 1:
                    c.fluid.val[fluid] = 1
                    c.fluid.val0[fluid] = 1

                    try:
                        if tmp_set[fluid]:
                            c.fluid.val_set[fluid] = True
                        else:
                            c.fluid.val_set[fluid] = False
                    except:
                        c.fluid.val_set[fluid] = False

                else:

                    if fluid in tmp.keys():
                        # if fluid in keys and is_set
                        if tmp_set[fluid]:
                            c.fluid.val[fluid] = tmp[fluid]
                            c.fluid.val0[fluid] = tmp[fluid]
                            c.fluid.val_set[fluid] = True
                        # if fluid in keys and not is_set
                        else:
                            c.fluid.val[fluid] = 0
                            c.fluid.val0[fluid] = 0
                            c.fluid.val_set[fluid] = False

                    # if there is a starting value
                    elif fluid in tmp0.keys():
                        if not tmp_set[fluid]:
                            c.fluid.val[fluid] = tmp0[fluid]
                            c.fluid.val0[fluid] = tmp0[fluid]
                            c.fluid.val_set[fluid] = False

                    # if fluid not in keys
                    else:
                        c.fluid.val[fluid] = 0
                        c.fluid.val0[fluid] = 0
                        c.fluid.val_set[fluid] = False

        if len(self.fluids) == 1:
            return

        for cp in self.comps.index:
            if isinstance(cp, cmp.combustion_chamber):
                cp.initialise_fluids(self)
                for c in self.comps.loc[cp].o:
                    self.init_target(c, c.t)

        for c in self.conns.index:
            if any(c.fluid.val_set.values()):
                self.init_target(c, c.t)
                self.init_source(c, c.s)

    def init_target(self, c, start):
        """
        propagates the fluids towards connections target,
        ends when reaching sink, merge or combustion chamber

        :param c: connection to initialise
        :type c: tespy.connections.connection
        :param start: fluid propagation startingpoint, in some cases needed
            to exit the recursion
        :type start: tespy.connections.connection
        :returns: no return value
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

        if isinstance(c.t, cmp.drum) and c.t != start:
            start = c.t
            for outconn in self.comps.loc[c.t].o:
                for fluid, x in c.fluid.val.items():
                    if not outconn.fluid.val_set[fluid]:
                        outconn.fluid.val[fluid] = x

                self.init_target(outconn, start)

    def init_source(self, c, start):
        """
        propagates the fluids towards connections source,
        ends when reaching source, merge or combustion chamber

        :param c: connection to initialise
        :type c: tespy.connections.connection
        :param start: fluid propagation startingpoint, in some cases needed
            to exit the recursion
        :type start: tespy.connections.connection
        :returns: no return value
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

        if isinstance(c.s, cmp.drum) and c.s != start:
            start = c.s
            for inconn in self.comps.loc[c.s].i:
                for fluid, x in c.fluid.val.items():
                    if not inconn.fluid.val_set[fluid]:
                        inconn.fluid.val[fluid] = x

                self.init_source(inconn, start)

    def init_properties(self):
        """
        initialises the fluid properties on every connection of the network

        - sets standard values for :code:`m0, p0, h0` if not user specified
        - sets :code:`var = var0` if var_set is False
        - initialises reference objects
        - performs target fluid propagation from merges
        - sets initial values for enthalpy at given vapour mass fraction or
          temperature

        :returns: no return value
        """
        # fluid properties
        for c in self.conns.index:
            for key in ['m', 'p', 'h', 'T', 'x']:
                if not c.get_attr(key).unit_set and key != 'x':
                    c.get_attr(key).unit = self.get_attr(key + '_unit')
                if key not in ['T', 'x'] and not c.get_attr(key).val_set:
                    self.init_val0(c, key)
                    c.get_attr(key).val_SI = (
                        c.get_attr(key).val0 *
                        self.get_attr(key)[c.get_attr(key).unit])
                elif key not in ['T', 'x'] and c.get_attr(key).val_set:
                    c.get_attr(key).val_SI = (
                        c.get_attr(key).val *
                        self.get_attr(key)[c.get_attr(key).unit])
                elif key == 'T' and c.T.val_set:
                    c.T.val_SI = ((c.T.val + self.T[c.T.unit][0]) *
                                  self.T[c.T.unit][1])
                elif key == 'x' and c.x.val_set:
                    c.x.val_SI = c.x.val
                else:
                    continue

        # fluid properties with referenced objects
        for c in self.conns.index:
            for key in ['m', 'p', 'h', 'T']:
                if c.get_attr(key).ref_set and not c.get_attr(key).val_set:
                    c.get_attr(key).val_SI = (
                            c.get_attr(key).ref.obj.get_attr(key).val_SI *
                            c.get_attr(key).ref.f + c.get_attr(key).ref.d)

        # propate fluids towards merges targets
        # add merge cp to list redo, if number of fluids at merges outlet is
        # still zero
        redo = []
        if len(self.fluids) >= 1:
            for cp in self.comps.index:
                if isinstance(cp, cmp.merge):
                    cp.initialise_fluids(self)
                    for c in self.comps.loc[cp].o:
                        if hlp.num_fluids(c.fluid.val) != 0:
                            self.init_target(c, c.t)
                        else:
                            redo += [cp]

        # propagete fluids towards merges targets of the redo list
        # do this, while the list is not empty
        # if the iteration number is over 50, stop calculation
        i = 0
        while len(redo) != 0:
            for cp in redo:
                cp.initialise_fluids(self)
                for c in self.comps.loc[cp].o:
                    if hlp.num_fluids(c.fluid.val) != 0:
                        self.init_target(c, c.t)
                        redo.remove(cp)

            if i > 50:
                msg = 'Media propagation failed.'
                raise hlp.MyNetworkError(msg)
            i += 1

        for c in self.conns.index:
            if c.x.val_set and not c.h.val_set:
                c.h.val_SI = hlp.h_mix_pQ(c.to_flow(), c.x.val_SI)

            if c.T.val_set and not c.h.val_set:
                c.h.val_SI = hlp.h_mix_pT(c.to_flow(), c.T.val_SI)

    def init_val0(self, c, key):
        """
        sets standard initialisation values for pressure
        values for pressure deriven by

        - attached components or
        - unspecific value (1e5 for pressure)

        :param c: connection to initialise
        :type c: tespy.connections.connection
        :returns: no return value
        """
        if key == 'x' or key == 'T':
            return

        # starting value for mass flow
        if math.isnan(c.get_attr(key).val0) and key == 'm':
            c.get_attr(key).val0 = 1
            return

        # generic starting values
        if math.isnan(c.get_attr(key).val0):
            val_s = c.s.initialise_source(c, key)
            val_t = c.t.initialise_target(c, key)
#            print(key, c.s.label, val_s, c.t.label, val_t)
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
            c.get_attr(key).val0 = (
                c.get_attr(key).val0 /
                self.get_attr(key)[self.get_attr(key + '_unit')])

    def init_csv(self):
        """
        initialise network from .csv file, used for

        - preprocessing before offdesign-calculations (design_file)
        - fluid properties and fluid initialisation (init_file)

        :returns: no return value
        """

        if self.mode == 'offdesign':
            for c in self.conns.index:
                c.m_tmp = c.m.val_SI
                c.p_tmp = c.p.val_SI
                c.h_tmp = c.h.val_SI
                c.fluid_tmp = c.fluid.val.copy()

            df = pd.read_csv(self.design_file, index_col=0, delimiter=';',
                             decimal=self.dec)
            self.conns.apply(network.init_design_file, axis=1,
                             args=(self, df, ))

            # component characteristics creation for offdesign calculation
            self.processing('pre')

            for c in self.conns.index:
                c.m.val_SI = c.m_tmp
                c.p.val_SI = c.p_tmp
                c.h.val_SI = c.h_tmp
                c.fluid.val = c.fluid_tmp

        if self.init_file is not None:
            df = pd.read_csv(self.init_file, index_col=0, delimiter=';',
                             decimal=self.dec)
            self.conns.apply(network.init_init_file, axis=1,
                             args=(self, df, ))

        for c in self.conns.index:
            c.m.val0 = c.m.val_SI / self.m[c.m.unit]
            c.p.val0 = c.p.val_SI / self.p[c.p.unit]
            c.h.val0 = c.h.val_SI / self.h[c.h.unit]
            c.fluid.val0 = c.fluid.val.copy()

    def init_design_file(c, nw, df):
        """
        overwrite variables with values from design file

        :param c: c are the connections of the network
        :type c: landas dataframe index object
        :param nw: tespy network
        :type nw: tespy.networks.network
        :param df: data from csv file
        :type df: pandas.DataFrame
        :returns: no return value
        """
        # match connection (source, source_id, target, target_id) on
        # connection objects of design file
        df_tmp = (df.s == c.s.label).to_frame()
        df_tmp.loc[:, 's_id'] = (df.s_id == c.s_id)
        df_tmp.loc[:, 't'] = (df.t == c.t.label)
        df_tmp.loc[:, 't_id'] = (df.t_id == c.t_id)
        # is True does not work the intended way here!
        s = df_tmp['s'] == True
        s_id = df_tmp['s_id'] == True
        t = df_tmp['t'] == True
        t_id = df_tmp['t_id'] == True
        # overwrite all properties with design file
        conn = df_tmp.index[s & s_id & t & t_id][0]
        c.name.m.val_SI = df.loc[conn].m * nw.m[df.loc[conn].m_unit]
        c.name.p.val_SI = df.loc[conn].p * nw.p[df.loc[conn].p_unit]
        c.name.h.val_SI = df.loc[conn].h * nw.h[df.loc[conn].h_unit]
        for fluid in nw.fluids:
            c.name.fluid.val[fluid] = df.loc[conn][fluid]

    def init_init_file(c, nw, df):
        """
        overwrite non set variables with values from initialisation file

        :param c: c are the connections of the network
        :type c: landas dataframe index object
        :param nw: tespy network
        :type nw: tespy.networks.network
        :param df: data from csv file
        :type df: pandas.DataFrame
        :returns: no return value
        """
        # match connection (source, source_id, target, target_id) on
        # connection objects of design file
        df_tmp = (df.s == c.s.label).to_frame()
        df_tmp.loc[:, 's_id'] = (df.s_id == c.s_id)
        df_tmp.loc[:, 't'] = (df.t == c.t.label)
        df_tmp.loc[:, 't_id'] = (df.t_id == c.t_id)
        # is True does not work the intended way here!
        s = df_tmp['s'] == True
        s_id = df_tmp['s_id'] == True
        t = df_tmp['t'] == True
        t_id = df_tmp['t_id'] == True
        if len(df_tmp.index[s & s_id & t & t_id]) > 0:
            conn = df_tmp.index[s & s_id & t & t_id][0]
            if not c.name.m.val_set:
                c.name.m.val_SI = df.loc[conn].m * nw.m[df.loc[conn].m_unit]
            if not c.name.p.val_set:
                c.name.p.val_SI = df.loc[conn].p * nw.p[df.loc[conn].p_unit]
            if not c.name.h.val_set:
                c.name.h.val_SI = df.loc[conn].h * nw.h[df.loc[conn].h_unit]
            for fluid in nw.fluids:
                if not c.name.fluid.val_set[fluid]:
                    c.name.fluid.val[fluid] = df.loc[conn][fluid]

        if c.name.T.val_set and not c.name.h.val_set:
            c.name.h.val_SI = hlp.h_mix_pT(c.name.to_flow(), c.name.T.val_SI)
        if c.name.x.val_set and not c.name.h.val_set:
            c.name.h.val_SI = hlp.h_mix_pQ(c.name.to_flow(), c.name.x.val_SI)

    def init_offdesign(self):
        """
        auto switches components and connections from design to offdesign mode.

        **components**

        If :code:`cp.mode == 'auto'` all parameters stated in the components
        attribute :code:`cp.design` will be unset and all parameters stated in
        the components attribute :code:`cp.offdesign` will be set instead.

        The auto-switch can be deactivated by using
        :code:`your_component.set_attr(mode='man')`

        **connections**

        All parameters given in the connections attribute :code:`c.design`
        will be unset.

        :returns: no return value
        """
        for cp in self.comps.index:
            if cp.mode == 'auto':
                for var in cp.design:
                    if cp.get_attr(var).is_set:
                        cp.get_attr(var).set_attr(is_set=False)
                for var in cp.offdesign:
                    if not cp.get_attr(var).is_set:
                        cp.get_attr(var).set_attr(is_set=True)

        for c in self.conns.index:
            for var in c.design:
                if c.get_attr(var).val_set:
                    c.get_attr(var).set_attr(val_set=False)
                if c.get_attr(var).ref_set:
                    c.get_attr(var).set_attr(ref_set=False)

            for var in c.offdesign:
                c.get_attr(var).set_attr(val_set=True)

    def solve(self, mode, init_file=None, design_file=None, dec='.',
              max_iter=50, parallel=False, init_only=False):
        """
        solves the network:

        - checks network consistency
        - initialises network
        - starts calculation

        :param mode: calculation mode (design, offdesign)
        :type mode: str
        :param init_file: .csv-file to use for initialisation
        :type init_file: str
        :param design_file: .csv-file containing network design point
        :type design_file: str
        :returns: no return value
        """
        self.init_file = init_file
        self.design_file = design_file
        self.dec = dec
        self.max_iter = max_iter
        self.parallel = parallel

        if mode != 'offdesign' and mode != 'design':
            msg = 'Mode must be \'design\' or \'offdesign\'.'
            raise hlp.MyNetworkError(msg)
        else:
            self.mode = mode

        if not self.checked:
            self.check_network()

        self.initialise()

        print('Network initialised.')

        if init_only:
            return

        # vectors for convergence history (massflow, pressure, enthalpy)
        self.convergence[0] = np.zeros((len(self.conns), 0))
        self.convergence[1] = np.zeros((len(self.conns), 0))
        self.convergence[2] = np.zeros((len(self.conns), 0))
        self.res = np.array([])

        print('Solving network.')

        self.vec_res = []
        self.iter = 0
        self.relax = 1
        self.num_restart = 0
        # number of variables
        self.num_vars = len(self.fluids) + 3
        self.solve_determination()

        # parameters for code parallelisation
        if self.parallel:
            self.cores = cpu_count()
            self.partit = self.cores
            self.comps_split = []
            self.conns_split = []
            self.pool = Pool(self.cores)

            for g, df in self.comps.groupby(np.arange(len(self.comps)) //
                                            (len(self.comps) / self.partit)):
                self.comps_split += [df]

            for g, df in self.conns.groupby(np.arange(len(self.conns)) //
                                            (len(self.conns) / self.partit)):
                self.conns_split += [df]

        else:
            self.partit = 1
            self.comps_split = [self.comps]
            self.conns_split = [self.conns]

        start_time = time.time()
        self.solve_loop()
        end_time = time.time()

        self.processing('post')

        if self.parallel:
            self.pool.close()
            self.pool.join()

        for c in self.conns.index:
            c.T.val_SI = hlp.T_mix_ph(c.to_flow())
            c.T.val = (c.T.val_SI /
                       self.T[c.T.unit][1] - self.T[c.T.unit][0])
            c.m.val = c.m.val_SI / self.m[c.m.unit]
            c.p.val = c.p.val_SI / self.p[c.p.unit]
            c.h.val = c.h.val_SI / self.h[c.h.unit]
            c.T.val0 = c.T.val
            c.m.val0 = c.m.val
            c.p.val0 = c.p.val
            c.h.val0 = c.h.val
            c.fluid.val0 = c.fluid.val.copy()

        print('Calculation complete.')
        msg = ('Total iterations:' + str(self.iter) + ' - '
               'Calculation time:' + str(round(end_time - start_time, 1)) +
               's - Iterations per second:' +
               str(round(self.iter / (end_time - start_time), 2)))
        print(msg)

    def solve_loop(self):
        """
        loop of the newton algorithm

        **Improvememts**
        """

        print('iter\t| residual')
        for self.iter in range(self.max_iter):

            self.convergence[0] = np.column_stack((
                    self.convergence[0], [0] * len(self.conns)))
            self.convergence[1] = np.column_stack((
                    self.convergence[1], [0] * len(self.conns)))
            self.convergence[2] = np.column_stack((
                    self.convergence[2], [0] * len(self.conns)))

            self.solve_control()
            self.res = np.append(self.res, norm(self.vec_res))

            print(self.iter + 1, '\t|', '{:.2e}'.format(norm(self.vec_res)))

            k = 0
            for c in self.conns.index:
                self.convergence[0][k][self.iter] = c.m.val_SI
                self.convergence[1][k][self.iter] = c.p.val_SI
                self.convergence[2][k][self.iter] = c.h.val_SI
                k += 1

            if self.iter > 3:
                self.relax = 1

            self.iter += 1

            # stop calculation after rediculous amount of iterations
            if self.iter > 3:
                if self.res[-1] < hlp.err ** (1 / 2):
                    if self.num_restart > 0:
                        self.convergence[0] = np.column_stack((
                                self.convergence[0], [0] * len(self.conns)))
                        self.convergence[1] = np.column_stack((
                                self.convergence[1], [0] * len(self.conns)))
                        self.convergence[2] = np.column_stack((
                                self.convergence[2], [0] * len(self.conns)))
                        self.vec_z = np.zeros((self.num_vars *
                                               len(self.conns),))
                    break

            if self.iter > 15:
                if (all(self.res[(self.iter - 3):] >= self.res[-2]) and
                        self.res[-1] >= self.res[-2]):
                    print('ERROR: Convergence is making no progress, '
                          'calculation stopped.')
                    if self.num_restart > 0:
                        self.convergence[0] = np.column_stack((
                                self.convergence[0], [0] * len(self.conns)))
                        self.convergence[1] = np.column_stack((
                                self.convergence[1], [0] * len(self.conns)))
                        self.convergence[2] = np.column_stack((
                                self.convergence[2], [0] * len(self.conns)))
                        self.vec_z = np.zeros((self.num_vars *
                                               len(self.conns),))
                    break

    def solve_control(self):
        """
        calculation step of newton algorithm

        - calculate the residual value for each equation
        - calculate the jacobian matrix
        - calculate new values for variables
        - restrict fluid properties to predefined range
        - check component parameters for consistency
        - restart calculation of network with adjusted relaxation factors,
          if linear dependency is detected after first successfull iteration

        :returns: no return value
        :raises: :code:`hlp.MyNetworkError` if network is under-determined.

        **Improvememts**
        """
        self.vec_res = []
        self.solve_components()
        self.solve_connections()
        self.solve_busses()

        lin_dep = True
        try:
            self.vec_z = inv(self.mat_deriv).dot(-np.asarray(self.vec_res))
            lin_dep = False
        except:
            pass

        # check for linear dependency
        if lin_dep:
            if self.iter > 0 and self.num_restart < 3:
                for c in self.conns.index:
                    c.m.val_SI = c.m.val0
                    c.p.val_SI = c.p.val0
                    c.h.val_SI = c.h.val0
                    c.fluid.val = c.fluid.val0.copy()

                self.relax *= 0.8
                self.num_restart += 1
                self.iter = 0
                self.convergence[0] = np.zeros((len(self.conns), 0))
                self.convergence[1] = np.zeros((len(self.conns), 0))
                self.convergence[2] = np.zeros((len(self.conns), 0))
                self.res = np.array([])
                self.vec_res = []
                print('Adjusting relaxation factors and'
                      'restarting calculation.')
                self.solve_loop()
            else:
                msg = ('error calculating the network:\n'
                       'singularity in jacobian matrix, possible reasons are\n'
                       '-> given Temperature with given pressure in two phase '
                       'region, try setting enthalpy instead or '
                       'provide accurate starting value for pressure.\n'
                       '-> given logarithmic temperature differences '
                       'or kA-values for heat exchangers, \n'
                       '-> support better starting values.\n'
                       '-> bad starting value for fuel mass flow of '
                       'combustion chamber, provide small (near to zero, '
                       'but not zero) starting value.')
                raise hlp.MyNetworkError(msg)

        # add increment
        i = 0
        for c in self.conns.index:
            if not c.m.val_set:
                c.m.val_SI += self.vec_z[i * (self.num_vars)] * self.relax
            if not c.p.val_set:
                c.p.val_SI += self.vec_z[i * (self.num_vars) + 1] * self.relax
            if not c.h.val_set:
                c.h.val_SI += self.vec_z[i * (self.num_vars) + 2] * self.relax

            j = 0
            for fluid in self.fluids:
                # add increment
                if not c.fluid.val_set[fluid]:
                    c.fluid.val[fluid] += (
                            self.vec_z[i * (self.num_vars) + 3 + j])

                # prevent bad changes within solution process
                if c.fluid.val[fluid] < hlp.err:
                    c.fluid.val[fluid] = 0
                if c.fluid.val[fluid] > 1 - hlp.err:
                    c.fluid.val[fluid] = 1

                j += 1
            i += 1

        # check properties for consistency
        if self.iter < 3:
            for cp in self.comps.index:
                cp.convergence_check(self)

            for c in self.conns.index:
                self.solve_check_properties(c)

    def solve_check_properties(self, c):
        """
        checks for invalid fluid properties in solution progress and adjusts
        values if necessary

        - check pressure
        - check enthalpy
        - check temperature

        :param c: connection object to check
        :type c: tespy.connections.connection
        :returns: no return value
        """
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

        # make sure, that at given temperatures values stay within feasible
        # range:
        # for mixtures: calculate maximum enthalpy and compare with
        # acutal value
        # for pure fluids:
        # obtain maximum temperature from fluid properties directly
        ############ this part seems to cause bad convergence ############
        if c.T.val_set and not c.h.val_set and not c.p.val_set:
            self.solve_check_temperature(c, 'min')
            self.solve_check_temperature(c, 'max')

    def solve_check_temperature(self, c, pos):
        """
        checks for invalid fluid temperatures in solution progress and adjusts
        values if necessary

        - check if feasible temperatures are within user specified limits and
          adjust limits if necessary

        :param c: connection object to check
        :type c: tespy.connections.connection
        :param pos: check at upper or lower boundary
        :type pos: str
        :returns: no return value
        """

        val = 'T' + pos
        if pos == 'min':
            fac = 1.1
            idx = 0
        else:
            fac = 0.9
            idx = 1

        try:
            T = CPPSI(val, 'T', 0, 'P', 0, hlp.single_fluid(c.fluid.val)) * fac
        except:
            T = self.T_range_SI[idx]

        if pos == 'min':
            if T < self.T_range_SI[idx]:
                T = self.T_range_SI[idx]
        else:
            if T > self.T_range_SI[idx]:
                T = self.T_range_SI[idx]

        p_temp = c.p.val_SI
        if 'INCOMP::' in hlp.single_fluid(c.fluid.val):
            c.p.val_SI = CPPSI('P', 'T', T, 'Q', 0,
                               hlp.single_fluid(c.fluid.val))

        h = hlp.h_mix_pT(c.to_flow(), T)
        c.p.val_SI = p_temp

        if pos == 'min':
            if c.h.val_SI < h:
                c.h.val_SI = h * fac
        else:
            if c.h.val_SI > h:
                c.h.val_SI = h * fac

    def solve_components(self):
        """
        calculates the equations and the partial derivatives for the networks
        components.

        - iterate through components in network to get residuals and
          derivatives
        - append residuals to residual vector
        - place partial derivatives in jacobian matrix

        :returns: no return value

        **Improvements**

        - search a way to speed up locating the data within the matrix
        """
        self.mat_deriv = np.zeros((len(self.conns) * (self.num_vars),
                                   len(self.conns) * (self.num_vars)))

        if self.parallel:
            data = self.solve_parallelize(network.solve_comp, self.comps_split)

        else:
            data = [network.solve_comp(args=(self, self.comps_split, ))]

        sum_eq = 0
        for part in range(self.partit):

            self.vec_res += [it for ls in data[part][0].tolist()
                             for it in ls]
            k = 0
            for cp in self.comps_split[part].index:

                if (not isinstance(cp, cmp.source) and
                        not isinstance(cp, cmp.sink)):

                    i = 0
                    num_eq = len(data[part][1].iloc[k][0])
                    inlets = self.comps.loc[cp].i.tolist()
                    outlets = self.comps.loc[cp].o.tolist()
                    for c in inlets + outlets:

                        loc = self.conns.index.get_loc(c)
                        self.mat_deriv[sum_eq:sum_eq + num_eq,
                                       loc * (self.num_vars):
                                       (loc + 1) * self.num_vars] = (
                                           data[part][1].iloc[k][0][:, i])
                        i += 1
                    sum_eq += num_eq
                k += 1

    def solve_parallelize(self, func, data):
        return self.pool.map(func, [(self, [i],) for i in data])

    def solve_comp(args):
        nw, data = args
        return [
                data[0].apply(network.solve_comp_eq, axis=1, args=(nw,)),
                data[0].apply(network.solve_comp_deriv, axis=1, args=(nw,))
        ]

    def solve_comp_eq(cp, nw):
        return cp.name.equations(nw)

    def solve_comp_deriv(cp, nw):
        return [cp.name.derivatives(nw)]

    def solve_connections(self):
        """
        calculates the equations and the partial derivatives for the networks
        connectons.

        - iterate through connections in network to get residuals and
          derivatives
        - append residuals to residual vector
        - place partial derivatives in jacobian matrix

        :returns: no return value

        **improvements**

        - parallelize fluid calculation
        """
        # parallelization does not seem to speed up calculation for connection
        # properties!!

        if self.parallel:
            df = self.conns_split
            data = self.solve_parallelize(network.solve_conn, df)

        else:
            df = [self.conns]
            data = [network.solve_conn(args=(self, df, ))]

        # write data in residual vector and jacobian matrix
        sum_eq = len(self.vec_res)
        var = {0: 'm', 1: 'p', 2: 'h', 3: 'T', 4: 'x',
               5: 'm', 6: 'p', 7: 'h', 8: 'T'}
        for part in range(self.partit):

            self.vec_res += [it for ls in data[part][0].tolist()
                             for it in ls if it is not None]
            k = 0
            for c in self.conns_split[part].index:

                # variable counter
                i = 0
                loc = self.conns.index.get_loc(c)
                for it in data[part][1].iloc[k]:

                    if it is not None:

                        self.mat_deriv[sum_eq:sum_eq + 1, loc *
                                       self.num_vars: (loc + 1) *
                                       self.num_vars] = it[0, 0]
                        if it[0].shape[0] == 2:

                            c_ref = c.get_attr(var[i]).get_attr('ref')
                            loc_ref = self.conns.index.get_loc(c_ref.obj)
                            self.mat_deriv[sum_eq:sum_eq + 1, loc_ref *
                                           self.num_vars: (loc_ref + 1) *
                                           self.num_vars] = it[0, 1]

                        sum_eq += 1

                    i += 1

                k += 1

        # fluids, no parallelization available yet
        row = sum_eq
        for c in self.conns.index:

            col = self.conns.index.get_loc(c) * (self.num_vars)
            j = 0
            for f in self.fluids:

                if c.fluid.val_set[f]:
                    self.mat_deriv[row, col + 3 + j] = 1
                    self.vec_res += [0]
                    row += 1

                j += 1

            if c.fluid.balance:

                j = 0
                res = 1
                for f in self.fluids:

                    res -= c.fluid.val[f]
                    self.mat_deriv[row, col + 3 + j] = -1
                    j += 1

                self.vec_res += [res]
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
                nw.solve_prop_ref_deriv(c.name, 'm'),
                nw.solve_prop_ref_deriv(c.name, 'p'),
                nw.solve_prop_ref_deriv(c.name, 'h'),
                nw.solve_prop_ref_deriv(c.name, 'T')]

    def solve_busses(self):
        """
        calculates the equations and the partial derivatives for the networks
        busses.

        - iterate through busses in network to get residuals and
          derivatives
        - append residuals to residual vector
        - place partial derivatives in jacobian matrix

        :returns: no return value
        """
        row = len(self.vec_res)
        for b in self.busses:
            if b.P_set:
                P_res = 0
                for cp in b.comps.index:
                    i = self.comps.loc[cp].i.tolist()
                    o = self.comps.loc[cp].o.tolist()

                    P_res += cp.bus_func(i, o) * b.comps.loc[cp].factor
                    deriv = -cp.bus_deriv(i, o)

                    j = 0
                    for c in i + o:
                        loc = self.conns.index.get_loc(c)
                        self.mat_deriv[row, loc * (self.num_vars):
                                       (loc + 1) * self.num_vars] = (
                            deriv[:, j] * b.comps.loc[cp].factor)
                        j += 1

                self.vec_res += [b.P - P_res]

                row += 1

    def solve_prop_eq(self, c, var):
        r"""
        calculate residuals for given mass flow,
        pressure, enthalpy, temperature and vapour mass fraction

        :param c: connections object to apply calculations on
        :type c: tespy.connections.connection
        :param var: variable to perform calculation
        :type var: str
        :returns: no return value

        **mass flow, pressure and enthalpy**

        **equation for numeric values**

        .. math::
            0 = 0

        **derivative for numeric values**

        .. math::
            J\left(\frac{\partial f_{i}}{\partial m_{j}}\right) = 1\\
            \text{for equation i, connection j}\\
            \text{pressure and enthalpy analogously}

        **temperatures**

        **equation for numeric values**

        .. math::
            0 = T_{j} - T \left( p_{j}, h_{j}, fluid_{j} \right)

        **derivative for numeric values**

        .. math::
            J\left(\frac{\partial f_{i}}{\partial p_{j}}\right) =
            -\frac{dT_{j}}{dp_{j}}\\
            J(\left(\frac{\partial f_{i}}{\partial h_{j}}\right) =
            -\frac{dT_{j}}{dh_{j}}\\
            J\left(\frac{\partial f_{i}}{\partial fluid_{j,k}}\right) =
            - \frac{dT_{j}}{dfluid_{j,k}}
            \; , \forall k \in \text{fluid components}\\
            \text{for equation i, connection j}

        **vapour mass fraction**

        .. note::
            works with pure fluids only!

        **equation for numeric values**

        .. math::
            0 = h_{j} - h \left( p_{j}, x_{j}, fluid_{j} \right)

        **derivative for numeric values**

        .. math::
            J\left(\frac{\partial f_{i}}{\partial p_{j}}\right) =
            -\frac{\partial h \left( p_{j}, x_{j}, fluid_{j} \right)}
            {\partial p_{j}}\\
            J(\left(\frac{\partial f_{i}}{\partial h_{j}}\right) = 1\\
            \text{for equation i, connection j, x: vapour mass fraction}
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

        else:
            if c.x.val_set:
                flow = c.to_flow()
                return c.h.val_SI - hlp.h_mix_pQ(flow, c.x.val_SI)
            else:
                return None

    def solve_prop_ref_eq(self, c, var):
        r"""
        calculate residuals for referenced mass flow,
        pressure, enthalpy and temperature

        :param c: connections object to apply calculations on
        :type c: tespy.connections.connection
        :param var: variable to perform calculation
        :type var: str
        :returns: no return value

        **mass flow, pressure and enthalpy**

        .. math::
            0 = m_{j} - m_{j,ref} \cdot a + b

        **temperatures**

        .. math::
            0 = T \left( p_{j}, h_{j}, fluid_{j} \right) -
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
        calculate derivatives for given mass flow,
        pressure, enthalpy, temperature and vapour mass fraction

        :param c: connections object to apply calculations on
        :type c: tespy.connections.connection
        :param var: variable to perform calculation
        :type var: str
        :returns: no return value

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

        **vapour mass fraction**

        .. note::
            works with pure fluids only!

        .. math::
            J\left(\frac{\partial f_{i}}{\partial p_{j}}\right) =
            -\frac{\partial h \left( p_{j}, x_{j}, fluid_{j} \right)}
            {\partial p_{j}}\\
            J(\left(\frac{\partial f_{i}}{\partial h_{j}}\right) = 1\\
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

        else:

            if c.x.val_set:

                flow = c.to_flow()
                deriv = np.zeros((1, 1, self.num_vars))
                deriv[0, 0, 1] = -hlp.dh_mix_dpQ(flow, c.x.val_SI)
                deriv[0, 0, 2] = 1
                return deriv

            else:
                return None

    def solve_prop_ref_deriv(self, c, var):
        r"""
        calculate residuals for referenced mass flow,
        pressure, enthalpy and temperature

        :param c: connections object to apply calculations on
        :type c: tespy.connections.connection
        :param var: variable to perform calculation
        :type var: str
        :returns: no return value

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
        calculates the number of given parameters

        :returns: no return value
        :raises: :code:`MyNetworkError`
        """
        n = 0
        for cp in self.comps.index:
            n += len(cp.equations(self))

        for c in self.conns.index:
            n += [c.m.val_set, c.p.val_set, c.h.val_set,
                  c.T.val_set, c.x.val_set].count(True)
            n += [c.m.ref_set, c.p.ref_set, c.h.ref_set,
                  c.T.ref_set].count(True)
            n += list(c.fluid.val_set.values()).count(True)
            n += [c.fluid.balance].count(True)

        for b in self.busses:
            n += [b.P_set].count(True)

        if n > self.num_vars * len(self.conns.index):
            msg = ('You have provided too many parameters:',
                   self.num_vars * len(self.conns.index), 'required, ',
                   n, 'given.')
            raise hlp.MyNetworkError(msg)
        elif n < self.num_vars * len(self.conns.index):
            msg = ('You have not provided enough parameters:',
                   self.num_vars * len(self.conns.index), 'required, ',
                   n, 'given.')
            raise hlp.MyNetworkError(msg)
        else:
            return

# %% pre and post processing

    def processing(self, mode):
        """
        preprocessing or postprocessing for components: calculation of
        components attributes

        :param mode: mode selection (pre/post) for pre- or postprocessing
        :type mode: str
        :returns: no return value
        """
        modes = ['post', 'pre']
        if mode not in modes:
            msg = ('Processing mode must be \'pre\' for offdesign preparation '
                   'or \'post\'.')
            raise hlp.MyNetworkError(msg)

        self.comps.apply(network.process_components, axis=1,
                         args=(self, mode,))

        if mode == 'pre':
            print('Preprocessing done.')
        else:
            print('Postprocessing.')
            # clear fluid property memory
            hlp.memorise.del_memory(self.fluids)
            self.process_busses()
            print('Done.')

    def process_busses(self):
        """
        processing the networks busses

        :returns: no return value
        """
        for b in self.busses:
            b.P = 0
            for cp in b.comps.index:
                i = self.comps.loc[cp].i.tolist()
                o = self.comps.loc[cp].o.tolist()

                b.P += cp.bus_func(i, o) * b.comps.loc[cp].factor

    def process_components(cols, nw, mode):
        """
        postprocessing: calculate components attributes

        :param cols: cols are the components of the network
        :type cols: landas dataframe index object
        :returns: no return value
        """
        cols.name.calc_parameters(nw, mode)

# %% printing and plotting

    def print_results(self):
        """
        prints the calculations results for components and connections

        :returns: no return value

        **Improvements**

        adjust number of decimal places according to specified units
        """

        P_res = [x.P for x in self.busses if x.label == 'P_res']
        Q_diss = [x.P for x in self.busses if x.label == 'Q_diss']

        if len(P_res) != 0 and len(Q_diss) != 0:
            print('process key figures')
            print('eta_th = ' + str(1 - sum(Q_diss) /
                  (sum(P_res) + sum(Q_diss))))
            print('eps_hp = ' + str(abs(sum(Q_diss)) / sum(P_res)))
            print('eps_cm = ' + str(abs(sum(Q_diss)) / sum(P_res) - 1))

        msg = 'Do you want to print the components parammeters?'
        if hlp.query_yes_no(msg):
            self.comps.apply(network.print_components, axis=1,
                             args=(self,))

        msg = 'Do you want to print the connections parammeters?'
        if hlp.query_yes_no(msg):
            df = pd.DataFrame(columns=['m / (' + self.m_unit + ')',
                                       'p / (' + self.p_unit + ')',
                                       'h / (' + self.h_unit + ')',
                                       'T / (' + self.T_unit + ')'])
            for c in self.conns.index:
                df.loc[c.s.label + ' -> ' + c.t.label] = (
                        ['{:.2e}'.format(c.m.val_SI / self.m[self.m_unit]),
                         '{:.2e}'.format(c.p.val_SI / self.p[self.p_unit]),
                         '{:.4e}'.format(c.h.val_SI / self.h[self.h_unit]),
                         '{:.2e}'.format(c.T.val_SI / self.T[self.T_unit][1] -
                                         self.T[self.T_unit][0])]
                        )
            print(df)

    def print_components(cols, nw):
        """
        postprocessing: calculate components attributes and print them to
        prompt

        :param cols: cols are the components of the network
        :type cols: landas dataframe index object
        :returns: no return value
        """

        cols.name.print_parameters(nw)

    def plot_convergence(self):
        """
        plots the convergence history of all mass flows, pressures and
        enthalpies as absolute values

        :returns: no return value
        """

        num_flows = len(self.conns.index)
        cm = plt.get_cmap('autumn')
        cNorm = colors.Normalize(vmin=0, vmax=num_flows - 1)
        scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
        color = [scalarMap.to_rgba(i) for i in range(num_flows)]

        num_steps = len(self.convergence[0][0])
        x = np.linspace(1, num_steps, num_steps)

        i = 0
        subplt_label = ['massflow', 'pressure', 'enthalpy']
        f, axarr = plt.subplots(3, sharex=True)
        f.suptitle('convergence history', fontsize=16)
        for subplt in axarr:
            subplt.grid()
            subplt.title.set_text(subplt_label[i])
            i += 1

        k = 0
        for c in self.conns.index:
            i = 0
            for prop in self.convergence:
                    axarr[i].plot(x, prop[k][:],
                                  color=color[k],
                                  label=c.s.label + ' -> ' + c.t.label)
                    i += 1
            k += 1

        axarr[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
        f.subplots_adjust(right=0.8, hspace=0.2)
        plt.show()

# %% saving

    def save(self, filename, **kwargs):
        """
        saves the results in two files:

        - results file and
        - components file

        :param filename: suffix for the .csv-file
        :type filename: str
        :returns: no return value
        """

        path = './' + filename + '/'

        if not os.path.exists(path):
            os.makedirs(path)

        self.save_connections(path + 'results.csv')

        if kwargs.get('structure', False):

            if not os.path.exists(path + 'comps/'):
                os.makedirs(path + 'comps/')

            self.save_network(path + 'netw.csv')
            self.save_connections(path + 'conn.csv', structure=True)
            self.save_components(path + 'comps/')
            self.save_busses(path + 'comps/bus.csv')
            self.save_characteristics(path + 'comps/char.csv')

    def save_network(self, fn):
        """
        saves basic network configuration

        :param fn: filename
        :type fn: str
        :returns: no return value
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
        """
        saves connections to fn, saves network structure data if structure is
        True

        - uses connections object id as row identifier and saves
            * connections source and target as well as
            * properties with references and
            * fluid vector (including user specification if structure is True)
        - connections source and target are identified by its labels

        :param fn: filename
        :type fn: str
        :returns: no return value
        """

        df = pd.DataFrame()
        df['id'] = self.conns.apply(network.get_id, axis=1)

        df['s'] = self.conns.apply(network.get_props, axis=1,
                                   args=('s', 'label'))
        df['s_id'] = self.conns.apply(network.get_props, axis=1,
                                      args=('s_id',))

        df['t'] = self.conns.apply(network.get_props, axis=1,
                                   args=('t', 'label'))
        df['t_id'] = self.conns.apply(network.get_props, axis=1,
                                      args=('t_id',))

        if structure:
            df['design'] = self.conns.apply(network.get_props, axis=1,
                                            args=('design',))
            df['offdesign'] = self.conns.apply(network.get_props, axis=1,
                                               args=('offdesign',))

        cols = ['m', 'p', 'h', 'T', 'x']
        for key in cols:
            df[key] = self.conns.apply(network.get_props, axis=1,
                                       args=(key, 'val'))
            df[key + '_unit'] = self.conns.apply(network.get_props, axis=1,
                                                 args=(key, 'unit'))

            if structure:
                df[key + '_unit_set'] = self.conns.apply(
                        network.get_props, axis=1, args=(key, 'unit_set'))
                df[key + '0'] = self.conns.apply(network.get_props, axis=1,
                                                 args=(key, 'val0'))
                df[key + '_set'] = self.conns.apply(network.get_props, axis=1,
                                                    args=(key, 'val_set'))
                df[key + '_ref'] = self.conns.apply(
                        network.get_props, axis=1,
                        args=(key, 'ref', 'obj',)).astype(str)
                df[key + '_ref'] = df[key + '_ref'].str.extract(r' at (.*?)>',
                                                                expand=False)
                df[key + '_ref_f'] = self.conns.apply(
                        network.get_props, axis=1, args=(key, 'ref', 'f',))
                df[key + '_ref_d'] = self.conns.apply(
                        network.get_props, axis=1, args=(key, 'ref', 'd',))
                df[key + '_ref_set'] = self.conns.apply(
                        network.get_props, axis=1, args=(key, 'ref_set',))

        for val in sorted(self.fluids):
            df[val] = self.conns.apply(network.get_props, axis=1,
                                       args=('fluid', 'val', val))

            if structure:
                df[val + '0'] = self.conns.apply(network.get_props, axis=1,
                                                 args=('fluid', 'val0', val))
                df[val + '_set'] = self.conns.apply(
                        network.get_props, axis=1,
                        args=('fluid', 'val_set', val))

        if structure:
            df['balance'] = self.conns.apply(network.get_props, axis=1,
                                             args=('fluid', 'balance'))

        df.to_csv(fn, sep=';', decimal='.', index=False, na_rep='nan')

    def save_components(self, path):
        """
        saves the components to filename/comps/*.csv

        - uses components labels as row identifier
        - writes:
            * components incomming and outgoing connections (object id)
            * components parametrisation

        :param path: path to the files
        :type path: str
        :returns: no return value
        """

        # create / overwrite csv file
        cp_sort = self.comps
        cp_sort['cp'] = cp_sort.apply(network.get_class_base, axis=1)
        cp_sort['busses'] = cp_sort.apply(network.get_busses, axis=1,
                                          args=(self.busses,))
        cp_sort['bus_factors'] = cp_sort.apply(network.get_bus_factors,
                                               axis=1,
                                               args=(self.busses,))

        pd.options.mode.chained_assignment = None
        for c in cp_sort.cp.unique():
            df = cp_sort[cp_sort['cp'] == c]

            cols = ['label', 'mode', 'design', 'offdesign']
            for col in cols:
                df[col] = df.apply(network.get_props, axis=1,
                                   args=(col,))

            for col, dc in df.index[0].attr_prop().items():
                if isinstance(dc, hlp.dc_cc):
                    df[col] = df.apply(network.get_props, axis=1,
                                       args=(col, 'func')).astype(str)
                    df[col] = df[col].str.extract(r' at (.*?)>', expand=False)
                    df[col + '_set'] = df.apply(network.get_props, axis=1,
                                                args=(col, 'is_set'))
                    df[col + '_method'] = df.apply(network.get_props, axis=1,
                                                   args=(col, 'method'))
                    df[col + '_param'] = df.apply(network.get_props, axis=1,
                                                  args=(col, 'param'))

                elif isinstance(dc, hlp.dc_cp):
                    df[col] = df.apply(network.get_props, axis=1,
                                       args=(col, 'val'))
                    df[col + '_set'] = df.apply(network.get_props, axis=1,
                                                args=(col, 'is_set'))
                    df[col + '_var'] = df.apply(network.get_props, axis=1,
                                                args=(col, 'is_var'))

                else:
                    continue

            df.set_index('label', inplace=True)
            df.drop('i', axis=1, inplace=True)
            df.drop('o', axis=1, inplace=True)
            df.to_csv(path + c + '.csv', sep=';', decimal='.',
                      index=True, na_rep='nan')

    def save_busses(self, fn):
        """
        saves the busses parametrisation

        :param fn: filename
        :type fn: str
        :returns: no return value
        """

        df = pd.DataFrame({'id': self.busses}, index=self.busses)
        df['id'] = df.apply(network.get_id, axis=1)

        df['label'] = df.apply(network.get_props, axis=1,
                               args=('label',))
        df['P'] = df.apply(network.get_props, axis=1,
                           args=('P',))
        df['P_set'] = df.apply(network.get_props, axis=1,
                               args=('P_set',))

        df.to_csv(fn, sep=';', decimal='.', index=False, na_rep='nan')

    def save_characteristics(self, fn):
        """
        saves the busses parametrisation to filename_bus.csv

        - uses connections object id as row identifier
            * properties and property_set (False/True)
            * referenced objects
            * fluids and fluid_set vector
            * connections source and target
        - connections source and target are identified by its labels

        :param filename: suffix for the .csv-file
        :type filename: str
        :returns: no return value
        """

        cp_sort = self.comps
        cp_sort['cp'] = cp_sort.apply(network.get_class_base, axis=1)

        chars = []
        for c in cp_sort.cp.unique():
            df = cp_sort[cp_sort['cp'] == c]

            for col, dc in df.index[0].attr_prop().items():
                if isinstance(dc, hlp.dc_cc):
                    chars += df.apply(network.get_props, axis=1,
                                      args=(col, 'func')).tolist()
                else:
                    continue

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

                if args[0] == 'fluid' and args[1] != 'balance':
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

    def get_bus_factors(c, *args):
        factors = []
        for bus in args[0]:
            if c.name in bus.comps.index:
                factors += [bus.comps.loc[c.name].factor]
        return factors
