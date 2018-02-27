"""
.. module:: networks
    :synopsis: contains logic of the network

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>
"""

import math

import csv

import pandas as pd
from multiprocessing import cpu_count, Pool

import numpy as np
from numpy.linalg import inv
from numpy.linalg import norm

from tespy.components import components as cmp
from tespy import connections as con
from tespy import helpers as hlp

import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as colors

import collections

import time
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
    - memo (*bool*), initialise fluid property memorisation?

    **example**

    .. code-block:: python

        from tespy import nwk

        fluid_list = ['Air', 'water']
        nw = nwk.network(fluid_list, p_unit='bar')

    **improvements**

    """

    # unit systems, calculation is alsways performed with SI-units
    # zu conversion/converter umbenennen?
    m = {
        'kg / s': 1,
        't / h': 3.6
    }
    p = {
        'Pa': 1,
        'psi': 6.8948e3,
        'bar': 1e5,
        'MPa': 1e6
    }
    h = {
        'J / kg': 1,
        'kJ / kg': 1e3,
        'MJ / kg': 1e6
    }
    T = {
        'C': [273.15, 1],
        'F': [459.67, 5 / 9],
        'K': [0, 1]
    }
    SI_units = {
          'm': 'kg / s',
          'p': 'Pa',
          'h': 'J / kg',
          'T': 'K'
          }

    def __init__(self, fluids, **kwargs):

        self.checked = False
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

        # standard unit set
        self.m_unit = network.SI_units['m']
        self.p_unit = network.SI_units['p']
        self.h_unit = network.SI_units['h']
        self.T_unit = network.SI_units['T']

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
        if self.m_unit not in network.m.keys():
            msg = ('Allowed units for mass flow are: ' +
                   str(network.m.keys()))
            raise hlp.MyNetworkError(msg)

        if self.p_unit not in network.p.keys():
            msg = ('Allowed units for pressure are: ' +
                   str(network.p.keys()))
            raise hlp.MyNetworkError(msg)

        if self.h_unit not in network.h.keys():
            msg = ('Allowed units for enthalpy are: ' +
                   str(network.h.keys()))
            raise hlp.MyNetworkError(msg)

        if self.T_unit not in network.T.keys():
            msg = ('Allowed units for temperature are: ' +
                   str(network.T.keys()))
            raise hlp.MyNetworkError(msg)

        # value ranges
        if not isinstance(self.p_range, list):
            msg = ('Specify the value range as list: [p_min, p_max]')
            raise TypeError(msg)
        else:
            self.p_range_SI[0] = self.p_range[0] * network.p[self.p_unit]
            self.p_range_SI[1] = self.p_range[1] * network.p[self.p_unit]

        if not isinstance(self.h_range, list):
            msg = ('Specify the value range as list: [h_min, h_max]')
            raise TypeError(msg)
        else:
            self.h_range_SI[0] = self.h_range[0] * network.h[self.h_unit]
            self.h_range_SI[1] = self.h_range[1] * network.h[self.h_unit]

        if not isinstance(self.T_range, list):
            msg = ('Specify the value range as list: [T_min, T_max]')
            raise TypeError(msg)
        else:
            self.T_range_SI[0] = ((self.T_range[0] +
                                   network.T[self.T_unit][0]) *
                                  network.T[self.T_unit][1])
            self.T_range_SI[1] = ((self.T_range[1] +
                                   network.T[self.T_unit][0]) *
                                  network.T[self.T_unit][1])

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
        if self.m_unit not in network.m.keys():
            msg = ('Allowed units for mass flow are: ' +
                   str(network.m.keys()))
            raise hlp.MyNetworkError(msg)

        if self.p_unit not in network.p.keys():
            msg = ('Allowed units for pressure are: ' +
                   str(network.p.keys()))
            raise hlp.MyNetworkError(msg)

        if self.h_unit not in network.h.keys():
            msg = ('Allowed units for enthalpy are: ' +
                   str(network.h.keys()))
            raise hlp.MyNetworkError(msg)

        if self.T_unit not in network.T.keys():
            msg = ('Allowed units for temperature are: ' +
                   str(network.T.keys()))
            raise hlp.MyNetworkError(msg)

        # value ranges
        if not isinstance(self.p_range, list):
            msg = ('Specify the value range as list: [p_min, p_max]')
            raise TypeError(msg)
        else:
            self.p_range_SI[0] = self.p_range[0] * network.p[self.p_unit]
            self.p_range_SI[1] = self.p_range[1] * network.p[self.p_unit]

        if not isinstance(self.h_range, list):
            msg = ('Specify the value range as list: [h_min, h_max]')
            raise TypeError(msg)
        else:
            self.h_range_SI[0] = self.h_range[0] * network.h[self.h_unit]
            self.h_range_SI[1] = self.h_range[1] * network.h[self.h_unit]

        if not isinstance(self.T_range, list):
            msg = ('Specify the value range as list: [T_min, T_max]')
            raise TypeError(msg)
        else:
            self.T_range_SI[0] = ((self.T_range[0] +
                                   network.T[self.T_unit][0]) *
                                  network.T[self.T_unit][1])
            self.T_range_SI[1] = ((self.T_range[1] +
                                   network.T[self.T_unit][0]) *
                                  network.T[self.T_unit][1])

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
        for c in self.conns.index:
            if any(c.fluid_set.values()):
                fluid_tmp = c.fluid.copy()
                fluid_set_tmp = c.fluid_set.copy()
                c.fluid = {}
                c.fluid_set = {}
            else:
                fluid_tmp = {}
                fluid_set_tmp = {}
            for fluid in self.fluids:
                if fluid in fluid_tmp.keys():
                    c.fluid[fluid] = fluid_tmp[fluid]
                    if fluid_set_tmp[fluid]:
                        c.fluid_set[fluid] = True
                    else:
                        c.fluid_set[fluid] = False
                else:
                    if len(self.fluids) == 1:
                        c.fluid[fluid] = 1
                    else:
                        c.fluid[fluid] = 0

                    c.fluid_set[fluid] = False

        if len(self.fluids) == 1:
            return

        for cp in self.comps.index:
            cp.comp_init(self)
            if isinstance(cp, cmp.combustion_chamber):
                cp.initialise_fluids(self)
                for c in self.comps.loc[cp].o:
                    self.init_target(c, c.t)

        for c in self.conns.index:
            if any(c.fluid_set.values()):
                self.init_target(c, c.t)
                self.init_source(c, c.s)

        for c in self.conns.index:
            if hlp.num_fluids(c.fluid) != 0:
                self.init_target(c, c.t)
                self.init_source(c, c.s)

        for c in self.conns.index[::-1]:
            if hlp.num_fluids(c.fluid) != 0:
                self.init_source(c, c.s)
                self.init_target(c, c.t)

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

            for fluid, x in c.fluid.items():
                if not outc.fluid_set[fluid]:
                    outc.fluid[fluid] = x

            self.init_target(outc, start)

        if isinstance(c.t, cmp.splitter):
            for outconn in self.comps.loc[c.t].o:
                for fluid, x in c.fluid.items():
                    if not outconn.fluid_set[fluid]:
                        outconn.fluid[fluid] = x

                self.init_target(outconn, start)

        if isinstance(c.t, cmp.drum) and c.t != start:
            start = c.t
            for outconn in self.comps.loc[c.t].o:
                for fluid, x in c.fluid.items():
                    if not outconn.fluid_set[fluid]:
                        outconn.fluid[fluid] = x

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

            for fluid, x in c.fluid.items():
                if not inc.fluid_set[fluid]:
                    inc.fluid[fluid] = x

            self.init_source(inc, start)

        if isinstance(c.s, cmp.splitter):
            for inconn in self.comps.loc[c.s].i:
                for fluid, x in c.fluid.items():
                    if not inconn.fluid_set[fluid]:
                        inconn.fluid[fluid] = x

                self.init_source(inconn, start)

        if isinstance(c.s, cmp.drum) and c.s != start:
            start = c.s
            for inconn in self.comps.loc[c.s].i:
                for fluid, x in c.fluid.items():
                    if not inconn.fluid_set[fluid]:
                        inconn.fluid[fluid] = x

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
        for c in self.conns.index:
            if not isinstance(c.m, con.ref):
                c.m *= network.m[self.m_unit]
            if not c.m_set:
                c.m = c.m0

            self.init_p0(c)
            if not isinstance(c.p, con.ref):
                c.p *= network.p[self.p_unit]
            if not c.p_set:
                c.p = c.p0

            self.init_h0(c)
            if not isinstance(c.h, con.ref):
                c.h *= network.h[self.h_unit]
            if not c.h_set:
                c.h = c.h0

        self.init_refobj()

        # propate fluids towards merges targets
        # add merge cp to list redo, if number of fluids at merges outlet is
        # still zero
        redo = []
        if len(self.fluids) >= 1:
            for cp in self.comps.index:
                if isinstance(cp, cmp.merge):
                    cp.initialise_fluids(self)
                    for c in self.comps.loc[cp].o:
                        if hlp.num_fluids(c.fluid) != 0:
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
                    if hlp.num_fluids(c.fluid) != 0:
                        self.init_target(c, c.t)
                        redo.remove(cp)

            if i > 50:
                msg = 'Media propagation failed.'
                raise hlp.MyNetworkError(msg)
            i += 1

        for c in self.conns.index:
            if c.x_set:
                flow = [c.m, c.p, c.h, c.fluid]
                c.h = hlp.h_mix_pQ(flow, c.x)

            if c.T_set and not isinstance(c.T, con.ref):
                c.T = ((c.T + network.T[self.T_unit][0]) *
                       network.T[self.T_unit][1])
                flow = [c.m, c.p, c.h, c.fluid]
                c.h = hlp.h_mix_pT(flow, c.T)

    def init_p0(self, c):
        """
        sets standard initialisation values for pressure
        values for pressure deriven by

        - user specification,
        - attached components or
        - unspecific value (1e5 for pressure)

        :param c: connection to initialise
        :type c: tespy.connections.connection
        :returns: no return value
        """
        if math.isnan(c.p0):
            source_p = c.s.initialise_source_p(c)
            target_p = c.t.initialise_target_p(c)
            if source_p == 0 and target_p == 0:
                c.p0 = 1e5
            elif source_p == 0:
                c.p0 = target_p
            elif target_p == 0:
                c.p0 = source_p
            else:
                c.p0 = (source_p + target_p) / 2

        else:
            c.p0 *= network.p[self.p_unit]

    def init_h0(self, c):
        """
        sets standard initialisation values for enthalpy
        values for enthalpy deriven by

        - user specification,
        - attached components or
        - unspecific value (1e6 for enthalpy)

        :param c: connection to initialise
        :type c: tespy.connections.connection
        :returns: no return value
        """
        if math.isnan(c.h0):
            source_h = c.s.initialise_source_h(c)
            target_h = c.t.initialise_target_h(c)
            if source_h == 0 and target_h == 0:
                c.h0 = 1e6
            elif source_h == 0:
                c.h0 = target_h
            elif target_h == 0:
                c.h0 = source_h
            else:
                c.h0 = (source_h + target_h) / 2

        else:
            c.h0 *= network.h[self.h_unit]

    def init_refobj(self):
        """
        initialise reference objects as fluid properties

        - adds a ref-attributes to connection (of type connections)
        - sets starting value for the variable

        :returns: no return value
        """
        for c in self.conns.index:
            if isinstance(c.m, con.ref):
                c.m_ref = c.m
                c.m = c.m_ref.obj.m * c.m.f + c.m.d
            if isinstance(c.p, con.ref):
                c.p_ref = c.p
                c.p = c.p_ref.obj.p * c.p.f + c.p.d
            if isinstance(c.h, con.ref):
                c.h_ref = c.h
                c.h = c.h_ref.obj.h * c.h.f + c.h.d
            if isinstance(c.T, con.ref):
                c.T_ref = c.T

    def init_csv(self):
        """
        initialise network from .csv file, used for

        - preprocessing before offdesign-calculations (design_file)
        - fluid properties and fluid initialisation (init_file)

        :returns: no return value
        """

        if self.mode == 'offdesign':
            for c in self.conns.index:
                c.m_tmp = c.m
                c.p_tmp = c.p
                c.h_tmp = c.h
                c.fluid_tmp = c.fluid.copy()

            df = pd.read_csv(self.design_file, index_col=0, delimiter=';',
                             decimal=self.dec)
            self.conns.apply(network.init_design_file, axis=1,
                             args=(self, df, ))

            # component characteristics creation for offdesign calculation
            self.processing('pre')

            for c in self.conns.index:
                c.m = c.m_tmp
                c.p = c.p_tmp
                c.h = c.h_tmp
                c.fluid = c.fluid_tmp

        if self.init_file is not None:
            df = pd.read_csv(self.init_file, index_col=0, delimiter=';',
                             decimal=self.dec)
            self.conns.apply(network.init_init_file, axis=1,
                             args=(self, df, ))

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
        c.name.m = df.loc[conn].m * network.m[nw.m_unit]
        c.name.p = df.loc[conn].p * network.p[nw.p_unit]
        c.name.h = df.loc[conn].h * network.h[nw.h_unit]
        for fluid in nw.fluids:
            c.name.fluid[fluid] = df.loc[conn][fluid]

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
            if not c.name.m_set:
                c.name.m = df.loc[conn].m * network.m[nw.m_unit]
            if not c.name.p_set:
                c.name.p = df.loc[conn].p * network.p[nw.p_unit]
            if not c.name.h_set:
                c.name.h = df.loc[conn].h * network.h[nw.h_unit]
            for fluid in nw.fluids:
                if not c.name.fluid_set[fluid]:
                    c.name.fluid[fluid] = df.loc[conn][fluid]

        if c.name.T_set and not isinstance(c.name.T, con.ref):
            c.name.h = hlp.h_mix_pT(c.name.as_list(), c.name.T)
        if c.name.x_set:
            c.name.h = hlp.h_mix_pQ(c.name.as_list(), c.name.x)

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
                    if cp.__dict__[var + '_set']:
                        cp.__dict__[var + '_set'] = False
                for var in cp.offdesign:
                    if not cp.__dict__[var + '_set']:
                        cp.__dict__[var + '_set'] = True

        for c in self.conns.index:
            for var in c.design:
                if c.__dict__[var + '_set']:
                    c.__dict__[var + '_set'] = False

            for var in c.offdesign:
                c.__dict__[var + '_set'] = True

    def solve(self, mode, init_file=None, design_file=None, dec='.',
              max_iter=50, parallel=False):
        """
        solves the network:

        - checks network consistency
        - initialises network
        - calculates network

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
        self.dec = '.'
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

        # vectors for convergence history (massflow, pressure, enthalpy)
        self.convergence[0] = np.zeros((len(self.conns), 0))
        self.convergence[1] = np.zeros((len(self.conns), 0))
        self.convergence[2] = np.zeros((len(self.conns), 0))
        self.res = np.array([])

        print('Solving network.')

        start_time = time.time()

        self.vec_res = []
        self.iter = 0
        self.relax = 1
        # number of variables
        self.num_vars = len(self.fluids) + 3
        self.solve_determination()

        # parameters for code parallelisation
        self.cores = cpu_count()
        self.partit = self.cores
        self.comps_split = []
        self.pool = Pool(self.cores)

        for g, df in self.comps.groupby(np.arange(len(self.comps)) //
                                        (len(self.comps) / self.partit)):
            self.comps_split += [df]

        print('iter\t| residual')
        for self.iter in range(max_iter):
            self.convergence[0] = np.column_stack((self.convergence[0],
                                                   [0] * len(self.conns)))
            self.convergence[1] = np.column_stack((self.convergence[1],
                                                   [0] * len(self.conns)))
            self.convergence[2] = np.column_stack((self.convergence[2],
                                                   [0] * len(self.conns)))
            self.solve_loop()
            self.res = np.append(self.res, norm(self.vec_res))

            print(self.iter + 1, '\t|', '{:.2e}'.format(norm(self.vec_res)))

            k = 0
            for c in self.conns.index:
                self.convergence[0][k][self.iter] = c.m
                self.convergence[1][k][self.iter] = c.p
                self.convergence[2][k][self.iter] = c.h
                k += 1

            self.iter += 1

            # stop calculation after rediculous amount of iterations
            if self.iter > 3:
                if (all(self.res[(self.iter - 4):] <
                        hlp.err ** (1 / 2))):
                    break

            if self.iter > 20:
                if (all(self.res[(self.iter - 5):] >= self.res[-4]) and
                        self.res[-1] > self.res[-2]):
                    print('Convergence is making no progress, '
                          'calculation stopped.')
                    break
        end_time = time.time()

        self.processing('post')

        self.pool.close()
        self.pool.join()

        for c in self.conns.index:
            c.T = (hlp.T_mix_ph([c.m, c.p, c.h, c.fluid]) /
                   network.T[self.T_unit][1] - network.T[self.T_unit][0])
            c.m /= network.m[self.m_unit]
            c.p /= network.p[self.p_unit]
            c.h /= network.h[self.h_unit]
            c.m0 = c.m
            c.p0 = c.p
            c.h0 = c.h

        print('Calculation complete.')
        msg = ('Total iterations:' + str(self.iter) + ' - '
               'Calculation time:' + str(round(end_time - start_time, 1)) +
               's - Iterations per second:' +
               str(round(self.iter / (end_time - start_time), 2)))
        print(msg)

    def solve_loop(self):
        """
        calculates the solution of the network with n-dimensional newton
        algorithm:

        - calculate the residual value for each equation
        - calculate the jacobian matrix
        - calculate new values for variables
        - restrict fluid properties to predefined range
        - check component parameters for consistency

        :returns: no return value
        :raises: :code:`hlp.MyNetworkError` if network is under-determined.

        **Improvememts**

        - add the possibility for user specified property ranges:
          min and max pressure, enthalpy and temperature
        """
        self.vec_res = []
        self.solve_components()
        self.solve_connections()
        self.solve_busses()

        try:
            vec_z = inv(self.mat_deriv).dot(-np.asarray(self.vec_res))
        except:
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
            if not c.m_set or hasattr(c, 'm_ref'):
                c.m += vec_z[i * (self.num_vars)]
            if not c.p_set or hasattr(c, 'p_ref'):
                c.p += vec_z[i * (self.num_vars) + 1] * self.relax
            if not c.h_set or hasattr(c, 'h_ref'):
                c.h += vec_z[i * (self.num_vars) + 2]

            l = 0
            for fluid in c.fluid.keys():
                # add increment
                if not c.fluid_set[fluid]:
                    c.fluid[fluid] += vec_z[i * (self.num_vars) + 3 + l]

                # prevent bad changes within solution process
                if c.fluid[fluid] < hlp.err:
                    c.fluid[fluid] = 0
                if c.fluid[fluid] > 1 - hlp.err:
                    c.fluid[fluid] = 1

                l += 1
            i += 1

            self.solve_check_properties(c)

        # check properties for consistency
        if self.init_file is None and self.iter < 3:
            for cp in self.comps.index:
                cp.convergence_check(self)

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
        if c.p <= self.p_range_SI[0]:
            c.p = self.p_range_SI[0]
        if c.p >= self.p_range_SI[1]:
            c.p = self.p_range_SI[1]

        # enthalpy
        if c.h < self.h_range_SI[0]:
            c.h = self.h_range_SI[0]
        if c.h > self.h_range_SI[1]:
            c.h = self.h_range_SI[1]

        # make sure, that at given temperatures values stay within feasible
        # range:
        # for mixtures: calculate maximum enthalpy and compare with
        # acutal value
        # for pure fluids:
        # obtain maximum temperature from fluid properties directly
        if self.iter < 3 and c.T_set and not c.h_set:
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
            fac = 1.01
            idx = 0
        else:
            fac = 0.99
            idx = 1

        try:
            T = CPPSI(val, 'T', 0, 'P', 0, hlp.single_fluid(c.fluid)) * fac
        except:
            T = self.T_range_SI[idx]

        if pos == 'min':
            if T < self.T_range_SI[idx]:
                T = self.T_range_SI[idx]
        else:
            if T > self.T_range_SI[idx]:
                T = self.T_range_SI[idx]

        p_temp = c.p
        if 'INCOMP::' in hlp.single_fluid(c.fluid):
            c.p = CPPSI('P', 'T', T, 'Q', 0, hlp.single_fluid(c.fluid))

        h = hlp.h_mix_pT(c.as_list(), T)
        c.p = p_temp

        if pos == 'min':
            if c.h < h:
                c.h = h * fac
        else:
            if c.h > h:
                c.h = h * fac

    def solve_components(self):
        """
        calculates the equations and the partial derivatives for the networks
        components.

        - iterate through components in network to get residuals and
          derivatives
        - append residuals to residual vector
        - place partial derivatives in jacobian matrix

        :returns: no return value
        """
        self.debug = []
        self.mat_deriv = np.zeros((len(self.conns) * (self.num_vars),
                                   len(self.conns) * (self.num_vars)))

        if self.parallel:

            data = self.solve_parallelize(network.solve_eq)

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

        else:

            sum_eq = 0
            for cp in self.comps.index.values:

                self.vec_res += cp.equations(self)
                vec_deriv = cp.derivatives(self)
                if (not isinstance(cp, cmp.source) and
                        not isinstance(cp, cmp.sink)):

                    i = 0
                    num_eq = len(vec_deriv)
                    inlets = self.comps.loc[cp].i.tolist()
                    outlets = self.comps.loc[cp].o.tolist()
                    # insert derivatives row - wise
                    # loc is the location of c in the jacobian matrix
                    for c in inlets + outlets:

                        loc = self.conns.index.get_loc(c)
                        self.mat_deriv[sum_eq:sum_eq + num_eq,
                                       loc * (self.num_vars):
                                       (loc + 1) * self.num_vars] = (
                                           vec_deriv[:, i])
                        i += 1
                    sum_eq += num_eq

    def solve_parallelize(self, func):
        return self.pool.map(func, [(self, [i],) for i in self.comps_split])

    def solve_eq(args):
        nw, data = args
        return [
                data[0].apply(network.eq, axis=1, args=(nw,)),
                data[0].apply(network.deriv, axis=1, args=(nw,))
        ]

    def eq(cp, nw):
        return cp.name.equations(nw)

    def deriv(cp, nw):
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
        """
        row = len(self.vec_res)

        for c in self.conns.index:
            col = self.conns.index.get_loc(c) * (self.num_vars)

            if c.m_set:
                self.solve_m_p_h(c, row, col, 'm')
                row += 1
            if c.p_set:
                self.solve_m_p_h(c, row, col, 'p')
                row += 1
            if c.h_set:
                self.solve_m_p_h(c, row, col, 'h')
                row += 1
            if c.T_set:
                self.solve_T(c, row, col)
                row += 1
            if c.x_set:
                self.solve_x(c, row, col)
                row += 1

            l = 0
            for f in c.fluid.keys():
                if c.fluid_set[f]:
                    self.mat_deriv[row, col + 3 + l] = 1
                    self.vec_res += [0]
                    row += 1
                    self.debug += [c]
                l += 1

            if c.fluid_balance:
                l = 0
                res = 1
                for f, x in c.fluid.items():
                    res -= x
                    self.mat_deriv[row, col + 3 + l] = -1
                    l += 1

                self.vec_res += [res]
                row += 1
                self.debug += [c]

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
                for c in b.comps.index:
                    i = self.comps.loc[c].i[0]
                    o = self.comps.loc[c].o[0]

                    factor = b.comps.loc[c].factor
                    self.comps.loc[c].i.tolist()

                    P_res += i.m * (o.h - i.h) * factor

                    col = self.conns.index.get_loc(i) * self.num_vars
                    self.mat_deriv[row, col] = -(o.h - i.h) * factor
                    self.mat_deriv[row, col + 2] = i.m * factor

                    col = self.conns.index.get_loc(o) * self.num_vars
                    self.mat_deriv[row, col + 2] = -i.m * factor

                self.vec_res += [b.P - P_res]

                row += 1

    def solve_m_p_h(self, c, row, col, var):
        r"""
        calculate residuals and partial derivatives for given mass flow,
        pressure and/or enthalpy

        :param c: connections object to apply calculations on
        :type c: tespy.connections.connection
        :param row: index of row to insert into jacobian matrix
        :type row: int
        :param col: index of column for connection c in jacobian matrix
        :type col: int
        :param var: variable to perform calculation (m, p, h)
        :type var: str
        :returns: no return value

        **equation for numeric values**

        .. math::
            0 = 0

        **derivative for numeric values**

        .. math::
            J\left(\frac{\partial f_{i}}{\partial m_{j}}\right) = 1\\
            \text{for equation i, connection j}\\
            \text{pressure and enthalpy analogously}

        **equation for referenced values**

        .. math::
            0 = m_{j} - m_{j,ref} \cdot a + b

        **derivatives for referenced values**

        .. math::
            J\left(\frac{\partial f_{i}}{\partial m_{j}}\right) = 1\\
            J\left(\frac{\partial f_{i}}{\partial m_{j,ref}}\right) = - a\\
            \text{for equation i, connection j}\\
            \text{pressure and enthalpy analogously}
        """
        pos = {'m': 0, 'p': 1, 'h': 2}
        if hasattr(c, (var + '_ref')):
            c_ref = c.get_attr(var + '_ref')
            ref_col = self.conns.index.get_loc(c_ref.obj) * (self.num_vars)
            self.mat_deriv[row, col + pos[var]] = 1
            self.mat_deriv[row, ref_col + pos[var]] = -c_ref.f
            self.vec_res += [c.get_attr(var) -
                             (c_ref.obj.get_attr(var) * c_ref.f + c_ref.d)]
        else:
            self.mat_deriv[row, col + pos[var]] = 1
            self.vec_res += [0]

        self.debug += [c]

    def solve_T(self, c, row, col):
        r"""
        calculate residuals and partial derivatives for given temperature

        - the calculation of partial derivatives of temperature to fluids
          is skipped if residuals become very small (speeds up calculation)

        :param c: connections object to apply calculations on
        :type c: tespy.connections.connection
        :param row: index of row to insert into jacobian matrix
        :type row: int
        :param col: index of column for connection c in jacobian matrix
        :type col: int
        :returns: no return value

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

        **equation for numeric values**

        .. math::
            0 = T \left( p_{j}, h_{j}, fluid_{j} \right) -
            T \left( p_{j}, h_{j}, fluid_{j} \right) \cdot a + b

        **derivative for numeric values**

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
        flow = c.as_list()
        # if is reference to other connection
        if hasattr(c, 'T_ref'):
            ref_col = self.conns.index.get_loc(c.T_ref.obj) * (self.num_vars)
            flow_ref = c.T_ref.obj.as_list()
            # residual
            res = hlp.T_mix_ph(flow) - (hlp.T_mix_ph(flow_ref) *
                                        c.T_ref.f + c.T_ref.d)
            self.vec_res += [res]

            # derivatives
            # dT / dp
            self.mat_deriv[row, col + 1] = hlp.dT_mix_dph(flow)
            self.mat_deriv[row, ref_col + 1] = (
                -hlp.dT_mix_dph(flow_ref) * c.T_ref.f)
            # dT / dh
            self.mat_deriv[row, col + 2] = hlp.dT_mix_pdh(flow)
            self.mat_deriv[row, ref_col + 2] = (
                -hlp.dT_mix_pdh(flow_ref) * c.T_ref.f)
            # dT / dFluid
            if len(self.fluids) != 1:
                dT_dfluid = hlp.dT_mix_ph_dfluid(flow)
                dT_dfluid_ref = hlp.dT_mix_ph_dfluid(flow_ref)
                for i in range(self.num_vars - 3):
                    self.mat_deriv[row, col + 3 + i] = dT_dfluid[i]
                    self.mat_deriv[row, ref_col + 3 + i] = (
                        -dT_dfluid_ref[i] * c.T_ref.f)

            # this part is speeding up calculation, but might lead to
            # singularities
#            if abs(res) > hlp.err ** 3:
                # dT / dp
#                self.mat_deriv[row, col + 1] = hlp.dT_mix_dph(flow)
#                self.mat_deriv[row, ref_col + 1] = (
#                    -hlp.dT_mix_dph(flow_ref) * c.T_ref.f)
                # dT / dh
#                self.mat_deriv[row, col + 2] = hlp.dT_mix_pdh(flow)
#                self.mat_deriv[row, ref_col + 2] = (
#                    -hlp.dT_mix_pdh(flow_ref) * c.T_ref.f)
                # dT / dFluid
                # stop calculating temperature derivatives if residuals in
                # temperature are small (very time consuming)
#                dT_dfluid = hlp.dT_mix_ph_dfluid(flow)
#                dT_dfluid_ref = hlp.dT_mix_ph_dfluid(flow_ref)
#                for i in range(self.num_vars - 3):
#                    self.mat_deriv[row, col + 3 + i] = dT_dfluid[i]
#                    self.mat_deriv[row, ref_col + 3 + i] = (
#                        -dT_dfluid_ref[i] * c.T_ref.f)
#            else:
                # values do not matter, just set differently to prevent
                # singularity
#                self.mat_deriv[row, col + 1] = 0.5
#                self.mat_deriv[row, ref_col + 1] = -c.T_ref.f * 0.7
#                self.mat_deriv[row, col + 2] = 2
#                self.mat_deriv[row, ref_col + 2] = -c.T_ref.f * 1.6
#                for i in range(self.num_vars - 3):
#                    self.mat_deriv[row, col + 3 + i] = 1.2
#                    self.mat_deriv[row, ref_col + 3 + i] = -0.8

        # if value is set for temperature
        else:
            # residual
            res = c.T - hlp.T_mix_ph(flow)
            self.vec_res += [res]
            # derivatives
            self.mat_deriv[row, col + 1] = -hlp.dT_mix_dph(flow)
            self.mat_deriv[row, col + 2] = -hlp.dT_mix_pdh(flow)

            if len(self.fluids) != 1:
                dT_dfluid = hlp.dT_mix_ph_dfluid(flow)
                self.mat_deriv[row, col + 3:col + self.num_vars] = -dT_dfluid

            # this part is speeding up calculation, but might lead to
            # singularities
#            if abs(res) > hlp.err ** 3:
#                self.mat_deriv[row, col + 1] = -hlp.dT_mix_dph(flow)
#                self.mat_deriv[row, col + 2] = -hlp.dT_mix_pdh(flow)
#
#                dT_dfluid = hlp.dT_mix_ph_dfluid(flow)
#                self.mat_deriv[row, col + 3:col + self.num_vars] = -dT_dfluid
#            else:
                # values do not matter, just set differently to prevent
                # singularity
#                self.mat_deriv[row, col + 1] = -0.1
#                self.mat_deriv[row, col + 2] = -2
#                self.mat_deriv[row, col + 3:col + self.num_vars] = -0.7

        self.debug += [c]

    def solve_x(self, c, row, col):
        r"""
        calculate residuals and partial derivatives for given vapour mass
        fraction

        :param c: connections object to apply calculations on
        :type c: tespy.connections.connection
        :param row: index of row to insert into jacobian matrix
        :type row: int
        :param col: index of column for connection c in jacobian matrix
        :type col: int
        :returns: no return value

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
        flow = c.as_list()
        self.mat_deriv[row, col + 1] = -hlp.dh_mix_dpQ(flow, c.x)
        self.mat_deriv[row, col + 2] = 1
        self.vec_res += [(c.h - hlp.h_mix_pQ(flow, c.x))]

        self.debug += [c]

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
            n += [c.m_set, c.p_set, c.h_set, c.T_set, c.x_set].count(True)
            n += list(c.fluid_set.values()).count(True)
            n += [c.fluid_balance].count(True)

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

            self.comps.apply(network.process_components, axis=1,
                             args=(self, mode,))
            self.process_busses()
            print('Done.')

    def process_busses(self):
        """
        processing the networks busses

        :returns: no return value
        """
        for b in self.busses:
            b.P = 0
            for c in b.comps.index:
                i = self.comps.loc[c].i[0]
                o = self.comps.loc[c].o[0]

                factor = b.comps.loc[c].factor

                b.P += i.m * (o.h - i.h) * factor

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

        # convert to SI-units
        for c in self.conns.index:
            c.m *= network.m[self.m_unit]
            c.p *= network.p[self.p_unit]
            c.h *= network.h[self.h_unit]
            c.T = ((c.T + network.T[self.T_unit][0]) *
                   network.T[self.T_unit][1])

        P_res = [x.P for x in self.busses if x.label == 'P_res']
        Q_diss = [x.P for x in self.busses if x.label == 'Q_diss']

        if len(P_res) != 0 and len(Q_diss) != 0:
            print('cycle process key figures')
            print('eta_th = ' + str(1 - sum(Q_diss) /
                  (sum(P_res) + sum(Q_diss))))

        msg = 'Do you want to print the components parammeters?'
        if hlp.query_yes_no(msg):
            self.comps.apply(network.print_components, axis=1,
                             args=(self,))

        # convert back to specified units
        for c in self.conns.index:
            c.m /= network.m[self.m_unit]
            c.p /= network.p[self.p_unit]
            c.h /= network.h[self.h_unit]
            c.T = (c.T / network.T[self.T_unit][1] -
                   network.T[self.T_unit][0])

        msg = 'Do you want to print the connections parammeters?'
        if hlp.query_yes_no(msg):
            df = pd.DataFrame(columns=['m / (' + self.m_unit + ')',
                                       'p / (' + self.p_unit + ')',
                                       'h / (' + self.h_unit + ')',
                                       'T / (' + self.T_unit + ')'])
            for c in self.conns.index:
                df.loc[c.s.label + ' -> ' + c.t.label] = (
                        ['{:.2e}'.format(c.m),
                         '{:.2e}'.format(c.p),
                         '{:.4e}'.format(c.h),
                         '{:.2e}'.format(c.T)]
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
                    axarr[i].plot(x, prop[k][:] / prop[k][-1],
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

        dec = kwargs.get('dec', '.')

        self.save_connections(filename, dec)
        self.save_components(filename, dec)
        self.save_busses(filename, dec)
        self.save_references(filename, dec)

    def save_connections(self, filename, dec):
        """
        saves the connection logic and parametrisation to filename_conn.csv

        - uses connections object id as row identifier
            * properties and property_set (False/True)
            * referenced objects
            * fluids and fluid_set vector
            * connections source and target
        - connections source and target are identified by its labels

        :param filename: suffix for the .csv-file
        :type filename: str
        :param dec: decimal separator
        :type dec: str
        :returns: no return value
        """

        # save fluids and property information of connections from dataframe
        fn = filename + '_conn.csv'

        df = pd.DataFrame()
        df['id'] = self.conns.apply(network.save_id, axis=1)

        cols = ['m', 'p', 'h', 'T', 'x',
                'm_set', 'p_set', 'h_set', 'T_set', 'x_set',
                'm_ref', 'p_ref', 'h_ref', 'T_ref']
        for val in cols:
            df[val] = self.conns.apply(network.save_props, axis=1,
                                       args=(val,))
        for val in sorted(self.fluids):
            df[val] = self.conns.apply(network.save_fluids, axis=1,
                                       args=(val,))
            df[val + '_set'] = self.conns.apply(network.save_fluids_set,
                                                axis=1, args=(val,))

        df['s'] = self.conns.apply(network.save_comps_label, axis=1,
                                   args=('s',))
        df['s_id'] = self.conns.apply(network.save_props, axis=1,
                                      args=('s_id',))

        df['t'] = self.conns.apply(network.save_comps_label, axis=1,
                                   args=('t',))
        df['t_id'] = self.conns.apply(network.save_props, axis=1,
                                      args=('t_id',))

        df.to_csv(fn, sep=';', decimal=dec, index=False, na_rep='nan')

    def save_components(self, filename, dec):
        """
        saves the components to filename_comp.csv

        - uses components labels as row identifier
        - writes:
            * components incomming and outgoing connections (object id)
            * components parametrisation

        :param filename: suffix for the .csv-file
        :type filename: str
        :returns: no return value
        """

        # create / overwrite csv file
        fn = filename + '_comp.csv'
        with open(fn, 'w') as csvfile:
            f = csv.writer(csvfile, delimiter=';')

            # write collumn heading
            f.writerow(['id', 'comp', 'i', 'o', 'mode', 'busses'])

            # write data in csv - file
            for i in range(len(self.comps)):
                c = self.comps.iloc[i]
                cp = self.comps.index[i]
                busses = []
                for b in self.busses:
                    if cp in b.comps.index.tolist():
                        busses += [str(b)[str(b).find(' at ') + 4:-1],
                                   str(b.comps.loc[cp][0]).replace('.', dec)]
                parset = []
                for var in cp.attr():
                    if (var != 'label' and var != 'mode' and
                            var != 'design' and var != 'offdesign'):
                        val = str(cp.get_attr(var)).replace('.', dec)
                        parset += [var, val]
                        if cp.get_attr(var + '_set'):
                            parset += [True]
                        else:
                            parset += [False]

                f.writerow([cp.label,
                            cp.__class__.__name__,
                            [str(x)[str(x).find(' at ') + 4:-1] for x in c.i],
                            [str(x)[str(x).find(' at ') + 4:-1] for x in c.o],
                            cp.mode,
                            busses,
                            *parset])

    def save_busses(self, filename, dec):
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
        :param dec: decimal separator
        :type dec: str
        :returns: no return value
        """

        # save fluids and property information of connections from dataframe
        fn = filename + '_bus.csv'

        df = pd.DataFrame({'id': self.busses}, index=self.busses)
        df['id'] = df.apply(network.save_id, axis=1)

        cols = ['label', 'P']
        for val in cols:
            df[val] = df.apply(network.save_props, axis=1,
                               args=(val,))

        df.to_csv(fn, sep=';', decimal=dec, index=False, na_rep='nan')

    def save_references(self, filename, dec):
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
        :param dec: decimal separator
        :type dec: str
        :returns: no return value
        """

        # save fluids and property information of connections from dataframe
        fn = filename + '_bus.csv'

        refs = []

        for c in self.conns.index:
            for attr in ['m_ref', 'p_ref', 'h_ref', 'T_ref']:
                if hasattr(c, attr):
                    refs += [c.get_attr(attr)]

        df = pd.DataFrame({'id': refs}, index=refs)
        df['id'] = df.apply(network.save_id, axis=1)

        cols = ['label', 'P']
        for val in cols:
            df[val] = df.apply(network.save_props, axis=1,
                               args=(val,))

        df.to_csv(fn, sep=';', decimal=dec, index=False, na_rep='nan')

    def save_id(c):
        return str(c.name)[str(c.name).find(' at ') + 4:-1]

    def save_props(c, *args):
        if hasattr(c.name, args[0]):
            if '_ref' in args[0]:
                ref = c.name.get_attr(args[0])
                return str(ref)[str(ref).find(' at ') + 4:-1]
            else:
                return c.name.get_attr(args[0])
        else:
            return ''

    def save_comps_label(c, *args):
        return c.name.get_attr(args[0]).label

    def save_fluids(c, *args):
        return c.name.fluid[args[0]]

    def save_fluids_set(c, *args):
        return c.name.fluid_set[args[0]]
