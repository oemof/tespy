"""
.. module:: networks
    :synopsis: contains logic of the network

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>
"""

import math

import csv

import pandas as pd

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
    """
    class network aggregates information on components, connections and busses
    """

    # unit systems, calculation is alsways performed with SI-units
    m_unit = {
        'kg / s': 1,
        't / h': 3.6
    }
    p_unit = {
        'Pa': 1,
        'psi': 6.8948e3,
        'bar': 1e5,
        'MPa': 1e6
    }
    h_unit = {
        'J / kg': 1,
        'kJ / kg': 1e3,
        'MJ / kg': 1e6
    }
    T_unit = {
        'C': [273.15, 1],
        'F': [459.67, 5 / 9],
        'K': [0, 1]
    }

    def __init__(self, **kwargs):
        self.conns = pd.DataFrame(columns=['s', 's_id', 't', 't_id'])
        self.fluids = []

        self.convergence = np.array([0, 0, 0], dtype=object)
        self.busses = []

        # standard unit set
        self.m = 'kg / s'
        self.p = 'bar'
        self.h = 'J / kg'
        self.T = 'K'

        # standard value range
        self.p_range = [0.02, 300]
        self.h_range = [1e3, 7e6]
        self.T_range = [273.15, 1773.15]

        for key in kwargs:
            if key in self.var():
                self.__dict__.update({key: kwargs[key]})

        if not isinstance(self.fluids, list):
            raise hlp.MyNetworkError('Fluids must be stated as list.')
        else:
            self.fluids = sorted(self.fluids)

        if self.m not in network.m_unit.keys():
            msg = ('Allowed units for mass flow are: ' +
                   str(network.m_unit.keys()))
            raise hlp.MyNetworkError(msg)

        if self.p not in network.p_unit.keys():
            msg = ('Allowed units for pressure are: ' +
                   str(network.p_unit.keys()))
            raise hlp.MyNetworkError(msg)

        if self.h not in network.h_unit.keys():
            msg = ('Allowed units for enthalpy are: ' +
                   str(network.h_unit.keys()))
            raise hlp.MyNetworkError(msg)

        if self.T not in network.T_unit.keys():
            msg = ('Allowed units for temperature are: ' +
                   str(network.T_unit.keys()))
            raise hlp.MyNetworkError(msg)

        if not isinstance(self.p_range, list):
            msg = ('Specify the value range as list: [p_min, p_max]')
            raise hlp.MyNetworkError(msg)
        else:
            self.p_range[0] *= network.p_unit[self.p]
            self.p_range[1] *= network.p_unit[self.p]

        if not isinstance(self.h_range, list):
            msg = ('Specify the value range as list: [h_min, h_max]')
            raise hlp.MyNetworkError(msg)
        else:
            self.h_range[0] *= network.h_unit[self.h]
            self.h_range[1] *= network.h_unit[self.h]

        if not isinstance(self.T_range, list):
            msg = ('Specify the value range as list: [T_min, T_max]')
            raise hlp.MyNetworkError(msg)
        else:
            self.T_range[0] = ((self.T_range[0] + network.T_unit[self.T][0]) *
                               network.T_unit[self.T][1])
            self.T_range[1] = ((self.T_range[1] + network.T_unit[self.T][0]) *
                               network.T_unit[self.T][1])

        for f in self.fluids:
            hlp.molar_masses[f] = CPPSI('M', f)
            hlp.gas_constants[f] = CPPSI('GAS_CONSTANT', f)

        hlp.memorise(len(self.fluids))

    def var(self):
        return ['m', 'p', 'h', 'T', 'p_range', 'h_range', 'T_range',  'fluids']

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

    def del_conns(self, c):
        """
        delets connections from a network

        :param c: connections object to delete
        :type c: tespy.connection
        :returns: no return value
        :raises: :code:`KeyError` if connections object c is not in the network
        """
        self.conns.drop(self.conns.index(c))

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

        print('Have you adjusted the value ranges for pressure, enthalpy and '
              'temperature according to the specified unit system?')

        if len(self.fluids) == 0:
            msg = ('Network has no fluids, please specify a list with fluids '
                   'on network creation.')
            raise hlp.MyNetworkError(msg)
        self.init_components()  # build the dataframe for components
        self.init_fluids()  # start standard fluid initialisation
        self.init_properties()  # start standard property initialisation

        if self.mode == 'offdesign' and self.design_file is None:
            msg = ('Please provide \'design_file\' for every offdesign '
                   'calculation.')
            raise hlp.MyNetworkError(msg)  # must provide design_file
        else:
            self.init_csv()  # initialisation from csv

        if self.mode == 'offdesign':
            self.init_offdesign()  # characteristics for offdesign

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
                    c.fluid[fluid] = 0
                    c.fluid_set[fluid] = False

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
            inconn = [x for x in self.comps.loc[c.s].o if
                      x in self.comps.loc[c.t].i]
            inconn_id = self.comps.loc[c.t].i.tolist().index(inconn[0])
            outconn = self.comps.loc[c.t].o.tolist()[inconn_id]
            for fluid, x in c.fluid.items():
                if not outconn.fluid_set[fluid]:
                    outconn.fluid[fluid] = x

            self.init_target(outconn, start)

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
            outconn = [x for x in self.comps.loc[c.t].i if
                       x in self.comps.loc[c.s].o]
            outconn_id = self.comps.loc[c.s].o.tolist().index(outconn[0])
            inconn = self.comps.loc[c.s].i.tolist()[outconn_id]
            for fluid, x in c.fluid.items():
                if not inconn.fluid_set[fluid]:
                    inconn.fluid[fluid] = x

            self.init_source(inconn, start)

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
                c.m *= network.m_unit[self.m]
            if not c.m_set:
                c.m = c.m0

            self.init_p0(c)
            if not isinstance(c.p, con.ref):
                c.p *= network.p_unit[self.p]
            if not c.p_set:
                c.p = c.p0

            self.init_h0(c)
            if not isinstance(c.h, con.ref):
                c.h *= network.h_unit[self.h]
            if not c.h_set:
                c.h = c.h0

        self.init_refobj()

        # propate fluids towards merges targets
        # add merge cp to list redo, if number of fluids at merges outlet is
        # still zero
        redo = []
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
                c.T = ((c.T + network.T_unit[self.T][0]) *
                       network.T_unit[self.T][1])
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
            c.p0 *= network.p_unit[self.p]

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
            c.h0 *= network.h_unit[self.h]

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

            df = pd.read_csv(self.design_file, index_col=0, delimiter=';')
            for c in self.conns.index:
                # match connection (source, source_id, target, target_id) on connection objects
                # of design file
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
                for conn in df_tmp.index[s & s_id & t & t_id]:
                    c.m = df.loc[conn].m * network.m_unit[self.m]
                    c.p = df.loc[conn].p * network.p_unit[self.p]
                    c.h = df.loc[conn].h * network.h_unit[self.h]
                    for fluid in self.fluids:
                        c.fluid[fluid] = df.loc[conn][fluid]

            # component characteristics creation for offdesign calculation
            self.process_components('pre')

            for c in self.conns.index:
                c.m = c.m_tmp
                c.p = c.p_tmp
                c.h = c.h_tmp
                c.fluid = c.fluid_tmp

        if self.init_file is not None:

            df = pd.read_csv(self.init_file, index_col=0, delimiter=';')
            for c in self.conns.index:
                df_tmp = (df.s == c.s.label).to_frame()
                df_tmp.loc[:, 's_id'] = (df.s_id == c.s_id)
                df_tmp.loc[:, 't'] = (df.t == c.t.label)
                df_tmp.loc[:, 't_id'] = (df.t_id == c.t_id)
                s = df_tmp['s'] == True
                s_id = df_tmp['s_id'] == True
                t = df_tmp['t'] == True
                t_id = df_tmp['t_id'] == True
                for conn in df_tmp.index[s & s_id & t & t_id]:
                    if not c.m_set:
                        c.m = df.loc[conn].m * network.m_unit[self.m]
                    if not c.p_set:
                        c.p = df.loc[conn].p * network.p_unit[self.p]
                    if not c.h_set:
                        c.h = df.loc[conn].h * network.h_unit[self.h]
                    for fluid in self.fluids:
                        if not c.fluid_set[fluid]:
                            c.fluid[fluid] = df.loc[conn][fluid]

            for c in self.conns.index:
                if c.T_set and not isinstance(c.T, con.ref):
                    c.h = hlp.h_mix_pT(c.as_list(), c.T)
                if c.x_set:
                    c.h = hlp.h_mix_pQ(c.as_list(), c.x)

    def init_offdesign(self):
        """
        auto switches components and connections from design to offdesign mode.

        **components**

        If :code:`cp.mode == 'auto'` all parameters stated the components
        method :code:`cp.design()` will be unset and all parameters stated in
        the components method :code:`cp.offdesign()` will be set instead.

        The auto-switch can be deactivated by using
        :code:`your_component.set_attr(mode='man')`

        **connections**

        All parameters given in the connections attribute :code:`c.design`
        will be unset.

        :returns: no return value
        """
        for cp in self.comps.index:
            if cp.mode == 'auto':
                for var in cp.design():
                    if cp.__dict__[var + '_set']:
                        cp.__dict__[var + '_set'] = False
                for var in cp.offdesign():
                    if not cp.__dict__[var + '_set']:
                        cp.__dict__[var + '_set'] = True

        for c in self.conns.index:
            for var in c.design:
                if c.__dict__[var + '_set']:
                    c.__dict__[var + '_set'] = False

    def solve(self, mode, init_file=None, design_file=None):
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

        if mode != 'offdesign' and mode != 'design':
            msg = 'Mode must be \'design\' or \'offdesign\'.'
            raise hlp.MyNetworkError(msg)
        else:
            self.mode = mode

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

        print('iter\t| residual')
        for self.iter in range(250):
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

            if self.iter > 10:
                if all(self.res[(self.iter - 5):] >= self.res[-4]):
                    print('Convergence is making no progress, '
                          'calculation stopped.')
                    break
        end_time = time.time()

        for c in self.conns.index:
            c.T = (hlp.T_mix_ph([c.m, c.p, c.h, c.fluid]) /
                   network.T_unit[self.T][1] - network.T_unit[self.T][0])
            c.m /= network.m_unit[self.m]
            c.p /= network.p_unit[self.p]
            c.h /= network.h_unit[self.h]
            c.m0 = c.m
            c.p0 = c.p
            c.h0 = c.h

        print('Calculation complete.')
        print('Total iterations:', self.iter, ' -  '
              'Calculation time:', round(end_time - start_time, 1), 's - '
              'Iterations per second:',
              round(self.iter / (end_time - start_time), 2))

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
            if self.num_vars * len(self.conns.index) > len(self.vec_res):
                msg = ('You have not provided enough parameters:',
                       self.num_vars * len(self.conns.index), 'required, ',
                       len(self.vec_res), 'given.')
                raise hlp.MyNetworkError(msg)
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

            # prevent bad changes within solution process
            if c.p <= self.p_range[0]:
                c.p = self.p_range[0]
            if c.p >= self.p_range[1]:
                c.p = self.p_range[1]
            if c.h < self.h_range[0]:
                c.h = self.h_range[0]
            if c.h > self.h_range[1]:
                c.h = self.h_range[1]

            # make sure, that at given temperatures values stay within feasible
            # enthalpy range: calculate maximum enthalpy and compare with
            # acutal value
            if self.iter < 5 and c.T_set:
                h_max = hlp.h_mix_pT(c.as_list(), self.T_range[1])
                if c.h > h_max:
                    c.h = h_max * 0.98
                if hasattr(c, 'T_ref'):
                    c = c.T_ref.obj
                    h_max = hlp.h_mix_pT(c.as_list(), self.T_range[1])
                    if c.h > h_max:
                        c.h = h_max * 0.98

                h_min = hlp.h_mix_pT(c.as_list(), self.T_range[0])
                if c.h < h_min:
                    c.h = h_min * 1.02
                if hasattr(c, 'T_ref'):
                    c = c.T_ref.obj
                    h_min = hlp.h_mix_pT(c.as_list(), self.T_range[0])
                    if c.h < h_min:
                        c.h = h_min * 1.02

        # check properties for consistency
        if self.init_file is None and self.iter < 5:
            for cp in self.comps.index:
                cp.convergence_check(self)

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
        num_eq = 0
        for cp in self.comps.index.values:
            self.vec_res += cp.equations(self)
            vec_deriv = cp.derivatives(self)
            if not isinstance(cp, cmp.source) and not isinstance(cp, cmp.sink):
                i = 0
                inlets = self.comps.loc[cp].i.tolist()
                outlets = self.comps.loc[cp].o.tolist()
                self.debug += [cp] * (len(vec_deriv))
                # insert derivatives row - wise
                # loc is the location of c in the jacobian matrix
                for c in inlets + outlets:
                    loc = self.conns.index.get_loc(c)
                    self.mat_deriv[num_eq:len(self.vec_res),
                                   loc * (self.num_vars):
                                   (loc + 1) * self.num_vars] = (
                                       vec_deriv[:, i])
                    i += 1

            num_eq = len(self.vec_res)

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
                row = self.solve_check_row(row)
            if c.p_set:
                self.solve_m_p_h(c, row, col, 'p')
                row = self.solve_check_row(row)
            if c.h_set:
                self.solve_m_p_h(c, row, col, 'h')
                row = self.solve_check_row(row)
            if c.T_set:
                self.solve_T(c, row, col)
                row = self.solve_check_row(row)
            if c.x_set:
                self.solve_x(c, row, col)
                row = self.solve_check_row(row)

            l = 0
            for f in c.fluid.keys():
                if c.fluid_set[f]:
                    self.mat_deriv[row, col + 3 + l] = 1
                    self.vec_res += [0]
                    row = self.solve_check_row(row)
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
                row = self.solve_check_row(row)
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

                row = self.solve_check_row(row)

    def solve_check_row(self, row):
        """
        check row index increment for jacobian matrix

        :param row: index of last row inserted into jacobian matrix
        :type row: int
        :returns: row + 1 (int) - row increment
        :raises: :code:`hlp.MyNetworkError`, if the row count is higher than
                 the number of total variables in the network (overdetermined).
        """
        if row < self.num_vars * len(self.conns.index):
            return row + 1
        else:
            raise hlp.MyNetworkError('You provide too many parameters!')

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

    def process_components(self, mode):
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

        if mode == 'pre':
            self.comps.apply(network.process_pre, axis=1)
            print('Preprocessing done.')

        if mode == 'post':
            # convert to SI-units
            for c in self.conns.index:
                c.m *= network.m_unit[self.m]
                c.p *= network.p_unit[self.p]
                c.h *= network.h_unit[self.h]
                c.T = ((c.T + network.T_unit[self.T][0]) *
                       network.T_unit[self.T][1])

# process key figures
            for b in self.busses:
                b.P = 0
                for c in b.comps.index:
                    i = self.comps.loc[c].i[0]
                    o = self.comps.loc[c].o[0]

                    factor = b.comps.loc[c].factor

                    b.P += i.m * (o.h - i.h) * factor

            P = [x.P for x in self.busses if x.label == 'P_res']
            Q_diss = [x.P for x in self.busses if x.label == 'Q_diss']

            print('Postprocessing.')

            print('##### Kennzahlen des Kreisprozesses #####')
            if len(P) != 0 and len(Q_diss) != 0:
                print('eta_th = ', 1 - Q_diss / (P + Q_diss))
            msg = 'Do you want to print the components parammeters?'
            if hlp.query_yes_no(msg):
                self.comps.apply(network.process_post_print, axis=1)
            else:
                self.comps.apply(network.process_post_calc, axis=1)

            # convert back to specified units
            for c in self.conns.index:
                c.m /= network.m_unit[self.m]
                c.p /= network.p_unit[self.p]
                c.h /= network.h_unit[self.h]
                c.T = (c.T / network.T_unit[self.T][1] -
                       network.T_unit[self.T][0])

    def process_post_calc(cols):
        """
        postprocessing: calculate components attributes

        :param cols: cols are the components of the network
        :type cols: landas dataframe index object
        :returns: no return value
        """
        inlets, outlets = cols.i.tolist(), cols.o.tolist()
        cols.name.calc_parameters(inlets, outlets, 'post')

    def process_post_print(cols):
        """
        postprocessing: calculate components attributes and print them to
        prompt

        :param cols: cols are the components of the network
        :type cols: landas dataframe index object
        :returns: no return value
        """
        inlets, outlets = cols.i.tolist(), cols.o.tolist()
        cols.name.calc_parameters(inlets, outlets, 'post')
        cols.name.print_parameters(inlets, outlets)

    def process_pre(cols):
        """
        preprocessing: calculate components attributes

        :param cols: cols are the components of the network
        :type cols: landas dataframe index object
        :returns: no return value
        """
        inlets, outlets = cols.i.tolist(), cols.o.tolist()
        cols.name.calc_parameters(inlets, outlets, 'pre')

    def plot_convergence(self, mode):
        """
        plots the convergence history of all mass flows, pressures and
        enthalpies as absolute values or relative to last value

        :param mode: absolute values or relative to last value (abs, rel)
        :type mode: str
        :returns: no return value
        :raises: :code:`ValueError` if mode is neither 'abs' nor 'rel'
        """
        if mode not in ['abs', 'rel']:
            raise ValueError(
                'Plotting mode must be either \'abs\' or \'rel\'.')

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
        if mode == 'rel':
            for c in self.conns.index:
                i = 0
                for prop in self.convergence:
                        axarr[i].semilogy(x, prop[k][:] / prop[k][-1],
                                          color=color[k],
                                          label=c.s.label + ' -> ' + c.t.label)
                        i += 1
                k += 1

        else:
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

    def save(self, filename):
        """
        saves the results in two files:

        - results file and
        - components file

        :param filename: suffix for the .csv-file
        :type filename: str
        :returns: no return value
        """
        self.save_connections(filename)
        self.save_components(filename)

    def save_connections(self, filename):
        """
        saves the results to filename_results.csv

        - uses connections object id as row identifier
        - writes:
            * properties (value and property_set: False/True)
            * fluid vector
            * fluid_set vector
            * connections source and target
            * connections parametrisation
        - connections source and target identified by its labels

        :param filename: suffix for the .csv-file
        :type filename: str
        :returns: no return value
        """

# create / overwrite csv file
        fn = filename + '_results.csv'
        with open(fn, 'w') as csvfile:
            f = csv.writer(csvfile, delimiter=';')

# write collumn headings
            cols = (['id', 'm', 'p', 'h', 'T', 'x'] + sorted(self.fluids) +
                    ['m_set', 'p_set', 'h_set', 'T_set', 'x_set'] +
                    ['{}_{}'.format(fluid, 'set')
                        for fluid in sorted(self.fluids)] +
                    ['s', 's_id', 't', 't_id'])
            f.writerow(cols)

# write data in csv - file
            for i in range(len(self.conns)):
                c = self.conns.index[i]
                cp = self.conns.iloc[i]
                f.writerow([str(c)[str(c).find(' at ') + 4:-1],
                            c.m, c.p, c.h, c.T, c.x] +
                           [x for fluid, x in sorted(c.fluid.items())] +
                           [c.m_set, c.p_set, c.h_set, c.T_set, c.x_set] +
                           [fluid_set for fluid_set in c.fluid_set.values()] +
                           [cp.s.label, cp.s_id, cp.t.label, cp.t_id])

    def save_components(self, filename):
        """
        saves the components to filename_comps.csv

        - uses components labels as row identifier
        - writes:
            * components incomming and outgoing connections (object id)
            * components parametrisation

        :param filename: suffix for the .csv-file
        :type filename: str
        :returns: no return value
        """

# create / overwrite csv file
        fn = filename + '_comps.csv'
        with open(fn, 'w') as csvfile:
            f = csv.writer(csvfile, delimiter=';')

# write collumn heading
            f.writerow(['id', 'comp', 'i', 'o', 'mode'])

# write data in csv - file
            for i in range(len(self.comps)):
                c = self.comps.iloc[i]
                cp = self.comps.index[i]
                parset = []
                for var in cp.attr():
                    if var != 'label' and var != 'mode':
                        if cp.get_attr(var + '_set'):
                            parset += [var, cp.get_attr(var)]

                f.writerow([cp.label,
                            cp.__class__.__name__,
                            [str(x)[str(x).find(' at ') + 4:-1] for x in c.i],
                            [str(x)[str(x).find(' at ') + 4:-1] for x in c.o],
                            cp.mode,
                            *parset])
