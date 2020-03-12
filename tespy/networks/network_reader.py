# -*- coding: utf-8

"""Module for loading a tespy network from saved state.

Use the method :func:`tespy.networks.network_reader.load_network` for importing
a network from a saved state.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/networks/network_reader.py

SPDX-License-Identifier: MIT
"""

import pandas as pd


from tespy.connections import connection, bus, ref

from tespy.components import (basics, combustion, heat_exchangers, nodes,
                              piping, reactors, turbomachinery)

from tespy.networks.networks import network

from tespy.tools.data_containers import (dc_cc, dc_cm, dc_cp, dc_flu, dc_gcp,
                                         dc_prop, dc_simple)
from tespy.tools.characteristics import char_line, char_map, compressor_map
from tespy.tools.helpers import modify_path_os
import os
import ast
import logging


global comp_target_classes
comp_target_classes = {
    'cycle_closer': basics.cycle_closer,
    'sink': basics.sink,
    'source': basics.source,
    'subsystem_interface': basics.subsystem_interface,
    'combustion_chamber': combustion.combustion_chamber,
    'combustion_chamber_stoich': combustion.combustion_chamber_stoich,
    'combustion_engine': combustion.combustion_engine,
    'condenser': heat_exchangers.condenser,
    'desuperheater': heat_exchangers.desuperheater,
    'heat_exchanger': heat_exchangers.heat_exchanger,
    'heat_exchanger_simple': heat_exchangers.heat_exchanger_simple,
    'solar_collector': heat_exchangers.solar_collector,
    'drum': nodes.drum,
    'merge': nodes.merge,
    'node': nodes.node,
    'separator': nodes.separator,
    'splitter': nodes.splitter,
    'pipe': piping.pipe,
    'valve': piping.valve,
    'water_electrolyzer': reactors.water_electrolyzer,
    'compressor': turbomachinery.compressor,
    'pump': turbomachinery.pump,
    'turbine': turbomachinery.turbine
}


global map_target_classes
map_target_classes = {
    'char_map': char_map,
    'compressor_map': compressor_map
}

# %% network loading


def load_network(path):
    r"""
    Load a network from a base path.

    Parameters
    ----------
    path : str
        The path to the network data.

    Returns
    -------
    nw : tespy.networks.networks.network
        TESPy networks object.

    Note
    ----
    If you save the network structure of an existing TESPy network, it will be
    saved to the path you specified. The structure of the saved data in that
    path is the structure you need to provide in the path for loading the
    network.

    The structure of the path must be as follows:

    - Folder: path (e. g. 'mynetwork')
    - Subfolder: comps (e. g. 'mynetwork/comps') containing

        - bus.csv*
        - char.csv*
        - char_map.csv*
        - component_class_name.csv (e. g. heat_exchanger.csv)

    - conns.csv
    - netw.csv

    The imported network has the following additional features:

    - Connections are accessible by their target's label and id, e. g. for a
      connection going to 'condenser' at inlet 'in2' use
      :code:`myimportednetwork.imp_conns['condenser:in2']`.
    - Components are accessible by label, e. g. for a component
      'heat exchanger' :code:`myimportednetwork.imp_comps['heat exchanger']`.
    - Busses are accessible by label, e. g. for a bus 'power input'
      :code:`myimportednetwork.imp_busses['power input']`.

    Example
    -------
    Create a network and export it. This is followed by loading the network
    with the network_reader module. All network information stored will be
    passed to a new network object. Components and busses will be accessible
    by label, connections by
    :code:`'source.label:source.id_target.label:target.id'`. The following
    example setup is simple gas turbine setup with compressor, combustion
    chamber and turbine.

    >>> from tespy.components import (sink, source, combustion_chamber,
    ... compressor, turbine)
    >>> from tespy.connections import connection, ref, bus
    >>> from tespy.networks import load_network, network
    >>> import shutil
    >>> fluid_list = ['CH4', 'O2', 'N2', 'CO2', 'H2O', 'Ar']
    >>> nw = network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ... h_unit='kJ / kg', T_range=[250, 1300], iterinfo=False)
    >>> air = source('air')
    >>> f = source('fuel')
    >>> c = compressor('compressor')
    >>> comb = combustion_chamber('combustion')
    >>> t = turbine('turbine')
    >>> si = sink('sink')
    >>> inc = connection(air, 'out1', c, 'in1')
    >>> cc = connection(c, 'out1', comb, 'in1')
    >>> fc = connection(f, 'out1', comb, 'in2')
    >>> ct = connection(comb, 'out1', t, 'in1')
    >>> outg = connection(t, 'out1', si, 'in1')
    >>> nw.add_conns(inc, cc, fc, ct, outg)

    Specify component and connection properties. The intlet pressure at the
    compressor and the outlet pressure after the turbine are identical. For the
    compressor, the pressure ratio and isentropic efficiency are design
    parameters. A compressor map (efficiency vs. mass flow and pressure rise
    vs. mass flow) is selected for the compressor. Fuel is Methane.

    >>> c.set_attr(pr=10, eta_s=0.88, design=['eta_s', 'pr'],
    ... offdesign=['char_map'])
    >>> t.set_attr(eta_s=0.9, design=['eta_s'],
    ... offdesign=['eta_s_char', 'cone'])
    >>> inc.set_attr(fluid={'N2': 0.7556, 'O2': 0.2315, 'Ar': 0.0129, 'CH4': 0,
    ... 'H2O': 0}, fluid_balance=True, T=25, p=1)
    >>> fc.set_attr(fluid={'N2': 0, 'O2': 0, 'Ar': 0, 'CH4': 0.96, 'H2O': 0,
    ... 'CO2': 0.04}, T=25)
    >>> ct.set_attr(T=1100)
    >>> outg.set_attr(p=ref(inc, 1, 0))

    The total power output is set to 1 MW, electrical or mechanical
    efficiencies are not considered in this example. The documentation
    example in class :func:`tespy.connections.bus` provides more information
    on efficiencies of generators, for instance.

    >>> power = bus('total power output', P=-1e6)
    >>> power.add_comps({'c': c}, {'c': t})
    >>> nw.add_busses(power)
    >>> nw.solve('design')
    >>> nw.save('exported_nwk')
    >>> c.set_attr(igva='var')
    >>> nw.solve('offdesign', design_path='exported_nwk',
    ... init_path='exported_nwk')
    >>> round(t.eta_s.val, 1)
    0.9
    >>> power.set_attr(P=-0.75e6)
    >>> nw.solve('offdesign', design_path='exported_nwk',
    ... init_path='exported_nwk')
    >>> eta_s_t = round(t.eta_s.val, 3)
    >>> igva = round(c.igva.val, 3)
    >>> eta_s_t
    0.898
    >>> igva
    20.139

    The designed network is exported to the path 'exported_nwk'. Now import the
    network and recalculate. Check if the results match with the previous
    calculation in design and offdesign case.

    >>> imported_nwk = load_network('exported_nwk')
    >>> imported_nwk.set_attr(iterinfo=False)
    >>> imported_nwk.solve('design')
    >>> round(imported_nwk.imp_comps['turbine'].eta_s.val, 3)
    0.9
    >>> imported_nwk.imp_comps['compressor'].set_attr(igva='var')
    >>> imported_nwk.solve('offdesign', design_path='exported_nwk',
    ... init_path='exported_nwk')
    >>> round(imported_nwk.imp_comps['turbine'].eta_s.val, 3)
    0.9
    >>> imported_nwk.imp_busses['total power output'].set_attr(P=-0.75e6)
    >>> imported_nwk.solve('offdesign', design_path='exported_nwk',
    ... init_path='exported_nwk')
    >>> round(imported_nwk.imp_comps['turbine'].eta_s.val, 3) == eta_s_t
    True
    >>> round(imported_nwk.imp_comps['compressor'].igva.val, 3) == igva
    True
    >>> shutil.rmtree('./exported_nwk', ignore_errors=True)
    """
    if path[-1] != '/' and path[-1] != '\\':
        path += '/'

    path_comps = modify_path_os(path + 'comps/')
    path = modify_path_os(path)

    msg = 'Reading network data from base path ' + path + '.'
    logging.info(msg)

    # load characteristics
    fn = path_comps + 'char_line.csv'
    try:
        char_lines = pd.read_csv(fn, sep=';', decimal='.',
                                 converters={'x': ast.literal_eval,
                                             'y': ast.literal_eval})
        msg = 'Reading characteristic lines data from ' + fn + '.'
        logging.debug(msg)

    except FileNotFoundError:
        char_lines = pd.DataFrame(columns=['id', 'type', 'x', 'y'])

    # load characteristic maps
    fn = path_comps + 'char_map.csv'
    try:
        msg = 'Reading characteristic maps data from ' + fn + '.'
        logging.debug(msg)
        char_maps = pd.read_csv(fn, sep=';', decimal='.',
                                converters={'x': ast.literal_eval,
                                            'y': ast.literal_eval,
                                            'z1': ast.literal_eval,
                                            'z2': ast.literal_eval})

    except FileNotFoundError:
        char_maps = pd.DataFrame(columns=['id', 'type', 'x', 'y', 'z1', 'z2'])

    # load components
    comps = pd.DataFrame()

    files = os.listdir(path_comps)
    for f in files:
        if f != 'bus.csv' and f != 'char_line.csv' and f != 'char_map.csv':
            fn = path_comps + f
            df = pd.read_csv(fn, sep=';', decimal='.',
                             converters={'design': ast.literal_eval,
                                         'offdesign': ast.literal_eval,
                                         'busses': ast.literal_eval,
                                         'bus_param': ast.literal_eval,
                                         'bus_P_ref': ast.literal_eval,
                                         'bus_char': ast.literal_eval})

            # create components
            df['instance'] = df.apply(construct_comps, axis=1,
                                      args=(char_lines, char_maps, ))
            comps = pd.concat((comps, df[['instance', 'label', 'busses',
                                          'bus_param', 'bus_P_ref',
                                          'bus_char']]), axis=0)

            msg = 'Reading component data (' + f[:-4] + ') from ' + fn + '.'
            logging.debug(msg)

    comps = comps.set_index('label')
    msg = 'Created network components.'
    logging.info(msg)

    # create network
    nw = construct_network(path)

    # make components accessible by labels
    nw.imp_comps = comps.to_dict()['instance']

    # load connections
    fn = path + 'conn.csv'
    conns = pd.read_csv(fn, sep=';', decimal='.',
                        converters={'design': ast.literal_eval,
                                    'offdesign': ast.literal_eval})

    msg = 'Reading connection data from ' + fn + '.'
    logging.debug(msg)

    # create connections
    conns['instance'] = conns.apply(construct_conns, axis=1, args=(comps, nw,))
    conns.apply(conns_set_ref, axis=1, args=(conns,))
    conns = conns.set_index('id')

    nw.imp_conns = {}
    # add connections to network
    for c in conns['instance']:
        nw.add_conns(c)
        nw.imp_conns[c.s.label + ':' + c.s_id + '_'
                     + c.t.label + ':' + c.t_id] = c

    msg = 'Created connections.'
    logging.info(msg)

    # load busses
    try:
        fn = path_comps + 'bus.csv'
        busses = pd.read_csv(fn, sep=';', decimal='.')
        msg = 'Reading bus data from ' + fn + '.'
        logging.debug(msg)

    except FileNotFoundError:
        busses = pd.DataFrame()
        msg = 'No bus data found!'
        logging.debug(msg)

    # create busses
    nw.imp_busses = {}
    if len(busses) > 0:
        busses['instance'] = busses.apply(construct_busses, axis=1)

        # add components to busses
        comps.apply(busses_add_comps, axis=1, args=(busses, char_lines,))

        # add busses to network
        for b in busses['instance']:
            nw.add_busses(b)
            nw.imp_busses[b.label] = b

        msg = 'Created busses.'
        logging.info(msg)

    msg = 'Created network.'
    logging.info(msg)
    return nw


# %% create components


def construct_comps(c, *args):
    r"""
    Create TESPy component from class name and set parameters.

    Parameters
    ----------
    c : pandas.core.series.Series
        Component information from .csv-file.

    args[0] : pandas.core.frame.DataFrame
        DataFrame containing the data of characteristic lines.

    args[1] : pandas.core.frame.DataFrame
        DataFrame containing the data of characteristic maps.

    Returns
    -------
    instance : tespy.components.components.component
        TESPy component object.
    """
    target_class = comp_target_classes[c.cp]
    instance = target_class(c.label)
    kwargs = {}

    # basic properties
    for key in ['design', 'offdesign', 'design_path', 'local_design',
                'local_offdesign']:
        if key in c:
            kwargs[key] = c[key]

    for key, value in instance.attr().items():
        if key in c:
            # component parameters
            if isinstance(value, dc_cp):
                kwargs[key] = dc_cp(val=c[key], is_set=c[key + '_set'],
                                    is_var=c[key + '_var'])

            # component parameters
            elif isinstance(value, dc_simple):
                kwargs[key] = dc_simple(val=c[key], is_set=c[key + '_set'])

            # component characteristics
            elif isinstance(value, dc_cc):
                # finding x and y values of the characteristic function
                values = args[0]['id'] == c[key]

                try:
                    x = args[0][values].x.values[0]
                    y = args[0][values].y.values[0]
                    if 'extrapolate' in args[0].columns:
                        extrapolate = args[0][values].extrapolate.values[0]
                    else:
                        extrapolate = False
                    char = char_line(x=x, y=y, extrapolate=extrapolate)

                except IndexError:

                    char = None
                    msg = ('Could not find x and y values for characteristic '
                           'line, using defaults instead for function ' + key +
                           ' at component ' + c.label + '.')
                    logging.warning(msg)

                kwargs[key] = dc_cc(is_set=c[key + '_set'],
                                    param=c[key + '_param'], func=char)

            # component characteristics
            elif isinstance(value, dc_cm):
                # finding x and y values of the characteristic function
                values = args[1]['id'] == c[key]

                try:
                    x = list(args[1][values].x.values[0])
                    y = list(args[1][values].y.values[0])
                    z1 = list(args[1][values].z1.values[0])
                    z2 = list(args[1][values].z2.values[0])
                    target_class = map_target_classes[
                        args[1][values].type.values[0]]
                    char = target_class(x=x, y=y, z1=z1, z2=z2)

                except IndexError:
                    char = None
                    msg = ('Could not find x, y, z1 and z2 values for '
                           'characteristic map of component ' + c.label + '!')
                    logging.warning(msg)

                kwargs[key] = dc_cm(is_set=c[key + '_set'],
                                    param=c[key + '_param'],
                                    func=char)

            # grouped component parameters
            elif isinstance(value, dc_gcp):
                kwargs[key] = dc_gcp(method=c[key])

    instance.set_attr(**kwargs)
    return instance

# %% create network object


def construct_network(path):
    r"""
    Create TESPy network from the data provided in the netw.csv-file.

    Parameters
    ----------
    path : str
        Base-path to stored network data.

    Returns
    -------
    nw : tespy.networks.networks.network
        TESPy network object.
    """
    # read network .csv-file
    netw = pd.read_csv(path + 'netw.csv', sep=';', decimal='.',
                       converters={'fluids': ast.literal_eval})
    f_list = netw['fluids'][0]

    kwargs = {}
    kwargs['m_unit'] = netw['m_unit'][0]
    kwargs['p_unit'] = netw['p_unit'][0]
    kwargs['p_range'] = [netw['p_min'][0], netw['p_max'][0]]
    kwargs['h_unit'] = netw['h_unit'][0]
    kwargs['h_range'] = [netw['h_min'][0], netw['h_max'][0]]
    kwargs['T_unit'] = netw['T_unit'][0]
    kwargs['T_range'] = [netw['T_min'][0], netw['T_max'][0]]

    # create network object with its properties
    nw = network(fluids=f_list, **kwargs)

    return nw

# %% create connections


def construct_conns(c, *args):
    r"""
    Create TESPy connection from data in the .csv-file and its parameters.

    Parameters
    ----------
    c : pandas.core.series.Series
        Connection information from .csv-file.

    args[0] : pandas.core.frame.DataFrame
        DataFrame containing all created components.

    Returns
    -------
    conn : tespy.connections.connection
        TESPy connection object.
    """
    # create connection
    conn = connection(args[0].instance[c.s], c.s_id,
                      args[0].instance[c.t], c.t_id)

    kwargs = {}
    # read basic properties
    for key in ['design', 'offdesign', 'design_path', 'local_design',
                'local_offdesign']:
        if key in c:
            kwargs[key] = c[key]

    # read fluid properties
    for key in ['m', 'p', 'h', 'T', 'x', 'v', 'Td_bp']:
        if key in c:
            if key in c:
                kwargs[key] = dc_prop(val=c[key], val0=c[key + '0'],
                                      val_set=c[key + '_set'],
                                      unit=c[key + '_unit'],
                                      unit_set=c[key + '_unit_set'],
                                      ref=None, ref_set=c[key + '_ref_set'])

    key = 'state'
    if key in c:
        kwargs[key] = dc_simple(val=c[key], is_set=c[key + '_set'])

    # read fluid vector
    val = {}
    val0 = {}
    val_set = {}
    for key in args[1].fluids:
        if key in c:
            val[key] = c[key]
            val0[key] = c[key + '0']
            val_set[key] = c[key + '_set']

    kwargs['fluid'] = dc_flu(val=val, val0=val0, val_set=val_set,
                             balance=c['balance'])

    # write properties to connection and return connection object
    conn.set_attr(**kwargs)
    return conn

# %% set references on connections


def conns_set_ref(c, *args):
    r"""
    Set references on connections as specified in connection data.

    Parameters
    ----------
    c : pandas.core.series.Series
        Connection information from .csv-file.

    args[0] : pandas.core.frame.DataFrame
        DataFrame containing all created connections.

    Returns
    -------
    instance : tespy.connections.ref
        TESPy reference object.
    """
    for col in ['m', 'p', 'h', 'T']:
        # search for referenced connections
        if isinstance(c[col + '_ref'], str):
            # create reference object
            instance = args[0].instance[c[col + '_ref'] ==
                                        args[0]['id']].values[0]
            # write to connection properties
            c['instance'].get_attr(col).ref = ref(instance,
                                                  c[col + '_ref_f'],
                                                  c[col + '_ref_d'])

# %% create busses


def construct_busses(c, *args):
    r"""
    Create busses of the network.

    Parameters
    ----------
    c : pandas.core.series.Series
        Bus information from .csv-file.

    Returns
    -------
    b : tespy.connections.bus
        TESPy bus object.
    """
    # set up bus with label and specify value for power
    b = bus(c.label, P=c.P)
    b.P.is_set = c.P_set
    return b

# %% add components to busses


def busses_add_comps(c, *args):
    r"""
    Add components to busses according to data from .csv file.

    Parameters
    ----------
    c : pandas.core.series.Series
        Component information from .csv-file.

    args[0] : pandas.core.frame.DataFrame
        DataFrame containing all created busses.

    args[1] : pandas.core.frame.DataFrame
        DataFrame containing all created characteristic lines.
    """
    i = 0
    for b in c.busses:
        p, P_ref, char = c.bus_param[i], c.bus_P_ref[i], c.bus_char[i]

        values = char == args[1]['id']
        char = char_line(x=args[1][values].x.values[0],
                         y=args[1][values].y.values[0])

        # add component with corresponding details to bus
        args[0].instance[b == args[0]['label']].values[0].add_comps(
                {'c': c.instance, 'p': p, 'P_ref': P_ref, 'char': char})
        i += 1
