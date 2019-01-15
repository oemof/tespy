"""
.. module:: network_reader
    :synopsis:

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>
"""

import pandas as pd
import numpy as np
from tespy import cmp, con, nwk, hlp, cmp_char
import os
import ast

# %% network loading


def load_nwk(path):
    r"""
    Loads a network from a base path.

    Parameters
    ----------
    path : String
        The path to the network data.

    Returns
    -------
    nw : tespy.networks.network
        TESPy networks object.

    Note
    ----
    If you save the network structure of an existing TESPy network, it will be saved to the path you specified.
    The structure of the saved data in that path is the structure you need to provide in the path for loading the network.

    The structure of the path must be as follows:

    - Folder: path (e. g. 'mynetwork')
    - Subfolder: comps (e. g. 'mynetwork/comps') containing
        - bus.csv*
        - char.csv*
        - component_class_name.csv (e. g. heat_exchanger.csv)
    - conns.csv
    - netw.csv

    Example
    -------
    >>> from tespy import cmp, con, nwk, hlp, nwkr
    >>> fluid_list = ['water']
    >>> nw = nwk.network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ...     h_unit='kJ / kg')
    >>> nw.set_printoptions(print_level='err')
    >>> si = cmp.sink('sink')
    >>> so = cmp.source('source')
    >>> t = cmp.turbine('turbine')
    >>> inc = con.connection(so, 'out1', t, 'in1')
    >>> outg = con.connection(t, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)
    >>> t.set_attr(pr=0.02, eta_s=0.8, P=-1e5, design=['eta_s', 'pr'],
    ...     offdesign=['eta_s_char', 'cone'])
    >>> inc.set_attr(fluid={'water': 1}, T=600)
    >>> outg.set_attr(p=0.5)
    >>> nw.solve('design')
    >>> nw.save('tmp', structure=True)
    >>> t.set_attr(P=-9e4)
    >>> nw.solve('offdesign', design_file='tmp/results.csv')
    >>> round(t.eta_s.val, 3)
    0.798
    >>> nw2 = nwkr.load_nwk('tmp')
    Reading network data.
    Created components.
    Created connections.
    >>> nw2.set_printoptions(print_level='err')
    >>> nw2.solve('design')
    >>> nw2.solve('offdesign', design_file='tmp/results.csv')
    """
    print('Reading network data.')

    # load characteristics
    chars = pd.read_csv(path + '/comps/char.csv', sep=';', decimal='.',
                        converters={'x': ast.literal_eval,
                                    'y': ast.literal_eval})

    # load characteristic maps
    char_maps = pd.read_csv(path + '/comps/char_map.csv', sep=';', decimal='.',
                            converters={'x': ast.literal_eval,
                                        'y': ast.literal_eval,
                                        'z1': ast.literal_eval,
                                        'z2': ast.literal_eval})

    # load components
    comps = pd.DataFrame()
    inter = pd.DataFrame()

    files = os.listdir(path + '/comps/')
    for f in files:
        if f != 'bus.csv' and f != 'char.csv' and f != 'char_map.csv':
            df = pd.read_csv(path + '/comps/' + f, sep=';', decimal='.',
                             converters={'design': ast.literal_eval,
                                         'offdesign': ast.literal_eval,
                                         'busses': ast.literal_eval,
                                         'bus_param': ast.literal_eval,
                                         'bus_P_ref': ast.literal_eval,
                                         'bus_char': ast.literal_eval})

            # create components
            df['instance'] = df.apply(construct_comps, axis=1, args=(chars, char_maps, ))
            comps = pd.concat((comps, df[['instance', 'label', 'busses',
                                          'bus_param', 'bus_P_ref',
                                          'bus_char']]), axis=0)

            df['inter'] = df.apply(get_interface, axis=1)
            inter = pd.concat((inter, df[['instance', 'label', 'inter']]), axis=0)

    comps = comps.set_index('label')
    print('Created components.')

    # create network
    nw = construct_network(path)
    inter = inter[inter['inter'] == True].drop('inter', axis=1)

    # make interfaces and components accessible by labels
    nw.imp_comps = comps.to_dict()['instance']
    nw.inter = inter.set_index('label').to_dict()['instance']

    # load connections
    conns = pd.read_csv(path + '/conn.csv', sep=';', decimal='.',
                        converters={'design': ast.literal_eval,
                                    'offdesign': ast.literal_eval})

    # create connections
    conns['instance'] = conns.apply(construct_conns, axis=1, args=(comps, nw,))
    conns.apply(conns_set_ref, axis=1, args=(conns,))
    conns = conns.set_index('id')

    nw.imp_conns = {}
    # add connections to network
    for c in conns['instance']:
        nw.add_conns(c)
        nw.imp_conns[c.t.label + ':' + c.t_id] = c

    print('Created connections.')

    # load busses
    busses = pd.DataFrame()
    busses = pd.read_csv(path + '/comps/bus.csv', sep=';', decimal='.')
    # create busses
    nw.imp_busses = {}
    if len(busses) > 0:
        busses['instance'] = busses.apply(construct_busses, axis=1)

        # add components to busses
        comps.apply(busses_add_comps, axis=1, args=(busses, chars,))

        # add busses to network
        for b in busses['instance']:
            nw.add_busses(b)
            nw.imp_busses[b.label] = b

        print('Created busses.')

    return nw


# %% create components


def construct_comps(c, *args):
    r"""
    Creates TESPy component from class name provided in the .csv-file and specifies its parameters.

    Parameters
    ----------
    c : pandas.core.series.Series
        Component information from .csv-file.

    args[0] : pandas.core.frame.DataFrame
        DataFrame containing the x and y data of characteristic functions.

    args[1] : pandas.core.frame.DataFrame
        DataFrame containing the x, y, z1 and z2 data of characteristic maps.

    Returns
    -------
    instance : tespy.components.components.component
        TESPy component object.
    """
    if c.interface:
        instance = cmp.subsys_interface(c.label, num_inter=1)
    else:
        target_class = getattr(cmp, c.cp)
        instance = target_class(c.label)
    kwargs = {}

    # basic properties
    for key in ['mode', 'design', 'offdesign']:
        kwargs[key] = c[key]

    for key, value in instance.attr().items():
        if key in c:
            # component parameters
            if isinstance(value, hlp.dc_cp):
                dc = hlp.dc_cp(val=c[key],
                               is_set=c[key + '_set'],
                               is_var=c[key + '_var'])
                kwargs[key] = dc
            # component characteristics
            elif isinstance(value, hlp.dc_cc):
                # finding x and y values of the characteristic function
                values = args[0]['id'] == c[key]

                try:
                    x = args[0][values].x.values[0]
                    y = args[0][values].y.values[0]
                    char = cmp_char.characteristics(x=x, y=y, method=c[key + '_method'])
                except IndexError:
                    # if characteristics are missing (for compressor map atm)
                    x = cmp_char.characteristics().x
                    y = cmp_char.characteristics().y

                dc = hlp.dc_cc(is_set=c[key + '_set'],
                               method=c[key + '_method'],
                               param=c[key + '_param'],
                               func=char,
                               x=x, y=y)
                kwargs[key] = dc
            # component characteristics
            elif isinstance(value, hlp.dc_cm):
                # finding x and y values of the characteristic function
                values = args[1]['id'] == c[key]

                try:
                    x = list(args[1][values].x.values[0])
                    y = list(args[1][values].y.values[0])
                    z1 = list(args[1][values].z1.values[0])
                    z2 = list(args[1][values].z2.values[0])
                    char_map = cmp_char.char_map(x=x, y=y, z1=z1, z2=z2, method=c[key + '_method'])
                except IndexError:
                    # if characteristics are missing (for compressor map atm)
                    x = cmp_char.char_map().x
                    y = cmp_char.char_map().y
                    z1 = cmp_char.char_map().z1
                    z2 = cmp_char.char_map().z2

                dc = hlp.dc_cc(is_set=c[key + '_set'],
                               method=c[key + '_method'],
                               param=c[key + '_param'],
                               func=char_map,
                               x=x, y=y, z1=z1, z2=z2)
                kwargs[key] = dc
            # grouped component parameters
            elif isinstance(value, hlp.dc_gcp):
                dc = hlp.dc_cp(method=c[key])
                kwargs[key] = dc
            else:
                continue

    instance.set_attr(**kwargs)
    return instance


def get_interface(c):
    r"""
    Checks, if a component is marked as interface.

    Parameters
    ----------
    c : pandas.core.series.Series
        Component information from .csv-file.

    Returns
    -------
    is_interface : bool
        Returns True, if component is marked as interface.
    """
    if c.interface:
        return True
    else:
        return False

# %% create network object


def construct_network(path):
    r"""
    Creates TESPy network from the data provided in the netw.csv-file.

    Parameters
    ----------
    path : String
        Base-path to stored network data.

    Returns
    -------
    nw : tespy.networks.network
        TESPy network object.
    """
    # read network .csv-file
    netw = pd.read_csv(path + '/netw.csv', sep=';', decimal='.',
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
    nw = nwk.network(fluids=f_list, **kwargs)

    return nw

# %% create connections


def construct_conns(c, *args):
    r"""
    Creates TESPy connection from data in the .csv-file and specifies its parameters.

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
    conn = con.connection(args[0].instance[c.s], c.s_id, args[0].instance[c.t], c.t_id)

    kwargs = {}
    # read basic properties
    for key in ['design', 'offdesign']:
        if key in c:
            kwargs[key] = c[key]

    # read fluid properties
    for key in ['m', 'p', 'h', 'T', 'x', 'v']:
        if key in c:
            dc = hlp.dc_prop(val=c[key], val0=c[key + '0'], val_set=c[key + '_set'],
                             unit=c[key + '_unit'], unit_set=c[key + '_unit_set'],
                             ref=None, ref_set=c[key + '_ref_set'])
            kwargs[key] = dc

    # read fluid vector
    val = {}
    val0 = {}
    val_set = {}
    for key in args[1].fluids:
        if key in c:
            val[key] = c[key]
            val0[key] = c[key + '0']
            val_set[key] = c[key + '_set']

    kwargs['fluid'] = hlp.dc_flu(val=val, val0=val0, val_set=val_set, balance=c['balance'])

    # write properties to connection and return connection object
    conn.set_attr(**kwargs)
    return conn

# %% set references on connections


def conns_set_ref(c, *args):
    r"""
    Sets references on connections as specified in connection data.

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
            c['instance'].get_attr(col).ref = con.ref(instance,
                                                      c[col + '_ref_f'],
                                                      c[col + '_ref_d'])

# %% create busses


def construct_busses(c, *args):
    r"""
    Creates busses of the network.

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
    b = con.bus(c.label, P=c.P)
    b.P.val_set = c.P_set
    return b

# %% add components to busses


def busses_add_comps(c, *args):
    r"""
    Adds components to busses according to data from .csv file.

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
        char = cmp_char.characteristics(x=args[1][values].x.values[0], y=args[1][values].y.values[0])

        # add component with corresponding details to bus
        args[0].instance[b == args[0]['id']].values[0].add_comps(
                {'c': c.instance, 'p': p, 'P_ref': P_ref, 'char': char})
        i += 1
