"""
.. module:: network_reader
    :synopsis:

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>
"""

import pandas as pd
from tespy import cmp, con, nwk, hlp, cmp_char
import os
import ast

# %% network loading


def load_nwk(path):
    r"""
    loads a network from a path given the following structure:

    - path (e. g. 'mynetwork')
        - comps
            - bus.csv
            - char.csv
            - heat_exchanger.csv
            - ...
        - conns.csv
        - netw.csv

    .. note::
        If you save the network structure of an existing TESPy network, it will
        be stored this way for you. The returned network object is ready for
        calculation (given good parametrisation).

    :param path: path to stored network data
    :type path: str
    :returns: nw (*tespy.networks.network*) - TESPy network object
    """
    print('Reading network data...')

    # load characteristics
    chars = pd.read_csv(path + '/comps/char.csv', sep=';', decimal='.',
                        converters={'x': ast.literal_eval,
                                    'y': ast.literal_eval})
    chars['char'] = chars.apply(construct_chars, axis=1)

    # load components
    comps = pd.DataFrame()

    files = os.listdir(path + '/comps/')
    for f in files:
        if f != 'bus.csv' and f != 'char.csv':
            df = pd.read_csv(path + '/comps/' + f, sep=';', decimal='.',
                             converters={'design': ast.literal_eval,
                                         'offdesign': ast.literal_eval,
                                         'busses': ast.literal_eval,
                                         'bus_param': ast.literal_eval,
                                         'bus_P_ref': ast.literal_eval,
                                         'bus_char': ast.literal_eval})

            # create components
            df['instance'] = df.apply(construct_comps, axis=1, args=(chars,))
            comps = pd.concat((comps, df[['instance', 'label', 'busses',
                                          'bus_param', 'bus_P_ref',
                                          'bus_char']]),
                              axis=0)

    comps = comps.set_index('label')
    print('Created components')

    # create network
    nw = construct_network(path)

    # load connections
    conns = pd.read_csv(path + '/conn.csv', sep=';', decimal='.',
                        converters={'design': ast.literal_eval,
                                    'offdesign': ast.literal_eval})

    # create connections
    conns['instance'] = conns.apply(construct_conns, axis=1, args=(comps, nw,))
    conns.apply(conns_set_ref, axis=1, args=(conns,))
    conns = conns.set_index('id')

    # add connections to network
    for c in conns['instance']:
        nw.add_conns(c)

    print('Created connections')
    nw.check_network()

    # load busses
    busses = pd.DataFrame()
    busses = pd.read_csv(path + '/comps/bus.csv', sep=';', decimal='.')
    # create busses
    if len(busses) > 0:
        busses['instance'] = busses.apply(construct_busses, axis=1)

        # add components to busses
        comps.apply(busses_add_comps, axis=1, args=(busses, chars,))

        # add busses to network
        for b in busses['instance']:
            nw.add_busses(b)

        print('Created busses')

    return nw


# %% create components


def construct_comps(c, *args):
    r"""
    creates TESPy component from class name provided in the .csv-file and
    specifies its parameter

    :param c: component information
    :type c: pandas.core.series.Series
    :returns: instance (*tespy.components.component*) - TESPy component object

    **additional arguments in args**

    - args[0]: char (*pandas.core.frame.DataFrame*) - DataFrame containing the
      x and y data of characteristic functions
    """

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
                x = args[0][values]['x'].values[0]
                y = args[0][values]['y'].values[0]
                dc = hlp.dc_cc(is_set=c[key + '_set'],
                               method=c[key + '_method'],
                               param=c[key + '_param'],
                               x=x, y=y)
                kwargs[key] = dc
            else:
                continue

    instance.set_attr(**kwargs)
    return instance

# %% create network object


def construct_network(path):
    r"""
    creates TESPy network from the data provided in the .csv-file

    :param path: path to stored network data
    :type path: str
    :returns: nw (*tespy.networks.network*) - TESPy network object
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

# %% create network object


def construct_chars(c):
    r"""
    creates TESPy characteristic functions

    :param c: connection information
    :type c: pandas.core.series.Series
    :returns: instance (*tespy.components.characteristics*) - TESPy
              characteristics object
    """

    char = cmp_char.characteristics(x=c.x, y=c.y)
    return char

# %% create connections


def construct_conns(c, *args):
    r"""
    creates TESPy component from class name provided in the .csv-file and
    specifies its parameters

    :param c: connection information
    :type c: pandas.core.series.Series
    :returns: instance (*tespy.components.component*) - TESPy component object

    **additional arguments in args**

    - args[0]: comps (*pandas.core.frame.DataFrame*) - DataFrame containing all
      created components
    """

    # create connection
    conn = con.connection(args[0].instance[c.s], c.s_id,
                          args[0].instance[c.t], c.t_id)

    kwargs = {}
    # read basic properties
    for key in ['design', 'offdesign']:
        if key in c:
            kwargs[key] = c[key]

    # read fluid properties
    for key in ['m', 'p', 'h', 'T', 'x']:
        if key in c:
            dc = hlp.dc_prop(val=c[key],
                             val0=c[key + '0'],
                             val_set=c[key + '_set'],
                             unit=c[key + '_unit'],
                             unit_set=c[key + '_unit_set'],
                             ref=None,
                             ref_set=c[key + '_ref_set'])

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

    kwargs['fluid'] = hlp.dc_flu(val=val, val0=val0, val_set=val_set,
                                 balance=c['balance'])

    # write properties to connection and return connection object
    conn.set_attr(**kwargs)
    return conn

# %% set references on connections


def conns_set_ref(c, *args):
    r"""
    sets references on the created connections

    :param c: connection information
    :type c: pandas.core.series.Series
    :returns: instance (*tespy.components.component*) - TESPy component object

    **additional arguments in args**

    - args[0]: conns (*pandas.core.frame.DataFrame*) - DataFrame containing all
      created connections
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
    creates busses

    :param c: bus information
    :type c: pandas.core.series.Series
    :returns: b (*tespy.connection.bus*) - TESPy bus object
    """

    # set up bus with label and specify value for power
    b = con.bus(c.label, P=c.P)
    b.P.val_set = c.P_set
    return b

# %% add components to busses


def busses_add_comps(c, *args):
    r"""
    adds components to the busses

    :param c: component information
    :type c: pandas.core.series.Series
    :returns: instance (*tespy.components.component*) - TESPy component object

    **additional arguments in args**

    - args[0]: busses (*pandas.core.frame.DataFrame*) - DataFrame containing
      all created busses
    """

    i = 0
    for b in c.busses:
        p, P_ref, char = c.bus_param[i], c.bus_P_ref[i], c.bus_char[i]

        values = char == args[1]['id']
        char = args[1]['char'][values[values == True].index[0]]

        # add component with corresponding details to bus
        args[0].instance[b == args[0]['id']
                         ].values[0].add_comps({'c': c.instance,
                                                'p': p, 'P_ref': P_ref,
                                                'char': char})
        i += 1
