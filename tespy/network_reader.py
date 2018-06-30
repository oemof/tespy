"""
.. module:: network_reader
    :synopsis:

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>

TODO: add documentation
"""

import pandas as pd
from tespy import cmp, con, nwk, hlp
import os
import ast

# %% network loading


def load_nw(path):
    """
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

    # load components
    comps = pd.DataFrame()

    files = os.listdir(path + '/comps/')
    for f in files:
        if f != 'bus.csv' and f != 'char.csv':
            df = pd.read_csv(path + '/comps/' + f, sep=';', decimal='.',
                             converters={'design': ast.literal_eval,
                                         'offdesign': ast.literal_eval,
                                         'busses': ast.literal_eval,
                                         'bus_factors': ast.literal_eval})

            # create components
            df['instance'] = df.apply(construct_comps, axis=1, args=(chars,))
            comps = pd.concat((comps, df[['instance', 'label', 'busses',
                                          'bus_factors']]),
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
        comps.apply(busses_add_comps, axis=1, args=(busses,))

        # add busses to network
        for b in busses['instance']:
            nw.add_busses(b)

        print('Created busses')

    return nw


# %% create components


def construct_comps(c, *args):

    target_class = getattr(cmp, c.cp)
    instance = target_class(c.label)
    kwargs = {}
    for key in ['mode', 'design', 'offdesign']:
        kwargs[key] = c[key]

    for key, value in instance.attr_prop().items():
        if isinstance(value, hlp.dc_cp):
            dc = hlp.dc_cp(val=c[key],
                           is_set=c[key + '_set'],
                           is_var=c[key + '_var'])
            kwargs[key] = dc
        elif isinstance(value, hlp.dc_cc):
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

    nw = nwk.network(fluids=f_list, **kwargs)

    return nw

# %% create connections


def construct_conns(c, *args):

    conn = con.connection(args[0].instance[c.s], c.s_id,
                          args[0].instance[c.t], c.t_id)

    kwargs = {}
    for key in ['design', 'offdesign']:
        kwargs[key] = c[key]

    for key in ['m', 'p', 'h', 'T', 'x']:
        dc = hlp.dc_prop(val=c[key],
                         val0=c[key + '0'],
                         val_set=c[key + '_set'],
                         unit=c[key + '_unit'],
                         unit_set=c[key + '_unit_set'],
                         ref=None,
                         ref_set=c[key + '_ref_set'])

        kwargs[key] = dc

    val = {}
    val0 = {}
    val_set = {}
    for key in args[1].fluids:
        val[key] = c[key]
        val0[key] = c[key + '0']
        val_set[key] = c[key + '_set']

    kwargs['fluid'] = hlp.dc_flu(val=val, val0=val0, val_set=val_set,
                                 balance=c['balance'])

    conn.set_attr(**kwargs)

    return conn

# %% set references on connections


def conns_set_ref(c, *args):

    for col in ['m', 'p', 'h', 'T']:
        if isinstance(c[col + '_ref'], str):
            instance = args[0].instance[c[col + '_ref'] ==
                                        args[0]['id']].values[0]
            c['instance'].get_attr(col).ref = con.ref(instance,
                                                      c[col + '_ref_f'],
                                                      c[col + '_ref_d'])

# %% create busses


def construct_busses(c, *args):

    b = con.bus(c.label, P=c.P)
    b.P_set = c.P_set
    return b

# %% add components to busses


def busses_add_comps(c, *args):

    i = 0
    for b in c.busses:
        f = c.bus_factors[i]
        args[0].instance[b == args[0]['id']
                         ].values[0].add_comps([c.instance, f])
        i += 1
