"""
.. module:: csv_reader
    :platforms: all
    :synopsis:

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>
"""

import pandas as pd
from tespy import cmp, con, nwk
import csv


class process_director:
    """
    documentation
    """
    def __init__(self):
        self.imported_comps = pd.DataFrame(columns=['comp'])

    def construct_component(self, component, label):
        target_class = getattr(cmp, component)
        instance = target_class(label)
        self.imported_comps.loc[label] = [instance]


def load_nw(filename):
    """
    generate network class providing result files
    Required files:

    - ..._comps.csv
    - ..._results.csv

    TODO: REDO!
    """

# open csv component file and generate component objects
    director = process_director()
    with open(filename + '_comps.csv', 'r') as csvfile:
        f = csv.reader(csvfile, delimiter=';')
        next(f, None)
        i = 0
        for row in f:
            numel = (len(row) - 7) / 2
            director.construct_component(row[1], row[0])
# set attributes (P,Q,eta, ...) of the component
            if numel != 0:
                for k in range(int(numel)):
                    director.imported_comps['comp'].iloc[i].set_attr(
                         **{row[7 + 2 * k]: float(row[7 + 2 * k + 1])})
            i += 1

# create network object, open csv file with results
    conns = pd.read_csv(filename + '_results.csv', sep = ';', decimal = '.')

    ix = 6

    fluids = list(conns.columns[ix:conns.columns.get_loc('m_set')])
    print(fluids)

    imported_conns = []
    for i in range(len(conns)):
# create connection
        imported_conns += [con.connection(
                        director.imported_comps.loc[conns['s'].iloc[i]].comp,
                        conns['s_id'].iloc[i],
                        director.imported_comps.loc[conns['t'].iloc[i]].comp,
                        conns['t_id'].iloc[i])]
# set properties
        if conns['m_set'].iloc[i]:
            imported_conns[i].set_attr(m = conns['m'].iloc[i])
        if conns['p_set'].iloc[i]:
            imported_conns[i].set_attr(p = conns['p'].iloc[i])
        if conns['h_set'].iloc[i]:
            imported_conns[i].set_attr(h = conns['h'].iloc[i])
        if conns['T_set'].iloc[i]:
            imported_conns[i].set_attr(T = conns['T'].iloc[i])
        if conns['x_set'].iloc[i]:
            imported_conns[i].set_attr(x = conns['x'].iloc[i])

    nw = nwk.network(fluids=fluids)
# add connections to network
    for c in imported_conns:
        nw.add_conns(c)

    return nw
