# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 10:37:36 2017

@author: witte
"""

from tespy import con, cmp, nwk

# %% network

# define full fluid list for the network's variable space
fluid_list = ['Ar', 'N2', 'O2', 'CO2', 'CH4', 'H2O']
# define unit systems and fluid property ranges
nw = nwk.network(fluids=fluid_list, p_unit='bar', T_unit='C',
                 p_range=[0.5, 10], T_range=[10, 1200])

# %% components

# sinks & sources
amb = cmp.source('ambient')
sf = cmp.source('fuel')
fg = cmp.sink('flue gas outlet')

# combustion chamber
comb = cmp.combustion_chamber(label='combustion chamber')

# %% connections

amb_comb = con.connection(amb, 'out1', comb, 'in1')
sf_comb = con.connection(sf, 'out1', comb, 'in2')
comb_fg = con.connection(comb, 'out1', fg, 'in1')

nw.add_conns(sf_comb, amb_comb, comb_fg)

# %% component parameters

# set combustion chamber fuel, air to stoichometric air ratio and thermal input
comb.set_attr(fuel='CH4', lamb=3, ti=20000)

# %% connection parameters

# air from abient (ambient pressure and temperature), air composition must be
# stated component wise.
amb_comb.set_attr(p=1, T=20,
                  fluid={'Ar': 0.0129, 'N2': 0.7553, 'H2O': 0,
                         'CH4': 0, 'CO2': 0.0004, 'O2': 0.2314})

# fuel, pressure must not be stated, as pressure is the same at all inlets and
# outlets of the combustion chamber
sf_comb.set_attr(T=25,
                 fluid={'CO2': 0.04, 'Ar': 0, 'N2': 0,
                        'O2': 0, 'H2O': 0, 'CH4': 0.96})

# state desired flue gas temperature
comb_fg.set_attr()

# %% solving

mode = 'design'
nw.solve(mode=mode)
nw.print_results()
nw.save('combustion')
