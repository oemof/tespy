# -*- coding: utf-8 -*-
"""
@author: witte
"""
from tespy import con, nwk, cmp, hlp
import math

# %% network

fluids = ['water']

nw = nwk.network(fluids=fluids)
nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg',
            p_range=[0.01, 150], T_range=[5, 800], h_range=[10, 5000])

# %% components

# main components
turbine = cmp.turbine('turbine')
condenser = cmp.condenser('condenser')
pump = cmp.pump('pump')
steam_generator = cmp.heat_exchanger_simple('steam generator')
source = cmp.source('cycle opener')
sink = cmp.sink('cycle closer')

# cooling water
source_cw = cmp.source('cooling water inlet')
sink_cw = cmp.sink('cooling water outlet')

# %% connections

# main cycle
fs_in = con.connection(source, 'out1', turbine, 'in1')
ws = con.connection(turbine, 'out1', condenser, 'in1')
cond = con.connection(condenser, 'out1', pump, 'in1')
fw = con.connection(pump, 'out1', steam_generator, 'in1')
fs_out = con.connection(steam_generator, 'out1', sink, 'in1')
nw.add_conns(fs_in, ws, cond, fw, fs_out)

# cooling water
cw_in = con.connection(source_cw, 'out1', condenser, 'in2')
cw_out = con.connection(condenser, 'out2', sink_cw, 'in1')
nw.add_conns(cw_in, cw_out)

# %% busses

power = con.bus('total output power')
power.add_comps([turbine, 1], [pump, 1])
nw.add_busses(power)

# %% parametrization of components

turbine.set_attr(eta_s=0.9)
condenser.set_attr(pr1=0.98, pr2=0.98, ttd_u=5)
pump.set_attr(eta_s=0.8)
steam_generator.set_attr(pr=0.95, mode='man')

# %% parametrization of connections

# offdesign calculation: use parameter design for auto deactivation
# turbine inlet pressure is deriven by stodolas law, outlet pressure by
# characteristic of condenser
fs_in.set_attr(p=100, T=500, fluid={'water': 1}, design=['p'])

# closing the cycle: fluid properties at sink must be identical to fluid
# properties at the source
fs_out.set_attr(p=con.ref(fs_in, 1, 0), h=con.ref(fs_in, 1, 0))

cw_in.set_attr(T=20, p=5, fluid={'water': 1})
cw_out.set_attr(T=30)

# total output power as input parameter
power.set_attr(P=-10e6)

# %% solving

mode = 'design'

file = 'cr'

# solve the network, print the results to prompt and save
nw.solve(mode=mode)
nw.print_results()
nw.save(file)

# change to offdesign mode
mode = 'offdesign'

# reset power input
power.set_attr(P=-9e6)

# the design file holds the information on the design case
# initialisation from previously design process
nw.solve(mode=mode, design_file=file + '/results.csv',
         init_file=file + '/results.csv')
nw.print_results()
nw.save('cr_OD')
