# -*- coding: utf-8 -*-
"""
@author: witte
"""
from tespy import con, nwk, cmp
from CoolProp.CoolProp import PropsSI
import math
import matplotlib.pyplot as plt

# %% components

# main components
turbine=cmp.turbine('turbine',eta_s=0.9)
condenser=cmp.condenser('condenser',dp1=1,dp2=1)
pump=cmp.pump('pump',eta_s=0.8)
steam_generator=cmp.heat_exchanger_simple('steam generator',dp=0.95,mode='man')
source=cmp.source('source1')
sink=cmp.sink('sink1')

# cooling water
source_cw=cmp.source('source2')
sink_cw=cmp.sink('sink2')

# %% connections

# main cycle
fs=con.connection(source,'out1',turbine,'in1',p=100,T=500,m=100,fluid={'water':1})
ws=con.connection(turbine,'out1',condenser,'in1',p=1)
cond=con.connection(condenser,'out1',pump,'in1')
fw=con.connection(pump,'out1',steam_generator,'in1')
# closing the cylce:
# fluid properties are the same as at the turbine inlet
out=con.connection(steam_generator,'out1',sink,'in1',p=con.ref(fs,1,0),h=con.ref(fs,1,0))

# cooling water
cw_in=con.connection(source_cw,'out1',condenser,'in2',T=20,p=5,fluid={'water':1})
cw_out=con.connection(condenser,'out2',sink_cw,'in1',T=30)

# %% network

fluids=['water']

nw=nwk.network(fluids=fluids,p='bar',T='C',
               p_range=[0.01, 150], T_range=[5, 800])
nw.add_conns(fs,ws,cond,fw,cw_in,cw_out,out)
# %% solving

mode='design'

file='simple_cr_'+mode+'_conn.csv'

# solve the network and process the results
nw.solve(mode=mode)
nw.process_components(mode='post')

# save the results
nw.save('simple_cr_'+mode)

# change to offdesign mode
mode='offdesign'

# unset pressure at turbine inlet and turbine outlet
# reset mass fresh steam mass flow
# turbine inlet pressure is deriven by stodolas law, outlet pressure by
# characteristic of condenser
fs.set_attr(p=math.nan,m=90)
ws.set_attr(p=math.nan)

nw.solve(mode=mode,design_file=file,init_file=file)
nw.process_components(mode='post')

nw.save('simple_cr_'+mode)
