# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 12:00:18 2022

@author: mrk

This model includes
- Boiler
- Press
- Decanter (no steam is mixed into product in the pressWater)
- Centrifuge

"""
from tespy.components import Sink, Source, SimpleHeatExchanger
from tespy.connections import Connection
from tespy.networks import Network


from tespy.components import Separator,Merge,CycleCloser,Valve
from tespy.components.newcomponents import DiabaticSimpleHeatExchanger,MergeWithPressureLoss,SeparatorWithSpeciesSplits

import logging
#logging.basicConfig(level=logging.DEBUG)

network = Network(m_unit='kg / s', p_unit='bar', T_unit='C',h_unit='kJ / kg', h_range=[-1e2,4e2], iterinfo=False)

# Objects
source              = Source('source')
boiler              = SimpleHeatExchanger('boiler')
press               = SeparatorWithSpeciesSplits('press', num_out=2)
#pressWater          = Sink('pressWater')
presscake           = Sink('presscake')
decanter            = SeparatorWithSpeciesSplits('decanter', num_out=2)
grax                = Sink('grax')
oil                = Sink('oil')
centrifuge          = SeparatorWithSpeciesSplits('centrifuge',num_out=2)
stickWater          = Sink('stickWater')

# Connections
c1 = Connection(source, 'out1', boiler, 'in1')
c2 = Connection(boiler, 'out1', press, 'in1')
c3 = Connection(press, 'out1', presscake, 'in1')
c4 = Connection(press, 'out2', decanter, 'in1')
c5 = Connection(decanter, 'out1', grax, 'in1')
c6 = Connection(decanter, 'out2', centrifuge, 'in1')
c7 = Connection(centrifuge, 'out1', stickWater, 'in1')
c8 = Connection(centrifuge, 'out2', oil, 'in1')

network.add_conns(c1,c2,c3,c4,c5,c6,c7,c8)

# set global guess values
m0 = 100       # transform unit at some point
h0 = 1e2       # global guess value in kJ/kg
p0 = 2         # global guess value in bar

for c in network.conns['object']:
#     n_fl = len(network.fluids)
    c.set_attr(m0=m0,h0=h0,p0=p0)#,fluid0={'INCOMP::Water': 1/3, 'INCOMP::PHE': 1/3, 'INCOMP::S800': 1/3})

# set conditions around boiler
c1.set_attr(fluid={'INCOMP::Water': 0.8,'INCOMP::PHE': 0.15,'INCOMP::S800': 0.05}, m=m0, h=h0, p=p0, mixing_rule="incompressible")
c2.set_attr(h=h0,p=p0)

# set conditions around press
press.set_attr(SFS={
    'val': 0.7, 'is_set': True,
    'split_fluid' : 'S800', 'split_outlet' : "out1"})
c3.set_attr(fluid={'INCOMP::Water': 0.50, 'INCOMP::PHE': 0.05, 'INCOMP::S800': 0.45})

# set conditions around decanter
decanter.set_attr(SFS={
    'val': 0.3, 'is_set': True,
    'split_fluid' : 'S800', 'split_outlet' : "out1"})
c5.set_attr(fluid={'INCOMP::Water': 0.6, 'INCOMP::PHE': 0.05, 'INCOMP::S800': 0.35})

# set conditions around centrifuge
centrifuge.set_attr(SFS={
    'val': 0.8, 'is_set': True,
    'split_fluid' : 'PHE', 'split_outlet' : "out2"})
c8.set_attr(fluid={'INCOMP::PHE': 0.99, 'INCOMP::S800': 0.01})

# solve and print results
network.solve('design')

# network.print_results()
print(network.results['Connection'].loc[:, [c for c in network.results["Connection"] if "unit" not in c]])

print(0.7 * c1.m.val_SI * c1.fluid.val["S800"] / c3.fluid.val["S800"] == network.results['Connection'].loc["press:out1_presscake:in1", "m"])

network.results["Connection"].to_csv(f"{__file__.replace('.py', '')}tespy070.csv")