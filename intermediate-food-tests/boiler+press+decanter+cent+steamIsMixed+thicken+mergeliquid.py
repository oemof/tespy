# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 12:00:18 2022

@author: mrk


This model includes 
- Boiler 
- Boiler heat loss
- Press "modelled as separator"
- Presswater Boiler 
- Decanter
- centrifuge
- steam is mixed into product in the presswater


"""
# %%

from tespy.components import Sink, Source, HeatExchangerSimple 
from tespy.connections import Connection
from tespy.networks import Network
import shutil
import numpy as np
import matplotlib.pyplot as plt

from tespy.components import Separator,Merge,CycleCloser,Valve

import logging
logging.basicConfig(level=logging.DEBUG)

# %% Separating the flows (the press machine) 


fluid_list = ['INCOMP::FoodWater','INCOMP::FoodProtein','INCOMP::FoodFat']
network = Network(fluids=fluid_list, m_unit='kg / s', p_unit='bar', T_unit='C',h_unit='kJ / kg', h_range=[-1e8,1e8], iterinfo=True)

source           = Source('source')
boiler           = HeatExchangerSimple('boiler')
boilerheatloss   = HeatExchangerSimple('boilerheatloss')
press            = Separator('press', num_out=2)
#presscake        = Sink('presscake')
presswaterboiler = HeatExchangerSimple('presswaterboiler')
decanter         = Separator('decanter', num_out=2)
#grax             = Sink('grax')
centrifuge       = Separator('centrifuge',num_out=2)

oil              = Sink('oil')
steamsource      = Source('steamsource')
steammerge       = Merge('steammerge', num_in = 2)
thickener        = Separator('thickener',num_out=2)
vapourextract    = Sink('vapourextract')

liquidmerge      = Merge('liquidmerge', num_in = 3)
thickenedliquid  = Sink('thickenedliquid')

valvemerge2 = Valve('valvemerge2')
valvemerge3 = Valve('valvemerge3')

c1 = Connection(source, 'out1', boiler, 'in1')
c2 = Connection(boiler, 'out1', boilerheatloss, 'in1')
c3 = Connection(boilerheatloss, 'out1', press, 'in1')
c4 = Connection(press, 'out1', liquidmerge, 'in1')
c5 = Connection(press, 'out2', presswaterboiler, 'in1')
c6 = Connection(presswaterboiler, 'out1', steammerge, 'in1') 
c6_1 = Connection(steamsource, 'out1', steammerge, 'in2') 
c6_2 = Connection(steammerge, 'out1', decanter, 'in1')
c7  = Connection(decanter, 'out1', valvemerge2, 'in1')
c7o = Connection(valvemerge2, 'out1', liquidmerge, 'in2')
c8 = Connection(decanter, 'out2', centrifuge, 'in1')
c10 = Connection(centrifuge, 'out2', oil, 'in1')
c9 = Connection(centrifuge, 'out1', thickener, 'in1')
c11 = Connection(thickener, 'out1', vapourextract, 'in1')
c12 = Connection(thickener, 'out2', valvemerge3, 'in1')
c12o = Connection(valvemerge3, 'out1', liquidmerge, 'in3')
c13 = Connection(liquidmerge, 'out1', thickenedliquid, 'in1')

network.add_conns(c1,c2,c3,c4,c5,c6,c6_1,c6_2,c7,c7o,c8,c9,c10,c11,c12,c12o,c13)

c1.set_attr(fluid={'FoodWater': 0.81,'FoodProtein': 0.17,'FoodFat': 0.02}, m=50000/3600, T=0, p=1)
c2.set_attr(T=95,p=1)
c3.set_attr(T=90,p=1)

# TS split (press)
total_flow_in      = 50000/3600
protein_flow_in    = (total_flow_in * 0.17)
fat_flow_in        = (total_flow_in * 0.02)
presscake_protein_flow  = protein_flow_in * 0.65
presswater_protein_flow = protein_flow_in * 0.35
presscake_flow  = presscake_protein_flow/0.435
presswater_flow = total_flow_in - presscake_flow


c4.set_attr(fluid={'FoodWater': 0.54, 'FoodFat': 0.025, 'FoodProtein': 0.435},m=presscake_flow)
logging.basicConfig(level=logging.DEBUG)
#c5.set_attr(fluid0={'FoodWater': 0.90, 'FoodFat': 0.018, 'FoodProtein': 0.08})

# initial conditions
# concerning initial conditions, fluid composition must always be specified
# hmmm. (why can T0 not be used to set h0 here in the separator outlets)
#c4.set_attr(h0=1e2,T0=90,p0=1,m0=presscake_flow,fluid0={'FoodWater': 0.33, 'FoodFat': 0.33, 'FoodProtein': 0.33})
#c5.set_attr(h0=1e2,T0=90,p0=1,m0=presswater_flow,fluid0={'FoodWater': 0.33, 'FoodFat': 0.33, 'FoodProtein': 0.33})

c4.set_attr(h0=1e2,p0=1,fluid0={'FoodWater': 0.33, 'FoodFat': 0.33, 'FoodProtein': 0.33})
c5.set_attr(h0=1e2,p0=1,fluid0={'FoodWater': 0.33, 'FoodFat': 0.33, 'FoodProtein': 0.33})
c6.set_attr(T=100,p=1,fluid0={'FoodWater': 0.33, 'FoodFat': 0.33, 'FoodProtein': 0.33})

# TS split (decanter)
grax_protein_flow     = presswater_protein_flow * 0.35
decanter_protein_flow = presswater_protein_flow * 0.65
grax_flow             = grax_protein_flow/0.33

c7.set_attr(fluid={'FoodWater': 0.648, 'FoodFat': 0.022, 'FoodProtein': 0.330},m=grax_flow)
c7.set_attr(h0=1e2,p0=1,fluid0={'FoodWater': 0.33, 'FoodFat': 0.33, 'FoodProtein': 0.33})
c8.set_attr(h0=1e2,p0=1,fluid0={'FoodWater': 0.33, 'FoodFat': 0.33, 'FoodProtein': 0.33})

# TF split (decanter + centrifuge)
presscake_fat_flow = presscake_flow * 0.025
fat_decanter_centrifuge_flow = fat_flow_in-presscake_fat_flow
fat_oil_flow = fat_decanter_centrifuge_flow * 0.85
oil_flow = fat_oil_flow/0.999

c10.set_attr(fluid={'FoodWater': 0, 'FoodFat': 0.999, 'FoodProtein': 0.001},m=oil_flow)
c9.set_attr(h0=1e2,p0=1,fluid0={'FoodWater': 0.33, 'FoodFat': 0.33, 'FoodProtein': 0.33})
c10.set_attr(h0=1e2,p0=1,fluid0={'FoodWater': 0.33, 'FoodFat': 0.33, 'FoodProtein': 0.33})

# Adding the stream mass manually (should be calculated inside the presswaterboiler)
# pressure is propagated through separator..
# c6_1.set_attr(T=100,p=1,m=703.0669512402173,fluid={'FoodWater': 1, 'FoodFat': 0, 'FoodProtein': 0})
# c6_1.set_attr(T=100,m=703.0669512402173/3600,fluid={'FoodWater': 1, 'FoodFat': 0, 'FoodProtein': 0})
c6_1.set_attr(m=703.0669512402173/3600,fluid={'FoodWater': 1, 'FoodFat': 0, 'FoodProtein': 0})
c6_2.set_attr(T=100)

c11.set_attr(fluid={'FoodWater': 1, 'FoodFat': 0, 'FoodProtein': 0})
c12.set_attr(fluid={'FoodProtein': 0.3})
c12.set_attr(h0=1e2,p0=1,fluid0={'FoodWater': 0.33, 'FoodFat': 0.33, 'FoodProtein': 0.33})
c13.set_attr(h0=1e2,p0=1,fluid0={'FoodWater': 0.33, 'FoodFat': 0.33, 'FoodProtein': 0.33})

for c in network.conns['object']:
    c.set_attr(h0=1e2,p0=1,fluid0={'FoodWater': 0.33, 'FoodFat': 0.33, 'FoodProtein': 0.33})


network.solve('design')

print("")
print(network.results['Connection'])

print("")
print(f"boiler heat should be close to = {0.900*95*50000} kcal/h")
print(f"boiler heat should be close to = {900*4.184*95*50000/3600} J/s")
print(f"boiler.Q.val = {boiler.Q.val}")
print(f"necessary steam flow = {boiler.Q.val/(503*4184)*3600} kg/h")
print(f"boilerheatloss.Q.val = {boilerheatloss.Q.val}")
print("")

#print(f"PressCake composition is = {c4.fluid.val}")
#print(f"PressWater composition is  = {c5.fluid.val}")
#print(f"PressWater composition is  = {c6.fluid.val}")
print(f"PressWater heat is   = {presswaterboiler.Q.val}")
print(f"necessary steam flow = {presswaterboiler.Q.val/(503*4184)*3600} kg/h")

print(f"necessary water flow = {(c11.m.val*3600/3*540)/10/1000*1.1}")

print("")



print(" close ")



"""

RM = [[]]*10




RM[0] = {'m' : 50000,
         'Wpct' : 81,
         'Fpct' : 2,
         'T' : 0}

RM[0]['W'] = RM[0]['m']*RM[0]['Wpct']/100
RM[0]['F'] = RM[0]['m']*RM[0]['Fpct']/100
RM[0]['Spct'] = 100 - RM[0]['Wpct'] - RM[0]['Fpct']
RM[0]['S'] = RM[0]['m']*RM[0]['Spct']/100
# specific heat in kcal/kgC and heat of evaporation in kcal/kg
RM[0]['cp'] = (RM[0]['Wpct']/100+RM[0]['Fpct']/100*0.58+RM[0]['Spct']/100*0.44) 
RM[0]['h_fg'] =	503	

RM[1] = RM[0].copy()
RM[1]['T'] = 95
QdotBoil = RM[0]['m'] * RM[0]['cp'] * (RM[1]['T']-RM[0]['T'])
mdotBoil= QdotBoil/RM[0]['h_fg']
print(mdotBoil)

RM[2] = RM[1].copy()
RM[2]['Fpct'] = 2.5
RM[2]['Spct'] = 43.5
RM[2]['Wpct'] = 100 - RM[2]['Fpct'] - RM[2]['Spct']
TS_split = 65
RM[2]['S'] = RM[1]['S']*TS_split/100
RM[2]['m'] = RM[2]['S']/(RM[2]['Spct']/100)
RM[2]['W'] = RM[2]['m']*RM[2]['Wpct']/100
RM[2]['F'] = RM[2]['m']*RM[2]['Fpct']/100

RM[3] = RM[1].copy()
RM[3]['m'] = RM[1]['m'] - RM[2]['m']
RM[3]['W'] = RM[1]['W'] - RM[2]['W']
RM[3]['F'] = RM[1]['F'] - RM[2]['F']
RM[3]['S'] = RM[1]['S'] - RM[2]['S']
RM[3]['Wpct'] = RM[1]['W']/RM[3]['m']
RM[3]['Fpct'] = RM[1]['F']/RM[3]['m']
RM[3]['Spct'] = RM[1]['S']/RM[3]['m']
RM[3]['T'] = RM[1]['T']-5

RM[4] = RM[3].copy()
RM[4]['T'] = 100
QdotPress = RM[3]['m'] * RM[3]['cp'] * (RM[4]['T']-RM[3]['T'])
mdotPress = QdotPress/RM[4]['h_fg']
print(mdotPress)

"""


# %%
