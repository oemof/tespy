# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 12:00:18 2022

@author: mrk

This model includes
- Boiler


"""

from tespy.components import Sink, Source, HeatExchangerSimple
from tespy.connections import Connection
from tespy.networks import Network
import shutil
import numpy as np
import matplotlib.pyplot as plt

from tespy.components import Separator,Merge,CycleCloser,Valve
from tespy.components.newcomponents import DiabaticSimpleHeatExchanger,MergeWithPressureLoss,SeparatorWithSpeciesSplits

import logging
#logging.basicConfig(level=logging.DEBUG)


network = Network(m_unit='kg / s', p_unit='bar', T_unit='C',h_unit='kJ / kg', h_range=[-1e2,4e2], iterinfo=True)

# Objects
source           = Source('source')
boiler           = HeatExchangerSimple('boiler')
sink             = Sink('sink')

# Connections
c1 = Connection(source, 'out1', boiler, 'in1')
c2 = Connection(boiler, 'out1', sink, 'in1')

network.add_conns(c1,c2)

# set global guess values
m0 = 100      # transform unit at some point
h0 = 1e2      # global guess value in kJ/kg
p0 = 2        # global guess value in bar

# set conditions around boiler
# fluid_back_ends={'Water': "INCOMP", "PHE": "INCOMP", "S800": "INCOMP"}
c1.set_attr(fluid={'INCOMP::Water': 0.80,'INCOMP::PHE': 0.15,'INCOMP::S800': 0.05}, m=100, h=h0, p=p0, mixing_rule="incompressible")
c2.set_attr(h=h0,p=p0)

# solve and print results
network.solve('design')
# c1.set_attr(T=40, h=None)
# c2.set_attr(T=60, h=None)
# network.solve('design')

network.print_results()


print(network.results['Connection'])

network.results["Connection"].to_csv(f"{__file__.replace('.py', '')}tespy070.csv")