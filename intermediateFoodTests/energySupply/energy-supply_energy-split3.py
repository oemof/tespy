from tespy.components import Sink, Source, HeatExchangerSimple, Splitter
from tespy.connections import Connection, Ref, Bus
from tespy.networks import Network
import shutil
import numpy as np
import matplotlib.pyplot as plt

from tespy.components import Separator,Merge,CycleCloser,Valve,Splitter
from tespy.components.newcomponents import *

import logging
logging.basicConfig(level=logging.DEBUG)

from CoolProp.CoolProp import PropsSI




fluid_list = ['INCOMP::FoodWater']
network = Network(fluids=fluid_list, m_unit='kg / s', p_unit='bar', T_unit='C',h_unit='kJ / kg', h_range=[-1e2,4e3], iterinfo=True)

# Objects
PurchasedNaturalGas   = Source('PurchasedNaturalGas')

Splitter1              = SplitterEnergySupply('Splitter1',num_out = 2)

Consumer1             = Sink('Consumer1')
Consumer2             = Sink('Consumer2')



# Connections
c1 = Connection(PurchasedNaturalGas, 'out1', Splitter1, 'in1')
c2 = Connection(Splitter1, 'out1', Consumer1, 'in1')
c3 = Connection(Splitter1, 'out2', Consumer2, 'in1')

network.add_conns(c1,c2,c3)

m0 = 4.428

c1.set_attr(m=m0)
c2.set_attr(m=m0/4)

# guess
for c in network.conns['object']:
    c.set_attr(m0=m0,h0=0,p0=1,fluid0={'FoodWater': 1})

# arbitray values
c1.set_attr(h=0,p=1,fluid={'FoodWater': 1})
c2.set_attr(h=0,p=1)
c3.set_attr(h=0,p=1)

network.solve('design',init_only=True)

for c in network.conns['object']:
    print(c.p.val_SI)
for c in network.conns['object']:
    print(c.h.val_SI)
for c in network.conns['object']:
    print(c.T.val_SI)

# solve and print results
network.solve('design')

network.print_results()
print(network.results['Connection'])

