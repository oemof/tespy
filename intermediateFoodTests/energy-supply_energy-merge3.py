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
PurchasedNaturalGas1   = Source('PurchasedNaturalGas1')
PurchasedNaturalGas2   = Source('PurchasedNaturalGas2')

Merge                  = MergeEnergySupply('Merge',num_in = 2)

Consumer               = Sink('Consumer')


# Connections
c1 = Connection(PurchasedNaturalGas1, 'out1', Merge, 'in1')
c2 = Connection(PurchasedNaturalGas2, 'out1', Merge, 'in2')
c3 = Connection(Merge, 'out1', Consumer, 'in1')

network.add_conns(c1,c2,c3)

m0 = 4.428

c1.set_attr(m=m0)
c2.set_attr(m=m0/4)

# guess
for c in network.conns['object']:
    c.set_attr(m0=m0,h0=0,p0=1,fluid0={'FoodWater': 1})

# arbitray values
c1.set_attr(h=0,p=1,fluid={'FoodWater': 1})
c2.set_attr(h=0,p=1,fluid={'FoodWater': 1})
c3.set_attr(h=0,p=1)
# merge already propergate h, p and fluid

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

