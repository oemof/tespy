from tespy.components import Sink, Source, HeatExchangerSimple, Splitter
from tespy.connections import Connection, Ref, Bus
from tespy.networks import Network
import shutil
import numpy as np
import matplotlib.pyplot as plt

from tespy.components import Separator,Merge,CycleCloser,Valve,Splitter
from tespy.components.energySupplyComponents import *

import logging
logging.basicConfig(level=logging.DEBUG)

from CoolProp.CoolProp import PropsSI




fluid_list = ['INCOMP::Water']
network = Network(fluids=fluid_list, m_unit='kg / s', p_unit='bar', T_unit='C',h_unit='kJ / kg', h_range=[-1e2,4e3], iterinfo=True)

# Objects
PurchasedElectricity  = Source('PurchasedElectricity')
HeatPump              = MassFactorVCCEnergySupply('HeatPump')
Heating             = Sink('Consumer1')
Cooling             = Sink('Consumer2')

# Connections
c1 = Connection(PurchasedElectricity, 'out1', HeatPump, 'in1')
c2 = Connection(HeatPump, 'out1', Heating, 'in1')
c3 = Connection(HeatPump, 'out2', Cooling, 'in1')

network.add_conns(c1,c2,c3)

m0 = 4.428
c1.set_attr(m=m0)
#c2.set_attr(m=m0*3.5)

HeatPump.set_attr(COP=3)


# guess
for c in network.conns['object']:
    c.set_attr(m0=m0,h0=0,p0=1,fluid0={'INCOMP::Water': 1})

# arbitray values
c1.set_attr(h=0,p=1,fluid={'INCOMP::Water': 1}, mixing_rule="incompressible")
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
