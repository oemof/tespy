# %%[sec_1]
from tespy.networks import Network

# create a network object with R134a as fluid
my_plant = Network()
# %%[sec_2]
# set the unitsystem for temperatures to °C and for pressure to bar
my_plant.units.set_defaults(
    temperature="degC", pressure="bar", enthalpy="kJ/kg"
)
# %%[sec_3]
from tespy.components import (
    CycleCloser, Compressor, Valve, SimpleHeatExchanger
)

cc = CycleCloser('cycle closer')

# heat sink
co = SimpleHeatExchanger('condenser')
# heat source
ev = SimpleHeatExchanger('evaporator')

va = Valve('expansion valve')
cp = Compressor('compressor')
# %%[sec_4]
from tespy.connections import Connection

# connections of heat pump
c1 = Connection(cc, 'out1', ev, 'in1', label='1')
c2 = Connection(ev, 'out1', cp, 'in1', label='2')
c3 = Connection(cp, 'out1', co, 'in1', label='3')
c4 = Connection(co, 'out1', va, 'in1', label='4')
c0 = Connection(va, 'out1', cc, 'in1', label='0')

# this line is crutial: you have to add all connections to your network
my_plant.add_conns(c1, c2, c3, c4, c0)
# %%[sec_5]
co.set_attr(pr=0.98, Q=-1e6)
ev.set_attr(pr=0.98)
cp.set_attr(eta_s=0.85)

c2.set_attr(T=20, x=1, fluid={'R134a': 1})
c4.set_attr(T=80, x=0)
# %%[sec_6]
my_plant.solve(mode='design')
my_plant.print_results()

print(f'COP = {abs(co.Q.val) / cp.P.val}')
# %%[sec_7]
co.set_attr(Q=None)
c1.set_attr(m=5)

my_plant.solve('design')
my_plant.print_results()
# %%[sec_8]
cp.set_attr(pr=4)
c4.set_attr(T=None)

my_plant.solve('design')
my_plant.print_results()
# %%[sec_9]
cp.set_attr(pr=None, eta_s=None)
c3.set_attr(T=97.3)
c4.set_attr(T=80)

my_plant.solve('design')
my_plant.print_results()
# %%[sec_10]
# first go back to the original state of the specifications
co.set_attr(Q=-1e6)
cp.set_attr(pr=None, eta_s=0.85)
c1.set_attr(m=None)
c3.set_attr(T=None)
c4.set_attr(T=80)

import matplotlib.pyplot as plt
import numpy as np

# make text reasonably sized
plt.rc('font', **{'size': 18})


data = {
    'T_source': np.linspace(0, 40, 11),
    'T_sink': np.linspace(60, 100, 11),
    'eta_s': np.linspace(0.75, 0.95, 11) * 100
}
COP = {
    'T_source': [],
    'T_sink': [],
    'eta_s': []
}
description = {
    'T_source': 'Evaporation temperature in °C',
    'T_sink': 'Condensation temperature in °C',
    'eta_s': 'Isentropic efficiency in %'
}

for T in data['T_source']:
    c2.set_attr(T=T)
    my_plant.solve('design')
    COP['T_source'] += [abs(co.Q.val) / cp.P.val]

# reset to base temperature
c2.set_attr(T=20)

for T in data['T_sink']:
    c4.set_attr(T=T)
    my_plant.solve('design')
    COP['T_sink'] += [abs(co.Q.val) / cp.P.val]

# reset to base temperature
c4.set_attr(T=80)

for eta_s in data['eta_s']:
    cp.set_attr(eta_s=eta_s / 100)
    my_plant.solve('design')
    COP['eta_s'] += [abs(co.Q.val) / cp.P.val]

fig, ax = plt.subplots(1, 3, sharey=True, figsize=(16, 8))

[a.grid() for a in ax]

i = 0
for key in data:
    ax[i].scatter(data[key], COP[key], s=100, color="#1f567d")
    ax[i].set_xlabel(description[key])
    i += 1

ax[0].set_ylabel('COP of the heat pump')

plt.tight_layout()

fig.savefig('heat_pump_parametric.svg')
# %%[sec_11]
