# %%[sec_1]
from tespy.components.basics.cycle_closer import CycleCloser
from tespy.networks import Network
from tespy.components import (
    CycleCloser, Pipe, Pump, Valve, SimpleHeatExchanger
)
from tespy.connections import Connection

nw = Network()
nw.set_attr(T_unit='C', p_unit='bar', h_unit='kJ / kg')

# central heating plant
hs = SimpleHeatExchanger('heat source')
cc = CycleCloser('cycle closer')
pu = Pump('feed pump')

# consumer
cons = SimpleHeatExchanger('consumer')
val = Valve('control valve')

# pipes
pipe_feed = Pipe('feed pipe')
pipe_return = Pipe('return pipe')

# connections
c0 = Connection(cc, "out1", hs, "in1", label="0")
c1 = Connection(hs, "out1", pu, "in1", label="1")
c2 = Connection(pu, "out1", pipe_feed, "in1", label="2")
c3 = Connection(pipe_feed, "out1", cons, "in1", label="3")
c4 = Connection(cons, "out1", val, "in1", label="4")
c5 = Connection(val, "out1", pipe_return, "in1", label="5")
c6 = Connection(pipe_return, "out1", cc, "in1", label="6")
nw.add_conns(c0, c1, c2, c3, c4, c5, c6)
# %%[sec_2]
cons.set_attr(Q=-10000, pr=0.98)
hs.set_attr(pr=1)
pu.set_attr(eta_s=0.75)
pipe_feed.set_attr(Q=-250, pr=0.98)
pipe_return.set_attr(Q=-200, pr=0.98)

c1.set_attr(T=90, p=10, fluid={'INCOMP::Water': 1})
c2.set_attr(p=13)
c4.set_attr(T=65)

nw.solve(mode="design")
nw.print_results()
# %%[sec_3]
pipe_feed.set_attr(
    ks=0.0005,  # pipe's roughness in meters
    L=100,  # length in m
    D="var",  # diameter in m
)
pipe_return.set_attr(
    ks=0.0005,  # pipe's roughness in meters
    L=100,  # length in m
    D="var",  # diameter in m
)
nw.solve(mode="design")
nw.print_results()
# %%[sec_4]
pipe_feed.set_attr(D=pipe_feed.D.val, pr=None)
pipe_return.set_attr(D=pipe_return.D.val, pr=None)
pipe_feed.set_attr(
    Tamb=0,  # ambient temperature level in network's temperature unit
    kA="var"  # area independent heat transfer coefficient
)
pipe_return.set_attr(
    Tamb=0,  # ambient temperature level in network's temperature unit
    kA="var"  # area independent heat transfer coefficient
)
nw.solve(mode="design")
nw.print_results()
# %%[sec_5]
nw.set_attr(iterinfo=False)
pipe_feed.set_attr(Tamb=0, kA=pipe_feed.kA.val, Q=None)
pipe_return.set_attr(Tamb=0, kA=pipe_return.kA.val, Q=None)

import matplotlib.pyplot as plt
import numpy as np

# make text reasonably sized
plt.rc('font', **{'size': 18})

data = {
    'T_ambient': np.linspace(-10, 20, 7),
    'heat_load': np.linspace(3, 12, 10),
    'T_level': np.linspace(90, 60, 7)
}
eta = {
    'T_ambient': [],
    'heat_load': [],
    'T_level': []
}
heat_loss = {
    'T_ambient': [],
    'heat_load': [],
    'T_level': []
}

for T in data['T_ambient']:
    pipe_feed.set_attr(Tamb=T)
    pipe_return.set_attr(Tamb=T)
    nw.solve('design')
    eta['T_ambient'] += [abs(cons.Q.val) / (hs.Q.val + pu.P.val) * 100]
    heat_loss['T_ambient'] += [abs(pipe_feed.Q.val + pipe_return.Q.val)]

# reset to base temperature
pipe_feed.set_attr(Tamb=0)
pipe_return.set_attr(Tamb=0)

for Q in data['heat_load']:
    cons.set_attr(Q=-1e3 * Q)
    nw.solve('design')
    eta['heat_load'] += [abs(cons.Q.val) / (hs.Q.val + pu.P.val) * 100]
    heat_loss['heat_load'] += [abs(pipe_feed.Q.val + pipe_return.Q.val)]

# reset to base temperature
cons.set_attr(Q=-10e3)

for T in data['T_level']:
    c1.set_attr(T=T)
    c4.set_attr(T=T - 20)  # return flow temperature assumed 20 °C lower than feed
    nw.solve('design')
    eta['T_level'] += [abs(cons.Q.val) / (hs.Q.val + pu.P.val) * 100]
    heat_loss['T_level'] += [abs(pipe_feed.Q.val + pipe_return.Q.val)]

fig, ax = plt.subplots(2, 3, figsize=(16, 8), sharex='col', sharey='row')

ax = ax.flatten()
[a.grid() for a in ax]

i = 0
for key in data:
    ax[i].scatter(data[key], eta[key], s=100, color="#1f567d")
    ax[i + 3].scatter(data[key], heat_loss[key], s=100, color="#18a999")
    i += 1

ax[0].set_ylabel('Efficiency in %')
ax[3].set_ylabel('Heat losses in W')
ax[3].set_xlabel('Ambient temperature in °C')
ax[4].set_xlabel('Consumer heat load in kW')
ax[5].set_xlabel('District heating temperature level in °C')
plt.tight_layout()
fig.savefig('district_heating_partload.svg')
plt.close()
# %%[sec_6]
