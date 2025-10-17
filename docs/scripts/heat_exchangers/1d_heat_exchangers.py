# %% [markdown]
# # MovingBoundaryHeatExchanger (and SectionedHeatExchanger) class
#
# - 3 different categories of HeatExchangers
#
# 1. 0D models without secondary side
#    - SimpleHeatExchanger
#    - (ParabolicTrough, SolarCollector, Pipe)
# 2. 0D models with secondary side (no internal states know)
#    - HeatExchanger
#    - Desuperheater (x=1 at hot side outlet)
#    - Condenser (x=0 at hot side outlet, upper terminal temperature difference
#      and kA refer to hot side condensation temperature level)
# 3. 1D models (discretized internal states)
#    - SectionedHeatExchanger (cut into specified number of linear sections)
#    - MovingBoundaryHeatExchanger (cut into sections at phase boundaries)
#
# ## MovingBoundaryHeatExchanger
#
# Workflow
#
# 1. For hot and for cold side:
#    1. Check if phase changes between inlet and outlet
#    2. If yes, calculate enthalpy of phase change point (assumes pressure does
#       not change). Idea to include pressure drop at end.
#    3. Use enthalpies at inlets and outlets and phase change points of each side
#       to calculate heat exchanged from start to that point ($Q_i$)
#    4. Merge both heat exchanged lists from hot and cold side
#    5. Cut the heat exchanger into sections at each heat exchange step
# 2. Calculate internal temperatures
#    1. Calculate enthalpies at each step
#    $h_i = h_\text{in} + \frac{\dot Q_i}{\dot m}$, $h_i = h_\text{out} - \frac{\dot Q_i}{\dot m}$
#    1. Calculate temperature at each step $T_i=T(p, h_i)$
#    2. Calculate temperature difference at each step $\Delta T_i=T_{\text{hot,}i} - T_{\text{cold,}i}$
# 3. For UA only:
#    1. Calculate logarithmic temperature difference per section between adjacent steps $\Delta T_{\text{log,}k}=\frac{\Delta T_i - \Delta T_{i+1}}{\ln\frac{\Delta T_i}{\Delta T_{i+1}}}$
#    2. Calculate UA per section $UA_k=\frac{\dot Q_k}{\Delta T_{\text{log,}k}}$
#    3. Sum UA for overall $UA = \sum UA_k$
#

# %%
from tespy.components import MovingBoundaryHeatExchanger, Source, Sink
from tespy.connections import Connection
from tespy.networks import Network

# %%
nw = Network()
nw.units.set_defaults(
    temperature="°C",
    pressure="bar"
)

# %%
so1 = Source("source 1")
so2 = Source("source 2")

si1 = Sink("sink 1")
si2 = Sink("sink 2")

heatex = MovingBoundaryHeatExchanger("heatexchanger")

c1 = Connection(so1, "out1", heatex, "in1", label="c1")
c2 = Connection(heatex, "out1", si1, "in1", label="c2")
d1 = Connection(so2, "out1", heatex, "in2", label="d1")
d2 = Connection(heatex, "out2", si2, "in1", label="d2")

nw.add_conns(c1, c2, d1, d2)

# %% [markdown]
# - In principle, the MovingBoundaryHeatExchanger works exactly as the 0D model

# %%
c1.set_attr(fluid={"R290": 1}, td_dew=50, T_dew=60, m=5)
c2.set_attr(td_bubble=5)
d1.set_attr(fluid={"water": 1}, p=1, T=45)
d2.set_attr(T=55)

heatex.set_attr(dp1=0, dp2=0)

nw.solve("design")

# %% [markdown]
# - The component has internal pinch and overall UA based on identified sections
#   on top of final (upper ttd) and initial (lower ttd) pinch

# %%
nw.results["MovingBoundaryHeatExchanger"][["td_pinch", "UA", "ttd_u", "ttd_l", "kA"]]

# %% [markdown]
# - We can extract the sections: temperatures and heat transfer and make a QT plot

# %%
heat, T_hot, T_cold, heat_per_section, td_log_per_section = heatex.calc_sections()

# %%
from matplotlib import pyplot as plt


fig, ax = plt.subplots(1)

ax.plot((heat, heat),([i for (i,j) in zip(T_hot,list(T_cold))], [j for (i,j) in zip(T_hot,list(T_cold))]),color='grey', linestyle = "--")
ax.plot(heat, T_hot, "o-", color="red")
ax.plot(heat, T_cold, "o-", color="blue")
ax.set_xlabel(r'$\Delta\dot{H}$ [W]')
ax.set_ylabel(r'T [K]')


# %% [markdown]
# - pinches are T_hot - T_cold

# %%
T_hot - T_cold

# %%
c1.set_attr(T_dew=None)
heatex.set_attr(td_pinch=2.5)
nw.solve("design")

heat, T_hot, T_cold, heat_per_section, td_log_per_section = heatex.calc_sections()

fig, ax = plt.subplots(1)

ax.plot((heat, heat),([i for (i,j) in zip(T_hot,list(T_cold))], [j for (i,j) in zip(T_hot,list(T_cold))]),color='grey', linestyle = "--")
ax.plot(heat, T_hot, "o-", color="red")
ax.plot(heat, T_cold, "o-", color="blue")
ax.set_xlabel(r'$\Delta\dot{H}$ [W]')
ax.set_ylabel(r'T [K]')

T_hot - T_cold

# %% [markdown]
# - kA references the terminal temperatures and calculates td_log based on those
# - UA references the terminal temperatures of each section and sums the UA of
#   all sections to total UA

# %%
heatex.kA.val, heatex.UA.val

# %% [markdown]
# Fix UA for pinch calculation is also possible

heatex.set_attr(td_pinch=None, UA=350e3)
nw.solve("design")
nw.save("design.json")

heat, T_hot, T_cold, heat_per_section, td_log_per_section = heatex.calc_sections()

fig, ax = plt.subplots(1)

ax.plot((heat, heat),([i for (i,j) in zip(T_hot,list(T_cold))], [j for (i,j) in zip(T_hot,list(T_cold))]),color='grey', linestyle = "--")
ax.plot(heat, T_hot, "o-", color="red")
ax.plot(heat, T_cold, "o-", color="blue")
ax.set_xlabel(r'$\Delta\dot{H}$ [W]')
ax.set_ylabel(r'T [K]')

T_hot - T_cold

import numpy as np


heatex.set_attr(UA=None, alpha_ratio=1, area_ratio=1, refrigerant_index=0, re_exp_sf=0.55, re_exp_r=0.7)

m_range = np.linspace(5, 2, 11)
UA_offdesign = []
for m in m_range:
    c1.set_attr(m=m)
    nw.solve("offdesign", design_path="design.json")
    UA_offdesign.append(heatex.UA.val)


fig, ax = plt.subplots(1)

ax.scatter(m_range, UA_offdesign)

from tespy.components import SectionedHeatExchanger

nw = Network()
nw.units.set_defaults(
    temperature="°C",
    pressure="bar"
)

so1 = Source("source 1")
so2 = Source("source 2")

si1 = Sink("sink 1")
si2 = Sink("sink 2")

heatex = SectionedHeatExchanger("heatexchanger")

c1 = Connection(so1, "out1", heatex, "in1", label="c1")
c2 = Connection(heatex, "out1", si1, "in1", label="c2")
d1 = Connection(so2, "out1", heatex, "in2", label="d1")
d2 = Connection(heatex, "out2", si2, "in1", label="d2")

nw.add_conns(c1, c2, d1, d2)

c1.set_attr(fluid={"R290": 1}, td_dew=50, T_dew=60, m=5)
c2.set_attr(td_bubble=5)
d1.set_attr(fluid={"water": 1}, p=1, T=45)
d2.set_attr(T=55)

heatex.set_attr(dp1=1, dp2=0)

nw.solve("design")

heat, T_hot, T_cold, heat_per_section, td_log_per_section = heatex.calc_sections()

fig, ax = plt.subplots(1)

ax.plot((heat, heat),([i for (i,j) in zip(T_hot,list(T_cold))], [j for (i,j) in zip(T_hot,list(T_cold))]),color='grey', linestyle = "-", linewidth=0.5)
ax.plot(heat, T_hot, "o-", color="red")
ax.plot(heat, T_cold, "o-", color="blue")
ax.set_xlabel(r'$\Delta\dot{H}$ [W]')
ax.set_ylabel(r'T [K]')

heat_defaultsteps = heat # for later comparison
T_hot_defaultsteps = T_hot

heatex.set_attr(num_sections=1)
nw.solve("design")

heat, T_hot, T_cold, heat_per_section, td_log_per_section = heatex.calc_sections()

fig, ax = plt.subplots(1)

ax.plot((heat, heat),([i for (i,j) in zip(T_hot,list(T_cold))], [j for (i,j) in zip(T_hot,list(T_cold))]),color='grey', linestyle = "-", linewidth=0.5)
ax.plot(heat_defaultsteps, T_hot_defaultsteps, "o-", color="grey", markersize=2, linewidth = 0.5, label = "sectioned heat exchanger (50 sections)")
ax.plot(heat, T_hot, "o-", color="red", label = "0D Heat Exchanger")
ax.plot(heat, T_cold, "o-", color="blue")
ax.set_xlabel(r'$\Delta\dot{H}$ [W]')
ax.set_ylabel(r'T [K]')
ax.legend()

plt.show()


heatex.set_attr(num_sections=6)
nw.solve("design")

heat, T_hot, T_cold, heat_per_section, td_log_per_section = heatex.calc_sections()

fig, ax = plt.subplots(1)

ax.plot((heat, heat),([i for (i,j) in zip(T_hot,list(T_cold))], [j for (i,j) in zip(T_hot,list(T_cold))]),color='grey', linestyle = "-", linewidth=0.5)
ax.plot(heat_defaultsteps, T_hot_defaultsteps, "o-", color="grey", markersize=2, linewidth = 0.5, label = "sectioned heat exchanger (50 sections)")
ax.plot(heat, T_hot, "o-", color="red", label = "sectioned heat exchanger (6 sections)")
ax.plot(heat, T_cold, "o-", color="blue")
ax.set_xlabel(r'$\Delta\dot{H}$ [W]')
ax.set_ylabel(r'T [K]')
ax.legend()

plt.show()
