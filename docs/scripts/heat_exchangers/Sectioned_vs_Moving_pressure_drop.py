from matplotlib import pyplot as plt
from tespy.components import SectionedHeatExchanger, MovingBoundaryHeatExchanger, Source, Sink
from tespy.connections import Connection
from tespy.networks import Network


nw = Network()
nw.units.set_defaults(
    temperature="°C",
    pressure="bar",
    heat="MW"
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

heatex.set_attr(dp1=2, dp2=0.5)

nw.solve("design")

heat, T_hot, T_cold, heat_per_section, td_log_per_section = heatex.calc_sections()
heat /= 1e6

fig, ax = plt.subplots(1, figsize=(10, 6))

annotation_color = "black"

heat_defaultsteps = heat # for later comparison
T_hot_defaultsteps = T_hot

# new Network with MovingBoundary

nw = Network()
nw.units.set_defaults(
    temperature="°C",
    pressure="bar",
    heat="MW"
)

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

c1.set_attr(fluid={"R290": 1}, td_dew=50, T_dew=60, m=5)
c2.set_attr(td_bubble=5)
d1.set_attr(fluid={"water": 1}, p=1, T=45)
d2.set_attr(T=55)

heatex.set_attr(dp1=2, dp2=0.5)

nw.solve("design")

heat, T_hot, T_cold, heat_per_section, td_log_per_section = heatex.calc_sections()
heat /= 1e6

fig, ax = plt.subplots(1, figsize=(10, 6))

ax.plot((heat, heat), ([T for T in T_hot], [T for T in T_cold]), color=annotation_color, linestyle="-", linewidth=0.5)
ax.plot(heat_defaultsteps, T_hot_defaultsteps, "o-", color=annotation_color, markersize=2, linewidth=0.5, label="SectionedHeatExchanger (50 sections)")
ax.plot(heat, T_hot, "o-", color="red", label="MovingBoundaryHeatExchanger")
ax.plot(heat, T_cold, "o-", color="blue")
ax.set_ylabel("temperature in K")
ax.set_xlabel("heat transferred in MW")
ax.legend()

fig.savefig("Sectioned_vs_Moving_pressure_drop.svg", bbox_inches="tight")
