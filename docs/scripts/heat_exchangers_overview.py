from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from tespy.components import HeatExchanger, Source, Sink
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

heatex = HeatExchanger("heatexchanger")

c1 = Connection(so1, "out1", heatex, "in1", label="c1")
c2 = Connection(heatex, "out1", si1, "in1", label="c2")
d1 = Connection(so2, "out1", heatex, "in2", label="d1")
d2 = Connection(heatex, "out2", si2, "in1", label="d2")

nw.add_conns(c1, c2, d1, d2)

c1.set_attr(fluid={"R290": 1}, td_dew=50, T_dew=60, m=5)
c2.set_attr(td_bubble=5)
d1.set_attr(fluid={"water": 1}, p=1, T=45)
d2.set_attr(T=55)

heatex.set_attr(dp1=0, dp2=0)

nw.solve("design")


heat = [0, abs(heatex.Q.val)]
T_hot = [c2.T.val_SI, c1.T.val_SI]
T_cold = [d1.T.val_SI, d2.T.val_SI]


fig, ax = plt.subplots(1)

ax.plot(heat, T_hot, "o-", color="red")
ax.plot(heat, T_cold, "o-", color="blue")

annotation_color = "black"

ttd_u_location = (heat[1] * (1 + 1 / 15), T_cold[1]), (heat[1] * (1 + 1 / 15), T_hot[1])
bracket = mpatches.FancyArrowPatch(
    *ttd_u_location,
    connectionstyle=f"bar,angle=0,fraction=0",
    arrowstyle='-',                      # no arrowhead
    linewidth=2, color=annotation_color
)
ax.add_patch(bracket)

offset = heat[1] / 15     # absolute horizontal offset from the line
bar_len = heat[1] / 15    # absolute size of the bracket ends (in data units)
ax.plot([heat[1] + offset - bar_len/2, heat[1] + offset + bar_len/2], [T_hot[1], T_hot[1]], color=annotation_color, lw=2)
ax.plot([heat[1] + offset - bar_len/2, heat[1] + offset + bar_len/2], [T_cold[1], T_cold[1]], color=annotation_color, lw=2)

ax.text(heat[1] * 1.1, (T_hot[1] + T_cold[1]) / 2, r'$\Delta T = \text{ttd_u}$', va='center', color=annotation_color, rotation=0)

ttd_l_location = (0 - offset, T_cold[0]), (0 - offset, T_hot[0])
bracket = mpatches.FancyArrowPatch(
    *ttd_l_location,
    connectionstyle=f"bar,angle=0,fraction=0",
    arrowstyle='-',                      # no arrowhead
    linewidth=2, color=annotation_color
)
ax.add_patch(bracket)

ax.text(heat[0] - offset * 5, (T_hot[0] + T_cold[0]) / 2, r'$\Delta T = \text{ttd_l}$', va='center', color=annotation_color, rotation=0)
ax.plot([heat[0] - offset - bar_len/2, heat[0] - offset + bar_len/2], [T_hot[0], T_hot[0]], color=annotation_color, lw=2)
ax.plot([heat[0] - offset - bar_len/2, heat[0] - offset + bar_len/2], [T_cold[0], T_cold[0]], color=annotation_color, lw=2)
ax.set_xbound([- 6 * offset, heat[1] + 6 * offset])
ax.set_ylabel("temperature in K")
ax.set_xlabel("heat transferred in MW")

fig.savefig("HeatExchanger.svg", bbox_inches="tight")
