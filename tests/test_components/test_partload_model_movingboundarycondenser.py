from tespy.components import HeatExchanger, Source, Sink, Compressor, MovingBoundaryCondenser
from tespy.connections import Connection
from tespy.networks import Network
import numpy as np


nw = Network(T_unit="C", p_unit="bar")


so1 = Source("cw source")
so2 = Source("wf source")

# multiple-boundary heat exchanger
# allow for U value change as function of volumetric flow/mass flow/...
cd = MovingBoundaryCondenser("Condenser")
cp = Compressor("compressor")

si1 = Sink("cw sink")
si2 = Sink("wf sink")


c1 = Connection(so1, "out1", cd, "in2", label="1")
c2 = Connection(cd, "out2", si1, "in1", label="2")

c10 = Connection(so2, "out1", cp, "in1", label="10")
c11 = Connection(cp, "out1", cd, "in1", label="11")
c12 = Connection(cd, "out1", si2, "in1", label="12")

nw.add_conns(c1, c2, c10, c11, c12)

cd.set_attr(pr1=1, pr2=1)
cp.set_attr(eta_s=0.8)

c10.set_attr(T=20, x=1)
c11.set_attr(fluid={"NH3": 1})
c12.set_attr(x=0, T=80)

c1.set_attr(fluid={"INCOMP::Water": 1}, m=10, T=70, p=1, h0=1e5)
c2.set_attr(h0=1e5, T=80)

cd.set_attr(U_desup=4000, U_cond=16000, U_subcool=5000)

nw.solve("design")
nw.save("design")
nw.print_results()


# test pinch specification
c12.set_attr(T=None)
cd.set_attr(td_pinch=3)
nw.solve("design")
nw.print_results()
# exit()
# print(de.A_desup, de.A_cond)

# c12.set_attr(T=None)
cd.set_attr(A=cd.A.val, td_pinch=None)

# Alternative: fix the input temperatures and mass flows
# outlet conditions (condensing temperature and water outlet are unknows)
# c2.set_attr(T=None)
# c10.set_attr(m=c10.m.val)

#

# get rid of warnings
cd.zeta1.set_attr(min_val=-2)

nw.solve("design")
nw.print_results()
# exit()
# print(c2.T.val)


Q = []
T_cond = []
m_refrig = []
dT_pinch = []
Q_in_sections = []
for m in np.linspace(12, 5, 20):
    c1.set_attr(m=m)
    nw.solve("design")
    m_refrig += [c12.m.val]
    T_cond += [c12.T.val]
    Q += [abs(cd.Q.val)]
    dT_pinch += [cd.td_pinch.val]
    _, h_at_steps = cd.get_U_sections_and_h_steps()
    _, Q_in_section = cd.calc_td_log_and_Q_in_sections(h_at_steps)
    Q_in_sections += [Q_in_section]


from matplotlib import pyplot as plt


fig, ax = plt.subplots(4, sharex=True, figsize=(12, 8))

ax[0].scatter(Q, m_refrig)
ax[0].set_ylabel("refrigerant mass flow")
ax[1].scatter(Q, T_cond)
ax[1].set_ylabel("condensation temperature")
ax[2].scatter(Q, dT_pinch)
ax[2].set_ylabel("pinch temperature difference")
ax[3].scatter(Q, [abs(q[0] / Q) for q, Q in zip(Q_in_sections, Q)], label="desup")
ax[3].scatter(Q, [abs(q[1] / Q) for q, Q in zip(Q_in_sections, Q)], label="cond")
ax[3].legend()
ax[3].set_ylabel("heat transfer shares of total heat transfer")

ax[3].set_xlabel("total heat transfer")

[_.grid() for _ in ax]

plt.tight_layout()
fig.savefig("mb_U_sections_diff_m.png")

plt.close()

c1.set_attr(m=10)
nw.solve("design")
nw.print_results()

cd.set_attr(A=None, UA=cd.UA.val)

Q = []
T_cond = []
m_refrig = []
dT_pinch = []
Q_in_sections = []
for m in np.linspace(12, 5, 20):
    c1.set_attr(m=m)
    nw.solve("design")
    m_refrig += [c12.m.val]
    T_cond += [c12.T.val]
    Q += [abs(cd.Q.val)]
    # Q_cond += [cd.Q_cond]
    dT_pinch += [cd.td_pinch.val]
    _, h_at_steps = cd.get_U_sections_and_h_steps()
    _, Q_in_section = cd.calc_td_log_and_Q_in_sections(h_at_steps)
    Q_in_sections += [Q_in_section]

print(Q_in_sections)

fig, ax = plt.subplots(4, sharex=True, figsize=(12, 8))

ax[0].scatter(Q, m_refrig)
ax[0].set_ylabel("refrigerant mass flow")
ax[1].scatter(Q, T_cond)
ax[1].set_ylabel("condensation temperature")
ax[2].scatter(Q, dT_pinch)
ax[2].set_ylabel("pinch temperature difference")
ax[3].scatter(Q, [abs(q[0] / Q) for q, Q in zip(Q_in_sections, Q)], label="desup")
ax[3].scatter(Q, [abs(q[1] / Q) for q, Q in zip(Q_in_sections, Q)], label="cond")
ax[3].legend()
ax[3].set_ylabel("heat transfer shares of total heat transfer")

ax[3].set_xlabel("total heat transfer")

[_.grid() for _ in ax]

plt.tight_layout()
fig.savefig("mb_UA_diff_m.png")