from tespy.components import HeatExchanger, Source, Sink, Compressor, MovingBoundaryHeatExchanger
from tespy.connections import Connection
from tespy.networks import Network
import numpy as np


nw = Network(T_unit="C", p_unit="bar")


so1 = Source("cw source")
so2 = Source("wf source")

# multiple-boundary heat exchanger
# allow for U value change as function of volumetric flow/mass flow/...
cd = MovingBoundaryHeatExchanger("Condenser")
cp = Compressor("compressor")

si1 = Sink("cw sink")
si2 = Sink("wf sink")


c1 = Connection(so1, "out1", cd, "in2", label="1")
c2 = Connection(cd, "out2", si1, "in1", label="2")

# c10 = Connection(so2, "out1", cp, "in1", label="10")
c11 = Connection(so2, "out1", cd, "in1", label="11")
c12 = Connection(cd, "out1", si2, "in1", label="12")

nw.add_conns(c1, c2, c11, c12)

cp.set_attr(eta_s=0.8)

# c10.set_attr(T=20, x=1)
c11.set_attr(fluid={"NH3": 1}, m=10, Td_bp=40, p=40)
c12.set_attr(x=0, p=40)

c1.set_attr(fluid={"Water": 1}, p=1, Td_bp=-5, h0=1e5)
c2.set_attr(p=1, Td_bp=20)


# cd.set_attr(U_desup=4000, U_cond=16000, U_subcool=5000)

nw.set_attr(iterinfo=False)
nw.solve("design")

c11.set_attr(p=None)
c12.set_attr(p=None)
cd.set_attr(td_pinch=5, pr1=1)
nw.solve("design")
c11.set_attr(Td_bp=30.4)
# c12.set_attr(Td_bp=-Td_bp)
nw.solve("design")

# Q_sections = cd._assign_sections()
# T_steps_hot, T_steps_cold = cd._get_T_at_steps(Q_sections)
# print(cd.UA.val, cd.td_pinch.val, Q_sections, T_steps_hot, T_steps_cold)
# print(cd._calc_UA_in_sections(T_steps_hot, T_steps_cold, np.diff(Q_sections)))

# c11.set_attr(Td_bp=30.5)
# # c12.set_attr(Td_bp=-Td_bp)
# nw.solve("design")

# # Q_sections = cd._assign_sections()
# # print(cd.UA.val, cd.td_pinch.val, Q_sections, cd._get_T_at_steps(Q_sections))

# c11.set_attr(Td_bp=30.4)
# # c12.set_attr(Td_bp=-Td_bp)
# nw.solve("design")

# Q_sections2 = cd._assign_sections()
# T_steps_hot2, T_steps_cold2 = cd._get_T_at_steps(Q_sections2)
# print(cd.UA.val, cd.td_pinch.val, Q_sections, T_steps_hot2, T_steps_cold2)
# print(cd._calc_UA_in_sections(T_steps_hot2, T_steps_cold2, np.diff(Q_sections2)))

# exit()

for Td_bp in np.linspace(40, 2):
    c11.set_attr(Td_bp=Td_bp)
    # c12.set_attr(Td_bp=-Td_bp)
    nw.solve("design")
    # print(cd.td_pinch.val, cd.UA.val)
    Q_sections = cd._assign_sections()
    print(Td_bp, cd.UA.val)
