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
c11.set_attr(fluid={"NH3": 1}, Td_bp=5, p=46)
c12.set_attr(Td_bp=-5, p=46)

c1.set_attr(fluid={"Water": 1}, m=1, p=0.3, Td_bp=-5, h0=1e5)
c2.set_attr(Td_bp=5, p=0.3)

# cd.set_attr(U_desup=4000, U_cond=16000, U_subcool=5000)

nw.set_attr(iterinfo=False)
nw.solve("design")

for Td_bp in np.linspace(5, 0.2):
    c11.set_attr(Td_bp=Td_bp)
    c12.set_attr(Td_bp=-Td_bp)
    nw.solve("design")
