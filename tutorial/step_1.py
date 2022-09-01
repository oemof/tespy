# %%[sec_1]
from tespy.networks import Network

# %% network

nw = Network(
    fluids=["water", "NH3"],
    T_unit="C", p_unit="bar", h_unit="kJ / kg", m_unit="kg / s"
)
# %%[sec_2]
from tespy.components import Condenser
from tespy.components import CycleCloser
from tespy.components import HeatExchangerSimple
from tespy.components import Pump
from tespy.components import Sink
from tespy.components import Source

# sources & sinks
c_in = Source("coolant in")
cons_closer = CycleCloser("consumer cycle closer")
va = Sink("valve")

# consumer system
cd = Condenser("condenser")
rp = Pump("recirculation pump")
cons = HeatExchangerSimple("consumer")
# %%[sec_3]
from tespy.connections import Connection

c0 = Connection(c_in, "out1", cd, "in1", label="0")
c1 = Connection(cd, "out1", va, "in1", label="1")

c20 = Connection(cons_closer, "out1", rp, "in1", label="20")
c21 = Connection(rp, "out1", cd, "in2", label="21")
c22 = Connection(cd, "out2", cons, "in1", label="22")
c23 = Connection(cons, "out1", cons_closer, "in1", label="23")

nw.add_conns(c0, c1, c20, c21, c22, c23)
# %%[sec_4]
cd.set_attr(
    pr1=0.99, pr2=0.99, ttd_u=5,
    design=["pr2", "ttd_u"], offdesign=["zeta2", "kA_char"]
)
rp.set_attr(eta_s=0.8, design=["eta_s"], offdesign=["eta_s_char"])
cons.set_attr(pr=0.99, design=["pr"], offdesign=["zeta"])
# %%[sec_5]
c0.set_attr(T=170, fluid={"water": 0, "NH3": 1})
c20.set_attr(T=60, p=2, fluid={"water": 1, "NH3": 0})
c22.set_attr(T=90)

# key design paramter
cons.set_attr(Q=-230e3)
# %%[sec_6]
nw.solve("design")
nw.print_results()
nw.save("condenser")

cons.set_attr(Q=-200e3)

nw.solve("offdesign", design_path="condenser")
nw.print_results()
# %%[sec_7]
