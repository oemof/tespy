# %%[sec_1]
from tespy.networks import Network

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
nw.save("consumer_part")

# run an offdesign test: without modifying any value it should give you
# the same results as the design calculation
nw.solve("offdesign", design_path="consumer_part")
nw.print_results()
# %%[sec_7]
from tespy.components import Valve, Drum, HeatExchanger

# ambient heat source
amb_in = Source("source ambient")
amb_out = Sink("sink ambient")

# evaporator system
va = Valve("valve")
dr = Drum("drum")
ev = HeatExchanger("evaporator")
su = HeatExchanger("superheater")

# virtual source
cp1 = Sink("compressor 1")
# %%[sec_8]
nw.del_conns(c1)

# evaporator system
c1 = Connection(cd, "out1", va, "in1", label="1")
c2 = Connection(va, "out1", dr, "in1", label="2")
c3 = Connection(dr, "out1", ev, "in2", label="3")
c4 = Connection(ev, "out2", dr, "in2", label="4")
c5 = Connection(dr, "out2", su, "in2", label="5")
c6 = Connection(su, "out2", cp1, "in1", label="6")

nw.add_conns(c1, c2, c3, c4, c5, c6)

c13 = Connection(amb_in, "out1", su, "in1", label="13")
c14 = Connection(su, "out1", ev, "in1", label="14")
c15 = Connection(ev, "out1", amb_out, "in1", label="15")

nw.add_conns(c13, c14, c15)
# %%[sec_9]
from tespy.tools.characteristics import CharLine
from tespy.tools.characteristics import load_default_char as ldc
# evaporator system

kA_char1 = ldc('heat exchanger', 'kA_char1', 'DEFAULT', CharLine)
kA_char2 = ldc('heat exchanger', 'kA_char2', 'EVAPORATING FLUID', CharLine)

ev.set_attr(
    pr1=0.99, ttd_l=5, kA_char1=kA_char1, kA_char2=kA_char2,
    design=['pr1', 'ttd_l'], offdesign=['zeta1', 'kA_char']
)
su.set_attr(
    pr1=0.99, pr2=0.99, ttd_u=2,
    design=['pr1', 'pr2', 'ttd_u'], offdesign=['zeta1', 'zeta2', 'kA_char']
)
# %%[sec_10]
from tespy.connections import Ref
# evaporator system cold side
c3.set_attr(m=Ref(c2, 0.75, 0))
# we provide this keyword for numerical stability
c6.set_attr(state="g")

# evaporator system hot side
c13.set_attr(T=15, fluid={'water': 1, 'NH3': 0})
c15.set_attr(T=9, p=1.013)
# %%[sec_11]
nw.solve("design")
nw.print_results()
nw.save("evaporator_part")

# run the offdesign test
nw.solve("offdesign", design_path="evaporator_part")
nw.print_results()
# %%[sec_12]

