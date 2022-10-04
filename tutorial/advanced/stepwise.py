# %%[sec_1]
from tespy.networks import Network
working_fluid = "NH3"

nw = Network(
    fluids=["water", working_fluid],
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
c_in = Source("refrigerant in")
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
cd.set_attr(pr1=0.99, pr2=0.99)
rp.set_attr(eta_s=0.75)
cons.set_attr(pr=0.99)
# %%[sec_5]
from CoolProp.CoolProp import PropsSI as PSI
p_cond = PSI("P", "Q", 1, "T", 273.15 + 95, working_fluid) / 1e5
c0.set_attr(T=170, p=p_cond, fluid={"water": 0, working_fluid: 1})
c20.set_attr(T=60, p=2, fluid={"water": 1, working_fluid: 0})
c22.set_attr(T=90)

# key design paramter
cons.set_attr(Q=-230e3)
# %%[sec_6]
nw.solve("design")
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

c17 = Connection(amb_in, "out1", su, "in1", label="17")
c18 = Connection(su, "out1", ev, "in1", label="18")
c19 = Connection(ev, "out1", amb_out, "in1", label="19")

nw.add_conns(c17, c18, c19)
# %%[sec_9]
ev.set_attr(pr1=0.99)
su.set_attr(pr1=0.99, pr2=0.99)
# %%[sec_10]
# evaporator system cold side
p_evap = PSI("P", "Q", 1, "T", 273.15 + 5, working_fluid) / 1e5
c4.set_attr(x=0.9, p=p_evap)

h_sat = PSI("H", "Q", 1, "T", 273.15 + 15, working_fluid) / 1e3
c6.set_attr(h=h_sat)

# evaporator system hot side
c17.set_attr(T=15, fluid={"water": 1, working_fluid: 0})
c19.set_attr(T=9, p=1.013)
# %%[sec_11]
nw.solve("design")
nw.print_results()
# %%[sec_12]
from tespy.components import Compressor, Splitter, Merge

cp1 = Compressor("compressor 1")
cp2 = Compressor("compressor 2")

ic = HeatExchanger("intermittent cooling")
hsp = Pump("heat source pump")

sp = Splitter("splitter")
me = Merge("merge")
cv = Valve("control valve")

hs = Source("ambient intake")
cc = CycleCloser("heat pump cycle closer")
# %%[sec_13]
nw.del_conns(c0, c6, c17)

c6 = Connection(su, "out2", cp1, "in1", label="6")
c7 = Connection(cp1, "out1", ic, "in1", label="7")
c8 = Connection(ic, "out1", cp2, "in1", label="8")
c9 = Connection(cp2, "out1", cc, "in1", label="9")
c0 = Connection(cc, "out1", cd, "in1", label="0")

c11 = Connection(hs, "out1", hsp, "in1", label="11")
c12 = Connection(hsp, "out1", sp, "in1", label="12")
c13 = Connection(sp, "out1", ic, "in2", label="13")
c14 = Connection(ic, "out2", me, "in1", label="14")
c15 = Connection(sp, "out2", cv, "in1", label="15")
c16 = Connection(cv, "out1", me, "in2", label="16")
c17 = Connection(me, "out1", su, "in1", label="17")

nw.add_conns(c6, c7, c8, c9, c0, c11, c12, c13, c14, c15, c16, c17)
# %%[sec_14]
pr = (c1.p.val / c5.p.val) ** 0.5
cp1.set_attr(pr=pr)
ic.set_attr(pr1=0.99, pr2=0.98)
hsp.set_attr(eta_s=0.75)
# %%[sec_15]
c0.set_attr(p=p_cond, fluid={"water": 0, working_fluid: 1})

c6.set_attr(h=c5.h.val + 10)
c8.set_attr(h=c5.h.val + 10)

c7.set_attr(h=c5.h.val * 1.2)
c9.set_attr(h=c5.h.val * 1.2)

c11.set_attr(p=1.013, T=15, fluid={"water": 1, working_fluid: 0})
c14.set_attr(T=30)
# %% [sec_16]
nw.solve("design")

c0.set_attr(p=None)
cd.set_attr(ttd_u=5)

c4.set_attr(p=None)
ev.set_attr(ttd_l=5)

c6.set_attr(h=None)
su.set_attr(ttd_u=5)

c7.set_attr(h=None)
cp1.set_attr(eta_s=0.8)

c9.set_attr(h=None)
cp2.set_attr(eta_s=0.8)

c8.set_attr(h=None, Td_bp=4)
nw.solve("design")
nw.save("system_design")
# %% [sec_17]
cp1.set_attr(design=["eta_s"], offdesign=["eta_s_char"])
cp2.set_attr(design=["eta_s"], offdesign=["eta_s_char"])
rp.set_attr(design=["eta_s"], offdesign=["eta_s_char"])
hsp.set_attr(design=["eta_s"], offdesign=["eta_s_char"])

cons.set_attr(design=["pr"], offdesign=["zeta"])

cd.set_attr(
    design=["pr2", "ttd_u"], offdesign=["zeta2", "kA_char"]
)

from tespy.tools.characteristics import CharLine
from tespy.tools.characteristics import load_default_char as ldc

kA_char1 = ldc("heat exchanger", "kA_char1", "DEFAULT", CharLine)
kA_char2 = ldc("heat exchanger", "kA_char2", "EVAPORATING FLUID", CharLine)
ev.set_attr(
    kA_char1=kA_char1, kA_char2=kA_char2,
    design=["pr1", "ttd_l"], offdesign=["zeta1", "kA_char"]
)

su.set_attr(
    design=["pr1", "pr2", "ttd_u"], offdesign=["zeta1", "zeta2", "kA_char"]
)

ic.set_attr(
    design=["pr1", "pr2"], offdesign=["zeta1", "zeta2", "kA_char"]
)
c14.set_attr(design=["T"])
nw.solve("offdesign", design_path="system_design")
nw.print_results()
# %% [sec_18]
import numpy as np
nw.set_attr(iterinfo=False)

for Q in np.linspace(1, 0.6, 5) * cons.Q.val:
    cons.set_attr(Q=Q)
    nw.solve("offdesign", design_path="system_design")
    print(
        "COP:",
        abs(cons.Q.val) / (cp1.P.val + cp2.P.val + hsp.P.val + rp.P.val)
    )
# %% [sec_19]
