# %%[sec_1]
from tespy.networks import Network

nw = Network(
    T_unit="C", p_unit="bar", h_unit="kJ / kg", m_unit="t / h"
)
# %%[sec_2]
from tespy.components import Merge
from tespy.components import Pump
from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Splitter
from tespy.components import Turbine
from tespy.components import Valve

# sources & sinks
live_steam = Source("live steam")
makeup_water = Source("makeup water")
process_steam = Sink("process steam")
process_condensate = Source("process condensate")
pre_feed_water = Sink("pre_feed_water")

# other components
extraction_turbine = Turbine("extraction turbine")
condensing_turbine = Turbine("condensing turbine")
condenser = SimpleHeatExchanger("condenser")
steam_cond_pump = Pump("steam conditioner pump")
condensate_pump = Pump("condensate pump")

merge_1 = Merge("merge 1")
merge_2 = Merge("merge 2", num_in=3)
conditioning_valve = Valve("conditioning valve")
splitter_1 = Splitter("steam extraction")
splitter_2 = Splitter("condensate split")
# %%[sec_3]
from tespy.connections import Connection


c1 = Connection(live_steam, "out1", extraction_turbine, "in1", label="1")
c2 = Connection(extraction_turbine, "out1", splitter_1, "in1", label="2")
c3 = Connection(splitter_1, "out2", merge_1, "in2", label="3")

c4 = Connection(splitter_1, "out1", conditioning_valve, "in1", label="4")
c5 = Connection(conditioning_valve, "out1", condensing_turbine, "in1", label="5")
c6 = Connection(condensing_turbine, "out1", condenser, "in1", label="6")
c7 = Connection(condenser, "out1", splitter_2, "in1", label="7")
c8 = Connection(splitter_2, "out1", condensate_pump, "in1", label="8")

c9 = Connection(splitter_2, "out2", steam_cond_pump, "in1", label="9")
c10 = Connection(steam_cond_pump,"out1", merge_1,"in1", label="10")
c11 = Connection(merge_1, "out1",process_steam,"in1", label="11")

c12 = Connection(condensate_pump,"out1", merge_2, "in1", label="12")
c13 = Connection(process_condensate, "out1", merge_2, "in2", label="13")
c14 = Connection(makeup_water,"out1", merge_2, "in3", label="14")
c15 = Connection(merge_2, "out1", pre_feed_water, "in1", label="15")

nw.add_conns(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15)
# %%[sec_4]
from tespy.connections import Ref


Tamb = 15
pamb = 1.013

c1.set_attr(fluid={"water": 1}, T=480, p=88.5, m=200)

extraction_mass_flow = 160 # t/h
c3.set_attr(m=extraction_mass_flow)
c7.set_attr(x=0, T=Tamb + 18)

process_steam_pressure = 13
process_steam_superheating = 10
c11.set_attr(Td_bp=process_steam_superheating, p=process_steam_pressure)
c12.set_attr(p=pamb)

condensate_return_fraction = 0.75
c13.set_attr(fluid={"water": 1}, x=0, m=Ref(c11, condensate_return_fraction, 0))
c14.set_attr(fluid={"water": 1}, T=Tamb)
c15.set_attr(m=200)
# %%[sec_5]
extraction_turbine.set_attr(eta_s=0.9)
condensing_turbine.set_attr(eta_s=0.9)

condenser.set_attr(pr=1)
conditioning_valve.set_attr(pr=0.85)
condensate_pump.set_attr(eta_s=0.7)
steam_cond_pump.set_attr(eta_s=0.7)

nw.solve("design")

c3.set_attr(m=None)
c11.set_attr(m=extraction_mass_flow)

nw.solve("design")

nw.print_results()
# %%[sec_6]
from tespy.components import Drum
from tespy.components import HeatExchanger


nw.del_conns(c15)

steam_deaeration = Sink("steam deaeration")
to_hp_cond = Sink("to hp condenser")

lp_feed_pump = Pump("lp feed pump")
hp_feed_pump = Pump("hp feed pump")

lp_drum = Drum("low pressure drum")
splitter_3 = Splitter("splitter 3")
splitter_4 = Splitter("splitter 4")
merge_3 = Merge("merge 3")

lp_evaporator = SimpleHeatExchanger("low pressure evaporator")
feed_water_preheater = HeatExchanger("feed water preheater")
hp_eco_1 = SimpleHeatExchanger("high pressure economizer 1")
hp_eco_2 = SimpleHeatExchanger("high pressure economizer 2")

c15 = Connection(merge_2, "out1", lp_feed_pump, "in1", label="15")
c16 = Connection(lp_feed_pump, "out1", feed_water_preheater, "in2", label="16")
c17 = Connection(feed_water_preheater, "out2", lp_drum, "in1", label="17")
c18 = Connection(lp_drum, "out1", splitter_3, "in1", label="18")

c19 = Connection(splitter_3, "out1", lp_evaporator, "in1", label="19")
c20 = Connection(lp_evaporator, "out1", lp_drum, "in2", label="20")
c21 = Connection(lp_drum, "out2", steam_deaeration, "in1", label="21")

c22 = Connection(splitter_3, "out2", feed_water_preheater, "in1", label="22")
c23 = Connection(feed_water_preheater, "out1", hp_feed_pump, "in1", label="23")
c24 = Connection(hp_feed_pump, "out1", splitter_4, "in1", label="24")

c25 = Connection(splitter_4, "out1", hp_eco_1, "in1", label="25")
c26 = Connection(hp_eco_1, "out1", hp_eco_2, "in1", label="26")
c27 = Connection(hp_eco_2, "out1", merge_3, "in1", label="27")

c28 = Connection(splitter_4, "out2", merge_3, "in2", label="28")
c29 = Connection(merge_3, "out1", to_hp_cond, "in1", label="29")

nw.add_conns(c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27, c28, c29)
# %%[sec_7]
c16.set_attr(p=2.3, T=55)
c17.set_attr(T=100, p0=2.3)
c20.set_attr(x=.05)
c21.set_attr(m=0.1)

c24.set_attr(p=88.5)
c26.set_attr(T=116, p0=88.5)
c27.set_attr(T=231, p0=88.5)
c28.set_attr(m=0)

lp_feed_pump.set_attr(eta_s=0.7)
hp_feed_pump.set_attr(eta_s=0.7)

feed_water_preheater.set_attr(pr1=0.9, pr2=0.9)

hp_eco_2.set_attr(pr=1)

nw.solve("design")
nw.print_results()
# %%[sec_8]
nw.del_conns(c29, c1)

hp_blowdown = Sink("hp blowdown")

hp_drum = Drum("high pressure drum")
splitter_5 = Splitter("splitter 5")
splitter_6 = Splitter("splitter 6")
merge_4 = Merge("merge 4")

hp_condenser = HeatExchanger("hp condenser")
hp_eco_3 = SimpleHeatExchanger("high pressure economizer 3")
hp_evaporator = SimpleHeatExchanger("high pressure evaporator")
hp_s_1 = SimpleHeatExchanger("high pressure superheater 1")
hp_s_2 = SimpleHeatExchanger("high pressure superheater 2")

c29 = Connection(merge_3, "out1", hp_condenser, "in2", label="29")
c30 = Connection(hp_condenser, "out2", hp_eco_3, "in1", label="30")
c31 = Connection(hp_eco_3, "out1", hp_drum, "in1", label="31")
c32 = Connection(hp_drum, "out1", splitter_5, "in1", label="32")

c33 = Connection(splitter_5, "out1", hp_blowdown, "in1", label="33")

c34 = Connection(splitter_5, "out2", hp_evaporator, "in1", label="34")
c35 = Connection(hp_evaporator, "out1", hp_drum, "in2", label="35")

c36 = Connection(hp_drum, "out2", splitter_6, "in1", label="36")

c37 = Connection(splitter_6, "out1", hp_condenser, "in1", label="37")
c38 = Connection(hp_condenser, "out1", merge_4, "in1", label="38")

c39 = Connection(splitter_6, "out2", hp_s_1, "in1", label="39")
c40 = Connection(hp_s_1, "out1", merge_4, "in2", label="40")
c41 = Connection(merge_4, "out1", hp_s_2, "in1", label="41")

c42 = Connection(hp_s_2, "out1", extraction_turbine, "in1", label="42")

nw.add_conns(c29, c30, c31, c32, c33, c34, c35, c36, c37, c38, c39, c40, c41, c42)
# %%[sec_9]
from tespy.connections import Bus


c31.set_attr(T=290, p0=88.5)
c33.set_attr(m=2)
c35.set_attr(x=0.05)
c38.set_attr(x=0)
c40.set_attr(T=529)
c41.set_attr(T=429)
c42.set_attr(fluid={"water": 1}, T=480)

hp_condenser.set_attr(pr1=1, pr2=1)
hp_eco_3.set_attr(pr=1)
hp_s_2.set_attr(pr=1)

power_output = Bus("power output")
power_output.add_comps(
    {"comp": extraction_turbine, "base": "component", "char": 0.97},
    {"comp": condensing_turbine, "base": "component", "char": 0.97},
    {"comp": condensate_pump, "base": "bus", "char": 0.97},
    {"comp": steam_cond_pump, "base": "bus", "char": 0.97},
    {"comp": lp_feed_pump, "base": "bus", "char": 0.97},
    {"comp": hp_feed_pump, "base": "bus", "char": 0.97}
)
nw.add_busses(power_output)

nw.solve("design")
nw.print_results()
# %% [sec_10]
nw.del_conns(c19, c20, c25, c26, c27, c30, c31, c34, c35, c39, c40, c41, c42)

lp_evaporator = HeatExchanger("low pressure evaporator")
hp_eco_1 = HeatExchanger("high pressure economizer 1")
hp_eco_2 = HeatExchanger("high pressure economizer 2")
hp_eco_3 = HeatExchanger("high pressure economizer 3")
hp_evaporator = HeatExchanger("high pressure evaporator")
hp_s_1 = HeatExchanger("high pressure superheater 1")
hp_s_2 = HeatExchanger("high pressure superheater 2")

c19 = Connection(splitter_3, "out1", lp_evaporator, "in2", label="19")
c20 = Connection(lp_evaporator, "out2", lp_drum, "in2", label="20")

c25 = Connection(splitter_4, "out1", hp_eco_1, "in2", label="25")
c26 = Connection(hp_eco_1, "out2", hp_eco_2, "in2", label="26")
c27 = Connection(hp_eco_2, "out2", merge_3, "in1", label="27")

c30 = Connection(hp_condenser, "out2", hp_eco_3, "in2", label="30")
c31 = Connection(hp_eco_3, "out2", hp_drum, "in1", label="31")

c34 = Connection(splitter_5, "out2", hp_evaporator, "in2", label="34")
c35 = Connection(hp_evaporator, "out2", hp_drum, "in2", label="35")

c39 = Connection(splitter_6, "out2", hp_s_1, "in2", label="39")
c40 = Connection(hp_s_1, "out2", merge_4, "in2", label="40")
c41 = Connection(merge_4, "out1", hp_s_2, "in2", label="41")

c42 = Connection(hp_s_2, "out2", extraction_turbine, "in1", label="42")
nw.add_conns(c19, c20, c25, c26, c27, c30, c31, c34, c35, c39, c40, c41, c42)

c20.set_attr(x=.05)
c26.set_attr(T=116)
c27.set_attr(T=231, p0=88.5)
c31.set_attr(T=290, p0=88.5)
c33.set_attr(m=2)
c35.set_attr(x=0.05)
c38.set_attr(x=0)
c40.set_attr(T=529)
c41.set_attr(T=429)
c42.set_attr(fluid={"water": 1}, T=480)

gt_exhaust = Source("gas turbine exhaust")
chimney = Sink("chimney")

c43 = Connection(gt_exhaust, "out1", hp_s_2, "in1", label="43")
c44 = Connection(hp_s_2, "out1", hp_s_1, "in1", label="44")
c45 = Connection(hp_s_1, "out1", hp_evaporator, "in1", label="45")
c46 = Connection(hp_evaporator, "out1", hp_eco_3, "in1", label="46")
c47 = Connection(hp_eco_3, "out1", hp_eco_2, "in1", label="47")
c48 = Connection(hp_eco_2, "out1", lp_evaporator, "in1", label="48")
c49 = Connection(lp_evaporator, "out1", hp_eco_1, "in1", label="49")
c50 = Connection(hp_eco_1, "out1", chimney, "in1", label="50")
nw.add_conns(c43, c44, c45, c46, c47, c48, c49, c50)

x_gh_CO2 = 0.061
x_gh_H2O = 0.07768
x_gh_N2 = 0.7364
x_gh_O2 = 0.1126
x_gh_Ar = 1 - x_gh_CO2 - x_gh_H2O - x_gh_N2 - x_gh_O2

T_gh = 713
G_gh = 4 * 482
c43.set_attr(
    fluid={"O2": x_gh_O2, "N2": x_gh_N2, "CO2": x_gh_CO2, "water": x_gh_H2O, "Ar": x_gh_Ar},
    m=G_gh,
    T=T_gh,
    p=pamb,
)

lp_evaporator.set_attr(pr1=1)
hp_eco_1.set_attr(pr1=1)
hp_eco_2.set_attr(pr1=1, pr2=1)
hp_eco_3.set_attr(pr1=1, pr2=1)

hp_evaporator.set_attr(pr1=1)
hp_s_1.set_attr(pr1=1)
hp_s_2.set_attr(pr1=1, pr2=1)

nw.solve("design")
nw.print_results()
# %% [sec_11]
