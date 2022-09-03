# %%[sec_1]
from tespy.networks import Network

from tespy.components import (
    Condenser, Compressor, CycleCloser,  HeatExchanger,
    HeatExchangerSimple, Pump, Sink, Source, Valve
    )

from tespy.connections import Connection, Bus
# %%[sec_2]
wf = "NH3"

# network
nw = Network(
    fluids=["water", wf],
    T_unit="C", p_unit="bar", h_unit="kJ / kg", m_unit="kg / s"
    )

# components
cycle_closer = CycleCloser("Refrigerant Cycle Closer")

# heat source
heatsource_feedflow = Source("Heat Source Feed Flow")
heatsource_pump = Pump("Heat Source Recirculation Pump")
heatsource_evaporator = HeatExchanger("Heat Source Evaporator")
heatsource_backflow = Sink("Heat Source Back Flow")

# compression
compressor = Compressor("Compressor")

# heat sink
cons_pump = Pump("Heat Sink Recirculation Pump")
condenser = Condenser("Heat Sink Condenser")
cons_heatsink = HeatExchangerSimple("Heat Consumer")
cons_cycle_closer = CycleCloser("Consumer Feed Flow")

# internal heat exchange
int_heatex = HeatExchanger("Internal Heat Exchanger")

# expansion
valve = Valve("Expansion Valve")

# connections
# main cycle
c0 = Connection(cycle_closer, "out1", heatsource_evaporator, "in2", label="0")
c1 = Connection(heatsource_evaporator, "out2", int_heatex, "in2", label="1")
c2 = Connection(int_heatex, "out2", compressor, "in1", label="2")
c3 = Connection(compressor, "out1", condenser, "in1", label="3")
c4 = Connection(condenser, "out1", int_heatex, "in1", label="4")
c5 = Connection(int_heatex, "out1", valve, "in1", label="5")
c6 = Connection(valve, "out1", cycle_closer, "in1", label="6")

nw.add_conns(
    c0, c1, c2, c3, c4,
    c5, c6
    )

# heat source
c11 = Connection(heatsource_feedflow, "out1", heatsource_pump, "in1", label="11")
c12 = Connection(heatsource_pump, "out1", heatsource_evaporator, "in1", label="12")
c13 = Connection(heatsource_evaporator, "out1", heatsource_backflow, "in1", label="13")

nw.add_conns(c11, c12, c13)

# heat sink
c21 = Connection(cons_cycle_closer, "out1", cons_pump, "in1", label="21")
c22 = Connection(cons_pump, "out1", condenser, "in2", label="22")
c23 = Connection(condenser, "out2", cons_heatsink, "in1", label="23")
c24 = Connection(cons_heatsink, "out1", cons_cycle_closer, "in1", label="24")

nw.add_conns(c21, c22, c23, c24)
# %%[sec_3]
# parametrization connections
# set feedflow and backflow temperature of heat source and consumer
T_hs_bf = 10
T_hs_ff = 15
T_cons_bf = 50
T_cons_ff = 90

# consumer cycle
c23.set_attr(T=T_cons_ff, p=10, fluid={"water": 1, wf: 0})
c24.set_attr(T=T_cons_bf)

# heat source cycle
c11.set_attr(T=T_hs_ff, p=1, fluid={"water": 1, wf: 0})
c13.set_attr(T=T_hs_bf, p=1)

# evaporation to fully saturated gas
c1.set_attr(x=1, fluid={"water": 0, wf: 1})
# degree of overheating after internal heat exchanger (evaporation side)
c2.set_attr(Td_bp=10)

# parametrization components
# isentropic efficiency
cons_pump.set_attr(eta_s=0.8)
heatsource_pump.set_attr(eta_s=0.8)
compressor.set_attr(eta_s=0.85)

# pressure ratios
condenser.set_attr(pr1=0.98, pr2=0.98)
heatsource_evaporator.set_attr(pr1=0.98, pr2=0.98)
cons_heatsink.set_attr(pr=0.99)
int_heatex.set_attr(pr1=0.98, pr2=0.98)

# temperature differences
heatsource_evaporator.set_attr(ttd_l=5)
condenser.set_attr(ttd_u=5)

# consumer heat demand
cons_heatsink.set_attr(Q=-1e6)

try:
    nw.solve("design")
except ValueError as e:
    print(e)
# %%[sec_4]
import CoolProp.CoolProp as CP

# evaporation point
p_eva = CP.PropsSI("P", "Q", 1, "T", T_hs_bf - 5 + 273.15, wf) * 1e-5
c1.set_attr(p=p_eva)
heatsource_evaporator.set_attr(ttd_l=None)

# condensation point
p_cond = CP.PropsSI("P", "Q", 0, "T", T_cons_ff + 5 + 273.15, wf) * 1e-5
c4.set_attr(p=p_cond)
condenser.set_attr(ttd_u=None)

# internal heat exchanger to compressor enthalpy
h_evap = CP.PropsSI("H", "Q", 1, "T", T_hs_bf - 5 + 273.15, wf) * 1e-3
c2.set_attr(Td_bp=None, h=h_evap * 1.01)

# solve the network again
nw.solve("design")
# %%[sec_5]
# evaporation point
c1.set_attr(p=None)
heatsource_evaporator.set_attr(ttd_l=5)

# condensation point
c4.set_attr(p=None)
condenser.set_attr(ttd_u=5)

# internal heat exchanger superheating
c2.set_attr(Td_bp=5, h=None)

# solve the network again
nw.solve("design")

# calculate the COP
cop = abs(
    cons_heatsink.Q.val
    / (cons_pump.P.val + heatsource_pump.P.val + compressor.P.val)
)
print(cop)
# %%[sec_6]
def generate_network_with_starting_values(wf):
    # network
    nw = Network(
        fluids=["water", wf],
        T_unit="C", p_unit="bar", h_unit="kJ / kg", m_unit="kg / s",
        iterinfo=False
    )

    # components
    cycle_closer = CycleCloser("Refrigerant Cycle Closer")

    # heat source
    heatsource_feedflow = Source("Heat Source Feed Flow")
    heatsource_pump = Pump("Heat Source Recirculation Pump")
    heatsource_evaporator = HeatExchanger("Heat Source Evaporator")
    heatsource_backflow = Sink("Heat Source Back Flow")

    # compression
    compressor = Compressor("Compressor")

    # heat sink
    cons_pump = Pump("Heat Sink Recirculation Pump")
    condenser = Condenser("Heat Sink Condenser")
    cons_heatsink = HeatExchangerSimple("Heat Consumer")
    cons_cycle_closer = CycleCloser("Consumer Feed Flow")

    # internal heat exchange
    int_heatex = HeatExchanger("Internal Heat Exchanger")

    # expansion
    valve = Valve("Expansion Valve")

    # connections
    # main cycle
    c0 = Connection(cycle_closer, "out1", heatsource_evaporator, "in2", label="0")
    c1 = Connection(heatsource_evaporator, "out2", int_heatex, "in2", label="1")
    c2 = Connection(int_heatex, "out2", compressor, "in1", label="2")
    c3 = Connection(compressor, "out1", condenser, "in1", label="3")
    c4 = Connection(condenser, "out1", int_heatex, "in1", label="4")
    c5 = Connection(int_heatex, "out1", valve, "in1", label="5")
    c6 = Connection(valve, "out1", cycle_closer, "in1", label="6")

    nw.add_conns(
        c0, c1, c2, c3, c4,
        c5, c6
        )

    # heat source
    c11 = Connection(heatsource_feedflow, "out1", heatsource_pump, "in1", label="11")
    c12 = Connection(heatsource_pump, "out1", heatsource_evaporator, "in1", label="12")
    c13 = Connection(heatsource_evaporator, "out1", heatsource_backflow, "in1", label="13")

    nw.add_conns(c11, c12, c13)

    # heat sink
    c21 = Connection(cons_cycle_closer, "out1", cons_pump, "in1", label="20")
    c22 = Connection(cons_pump, "out1", condenser, "in2", label="21")
    c23 = Connection(condenser, "out2", cons_heatsink, "in1", label="22")
    c24 = Connection(cons_heatsink, "out1", cons_cycle_closer, "in1", label="23")

    nw.add_conns(c21, c22, c23, c24)

    # set feedflow and backflow temperature of heat source and consumer
    T_hs_bf = 10
    T_hs_ff = 15
    T_cons_bf = 50
    T_cons_ff = 90

    # consumer cycle
    c23.set_attr(T=T_cons_ff, p=10, fluid={"water": 1, wf: 0})
    c24.set_attr(T=T_cons_bf)

    # heat source cycle
    c11.set_attr(T=T_hs_ff, p=1, fluid={"water": 1, wf: 0})
    c13.set_attr(T=T_hs_bf, p=1)

    # evaporation to fully saturated gas
    c1.set_attr(x=1, fluid={"water": 0, wf: 1})

    # parametrization components
    # isentropic efficiency
    cons_pump.set_attr(eta_s=0.8)
    heatsource_pump.set_attr(eta_s=0.8)
    compressor.set_attr(eta_s=0.85)

    # pressure ratios
    condenser.set_attr(pr1=0.98, pr2=0.98)
    heatsource_evaporator.set_attr(pr1=0.98, pr2=0.98)
    cons_heatsink.set_attr(pr=0.99)
    int_heatex.set_attr(pr1=0.98, pr2=0.98)

    # evaporation point
    p_eva = CP.PropsSI("P", "Q", 1, "T", T_hs_bf - 5 + 273.15, wf) * 1e-5
    c1.set_attr(p=p_eva)

    # condensation point
    p_cond = CP.PropsSI("P", "Q", 0, "T", T_cons_ff + 5 + 273.15, wf) * 1e-5
    c4.set_attr(p=p_cond)

    # internal heat exchanger to compressor enthalpy
    h_evap = CP.PropsSI("H", "Q", 1, "T", T_hs_bf - 5 + 273.15, wf) * 1e-3
    c2.set_attr(h=h_evap * 1.01)

    # consumer heat demand
    cons_heatsink.set_attr(Q=-1e6)

    power_bus = Bus("Total power input")
    heat_bus = Bus("Total heat production")
    power_bus.add_comps(
        {"comp": compressor, "base": "bus"},
        {"comp": cons_pump, "base": "bus"},
        {"comp": heatsource_pump, "base": "bus"},
    )
    heat_bus.add_comps({"comp": cons_heatsink})

    nw.add_busses(power_bus, heat_bus)

    nw.solve("design")

        # evaporation point
    c1.set_attr(p=None)
    heatsource_evaporator.set_attr(ttd_l=5)

    # condensation point
    c4.set_attr(p=None)
    condenser.set_attr(ttd_u=5)

    # internal heat exchanger superheating
    c2.set_attr(Td_bp=5, h=None)

    # solve the network again
    nw.solve("design")

    return nw
# %%[sec_7]
import matplotlib.pyplot as plt
import pandas as pd


# make text reasonably sized
plt.rc("font", **{"size": 18})

cop = pd.DataFrame(columns=["COP"])

for wf in ["NH3", "R22", "R134a", "R152a", "R290", "R718"]:
    nw = generate_network_with_starting_values(wf)

    power = nw.busses["Total power input"].P.val
    heat = abs(nw.busses["Total heat production"].P.val)
    cop.loc[wf] = heat / power

fig, ax = plt.subplots(1, figsize=(16, 8))

cop.plot.bar(ax=ax, legend=False)

ax.set_axisbelow(True)
ax.yaxis.grid(linestyle="dashed")
ax.set_xlabel("Name of working fluid")
ax.set_ylabel("Coefficicent of performance")
plt.tight_layout()

fig.savefig("COP_by_wf.svg")
# %%[sec_8]
