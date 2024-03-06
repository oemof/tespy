from tespy.networks import Network
from tespy.components import (Sink, Source, CombustionChamber, HeatExchanger, Drum, Pump)
from tespy.connections import Connection, Bus
from tespy.tools import ExergyAnalysis
from tespy.tools.helpers import get_chem_ex_lib
chemexlib = get_chem_ex_lib("Ahrendts")
from CoolProp.CoolProp import PropsSI as CP

# fluids
air = {'O2': 0.21, 'N2': 0.79}
fuel = {'CH4': 1}
water = {'H2O': 1}

# network
steamboiler = Network(T_unit='C', p_unit='bar', h_unit='kJ / kg', m_unit='kg / s')

# components
src_air = Source('air-source')
src_fuel = Source('fuel-source')
snk_fluegas = Sink('flue-gas-sink')
src_water = Source('water')
snk_steam = Sink('steam')
cmp_cc = CombustionChamber('combustion chamber')
cmp_pump = Pump('pump')
cmp_econ = HeatExchanger('economizer')
cmp_drum = Drum('drum')
cmp_evap = HeatExchanger('evaporator')

# connections
c11 = Connection(src_air, 'out1', cmp_cc, 'in1', label='11 (air)')
c12 = Connection(src_fuel, 'out1', cmp_cc, 'in2', label='12 (fuel)')
c13 = Connection(cmp_cc, 'out1', cmp_evap, 'in1', label='13 (fluegas)')
c14 = Connection(cmp_evap, 'out1', cmp_econ, 'in1', label='14 (fluegas)')
c15 = Connection(cmp_econ, 'out1', snk_fluegas, 'in1', label='15 (fluegas)')
c21 = Connection(src_water, 'out1', cmp_pump, 'in1', label='21 (water, lq.)')
c22 = Connection(cmp_pump, 'out1', cmp_econ, 'in2', label='22 (water, lq.)')
c23 = Connection(cmp_econ, 'out2', cmp_drum, 'in1', label='23 (water, x=0)')
c24 = Connection(cmp_drum, 'out1', cmp_evap, 'in2', label='24 (drum circulation)')
c25 = Connection(cmp_evap, 'out2', cmp_drum, 'in2', label='25 (drum circulation)')
c26 = Connection(cmp_drum, 'out2', snk_steam, 'in1', label='26 (steam, x=1)')

steamboiler.add_conns(c11,c12,c13,c14,c15,c21,c22,c23,c24,c25,c26)

# parameters

# components
cmp_cc.set_attr(lamb=1.1)
cmp_pump.set_attr(eta_s=0.8)
cmp_econ.set_attr(pr1=1, pr2=1)
cmp_evap.set_attr(pr1=1)

# connections
c11.set_attr(p=1, T=25, fluid=air)
c12.set_attr(m=1.5, T=25, fluid=fuel)
c21.set_attr(m=25, p=1, T=60, fluid=water)
c22.set_attr(p=40)
c23.set_attr(h=CP('H', 'P', 40*1E5, 'Q', 0, 'water')/1E3)
c25.set_attr(x=0.6)

# solve
steamboiler.solve(mode='design')
steamboiler.print_results()

c12.set_attr(m=None)
cmp_econ.set_attr(ttd_l=50)

steamboiler.solve(mode='design')
steamboiler.print_results()


""" +++ exergy analysis +++ """
# define ambient
T_amb = 25
p_amb = 1.02

# define busses
fuel_bus = Bus("fuel")
fuel_bus.add_comps({"comp": src_fuel, "base": "bus",
                    "comp": src_air, "base": "bus"})

loss_bus = Bus("loss")
loss_bus.add_comps({"comp": snk_fluegas, "base": "component"})

product_bus = Bus("product")
product_bus.add_comps({"comp": src_water, "base": "bus",
                     "comp": snk_steam, "base": "component"})


# exergy and exergoeconomic analysis
exe_eco_input = {'pump_Z': 5, 'economizer_Z': 2, 'combustion-chamber_Z': 2, 'drum_Z': 4, 'evaporator_Z': 3,
                 'air-source_c': 0.02, 'fuel-source_c': 0.01, 'water_c': 0.01}
ean = ExergyAnalysis(steamboiler, E_P=[product_bus], E_F=[fuel_bus], E_L=[loss_bus])
ean.analyse(pamb=p_amb, Tamb=T_amb, Chem_Ex=chemexlib)
ean.evaluate_exergoeconomics(Exe_Eco_Costs=exe_eco_input, Tamb=T_amb)
ean.print_results(Exe_Eco_An=True)
