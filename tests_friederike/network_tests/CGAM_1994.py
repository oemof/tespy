from tespy.networks import Network
from tespy.components import (Sink, Source, Compressor, DiabaticCombustionChamber, HeatExchanger, Turbine, Drum, Pump)
from tespy.connections import Connection, Bus
from tespy.tools import ExergyAnalysis
from tespy.tools.helpers import get_chem_ex_lib
chemexlib = get_chem_ex_lib("Ahrendts")
from CoolProp.CoolProp import PropsSI as CP

# fluids
air_molar = {'O2': 0.2059, 'N2': 0.7748, 'CO2': 0.0003, 'H2O': 0.019}
molar_masses = {key: CP('M', key) * 1000 for key in air_molar}
M_air = sum([air_molar[key] * molar_masses[key] for key in air_molar])

air = {key: value / M_air * molar_masses[key] for key, value in air_molar.items()}
fuel = {'CH4': 1}

# network
cgam = Network(T_unit='C', p_unit='bar', h_unit='kJ / kg', m_unit='kg / s')

# components
src_air = Source('air-source')
src_fuel = Source('fuel-source')
snk_fluegas = Sink('flue-gas-sink')
cmp_AC = Compressor('air-compressor')
cmp_APH = HeatExchanger('air-pre-heater')
cmp_CC = DiabaticCombustionChamber('combustion-chamber')
cmp_EXP = Turbine('expander')

# connections
c1 = Connection(src_air, 'out1', cmp_AC, 'in1', label='1')
c2 = Connection(cmp_AC, 'out1', cmp_APH, 'in2', label='2')
c3 = Connection(cmp_APH, 'out2', cmp_CC, 'in1', label='3')
c4 = Connection(cmp_CC, 'out1', cmp_EXP, 'in1', label='4')
c5 = Connection(cmp_EXP, 'out1', cmp_APH, 'in1', label='5')
c6 = Connection(cmp_APH, 'out1', snk_fluegas, 'in1', label='6')
c10 = Connection(src_fuel, 'out1', cmp_CC, 'in2', label='10')

cgam.add_conns(c1,c2,c3,c4,c5,c6,c10)

# parameters
# connections
c1.set_attr(p=1.013, T=25, fluid=air, m=95.9185)
c2.set_attr(p=10.13, T=347.66)
c3.set_attr(p=9.62, T=576.85)
c4.set_attr(p=9.14, T=1246.85)
c5.set_attr(p=1.1, T=712.48)
c6.set_attr(p=1.07)
c10.set_attr(p=12, T=25, fluid=fuel, m=1.7653)

# solve network
cgam.solve('design')
cgam.print_results()


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

generator = Bus("product")
generator.add_comps({"comp": cmp_EXP, "base": "component",
                     "comp": cmp_AC, "base": "bus"})


# exergy and exergoeconomic analysis
exe_eco_input = {'air-compressor_Z': 5, 'air-pre-heater_Z': 2, 'combustion-chamber_Z': 2, 'expander_Z': 4,
                 'air-source_c': 0.02, 'fuel-source_c': 0.01}
ean = ExergyAnalysis(cgam, E_P=[generator], E_F=[fuel_bus], E_L=[loss_bus])
ean.analyse(pamb=p_amb, Tamb=T_amb, Chem_Ex=chemexlib)
ean.evaluate_exergoeconomics(Exe_Eco_Costs=exe_eco_input, Tamb=T_amb)
ean.print_results(Exe_Eco_An=True)
