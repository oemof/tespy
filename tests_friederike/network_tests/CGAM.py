from CoolProp.CoolProp import PropsSI as CPSI
from tespy.networks import Network
from tespy.components import (
    HeatExchanger, Turbine, Compressor, Drum,
    DiabaticCombustionChamber, Sink, Source
)
from tespy.connections import Connection, Bus
from tespy.tools.helpers import get_chem_ex_lib
chemexlib = get_chem_ex_lib("Ahrendts")

nw = Network(p_unit='bar', T_unit='C')

air_molar = {'O2': 0.2059, 'N2': 0.7748, 'CO2': 0.0003, 'H2O': 0.019, 'CH4': 0}
molar_masses = {key: CPSI('M', key) * 1000 for key in air_molar}
M_air = sum([air_molar[key] * molar_masses[key] for key in air_molar])

air = {key: value / M_air * molar_masses[key] for key, value in air_molar.items()}

amb = Source('ambient air')
ch4 = Source('methane')
fw = Source('feed water')

ch = Sink('chimney')
ls = Sink('live steam')

cmp = Compressor('compressor')
aph = HeatExchanger('air preheater')
cb = DiabaticCombustionChamber('combustion chamber')
tur = Turbine('gas turbine')

eva = HeatExchanger('evaporator')
eco = HeatExchanger('economizer')
dr = Drum('drum')

c1 = Connection(amb, 'out1', cmp, 'in1', label='1')
c2 = Connection(cmp, 'out1', aph, 'in2', label='2')
c3 = Connection(aph, 'out2', cb, 'in1', label='3')
c10 = Connection(ch4, 'out1', cb, 'in2', label='10')

nw.add_conns(c1, c2, c3, c10)

c4 = Connection(cb, 'out1', tur, 'in1', label='4')
c5 = Connection(tur, 'out1', aph, 'in1', label='5')
c6 = Connection(aph, 'out1', eva, 'in1', label='6')
c6p = Connection(eva, 'out1', eco, 'in1', label='6p')
c7 = Connection(eco, 'out1', ch, 'in1', label='7')

nw.add_conns(c4, c5, c6, c6p, c7)

c8 = Connection(fw, 'out1', eco, 'in2', label='8')
c8p = Connection(eco, 'out2', dr, 'in1', label='8p')
c11 = Connection(dr, 'out1', eva, 'in2', label='11')
c11p = Connection(eva, 'out2', dr, 'in2', label='11p')
c9 = Connection(dr, 'out2', ls, 'in1', label='9')

nw.add_conns(c8, c8p, c11, c11p, c9)

c8.set_attr(p=20, T=25, m=14, fluid={'water': 1})
c1.set_attr(p=1.013, T=25, fluid=air, m=91.753028)
c10.set_attr(T=25, fluid={'CH4': 1}, p=12)
c7.set_attr(p=1.013)
c3.set_attr(T=850 - 273.15)
c4.set_attr(T=1520 - 273.15)
c8p.set_attr(Td_bp=-15)
c11p.set_attr(x=0.5)

cmp.set_attr(pr=10, eta_s=0.86)
cb.set_attr(eta=0.98, pr=0.95)
tur.set_attr(eta_s=0.86)
aph.set_attr(pr1=0.97, pr2=0.95)
eva.set_attr(pr1=0.95 ** 0.5)
eco.set_attr(pr1=0.95 ** 0.5, pr2=1)

power = Bus('total power')
power.add_comps({'comp': cmp, 'base': 'bus'}, {'comp': tur})

nw.add_busses(power)

heat_output = Bus('heat output')
power_output = Bus('power output')
fuel_input = Bus('fuel input')

heat_output.add_comps(
    {'comp': eco, 'char': -1},
    {'comp': eva, 'char': -1})
power_output.add_comps(
    {'comp': cmp, 'base': 'bus', 'char': 1},
    {'comp': tur, 'char': 1})
fuel_input.add_comps({'comp': cb, 'base': 'bus'})
nw.add_busses(heat_output, power_output, fuel_input)

nw.solve('design')
nw.print_results()

'''
Singularity in jacobian matrix, calculation aborted! 
Make sure your network does not have any linear dependencies in the parametrisation
'''
