# -*- coding: utf-8 -*-
"""
@author: Karim Shawky
"""
from CoolProp.CoolProp import PropsSI as CPSI
from tespy.networks import Network
from tespy.connections import Connection, Bus
from tespy.components import (
    Source, Sink, HeatExchanger,
    DiabaticCombustionChamber, Turbine, Compressor, Drum)
from tespy.tools.helpers import get_chem_ex_lib
from tespy.tools.analyses import ExergyAnalysis
import plotly.graph_objects as go


Chem_Ex = get_chem_ex_lib("Ahrendts")

fluid_list = ['O2', 'H2O', 'N2', 'CO2', 'CH4']
nwk = Network(fluids=fluid_list, p_unit='bar', T_unit='K')

air_molar = {
    'O2': 0.2059, 'N2': 0.7748, 'CO2': 0.0003, 'H2O': 0.019, 'CH4': 0
}
molar_masses = {key: CPSI('M', key) * 1000 for key in air_molar}
M_air = sum([air_molar[key] * molar_masses[key] for key in air_molar])

air = {key: value / M_air * molar_masses[key] for key, value in air_molar.items()}

water = {f: (0 if f != 'H2O' else 1) for f in air}
fuel = {f: (0 if f != 'CH4' else 1) for f in air}

amb = Source('ambient air')
ch4 = Source('methane')
h2o= Source('water')

flueout= Sink('Fluegas out')
steam= Sink('steam')


aph=HeatExchanger('Air Prehater')
cc=DiabaticCombustionChamber('Combustion Chamber')
exp= Turbine('Expander')
ac=Compressor('Compressor')

eva = HeatExchanger('evaporator', fkt_group='HRSG')
eco = HeatExchanger('economizer',  fkt_group='HRSG')
dr = Drum('drum',  fkt_group='HRSG')


amb_in_ac=Connection(amb, 'out1', ac, 'in1')
ac_in_aph=Connection(ac, 'out1', aph, 'in2')
aph_in_cc=Connection(aph, 'out2', cc, 'in1')

ch4_in_cc=Connection(ch4, 'out1', cc, 'in2')

cc_in_exp=Connection(cc, 'out1', exp, 'in1')
exp_in_aph=Connection(exp, 'out1', aph, 'in1')

#aph_in_out=Connection(aph, 'out1', flueout, 'in1')

aph_in_eva = Connection(aph, 'out1', eva, 'in1')
eva_in_eco = Connection(eva, 'out1', eco, 'in1')
eco_in_out = Connection(eco, 'out1', flueout, 'in1')


water_in_eco = Connection(h2o, 'out1', eco, 'in2')
eco_in_dr = Connection(eco, 'out2', dr, 'in1')
dr_in_eva = Connection(dr, 'out1', eva, 'in2')
eva_in_dr = Connection(eva, 'out2', dr, 'in2')
dr_in_steam = Connection(dr, 'out2', steam, 'in1')



nwk.add_conns(amb_in_ac,ac_in_aph,
              aph_in_cc, ch4_in_cc, cc_in_exp, exp_in_aph,
              aph_in_eva, eva_in_eco, eco_in_out, water_in_eco, eco_in_dr,
              dr_in_eva, eva_in_dr, dr_in_steam)

power = Bus('total power')
power.add_comps({'comp': ac, 'base': 'bus'}, {'comp': exp})

E_F= Bus('exergy fuel')
E_F.add_comps({'comp':ch4, 'base':'bus'},{'comp':amb, 'base':'bus'})
Heat_water=Bus('Heat product')
Heat_water.add_comps({'comp':h2o, 'base':'bus'},{'comp':steam, 'base':'component'})

E_Loss=Bus('exergy loss')
E_Loss.add_comps({'comp':flueout})

nwk.add_busses(E_F, E_Loss, power, Heat_water)


amb_in_ac.set_attr(fluid=air, T=298.15, p=1.013, m=91.28)

aph_in_cc.set_attr(T=850)
ch4_in_cc.set_attr(fluid=fuel, T=298.15, p=12)

cc_in_exp.set_attr(T=1520)

water_in_eco.set_attr(fluid=water ,p=20, T=298.15, m=14)
eco_in_dr.set_attr(Td_bp=-15)
eva_in_dr.set_attr(x=0.5)


cc.set_attr(eta=0.98, pr=0.95)
exp.set_attr(eta_s=0.86)
aph.set_attr(pr1=0.97, pr2=0.95)
ac.set_attr(eta_s=0.86, pr=10)

eva.set_attr(pr1=0.95 ** 0.5)
eco.set_attr(pr1=0.95 ** 0.5, pr2=1)

power.set_attr(P=-30e6)
nwk.solve('design')
nwk.print_results()



ean=ExergyAnalysis(nwk, E_P=[power, Heat_water], E_F=[E_F], E_L=[E_Loss])
ean.analyse(pamb=1.013, Tamb=298.15)


ean.print_results()

network_result = ean.network_data.to_frame().transpose()



for c in nwk.conns["object"]:
    c.get_physical_exergy(1.013e5, 298.15)
    c.get_chemical_exergy(1.013e5, 298.15, Chem_Ex)

##Following Tstsaronis definition
##Exergy Balance Combustion Chamber without Split

# E_P=cc_in_exp.Ex_physical+cc_in_exp.Ex_chemical-(amb_in_cc.Ex_physical+amb_in_cc.Ex_chemical)
# E_F= (ch4_in_cc.Ex_chemical+ch4_in_cc.Ex_physical)

# 70.17%

# # Exergy Balance Combustion Chamber Chemical/Physical Split
# E_P= cc_in_exp.Ex_physical-aph_in_cc.Ex_physical-ch4_in_cc.Ex_physical
# E_F= ch4_in_cc.Ex_chemical-cc_in_exp.Ex_chemical+aph_in_cc.Ex_chemical

#print(E_P) # 59MW
#print(E_F-E_P) #25.5MW
# 69.81%


#Bejan
#E_P=cc_in_exp.Ex_physical+cc_in_exp.Ex_chemical
#E_F= aph_in_cc.Ex_chemical+ aph_in_cc.Ex_physical + ch4_in_cc.Ex_chemical+ ch4_in_cc.Ex_physical


# Always when there is a combustion chamber it should calculate the chemical exergy
#79.9%


# E_P= (aph_in_cc.m.val* (cc_in_exp.ex_physical+cc_in_exp.ex_chemical)-
#       (aph_in_cc.Ex_chemical+ aph_in_cc.Ex_physical) +
#       ch4_in_cc.m.val*(-ch4_in_cc.ex_physical+cc_in_exp.ex_physical))
# E_F= ch4_in_cc.m.val *(ch4_in_cc.ex_chemical-cc_in_exp.ex_chemical )
#
# 69.96 %
E_P=dr_in_steam.Ex_chemical+dr_in_steam.Ex_physical-(water_in_eco.Ex_chemical
                                                     +water_in_eco.Ex_physical)
E_F=aph_in_eva.Ex_physical-eco_in_out.Ex_physical


epsilon=E_P/E_F
print(epsilon)
links, nodes = ean.generate_plotly_sankey_input()
links['value'] = [val / links['value'][0] for val in links['value']]

fig = go.Figure(data=[go.Sankey(
    arrangement="snap",
    textfont={"family": "Linux Libertine O"},
    node={
        "label": nodes,
        'pad':11,
        'color': 'orange'},
    link=links)])
fig.show()

#print(ebsilon)
#plot(fig, filename='NH3_sankey')