# -*- coding: utf-8 -*-
"""
@author: Karim Shawky
"""
from CoolProp.CoolProp import PropsSI as CPSI
from tespy.networks import Network
from tespy.connections import Connection, Bus
from tespy.components import (
    Source, Sink, HeatExchanger,
    DiabaticCombustionChamber, Turbine, Compressor)
from tespy.tools.Chem_Ex_libs.libChemExAhrendts import Chem_Ex
from tespy.tools.analyses import ExergyAnalysis

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


flueout= Sink('Fluegas out')

#aph=HeatExchanger('Air Prehater')
cc=DiabaticCombustionChamber('Combustion Chamber')
#exp= Turbine('Expander')
#ac=Compressor('Compressor')


# amb_in_ac=Connection(amb, 'out1', ac, 'in1')
# ac_in_aph=Connection(ac, 'out1', aph, 'in2')

aph_in_cc=Connection(amb, 'out1', cc, 'in1')
ch4_in_cc=Connection(ch4, 'out1', cc, 'in2')

cc_in_exp=Connection(cc, 'out1', flueout, 'in1')

#exp_in_aph=Connection(exp, 'out1', aph, 'in1')

#aph_in_out=Connection(aph, 'out1', flueout, 'in1')
nwk.add_conns(#amb_in_ac,ac_in_aph,
              aph_in_cc, ch4_in_cc, cc_in_exp)#, exp_in_aph, aph_in_out)

# power = Bus('total power')
# power.add_comps({'comp': ac, 'base': 'bus'}, {'comp': exp})

E_F= Bus('exergy fuel')
E_F.add_comps({'comp':ch4, 'base':'bus'})
E_P=Bus('exergy product')
E_P.add_comps({'comp':flueout},{'comp':amb, 'base':'bus'})

nwk.add_busses(E_F, E_P)
#amb_in_ac.set_attr(fluid=air, T=298.15, p=1.013)

aph_in_cc.set_attr( fluid=air , T=850, p=9.623, m=91.28)
ch4_in_cc.set_attr(fluid=fuel, T=298.15, p=12)
#aph_in_out.set_attr(T=779.75)
cc_in_exp.set_attr(T=1520)

cc.set_attr(eta=0.98, pr=0.95)
# exp.set_attr(eta_s=0.86, P=-30e6)
# aph.set_attr(pr1=0.97, pr2=0.95)
# ac.set_attr(eta_s=0.86)

nwk.solve('design')
nwk.print_results()



ean=ExergyAnalysis(nwk, E_P=[E_P], E_F=[E_F])
ean.analyse(pamb=1.013e5, Tamb=298.15)


ean.print_results()



for c in nwk.conns["object"]:
    c.get_physical_exergy(1.013e5, 298.15)
    c.get_chemical_exergy(1.013e5, 298.15, Chem_Ex)
    
##Following Tstsaronis definition
##Exergy Balance Combustion Chamber without Split

# E_P=cc_in_exp.Ex_physical+cc_in_exp.Ex_chemical-(amb_in_cc.Ex_physical+amb_in_cc.Ex_chemical)
# E_F= (ch4_in_cc.Ex_chemical+ch4_in_cc.Ex_physical) 

# 70.17%

# # Exergy Balance Combustion Chamber Chemical/Physical Split
E_P= cc_in_exp.Ex_physical-aph_in_cc.Ex_physical-ch4_in_cc.Ex_physical
E_F= ch4_in_cc.Ex_chemical-cc_in_exp.Ex_chemical+aph_in_cc.Ex_chemical

print(E_P) # 59MW
print(E_F-E_P) #25.5MW
# 69.81%


#Bejan
# E_P=cc_in_exp.Ex_physical+cc_in_exp.Ex_chemical
# E_F= aph_in_cc.Ex_chemical+ aph_in_cc.Ex_physical + ch4_in_cc.Ex_chemical+ ch4_in_cc.Ex_physical

#  
 
# Always when there is a combustion chamber it should calculate the chemical exergy
#79.9%

print(ch4_in_cc.Ex_physical)    
ebsilon= E_P/E_F
#print(ebsilon)  
