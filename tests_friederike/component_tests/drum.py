from tespy.networks import Network
from tespy.components import (Source, Sink, Drum)
from tespy.connections import Connection, Bus
from tespy.tools import ExergyAnalysis
from tespy.tools.helpers import get_chem_ex_lib
chemexlib = get_chem_ex_lib("Ahrendts")


# network
nw = Network(fluids=["Water"], T_unit="C", p_unit="bar", h_unit="kJ / kg")

# components
drum = Drum('Drum')
so1 = Source("source 1")
so2 = Source("source 2")
si_sat_l = Sink("sat liquid")
si_sat_g = Sink("sat gas")

# connections
liq_in_2_hx = Connection(so1, 'out1', drum, 'in1', label='inlet 1')
gas_in_2_hx = Connection(so2, 'out1', drum, 'in2', label='inlet 2')
sat_liq_out = Connection(drum, 'out1', si_sat_l, 'in1', label='saturated liquid')
sat_gas_out = Connection(drum, 'out2', si_sat_g, 'in1', label='saturated gas')

nw.add_conns(liq_in_2_hx, sat_liq_out)
nw.add_conns(gas_in_2_hx, sat_gas_out)


# define parameters
liq_in_2_hx.set_attr(T=60, p=1, fluid={'Water': 1}, m=5)
gas_in_2_hx.set_attr(T=200, fluid={'Water': 1}, m=50)


# solve
nw.solve(mode='design')
nw.print_results()


""" +++ exergy analysis +++ """
# define ambient
T_amb = 10
p_amp = 1

# define busses (no need to add them to system)
eco_out = Bus('liquid in')
eco_out.add_comps({'comp': so1, 'base': 'bus'})

liq_out = Bus('liquid out')
liq_out.add_comps({'comp': si_sat_l})

eva_out = Bus('gas in')
eva_out.add_comps({'comp': so2, 'base': 'bus'})

gas_out = Bus('gas out')
gas_out.add_comps({'comp': si_sat_g})



ean = ExergyAnalysis(nw, E_F=[eco_out, eva_out], E_P=[liq_out, gas_out], E_L=[])
exe_eco_input = {'Drum_Z': 50, 'source 1_c': 10/(10**9/3600), 'source 2_c': 10/(10**9/3600)}
ean.analyse(pamb=p_amp, Tamb=T_amb,  Chem_Ex= chemexlib)
ean.evaluate_exergoeconomics(Tamb=T_amb, Exe_Eco_Costs=exe_eco_input)
ean.print_results(Exe_Eco_An=True)
