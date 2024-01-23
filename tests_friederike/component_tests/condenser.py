from tespy.networks import Network
from tespy.components import (Source, Sink, Condenser)
from tespy.connections import Connection, Bus
from tespy.tools import ExergyAnalysis
from tespy.tools.helpers import get_chem_ex_lib
chemexlib = get_chem_ex_lib("Ahrendts")


nw = Network(T_unit="C", p_unit="bar", h_unit="kJ / kg")

# components
cond = Condenser('Heat Exchanger')
so_h = Source("Source Hot")
si_h = Sink("Sink Hot")
so_c = Source("Source Cold")
si_c = Sink("Sink Cold")

# connections
c1 = Connection(so_h, 'out1', cond, 'in1', label='Hot In')
c2 = Connection(so_c, 'out1', cond, 'in2', label='Cold In')
c3 = Connection(cond, 'out1', si_h, 'in1', label='Hot Out')
c4 = Connection(cond, 'out2', si_c, 'in1', label='Cold Out')

nw.add_conns(c1, c2, c3, c4)

# define parameters
cond.set_attr(pr1=1, pr2=1)
c1.set_attr(T=110, p=1.2, fluid={'Water': 1}, m=5)
c2.set_attr(T=20, p=1.2, fluid={'Water': 1})
c4.set_attr(T=30)

# solve
nw.solve(mode='design')
nw.print_results()


""" +++ exergy analysis +++ """
# define ambient
T_amb = 10
p_amp = 1

# define busses (no need to add them to system)
working_fluid_heat_out = Bus('Heat In')
working_fluid_heat_out.add_comps({'comp': so_h, 'base': 'bus'}, {'comp':si_h})

cooling_fluid_heat_in = Bus('Heated fluid')
cooling_fluid_heat_in.add_comps({'comp':so_c, 'base': 'bus'}, {'comp': si_c})

# exergy and exergoeconomic analysis
exe_eco_input = {'Heat Exchanger_Z': 5e3, 'Source Hot_c': 0.02, 'Source Cold_c': 0.05}
ean = ExergyAnalysis(nw, E_F=[working_fluid_heat_out], E_P=[cooling_fluid_heat_in], E_L=[])
ean.analyse(pamb=p_amp, Tamb=T_amb, Chem_Ex=chemexlib)
ean.evaluate_exergoeconomics(Tamb=T_amb, Exe_Eco_Costs=exe_eco_input)
ean.print_results(Exe_Eco_An=True)
