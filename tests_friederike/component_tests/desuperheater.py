from tespy.networks import Network
from tespy.components import (Source, Sink, Desuperheater)
from tespy.connections import Connection, Bus
from tespy.tools import ExergyAnalysis
from tespy.tools.helpers import get_chem_ex_lib
chemexlib = get_chem_ex_lib("Ahrendts")


nw = Network(T_unit="C", p_unit="bar", h_unit="kJ / kg")

# components
desu = Desuperheater('Desuperheater')      # cools fluid to saturated gas state
et_in = Source("Ethanol Inlet")
et_out = Sink("Ethanol Outlet")
cw_in = Source("Cooling Water Inlet")
cw_out = Sink("Cooling Water Outlet")

# connections
c1 = Connection(et_in, 'out1', desu, 'in1', label='Ethanol In')
c2 = Connection(cw_in, 'out1', desu, 'in2', label='Water In')
c3 = Connection(desu, 'out1', et_out, 'in1', label='Ethanol Out')
c4 = Connection(desu, 'out2', cw_out, 'in1', label='Water Out')

nw.add_conns(c1, c2, c3, c4)

# define parameters
desu.set_attr(pr1=1, pr2=1)
c1.set_attr(Td_bp=100, v=10, fluid={'ethanol': 1})
c3.set_attr(p=1)
c2.set_attr(T=15, v=1, fluid={'water': 1})
c4.set_attr(p=1)

# solve
nw.solve(mode='design')
nw.print_results()


""" +++ exergy analysis +++ """
# define ambient
T_amb = 20
p_amp = 1

# define busses (no need to add them to system)
ethanol = Bus('Ethanol')
ethanol.add_comps({'comp': et_in, 'base': 'bus'}, {'comp':et_out})

cooling_water = Bus('Cooling Water')
cooling_water.add_comps({'comp':cw_in, 'base': 'bus'}, {'comp': cw_out})

# exergy and exergoeconomic analysis
exe_eco_input = {'Desuperheater_Z': 5e3, 'Ethanol Inlet_c': 0.02, 'Cooling Water Inlet_c': 0.05}
ean = ExergyAnalysis(nw, E_F=[cooling_water], E_P=[ethanol], E_L=[])
ean.analyse(pamb=p_amp, Tamb=T_amb, Chem_Ex=chemexlib)
ean.evaluate_exergoeconomics(Tamb=T_amb, Exe_Eco_Costs=exe_eco_input)
ean.print_results(Exe_Eco_An=True)
