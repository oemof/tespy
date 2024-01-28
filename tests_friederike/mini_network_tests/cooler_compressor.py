from tespy.networks import Network
from tespy.components import (Source, Sink, HeatExchanger, Compressor)
from tespy.connections import Connection, Bus
from tespy.tools import ExergyAnalysis
from tespy.tools.helpers import get_chem_ex_lib
chemexlib = get_chem_ex_lib("Ahrendts")


nw = Network(T_unit="C", p_unit="bar", h_unit="kJ / kg")

# components
wf_in = Source("Working Fluid Inlet")
hx = HeatExchanger('Cooler')
cp = Compressor("Compressor")
wf_out = Sink("Working Fluid Outlet")
cw_in = Source("Cooling Water Inlet")
cw_out = Sink("Cooling Water Outlet")


# connections
c1 = Connection(wf_in, 'out1', hx, 'in1', label='Inlet')
c2 = Connection(hx, 'out1', cp, 'in1', label='Cooler to Compressor')
c3 = Connection(cp, 'out1', wf_out, 'in1', label='Outlet')
c4 = Connection(cw_in, 'out1', hx, 'in2', label='Cooling Water Inlet')
c5 = Connection(hx, 'out2', cw_out, 'in1', label='Cooling Water Outlet')

nw.add_conns(c1, c2, c3, c4, c5)


# define parameters
hx.set_attr(pr1=1, pr2=1)
cp.set_attr(pr=5, eta_s=0.8)
c1.set_attr(T=25, p=1.2, fluid={'N2': 1}, m=5)
c2.set_attr(T=20)
c4.set_attr(T=10, p=1.2, fluid={'H2O': 1}, m=2)


# solve
nw.solve(mode='design')
nw.print_results()


""" +++ exergy analysis +++ """
# define ambient
T_amb = 5
p_amp = 1

# define busses (no need to add them to system)
working_fluid = Bus('Working Fluid')
working_fluid.add_comps({'comp': wf_in, 'base': 'bus'}, {'comp':wf_out, 'base': 'component'})

cooling_water = Bus('Cooling Water')
cooling_water.add_comps({'comp':cw_in, 'base': 'bus'}, {'comp': cw_out, 'base': 'component'})

power = Bus('Power Input')
power.add_comps({'comp': cp, 'base': 'bus'})


# exergy and exergoeconomic analysis
hx.dissipative = True
hx.serving_component = cp
exe_eco_input = {'Cooler_Z': 5e1, 'Compressor_Z': 8e1, 'Working Fluid Inlet_c': 0.05, 'Cooling Water Inlet_c': 0.02, 'Power Input_c': 0.2}
ean = ExergyAnalysis(nw, E_F=[cooling_water, power], E_P=[working_fluid], E_L=[])
ean.analyse(pamb=p_amp, Tamb=T_amb, Chem_Ex=chemexlib)
ean.evaluate_exergoeconomics(Tamb=T_amb, Exe_Eco_Costs=exe_eco_input)
ean.print_results(Exe_Eco_An=True)
