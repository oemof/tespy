from tespy.networks import Network
from tespy.components import (Source, Sink, DiabaticCombustionChamber, Compressor)
from tespy.connections import Connection, Bus, Ref
from tespy.tools import ExergyAnalysis

nw = Network(T_unit="C", p_unit="bar", h_unit="kJ / kg")

# components
comp = Compressor('Compressor')
src = Source('Source')
snk = Sink('Sink')

# connections
c1 = Connection(src, 'out1', comp, 'in1', label='Source to Compressor')
c2 = Connection(comp, 'out1', snk, 'in1', label='Compressor to Sink')

nw.add_conns(c1, c2)

# define parameters
comp.set_attr(pr=5, eta_s=0.8)

c1.set_attr(fluid={'air': 1}, p=1, T=20, v=50)

# solve
nw.solve(mode='design')
nw.print_results()


""" +++ exergy analysis +++ """
# define ambient
p_amb = 0.1
T_amb = 25

# define busses (no need to add them to system)
power = Bus('Power Input')
power.add_comps({'comp': comp, 'base': 'bus'})

air_in = Bus('Air In')
air_in.add_comps({'comp': src, 'base': 'bus'})

air_out = Bus('Air Out')
air_out.add_comps({'comp': snk})


# exergy analysis
exe_eco_input = {'Compressor_Z': 500, 'Source_c': 0.0001, 'Power Input_c': 0.002}
ean = ExergyAnalysis(nw, E_P=[air_out], E_F=[air_in, power], E_L=[])
ean.analyse(pamb=p_amb, Tamb=T_amb)
ean.evaluate_exergoeconomics(Exe_Eco_Costs=exe_eco_input, Tamb=T_amb)
ean.print_results(Exe_Eco_An=True)
