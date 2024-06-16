from tespy.networks import Network
from tespy.components import (Source, Sink, Splitter)
from tespy.connections import Connection, Bus
from tespy.tools import ExergyAnalysis
from tespy.tools.helpers import get_chem_ex_lib
chemexlib = get_chem_ex_lib("Ahrendts")

# out3 is not available. Choose from out1, out2.

# network
nw = Network(T_unit="C", p_unit="bar", h_unit="kJ / kg")

# components
sp = Splitter('Splitter', num_out = 2)
so = Source("source_stream")
si = Sink("out_1")
si_2 = Sink("out_2")
#si_3 = Sink("out_3")

# connections
inl = Connection(so, 'out1', sp, 'in1', label='stream in splitter')
out1 = Connection(sp, 'out1', si, 'in1', label='out 1')
out2 = Connection(sp, 'out2', si_2, 'in1', label='out 2')
#out3 = Connection(sp, 'out3', si_3, 'in1', label='out 3')

nw.add_conns(inl, out1, out2)       # out3


# define parameters
inl.set_attr(T=60, p=2, fluid={'Water': 1}, m=5)
out1.set_attr(m=1)
#out1.set_attr(m=2)

# solve
nw.solve(mode='design')
nw.print_results()


""" +++ exergy analysis +++ """
# define ambient
T_amb = 20
p_amp = 1

# define busses
so_in_B = Bus('in')
so_in_B.add_comps({'comp': so, 'base': 'bus'})

si_out = Bus('out')
si_out.add_comps({'comp': si})

si_out2 = Bus('out2')
si_out2.add_comps({'comp': si_2})

# exergy balance not closed if si und si_2 in einem bus zsm

# exergoeconomic analysis
ean = ExergyAnalysis(nw, E_F=[so_in_B], E_P=[si_out2, si_out], E_L=[])
exe_eco_input = {'Splitter_Z': 5e3, 'source_stream_c': 10}
ean.analyse(pamb=p_amp, Tamb=T_amb, Chem_Ex=chemexlib)
ean.evaluate_exergoeconomics(Tamb=T_amb, Exe_Eco_Costs=exe_eco_input)
ean.print_results(Exe_Eco_An=True)
