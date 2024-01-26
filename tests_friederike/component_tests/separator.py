from tespy.networks import Network
from tespy.components import Sink, Source, Separator
from tespy.connections import Connection, Bus
from tespy.tools import ExergyAnalysis
from tespy.tools.helpers import get_chem_ex_lib
chemexlib = get_chem_ex_lib("Ahrendts")

# network
nw = Network(p_unit='bar', T_unit='C',iterinfo=False)

# components
so = Source('Source')
si1 = Sink('Sink1')
si2 = Sink('Sink2')
sep = Separator('Separator', num_out=2)

# connections
inc = Connection(so, 'out1', sep, 'in1')
outg1 = Connection(sep, 'out1', si1, 'in1')
outg2 = Connection(sep, 'out2', si2, 'in1')
nw.add_conns(inc, outg1, outg2)

# define parameters
inc.set_attr(fluid={'O2': 0.23, 'N2': 0.77}, p=1, T=20, m=5)
outg1.set_attr(fluid={'O2': 0.1, 'N2': 0.9}, m=1)
outg2.set_attr(fluid0={'O2': 0.5, 'N2': 0.5})

# solve
nw.solve('design')
nw.print_results()


""" +++ exergy analysis +++ """
# define ambient
T_amb = 20
p_amp = 1

# define busses
so_in_B = Bus('in')
so_in_B.add_comps({'comp': so, 'base': 'bus'})

si_out_B = Bus('out')
si_out_B.add_comps({'comp': si1})

si2_out_B = Bus('out2')
si2_out_B.add_comps({'comp': si2})



# exergy analysis
ean = ExergyAnalysis(nw, E_F=[so_in_B], E_P=[si2_out_B, si_out_B], E_L=[])
ean.analyse(pamb=p_amp, Tamb=T_amb, Chem_Ex=chemexlib)
ean.print_results()

"""
# exergoeconomic analysis
ean = ExergyAnalysis(nw, E_F=[so_in_B], E_P=[si2_out_B, si_out_B], E_L=[])
exe_eco_input = {'Separator_Z': 5e3, 'Source_c': 10}
ean.analyse(pamb=p_amp, Tamb=T_amb, Chem_Ex=chemexlib)
ean.evaluate_exergoeconomics(Tamb=T_amb, Exe_Eco_Costs=exe_eco_input)
ean.print_results(Exe_Eco_An=True)
"""
