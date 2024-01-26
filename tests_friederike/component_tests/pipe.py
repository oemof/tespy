from tespy.networks import Network
from tespy.components import Source, Sink, Pipe
from tespy.connections import Connection, Bus
from tespy.tools import ExergyAnalysis
from tespy.tools.helpers import get_chem_ex_lib
chemexlib = get_chem_ex_lib("Ahrendts")

# network
nw = Network(p_unit='bar', T_unit='C', iterinfo=False)

# components
pi = Pipe('Pipe')
so = Source('Inlet')
si = Sink('Outlet')

# connections
c1 = Connection(so, 'out1', pi, 'in1', label='in')
c2 = Connection(pi, 'out1', si, 'in1', label='out')
nw.add_conns(c1, c2)

# define parameters
pi.set_attr(pr=0.975, Q=0, design=['pr'], L=100, D='var', ks=5e-5)
c1.set_attr(fluid={'ethanol': 1}, m=10, T=30, p=3)


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

si_out = Bus('out')
si_out.add_comps({'comp': si})


# exergoeconomic analysis
ean = ExergyAnalysis(nw, E_F=[so_in_B], E_P=[si_out], E_L=[])
exe_eco_input = {'Pipe_Z': 5e3, 'Inlet_c': 10}
ean.analyse(pamb=p_amp, Tamb=T_amb, Chem_Ex=chemexlib)
ean.evaluate_exergoeconomics(Tamb=T_amb, Exe_Eco_Costs=exe_eco_input)
ean.print_results(Exe_Eco_An=True)