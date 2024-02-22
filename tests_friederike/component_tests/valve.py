from tespy.components import Sink, Source, Valve
from tespy.connections import Connection, Bus
from tespy.tools.analyses import ExergyAnalysis
from tespy.networks import Network
from tespy.tools.helpers import get_chem_ex_lib
chemexlib = get_chem_ex_lib("Ahrendts")

testCase = 3

# network
nw = Network(p_unit='bar', T_unit='C', iterinfo=False)

# components
so = Source('source')
si = Sink('sink')
v = Valve('valve')

# connections
so_v = Connection(so, 'out1', v, 'in1', label='in')
v_si = Connection(v, 'out1', si, 'in1', label='out')
nw.add_conns(so_v, v_si)

# define parameters
v.set_attr(offdesign=['zeta'])
so_v.set_attr(fluid={'CH4': 1}, m=1, design=['m'])
match testCase:
    case 1:         # Tin, Tout > T0  -> dissipativ
        so_v.set_attr(p=1.5)
        v_si.set_attr(T=70, p=1.2) # non-dissipative
    case 2:         # Tin > T0, Tout <= T0
        so_v.set_attr(p=1.5)
        v_si.set_attr(T=49.9, p=1.2) # non-dissipative
    case 3:         # Tin, Tout <= T0
        so_v.set_attr(p=1.5)
        v_si.set_attr(T=20, p=1.2) # non-dissipative

# solve
nw.solve('design')
nw.print_results()


""" +++ exergy analysis +++ """
# define ambient
p_amb = 1
T_amb = 50

# define busses
inlet = Bus('in')
inlet.add_comps({'comp': so, 'base': 'bus'})

outlet = Bus('out')
outlet.add_comps({'comp': si, 'base': 'component'})

# exergy and exergoeconomic analysis
exe_eco_input = {'source_c': 0.02, 'valve_Z': 5e3}
ean = ExergyAnalysis(nw, E_F=[inlet], E_P=[outlet], E_L=[])
ean.analyse(pamb=p_amb, Tamb=T_amb, Chem_Ex=chemexlib)
ean.evaluate_exergoeconomics(Exe_Eco_Costs=exe_eco_input, Tamb=T_amb)
ean.print_results(Exe_Eco_An=True)
