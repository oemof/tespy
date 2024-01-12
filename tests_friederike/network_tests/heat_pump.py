from tespy.networks import Network
from tespy.components import (CycleCloser, Compressor, Valve, SimpleHeatExchanger, HeatExchanger, Source, Sink)
from tespy.connections import Connection, Bus
from tespy.tools import ExergyAnalysis
from tespy.tools.helpers import get_chem_ex_lib
chemexlib = get_chem_ex_lib("Ahrendts")


nw = Network(T_unit='C', p_unit='bar', h_unit='kJ / kg')

# components
cc = CycleCloser('cycle closer')
ev = SimpleHeatExchanger('evaporator')
cp = Compressor('compressor')
co = HeatExchanger('condenser')
va = Valve('expansion valve')
con_so = Source('consumer in')
con_si = Sink('consumer out')

# connections
c1 = Connection(cc, 'out1', ev, 'in1', label='1')
c2 = Connection(ev, 'out1', cp, 'in1', label='2')
c3 = Connection(cp, 'out1', co, 'in1', label='3')
c4 = Connection(co, 'out1', va, 'in1', label='4')
c0 = Connection(va, 'out1', cc, 'in1', label='0')

c11 = Connection(con_so, 'out1', co, 'in2', label='11')
c12 = Connection(co, 'out2', con_si, 'in1', label='22')
nw.add_conns(c1, c2, c3, c4, c0, c11, c12)

# define parameters
co.set_attr(pr1=0.98, pr2=1)
ev.set_attr(pr=0.98)
cp.set_attr(eta_s=0.85)

c2.set_attr(T=20, x=1, fluid={'CO2': 1})            # war R134a aber das gibt es nicht im Ahrendts Stoffmodell
c4.set_attr(T=80, x=0)

c11.set_attr(T=20, p=1, m=1, fluid={'H2O': 1})
c12.set_attr(T=40)

# define necessary busses
power_in = Bus('power_in')      # compressor needs power (fuel of the heat pump)
power_in.add_comps({'comp': cp})
nw.add_busses(power_in)

heat_out = Bus('heat_out')      # condenser removes heat (product of the heat pump)
heat_out.add_comps({'comp': con_so, 'base': 'bus'}, {'comp': con_si})
nw.add_busses(heat_out)

# solve
nw.solve(mode='design')
nw.print_results()


""" +++ exergy analysis +++ """
# define ambient
T_amb = 28
p_amb = 1.013

# exergy analysis
ean = ExergyAnalysis(nw, E_F=[power_in], E_P=[heat_out], E_L=[])
ean.analyse(pamb=p_amb, Tamb=T_amb, Chem_Ex=chemexlib)
ean.print_results()


"""
# exergy and exergoeconomic analysis
exe_eco_input = {'condenser_Z': 50, 'evaporator_Z': 20, 'expansion valve_Z': 40, 'compressor_Z': 100, 'power in_c': 0.003}
ean = ExergyAnalysis(nw, E_F=[power_in], E_P=[heat_out], E_L=[])
ean.analyse(pamb=p_amb, Tamb=T_amb, Exe_Eco_An = True, Exe_Eco_Costs = exe_eco_input)
ean.print_results(Exe_Eco_An=True)
"""

# before changing condenser from SimpleHeatExchanger ot HeatExchanger in order
# to be able to assign an exergetic product to the system

"""
from tespy.networks import Network
from tespy.components import (CycleCloser, Compressor, Valve, SimpleHeatExchanger)
from tespy.connections import Connection, Bus
from tespy.tools import ExergyAnalysis


nw = Network(T_unit='C', p_unit='bar', h_unit='kJ / kg')

# components
cc = CycleCloser('cycle closer')
ev = SimpleHeatExchanger('evaporator')
cp = Compressor('compressor')
co = SimpleHeatExchanger('condenser')
va = Valve('expansion valve')

# connections
c1 = Connection(cc, 'out1', ev, 'in1', label='1')
c2 = Connection(ev, 'out1', cp, 'in1', label='2')
c3 = Connection(cp, 'out1', co, 'in1', label='3')
c4 = Connection(co, 'out1', va, 'in1', label='4')
c0 = Connection(va, 'out1', cc, 'in1', label='0')

nw.add_conns(c1, c2, c3, c4, c0)

# define parameters
co.set_attr(pr=0.98, Q=-1e6)
ev.set_attr(pr=0.98)
cp.set_attr(eta_s=0.85)

c2.set_attr(T=20, x=1, fluid={'R134a': 1})
c4.set_attr(T=80, x=0)

# define necessary busses
power_in = Bus('power_in')      # compressor needs power (fuel of the heat pump)
power_in.add_comps({'comp': cp})
nw.add_busses(power_in)

heat_out = Bus('heat_out')      # condenser removes heat (product of the heat pump)
heat_out.add_comps({'comp': co})
nw.add_busses(heat_out)

# solve
nw.solve(mode='design')
nw.print_results()


+++ exergy analysis +++ 
# define ambient
T_amb = 28
p_amb = 1.013

# exergy analysis
ean = ExergyAnalysis(nw, E_F=[power_in], E_P=[heat_out], E_L=[])
ean.analyse(pamb=p_amb, Tamb=T_amb)
ean.print_results()


# exergy and exergoeconomic analysis
exe_eco_input = {'condenser_Z': 50, 'evaporator_Z': 20, 'expansion valve_Z': 40, 'compressor_Z': 100, 'power in_c': 0.003}
ean = ExergyAnalysis(nw, E_F=[power_in], E_P=[heat_out], E_L=[])
ean.analyse(pamb=p_amb, Tamb=T_amb, Exe_Eco_An = True, Exe_Eco_Costs = exe_eco_input)
ean.print_results(Exe_Eco_An=True)

"""