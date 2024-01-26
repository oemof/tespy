from tespy.networks import Network
from tespy.components import (Source, Sink, SimpleHeatExchanger)
from tespy.connections import Connection, Bus
from tespy.tools import analyses, ExergyAnalysis
import numpy as np

# network
nw = Network(T_unit="C", p_unit="bar", h_unit="kJ / kg")

# components
hx = SimpleHeatExchanger("My Heat Exchanger")
so = Source("My Source")
si = Sink("My Sink")

# connections
c1 = Connection(so, 'out1', hx, 'in1', label="My Inlet")
c2 = Connection(hx, 'out1', si, 'in1', label="My Outlet")
nw.add_conns(c1, c2)

# define parameters
hx.set_attr(pr=0.95, Q=1000000)
c1.set_attr(fluid={'Water': 1.0}, T=350, p=10, m=40)

# solve network
nw.solve(mode='design')
nw.print_results()

""" +++ exergy analysis +++ """
# define ambient
T_amb = 20
p_amb = 1

# define busses
heat_in_B = Bus('Heat In')
heat_in_B.add_comps({'comp': hx, 'base': 'bus'})

mass_flow_B = Bus('water flow')
mass_flow_B.add_comps({'comp': so, 'base': 'bus'}, {'comp': si})

nw.add_busses(heat_in_B)

# exergy and exergoeconomic analysis
exe_eco_input = {'My Heat Exchanger_Z': 100, 'My Source_c': 2e-4, 'Heat In_c': 5e-4}
ean = ExergyAnalysis(nw, E_F=[heat_in_B], E_P=[mass_flow_B], E_L=[])
ean.analyse(pamb=p_amb, Tamb=T_amb)
ean.evaluate_exergoeconomics(Exe_Eco_Costs=exe_eco_input, Tamb=T_amb)
ean.print_results(Exe_Eco_An=True)