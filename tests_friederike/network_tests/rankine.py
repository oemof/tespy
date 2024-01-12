from tespy.networks import Network
from tespy.connections import Bus
from tespy.components import (
    CycleCloser, Pump, Condenser, Turbine, SimpleHeatExchanger, Source, Sink
)
from tespy.connections import Connection
from tespy.tools import ExergyAnalysis

nw = Network(T_unit='C', p_unit='bar', h_unit='kJ / kg')

# components
cc = CycleCloser('cycle closer')
sg = SimpleHeatExchanger('steam generator')
co = Condenser('main condenser')
tu = Turbine('steam turbine')
fp = Pump('feed pump')

cwso = Source('cooling water source')
cwsi = Sink('cooling water sink')

# connections
c1 = Connection(cc, 'out1', tu, 'in1', label='1')
c2 = Connection(tu, 'out1', co, 'in1', label='2')
c3 = Connection(co, 'out1', fp, 'in1', label='3')
c4 = Connection(fp, 'out1', sg, 'in1', label='4')
c0 = Connection(sg, 'out1', cc, 'in1', label='0')

nw.add_conns(c1, c2, c3, c4, c0)

c11 = Connection(cwso, 'out1', co, 'in2', label='11')
c12 = Connection(co, 'out2', cwsi, 'in1', label='12')

nw.add_conns(c11, c12)

# define parameters
co.set_attr(pr1=1, pr2=0.98)
sg.set_attr(pr=0.9)
tu.set_attr(eta_s=0.9)
fp.set_attr(eta_s=0.75)

c11.set_attr(T=20, p=1.2, fluid={'water': 1})
c12.set_attr(T=30)
c1.set_attr(T=600, p=150, m=10, fluid={'water': 1})
c2.set_attr(p=0.1)

# busses
power = Bus("electrical power output")
power.add_comps(
    {"comp": tu, "char": 0.97, "base": "component"},
    {"comp": fp, "char": 0.97, "base": "bus"},
)
nw.add_busses(power)

# solve
nw.solve(mode='design')
nw.print_results()

""" +++ exergy analysis +++ """
# define ambient
p_amb = 1.013
T_amb = 2.8

heat_source = Bus("heat effort")
heat_source.add_comps({"comp": sg})

heat_sink = Bus("heat loss")
heat_sink.add_comps({"comp": cwso, "base": "bus"}, {"comp": cwsi})

# exergy analysis
ean = ExergyAnalysis(network=nw, E_F=[heat_source], E_P=[power], E_L=[heat_sink])
ean.analyse(pamb=p_amb, Tamb=T_amb)
ean.print_results()

