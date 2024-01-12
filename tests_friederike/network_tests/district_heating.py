from tespy.networks import Network
from tespy.components import (CycleCloser, Pipe, Pump, Valve, SimpleHeatExchanger, HeatExchanger, Source, Sink)
from tespy.connections import Connection, Bus, Ref
from tespy.tools import ExergyAnalysis

nw = Network(T_unit='C', p_unit='bar', h_unit='kJ / kg')

# central heating plant
hs = SimpleHeatExchanger('heat source')
cc = CycleCloser('cycle closer')
pu = Pump('feed pump')

# consumer
cons = HeatExchanger('consumer')
val = Valve('control valve')
con_so = Source("consumer in")
con_si = Sink("consumer out")

# pipes
pipe_feed = Pipe('feed pipe')
pipe_return = Pipe('return pipe')

# connections
c0 = Connection(cc, "out1", hs, "in1", label="0")
c1 = Connection(hs, "out1", pu, "in1", label="1")
c2 = Connection(pu, "out1", pipe_feed, "in1", label="2")
c3 = Connection(pipe_feed, "out1", cons, "in1", label="3")
c4 = Connection(cons, "out1", val, "in1", label="4")
c5 = Connection(val, "out1", pipe_return, "in1", label="5")
c6 = Connection(pipe_return, "out1", cc, "in1", label="6")

c11 = Connection(con_so, "out1", cons, "in2", label="11")
c12 = Connection(cons, "out2", con_si, "in1", label="12")
nw.add_conns(c0, c1, c2, c3, c4, c5, c6, c11, c12)


# define parameters
cons.set_attr(Q=-10000)
hs.set_attr(pr=1)
pu.set_attr(eta_s=0.75)
pipe_feed.set_attr(Q=-250, pr=0.98)
pipe_return.set_attr(Q=-200, pr=0.98)

c1.set_attr(T=90, p=10, fluid={'water': 1.0})
c2.set_attr(p=13)
c4.set_attr(T=65, p=Ref(c3, 1, 0))

c11.set_attr(T=25, p=1, fluid={'water': 1.0})
c12.set_attr(T=40, p=1)

""" +++ exergy analysis +++ """
# define ambient
T_amb = 5
p_amp = 1

# define busses
consumer = Bus("consumer bus")
consumer.add_comps(
    {'comp': con_so, 'base': 'bus'},
    {'comp': con_si})
heat_source = Bus("heat source bus")
heat_source.add_comps({"comp": hs})
pump = Bus("feed pump bus")
pump.add_comps({"comp": pu})
# nw.add_busses(consumer, heat_source, pump)

# solve network
nw.solve(mode='design')
nw.print_results()

# exergy analysis
ean = ExergyAnalysis(nw, E_F=[heat_source, pump], E_P=[consumer], E_L=[])
ean.analyse(pamb=p_amp, Tamb=T_amb)
ean.print_results()


'''
# exergy and exergoeconomic analysis
exe_eco_input = {'heat source_Z': 100, 'consumer_Z': 80, 'feed pipe_Z': 70, 'return pipe_Z': 120}
ean = ExergyAnalysis(nw, E_F=[heat_source, pump], E_P=[consumer], E_L=[])
ean.analyse(pamb=p_amp, Tamb=T_amb, Exe_Eco_An=True, Exe_Eco_Costs=exe_eco_input)
ean.print_results(Exe_Eco_An=True)
'''


# following code is with consumer as SimpleHeatExchanger instead of HeatExchanger
# then exergy analysis doesn't work since product is removal of Exergy from system

'''
from tespy.networks import Network
from tespy.components import (CycleCloser, Pipe, Pump, Valve, SimpleHeatExchanger)
from tespy.connections import Connection, Bus
from tespy.tools import ExergyAnalysis

nw = Network(T_unit='C', p_unit='bar', h_unit='kJ / kg')

# central heating plant
hs = SimpleHeatExchanger('heat source')
cc = CycleCloser('cycle closer')
pu = Pump('feed pump')

# consumer
cons = SimpleHeatExchanger('consumer')
val = Valve('control valve')

# pipes
pipe_feed = Pipe('feed pipe')
pipe_return = Pipe('return pipe')

# connections
c0 = Connection(cc, "out1", hs, "in1", label="0")
c1 = Connection(hs, "out1", pu, "in1", label="1")
c2 = Connection(pu, "out1", pipe_feed, "in1", label="2")
c3 = Connection(pipe_feed, "out1", cons, "in1", label="3")
c4 = Connection(cons, "out1", val, "in1", label="4")
c5 = Connection(val, "out1", pipe_return, "in1", label="5")
c6 = Connection(pipe_return, "out1", cc, "in1", label="6")
nw.add_conns(c0, c1, c2, c3, c4, c5, c6)


# define parameters
cons.set_attr(Q=-10000, pr=0.98)
hs.set_attr(pr=1)
pu.set_attr(eta_s=0.75)
pipe_feed.set_attr(Q=-250, pr=0.98)
pipe_return.set_attr(Q=-200, pr=0.98)

c1.set_attr(T=90, p=10, fluid={'water': 1.0})
c2.set_attr(p=13)
c4.set_attr(T=65)


""" +++ exergy analysis +++ """
# define ambient
T_amb = 5
p_amp = 1

# define busses
consumer = Bus("consumer bus")
consumer.add_comps({"comp": cons})
heat_source = Bus("heat source bus")
heat_source.add_comps({"comp": hs})
pump = Bus("feed pump bus")
pump.add_comps({"comp": pu})
# nw.add_busses(consumer, heat_source, pump)

# solve network
nw.solve(mode='design')
nw.print_results()

# exergy analysis
ean = ExergyAnalysis(nw, E_F=[heat_source, pump], E_P=[consumer], E_L=[])
ean.analyse(pamb=p_amp, Tamb=T_amb)
ean.print_results()

# exergy and exergoeconomic analysis
exe_eco_input = {'heat source_Z': 100, 'consumer_Z': 80, 'feed pipe_Z': 70, 'return pipe_Z': 120}
ean = ExergyAnalysis(nw, E_F=[heat_source, pump], E_P=[consumer], E_L=[])
ean.analyse(pamb=p_amp, Tamb=T_amb, Exe_Eco_An=True, Exe_Eco_Costs=exe_eco_input)
ean.print_results(Exe_Eco_An=True)


'''
