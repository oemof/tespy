from tespy.networks import Network
from tespy.components import Source, Sink, Merge, Splitter, HeatExchangerSimple, Compressor, Pipe, Splitter, Turbine
from tespy.connections import Connection


nwk = Network(T_unit="C", p_unit="bar", h_unit='kJ / kg')

source = Source('source123')
source2 = Source('source1234')
merge = Merge('merge123')
component1 = HeatExchangerSimple('comp1', pr=1)
splitter = Splitter('splitter123')
component2 = HeatExchangerSimple('comp2')
sink = Sink('sink123')

c1 = Connection(source, 'out1', merge, 'in1', p=1, h=200, m=10, fluid={'R134a': 1})
c2 = Connection(merge, 'out1', component1, 'in1')
c3 = Connection(component1, 'out1', splitter, 'in1', h=180)
c4 = Connection(splitter, 'out1', component2, 'in1', m=1)
c5 = Connection(component2, 'out1', merge, 'in2', h=170)
c6 = Connection(splitter, 'out2', sink, 'in1')
nwk.add_conns(c1, c2, c3, c4, c5, c6)

nwk.solve('design')

c3.set_attr(h=None)
c5.set_attr(h=None)
component1.set_attr(Q=-1000)
component2.set_attr(Q=-500)

nwk.solve('design')

so1 = Source("air")
so2 = Source("Other gases")
m1 = Merge("gas mixing")
p1 = Pipe("test", pr=1, Q=0)
sp1 = Splitter("Splitter")
t1 = Turbine("Turbine", pr=.1, eta_s=.8)
cp1 = Compressor("Compressor", pr=10, eta_s=.8)
si1 = Sink("Sink1")
si2 = Sink("Sink2")

c21 = Connection(so1, "out1", m1, "in1", label="21", fluid={"N2": 0.76, "O2": 0.23, "Ar": 0.01}, m=10, T=400, p=1, mixing_rule="ideal-cond")
c22 = Connection(so2, "out1", m1, "in2", label="22", fluid={"H2O": 1}, m=.5, T=400)
c23 = Connection(m1, "out1", p1, "in1", label="23")
c24 = Connection(p1, "out1", sp1, "in1", label="24")
c25 = Connection(sp1, "out1", t1, "in1", label="25")
c26 = Connection(t1, "out1", si1, "in1", label="26")
c27 = Connection(sp1, "out2", cp1, "in1", label="27", m=4)
c28 = Connection(cp1, "out1", si2, "in1", label="28")

nwk.add_conns(c21, c22, c23, c24, c25, c26, c27, c28)

nwk.solve("design")