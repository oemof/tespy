import logging

from tespy.components import HeatExchangerSimple, Source, Sink, Merge, Separator 
from tespy.tools import ComponentProperties
from tespy.connections import Connection
from tespy.networks import Network
import numpy as np

from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp

from tespy.components.newComponents import \
    DiabaticSimpleHeatExchanger,MergeDeltaP,SeparatorWithSpeciesSplits, \
        SeparatorWithSpeciesSplitsDeltaT, SeparatorWithSpeciesSplitsDeltaTDeltaP

from tespy.components.newAdvancedComponents import TwoStreamEvaporator


logging.basicConfig(level=logging.DEBUG)

# %%

# caution, must write "Water" (capital W) in INCOMP backend -> CoolProp bug? Intentional?
fluids = ["HEOS::Water", "INCOMP::FoodFat", "INCOMP::FoodProtein"]
nw = Network(m_unit='kg / s', p_unit='bar', T_unit='C',h_unit='kJ / kg', h_range=[-1e2,4e3], iterinfo=True)

so1 = Source("Source 1")
so2 = Source("Source 2")
se = TwoStreamEvaporator("Separator",num_in=2,num_out=3)
si1 = Sink("Sink 1")
si2 = Sink("Sink 2")
si3 = Sink("Sink 3")

c1 = Connection(so1, "out1", se, "in1", label="1")
c2 = Connection(se, "out1", si1, "in1", label="2")
c3 = Connection(se, "out2", si2, "in1", label="3")

c4 = Connection(so2, "out1", se, "in2", label="4")
c5 = Connection(se, "out3", si3, "in1", label="5")

nw.add_conns(c1, c2, c3, c4, c5)

# set global guess values 
m0 = 74.6    # transform unit at some point [this is kt/yr]
h0 = 323.3        # global guess value in kJ/kg
p0 = 1        # global guess value in bar

# for c in nw.conns['object']:
#     n_fl = 3 # len(nw.fluids)
#     c.set_attr(m0=70,h0=h0,p0=p0,fluid0={'HEOS::Water': 0.92, 'INCOMP::FoodFat': 0.023, 'INCOMP::FoodProtein': 0.052})

# set some generic data for starting values
c1.set_attr(m=m0, p=p0, T=80, fluid={'HEOS::Water': 0.942, "INCOMP::FoodFat": 0.004, "INCOMP::FoodProtein": 0.054}, mixing_rule="incompressible")
c2.set_attr(p=p0, T=40, force_state='l', fluid={"INCOMP::FoodProtein": 0.075})
c3.set_attr(x=1, T=40, force_state='g', fluid={'HEOS::Water': 1.0, "INCOMP::FoodFat": 0.0, "INCOMP::FoodProtein": 0.0})

# Now it is possible to set the temperatures out of the separator differently
c4.set_attr(m=25, p=p0, x=1, fluid={'HEOS::Water': 1.0}, mixing_rule="incompressible")
c5.set_attr(p=p0)

nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
print(nw.results['Connection'])

c2.set_attr(fluid={"INCOMP::FoodProtein": None})
c5.set_attr(p=p0,x=0)
nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
print(nw.results['Connection'])


c5.set_attr(x=None)
se.set_attr(Q=5e7)
nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
print(nw.results['Connection'])


c5.set_attr(x=None)
se.set_attr(Q=None)
se.set_attr(kA=8e5)
nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
print(nw.results['Connection'])

c5.set_attr(x=None)
se.set_attr(Q=None)
se.set_attr(kA=None)
se.set_attr(KPI=6e5)
nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
print(nw.results['Connection'])


se.set_attr(KPI=None)
se.set_attr(kA=7e5)
se.set_attr(dTo=0)
c2.set_attr(fluid={"INCOMP::FoodProtein": 0.075})
c2.set_attr(T=None)
c3.set_attr(T=None)

nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
print(nw.results['Connection'])



# mass balance mode

c2.set_attr(fluid={"INCOMP::FoodProtein": 0.075})
c2.set_attr(p=p0, x=None, T=None, force_state=None)
c3.set_attr(p=p0, x=None, T=None, force_state=None)

c5.set_attr(x=None, force_state=None)
se.set_attr(Q=None)
se.set_attr(kA=None)
se.set_attr(dTo=None)
se.set_attr(KPI=None)
se.set_attr(deltaH=0)
nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
print(nw.results['Connection'])
