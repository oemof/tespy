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
        SeparatorWithSpeciesSplitsDeltaH,SeparatorWithSpeciesSplitsDeltaTDeltaPDeltaH

logging.basicConfig(level=logging.DEBUG)

# %%

# caution, must write "Water" (capital W) in INCOMP backend -> CoolProp bug? Intentional?
fluids = ["INCOMP::Water", "INCOMP::T66"]
nw = Network(fluids=fluids, m_unit='kg / s', p_unit='bar', T_unit='C',h_unit='kJ / kg', h_range=[-1e2,4e3], iterinfo=True)

so = Source("Source")
se = SeparatorWithSpeciesSplitsDeltaTDeltaPDeltaH("Separator",num_out=2)
#se = SeparatorWithSpeciesSplits("Separator") #,num_out=2)
si1 = Sink("Sink 1")
si2 = Sink("Sink 2")

c1 = Connection(so, "out1", se, "in1", label="1")
c2 = Connection(se, "out1", si1, "in1", label="2")
c3 = Connection(se, "out2", si2, "in1", label="3")

nw.add_conns(c1, c2, c3)

# set some generic data for starting values
c1.set_attr(m=1, p=1.0, h=1e2, fluid={"HEOS::Water": 0.9, "INCOMP::T66": 0.1}, mixing_rule="incompressible")
#c1.set_attr(m=1, p=1.2, T=50, fluid={"INCOMP::Water": 0.9, "INCOMP::T66": 0.1})
#c1.set_attr(T0=10) # it seems guess values are SI 

c2.set_attr(fluid={"HEOS::Water": 1, "INCOMP::T66": 0.0},force_state='g')
c3.set_attr(fluid={"HEOS::Water": 0.08},force_state='l')

c2.set_attr(p=1.0,T=100)
c3.set_attr(p=1.0,T=100)

nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
print(nw.results['Connection'])


se.set_attr(Q=2.0e6)
c3.set_attr(fluid={"HEOS::Water": None})
# c3.set_attr(fluid0={"INCOMP::Water": 0.08, "INCOMP::T66": 0.92})
# c3.set_attr(T0=100)
#c3.set_attr(T=None)
nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
print(nw.results['Connection'])





se.set_attr(Q=None)
se.set_attr(KPI=1.5e5)
nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
print(nw.results['Connection'])





import sys 
sys.exit()









print(nw.results['Connection'])
m_T66_c1 = c1.m.val * c1.fluid.val['T66']
m_T66_c2 = c2.m.val * c2.fluid.val['T66']
print(f"\n Species flow split is {m_T66_c2/m_T66_c1}")
print(f"\n")

se.set_attr(deltaH=0.0)
c2.set_attr(T=None)
c3.set_attr(T=None)


nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
print(nw.results['Connection'])
m_T66_c1 = c1.m.val * c1.fluid.val['T66']
m_T66_c2 = c2.m.val * c2.fluid.val['T66']
print(f"\n Species flow split is {m_T66_c2/m_T66_c1}")
print(f"\n")


