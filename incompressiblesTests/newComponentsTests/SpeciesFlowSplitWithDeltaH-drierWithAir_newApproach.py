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
        SeparatorWithSpeciesSplitsDeltaT

from tespy.components.newAdvancedComponents import TwoStreamDrier

logging.basicConfig(level=logging.DEBUG)

# %%

# caution, must write "Water" (capital W) in INCOMP backend -> CoolProp bug? Intentional?
#fluids = ["INCOMP::Water", "INCOMP::T66"]
nw = Network(m_unit='kg / s', p_unit='bar', T_unit='C',h_unit='kJ / kg', h_range=[-1e2,4e3], iterinfo=True)

so = Source("Source")
soAir = Source("SourceAir")
se = TwoStreamDrier("Separator",num_out=2,num_in=2)
#se = SeparatorWithSpeciesSplits("Separator") #,num_out=2)
si1 = Sink("Sink 1")
si2 = Sink("Sink 2")

c1 = Connection(so, "out1", se, "in1", label="1")
c2 = Connection(se, "out1", si1, "in1", label="2") # vapor
c3 = Connection(se, "out2", si2, "in1", label="3") # liquid
c4 = Connection(soAir, "out1", se, "in2", label="4")

nw.add_conns(c1, c2, c3, c4)

# for c in nw.conns['object']:
#     c.set_attr(m0=1, h0=100, p0=1.2)

# set some generic data for starting values
c1.set_attr(m=1, p=1.0, T=50, fluid={"HEOS::Water": 0.9, "INCOMP::T66": 0.1, "HEOS::Air": 0}, mixing_rule="incompressible")
c4.set_attr(m=50, p=1.0, T=80, fluid={"HEOS::Water": 0, "INCOMP::T66": 0, "HEOS::Air": 1}, mixing_rule="incompressible")

c3.set_attr(fluid={"HEOS::Water": 0.08, "HEOS::Air": 0})
c2.set_attr(fluid={"INCOMP::T66": 0})


#c3.set_attr(p=1.2,T=60,force_state='g')
c3.set_attr(p=1.0)
c2.set_attr(p=1.0)

se.set_attr(deltaH=0)

nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
print(nw.results['Connection'])


se.set_attr(deltaH=None)
c2.set_attr(p=1.0,T=None,force_state='g')
c3.set_attr(p=1.0,T=None,force_state='l') 
se.set_attr(dTwbProd=0)
nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
print(nw.results['Connection'])


c3.set_attr(fluid={"HEOS::Water": None})
se.set_attr(Q=se.Q.val)
se.set_attr(Q=1e6)

nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
print(nw.results['Connection'])


c3.set_attr(fluid={"HEOS::Water": None})
se.set_attr(Q=None)
se.set_attr(KPI=2.250e6)

nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
print(nw.results['Connection'])


c3.set_attr(fluid={"HEOS::Water": None})
se.set_attr(Q=None)
se.set_attr(KPI=None)
se.set_attr(kA=6.0e+04)

nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
print(nw.results['Connection'])

c3.set_attr(fluid={"HEOS::Water": None})
se.set_attr(Q=None)
se.set_attr(KPI=None)
se.set_attr(kA=None)
se.set_attr(WBeff=0.4)

nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
print(nw.results['Connection'])



