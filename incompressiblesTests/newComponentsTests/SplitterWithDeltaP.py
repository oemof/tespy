import logging

from tespy.components import HeatExchangerSimple, Source, Sink, Merge, Separator 
from tespy.tools import ComponentProperties
from tespy.connections import Connection
from tespy.networks import Network
import numpy as np

from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp

from tespy.components.newcomponents import DiabaticSimpleHeatExchanger,MergeDeltaP,SeparatorWithSpeciesSplits,SplitterDeltaP

logging.basicConfig(level=logging.DEBUG)

# %%

# caution, must write "Water" (capital W) in INCOMP backend -> CoolProp bug? Intentional?
fluids = ["INCOMP::Water", "INCOMP::T66"]
nw = Network(fluids=fluids, p_unit="bar", T_unit="C")

so = Source("Source")
se = SplitterDeltaP("Splitter",num_out=2)
si1 = Sink("Sink 1")
si2 = Sink("Sink 2")

c1 = Connection(so, "out1", se, "in1", label="1")
c2 = Connection(se, "out1", si1, "in1", label="2")
c3 = Connection(se, "out2", si2, "in1", label="3")

nw.add_conns(c1, c2, c3)

# set some generic data for starting values
c1.set_attr(m=1, p=1.2, T=30, fluid={"INCOMP::Water": 0.9, "INCOMP::T66": 0.1}, mixing_rule="incompressible")
c2.set_attr(m=0.6)

c2.set_attr(p=1.1)
c3.set_attr(p=1.0)

# add some guess values
#c2.set_attr(m0=0.5,p0=1.2,h0=1e5,T0=50,fluid0={"Water": 0.5, "T66": 0.5})
#c3.set_attr(m0=0.5,p0=1.2,h0=1e5,T0=50,fluid0={"Water": 0.5, "T66": 0.5})

nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
print(nw.results['Connection'])

# use delta P
c2.set_attr(p=None)
c3.set_attr(p=None)
se.set_attr(deltaP=0.25)

nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
print(nw.results['Connection'])
