import logging

from tespy.components import HeatExchangerSimple, Source, Sink, Merge, Separator 
from tespy.tools import ComponentProperties
from tespy.connections import Connection
from tespy.networks import Network
import numpy as np

from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp

from tespy.components.newcomponents import DiabaticSimpleHeatExchanger,MergeWithPressureLoss,SeparatorWithSpeciesSplits

# %%

# caution, must write "Water" (capital W) in INCOMP backend -> CoolProp bug? Intentional?
fluids = ["INCOMP::Water", "INCOMP::T66"]
nw = Network(fluids=fluids, p_unit="bar", T_unit="C")

so = Source("Source")
se = SeparatorWithSpeciesSplits("Separator")
si1 = Sink("Sink 1")
si2 = Sink("Sink 2")

c1 = Connection(so, "out1", se, "in1", label="1")
c2 = Connection(se, "out1", si1, "in1", label="2")
c3 = Connection(se, "out2", si2, "in1", label="3")

nw.add_conns(c1, c2, c3)

# set some generic data for starting values
c1.set_attr(m=1, p=1.2, T=30, fluid={"Water": 0.9, "T66": 0.1})

# set compositions 
c2.set_attr(fluid={"Water": 0.85, "T66": 0.15})
# or others 
#c3.set_attr(fluid={"Water": 0.85, "T66": 0.15})
# or others
#c2.set_attr(fluid={"Water": 0.85})
#c3.set_attr(fluid={"Water": 1})

# This one produce error because T66 is 0 and the equation cannot be solved.. 
#c2.set_attr(fluid={"Water": 1, "T66": 0.0})

# specify this one to avoid using the species flow split : SFS  
#c2.set_attr(m=0.5)
# set the species flow split, specify the fluid and the outlet too.. (we might need some checks of this)
se.set_attr(SFS={
    'val': 0.6, 'is_set': True, 
    'split_fluid' : 'T66', 'split_outlet' : "out1"})

# add some guess values
c2.set_attr(m0=0.5,p0=1.2,h0=1e5,T0=50,fluid0={"Water": 0.5, "T66": 0.5})
c3.set_attr(m0=0.5,p0=1.2,h0=1e5,T0=50,fluid0={"Water": 0.5, "T66": 0.5})

nw.solve("design")
nw.print_results()

print(nw.results['Connection'])

m_T66_c1 = c1.m.val * c1.fluid.val['T66']
m_T66_c2 = c2.m.val * c2.fluid.val['T66']

print(f"\n Species flow split is {m_T66_c2/m_T66_c1}")

print(f"\n heat flows are  {se.Q.val}")
print(f"\n")


