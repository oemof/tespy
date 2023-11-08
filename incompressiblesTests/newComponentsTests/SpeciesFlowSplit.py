import logging

from tespy.components import HeatExchangerSimple, Source, Sink, Merge, Separator 
from tespy.tools import ComponentProperties
from tespy.connections import Connection
from tespy.networks import Network
import numpy as np

from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp

from tespy.components.newcomponents import DiabaticSimpleHeatExchanger,MergeDeltaP,SeparatorWithSpeciesSplits

logging.basicConfig(level=logging.DEBUG)

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
c1.set_attr(m=1, p=1.2, T=30, fluid={"INCOMP::Water": 0.9, "INCOMP::T66": 0.1}, mixing_rule="incompressible")

# set compositions 
c2.set_attr(fluid={"INCOMP::Water": 0.85, "INCOMP::T66": 0.15})

se.set_attr(SFS={
    'val': 0.6, 'is_set': True, 
    'split_fluid' : 'INCOMP::T66', 'split_outlet' : "out1"})
nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
print(nw.results['Connection'])

se.set_attr(SFS={
    'val': 0.6, 'is_set': True, 
    'split_fluid' : 'INCOMP::Water', 'split_outlet' : "out1"})
nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
print(nw.results['Connection'])
