# %%
import logging

from tespy.components import HeatExchangerSimple, Source, Sink, Merge
from tespy.connections import Connection
from tespy.networks import Network
import numpy as np

from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp

from tespy.components.newcomponents import DiabaticSimpleHeatExchanger,MergeDeltaP,SeparatorWithSpeciesSplits

logging.basicConfig(level=logging.DEBUG)

# %%

# caution, must write "Water" (capital W) in INCOMP backend -> CoolProp bug? Intentional?

nw = Network(p_unit="bar", T_unit="C")

so = Source("Source")
so2 = Source("Source2")

me = MergeDeltaP("Merge")
si = Sink("Sink")

c1 = Connection(so, "out1", me, "in1", label="2")
c3 = Connection(so2, "out1", me, "in2", label="3")
c4 = Connection(me, "out1", si, "in1", label="4")

nw.add_conns(c1, c3, c4)

# set some generic data for starting values
c1.set_attr(m=1, p=2.1, h=0.5e5, fluid={"INCOMP::Water": 0.9, "INCOMP::T66": 0.1}, mixing_rule="incompressible")
# mix with pure water
c3.set_attr(m=0.05, p=2.2, h=0.5e5, fluid={"INCOMP::Water": 1, "INCOMP::T66": 0})

# set pressure ratios of heater and merge
me.set_attr(deltaP=1)
#c4.set_attr(p=1)

nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
print(nw.results['Connection'])


me.set_attr(deltaP=None)
c4.set_attr(p=1)

nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
print(nw.results['Connection'])