import logging

from tespy.components import HeatExchangerSimple, Source, Sink, Merge
from tespy.connections import Connection
from tespy.networks import Network
import numpy as np

from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp

from tespy.components.newcomponents import DiabaticSimpleHeatExchanger,MergeWithPressureLoss,SeparatorWithSpeciesSplits

logging.basicConfig(level=logging.DEBUG)

# caution, must write "Water" (capital W) in INCOMP backend -> CoolProp bug? Intentional?
fluids = ["INCOMP::Water", "INCOMP::T66"]


nw = Network(fluids=fluids, p_unit="bar", T_unit="C")

so = Source("Source")
so2 = Source("Source2")

#  Variant 2: Q is m (h_2 - h_1), Q_total is taking efficiency into account and represents the heat transfer over system
# boundary. For heat transfer into the system: Q = Q_total * eta, for heat transfer from the system: Q_total = Q * eta
he = DiabaticSimpleHeatExchanger("Heater")
me = MergeWithPressureLoss("Merge")
si = Sink("Sink")

c1 = Connection(so, "out1", he, "in1", label="1")
c2 = Connection(he, "out1", me, "in1", label="2")
c3 = Connection(so2, "out1", me, "in2", label="3")
c4 = Connection(me, "out1", si, "in1", label="4")

nw.add_conns(c1, c2, c3, c4)

# set some generic data for starting values
c1.set_attr(m=1, p=1.2, h=0.5e5, fluid={"INCOMP::Water": 0.9, "INCOMP::T66": 0.1})
c2.set_attr(h=2.2e5)
# mix with pure water
c3.set_attr(m=0.05, p=1.1, h=0.5e5, fluid={"INCOMP::Water": 1, "INCOMP::T66": 0})

# set pressure ratios of heater and merge
he.set_attr(pr=0.9)
me.set_attr(deltaP=0.15)
#c2.set_attr(p=2.2)
#c4.set_attr(p=2.2)

nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()

# use temperature to make it relatable
c1.set_attr(h=None, T=30)
c2.set_attr(h=None, T=50)

nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()

# add some heat
c2.set_attr(T=None)
# # efficiency is used for postprocessing here
he.set_attr(Q=1e5, eta=0.9)
nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()

he.set_attr(Q=1e5, Q_total=1.1e5, eta=None)
nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()

he.set_attr(Q=1e5, Q_total=1.5e5)
nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()

c2.set_attr(T=50)

# impose over system boundary heat transfer (cannot be lower than actual heat transfer, efficiency value cannot be > 1!)
# In this case, efficiency decreases
he.set_attr(Q=None, Q_total=1.1e5, eta=None)

nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()

# with set efficiency, temperature cannot be set anymore
c2.set_attr(T=None)
he.set_attr(Q_total=1.1e5, eta=.5)

nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()

# now cooling instead of heating, CoolProp or TESPy have issues with freezing temperatures, so > 0°C
c2.set_attr(T=5)
he.set_attr(Q_total=None, eta=None)

nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()

he.set_attr(Q_total=-.6e5, eta="var")

nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
