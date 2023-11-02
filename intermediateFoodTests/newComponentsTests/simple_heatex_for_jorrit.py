# %% two fluids
import logging

from tespy.components import HeatExchangerSimple, Source, Sink
from tespy.connections import Connection
from tespy.networks import Network

# fluid and network
fluids = ["INCOMP::Water", "INCOMP::T66"]
nw = Network(fluids=fluids, p_unit="bar", T_unit="C")

# components
so = Source("Source")
he = HeatExchangerSimple("Heater")
si = Sink("Sink")

# connections
c1 = Connection(so, "out1", he, "in1", label="1")
c2 = Connection(he, "out1", si, "in1", label="4")
nw.add_conns(c1, c2)

# set some conditions on connections
c1.set_attr(m=1, p=1.2, T=30, fluid={"INCOMP::Water": 0.9, "INCOMP::T66": 0.1})
c2.set_attr(T=50)

# set some conditions on component
he.set_attr(pr=1)

nw.solve("design")
nw.print_results()

# %% one fluid

import logging

from tespy.components import HeatExchangerSimple, Source, Sink
from tespy.connections import Connection
from tespy.networks import Network

# fluid and network
fluids = ["INCOMP::Water"]
nw = Network(fluids=fluids, p_unit="bar", T_unit="C")

# components
so = Source("Source")
he = HeatExchangerSimple("Heater")
si = Sink("Sink")

# connections
c1 = Connection(so, "out1", he, "in1", label="1")
c2 = Connection(he, "out1", si, "in1", label="2")
nw.add_conns(c1, c2)

# set some conditions on connections
c1.set_attr(m=1, p=1.2, T=30, fluid={"INCOMP::Water": 1})
c2.set_attr(T=50)

# set some conditions on component
he.set_attr(pr=1)

nw.solve("design")
nw.print_results()
