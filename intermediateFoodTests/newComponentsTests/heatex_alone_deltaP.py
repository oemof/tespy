# %%

import logging


from tespy.components import SimpleHeatExchanger, Source, Sink, Merge, Separator 
from tespy.tools import ComponentProperties
from tespy.connections import Connection
from tespy.networks import Network
import numpy as np

from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp

from tespy.components.newcomponents import SimpleHeatExchangerDeltaP, SimpleHeatExchangerDeltaPLossFactor,MergeWithPressureLoss,SeparatorWithSpeciesSplits

logging.basicConfig(level=logging.DEBUG)

# %%

# caution, must write "Water" (capital W) in INCOMP backend -> CoolProp bug? Intentional?
#fluids = ["INCOMP::FoodWater", "INCOMP::FoodProtein"]
#nw = Network(p_unit="bar", T_unit="C", iterinfo=True)
nw = Network(m_unit='kg / s', p_unit='bar', T_unit='C',h_unit='kJ / kg', h_range=[-1e2,4e2], iterinfo=True)


so = Source("Source")
#  Variant 2: Q is m (h_2 - h_1), Q_total is taking efficiency into account and represents the heat transfer over system
# boundary. For heat transfer into the system: Q = Q_total * eta, for heat transfer from the system: Q_total = Q * eta

he = SimpleHeatExchangerDeltaP("Heater")
#he = SimpleHeatExchangerDeltaPLossFactor("Heater")


si = Sink("Sink")

c1 = Connection(so, "out1", he, "in1", label="1")
c2 = Connection(he, "out1", si, "in1", label="4")

nw.add_conns(c1, c2)

# set some generic data for starting values
c1.set_attr(fluid={'INCOMP::Water': 0.80,'INCOMP::PHE': 0.15,'INCOMP::S800': 0.05}, m=1, T=30, p=1, mixing_rule="incompressible")
#c1.set_attr(m=1, p=2.2, T=30, fluid={"INCOMP::FoodWater": 0.9, "INCOMP::FoodProtein": 0.1}, mixing_rule="incompressible")
c2.set_attr(T=95) #,p=1)

# set pressure ratios of heater and merge
he.set_attr(deltaP=0.1)
#he.set_attr(pr=1.1)

#he.set_attr(Tamb=10)

#he.set_attr(L=10)
#he.set_attr(D=0.3)

#he.set_attr(Q=74000)
#he.set_attr(LF=0.1) # MRK so eta is (1-hlf) heat loss factor
#he.set_attr(Q_total=8.16e4) # MRK so eta is (1-hlf) heat loss factor
#he.set_attr(Q_loss=-7851.921461139966) # MRK so eta is (1-hlf) heat loss factor

nw.solve("design")
nw.print_results()

# print(nw.results['Connection'])
# he.Q.val
# he.Q_loss.val
# he.Q_total.val

# print(he.LF.val)
# print(he.Q_total.val)
# print(he.Q_loss.val)
