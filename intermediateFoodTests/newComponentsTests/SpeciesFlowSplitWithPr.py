import logging

from tespy.components import HeatExchangerSimple, Source, Sink, Merge, Separator 
from tespy.tools import ComponentProperties
from tespy.connections import Connection
from tespy.networks import Network
import numpy as np

from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp

from tespy.components.newcomponents import \
    DiabaticSimpleHeatExchanger,MergeWithPressureLoss,SeparatorWithSpeciesSplits, \
        SeparatorWithSpeciesSplitsAndDeltaT, SeparatorWithSpeciesSplitsAndDeltaTAndPr, SeparatorWithSpeciesSplitsAndPr

# %%

# caution, must write "Water" (capital W) in INCOMP backend -> CoolProp bug? Intentional?
fluids = ["INCOMP::FoodWater", "INCOMP::FoodProtein"]
nw = Network(fluids=fluids, m_unit='kg / s', p_unit='bar', T_unit='C',h_unit='kJ / kg', h_range=[-1e2,4e3], iterinfo=True)

so = Source("Source")
se = SeparatorWithSpeciesSplitsAndPr("Separator",num_out=2)
si1 = Sink("Sink 1")
si2 = Sink("Sink 2")

c1 = Connection(so, "out1", se, "in1", label="1")
c2 = Connection(se, "out1", si1, "in1", label="2")
c3 = Connection(se, "out2", si2, "in1", label="3")

nw.add_conns(c1, c2, c3)

# set global guess values 
m0 = 1    # transform unit at some point [this is kt/yr]
h0 = 1e2        # global guess value in kJ/kg
p0 = 5        # global guess value in bar

for c in nw.conns['object']:
    n_fl = 2 # len(nw.fluids)
    c.set_attr(m0=m0,h0=h0,p0=p0,fluid0={'INCOMP::FoodWater': 1/n_fl, 'INCOMP::FoodProtein': 1/n_fl})

# set some generic data for starting values
c1.set_attr(m=1, p=5, h=h0, fluid={"INCOMP::FoodWater": 0.9, "INCOMP::FoodProtein": 0.1}, mixing_rule="incompressible")
c2.set_attr(fluid={"INCOMP::FoodWater": 0.8, "INCOMP::FoodProtein": 0.2})
#c3.set_attr(fluid={"FoodProtein": 0.1})

se.set_attr(SFS={
    'val': 0.65, 'is_set': True, 
    'split_fluid' : 'FoodProtein', 'split_outlet' : "out1"})


# Now it is possible to set the temperatures out of the separator differently
#c2.set_attr(p=5)
#c3.set_attr(p=5)

se.set_attr(deltaP=0)

# Or to use a deltaT array instead
#se.set_attr(deltaT=[-10,-20])
#se.set_attr(deltaT=[0,0])

# # add some guess values
# c2.set_attr(m0=0.5,p0=1.2,T0=50,fluid0={"FoodWater": 0.5, "FoodProtein": 0.5})
# c3.set_attr(m0=0.5,p0=1.2,T0=50,fluid0={"FoodWater": 0.5, "FoodProtein": 0.5})

nw.solve("design")
nw.print_results()

print(nw.results['Connection'])

m_FoodProtein_c1 = c1.m.val * c1.fluid.val['FoodProtein']
m_FoodProtein_c2 = c2.m.val * c2.fluid.val['FoodProtein']


print(f"\n Species flow split is {m_FoodProtein_c2/m_FoodProtein_c1}")

#print(f"\n heat flows are  {se.Q.val}")
#print(se.Qout.val)

print(f"\n")



