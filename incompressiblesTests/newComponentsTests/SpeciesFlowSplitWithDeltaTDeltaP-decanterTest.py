import logging

from tespy.components import HeatExchangerSimple, Source, Sink, Merge, Separator 
from tespy.tools import ComponentProperties
from tespy.connections import Connection
from tespy.networks import Network
import numpy as np

from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp

from tespy.components.newcomponents import \
    DiabaticSimpleHeatExchanger,MergeDeltaP,SeparatorWithSpeciesSplits, \
        SeparatorWithSpeciesSplitsDeltaT, SeparatorWithSpeciesSplitsDeltaTDeltaP

logging.basicConfig(level=logging.DEBUG)

# %%

# caution, must write "Water" (capital W) in INCOMP backend -> CoolProp bug? Intentional?
fluids = ["INCOMP::FoodWater", "INCOMP::FoodFat", "INCOMP::FoodProtein"]
nw = Network(m_unit='kg / s', p_unit='bar', T_unit='C',h_unit='kJ / kg', h_range=[-1e2,4e3], iterinfo=True)

so = Source("Source")
se = SeparatorWithSpeciesSplitsDeltaTDeltaP("Separator",num_out=2)
si1 = Sink("Sink 1")
si2 = Sink("Sink 2")

c1 = Connection(so, "out1", se, "in1", label="1")
c2 = Connection(se, "out1", si1, "in1", label="2")
c3 = Connection(se, "out2", si2, "in1", label="3")

nw.add_conns(c1, c2, c3)

# set global guess values 
m0 = 76    # transform unit at some point [this is kt/yr]
h0 = 119        # global guess value in kJ/kg
p0 = 1        # global guess value in bar

water = 'INCOMP::FoodWater'
#water = 'Water' # wrapper gets HEOS backend

for c in nw.conns['object']:
    n_fl = 3 # len(nw.fluids)
    c.set_attr(m0=80,h0=h0,p0=p0,fluid0={'INCOMP::FoodWater': 0.92, 'INCOMP::FoodFat': 0.023, 'INCOMP::FoodProtein': 0.052})

# set some generic data for starting values
c1.set_attr(m=1, p=p0, T=95, fluid={'INCOMP::FoodWater': 0.902, "INCOMP::FoodFat": 0.023, "INCOMP::FoodProtein": 0.075}, mixing_rule="incompressible")
c2.set_attr(fluid={'INCOMP::FoodWater': 0.648, "INCOMP::FoodFat": 0.022, "INCOMP::FoodProtein": 0.33})

c3.set_attr(force_state='l')

se.set_attr(SFS={
    'val': 0.35, 'is_set': True, 
    'split_fluid' : 'INCOMP::FoodProtein', 'split_outlet' : "out1"})


# Now it is possible to set the temperatures out of the separator differently
c2.set_attr(T=90,p=p0)
c3.set_attr(T=90,p=p0)

nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
else: 
    nw.print_results()
    print(nw.results['Connection'])
    print(f"\n converged ")



import sys
sys.exit()

# deltaT
se.set_attr(deltaT=5)
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

# deltaT and different pressures
se.set_attr(deltaT=5)
c2.set_attr(p=4)
c3.set_attr(p=3)

nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
print(nw.results['Connection'])
m_T66_c1 = c1.m.val * c1.fluid.val['T66']
m_T66_c2 = c2.m.val * c2.fluid.val['T66']
print(f"\n Species flow split is {m_T66_c2/m_T66_c1}")
print(f"\n")

# deltaP
se.set_attr(deltaT=5,deltaP=1.2)
c2.set_attr(p=None)
c3.set_attr(p=None)

nw.solve("design")
if not nw.converged:
    raise Exception("not converged")
nw.print_results()
print(nw.results['Connection'])
m_T66_c1 = c1.m.val * c1.fluid.val['T66']
m_T66_c2 = c2.m.val * c2.fluid.val['T66']
print(f"\n Species flow split is {m_T66_c2/m_T66_c1}")
print(f"\n")


