
from tespy.components import Sink, Source, HeatExchangerSimple, Splitter
from tespy.connections import Connection
from tespy.networks import Network
import shutil
import numpy as np
import matplotlib.pyplot as plt

from tespy.components import Separator,Merge,CycleCloser,Valve
from tespy.components.newcomponents import DiabaticSimpleHeatExchanger,MergeWithPressureLoss,SeparatorWithSpeciesSplits,SplitWithFlowSplitter,SeparatorWithSpeciesSplitsAndDeltaT,SeparatorWithSpeciesSplitsAndDeltaTAndPr

import logging
#logging.basicConfig(level=logging.DEBUG)

fluid_list = ['HEOS::Water','INCOMP::PHE','INCOMP::S800']
network = Network(fluids=fluid_list, m_unit='kg / s', p_unit='bar', T_unit='C',h_unit='kJ / kg', h_range=[-1e2,4e2], iterinfo=True)

# Objects
source              = Source('source')
boiler              = HeatExchangerSimple('boiler')
press               = SeparatorWithSpeciesSplitsAndDeltaTAndPr('press', num_out=2)
#presswater          = Sink('presswater')
#presscake           = Sink('presscake')
decanter            = SeparatorWithSpeciesSplitsAndDeltaTAndPr('decanter', num_out=2)
#grax                = Sink('grax')
oil                = Sink('oil')
centrifuge          = SeparatorWithSpeciesSplitsAndDeltaTAndPr('centrifuge',num_out=2)
thickener        = SeparatorWithSpeciesSplitsAndDeltaTAndPr('thickener',num_out=2)
vapourextract1    = Sink('vapourextract1')
#solubles          = Sink('solubles')
liquidmerge      = MergeWithPressureLoss('liquidmerge', num_in = 3)
wetproduct     = Sink('wetproduct')
drier            = SeparatorWithSpeciesSplitsAndDeltaTAndPr('drier',num_out=2)
meal             = Sink('meal')
vapourextract2   = Sink('vapourextract2')

presswaterheater   = HeatExchangerSimple('presswaterheater')

# Connections
c1 = Connection(source, 'out1', boiler, 'in1')
c2 = Connection(boiler, 'out1', press, 'in1')
c3 = Connection(press, 'out1', liquidmerge, 'in1')
#c4 = Connection(press, 'out2', decanter, 'in1')
c4a = Connection(press, 'out2', presswaterheater, 'in1')
c4b= Connection(presswaterheater, 'out1', decanter, 'in1')
c5 = Connection(decanter, 'out1', liquidmerge, 'in2')
c6 = Connection(decanter, 'out2', centrifuge, 'in1')
c7 = Connection(centrifuge, 'out1', thickener, 'in1')
c8 = Connection(centrifuge, 'out2', oil, 'in1')
c9 = Connection(thickener, 'out1', liquidmerge, 'in3')
c10 = Connection(thickener, 'out2', vapourextract1, 'in1')
c11 = Connection(liquidmerge, 'out1', drier, 'in1')
c12 = Connection(drier, 'out1', meal, 'in1')
c13 = Connection(drier, 'out2', vapourextract2, 'in1')

network.add_conns(c1,c2,c3,c4a,c4b,c5,c6,c7,c8,c9,c10,c11,c12,c13)

# set global guess values
m0 = 100    # transform unit at some point 
h0 = 1e2    # global guess value in kJ/kg
p0 = 5      # global guess value in bar

for c in network.conns['object']:
    n_fl = len(network.fluids)
    c.set_attr(m0=m0,h0=h0,p0=p0,fluid0={'Water': 1/n_fl, 'S800': 1/n_fl, 'PHE': 1/n_fl})
    c.set_attr(p=p0)

# set conditions around boiler
c1.set_attr(fluid={'Water': 0.80,'PHE': 0.15,'S800': 0.05}, m=m0, T=5)
c2.set_attr(T=95)

# set conditions around press
press.set_attr(SFS={
    'val': 0.7, 'is_set': True,
    'split_fluid' : 'PHE', 'split_outlet' : "out1"})
c3.set_attr(fluid={'Water': 0.5, 'S800': 0.05, 'PHE': 0.45})
c3.set_attr(T=85)
c4a.set_attr(T=85)
c4b.set_attr(T=95)
#c4b.set_attr(p0=1)

# set conditions around decanter
decanter.set_attr(SFS={
    'val': 0.3, 'is_set': True,
    'split_fluid' : 'PHE', 'split_outlet' : "out1"})
c5.set_attr(fluid={'Water': 0.60, 'S800': 0.05, 'PHE': 0.35})
c5.set_attr(T=90)
c6.set_attr(T=90)

# set conditions around centrifuge
centrifuge.set_attr(SFS={
    'val': 0.8, 'is_set': True,
    'split_fluid' : 'S800', 'split_outlet' : "out2"})
c8.set_attr(fluid={'Water': 0, 'S800': 0.99, 'PHE': 0.01})
c7.set_attr(T=45)
c8.set_attr(T=80)

# set conditions around thickener
c10.set_attr(fluid={'Water': 1, 'S800': 0, 'PHE': 0})
c9.set_attr(fluid={'PHE': 0.25})
c10.set_attr(T=105)
c9.set_attr(T=105)

c10.set_attr(p=p0)

# set conditions around liquidMerge
#c11.set_attr(p=p0)

# set conditions around drier
c12.set_attr(fluid={'Water': 0.1})
c13.set_attr(fluid={'Water': 1, 'S800': 0, 'PHE': 0})
c12.set_attr(T=100)
c13.set_attr(T=100)

c13.set_attr(p=None,x=1)


# solve and print results
network.solve('design')

network.print_results()
print(network.results['Connection'])

oilmassflow = c8.m.val
print(f"oil mass flow is {oilmassflow}")
print(f"\n")
network.results["Connection"].to_csv(f"{__file__.replace('.py', '')}tespy070.csv")

# MJ to kwh
#
for o in network.comps['object']:
    if isinstance(o,SeparatorWithSpeciesSplitsAndDeltaTAndPr):
        print(f"heat exchange for {o.label} = {o.Q.val}")
print(f"\n")

for o in network.comps['object']:
    if isinstance(o,SeparatorWithSpeciesSplitsAndDeltaTAndPr):
        print(f"Total heat for {o.label} = {o.Q.val / (3.6*1e6)}")
print(f"\n")

print(f"Total heat for boiler is {boiler.Q.val/(3.6*1e6):.1f}")
print(f"Total heat for presswater heater is {presswaterheater.Q.val/(3.6*1e6):.1f}")

