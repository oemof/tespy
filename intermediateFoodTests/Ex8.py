
from tespy.components import Sink, Source, HeatExchangerSimple, Splitter
from tespy.connections import Connection
from tespy.networks import Network
import shutil
import numpy as np
import matplotlib.pyplot as plt

from tespy.components import Separator,Merge,CycleCloser,Valve
from tespy.components.newcomponents import DiabaticSimpleHeatExchanger,MergeWithPressureLoss,SeparatorWithSpeciesSplits,SplitWithFlowSplitter

import logging
#logging.basicConfig(level=logging.DEBUG)

fluid_list = ['INCOMP::Water','INCOMP::PHE','INCOMP::S800']
network = Network(fluids=fluid_list, m_unit='kg / s', p_unit='bar', T_unit='C',h_unit='kJ / kg', h_range=[-1e2,4e2], iterinfo=True)

# Objects
source              = Source('source')
boiler              = HeatExchangerSimple('boiler')
press               = SeparatorWithSpeciesSplits('press', num_out=2)
#presswater          = Sink('presswater')
#presscake           = Sink('presscake')
decanter            = SeparatorWithSpeciesSplits('decanter', num_out=2)
#grax                = Sink('grax')
oil                = Sink('oil')
centrifuge          = SeparatorWithSpeciesSplits('centrifuge',num_out=2)
thickener        = SeparatorWithSpeciesSplits('thickener',num_out=2)
vapourextract1    = Sink('vapourextract1')
#solubles          = Sink('solubles')
liquidmerge      = MergeWithPressureLoss('liquidmerge', num_in = 3)
wetproduct     = Sink('wetproduct')
drier            = SeparatorWithSpeciesSplits('drier',num_out=2)
meal             = Sink('meal')
vapourextract2   = Sink('vapourextract2')

# Connections
c1 = Connection(source, 'out1', boiler, 'in1')
c2 = Connection(boiler, 'out1', press, 'in1')
c3 = Connection(press, 'out1', liquidmerge, 'in1')
c4 = Connection(press, 'out2', decanter, 'in1')
c5 = Connection(decanter, 'out1', liquidmerge, 'in2')
c6 = Connection(decanter, 'out2', centrifuge, 'in1')
c7 = Connection(centrifuge, 'out1', thickener, 'in1')
c8 = Connection(centrifuge, 'out2', oil, 'in1')
c9 = Connection(thickener, 'out1', liquidmerge, 'in3')
c10 = Connection(thickener, 'out2', vapourextract1, 'in1')
c11 = Connection(liquidmerge, 'out1', drier, 'in1')
c12 = Connection(drier, 'out1', meal, 'in1')
c13 = Connection(drier, 'out2', vapourextract2, 'in1')

network.add_conns(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13)

# set global guess values
m0 = 100      # transform unit at some point
h0 = 1e2      # global guess value in kJ/kg
p0 = 2        # global guess value in bar

for c in network.conns['object']:
    c.set_attr(m0=m0,h0=h0,p0=p0,fluid0={'Water': 1/3, 'PHE': 1/3, 'S800': 1/3})

# set conditions around boiler
c1.set_attr(fluid={'Water': 0.80,'PHE': 0.15,'S800': 0.05}, m=m0, h=h0, p=p0)
c2.set_attr(h=h0,p=p0)

# set conditions around press
press.set_attr(SFS={
    'val': 0.7, 'is_set': True,
    'split_fluid' : 'S800', 'split_outlet' : "out1"})
c3.set_attr(fluid={'Water': 0.50, 'PHE': 0.05, 'S800': 0.45})

# set conditions around decanter
decanter.set_attr(SFS={
    'val': 0.3, 'is_set': True,
    'split_fluid' : 'S800', 'split_outlet' : "out1"})
c5.set_attr(fluid={'Water': 0.60, 'PHE': 0.05, 'S800': 0.35})

# set conditions around centrifuge
centrifuge.set_attr(SFS={
    'val': 0.8, 'is_set': True,
    'split_fluid' : 'PHE', 'split_outlet' : "out2"})
c8.set_attr(fluid={'Water': 0,'PHE': 0.99, 'S800': 0.01})

# set conditions around thickener
c10.set_attr(fluid={'Water': 1, 'S800': 0, 'PHE': 0})
c9.set_attr(fluid={'PHE': 0.25})

# set conditions around liquidMerge
c11.set_attr(p=p0)

# set conditions around drier
c12.set_attr(fluid={'Water': 0.1})
c13.set_attr(fluid={'Water': 1, 'S800': 0, 'PHE': 0})

# solve and print results
network.solve('design')

network.print_results()
print(network.results['Connection'])

oilmassflow = c8.m.val
print(f"oil mass flow is {oilmassflow}")
network.results["Connection"].to_csv(f"{__file__.replace('.py', '')}tespy070_1.csv")

# import sys
# sys.exit()

# %% oil cleansing

fluid_list2 = ['INCOMP::Water','INCOMP::S800']
network2 = Network(fluids=fluid_list2, m_unit='kg / s', p_unit='bar', T_unit='C',h_unit='kJ / kg', h_range=[-1e2,4e2], iterinfo=True)

# Objects
sourceFat           = Source('Fat')
sourceCitricAcid    = Source('CitricAcid')
centimerge          = MergeWithPressureLoss('centimerge', num_in=2)
T2                  = SplitWithFlowSplitter('T2', num_out=2)
Oil1                = Sink('Oil1')
stripper            = SplitWithFlowSplitter('stripper', num_out=2)
Oil2                = Sink('Oil2')
scrubOil            = Sink('scrubOil')

# Connections
c1 = Connection(sourceFat, 'out1', centimerge, 'in1')
c2 = Connection(sourceCitricAcid, 'out1', centimerge, 'in2')
c3 = Connection(centimerge, 'out1', T2, 'in1')
c4 = Connection(T2, 'out1', stripper, 'in1')
c5 = Connection(T2, 'out2', Oil1, 'in1')
c6 = Connection(stripper, 'out1', scrubOil, 'in1')
c7 = Connection(stripper, 'out2', Oil2, 'in1')

network2.add_conns(c1,c2,c3,c4,c5,c6,c7)

# set global guess values
m0 = oilmassflow    # transform unit at some point
h0 = 1e2            # global guess value in kJ/kg
p0 = 2              # global guess value in bar

for c in network2.conns['object']:
    c.set_attr(m0=m0,h0=h0,p0=p0)#,fluid0={'INCOMP::Water': 1/n_fl, 'INCOMP::S800': 1/n_fl})

# set conditions around merge
c1.set_attr(fluid={'Water': 0,'S800': 1}, m=m0, h=h0, p=p0)
c2.set_attr(fluid={'Water': 1,'S800': 0}, m=0.7*m0/1000, h=h0, p=p0)
c3.set_attr(p=p0)

T2.set_attr(FS={'val': 0.75, 'is_set': True, 'split_outlet' : "out1"})
stripper.set_attr(FS={'val': 0.05, 'is_set': True, 'split_outlet' : "out1"})

# solve and print results
network2.solve('design')
network2.print_results()
print(network2.results['Connection'])

network.results["Connection"].to_csv(f"{__file__.replace('.py', '')}tespy070_2.csv")
