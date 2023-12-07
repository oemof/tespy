
from tespy.components import Sink, Source, HeatExchangerSimple, Splitter
from tespy.connections import Connection
from tespy.networks import Network
import shutil
import numpy as np
import matplotlib.pyplot as plt

from tespy.components import Separator,Merge,CycleCloser,Valve
from tespy.components.newComponents import DiabaticSimpleHeatExchanger,MergeDeltaP,SeparatorWithSpeciesSplits,SplitterWithFlowSplitter,SeparatorWithSpeciesSplitsDeltaT,SeparatorWithSpeciesSplitsDeltaTDeltaP

import logging
#logging.basicConfig(level=logging.DEBUG)

network = Network(m_unit='kg / s', p_unit='bar', T_unit='C',h_unit='kJ / kg', h_range=[-1e2,4e2], iterinfo=True)

# Objects
source              = Source('source')
boiler              = HeatExchangerSimple('boiler')
press               = SeparatorWithSpeciesSplitsDeltaTDeltaP('press', num_out=2)
#presswater          = Sink('presswater')
#presscake           = Sink('presscake')
decanter            = SeparatorWithSpeciesSplitsDeltaTDeltaP('decanter', num_out=2)
#grax                = Sink('grax')
oil                = Sink('oil')
centrifuge          = SeparatorWithSpeciesSplitsDeltaTDeltaP('centrifuge',num_out=2)
thickener        = SeparatorWithSpeciesSplitsDeltaTDeltaP('thickener',num_out=2)
vapourextract1    = Sink('vapourextract1')
#solubles          = Sink('solubles')
liquidmerge      = MergeDeltaP('liquidmerge', num_in = 3)
wetproduct     = Sink('wetproduct')
drier            = SeparatorWithSpeciesSplitsDeltaTDeltaP('drier',num_out=2)
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
m0 = 100    # transform unit at some point [this is kt/yr]
h0 = 1e2    # global guess value in kJ/kg
p0 = 5      # global guess value in bar

for c in network.conns['object']:
    # n_fl = len(network.fluids)
    c.set_attr(m0=m0,h0=h0,p0=p0)#,fluid0={'Water': 1/n_fl, 'INCOMP::S800': 1/n_fl, 'INCOMP::PHE': 1/n_fl})
    c.set_attr(p=p0)

# set conditions around boiler
c1.set_attr(fluid={'Water': 0.8,'INCOMP::PHE': 0.15,'INCOMP::S800': 0.05}, m=m0, T=5, mixing_rule="incompressible")
c2.set_attr(T=95)

# set conditions around press
press.set_attr(SFS={
    'val': 0.7, 'is_set': True,
    'split_fluid' : 'PHE', 'split_outlet' : "out1"})
c3.set_attr(fluid={'Water': 0.50, 'INCOMP::S800': 0.05, 'INCOMP::PHE': 0.45})
c3.set_attr(T=85)
c4a.set_attr(T=85)
c4b.set_attr(T=95)
#c4b.set_attr(p0=1)

# set conditions around decanter
decanter.set_attr(SFS={
    'val': 0.3, 'is_set': True,
    'split_fluid' : 'PHE', 'split_outlet' : "out1"})
c5.set_attr(fluid={'Water': 0.60, 'INCOMP::S800': 0.05, 'INCOMP::PHE': 0.35})
c5.set_attr(T=90)
c6.set_attr(T=90)

# set conditions around centrifuge
centrifuge.set_attr(SFS={
    'val': 0.8, 'is_set': True,
    'split_fluid' : 'S800', 'split_outlet' : "out2"})
c8.set_attr(fluid={'Water': 0, 'INCOMP::S800': 0.99, 'INCOMP::PHE': 0.01})
c7.set_attr(T=45)
c8.set_attr(T=80)

# set conditions around thickener
c10.set_attr(fluid={'Water': 1, 'INCOMP::S800': 0, 'INCOMP::PHE': 0})
c9.set_attr(fluid={'INCOMP::PHE': 0.25})
c10.set_attr(T=105)
c9.set_attr(T=105)

c10.set_attr(p=p0)

# set conditions around liquidMerge
#c11.set_attr(p=p0)

# set conditions around drier
c12.set_attr(fluid={'Water': 0.1})
c13.set_attr(fluid={'Water': 1, 'INCOMP::S800': 0, 'INCOMP::PHE': 0})
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
    if isinstance(o,SeparatorWithSpeciesSplitsDeltaTDeltaP):
        print(f"heat exchange for {o.label} = {o.Q.val}")
print(f"\n")

for o in network.comps['object']:
    if isinstance(o,SeparatorWithSpeciesSplitsDeltaTDeltaP):
        print(f"Total heat for {o.label} = {o.Q.val / (3.6*1e6)}")
print(f"\n")

print(f"Total heat for boiler is {boiler.Q.val/(3.6*1e6):.1f}")
print(f"Total heat for presswater heater is {presswaterheater.Q.val/(3.6*1e6):.1f}")



def get_paths(self):
    def generate_paths(start, current_comp_path, current_conn_path):
        # If the component is not in the connections, append the current path, we are finished 
        if start not in comp_connections:
            comps_paths.append(current_comp_path)
            conns_paths.append(current_conn_path)
            return
        # Iterate through the connected lists
        for comp_list,conn_list in zip(comp_connections[start],conn_connections[start]):
            # Recursively generate paths
            generate_paths(comp_list[-1], current_comp_path + comp_list[1:], current_conn_path + conn_list)        

    # make lists of component and connection branches
    comp_branches = [b['components'] for b in self.nw.massflow_branches]
    conn_branches = [b['connections'] for b in self.nw.massflow_branches]
    
    # make connections dict to iterate over 
    comp_connections = {}
    conn_connections = {}
    for comps,conns in zip(comp_branches,conn_branches):
        if comps[0] not in comp_connections:
            comp_connections[comps[0]] = []
            conn_connections[comps[0]] = []
        comp_connections[comps[0]].append(comps)
        conn_connections[comps[0]].append(conns)

    # get sources to start from and iterate for each uniquie path
    sources = [k for k in comp_connections.keys() if type(k).__name__ == "Source"]
    comps_paths = []
    conns_paths = []
    for start in sources:
        generate_paths(start, [start], [])

    paths = {}
    paths['connections'] = sorted(conns_paths, key=len, reverse=True)
    paths['components'] = sorted(comps_paths, key=len, reverse=True)
    
    return paths

class SELF():
    pass
    
self = SELF()
self.nw = network

paths = get_paths(self)

print(paths['components'])

longest_path_components = max(paths['components'], key=len)
longest_path_components_label = [l.label for l in longest_path_components]

path = {}
path["T"] = []
path["source"] = []
path["source_id"] = []
path["target"] = []
for conns in paths['connections']:
    T = []
    source  = []
    source_id = []
    target  = []    
    for c in conns:
        T += [c.T.val]
        source += [c.source.label] 
        source_id += [longest_path_components_label.index(c.source.label)]
        target += [c.target.label] 
    path["T"].append(T)
    path["source"].append(source)
    path["source_id"].append(source_id)
    
    path["target"].append(target)

longest_path_components_label


import matplotlib.pyplot as plt

plt.figure()
for s,T in zip(path["source_id"],path["T"]):
    #x = np.array([i for i,_ in enumerate(T)])
    #y = np.array(T)
    x = [longest_path_components_label[i] for i in s]
    plt.plot(x, T)
    
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.legend()
plt.show()


print(x)

print(x)
print(x)
print(x)

import sys
sys.exit()




sankey_nodes = []
for name,comp in self.nw.comps.iterrows():
    sankey_nodes += [name]

#missing = [value for value in sankey_nodes if not value in longest_path_components_label]
#sankey_nodes = longest_path_components_label + missing


links = {
    'source': [],
    'target': [],
    'value': [],
    #'color': []
}


for i,c in self.nw.conns.iterrows():
    links["source"] += [sankey_nodes.index(c.source.label)]
    links["target"] += [sankey_nodes.index(c.target.label)]
    links["value"]  += [c['object'].m.val_SI]
    #links["color"]  
    #print(c)

from plotly.offline import plot
import plotly.graph_objects as go
import plotly as px



fig = go.Figure(go.Sankey(
    #arrangement="snap",
    node = dict(
      pad = 100,
      thickness = 20,
      line = dict(color = "grey", width = 0.5),
      label = sankey_nodes,
      color = "blue"
    ),
    # node={
    #     "label": sankey_nodes,
    #     'pad': 10,
    #     'thickness': 10,
    #     'color': 'orange'},
    link=links))
plot(fig, filename='R410A_sankey.html')


print("hey")


# "----------------------------------------------------------------------------------"
# "----------------------------------------------------------------------------------"
# "----------------------------------------------------------------------------------"
# "----------------------------------------------------------------------------------"




sankey_nodes = []
for name,comp in self.nw.comps.iterrows():
    sankey_nodes += [name]

#missing = [value for value in sankey_nodes if not value in longest_path_components_label]
#sankey_nodes = longest_path_components_label + missing

links = {
    'source': [],
    'target': [],
    'value': [],
    'color': []
}

#colors = px.colors.DEFAULT_PLOTLY_COLORS
colors = ["rgba(242, 142, 43, 0.75)",
          "rgba(118, 183, 178, 0.75)",
          "rgba(176, 122, 161, 0.75)",
          "rgba(156, 117, 95, 0.75)",
          "rgba(237, 201, 72, 0.75)",
          "rgba(186, 176, 172, 0.75)",
          "rgba(89, 161, 79, 0.75)",
          "rgba(255, 157, 167, 0.75)",
          "rgba(78, 121, 167, 0.75)",
          "rgba(225, 87, 89, 0.75)",
          "rgba(100, 100, 100, 1.00)"]

for fluid in self.nw.all_fluids:
    
    for i,c in self.nw.conns.iterrows():
        links["source"] += [sankey_nodes.index(c.source.label)]
        links["target"] += [sankey_nodes.index(c.target.label)]
        links["value"]  += [c['object'].m.val_SI * c['object'].fluid.val[fluid]]
        links["color"]  += [colors[list(self.nw.all_fluids).index(fluid)]]
        #links["color"]  
        #print(c)

from plotly.offline import plot
import plotly.graph_objects as go

fig = go.Figure(go.Sankey(
    #arrangement="snap",
    node = dict(
      pad = 200,
      thickness = 10,
      #line = dict(color = "grey", width = 0.5),
      label = sankey_nodes,
      color = "grey"
    ),
    # node={
    #     "label": sankey_nodes,
    #     'pad': 10,
    #     'thickness': 10,
    #     'color': 'orange'},
    link=links))
plot(fig, filename='R410A_sankey.html')


print("hey")







