# -*- coding: utf-8 -*-

from tespy.components import Compressor
from tespy.components import Condenser
from tespy.components import CycleCloser
from tespy.components import HeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Valve
from tespy.components import Pump
from tespy.connections import Connection
from tespy.connections import Bus
from tespy.networks import Network
from tespy.tools.characteristics import CharLine
from tespy.tools.characteristics import load_default_char as ldc
from tespy.tools import ExergyAnalysis
import numpy as np
from plotly.offline import plot
import plotly.graph_objects as go
from fluprodia import FluidPropertyDiagram
import pandas as pd

# %% network
pamb = 1.013  # ambient pressure
Tamb = 2.8  # ambient temperature

# mean geothermal temperature (mean value of ground feed and return flow)
Tgeo = 9.5

nw = Network(fluids=['water', 'NH3'], T_unit='C', p_unit='bar',
             h_unit='kJ / kg', m_unit='kg / s')

# %% components

cc = CycleCloser('cycle closer')

# heat pump system
cd = Condenser('condenser')
va = Valve('valve')
ev = HeatExchanger('evaporator')
cp = Compressor('compressor')

# geothermal heat collector
gh_in = Source('ground heat feed flow')
gh_out = Sink('ground heat return flow')
ghp = Pump('ground heat loop pump')

# heating system
hs_feed = Sink('heating system feed flow')
hs_ret = Source('heating system return flow')
hsp = Pump('heating system pump')

# %% connections

# heat pump system
cc_cd = Connection(cc, 'out1', cd, 'in1')
cd_va = Connection(cd, 'out1', va, 'in1')
va_ev = Connection(va, 'out1', ev, 'in2')
ev_cp = Connection(ev, 'out2', cp, 'in1')
cp_cc = Connection(cp, 'out1', cc, 'in1')
nw.add_conns(cc_cd, cd_va, va_ev, ev_cp, cp_cc)

# geothermal heat collector
gh_in_ghp = Connection(gh_in, 'out1', ghp, 'in1')
ghp_ev = Connection(ghp, 'out1', ev, 'in1')
ev_gh_out = Connection(ev, 'out1', gh_out, 'in1')
nw.add_conns(gh_in_ghp, ghp_ev, ev_gh_out)

# heating system
hs_ret_hsp = Connection(hs_ret, 'out1', hsp, 'in1')
hsp_cd = Connection(hsp, 'out1', cd, 'in2')
cd_hs_feed = Connection(cd, 'out2', hs_feed, 'in1')
nw.add_conns(hs_ret_hsp, hsp_cd, cd_hs_feed)


# %% component parametrization

# condenser
cd.set_attr(pr1=0.99, pr2=0.99, ttd_u=5, design=['pr2', 'ttd_u'],
            offdesign=['zeta2', 'kA_char'])
# evaporator
kA_char1 = ldc('heat exchanger', 'kA_char1', 'DEFAULT', CharLine)
kA_char2 = ldc('heat exchanger', 'kA_char2', 'EVAPORATING FLUID', CharLine)
ev.set_attr(pr1=0.99, pr2=0.99, ttd_l=5,
            kA_char1=kA_char1, kA_char2=kA_char2,
            design=['pr1', 'ttd_l'], offdesign=['zeta1', 'kA_char'])
# compressor
cp.set_attr(eta_s=0.8, design=['eta_s'], offdesign=['eta_s_char'])
# heating system pump
hsp.set_attr(eta_s=0.75, design=['eta_s'], offdesign=['eta_s_char'])
# ground heat loop pump
ghp.set_attr(eta_s=0.75, design=['eta_s'], offdesign=['eta_s_char'])


# %% connection parametrization

# heat pump system
cc_cd.set_attr(fluid={'water': 0, 'NH3': 1})
ev_cp.set_attr(Td_bp=3)

# geothermal heat collector
gh_in_ghp.set_attr(T=Tgeo + 1.5, p=1.5, fluid={'water': 1, 'NH3': 0},
                   )
ev_gh_out.set_attr(T=Tgeo - 1.5, p=1.5)

# heating system
cd_hs_feed.set_attr(T=40, p=2, fluid={'water': 1, 'NH3': 0})
hs_ret_hsp.set_attr(T=35, p=2)

# starting values
ev_cp.set_attr(p0=5)
cc_cd.set_attr(p0=18)

# %% create busses

# characteristic function for motor efficiency
x = np.array([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4])
y = np.array([0, 0.86, 0.9, 0.93, 0.95, 0.96, 0.95, 0.93])

# power bus
char = CharLine(x=x, y=y)
power = Bus('power input')
power.add_comps({'comp': cp, 'char': char, 'base': 'bus'},
                {'comp': ghp, 'char': char, 'base': 'bus'},
                {'comp': hsp, 'char': char, 'base': 'bus'})

# consumer heat bus
heat_cons = Bus('heating system')
heat_cons.add_comps({'comp': hs_ret, 'base': 'bus'}, {'comp': hs_feed})

# geothermal heat bus
heat_geo = Bus('geothermal heat')
heat_geo.add_comps({'comp': gh_in, 'base': 'bus'},
                   {'comp': gh_out})


nw.add_busses(power, heat_cons, heat_geo)

# %% key paramter

cd.set_attr(Q=-4e3)

# %% design calculation

path = 'NH3'
nw.solve('design')
# alternatively use:
# nw.solve('design', init_path=path)
print("\n##### DESIGN CALCULATION #####\n")
nw.print_results()
nw.save(path)

# %% plot h_log(p) diagram

# generate plotting data
result_dict = {}
result_dict.update({ev.label: ev.get_plotting_data()[2]})
result_dict.update({cp.label: cp.get_plotting_data()[1]})
result_dict.update({cd.label: cd.get_plotting_data()[1]})
result_dict.update({va.label: va.get_plotting_data()[1]})

# create plot
diagram = FluidPropertyDiagram('NH3')
diagram.set_unit_system(T='°C', p='bar', h='kJ/kg')

for key, data in result_dict.items():
    result_dict[key]['datapoints'] = diagram.calc_individual_isoline(**data)

diagram.set_limits(x_min=0, x_max=2100, y_min=1e0, y_max=2e2)
diagram.calc_isolines()
diagram.draw_isolines('logph')

for key in result_dict.keys():
    datapoints = result_dict[key]['datapoints']
    diagram.ax.plot(datapoints['h'], datapoints['p'], color='#ff0000')
    diagram.ax.scatter(datapoints['h'][0], datapoints['p'][0], color='#ff0000')

diagram.save('NH3_logph.svg')

# %% exergy analysis

ean = ExergyAnalysis(network=nw,
                     E_F=[power, heat_geo],
                     E_P=[heat_cons])
ean.analyse(pamb, Tamb)
print("\n##### EXERGY ANALYSIS #####\n")
ean.print_results()

# create sankey diagram
links, nodes = ean.generate_plotly_sankey_input()
fig = go.Figure(go.Sankey(
    arrangement="snap",
    node={
        "label": nodes,
        'pad': 11,
        'color': 'orange'},
    link=links))
plot(fig, filename='NH3_sankey.html')

# %% plot exergy destruction

# create data for bar chart
comps = ['E_F']
E_F = ean.network_data.E_F
# top bar
E_D = [0]  # no exergy destruction in the top bar
E_P = [E_F]  # add E_F as the top bar
for comp in ean.component_data.index:
    # only plot components with exergy destruction > 1 W
    if ean.component_data.E_D[comp] > 1:
        comps.append(comp)
        E_D.append(ean.component_data.E_D[comp])
        E_F = E_F-ean.component_data.E_D[comp]
        E_P.append(E_F)
comps.append("E_P")
E_D.append(0)
E_P.append(E_F)

# create data frame and save data
df_comps = pd.DataFrame(columns=comps)
df_comps.loc["E_D"] = E_D
df_comps.loc["E_P"] = E_P
df_comps.to_csv('NH3_E_D.csv')

# %% further calculations

print("\n#### FURTHER CALCULATIONS ####\n")
# switch off iterinfo
nw.set_attr(iterinfo=False)
# offdesign test
nw.solve('offdesign', design_path=path)

# %% calculate epsilon depending on:
#    - ambient temperature Tamb
#    - mean geothermal temperature Tgeo

Tamb_design = Tamb
Tgeo_design = Tgeo
i = 0

# create data ranges and frames
Tamb_range = [1, 4, 8, 12, 16, 20]
Tgeo_range = [11.5, 10.5, 9.5, 8.5, 7.5, 6.5]
df_eps_Tamb = pd.DataFrame(columns=Tamb_range)
df_eps_Tgeo = pd.DataFrame(columns=Tgeo_range)

# calculate epsilon depending on Tamb
eps_Tamb = []
print("Varying ambient temperature:\n")
for Tamb in Tamb_range:
    i += 1
    ean.analyse(pamb, Tamb)
    eps_Tamb.append(ean.network_data.epsilon)
    print("Case %d: Tamb = %.1f °C" % (i, Tamb))

# save to data frame
df_eps_Tamb.loc[Tgeo_design] = eps_Tamb
df_eps_Tamb.to_csv('NH3_eps_Tamb.csv')

# calculate epsilon depending on Tgeo
eps_Tgeo = []
print("\nVarying mean geothermal temperature:\n")
for Tgeo in Tgeo_range:
    i += 1
    # set feed and return flow temperatures around mean value Tgeo
    gh_in_ghp.set_attr(T=Tgeo + 1.5)
    ev_gh_out.set_attr(T=Tgeo - 1.5)
    nw.solve('offdesign', init_path=path, design_path=path)
    ean.analyse(pamb, Tamb_design)
    eps_Tgeo.append(ean.network_data.epsilon)
    print("Case %d: Tgeo = %.1f °C" % (i, Tgeo))

# save to data frame
df_eps_Tgeo.loc[Tamb_design] = eps_Tgeo
df_eps_Tgeo.to_csv('NH3_eps_Tgeo.csv')


# %% calculate epsilon and COP depending on:
#     - mean geothermal temperature Tgeo
#     - heating system Temperature Ths

# create data ranges and frames
Tgeo_range = [10.5, 8.5, 6.5]
Ths_range = [42.5, 37.5, 32.5]
df_eps_Tgeo_Ths = pd.DataFrame(columns=Ths_range)
df_cop_Tgeo_Ths = pd.DataFrame(columns=Ths_range)

# calculate epsilon and COP
print("\nVarying mean geothermal temperature and "
      "heating system temperature:\n")
for Tgeo in Tgeo_range:
    # set feed and return flow temperatures around mean value Tgeo
    gh_in_ghp.set_attr(T=Tgeo+1.5)
    ev_gh_out.set_attr(T=Tgeo-1.5)
    epsilon = []
    cop = []
    for Ths in Ths_range:
        i += 1
        # set feed and return flow temperatures around mean value Ths
        cd_hs_feed.set_attr(T=Ths+2.5)
        hs_ret_hsp.set_attr(T=Ths-2.5)
        if Ths == Ths_range[0]:
            nw.solve('offdesign', init_path=path, design_path=path)
        else:
            nw.solve('offdesign', design_path=path)
        ean.analyse(pamb, Tamb_design)
        epsilon.append(ean.network_data.epsilon)
        cop += [abs(cd.Q.val) / (cp.P.val + ghp.P.val + hsp.P.val)]
        print("Case %d: Tgeo = %.1f °C, Ths = %.1f °C" % (i, Tgeo, Ths))

    # save to data frame
    df_eps_Tgeo_Ths.loc[Tgeo] = epsilon
    df_cop_Tgeo_Ths.loc[Tgeo] = cop

df_eps_Tgeo_Ths.to_csv('NH3_eps_Tgeo_Ths.csv')
df_cop_Tgeo_Ths.to_csv('NH3_cop_Tgeo_Ths.csv')


# %% calculate epsilon and COP depending on:
#     - mean geothermal temperature Tgeo
#     - heating load Q_cond

# reset heating system temperature to design values
cd_hs_feed.set_attr(T=40)
hs_ret_hsp.set_attr(T=35)

# create data ranges and frames
Tgeo_range = [10.5, 8.5, 6.5]
Q_range = np.array([4.3e3, 4e3, 3.7e3, 3.4e3, 3.1e3, 2.8e3])
df_cop_Tgeo_Q = pd.DataFrame(columns=Q_range)
df_eps_Tgeo_Q = pd.DataFrame(columns=Q_range)

# calculate epsilon and COP
print("\nVarying mean geothermal temperature and "
      "heating load:\n")
for Tgeo in Tgeo_range:
    gh_in_ghp.set_attr(T=Tgeo+1.5)
    ev_gh_out.set_attr(T=Tgeo-1.5)
    cop = []
    epsilon = []
    for Q in Q_range:
        i += 1
        cd.set_attr(Q=-Q)
        if Q == Q_range[0]:
            nw.solve('offdesign', init_path=path, design_path=path)
        else:
            nw.solve('offdesign', design_path=path)
        ean.analyse(pamb, Tamb_design)
        cop += [abs(cd.Q.val) / (cp.P.val + ghp.P.val + hsp.P.val)]
        epsilon.append(ean.network_data.epsilon)
        print("Case %s: Tgeo = %.1f °C, Q = -%.1f kW" % (i, Tgeo, Q/1000))

    # save to data frame
    df_cop_Tgeo_Q.loc[Tgeo] = cop
    df_eps_Tgeo_Q.loc[Tgeo] = epsilon

df_cop_Tgeo_Q.to_csv('NH3_cop_Tgeo_Q.csv')
df_eps_Tgeo_Q.to_csv('NH3_eps_Tgeo_Q.csv')
