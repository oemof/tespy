# -*- coding: utf-8 -*-

from tespy import con, cmp, nwk
import matplotlib.pyplot as plt

# %% network

fluids = ['water']

nw = nwk.network(fluids=fluids, p_unit='bar', T_unit='C', h_unit='kJ / kg',
                 p_range=[0.02, 110], T_range=[20, 800], h_range=[15, 5000])

# %% components

# turbine part
vessel_turb = cmp.vessel(label='vessel_turb')
turbine_hp = cmp.turbine(label='turbine_hp')
split = cmp.splitter(label='splitter1')
turbine_lp = cmp.turbine(label='turbine_lp')

# condenser and preheater
condenser = cmp.condenser(label='condenser')
preheater = cmp.condenser(label='preheater')
vessel = cmp.vessel(label='vessel1')
merge = cmp.merge(label='merge1')

# feed water
pump = cmp.pump(label='pump')
steam_generator = cmp.heat_exchanger_simple(label='steam generator')

# sources and sinks
source = cmp.source(label='source')
sink = cmp.sink(label='sink')

# for cooling water
source_cw = cmp.source(label='source_cw')
sink_cw = cmp.sink(label='sink_cw')

# %% connections

# turbine part
fs_in = con.connection(source, 'out1', vessel_turb, 'in1')
fs = con.connection(vessel_turb, 'out1', turbine_hp, 'in1')
ext = con.connection(turbine_hp, 'out1', split, 'in1')
ext_pre = con.connection(split, 'out1', preheater, 'in1')
ext_turb = con.connection(split, 'out2', turbine_lp, 'in1')
nw.add_conns(fs_in, fs, ext, ext_pre, ext_turb)

# preheater and condenser
ext_cond = con.connection(preheater, 'out1', vessel, 'in1')
cond_ws = con.connection(vessel, 'out1', merge, 'in2')
turb_ws = con.connection(turbine_lp, 'out1', merge, 'in1')
ws = con.connection(merge, 'out1', condenser, 'in1')
nw.add_conns(ext_cond, cond_ws, turb_ws, ws)

# feed water
cond = con.connection(condenser, 'out1', pump, 'in1')
fw_c = con.connection(pump, 'out1', preheater, 'in2')
fw_w = con.connection(preheater, 'out2', steam_generator, 'in1')
fs_out = con.connection(steam_generator, 'out1', sink, 'in1')
nw.add_conns(cond, fw_c, fw_w, fs_out)

# cooling water
cw_in = con.connection(source_cw, 'out1', condenser, 'in2')
cw_out = con.connection(condenser, 'out2', sink_cw, 'in1')
nw.add_conns(cw_in, cw_out)

# %% busses

# power bus
power_bus = con.bus('power')
power_bus.add_comps([turbine_hp, -1], [turbine_lp, -1], [pump, -1])

# heating bus
heat_bus = con.bus('heat')
heat_bus.add_comps([condenser, -1])

nw.add_busses(power_bus, heat_bus)

# %% parametrization of components

vessel_turb.set_attr(mode='man')
turbine_hp.set_attr(eta_s=0.9)
turbine_lp.set_attr(eta_s=0.9)

condenser.set_attr(pr1=0.95, pr2=0.95, ttd_u=12)
preheater.set_attr(pr1=0.95, pr2=0.99, ttd_u=7)
vessel.set_attr(mode='man')

pump.set_attr(eta_s=0.8, mode='man')
steam_generator.set_attr(pr=0.95, mode='man')

# %% parametrization of connections

# fresh steam properties
fs_in.set_attr(p=100, T=550, fluid={'water': 1})

# pressure after turbine inlet vessel
fs.set_attr(p=100, design=['p'])

# pressure extraction steam
ext.set_attr(p=10, design=['p'])

# staring value for warm feed water
fw_w.set_attr(h0=310)

# closing the cycle: fluid properties at sink must be identical to fluid
# properties at the source
fs_out.set_attr(p=con.ref(fs_in, 1, 0), h=con.ref(fs_in, 1, 0))

# cooling water inlet
cw_in.set_attr(T=60, p=10, fluid={'water': 1})

# setting key parameters:
# Power of the plant
power_bus.set_attr(P=5e6)
#
cw_out.set_attr(T=110)


# %% solving

mode = 'design'

nw.solve(init_file=None, mode=mode, parallel=False)
nw.save('chp_' + mode)

file = 'chp_' + mode + '_conn.csv'
mode = 'offdesign'

# representation of part loads
P_range = [5.25e6, 5e6, 4.5e6, 4e6, 3.5e6, 3e6]

# temperatures for the heating system
T_fl = [80, 90, 100, 110, 120]

P = {}
Q = {}

# iterate over temperatures
for i in T_fl:
    cw_out.set_attr(T=i)
    P[i] = []
    Q[i] = []
    # iterate over power plant power output
    for j in P_range:
        power_bus.set_attr(P=j)

        # use an initialisation file with parameters similar to next
        # calculation
        if j == P_range[0]:
            init_file = file
        else:
            init_file = 'chp_' + mode + '_conn.csv'

        nw.solve(init_file=init_file, design_file=file, mode=mode)
        nw.save('chp_' + mode)
        P[i] += [power_bus.P]
        Q[i] += [heat_bus.P]

# plotting
colors = ['#00395b', '#74adc1', '#bfbfbf', '#b54036', '#ec6707']

fig, ax = plt.subplots()
j = 0
for i in T_fl:
    plt.plot(Q[i], P[i], '.-', Color=colors[j],
             label='$T_{VL}$ = '+str(i)+' °C', markersize=15, linewidth=2)
    j += 1
ax.set_ylabel('$P$ in W')
ax.set_xlabel('$\dot{Q}$ in W')
plt.title('P-Q diagram for CHP with backpressure steam turbine')
plt.legend(loc='lower left')

ax.set_ylim([0, 6e6])
ax.set_xlim([0, 14e6])

plt.show()

fig.savefig('PQ_diagram.svg')
