# %%[sec_1]
from tespy.networks import Network

# create a network object with R134a as fluid
my_plant = Network()
my_plant.units.set_defaults(
    temperature="degC", pressure="bar", enthalpy="kJ/kg"
)
# %%[sec_2]
from tespy.components import (
    CycleCloser, Pump, Condenser, Turbine, SimpleHeatExchanger, Source, Sink
)

cc = CycleCloser('cycle closer')
sg = SimpleHeatExchanger('steam generator')
mc = Condenser('main condenser')
tu = Turbine('steam turbine')
fp = Pump('feed pump')

cwso = Source('cooling water source')
cwsi = Sink('cooling water sink')

from tespy.connections import Connection

c1 = Connection(cc, 'out1', tu, 'in1', label='1')
c2 = Connection(tu, 'out1', mc, 'in1', label='2')
c3 = Connection(mc, 'out1', fp, 'in1', label='3')
c4 = Connection(fp, 'out1', sg, 'in1', label='4')
c0 = Connection(sg, 'out1', cc, 'in1', label='0')

my_plant.add_conns(c1, c2, c3, c4, c0)

c11 = Connection(cwso, 'out1', mc, 'in2', label='11')
c12 = Connection(mc, 'out2', cwsi, 'in1', label='12')

my_plant.add_conns(c11, c12)
# %%[sec_3]
mc.set_attr(pr1=1, pr2=0.98)
sg.set_attr(pr=0.9)
tu.set_attr(eta_s=0.9)
fp.set_attr(eta_s=0.75)

c11.set_attr(T=20, p=1.2, fluid={'water': 1})
c12.set_attr(T=30)
c1.set_attr(T=600, p=150, m=10, fluid={'water': 1})
c2.set_attr(p=0.1)

my_plant.solve(mode='design')
my_plant.print_results()
# %%[sec_4]
mc.set_attr(ttd_u=4)
c2.set_attr(p=None)

my_plant.solve(mode='design')
my_plant.print_results()

# %%[sec_5]
# Adding feature to plot the T-s Diagram using fluprodia library
# Importing necessary library
import matplotlib.pyplot as plt
import numpy as np
from fluprodia import FluidPropertyDiagram

# Initial Setup
diagram = FluidPropertyDiagram('water')
diagram.set_unit_system(T='°C', p='bar', h='kJ/kg')

# Storing the model result in the dictionary
result_dict = {}
result_dict.update(
    {cp.label: cp.get_plotting_data()[1] for cp in my_plant.comps['object']
     if cp.get_plotting_data() is not None})

# Iterate over the results obtained from TESPy simulation
for key, data in result_dict.items():
    # Calculate individual isolines for T-s diagram
    result_dict[key]['datapoints'] = diagram.calc_individual_isoline(**data)

# Create a figure and axis for plotting T-s diagram
fig, ax = plt.subplots(1, figsize=(20, 10))
isolines = {
    'Q': np.linspace(0, 1, 2),
    'p': np.array([1, 2, 5, 10, 20, 50, 100, 300]),
    'v': np.array([]),
    'h': np.arange(500, 3501, 500)
}

# Set isolines for T-s diagram
diagram.set_isolines(**isolines)
diagram.calc_isolines()

# Draw isolines on the T-s diagram
diagram.draw_isolines(fig, ax, 'Ts', x_min=0, x_max=7500, y_min=0, y_max=650)

# Adjust the font size of the isoline labels
for text in ax.texts:
    text.set_fontsize(10)

# Plot T-s curves for each component
for key in result_dict.keys():
    datapoints = result_dict[key]['datapoints']
    _ = ax.plot(datapoints['s'], datapoints['T'], color='#ff0000', linewidth=2)
    _ = ax.scatter(datapoints['s'][0], datapoints['T'][0], color='#ff0000')

# Set labels and title for the T-s diagram
ax.set_xlabel('Entropy, s in J/kgK', fontsize=16)
ax.set_ylabel('Temperature, T in °C', fontsize=16)
ax.set_title('T-s Diagram of Rankine Cycle', fontsize=20)

# Set font size for the x-axis and y-axis ticks
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
plt.tight_layout()

# Save the T-s diagram plot as an SVG file
fig.savefig('rankine_ts_diagram.svg')

# %%[sec_6]
from tespy.components import PowerBus, PowerSink, Motor, Generator
from tespy.connections import PowerConnection

electricity = PowerBus("electricity bus", num_in=1, num_out=2)
grid = PowerSink("grid")
motor = Motor("feed water pump motor")
generator = Generator("turbine generator")

e1 = PowerConnection(tu, "power", generator, "power_in", label="e1")
e2 = PowerConnection(generator, "power_out", electricity, "power_in1", label="e2")
e3 = PowerConnection(electricity, "power_out1", motor, "power_in", label="e3")
e4 = PowerConnection(motor, "power_out", fp, "power", label="e4")
e5 = PowerConnection(electricity, "power_out2", grid, "power", label="e5")

my_plant.add_conns(e1, e2, e3, e4, e5)

generator.set_attr(eta=0.97)
motor.set_attr(eta=0.97)

my_plant.solve(mode='design')
my_plant.print_results()
# %%[sec_7]
e5.set_attr(E=10e6)
c1.set_attr(m=None)

my_plant.solve(mode='design')
my_plant.print_results()
# %%[sec_8]
my_plant.set_attr(iterinfo=False)
c1.set_attr(m=20)
e5.set_attr(E=None)

# make text reasonably sized
plt.rc('font', **{'size': 18})

data = {
    'T_livesteam': np.linspace(450, 750, 7),
    'T_cooling': np.linspace(15, 45, 7),
    'p_livesteam': np.linspace(75, 225, 7)
}
eta = {
    'T_livesteam': [],
    'T_cooling': [],
    'p_livesteam': []
}
power = {
    'T_livesteam': [],
    'T_cooling': [],
    'p_livesteam': []
}

for T in data['T_livesteam']:
    c1.set_attr(T=T)
    my_plant.solve('design')
    eta['T_livesteam'] += [e5.E.val / sg.Q.val * 100]
    power['T_livesteam'] += [e5.E.val / 1e6]

# reset to base temperature
c1.set_attr(T=600)

for T in data['T_cooling']:
    c12.set_attr(T=T)
    c11.set_attr(T=T - 10)
    my_plant.solve('design')
    eta['T_cooling'] += [e5.E.val / sg.Q.val * 100]
    power['T_cooling'] += [e5.E.val / 1e6]

# reset to base temperature
c12.set_attr(T=30)
c11.set_attr(T=20)

for p in data['p_livesteam']:
    c1.set_attr(p=p)
    my_plant.solve('design')
    eta['p_livesteam'] += [e5.E.val / sg.Q.val * 100]
    power['p_livesteam'] += [e5.E.val / 1e6]

# reset to base pressure
c1.set_attr(p=150)


fig, ax = plt.subplots(2, 3, figsize=(16, 8), sharex='col', sharey='row')

ax = ax.flatten()
[a.grid() for a in ax]

i = 0
for key in data:
    ax[i].scatter(data[key], eta[key], s=100, color="#1f567d")
    ax[i + 3].scatter(data[key], power[key], s=100, color="#18a999")
    i += 1

ax[0].set_ylabel('Efficiency in %')
ax[3].set_ylabel('Power in MW')
ax[3].set_xlabel('Live steam temperature in °C')
ax[4].set_xlabel('Feed water temperature in °C')
ax[5].set_xlabel('Live steam pressure in bar')
plt.tight_layout()
fig.savefig('rankine_parametric-darkmode.svg')
plt.close()
# %%[sec_9]
mc.set_attr(design=["ttd_u"], offdesign=["kA"])
c11.set_attr(offdesign=["v"])
c12.set_attr(design=["T"])
c1.set_attr(design=["p"])
tu.set_attr(offdesign=["cone"])
# %%[sec_10]
my_plant.solve("design")
my_plant.save("rankine_design.json")
# %%[sec_11]
partload_efficiency = []
partload_m_range = np.linspace(20, 10, 11)

for m in partload_m_range:
    c1.set_attr(m=m)
    my_plant.solve("offdesign", design_path="rankine_design.json")
    partload_efficiency += [e5.E.val / sg.Q.val * 100]


fig, ax = plt.subplots(1, figsize=(16, 8))
ax.grid()
ax.scatter(partload_m_range, partload_efficiency, s=100, color="#1f567d")
ax.set_xlabel("Mass flow in kg/s")
ax.set_ylabel("Plant electrical efficiency in %")

plt.tight_layout()
fig.savefig('rankine_partload.svg')
plt.close()
# %%[sec_12]
