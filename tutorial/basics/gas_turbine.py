# %%[sec_1]
from tespy.networks import Network
from tespy.components import (
    DiabaticCombustionChamber, Turbine, Source, Sink, Compressor,
    Generator, PowerBus, PowerSink
)
from tespy.connections import Connection, Ref, PowerConnection

# define full fluid list for the network"s variable space
nw = Network()
nw.units.set_defaults(temperature="degC", pressure="bar", heat="MW", power="MW")

cp = Compressor("Compressor")
cc = DiabaticCombustionChamber("combustion chamber")
tu = Turbine("turbine")
air = Source("air source")
fuel = Source("fuel source")
fg = Sink("flue gas sink")
# %%[sec_2]
c2 = Connection(air, "out1", cc, "in1", label="2")
c3 = Connection(cc, "out1", fg, "in1", label="3")
c5 = Connection(fuel, "out1", cc, "in2", label="5")
nw.add_conns(c2, c3, c5)
# %%[sec_3]
cc.set_attr(pr=1, eta=1, lamb=1.5, ti=10)

c2.set_attr(
    p=1, T=20,
    fluid={"Ar": 0.0129, "N2": 0.7553, "CO2": 0.0004, "O2": 0.2314}
)
c5.set_attr(p=1, T=20, fluid={"CO2": 0.04, "CH4": 0.96, "H2": 0})

nw.solve(mode="design")
nw.print_results()
# %%[sec_4]
cc.set_attr(ti=None)
c5.set_attr(m=1)
nw.solve(mode="design")
# %%[sec_5]
cc.set_attr(lamb=None)
c3.set_attr(T=1400)
nw.solve(mode="design")
# %%[sec_6]
c5.set_attr(fluid={"CO2": 0.03, "CH4": 0.92, "H2": 0.05})
nw.solve(mode="design")
# %%[sec_7]
print(nw.results["Connection"])
# %%[sec_8]
nw.del_conns(c2, c3)
c1 = Connection(air, "out1", cp, "in1", label="1")
c2 = Connection(cp, "out1", cc, "in1", label="2")
c3 = Connection(cc, "out1", tu, "in1", label="3")
c4 = Connection(tu, "out1", fg, "in1", label="4")
nw.add_conns(c1, c2, c3, c4)

generator = Generator("generator")
grid = PowerSink("grid")
shaft = PowerBus("shaft", num_in=1, num_out=2)

e1 = PowerConnection(tu, "power", shaft, "power_in1", label="e1")
e2 = PowerConnection(shaft, "power_out1", cp, "power", label="e2")
e3 = PowerConnection(shaft, "power_out2", generator, "power_in", label="e3")
e4 = PowerConnection(generator, "power_out", grid, "power", label="e4")

nw.add_conns(e1, e2, e3, e4)

generator.set_attr(eta=0.98)
# %%[sec_9]
cp.set_attr(eta_s=0.85, pr=15)
tu.set_attr(eta_s=0.90)
c1.set_attr(
    p=1, T=20,
    fluid={"Ar": 0.0129, "N2": 0.7553, "CO2": 0.0004, "O2": 0.2314}
)
c3.set_attr(m=30)
c4.set_attr(p=Ref(c1, 1, 0))
nw.solve("design")
c3.set_attr(m=None, T=1200)
nw.solve("design")
nw.print_results()
# %%[sec_10]
# unset the value, set Referenced value instead
c5.set_attr(p=None)
c5.set_attr(p=Ref(c2, 1.05, 0))
nw.solve("design")
# %%[sec_11]
cc.set_attr(pr=0.97, eta=0.98)
nw.set_attr(iterinfo=False)
import matplotlib.pyplot as plt
import numpy as np

# make text reasonably sized
plt.rc('font', **{'size': 18})

data = {
    'T_3': np.linspace(900, 1400, 11),
    'pr': np.linspace(10, 30, 11)
}
power = {
    'T_3': [],
    'pr': []
}
eta = {
    'T_3': [],
    'pr': []
}

for T in data['T_3']:
    c3.set_attr(T=T)
    nw.solve('design')
    power['T_3'] += [abs(e4.E.val)]
    eta['T_3'] += [abs(e4.E.val) / cc.ti.val * 100]

# reset to base value
c3.set_attr(T=1200)

for pr in data['pr']:
    cp.set_attr(pr=pr)
    nw.solve('design')
    power['pr'] += [abs(e4.E.val)]
    eta['pr'] += [abs(e4.E.val) / cc.ti.val * 100]

# reset to base value
cp.set_attr(pr=15)

fig, ax = plt.subplots(2, 2, figsize=(16, 8), sharex='col', sharey='row')

ax = ax.flatten()
[(a.grid(), a.set_axisbelow(True)) for a in ax]

i = 0
for key in data:
    ax[i].scatter(data[key], eta[key], s=100, color="#1f567d")
    ax[i + 2].scatter(data[key], power[key], s=100, color="#18a999")
    i += 1

ax[0].set_ylabel('Efficiency in %')
ax[2].set_ylabel('Power in MW')
ax[2].set_xlabel('Turbine inlet temperature °C')
ax[3].set_xlabel('Compressure pressure ratio')

plt.tight_layout()
fig.savefig('gas_turbine_parametric.svg')
plt.close()
# %%[sec_12]
c3.set_attr(T=None)

data = np.linspace(0.1, 0.2, 6)
T3 = []

for oxy in data[::-1]:
    c3.set_attr(fluid={"O2": oxy})
    nw.solve('design')
    T3 += [c3.T.val]

T3 = T3[::-1]

# reset to base value
c3.fluid.is_set.remove("O2")
c3.set_attr(T=1200)

fig, ax = plt.subplots(1, figsize=(16, 8))

ax.scatter(data * 100, T3, s=100, color="#1f567d")
ax.grid()
ax.set_axisbelow(True)

ax.set_ylabel('Turbine inlet temperature in °C')
ax.set_xlabel('Oxygen mass fraction in flue gas in %')

plt.tight_layout()
fig.savefig('gas_turbine_oxygen.svg')
plt.close()

# %%[sec_13]
# fix mass fractions of all potential fluids except combustion gases
c5.set_attr(fluid={"CO2": 0.03, "O2": 0, "H2O": 0, "Ar": 0, "N2": 0, "CH4": None, "H2": None})
c5.set_attr(fluid_balance=True)


data = np.linspace(50, 60, 11)

CH4 = []
H2 = []

for ti in data:
    cc.set_attr(ti=ti)
    nw.solve('design')
    CH4 += [c5.fluid.val["CH4"] * 100]
    H2 += [c5.fluid.val["H2"] * 100]

nw.assert_convergence()

fig, ax = plt.subplots(1, figsize=(16, 8))

ax.scatter(data, CH4, s=100, color="#1f567d", label="CH4 mass fraction")
ax.scatter(data, H2, s=100, color="#18a999", label="H2 mass fraction")
ax.grid()
ax.set_axisbelow(True)
ax.legend()

ax.set_ylabel('Mass fraction of the fuel in %')
ax.set_xlabel('Thermal input in MW')
ax.set_ybound([0, 100])

plt.tight_layout()
fig.savefig('gas_turbine_fuel_composition.svg')
plt.close()
# %%[sec_14]
