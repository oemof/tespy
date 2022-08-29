# [sec_1]
from tespy.networks import Network
from tespy.components import (
    DiabaticCombustionChamber, Turbine, Source, Sink, Compressor
)
from tespy.connections import Connection, Ref, Bus

# define full fluid list for the network"s variable space
fluid_list = ["Ar", "N2", "O2", "CO2", "CH4", "H2O", "H2"]
nw = Network(fluids=fluid_list, p_unit="bar", T_unit="C")

cp = Compressor("Compressor")
cc = DiabaticCombustionChamber("combustion chamber")
tu = Turbine("turbine")
air = Source("air source")
fuel = Source("fuel source")
fg = Sink("flue gas sink")
# [sec_2]
c2 = Connection(air, "out1", cc, "in1", label="2")
c3 = Connection(cc, "out1", fg, "in1", label="3")
c5 = Connection(fuel, "out1", cc, "in2", label="5")
nw.add_conns(c2, c3, c5)
# [sec_3]
cc.set_attr(pr=1, eta=1, lamb=1.5, ti=10e6)

c2.set_attr(
    p=1, T=20,
    fluid={
        "Ar": 0.0129, "N2": 0.7553, "H2O": 0,
        "CH4": 0, "CO2": 0.0004, "O2": 0.2314, "H2": 0
    }
)
c5.set_attr(
    p=1, T=20,
    fluid={
        "CO2": 0.04, "Ar": 0, "N2": 0, "O2": 0,
        "H2O": 0, "CH4": 0.96, "H2": 0
    }
)

nw.solve(mode="design")
nw.print_results()
# [sec_4]
cc.set_attr(ti=None)
c5.set_attr(m=1)
nw.solve(mode="design")
# [sec_5]
cc.set_attr(lamb=None)
c3.set_attr(T=1400)
nw.solve(mode="design")
# [sec_6]
c5.set_attr(
    fluid={
        "CO2": 0.03, "Ar": 0, "N2": 0, "O2": 0,
        "H2O": 0, "CH4": 0.92, "H2": 0.05
    }
)
nw.solve(mode="design")
# [sec_7]
print(nw.results["Connection"])
# [sec_8]
nw.del_conns(c2, c3)
c1 = Connection(air, "out1", cp, "in1", label="1")
c2 = Connection(cp, "out1", cc, "in1", label="2")
c3 = Connection(cc, "out1", tu, "in1", label="3")
c4 = Connection(tu, "out1", fg, "in1", label="4")
nw.add_conns(c1, c2, c3, c4)

generator = Bus("generator")
generator.add_comps(
    {"comp": tu, "char": 0.98, "base": "component"},
    {"comp": cp, "char": 0.98, "base": "bus"},
)
nw.add_busses(generator)
# [sec_9]
cp.set_attr(eta_s=0.85)
tu.set_attr(eta_s=0.90)
c1.set_attr(
    p=1, T=20,
    fluid={
        "Ar": 0.0129, "N2": 0.7553, "H2O": 0,
        "CH4": 0, "CO2": 0.0004, "O2": 0.2314, "H2": 0
    }
)
c2.set_attr(p=10)
c3.set_attr(T=1400)
c4.set_attr(p=Ref(c1, 1, 0))
nw.solve("design")
nw.print_results()
# [sec_10]
c5.set_attr(p=10.1)
nw.solve("design")
# [sec_11]
nw.set_attr(iterinfo=False)
import matplotlib.pyplot as plt
import numpy as np

# make text reasonably sized
plt.rc('font', **{'size': 18})

data = {
    'eta': np.linspace(0.9, 1.0, 11) * 100,
    'pr': np.linspace(0.0, 0.1, 11) * 100
}
power = {
    'eta': [],
    'pr': []
}

for eta in data['eta']:
    cc.set_attr(eta=eta / 100)
    nw.solve('design')
    power['eta'] += [abs(generator.P.val) / 1e6]

# reset to base value
cc.set_attr(eta=0.98)

for pr in data['pr']:
    cc.set_attr(pr=1 - (pr / 100))
    nw.solve('design')
    power['pr'] += [abs(generator.P.val) / 1e6]

# reset to base value
cc.set_attr(pr=0.96)

fig, ax = plt.subplots(1, 2, figsize=(16, 8), sharey=True)

[a.grid() for a in ax]

i = 0
for key in data:
    ax[i].scatter(data[key], power[key], s=100, color="#1f567d")
    i += 1

ax[0].set_ylabel('Power in MW')
ax[0].set_xlabel('Combustion chamber efficiency in %')
ax[1].set_xlabel('Combustion chamber pressure loss in %')
plt.tight_layout()
fig.savefig('gas_turbine_power.svg')
plt.close()
# [sec_12]