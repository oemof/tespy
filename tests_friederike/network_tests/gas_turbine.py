from tespy.networks import Network
from tespy.components import (CombustionChamber, Turbine, Source, Sink, Compressor)
from tespy.connections import Connection, Ref, Bus
from tespy.tools import ExergyAnalysis
from tespy.tools.helpers import get_chem_ex_lib
chemexlib = get_chem_ex_lib("Ahrendts")

# entsprechend Tut3 aufgebaut

# Definition des Netwerks
nw = Network(T_unit='C', p_unit='bar', h_unit='kJ / kg', m_unit='kg / s')

# Definition der Komponenten
cp = Compressor("compressor")
cc = CombustionChamber("combustion chamber")
tu = Turbine("turbine")
air = Source("air source")
fuel = Source("fuel source")
fg = Sink("flue gas sink")

# Definition der Verbindungen
c1 = Connection(air, "out1", cp, "in1", label="1")
c2 = Connection(cp, "out1", cc, "in1", label="2")
c3 = Connection(cc, "out1", tu, "in1", label="3")
c4 = Connection(tu, "out1", fg, "in1", label="4")
c5 = Connection(fuel, "out1", cc, "in2", label="5")
nw.add_conns(c1, c2, c3, c4, c5)

# Definition der Parameter
c1.set_attr(
    m=151, p=4.5, T=100,
    fluid={"Ar": 0, "N2": 0.79, "CO2": 0, "O2": 0.21}
)
c2.set_attr(p=20, T=326.3)
c4.set_attr(p=1.1, T=609.86)
c5.set_attr(m=3.83, T=25, fluid={"CH4": 1})


# Bud
generator = Bus("generator")
generator.add_comps(
    {"comp": tu, "char": 0.95, "base": "component"},
    {"comp": cp, "char": 0.98, "base": "bus"},
)
nw.add_busses(generator)

# Lösen des Netzwerks
nw.solve(mode='design')
nw.print_results()


# Umgebung
T_amb = 25
p_amb = 1.02

# Busse, die für Exergieanalyse benötigt werden
fuel_bus = Bus("fuel")
fuel_bus.add_comps({"comp": fuel, "base": "bus"}, {"comp": air, "base": "bus"})

loss_bus = Bus("loss")
loss_bus.add_comps({"comp": fg, "base": "component"})

# Exergie- und exergoökonomische Analyse
exe_eco_input = {'turbine_Z': 50, 'compressor_Z': 40, 'combustion chamber_Z': 60, 'air source_c': 0, 'fuel source_c': 10}
ean = ExergyAnalysis(network=nw, E_F=[fuel_bus], E_P=[generator], E_L=[loss_bus])
ean.analyse(pamb=p_amb, Tamb=T_amb, Chem_Ex=chemexlib)
ean.evaluate_exergoeconomics(Exe_Eco_Costs=exe_eco_input, Tamb=T_amb)
ean.print_results(Exe_Eco_An=True)
