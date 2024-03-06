from tespy.networks import Network
from tespy.components import (CombustionChamber, Turbine, Source, Sink, Compressor)
from tespy.connections import Connection, Ref, Bus
from tespy.tools import ExergyAnalysis
from tespy.tools.helpers import get_chem_ex_lib
chemexlib = get_chem_ex_lib("Ahrendts")

# entsprechend Tut3 aufgebaut

nw = Network(p_unit="bar", T_unit="C")

# components
cp = Compressor("compressor")
cc = CombustionChamber("combustion chamber")
tu = Turbine("turbine")
air = Source("air source")
fuel = Source("fuel source")
fg = Sink("flue gas sink")

# connections
c3 = Connection(air, "out1", cp, "in1", label="3")
c4 = Connection(cp, "out1", cc, "in1", label="4")
c6 = Connection(cc, "out1", tu, "in1", label="6")
c7 = Connection(tu, "out1", fg, "in1", label="7")
c5 = Connection(fuel, "out1", cc, "in2", label="5")
nw.add_conns(c3, c4, c6, c7, c5)

# define parameters
c3.set_attr(
    m=151, p=4.5, T=100,
    fluid={"Ar": 0, "N2": 0.79, "CO2": 0, "O2": 0.21}
)
c4.set_attr(p=20, T=326.3)
#c6.set_attr(T=1308.87)
c7.set_attr(p=1.1, T=609.86)
c5.set_attr(m=3.83, T=25, fluid={"CH4": 1})
"""c3.set_attr(
    m=151, p=4.5, T=100,
    fluid={"Ar": 0, "N2": 0.79, "CO2": 0, "O2": 0.21}
)
c4.set_attr(p=20, T=326.3)
#c6.set_attr(T=1308.87)
c7.set_attr(p=1.1, T=609.86)
c5.set_attr(m=3.83, T=25, fluid={"CH4": 1})
"""


# bus representing generated electricity
generator = Bus("generator")
generator.add_comps(
    {"comp": tu, "char": 0.98, "base": "component"},
    {"comp": cp, "char": 0.98, "base": "bus"},
)
nw.add_busses(generator)

# solve
nw.solve(mode='design')
nw.print_results()


""" +++ exergy analysis +++ """
# define ambient
T_amb = 25
p_amb = 1.02

# define busses (no need to add them to system)
fuel_bus = Bus("fuel")
fuel_bus.add_comps({"comp": fuel, "base": "bus"})

loss_bus = Bus("loss")
loss_bus.add_comps({"comp": fg, "base": "component"})

# exergy analysis
ean = ExergyAnalysis(network=nw, E_F=[fuel_bus], E_P=[generator], E_L=[loss_bus])
ean.analyse(pamb=p_amb, Tamb=T_amb, Chem_Ex=chemexlib)
ean.print_results()
