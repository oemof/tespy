from tespy.networks import Network
from tespy.components import (Source, Sink, Turbine)
from tespy.connections import Connection, Bus
from tespy.tools import ExergyAnalysis
from tespy.tools.helpers import get_chem_ex_lib
chemexlib = get_chem_ex_lib("Ahrendts")

testCase = 3

# network
nw = Network(T_unit="C", p_unit="bar", h_unit="kJ / kg")

# components
tu = Turbine("Turbine")
so = Source("Source")
si = Sink("Sink")

# connections
c1 = Connection(so, 'out1', tu, 'in1', label="Inlet")
c2 = Connection(tu, 'out1', si, 'in1', label="Outlet")

nw.add_conns(c1,c2)

# define parameters
tu.set_attr(eta_s=1)

c1.set_attr(fluid={'Water': 1}, p=100, m=20)
c2.set_attr(p=1.2)

match testCase:
    case 1:
        c1.set_attr(T=600)
    case 2:
        c1.set_attr(T=25.1)
    case 3:
        c1.set_attr(T=24)


# solve
nw.solve(mode='design')
nw.print_results()


""" +++ exergy analysis +++ """
# define ambient
p_amb = 1
T_amb = 25

# define busses
power = Bus('power output')
power.add_comps({'comp': tu, 'char': 0.9, 'base': 'component'})

steam = Bus('fresh steam dif')
steam.add_comps({'comp': so, 'base': 'bus'}, {'comp': si})

# exergy and exergoeconomic analysis
exe_eco_input = {'Source_c': 0.02, 'Turbine_Z': 5e3}
ean = ExergyAnalysis(nw, E_F=[steam], E_P=[power], E_L=[])
ean.analyse(pamb=p_amb, Tamb=T_amb, Chem_Ex=chemexlib)
ean.evaluate_exergoeconomics(Exe_Eco_Costs=exe_eco_input, Tamb=T_amb)
ean.print_results(Exe_Eco_An=True)

# print(power.P.is_set)
# print(power.P)
# print(power.get_attr("P"))
# print(tu.E_bus)
# print(tu.E_bus["massless"])
# print(steam.comps.index[0].component())
# print(steam.comps.loc[steam.comps.index.str.contains("u")])
