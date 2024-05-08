from tespy.networks import Network
from tespy.components import (Source, Sink, CombustionChamber)
from tespy.connections import Connection, Bus
from tespy.tools import ExergyAnalysis
from tespy.tools.helpers import get_chem_ex_lib
chemexlib = get_chem_ex_lib("Ahrendts")

testCase = 1

nw = Network(T_unit="C", p_unit="bar", h_unit="kJ / kg")

# components
chamber = CombustionChamber('Combustion Chamber')
air = Source('Air Inlet')
fuelgas = Source('Fuel Gas')
exhaust = Sink('Exhaust')

# connections
c1 = Connection(air, 'out1', chamber, 'in1', label='Air to Chamber')
c2 = Connection(fuelgas, 'out1', chamber, 'in2', label='Fuel to Chamber')
c3 = Connection(chamber, 'out1', exhaust, 'in1', label='to Exhaust')

nw.add_conns(c1, c2, c3)

# define parameters
match testCase:
    case 1: # lambda > 1
        chamber.set_attr(ti=10e6, lamb=1.5)
        c1.set_attr(
            p=1.0, T=20,
            fluid={"Ar": 0.0129, "N2": 0.7553, "CO2": 0.0004, "O2": 0.2314})
        c2.set_attr(
            T=20,
            fluid={"CO2": 0.04, "CH4": 0.96})
    case 2: # lambda = 1
        # outlet temperature too high...
        chamber.set_attr(ti=10e1, lamb=1)
        c1.set_attr(
            p=1.0, T=20,
            fluid={"Ar": 0.0129, "N2": 0.7553, "CO2": 0.0004, "O2": 0.2314})
        c2.set_attr(
            T=20,
            fluid={"H2": 1})

# solve
nw.solve(mode='design')
nw.print_results()


""" +++ exergy analysis +++ """
# define ambient
p_amb = 1.0
T_amb = 25

# define busses (no need to add them to system)
heat_in = Bus('Heat In')
heat_in.add_comps({'comp': fuelgas, 'base': 'bus'})

air_in = Bus('Air In')
air_in.add_comps({'comp': air, 'base': 'bus'})

exhaust_out = Bus('Exhaust Gas')
exhaust_out.add_comps({'comp': exhaust})

# exergy analysis
exe_eco_input = {'Combustion Chamber_Z': 500, 'Air Inlet_c': 2, 'Fuel Gas_c': 5}
ean = ExergyAnalysis(nw, E_P=[exhaust_out], E_F=[heat_in, air_in], E_L=[])
ean.analyse(pamb=p_amb, Tamb=T_amb, Chem_Ex=chemexlib)
ean.evaluate_exergoeconomics(Tamb=T_amb, Exe_Eco_Costs=exe_eco_input)
ean.print_results(Exe_Eco_An=True)
