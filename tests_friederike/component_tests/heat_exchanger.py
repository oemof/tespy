from tespy.networks import Network
from tespy.components import (Source, Sink, HeatExchanger)
from tespy.connections import Connection, Bus
from tespy.tools import ExergyAnalysis
from tespy.tools.helpers import get_chem_ex_lib
chemexlib = get_chem_ex_lib("Ahrendts")

testCase = 2

nw = Network(T_unit="C", p_unit="bar", h_unit="kJ / kg")

# components
hx = HeatExchanger('Heat Exchanger')
so_h = Source("Source Hot")
si_h = Sink("Sink Hot")
so_c = Source("Source Cold")
si_c = Sink("Sink Cold")

# connections
c1 = Connection(so_h, 'out1', hx, 'in1', label='Hot In')
c2 = Connection(so_c, 'out1', hx, 'in2', label='Cold In')
c3 = Connection(hx, 'out1', si_h, 'in1', label='Hot Out')
c4 = Connection(hx, 'out2', si_c, 'in1', label='Cold Out')

nw.add_conns(c1, c2, c3, c4)

# define parameters
match testCase:
    case 1:
        hx.set_attr(pr1=1, pr2=1)
        c1.set_attr(T=110, p=1.2, fluid={'Water': 1}, m=5)
        c2.set_attr(T=25, p=1.2, fluid={'Water': 1}, m=6)
        c4.set_attr(T=30)
    case 2:
        hx.set_attr(pr1=1, pr2=1)
        c1.set_attr(T=15, p=1.2, fluid={'Water': 1}, m=5)
        c2.set_attr(T=2, p=1.2, fluid={'Water': 1}, m=6)
        c4.set_attr(T=10)
    case 3:
        hx.set_attr(pr1=1, pr2=1)
        c1.set_attr(T=35, p=1.2, fluid={'Water': 1}, m=5)
        c2.set_attr(T=15, p=1.2, fluid={'Water': 1}, m=6)
        c4.set_attr(T=30)
    case 4:
        hx.set_attr(pr1=1, pr2=1)
        c1.set_attr(T=22, p=1.2, fluid={'Water': 1}, m=2)
        c2.set_attr(T=15, p=1.2, fluid={'Water': 1}, m=6)
        c4.set_attr(T=19)
    case 5:
        hx.set_attr(pr1=1, pr2=1)
        c1.set_attr(T=30, p=1.2, fluid={'Water': 1}, m=5)
        c2.set_attr(T=10, p=1.2, fluid={'Water': 1}, m=5)
        c4.set_attr(T=15)
    case 6:
        hx.set_attr(pr1=1, pr2=1)
        c1.set_attr(T=40, p=1.2, fluid={'Water': 1}, m=5)
        c2.set_attr(T=18, p=1.2, fluid={'Water': 1}, m=5)
        c4.set_attr(T=21)


# solve
nw.solve(mode='design')
nw.print_results()


""" +++ exergy analysis +++ """
# define ambient
T_amb = 20
p_amp = 1

# define busses (no need to add them to system)
warm_stream = Bus('Warm Stream')
warm_stream.add_comps({'comp': so_h, 'base': 'bus'}, {'comp':si_h})

cold_stream = Bus('Cold Stream')
cold_stream.add_comps({'comp':so_c, 'base': 'bus'}, {'comp': si_c})

# exergy and exergoeconomic analysis
match testCase:
    case 1:
        # hot heats up cold
        exe_eco_input = {'Heat Exchanger_Z': 5e1, 'Source Hot_c': 0.05, 'Source Cold_c': 0.02}
        ean = ExergyAnalysis(nw, E_F=[warm_stream], E_P=[cold_stream], E_L=[])
    case 2:
        # cold cools hot
        exe_eco_input = {'Heat Exchanger_Z': 5e1, 'Source Hot_c': 0.02, 'Source Cold_c': 0.05}
        ean = ExergyAnalysis(nw, E_F=[cold_stream], E_P=[warm_stream], E_L=[])
    case _:
        exe_eco_input = {'Heat Exchanger_Z': 5e1, 'Source Hot_c': 0.05, 'Source Cold_c': 0.02}
        ean = ExergyAnalysis(nw, E_F=[warm_stream], E_P=[cold_stream], E_L=[])

ean.analyse(pamb=p_amp, Tamb=T_amb, Chem_Ex=chemexlib)
ean.evaluate_exergoeconomics(Tamb=T_amb, Exe_Eco_Costs=exe_eco_input)
ean.print_results(Exe_Eco_An=True)
