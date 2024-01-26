from tespy.networks import Network
from tespy.components import (Source, Sink, Merge)
from tespy.connections import Connection, Bus
from tespy.tools import ExergyAnalysis
from tespy.tools.helpers import get_chem_ex_lib
chemexlib = get_chem_ex_lib("Ahrendts")

testCase = 1            # (1: out > T0, 2: out = T0, 3: out < T0)
numInlets = 2           # 2 or 3

# network
nw = Network(T_unit="C", p_unit="bar", h_unit="kJ / kg")

# components
so_cold = Source("cold_stream")
si = Sink("out_mix")
so_hot = Source("hot_stream")
match numInlets:
    case 2:
        hx = Merge('Mixer', num_in = 2)
    case 3:
        hx = Merge('Mixer', num_in = 3)
        so_hot_2 = Source("hot_stream_2")


# connections
cold_in_mixer = Connection(so_cold, 'out1', hx, 'in1', label='cold_stream to Mixer')
hot_in_mixer = Connection(so_hot, 'out1', hx, 'in2', label='hot_stream to Mixer')
out_mix = Connection(hx, 'out1', si, 'in1', label='out as mix')

nw.add_conns(cold_in_mixer, out_mix, hot_in_mixer)

if numInlets > 2:
    hot2_in_mixer = Connection(so_hot_2, 'out1', hx, 'in3', label='hot_stream_2 to Mixer')
    nw.add_conns(hot2_in_mixer)


# define parameters
match testCase:
    case 1:
        out_mix.set_attr(T=80, p=2, fluid={'Water': 1}, m=10)
        cold_in_mixer.set_attr(T=70, fluid={'Water': 1}, m=3)
        if numInlets > 2:
            hot_in_mixer.set_attr(T=90, fluid={'Water': 1}, m=2)
    case 2:
        out_mix.set_attr(T=20, p=2, fluid={'Water': 1}, m=10)
        cold_in_mixer.set_attr(T=10, fluid={'Water': 1}, m=3)
        if numInlets > 2:
            hot_in_mixer.set_attr(T=30, fluid={'Water': 1}, m=2)
    case 3:
        out_mix.set_attr(T=15, p=2, fluid={'Water': 1}, m=10)
        cold_in_mixer.set_attr(T=5, fluid={'Water': 1}, m=3)
        if numInlets > 2:
            hot_in_mixer.set_attr(T=18, fluid={'Water': 1}, m=2)


# solve
nw.solve(mode='design')
nw.print_results()


""" +++ exergy analysis +++ """
# define ambient
T_amb = 20
p_amp = 1

# define busses
cold_in = Bus('cold')
cold_in.add_comps({'comp': so_cold, 'base': 'bus'})

hot_in = Bus('hot')
hot_in.add_comps({'comp': so_hot, 'base': 'bus'})

if numInlets > 2:
    hot2_in = Bus('hot2')
    hot2_in.add_comps({'comp': so_hot_2, 'base': 'bus'})

out_mix_B = Bus('out')
out_mix_B.add_comps({'comp': si})


match numInlets:
    case 2:
        exe_eco_input = {'Mixer_Z': 10, 'cold_stream_c': 1, 'hot_stream_c': 2}
        ean = ExergyAnalysis(nw, E_F=[cold_in, hot_in], E_P=[out_mix_B], E_L=[])
    case 3:
        exe_eco_input = {'Mixer_Z': 10, 'cold_stream_c': 1, 'hot_stream_c': 2, 'hot_stream_2_c': 3}
        ean = ExergyAnalysis(nw, E_F=[cold_in, hot_in, hot2_in], E_P=[out_mix_B], E_L=[])
ean.analyse(pamb=p_amp, Tamb=T_amb, Chem_Ex=chemexlib)
ean.evaluate_exergoeconomics(Tamb=T_amb, Exe_Eco_Costs=exe_eco_input)
ean.print_results(Exe_Eco_An=True)
