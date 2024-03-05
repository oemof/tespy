from tespy.networks import Network
from tespy.components import (Source, Compressor, DiabaticCombustionChamber, Turbine, Sink, Valve, Merge, Splitter)
from tespy.connections import Connection, Ref, Bus
from CoolProp.CoolProp import PropsSI as PSI
from tespy.tools.helpers import mass_flow

# von Jubran

# network
fluid_list = ["Ar", "N2", "O2", "CO2", "CH4", "H2O", "H2", "ethane", "C3H8", "n-Butane", "n-Pentane", "n-Hexane"]
nw = Network(fluids=fluid_list, p_unit="bar", T_unit="C")

# components
exp = Turbine('EXP')
ac = Compressor('AC')
spl = Splitter('Split')
tv1 = Valve('TV1')
# tv2 = Valve('TV2')
# m_lambda = Merge('lambda')
m_tit = Merge('TIT')
sc = DiabaticCombustionChamber('Stoichiometric combustion')
c1 = Source('Air')
c14 = Sink('Exaust')
c10 = Source('Fuel')
sink_abgas = Sink("Sink Abgas")
sink_cooling = Sink("Sink Cooling")
sink_spliter = Sink("Sink Splitter")
sink_merge_tit = Sink("Sink Mix tit")

# connections
so_to_comp = Connection(c1, 'out1', ac, 'in1', label='1')
comp_to_split = Connection(ac, 'out1', spl, 'in1', label='2')
split_to_sc = Connection(spl, 'out1', sc, 'in1', label='3')
split_to_tv1 = Connection(spl, 'out2', tv1, 'in1', label='7')

# later
sc_to_m_tit = Connection(sc, 'out1', m_tit, 'in1', label='12')
tv1_to_m_tit = Connection(tv1, 'out1', m_tit, 'in2', label='8')
#m_tit_to_sink = Connection(m_tit, 'out1', sink_merge_tit, 'in1', label='13')

# later 2
m_tit_to_exp = Connection(m_tit, 'out1', exp, 'in1', label='13')
exp_to_exit = Connection(exp, 'out1', c14, 'in1', label='14')
fuel_to_sc = Connection(c10, 'out1', sc, 'in2', label='10')


# nw.add_conns(so_to_comp,comp_to_split, split_to_sc, split_to_tv1, sc_to_m_tit, tv1_to_m_tit, m_tit_to_exp, exp_to_exit, fuel_to_sc)
nw.add_conns(so_to_comp,comp_to_split, split_to_sc, sc_to_m_tit, split_to_tv1, tv1_to_m_tit, fuel_to_sc, m_tit_to_exp, exp_to_exit)

#
generator = Bus("generator")
generator.add_comps(
    {"comp": exp, "char": 0.995, "base": "component"},  # 99.5 % mech. ele. efficiency
    {"comp": ac, "char": 0.995, "base": "bus"},  # 99.5 % mech. ele. efficiency
)
generator.set_attr(P=-4.3e8)
#nw.add_busses(generator)
# total 400 instead of mass flow as input

# parameters
# air source
air_0= {
        "CO2": 0.00027, "ethane": 0, "N2": 0.78092,
        "C3H8": 0, "CH4": 0,  "O2": 0.20947, "H2O":0,
        "n-Butane": 0, "n-Pentane": 0, "n-Hexane": 0,
        "Ar":0.00934, "H2": 0
    }

psatH2O = PSI('P', 'Q', 0, 'T', 5+273.15, 'water') * 0.00001
rh = 75
pH2O = float(rh) / 100 * psatH2O

xH2O = PSI('M', 'water') * (pH2O / 1.013) / (
    pH2O/1.013 * PSI('M', 'water') + (1 - pH2O/1.013) * PSI('M', 'air'))

air_0 = {key: value * (1-xH2O) for key, value in air_0.items()}
air_0['H2O'] = xH2O
air_0 = mass_flow(air_0)

so_to_comp.set_attr(m=841.75,
    p=1.013, T=5,
    fluid=air_0
)

# Fuel
fuel_n = {
        "CO2": 0.01800, "ethane": 0.03600, "N2": 0.10300,
        "C3H8": 0.00600, "CH4": 0.83380,  "O2": 0.000010, "H2O":0,
        "n-Butane": 0.00200, "n-Pentane": 0.00050, "n-Hexane": 0.00069,
        "Ar":0, "H2": 0
    }
fuel = mass_flow(fuel_n)
fuel_to_sc.set_attr(p=19.25,
                    T=195,
    fluid=fuel
)


# Chamber
sc.set_attr(pr=0.95, eta=0.98, lamb=1.9)
#sc_to_sink.set_attr(T=1700)

# # TIT
# m_tit_to_exp.set_attr(T=1340)

# compressor
ac.set_attr(eta_s=0.90, pr=19)

# # expander
exp.set_attr(eta_s=0.90)
exp_to_exit.set_attr(p=1.051)

# splitter
split_to_tv1.set_attr(m=148.83)
#m_tit_to_exp.set_attr(T=1340)


# valve
#tv1.set_attr(pr=1)

# starting values for out
m_tit_to_exp.set_attr(fluid0={"O2": 0.001})

# solve
nw.solve(mode="design")
nw.print_results()
