from tespy.networks import Network
from tespy.components import (Sink, Source, Turbine, HeatExchanger, CycleCloser, Compressor)
from tespy.connections import Connection, Bus
from tespy.tools import ExergyAnalysis
from tespy.tools.helpers import get_chem_ex_lib
chemexlib = get_chem_ex_lib("Ahrendts")

# network
nw = Network(T_unit='C', p_unit='bar', h_unit='kJ / kg', m_unit='kg / s', s_unit="kJ / kgK")

# components
closer = CycleCloser('Cycle closer')
cp = Compressor('Compressor')
turb = Turbine('Turbine')
ev = HeatExchanger('Cooling heat exchanger')
co = HeatExchanger('Heat sink heat exchanger')
water_in = Source('Water source')
water_out = Sink('Water sink')
air_in = Source('Air source')
air_out = Sink('Air sink')

# connections
c0 = Connection(ev, 'out2', closer, 'in1', label='0')
c1 = Connection(closer, 'out1', cp, 'in1', label='1')
c2 = Connection(cp, 'out1', co, 'in1', label='2')
c3 = Connection(co, 'out1', turb, 'in1', label='3')
c4 = Connection(turb, 'out1', ev, 'in2', label='4')

c11 = Connection(air_in, 'out1', ev, 'in1', label='11')
c12 = Connection(ev, 'out1', air_out, 'in1', label='12')

c21 = Connection(water_in, 'out1', co, 'in2', label='21')
c22 = Connection(co, 'out2', water_out, 'in1', label='22')

nw.add_conns(c0, c1, c2, c3, c4, c11, c12, c21, c22)

# define parameters
turb.set_attr(eta_s=0.8)
cp.set_attr(eta_s=0.8)
ev.set_attr(Q=-100e3)

c0.set_attr(T=-30, p=1, fluid={"Ar": 0.0129, "N2": 0.7553, "CO2": 0.0004, "O2": 0.2314})
c2.set_attr(p=5.25)
c3.set_attr(p=5, T=35)
c4.set_attr(p=1.05)

c11.set_attr(fluid={"Ar": 0.0129, "N2": 0.7553, "CO2": 0.0004, "O2": 0.2314}, T=-10, p=1)
c12.set_attr(p=1, T=-20)

c21.set_attr(fluid={'H2O': 1}, T=25, p=1.5)
c22.set_attr(p=1.5, T=40)


# busses
power = Bus('power input')
power.add_comps(
    {'comp': turb, 'char': 1, 'base': 'component'},
    {'comp': cp, 'char': 1, 'base': 'bus'})

cool_product_bus = Bus('cooling')
cool_product_bus.add_comps(
    {'comp': air_in, 'base': 'bus'},
    {'comp': air_out})

heat_loss_bus = Bus('heat sink')
heat_loss_bus.add_comps(
    {'comp': water_in, 'base': 'bus'},
    {'comp': water_out})

nw.add_busses(power, cool_product_bus, heat_loss_bus)

# solve
nw.solve(mode='design')
nw.print_results()


""" +++ exergy analysis +++ """
# define ambient
p_amb = 1
T_amb = 25

# exergy and exergoeconomic analysis
exe_eco_input = {'Turbine_Z': 7e3, 'Compressor_Z': 5e3, 'Cooling heat exchanger_Z': 1e3,
                 'Heat sink heat exchanger_Z': 1e3, 'Water source_c': 0.001, 'Air source_c': 0.001,
                 'power input_c': 0.01}
ean = ExergyAnalysis(network=nw, E_F=[power], E_P=[cool_product_bus], E_L=[heat_loss_bus], internal_busses=[])
# power is fuel but also internal bus -> how to implement this??
ean.analyse(pamb=p_amb, Tamb=T_amb, Chem_Ex=chemexlib)
ean.evaluate_exergoeconomics(Exe_Eco_Costs=exe_eco_input, Tamb=T_amb)
ean.print_results(Exe_Eco_An=True)
"""
ean = ExergyAnalysis(network=nw, E_F=[power], E_P=[cool_product_bus], E_L=[heat_loss_bus])
ean.analyse(pamb=p_amb, Tamb=T_amb, Chem_Ex=chemexlib)  # Exe_Eco_An=True, Exe_Eco_Costs=exe_eco_input
ean.print_results()                                     # Exe_Eco_An=True
"""
