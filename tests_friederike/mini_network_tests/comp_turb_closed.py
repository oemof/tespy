from tespy.networks import Network
from tespy.components import (Turbine, Compressor, CycleCloser)
from tespy.connections import Connection, Bus
from tespy.tools import ExergyAnalysis
from tespy.tools.helpers import get_chem_ex_lib
chemexlib = get_chem_ex_lib("Ahrendts")

# network
nw = Network(T_unit="C", p_unit="bar", h_unit="kJ / kg")

# components
turb = Turbine("Turbine")
compr = Compressor("Compressor")
cc = CycleCloser("Closer")

# Connections
closer_2_compr = Connection(cc, 'out1', compr, 'in1', label="Closer-Compr")
compr_2_turb = Connection(compr, 'out1', turb, 'in1', label="Compr-Turbine")
turb_2_closer = Connection(turb, 'out1', cc, 'in1', label="Turbine-Closer")

nw.add_conns(closer_2_compr, compr_2_turb, turb_2_closer)

# define parameters
turb.set_attr(eta_s=1)
compr.set_attr(pr=1.5, eta_s=1)
closer_2_compr.set_attr(fluid={'Water': 1}, m=2, p=1)
#compr_2_turb.set_attr(T=30)

# solve
nw.solve(mode='design')
nw.print_results()


""" +++ exergy analysis +++ """
# define ambient
pamb = 1
Tamb = 25

# define busses
power_out = Bus('power output')
power_out.add_comps(
    {'comp': turb, 'char': 0.9, 'base': 'component'})

power_in = Bus('power input')
power_in.add_comps(
    {'comp': compr, 'char': 0.8, 'base': 'bus'})



# can't define power input costs for compressor: need to add this part in code
exe_eco_input = {'power input_c': 10, 'Turbine_Z': 50, 'Compressor_Z': 50}
ean = ExergyAnalysis(nw, E_P=[power_out], E_F=[power_in], E_L=[], internal_busses=[])
ean.analyse(pamb=pamb, Tamb=Tamb, Chem_Ex=chemexlib)
ean.evaluate_exergoeconomics(Tamb=Tamb, Exe_Eco_Costs=exe_eco_input)
ean.print_results(Exe_Eco_An=True)
#ean.print_results()