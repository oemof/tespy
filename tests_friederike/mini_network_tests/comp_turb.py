from tespy.networks import Network
from tespy.components import (Source, Sink, Turbine, Compressor)
from tespy.connections import Connection, Bus
from tespy.tools import ExergyAnalysis
from tespy.tools.helpers import get_chem_ex_lib
chemexlib = get_chem_ex_lib("Ahrendts")

# network
nw = Network(T_unit="C", p_unit="bar", h_unit="kJ / kg")

# components
turb = Turbine("Turbine")
compr = Compressor("Compressor")
so = Source("Source")
si = Sink("Sink")

# Connections
so_2_compr = Connection(so, 'out1', compr, 'in1', label="Inlet")
compr_2_turb = Connection(compr, 'out1', turb, 'in1', label="Compr-Turbine")
turb_2_si = Connection(turb, 'out1', si, 'in1', label="Outlet")

nw.add_conns(so_2_compr, compr_2_turb, turb_2_si)

# define parameters
turb.set_attr(eta_s=0.8)
compr.set_attr(pr=2, eta_s=0.9)
so_2_compr.set_attr(fluid={'Water': 1}, T=600, p=100,  m=20)
turb_2_si.set_attr(p=1.1)

# solve
nw.solve(mode='design')
nw.print_results()


""" +++ exergy analysis +++ """
# define ambient
pamb = 1
Tamb = 25

# define busses
power = Bus('power output')
power.add_comps(
    {'comp': turb, 'char': 0.6, 'base': 'component'}, {'comp': compr, 'char': 1, 'base': 'bus'})

hot_steam = Bus('fresh steam dif')
hot_steam.add_comps(
    {'comp': so, 'base': 'bus'},
    {'comp': si})


# can't define power input costs for compressor: need to add this part in code
exe_eco_input = {'Source_c': 0.2, 'Turbine_Z': 5e3, 'Compressor_Z': 4e3}
ean = ExergyAnalysis(nw, E_P=[power], E_F=[hot_steam], E_L=[], internal_busses=[])
ean.analyse(pamb=pamb, Tamb=Tamb, Chem_Ex=chemexlib)
ean.print_results()
ean.evaluate_exergoeconomics(Tamb=Tamb, Exe_Eco_Costs=exe_eco_input)
ean.print_results(Exe_Eco_An=True)
