
from tespy.components import Sink, Source, DropletSeparator
from tespy.connections import Connection, Bus
from tespy.connections import Connection
from tespy.tools import ExergyAnalysis
from tespy.networks import Network
from tespy.tools.helpers import get_chem_ex_lib
chemexlib = get_chem_ex_lib("Ahrendts")


nw = Network(T_unit='C', p_unit='bar', h_unit='kJ / kg', iterinfo=False)

# components
so = Source('Two Phase Inflow')
si_gas = Sink('Gas Outflow')
si_liq = Sink('Liquid Outflow')
ds = DropletSeparator('Droplet Separator')

# connections
so_ds = Connection(so, 'out1', ds, 'in1', label='Two Phase Inflow')
ds_si_gas = Connection(ds, 'out2', si_gas, 'in1', label='Gas')
ds_si_liq = Connection(ds, 'out1', si_liq, 'in1', label='Liquid')

nw.add_conns(so_ds, ds_si_gas, ds_si_liq)

# define parameters
so_ds.set_attr(x=0.5, p=1.2, fluid={'Water': 1}, m=5)

# solve
nw.solve('design')
nw.print_results()


""" +++ exergy analysis +++ """
# define ambient
p_amb = 1
T_amb = 10

# define busses (no need to add them to system)
two_phase_in = Bus('Two Phase In')
two_phase_in.add_comps({'comp': so, 'base': 'bus'})

liquid_out = Bus('Liquid Out')
liquid_out.add_comps({'comp': si_liq})

gas_out = Bus('Gas Out')
gas_out.add_comps({'comp': si_gas})

# exergy analysis
exe_eco_input = {'Droplet Separator_Z': 50, 'Two Phase Inflow_c': 10}
ean = ExergyAnalysis(nw, E_P=[gas_out, liquid_out], E_F=[two_phase_in], E_L=[])
ean.analyse(pamb=p_amb, Tamb=T_amb, Chem_Ex=chemexlib)
ean.evaluate_exergoeconomics(Tamb=T_amb, Exe_Eco_Costs=exe_eco_input)
ean.print_results(Exe_Eco_An=True)
