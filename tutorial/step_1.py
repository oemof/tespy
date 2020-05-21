from tespy.components import condenser
from tespy.components import cycle_closer
from tespy.components import heat_exchanger_simple
from tespy.components import pump
from tespy.components import sink
from tespy.components import source
from tespy.connections import connection
from tespy.networks import network

# %% network

nw = network(fluids=['water', 'NH3'], T_unit='C', p_unit='bar',
             h_unit='kJ / kg', m_unit='kg / s')

# %% components

# sources & sinks

c_in = source('coolant in')
cons_closer = cycle_closer('consumer cycle closer')

va = sink('valve')

# consumer system

cd = condenser('condenser')
rp = pump('recirculation pump')
cons = heat_exchanger_simple('consumer')

# %% connections

# consumer system

c_in_cd = connection(c_in, 'out1', cd, 'in1')

close_rp = connection(cons_closer, 'out1', rp, 'in1')
rp_cd = connection(rp, 'out1', cd, 'in2')
cd_cons = connection(cd, 'out2', cons, 'in1')
cons_close = connection(cons, 'out1', cons_closer, 'in1')

nw.add_conns(c_in_cd, close_rp, rp_cd, cd_cons, cons_close)

# connection condenser - evaporator system

cd_va = connection(cd, 'out1', va, 'in1')

nw.add_conns(cd_va)

# %% component parametrization

cd.set_attr(pr1=0.99, pr2=0.99, ttd_u=5, design=['pr2', 'ttd_u'],
            offdesign=['zeta2', 'kA_char'])
rp.set_attr(eta_s=0.8, design=['eta_s'], offdesign=['eta_s_char'])
cons.set_attr(pr=0.99, design=['pr'], offdesign=['zeta'])

# %% connection parametrization

c_in_cd.set_attr(T=170, fluid={'water': 0, 'NH3': 1})
close_rp.set_attr(T=60, p=10, fluid={'water': 1, 'NH3': 0})
cd_cons.set_attr(T=90)

# %% key paramter

cons.set_attr(Q=-230e3)

# %% Calculation

nw.solve('design')
nw.print_results()
nw.save('condenser')

cons.set_attr(Q=-200e3)

nw.solve('offdesign', design_path='condenser')
nw.print_results()
