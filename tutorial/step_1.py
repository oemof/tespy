# -*- coding: utf-8 -*-

from tespy.components import Condenser
from tespy.components import CycleCloser
from tespy.components import HeatExchangerSimple
from tespy.components import Pump
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.networks import Network

# %% network

nw = Network(fluids=['water', 'NH3'], T_unit='C', p_unit='bar',
             h_unit='kJ / kg', m_unit='kg / s')

# %% components

# sources & sinks

c_in = Source('coolant in')
cons_closer = CycleCloser('consumer cycle closer')

va = Sink('valve')

# consumer system

cd = Condenser('condenser')
rp = Pump('recirculation pump')
cons = HeatExchangerSimple('consumer')

# %% connections

# consumer system

c_in_cd = Connection(c_in, 'out1', cd, 'in1')

close_rp = Connection(cons_closer, 'out1', rp, 'in1')
rp_cd = Connection(rp, 'out1', cd, 'in2')
cd_cons = Connection(cd, 'out2', cons, 'in1')
cons_close = Connection(cons, 'out1', cons_closer, 'in1')

nw.add_conns(c_in_cd, close_rp, rp_cd, cd_cons, cons_close)

# connection condenser - evaporator system

cd_va = Connection(cd, 'out1', va, 'in1')

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
