from tespy import cmp, con, nwk

# %% network

nw = nwk.network(fluids=['water', 'NH3'],
                 T_unit='C', p_unit='bar', h_unit='kJ / kg', m_unit='kg / s',
                 p_range=[0.1, 100], T_range=[1, 500], h_range=[15, 5000])

# %% components

# sources & sinks

c_in = cmp.source('coolant in')
cb = cmp.source('consumer back flow')
cf = cmp.sink('consumer feed flow')

va = cmp.sink('valve')

# consumer system

cd = cmp.condenser('condenser')
rp = cmp.pump('recirculation pump')
cons = cmp.heat_exchanger_simple('consumer')

# %% connections

# consumer system

c_in_cd = con.connection(c_in, 'out1', cd, 'in1')

cb_rp = con.connection(cb, 'out1', rp, 'in1')
rp_cd = con.connection(rp, 'out1', cd, 'in2')
cd_cons = con.connection(cd, 'out2', cons, 'in1')
cons_cf = con.connection(cons, 'out1', cf, 'in1')

nw.add_conns(c_in_cd, cb_rp, rp_cd, cd_cons, cons_cf)

# connection condenser - evaporator system

cd_va = con.connection(cd, 'out1', va, 'in1')

nw.add_conns(cd_va)

# %% component parametrization

cd.set_attr(pr1=0.99, pr2=0.99, ttd_u=5, design=['pr2', 'ttd_u'], offdesign=['zeta2', 'kA'])
rp.set_attr(eta_s=0.8, design=['eta_s'], offdesign=['eta_s_char'])
cons.set_attr(pr=0.99, design=['pr'], offdesign=['zeta'])

# %% connection parametrization

c_in_cd.set_attr(T=170, fluid={'water': 0, 'NH3': 1})
cb_rp.set_attr(T=60, p=10, fluid={'water': 1, 'NH3': 0})
cd_cons.set_attr(T=90)
cons_cf.set_attr(h=con.ref(cb_rp, 1, 0), p=con.ref(cb_rp, 1, 0))

# %% key paramter

cons.set_attr(Q=-230e3)

# %% Calculation

nw.solve('design')
nw.print_results()
nw.save('condenser')

cons.set_attr(Q=-200e3)

nw.solve('offdesign', design_path='condenser')
nw.print_results()
