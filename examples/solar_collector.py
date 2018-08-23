# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 10:37:36 2017

@author: witte
"""

from tespy import con, cmp, nwk
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

# %% network

fluid_list = ['H2O']
nw = nwk.network(fluids=fluid_list, p_unit='bar', T_unit='C',
                 p_range=[4, 10], T_range=[10, 200])

# %% components

# sinks & sources
back = cmp.source('to collector')
feed = cmp.sink('from collector')

# collector
coll = cmp.solar_collector(label='solar thermal collector')

# %% connections

b_c = con.connection(back, 'out1', coll, 'in1')
c_f = con.connection(coll, 'out1', feed, 'in1')

nw.add_conns(b_c, c_f)

# %% component parameters

# set combustion chamber fuel, air to stoichometric air ratio and thermal input
coll.set_attr(pr=0.99, Q=8e3, lkf_lin=1, lkf_quad=0.005, A=10, t_a=10)

# %% connection parameters

b_c.set_attr(p=5, T=20, fluid={'H2O': 1})
c_f.set_attr(p0=2, T=120)

# %% solving

# going through several parametrisation possibilities
mode = 'design'
nw.solve(mode=mode)
nw.print_results()

coll.set_attr(Q=7e3, E=9e3)
c_f.set_attr(T=np.nan)

nw.solve(mode=mode)
nw.print_results()

coll.set_attr(Q=np.nan, E=np.nan)
c_f.set_attr(T=100, m=1e-2)

nw.solve(mode=mode)
nw.print_results()

# looping over different temperature differences (assuming constant mass flow)
# and global radiation ()

c_f.set_attr(m=np.nan)
T_amb = np.linspace(0, 60, 14, dtype=float)
E_glob = np.linspace(100, 1000, 14, dtype=float)

df = pd.DataFrame(columns=(60 - T_amb))

for E in E_glob:
    eta = []
    coll.set_attr(E=E)
    for T in T_amb:
        coll.set_attr(t_a=T)
        nw.solve(mode=mode)
        eta += [coll.Q.val / (coll.E.val * coll.A.val)]
        if eta[-1] < 0:
            eta[-1] = np.nan

    df.loc[E] = eta

E, T = np.meshgrid(60 - T_amb, E_glob)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(E, T, df.as_matrix())
ax.set_xlabel('Temperaturdifferenz')
ax.set_ylabel('Globalstrahlung auf die schiefe Ebene')
ax.set_zlabel('Wirkungsgrad (nur thermische Verluste)')
plt.show()
