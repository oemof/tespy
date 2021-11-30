# -*- coding: utf-8 -*-

from tespy.components import HeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.networks import Network

# %% network

nw = Network(fluids=['INCOMP::MEG[0.2]', 'INCOMP::MPG[0.4]'], T_unit='C', p_unit='bar',
             h_unit='kJ / kg', m_unit='kg / s')

# %% components

meg_in = Source('meg source')
meg_out = Sink('meg sink')
hx = HeatExchanger('exchanger')
mpg_in = Source('mpg source')
mpg_out = Sink('mpg sink')

# %% connections

# hot side
meg_in_hx = Connection(meg_in, 'out1', hx, 'in2')
hx_meg_out = Connection(hx, 'out2', meg_out, 'in1')
nw.add_conns(meg_in_hx, hx_meg_out)

# cold side
mpg_in_hx = Connection(mpg_in, 'out1', hx, 'in1')
hx_mpg_out = Connection(hx, 'out1', mpg_out, 'in1')
nw.add_conns(mpg_in_hx, hx_mpg_out)


# %% component parametrization

hx.set_attr(Q=-1e3, pr1=1, pr2=1)

# %% connection parametrization

meg_in_hx.set_attr(fluid={'MEG[0.2]': 1, 'MPG[0.4]': 0}, m=0.4, T=10, p=1)
mpg_in_hx.set_attr(fluid={'MEG[0.2]': 0, 'MPG[0.4]': 1}, m=0.2, T=20, p=1)

path = 'binary_mixtures'
nw.solve('design')
nw.print_results()
nw.save(path)
