from tespy.components import Sink, Source, HeatExchanger
from tespy.components.newComponents import TwoStreamHeatExchanger

from tespy.connections import Connection
from tespy.networks import Network
from tespy.tools import document_model
import shutil

import logging
logging.basicConfig(level=logging.DEBUG)


nw = Network(T_unit='C', p_unit='bar', h_unit='kJ / kg', iterinfo=True)



hot_i = Source('hot in')
hot_o = Sink('hot out')

cold_i = Source('cold in')
cold_o = Sink('cold out')

he = TwoStreamHeatExchanger('heat exchanger')
#he.component()

hot_in_he = Connection(hot_i, 'out1', he, 'in1')
hot_he_out = Connection(he, 'out1', hot_o, 'in1')
cold_in_he = Connection(cold_i, 'out1', he, 'in2')
cold_he_out = Connection(he, 'out2', cold_o, 'in1')
nw.add_conns(hot_in_he, hot_he_out, cold_in_he, cold_he_out)

for c in nw.conns['object']:
    n_fl=3
    c.set_attr(m0=0.5,h0=1e2,p0=3, mixing_rule="incompressible")
    #c.set_attr(fluid0={"HEOS::Water": 1/n_fl, 'INCOMP::PHE': 1/n_fl, 'INCOMP::S800': 1/n_fl})
    #c.set_attr(fluid0={'INCOMP::Water': 0.80,'INCOMP::PHE': 0.15,'INCOMP::S800': 0.05})
    c.set_attr(force_state='l')
    c.set_attr(good_starting_values=True)

he.set_attr(pr1=1, pr2=1)
cold_in_he.set_attr(fluid={'HEOS::Water': 1}, T=1, m=1, p=3,mixing_rule="incompressible")
hot_in_he.set_attr(fluid={'HEOS::Water': 0.8,'INCOMP::PHE': 0.15,'INCOMP::S800': 0.05}, T=35, m=0.5, p=3)
hot_he_out.set_attr(T=30)

nw.solve("design",print_results=True)
if not nw.converged:
    raise Exception("not converged")
nw.print_results(colored=True, print_results=False) 
print(nw.results['Connection'][[c for c in nw.results['Connection'].columns if not (c[-4:]=="unit" or c in ["v","Td_bp","s"])]])

# play with terminal
he.set_attr(ttd_l = 3)
hot_he_out.set_attr(T=None)
nw.solve("design",print_results=True)
if not nw.converged:
    raise Exception("not converged")
nw.print_results(colored=True, print_results=False) 
print(nw.results['Connection'][[c for c in nw.results['Connection'].columns if not (c[-4:]=="unit" or c in ["v","Td_bp","s"])]])

# play with terminal
#he.set_attr(ttd_l = None,ttd_u = 3) # this is not possible because cold stream has lower capacitance flow
he.set_attr(ttd_l = None,ttd_min = 15)    # using new model with pinch
hot_he_out.set_attr(T=None)
nw.solve("design",print_results=True)
if not nw.converged:
    raise Exception("not converged")
nw.print_results(colored=True, print_results=False) 
print(nw.results['Connection'][[c for c in nw.results['Connection'].columns if not (c[-4:]=="unit" or c in ["v","Td_bp","s"])]])

Ti1 = nw.results['Connection']['T'][0]
To1 = nw.results['Connection']['T'][1]
Ti2 = nw.results['Connection']['T'][2]
To2 = nw.results['Connection']['T'][3]

dTA = (Ti1-To2)
dTB = (To1-Ti2)
import numpy as np
LMDT  = (dTA-dTB)/np.log(dTA/dTB)
UA = -he.Q.val/LMDT

if not he.kA.val == UA:
    raise Exception("UA did not compare")

he.set_attr(ttd_l = None,ttd_min = None,kA=10e3)    # using new model with pinch
nw.solve("design",print_results=True)
if not nw.converged:
    raise Exception("not converged")
nw.print_results(colored=True, print_results=False) 
print(nw.results['Connection'][[c for c in nw.results['Connection'].columns if not (c[-4:]=="unit" or c in ["v","Td_bp","s"])]])


# switch hot and cold 
cold_in_he.set_attr(T=50, m=1, p=3,mixing_rule="incompressible")
hot_in_he.set_attr(T=20, m=0.5, p=3)
hot_he_out.set_attr(T=25)
he.set_attr(ttd_l = None,ttd_min = None,kA=None)    # using new model with pinch
nw.solve("design",print_results=True)
if not nw.converged:
    raise Exception("not converged")
nw.print_results(colored=True, print_results=False) 
print(nw.results['Connection'][[c for c in nw.results['Connection'].columns if not (c[-4:]=="unit" or c in ["v","Td_bp","s"])]])

# why negative Q and ttd_u and ttd_l

cold_in_he.set_attr(T=50, m=1, p=3,mixing_rule="incompressible")
hot_in_he.set_attr(T=20, m=0.5, p=3)
hot_he_out.set_attr(T=None)
he.set_attr(ttd_l = None,ttd_min = 5)    # using new model with pinch
nw.solve("design",print_results=True)
if not nw.converged:
    raise Exception("not converged")
nw.print_results(colored=True, print_results=False) 
print(nw.results['Connection'][[c for c in nw.results['Connection'].columns if not (c[-4:]=="unit" or c in ["v","Td_bp","s"])]])




# nw.save('tmp')

# round(ex_he.T.val - he_cw.T.val, 0)
# ex_he.set_attr(v=0.075)
# nw.solve('offdesign', design_path='tmp')
# round(he_cw.T.val, 1)
# round(he_ex.T.val, 1)
# ex_he.set_attr(v=0.1, T=40)
# nw.solve('offdesign', design_path='tmp')
# document_model(nw)
# round(he_cw.T.val, 1)
# round(he_ex.T.val, 1)
# shutil.rmtree('./tmp', ignore_errors=True)
