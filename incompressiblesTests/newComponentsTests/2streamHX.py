from tespy.components import Sink, Source, HeatExchanger
from tespy.components.newcomponents import TwoStreamHeatExchanger

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


he.set_attr(pr1=1, pr2=1)
cold_in_he.set_attr(fluid={'water': 1}, T=10, m=1, p=3)
hot_in_he.set_attr(fluid={'water': 1}, T=35, m=0.5, p=3)
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
