import numpy as np
from tespy.tools.fluid_properties.wrappers import FluidPropertyWrapper
from tespy.tools.global_vars import gas_constants
from MyWrapper import MyWrapper
import logging
logging.basicConfig(level=logging.DEBUG)


# coefficients   a      b       c    d        
COEF = {
   "protein": {
       "unit" : "C",
       "cp": [2008.2,     1.2089, -1.3129*1e-3,    0.0],
       "d" : [1329.9,    -0.5184,          0.0,    0.0],
    }
}

myWrapper = MyWrapper("protein", Tref=298.15, coefs=COEF)  # same as in CoolProp
h = myWrapper.h_pT(1e5, 400)
T = myWrapper.T_ph(1e5, h)

# from tespy.tools.fluid_properties import CoolPropWrapper
# coolprop_water = CoolPropWrapper("H2O")
# h_cp = coolprop_water.h_pT(1e5, 400)
# T_cp = coolprop_water.T_ph(1e5, h_cp)


from tespy.components import Sink
from tespy.components import Source
from tespy.components import SimpleHeatExchanger
from tespy.connections import Connection
from tespy.networks import Network

nwk = Network(T_unit="C", p_unit="bar", iterinfo=True)

so = Source("Source")
hx = SimpleHeatExchanger("Heatex")
si = Sink("Sink")

c1 = Connection(so, "out1", hx, "in1", label="1")
c2 = Connection(hx, "out1", si, "in1", label="2")

nwk.add_conns(c1, c2)

c1.set_attr(
    m=1, p=1, T=20,
    fluid={"protein": 1}, fluid_engines={"protein": MyWrapper}, fluid_coefs = COEF
)
c2.set_attr(p=1, T=80)

nwk.solve("design")

hx.set_attr(Q=1.5e5)
c2.set_attr(T=None)
nwk.solve("design")

nwk.print_results()

print("hey")
