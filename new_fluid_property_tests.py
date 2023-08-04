from tespy.tools.fluid_properties import s_mix_ph, s_mix_pT, T_mix_ph, h_mix_pT, v_mix_ph, h_mix_pQ, T_sat_p, CoolPropWrapper, Q_mix_ph
from collections import OrderedDict
from tespy.tools.helpers import convert_from_SI
from tespy.tools.global_vars import ERR

from tespy.tools.fluid_properties.helpers import get_number_of_fluids
from tespy.tools.data_containers import FluidProperties as dc_prop, FluidComposition as dc_flu, SimpleDataContainer as dc_simple
from tespy.tools.fluid_properties.wrappers import FluidPropertyWrapper
from tespy.tools.fluid_properties.functions import dT_mix_pdh, dT_mix_dph, dv_mix_dph, dv_mix_pdh, dh_mix_dpQ, dT_sat_dp, isentropic
from tespy.tools.global_vars import fluid_property_data as fpd
from tespy.tools.logger import logger
from copy import deepcopy
import numpy as np


# coefficients    H+        S+     a        b       c        d       M
c = {'C':      [-2.101,  -6.540,  0.109,  38.940, -0.146, -17.385, 12.011],
    'S':      [-5.242, -59.014, 14.795,  24.075,  0.071,   0.0,   32.06],
    #  'N2':     [-7.069,  51.539, 24.229,  10.521,  0.180,  -2.315, 28.0134],
    'N2':     [-9.982,  16.203, 30.418,  2.544,  -0.238,  0.0, 28.0134],
    'O2':     [-9.589,  36.116, 29.154,   6.477, -0.184,  -1.017, 31.9988],
    'H2':     [-7.823, -22.966, 26.882,   3.586,  0.105,   0.0,    2.0158],
    'CO':   [-120.809,  18.937, 30.962,   2.439, -0.280,   0.0,   28.0104],
    'CO2':  [-413.886, -87.078, 51.128,   4.368, -1.469,   0.0,   44.0098],
    'H2O': [-253.871, -11.750, 34.376,   7.841, -0.423,   0.0,   18.0152],
    'H2Ol': [-289.932, -67.147, 20.355, 109.198,  2.033,   0.0,   18.0152],
    'CH4':   [-81.242,  96.731, 11.933,  77.647,  0.142, -18.414, 16.0426],
    'SO2':  [-315.422, -43.725, 49.936,   4.766, -1.046,   0.0,   64.0588],
    'H2S':   [-32.887,   1.142, 34.911,  10.686, -0.448,   0.0,   34.0758],
    'NH3':   [-60.244, -29.402, 37.321,  18.661, -0.649,   0.0,   17.0304]}


class KKHWrapper(FluidPropertyWrapper):

    def __init__(self, fluid, reference_temperature=298.15) -> None:
        self.fluid = fluid

        if fluid not in c:
            msg = "Fluid not available in KKH database"
            raise KeyError(msg)

        self.coefficients = c[fluid]
        self.h_ref = self._h_pT(None, reference_temperature)
        self.s_ref = self._s_pT(None, reference_temperature)

    def h_pT(self, p, T):
        return self._h_pT(p, T) - self.h_ref

    def _h_pT(self, p, T):
        y = T * 1e-3
        return 1e6 * (
            self.coefficients[0]
            + self.coefficients[2] * y
            + self.coefficients[3] / 2 * y ** 2
            - self.coefficients[4] / y
            + self.coefficients[5] / 3 * y ** 3
        ) / self.coefficients[6]

    def T_ph(self, p, h):
        # return newton(self.h_pT, self._T_ph_inverse_deriv, [], h)
        return _newton(self.h_pT, self._T_ph_inverse_deriv, h)

    def _T_ph_inverse_deriv(self, p, T):
        d = 0.001
        return (
            self.h_pT(p, T + d)
            - self.h_pT(p, T - d)
        ) / (2 * d)

    def s_pT(self, p, T):
        return self._s_pT(p, T) - self.s_ref

    def _s_pT(self, p, T):
        y = T * 1e-3
        return 1e3 * (
            self.coefficients[1]
            + self.coefficients[2] * np.log(T)
            + self.coefficients[3] * y
            - self.coefficients[4] / 2 * y ** -2
            + self.coefficients[5] / 2 * y ** 2
        ) / self.coefficients[6]

    def s_ph(self, p, h):
        return _newton(self.s_pT, self._T_ph_inverse_deriv, h)


def _newton(func, deriv, y, **kwargs):
    # default valaues
    x = kwargs.get('val0', 300)
    valmin = kwargs.get('valmin', 70)
    valmax = kwargs.get('valmax', 3000)
    max_iter = kwargs.get('max_iter', 10)
    tol_rel = kwargs.get('tol_rel', ERR)
    tol_abs = kwargs.get('tol_abs', ERR)
    tol_mode = kwargs.get('tol_mode', 'abs')

    # start newton loop
    expr = True
    i = 0
    while expr:
        # calculate function residual and new value
        res = y - func(None, x)
        x += res / deriv(None, x)

        # check for value ranges
        if x < valmin:
            x = valmin
        if x > valmax:
            x = valmax
        i += 1

        if i > max_iter:
            break
        if tol_mode == 'abs':
            expr = abs(res) >= tol_abs
        elif tol_mode == 'rel':
            expr = abs(res / y) >= tol_rel
        else:
            expr = abs(res / y) >= tol_rel or abs(res) >= tol_abs

    return x


# objects = []
# for i in range(1000):
#     fluid_data = OrderedDict()
#     fluid_data["water"] = {"wrapper": KKHWrapper("H2O"), "mass_fraction": 1.0}
#     fluid_data["air"] = {"wrapper": CoolPropWrapper("air", "HEOS"), "mass_fraction": 0.0}
#     objects += [fluid_data]

# objects[0]

# print(
#     h_mix_pT(1e5, 493.15, objects[0], "ideal-cond")
#     - h_mix_pT(1e5, 298.15, objects[0], "ideal-cond")
# )

# print(
#     T_mix_ph(1e5, 3.745e5, objects[0], "ideal-cond")
# )

# print(
#     s_mix_ph(1e5, 3.745e5, objects[0], "ideal-cond")
# )

# print(
#     s_mix_pT(1e5, 493.15, objects[0], "ideal-cond")
#     - s_mix_pT(1e5, 393.15, objects[0], "ideal-cond")
# )


# fluid_data = OrderedDict()
# fluid_data["water"] = {"wrapper": CoolPropWrapper("water", "HEOS"), "mass_fraction": 0.15}
# fluid_data["O2"] = {"wrapper": CoolPropWrapper("O2", "HEOS"), "mass_fraction": 0.05}
# fluid_data["N2"] = {"wrapper": CoolPropWrapper("N2", "HEOS"), "mass_fraction": 0.70}
# fluid_data["CO2"] = {"wrapper": CoolPropWrapper("CO2", "HEOS"), "mass_fraction": 0.1}


# print(
#     h_mix_pT(1e5, 293.15, fluid_data, "ideal")
# )

# print(
#     s_mix_pT(1e5, 293.15, fluid_data, "ideal")
# )

# h = h_mix_pT(1e5, 293.15, fluid_data, "ideal-cond")


# print(
#     s_mix_ph(1e5, h, fluid_data, "ideal-cond")
# )
# print(
#     s_mix_pT(1e5, 293.15, fluid_data, "ideal-cond")
# )
# s

from tespy.connections import Connection as NewConnection, Ref

# def conceptual_test():

from tespy.networks import Network
from tespy.components import Source, Sink, Pipe, Splitter, Turbine


nwk = Network(T_unit="C", p_unit="bar", h_unit='kJ / kg')

a = Source("source")
b = Sink("sink")
h = Turbine("turbine")

c1 = NewConnection(a, "out1", h, "in1", label="1")
c2 = NewConnection(h, "out1", b, "in1", label="2")

c = Source("source2")
d = Pipe("pipe")
e = Splitter("splitter")
f = Sink("sink2")
g = Sink("sink3")

c11 = NewConnection(c, "out1", d, "in1", label="11")
c12 = NewConnection(d, "out1", e, "in1", label="12")
c13 = NewConnection(e, "out1", f, "in1", label="13")
c14 = NewConnection(e, "out2", g, "in1", label="14")

# monkey patch new isentropic method
# import tespy.components.turbomachinery.turbine as t
# t.isentropic = isentropic


nwk.add_conns(c1, c2, c11, c12, c13, c14)

c1.set_attr(p=3, Td_bp=50, m=6, fluid={"Water": 1})
c2.set_attr(h=Ref(c1, 0.9, 0), fluid={"H2": 0})
h.set_attr(eta_s=0.9)

c11.set_attr(p=3, Td_bp=3, m=4)
d.set_attr(pr=1, Q="var")
c12.set_attr(Td_bp=5)
c13.set_attr(m=0.5, fluid={"Water": 1})

# nwk.solve("design")

# nwk.print_results()

# c11.set_attr(p=None, T=136)
# nwk.solve("design")
# nwk.print_results()


# conceptual_test()


from tespy.networks import Network
from tespy.components import Source, Sink, Merge, Splitter, HeatExchangerSimple, Compressor
# from tespy.connections import Connection

# network = Network(fluids=['R134a', 'Water'])

# network.set_attr(T_unit='C', p_unit='bar', h_unit='kJ / kg')

source = Source('source123')
source2 = Source('source1234')
merge = Merge('merge123')
component1 = HeatExchangerSimple('comp1', pr=1)
splitter = Splitter('splitter123')
component2 = HeatExchangerSimple('comp2')
sink = Sink('sink123')

c1 = NewConnection(source, 'out1', merge, 'in1', p=1, h=200, m=10, fluid={'R134a': 1})
# c3 = NewConnection(source2, 'out1', merge, 'in2', h=200, fluid={'Water': 1})
# c2 = NewConnection(merge, 'out1', sink, 'in1', fluid={'R134a': .2})
c2 = NewConnection(merge, 'out1', component1, 'in1')
c3 = NewConnection(component1, 'out1', splitter, 'in1', h=180)
c4 = NewConnection(splitter, 'out1', component2, 'in1', m=1)
c5 = NewConnection(component2, 'out1', merge, 'in2', h=170)
c6 = NewConnection(splitter, 'out2', sink, 'in1')
# c6 = NewConnection(splitter, 'out2', sink, 'in1')

# nwk.add_conns(c1, c2, c3)
nwk.add_conns(c1, c2, c3, c4, c5, c6)

nwk.solve('design')
nwk.print_results()

c3.set_attr(h=None)
c5.set_attr(h=None)
component1.set_attr(Q=-1000)
component2.set_attr(Q=-500)

nwk.solve('design')
# network.solve('design')
nwk.print_results()


so1 = Source("air")
so2 = Source("Other gases")
m1 = Merge("gas mixing")
p1 = Pipe("test", pr=1, Q=0)
sp1 = Splitter("Splitter")
t1 = Turbine("Turbine", pr=.1, eta_s=.8)
cp1 = Compressor("Compressor", pr=10, eta_s=.8)
si1 = Sink("Sink1")
si2 = Sink("Sink2")

c21 = NewConnection(so1, "out1", m1, "in1", label="21", fluid={"N2": 0.76, "O2": 0.23, "Ar": 0.01}, m=10, T=400, p=1, mixing_rule="ideal-cond")
c22 = NewConnection(so2, "out1", m1, "in2", label="22", fluid={"H2O": 1}, m=.5, T=400)
c23 = NewConnection(m1, "out1", p1, "in1", label="23")
c24 = NewConnection(p1, "out1", sp1, "in1", label="24")
c25 = NewConnection(sp1, "out1", t1, "in1", label="25")
c26 = NewConnection(t1, "out1", si1, "in1", label="26")
c27 = NewConnection(sp1, "out2", cp1, "in1", label="27", m=4)
c28 = NewConnection(cp1, "out1", si2, "in1", label="28")

nwk.add_conns(c21, c22, c23, c24, c25, c26, c27, c28)

from tespy.tools.helpers import TESPyNetworkError
try:
    nwk.solve("design")
except TESPyNetworkError as e:
    print(e)

for c in [c21, c22, c23]:#, c24, c25, c26, c27, c28]:
    print(c.fluid.val)