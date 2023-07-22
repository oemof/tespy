from fluid_properties import s_mix_ph, s_mix_pT, T_mix_ph, h_mix_pT, v_mix_ph, h_mix_pQ, T_sat_p, CoolPropWrapper
from collections import OrderedDict
from tespy.tools.helpers import newton
from tespy.tools.global_vars import ERR

from fluid_properties.helpers import newton_with_kwargs
from fluid_properties.wrappers import FluidPropertyWrapper
from fluid_properties.functions import dT_mix_pdh, dT_mix_dph, dv_mix_dph, dv_mix_pdh, dh_mix_dpQ, dT_sat_dp

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


fluid_data = OrderedDict()
fluid_data["water"] = {"wrapper": CoolPropWrapper("water", "HEOS"), "mass_fraction": 1.0}


# print(
#     h_mix_pT(1e5, 493.15, fluid_data, "ideal-cond")
#     - h_mix_pT(1e5, 393.15, fluid_data, "ideal-cond")
# )

# print(
#     s_mix_pT(1e5, 493.15, fluid_data, "ideal-cond")
#     - s_mix_pT(1e5, 393.15, fluid_data, "ideal-cond")
# )


from tespy.connections import Connection

class NewConnection(Connection):

    def __init__(self, comp1, outlet_id, comp2, inlet_id, label=None, **kwargs):
        super().__init__(comp1, outlet_id, comp2, inlet_id, label, **kwargs)

    def set_attr(self, **kwargs):
        super().set_attr(**kwargs)
        if "fluid_variable" in kwargs:
            self.fluid.is_var = kwargs.get("fluid_variable", False)
        if "fluid" in kwargs:
            for fluid in kwargs["fluid"]:
                engine = CoolPropWrapper
                back_end = "HEOS"
                self._create_fluid_wrapper(fluid, engine, back_end)

    def _create_fluid_wrapper(self, fluid, engine, back_end):
        self.fluid.back_end[fluid] = back_end
        self.fluid.engine[fluid] = engine
        self.fluid.wrapper[fluid] = engine(fluid, back_end)

    def preprocess(self):
        self.parameters = self.get_parameters().copy()

        self.var_pos = {}

        self.num_eq = 0
        for parameter in self.parameters:
            if self.get_attr(parameter).val_set:
                self.num_eq += 1

        self.num_vars = 0
        for variable in ["m", "p", "h"]:
            if self.get_attr(variable).is_var:
                self.var_pos[variable] = self.num_vars
                self.num_vars += 1

        if self.fluid.is_var:
            # if the fluid is not a variable, it does not count towards the
            # equations for this connection
            for fluid, x in self.fluid.val_set.items():
                if not x:
                    self.var_pos[fluid] = self.num_vars
                    self.num_vars += 1

        self.jacobian = np.zeros((self.num_eq, 1, self.num_vars))
        self.residual = np.zeros(self.num_eq)

    def get_parameters(self):
        return {
            "T": {
                "func": self.T_func, "deriv": self.T_deriv,
                "constant_deriv": False, "latex": None,
                "num_eq": 1
            },
            "v": {
                "func": self.v_func, "deriv": self.v_deriv,
                "constant_deriv": False, "latex": None,
                "num_eq": 1
            },
            "x": {
                "func": self.x_func, "deriv": self.x_deriv,
                "constant_deriv": False, "latex": None,
                "num_eq": 1
            },
            "Td_bp": {
                "func": self.Td_bp_func, "deriv": self.Td_bp_deriv,
                "constant_deriv": False, "latex": None,
                "num_eq": 1
            },
        }

    def build_fluid_data(self):
        self.fluid_data = {
            fluid: {
                "wrapper": self.fluid.wrapper[fluid],
                "mass_fraction": self.fluid.val[fluid]
            } for fluid in self.fluid.val
        }

    def T_func(self, k):
        self.residual[k] = self.T.val_SI - T_mix_ph(self.p.val_SI, self.h.val_SI, self.fluid_data)

    def T_deriv(self, k):
        if self.p.is_var:
            self.jacobian[k, 0, self.var_pos["p"]] = (
                dT_mix_dph(self.p.val_SI, self.h.val_SI, self.fluid_data)
            )
        if self.h.is_var:
            self.jacobian[k, 0, self.var_pos["h"]] = (
                dT_mix_pdh(self.p.val_SI, self.h.val_SI, self.fluid_data)
            )
        # if len(self.fluid.val) > 1:
        #     self.jacobian[0, 0, 3:] = dT_mix_ph_dfluid(
        #         c.p.val_SI, c.h.val_SI, self.fluid_data, T0=c.T.val_SI
        # )

    def v_func(self, k):
        self.residual[k] = (
            v_mix_ph(self.p.val_SI, self.h.val_SI, self.fluid_data)
            * self.m.val_SI - self.v.val_SI
        )

    def v_deriv(self, k):
        if self.m.is_var:
            self.jacobian[k, 0, self.var_pos["m"]] = v_mix_ph(self.p.val_SI, self.h.val_SI, self.fluid_data)
        if self.p.is_var:
            self.jacobian[k, 0, self.var_pos["p"]] = dv_mix_dph(self.p.val_SI, self.h.val_SI, self.fluid_data) * self.m.val_SI
        if self.h.is_var:
            self.jacobian[k, 0, self.var_pos["h"]] = dv_mix_pdh(self.p.val_SI, self.h.val_SI, self.fluid_data) * self.m.val_SI


    def x_func(self, k):
        # saturated steam fraction
        # how to check if condition?
        # if (np.absolute(self.residual[k]) > ERR ** 2 or
        #         self.iter % 2 == 0 or self.always_all_equations):
        self.residual[k] = self.h.val_SI - h_mix_pQ(self.p.val_SI, self.x.val_SI, self.fluid_data)

    def x_deriv(self, k):
        # if not self.increment_filter[col + 1]:
        if self.p.is_var:
            self.jacobian[k, 0, self.var_pos["p"]] = -dh_mix_dpQ(self.p.val_SI, self.x.val_SI, self.fluid_data)
        if self.h.is_var:
            self.jacobian[k, 0, self.var_pos["h"]] = 1

    def Td_bp_func(self, k):

        # temperature difference to boiling point
            # if (np.absolute(self.residual[k]) > ERR ** 2 or
            #         self.iter % 2 == 0 or self.always_all_equations):
        self.residual[k] = (
            T_mix_ph(self.p.val_SI, self.h.val_SI, self.fluid_data)
            - self.Td_bp.val_SI
            - T_sat_p(self.p.val_SI, self.fluid_data)
        )

    def Td_bp_deriv(self, k):
        # if not self.increment_filter[col + 1]:
        if self.p.is_var:
            self.jacobian[k, 0, self.var_pos["p"]] = (
                dT_mix_dph(self.p.val_SI, self.h.val_SI, self.fluid_data)
                - dT_sat_dp(self.p.val_SI, self.fluid_data)
            )
        if self.h.is_var:
            # if not self.increment_filter[col + 2]:
            self.jacobian[k, 0, self.var_pos["h"]] = dT_mix_pdh(
                self.p.val_SI, self.h.val_SI, self.fluid_data
            )

    def solve(self, increment_filter):
        k = 0
        for parameter in self.parameters:
            if self.get_attr(parameter).val_set:
                self.parameters[parameter]["func"](k)
                self.parameters[parameter]["deriv"](k)

                k += 1


from tespy.networks import Network
from tespy.components import Source, Sink, Pipe

nwk = Network(fluids=["water", "H2"], T_unit="C", p_unit="bar")

a = Source("source")
b = Sink("sink")

c1 = NewConnection(a, "out1", b, "in1", label="1")

c = Source("source2")
d = Pipe("pipe")
e = Sink("sink2")

c2 = NewConnection(c, "out1", d, "in1", label="2")
c3 = NewConnection(d, "out1", e, "in1", label="3")

nwk.add_conns(c1, c2, c3)

c1.set_attr(m=1, x=1, T=150, fluid={"water": 1, "H2": 0})

c2.set_attr(m=1, p=3, Td_bp=3)
d.set_attr(pr=1, Q="var")
c3.set_attr(Td_bp=5, fluid={"water": 1, "H2": 0})

nwk.solve("design")

nwk.print_results()

c2.set_attr(p=None, T=136)
nwk.solve("design")
nwk.print_results()