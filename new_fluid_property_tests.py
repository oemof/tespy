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

from tespy.connections import Connection, Ref

class NewConnection(Connection):

    def __init__(self, comp1, outlet_id, comp2, inlet_id, label=None, **kwargs):
        super().__init__(comp1, outlet_id, comp2, inlet_id, label, **kwargs)

    def set_attr(self, **kwargs):
        super().set_attr(**kwargs)
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
        self.num_eq = 0
        for parameter in self.parameters:
            if self.get_attr(parameter).val_set:
                self.num_eq += self.parameters[parameter].num_eq

        self.residual = np.zeros(self.num_eq)
        self.jacobian = OrderedDict()

    def get_parameters(self):
        return {
            'm': dc_prop(),
            'p': dc_prop(),
            'h': dc_prop(),
            'vol': dc_prop(),
            's': dc_prop(),
            'fluid': dc_flu(),
            'state': dc_simple(),
            "T": dc_prop(**{
                "func": self.T_func, "deriv": self.T_deriv, "num_eq": 1
            }),
            "v": dc_prop(**{
                "func": self.v_func, "deriv": self.v_deriv, "num_eq": 1
            }),
            "x": dc_prop(**{
                "func": self.x_func, "deriv": self.x_deriv, "num_eq": 1
            }),
            "Td_bp": dc_prop(**{
                "func": self.Td_bp_func, "deriv": self.Td_bp_deriv, "num_eq": 1
            }),
            "m_ref": dc_prop(**{
                "func": self.primary_ref_func, "deriv": self.primary_ref_deriv,
                "num_eq": 1, "func_params": {"variable": "m"}
            }),
            "p_ref": dc_prop(**{
                "func": self.primary_ref_func, "deriv": self.primary_ref_deriv,
                "num_eq": 1, "func_params": {"variable": "p"}
            }),
            "h_ref": dc_prop(**{
                "func": self.primary_ref_func, "deriv": self.primary_ref_deriv,
                "num_eq": 1, "func_params": {"variable": "h"}
            }),
        }

    def build_fluid_data(self):
        self.fluid_data = {
            fluid: {
                "wrapper": self.fluid.wrapper[fluid],
                "mass_fraction": self.fluid.val[fluid]
            } for fluid in self.fluid.val
        }

    def primary_ref_func(self, k, **kwargs):
        variable = kwargs['variable']
        self.get_attr(variable)
        ref = self.get_attr(f"{variable}_ref").val
        self.residual[k] = (
            self.get_attr(variable).val_SI
            - ref.obj.get_attr(variable).val_SI * ref.factor + ref.delta
        )

    def primary_ref_deriv(self, k, **kwargs):
        variable = kwargs['variable']
        ref = self.get_attr(f"{variable}_ref").val
        if self.get_attr(variable).is_var:
            self.jacobian[k, self.get_attr(variable).J_col] = 1

        if ref.obj.get_attr(variable).is_var:
            self.jacobian[k, ref.obj.get_attr(variable).J_col] = -ref.factor

    def calc_T(self):
        return T_mix_ph(self.p.val_SI, self.h.val_SI, self.fluid_data)

    def T_func(self, k, **kwargs):
        self.residual[k] = self.T.val_SI - self.calc_T()

    def T_deriv(self, k, **kwargs):
        if self.p.is_var:
            self.jacobian[k, self.p.J_col] = (
                -dT_mix_dph(self.p.val_SI, self.h.val_SI, self.fluid_data)
            )
        if self.h.is_var:
            self.jacobian[k, self.h.J_col] = (
                -dT_mix_pdh(self.p.val_SI, self.h.val_SI, self.fluid_data)
            )
        # if len(self.fluid.val) > 1:
        #     self.jacobian[0, 0, 3:] = dT_mix_ph_dfluid(
        #         c.p.val_SI, c.h.val_SI, self.fluid_data, T0=c.T.val_SI
        # )

    def calc_vol(self):
        try:
            vol = v_mix_ph(self.p.val_SI, self.h.val_SI, self.fluid_data, T0=self.T.val_SI)
        except NotImplementedError:
            vol = np.nan
        return vol

    def v_func(self, k, **kwargs):
        self.residual[k] = (
            v_mix_ph(self.p.val_SI, self.h.val_SI, self.fluid_data)
            * self.m.val_SI - self.v.val_SI
        )

    def v_deriv(self, k, **kwargs):
        if self.m.is_var:
            self.jacobian[k, self.m.J_col] = v_mix_ph(self.p.val_SI, self.h.val_SI, self.fluid_data)
        if self.p.is_var:
            self.jacobian[k, self.p.J_col] = dv_mix_dph(self.p.val_SI, self.h.val_SI, self.fluid_data) * self.m.val_SI
        if self.h.is_var:
            self.jacobian[k, self.h.J_col] = dv_mix_pdh(self.p.val_SI, self.h.val_SI, self.fluid_data) * self.m.val_SI

    def calc_x(self):
        try:
            x = Q_mix_ph(self.p.val_SI, self.h.val_SI, self.fluid_data)
        except NotImplementedError:
            x = np.nan
        return x

    def x_func(self, k, **kwargs):
        # saturated steam fraction
        # how to check if condition?
        # if (np.absolute(self.residual[k]) > ERR ** 2 or
        #         self.iter % 2 == 0 or self.always_all_equations):
        self.residual[k] = self.h.val_SI - h_mix_pQ(self.p.val_SI, self.x.val_SI, self.fluid_data)

    def x_deriv(self, k, **kwargs):
        # if not self.increment_filter[col + 1]:
        if self.p.is_var:
            self.jacobian[k, self.p.J_col] = -dh_mix_dpQ(self.p.val_SI, self.x.val_SI, self.fluid_data)
        if self.h.is_var: # and self.it == 0
            self.jacobian[k, self.h.J_col] = 1

    def calc_Td_bp(self):
        try:
            Td_bp = self.calc_T() - T_sat_p(self.p.val_SI, self.fluid_data)
        except NotImplementedError:
            Td_bp = np.nan
        return Td_bp

    def Td_bp_func(self, k, **kwargs):

        # temperature difference to boiling point
            # if (np.absolute(self.residual[k]) > ERR ** 2 or
            #         self.iter % 2 == 0 or self.always_all_equations):
        self.residual[k] = self.calc_Td_bp() - self.Td_bp.val_SI

    def Td_bp_deriv(self, k, **kwargs):
        # if not self.increment_filter[col + 1]:
        if self.p.is_var:
            self.jacobian[k, self.p.J_col] = (
                dT_mix_dph(self.p.val_SI, self.h.val_SI, self.fluid_data)
                - dT_sat_dp(self.p.val_SI, self.fluid_data)
            )
        if self.h.is_var:
            # if not self.increment_filter[col + 2]:
            self.jacobian[k, self.h.J_col] = dT_mix_pdh(
                self.p.val_SI, self.h.val_SI, self.fluid_data
            )

    def calc_s(self):
        return s_mix_ph(self.p.val_SI, self.h.val_SI, self.fluid_data, T0=self.T.val_SI)

    def solve(self, increment_filter):
        k = 0
        for parameter, data in self.parameters.items():
            if self.get_attr(parameter).val_set:
                data.func(k, **data.func_params)
                data.deriv(k, **data.func_params)

                k += 1

    def calc_results(self):
        self.T.val_SI = self.calc_T()
        number_fluids = get_number_of_fluids(self.fluid_data)
        if (
                number_fluids > 1
                and abs(
                    h_mix_pT(self.p.val_SI, self.T.val_SI, self.fluid_data)
                    - self.h.val_SI
                ) > ERR ** .5
        ):
            self.T.val_SI = np.nan
            self.vol.val_SI = np.nan
            self.v.val_SI = np.nan
            self.s.val_SI = np.nan
            msg = (
                'Could not find a feasible value for mixture temperature '
                'at connection ' + c.label + '. The values for '
                'temperature, specific volume, volumetric flow and '
                'entropy are set to nan.')
            logger.error(msg)
        else:
            self.vol.val_SI = self.calc_vol()
            self.v.val_SI = self.vol.val_SI * self.m.val_SI
            self.s.val_SI = self.calc_s()

        if number_fluids == 1:
            if not self.x.val_set:
                self.x.val_SI = self.calc_x()
            if not self.Td_bp.val_set:
                self.Td_bp.val_SI = self.calc_Td_bp()

            for prop in fpd.keys():
                self.get_attr(prop).val = convert_from_SI(
                    prop, self.get_attr(prop).val_SI, self.get_attr(prop).unit
                )

        self.m.val0 = self.m.val
        self.p.val0 = self.p.val
        self.h.val0 = self.h.val
        self.fluid.val0 = self.fluid.val.copy()

    def check_pressure_bounds(self, fluid):
        if self.p.val_SI < self.fluid.wrapper[fluid]._p_min:
            c.p.val_SI = self.fluid.wrapper[fluid]._p_min * 1.01
            logger.debug(self._property_range_message('p'))
        elif self.p.val_SI > self.fluid.wrapper[fluid]._p_max:
            self.p.val_SI = self.fluid.wrapper[fluid]._p_max
            logger.debug(self._property_range_message('p'))

    def check_enthalpy_bounds(self, fluid):
        # enthalpy
        try:
            hmin = self.fluid.wrapper[fluid].h_pT(
                self.p.val_SI, self.fluid.wrapper[fluid]._T_min
            )
        except ValueError:
            f = 1.05
            hmin = self.fluid.wrapper[fluid].h_pT(
                self.p.val_SI, self.fluid.wrapper[fluid]._T_min * f
            )

        T = self.fluid.wrapper[fluid]._T_max
        while True:
            try:
                hmax = self.fluid.wrapper[fluid].h_pT(self.p.val_SI, T)
                break
            except ValueError as e:
                T *= 0.99
                if T < self.fluid.wrapper[fluid]._T_min:
                    raise ValueError(e) from e

        if self.h.val_SI < hmin:
            if hmin < 0:
                self.h.val_SI = hmin * 0.9999
            else:
                self.h.val_SI = hmin * 1.0001
            logger.debug(self._property_range_message('h'))

        elif self.h.val_SI > hmax:
            self.h.val_SI = hmax * 0.9999
            logger.debug(self._property_range_message('h'))

    def check_two_phase_bounds(self, fluid):

        if (self.Td_bp.val_SI > 0 or (self.state.val == 'g' and self.state.is_set)):
            h = self.fluid.wrapper[fluid].h_pQ(self.p.val_SI, 1)
            if self.h.val_SI < h:
                self.h.val_SI = h * 1.01
                logger.debug(self._property_range_message('h'))
        elif (self.Td_bp.val_SI < 0 or (self.state.val == 'l' and self.state.is_set)):
            h = self.fluid.wrapper[fluid].h_pQ(self.p.val_SI, 0)
            if self.h.val_SI > h:
                self.h.val_SI = h * 0.99
                logger.debug(self._property_range_message('h'))

    def check_temperature_bounds(self):
        r"""
        Check if temperature is within user specified limits.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            Connection to check fluid properties.
        """
        Tmin = max(
            [w._T_min for f, w in self.fluid.wrapper.items() if self.fluid.val[f] > ERR]
        ) * 1.01
        Tmax = min(
            [w._T_max for f, w in self.fluid.wrapper.items() if self.fluid.val[f] > ERR]
        ) * 0.99
        hmin = h_mix_pT(self.p.val_SI, Tmin, self.fluid_data)
        hmax = h_mix_pT(self.p.val_SI, Tmax, self.fluid_data)

        if self.h.val_SI < hmin:
            self.h.val_SI = hmin
            logger.debug(self._property_range_message('h'))

        if self.h.val_SI > hmax:
            self.h.val_SI = hmax
            logger.debug(self._property_range_message('h'))

    def _property_range_message(self, prop):
        r"""
        Return debugging message for fluid property range adjustments.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            Connection to check fluid properties.

        prop : str
            Fluid property.

        Returns
        -------
        msg : str
            Debugging message.
        """
        msg = (
            f"{fpd[prop]['text'][0].upper()} {fpd[prop]['text'][1:]} out of "
            f"fluid property range at connection {self.label}, adjusting value "
            f"to {self.get_attr(prop).val_SI} {fpd[prop]['SI_unit']}."
        )
        return msg

# def conceptual_test():

from tespy.networks import Network
from tespy.components import Source, Sink, Pipe, Splitter, Turbine


nwk = Network(fluids=["Water", "H2", "R134a"], T_unit="C", p_unit="bar", h_unit='kJ / kg')

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
import tespy.components.turbomachinery.turbine as t
t.isentropic = isentropic


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
from tespy.components import Source, Sink, Merge, Splitter, HeatExchangerSimple
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