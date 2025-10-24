# %%[sec_1]
from tespy.components import PolynomialCompressor
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.networks import Network
# %%[sec_2]
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.helpers import TESPyComponentError


# %%[sec_3]
class PolynomialCompressorWithCooling(PolynomialCompressor):

    @staticmethod
    def inlets():
        return ['in1', 'in2']

    @staticmethod
    def outlets():
        return ['out1', 'out2']
# %%[sec_4]
    def get_mandatory_constraints(self) -> dict:
        constraints = super().get_mandatory_constraints()
        # this is a dictionary
        constraints["cooling_energy_balance_constraints"] = dc_cmc(
            func=self.cooling_energy_balance_func,
            dependents=self.cooling_energy_balance_dependents,
            num_eq_sets=1
        )
        return constraints
# %%[sec_5]
    def _preprocess(self, row_idx):
        if not self.eta_recovery.is_set:
            msg = (
                f"The component {self.label} of type {self.__class__.__name__}"
                "requires you to specify the share of heat recovery "
                "eta_recovery."
            )
            raise TESPyComponentError(msg)

        return super()._preprocess(row_idx)

    def cooling_energy_balance_func(self):
        residual = (
            self.inl[1].m.val_SI * (self.outl[1].h.val_SI - self.inl[1].h.val_SI)
            + self.inl[0].m.val_SI * (
                self.outl[0].h.val_SI
                - self.outl[0].h.val_SI / (1 - self.dissipation_ratio.val_SI)
                + self.inl[0].h.val_SI * (
                    self.dissipation_ratio.val_SI / (1 - self.dissipation_ratio.val_SI)
                )
            ) * self.eta_recovery.val_SI
        )
        return residual

    def cooling_energy_balance_dependents(self):
        return [
            self.inl[0].m, self.inl[1].m,
            self.inl[0].h, self.inl[1].h,
            self.outl[0].h, self.outl[1].h
        ]
# %%[sec_6]
    def get_parameters(self):
        params = super().get_parameters()
        params["eta_recovery"] = dc_cp(
            quantity="efficiency"
        )
        params["td_minimal"] = dc_cp(
            min_val=0,
            quantity="temperature_difference"
        )
        params["dp_cooling"] = dc_cp(
            min_val=0,
            structure_matrix=self.dp_structure_matrix,
            func_params={"inconn": 1, "outconn": 1, "dp": "dp_cooling"},
            quantity="pressure"
        )
        return params
# %%[sec_7]
    def calc_parameters(self):
        super().calc_parameters()

        i = self.inl[0]
        o = self.outl[0]
        h_2 = (
            (o.h.val_SI - i.h.val_SI * self.dissipation_ratio.val_SI)
            / (1 - self.dissipation_ratio.val_SI)
        )
        T_max_compressor_internal = T_mix_ph(
            self.outl[0].p.val_SI,
            h_2,
            self.outl[0].fluid_data,
            self.outl[0].mixing_rule,
            T0=self.outl[0].T.val_SI
        )
        self.td_minimal.val_SI = (
            T_max_compressor_internal
            - self.outl[1].T.val_SI
        )
# %%[sec_8]
        self.dp_cooling.val_SI = self.inl[1].p.val_SI - self.outl[1].p.val_SI
# %%[sec_9]
nw = Network()
nw.units.set_defaults(
    temperature="°C",
    pressure="bar"
)

gas_inlet = Source("gas inlet")
gas_outlet = Sink("gas outlet")

water_inlet = Source("water cold")
water_outlet = Sink("water hot")

compressor = PolynomialCompressorWithCooling("compressor")

c1 = Connection(gas_inlet, "out1", compressor, "in1", label="c1")
c2 = Connection(compressor, "out1", gas_outlet, "in1", label="c2")

b1 = Connection(water_inlet, "out1", compressor, "in2", label="b1")
b2 = Connection(compressor, "out2", water_outlet, "in1", label="b2")

nw.add_conns(c1, c2, b1, b2)
# %%[sec_10]
c1.set_attr(fluid={"R290": 1}, m=1, T_dew=10, td_dew=10)
c2.set_attr(T_dew=60, td_dew=50)

b1.set_attr(fluid={"water": 1}, m=0.5, T=15, p=1)
b2.set_attr(T=25, p=1)

compressor.set_attr(dissipation_ratio=0.1)

nw.solve("design")
nw.print_results()

b1.fluid.val, b2.fluid.val
# %%[sec_11]
c1.set_attr(fluid={"R290": 1}, m=1, T_dew=10, td_dew=10)
c2.set_attr(T_dew=60, td_dew=50)

b1.set_attr(fluid={"water": 1}, m=1, T=15, p=1)
b2.set_attr(T=25, p=1)
compressor.set_attr(dissipation_ratio=0.1)

nw.solve("design")
# %%[sec_12]
b1.set_attr(m=None)
nw.solve("design")
b1.m.val_SI
# %%[sec_13]
c1.set_attr(fluid={"R290": 1}, m=1, T_dew=10, td_dew=10)
c2.set_attr(T_dew=60, td_dew=25)

b1.set_attr(fluid={"water": 1}, T=15, p=1)
b2.set_attr(T=25, p=1)
compressor.set_attr(dissipation_ratio=0.1, eta_recovery=0.9)

nw.solve("design")
compressor.Q_diss.val

b1.m.val_SI * (b2.h.val_SI - b1.h.val_SI)
# %%[sec_14]
b1.set_attr(m=0.005)
b2.set_attr(T=None)
nw.solve("design")

b2.T.val, c2.T.val
h_2 = c1.h.val_SI + (c2.h.val_SI - c1.h.val_SI) / (1 - compressor.dissipation_ratio.val_SI)
c2.p.val_SI
# %%[sec_15]
T_mix_ph(c2.p.val_SI, h_2, c2.fluid_data) - 273.15
compressor.eta_s.val
# %%[sec_16]
c1.set_attr(fluid={"R290": 1}, m=1, T_dew=10, td_dew=10)
c2.set_attr(T_dew=60, td_dew=25)

b1.set_attr(fluid={"water": 1}, T=15, m=0.05, p=1)
b2.set_attr(p=1)
compressor.set_attr(dissipation_ratio=0.1, eta_recovery=0.9)

nw.solve("design")
# %%[sec_17]
c1.set_attr(fluid={"R290": 1}, m=1, T_dew=10, td_dew=10)
c2.set_attr(T_dew=60, td_dew=25)

b1.set_attr(fluid={"water": 1}, T=15, m=0.05, p=1)
b2.set_attr(p=0.9)
compressor.set_attr(dissipation_ratio=0.1, eta_recovery=0.9)

nw.solve("design")
nw.print_results()
# %%[sec_18]


