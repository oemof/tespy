# %%[sec_1]
from tespy.components import PolynomialCompressor
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.networks import Network


class PolynomialCompressorWithCooling(PolynomialCompressor):

    @staticmethod
    def inlets():
        return ['in1', 'in2']

    @staticmethod
    def outlets():
        return ['out1', 'out2']
# %%[sec_2]
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


c1.set_attr(fluid={"R290": 1}, m=1, T_dew=10, td_dew=10)
c2.set_attr(T_dew=60, td_dew=50)

b1.set_attr(fluid={"water": 1}, m=0.5, T=15, p=1)
b2.set_attr(T=25, p=1)

compressor.set_attr(dissipation_ratio=0.1)

nw.solve("design")
nw.print_results()

b1.fluid.val, b2.fluid.val
# %%[sec_3]
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc


def get_mandatory_constraints(self) -> dict:
    constraints = super().get_mandatory_constraints()
    # this is a dictionary
    constraints["cooling_energy_balance_constraints"] = dc_cmc(
        func=self.cooling_energy_balance_func,
        dependents=self.cooling_energy_balance_dependents,
        num_eq_sets=1
    )
    return constraints
# %%[sec_4]
def _preprocess(self, row_idx):
    if not self.dissipation_ratio.is_set:
        self.dissipation_ratio.is_set = True
        self.dissipation_ratio.val = 0
        self.dissipation_ratio.val_SI = 0
    return super()._preprocess(row_idx)

def cooling_energy_balance_func(self):
    usable_share = 0.9
    residual = (
        self.inl[1].m.val_SI * (self.outl[1].h.val_SI - self.inl[1].h.val_SI)
        + self.inl[0].m.val_SI * (
            self.outl[0].h.val_SI
            - self.outl[0].h.val_SI / (1 - self.dissipation_ratio.val_SI)
            + self.inl[0].h.val_SI * (self.dissipation_ratio.val_SI / (1 - self.dissipation_ratio.val_SI))
        ) * usable_share
    )
    return residual

def cooling_energy_balance_dependents(self):
    return [
        self.inl[0].m, self.inl[1].m,
        self.inl[0].h, self.inl[1].h,
        self.outl[0].h, self.outl[1].h
    ]
# %%[sec_5]


