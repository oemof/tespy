from tespy.components.displacementmachinery.polynomial_compressor import PolynomialCompressor
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.helpers import TESPyComponentError


class PolynomialCompressorWithCooling(PolynomialCompressor):

    def _preprocess(self, row_idx):
        if not self.eta_recovery.is_set:
            msg = (
                f"The component {self.label} of type {self.__class__.__name__}"
                "requires you to specify the share of heat recovery "
                "eta_recovery."
            )
            raise TESPyComponentError(msg)

        return super()._preprocess(row_idx)

    @staticmethod
    def inlets():
        return ['in1', 'in2']

    @staticmethod
    def outlets():
        return ['out1', 'out2']

    def get_mandatory_constraints(self) -> dict:
        constraints = super().get_mandatory_constraints()
        # this is a dictionary
        constraints["cooling_energy_balance_constraints"] = dc_cmc(
            func=self.cooling_energy_balance_func,
            dependents=self.cooling_energy_balance_dependents,
            num_eq_sets=1
        )
        return constraints

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

        self.dp_cooling.val_SI = self.inl[1].p.val_SI - self.outl[1].p.val_SI
