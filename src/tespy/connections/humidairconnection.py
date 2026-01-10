from CoolProp.CoolProp import HAPropsSI

from tespy.tools import fluid_properties as fp
from tespy.tools.data_containers import FluidComposition as dc_flu
from tespy.tools.data_containers import FluidProperties as dc_prop

from .connection import Connection
from .connection import connection_registry


@connection_registry
class HumidAirConnection(Connection):

    def get_variables(self):
        return {
            "m": self.m, "p": self.p, "h": self.h, "w": self.w
        }

    def get_parameters(self):
        return {
            "fluid": dc_flu(
                d=1e-5,
                description="mass fractions of the fluid composition"
            ),
            "m": dc_prop(
                quantity="mass_flow",
                description="mass flow of the fluid (system variable)"
            ),
            "p": dc_prop(
                quantity="pressure",
                description="absolute pressure of the fluid (system variable)"
            ),
            "h": dc_prop(
                quantity="enthalpy",
                description="mass specific enthalpy of the fluid (system variable)"
            ),
            "w": dc_prop(
                quantity="ratio",
                description="mass of water per mass of dry air"
            ),
            "T": dc_prop(
                func=self.T_func,
                dependents=self.T_dependents,
                num_eq=1,
                quantity="temperature",
                description="temperature of the fluid"
            ),
            "v": dc_prop(
                func=self.v_func,
                dependents=self.v_dependents,
                num_eq=1,
                quantity="volumetric_flow",
                description="volumetric flow of the fluid"
            ),
            "vol": dc_prop(
                quantity="specific_volume",
                description="specific volume of the fluid (output only)"
            ),
            "s": dc_prop(
                quantity="entropy",
                description="specific entropy of the fluid (output only)"
            ),
            # "r": dc_prop(
            #     func=self.r_func,
            #     dependents=self.r_dependents,
            #     quantity="ratio",
            #     description="relative humidity"
            # )
        }

    def calc_T(self, T0=None):
        return HAPropsSI(
            "T", "P", self.p.val_SI, "H", self.h.val_SI, "W", self.w.val_SI / 1e3
        )

    def T_dependents(self):
        return [self.p, self.h, self.w]

    def calc_vol(self, T0=None):
        return HAPropsSI(
            "V", "P", self.p.val_SI, "H", self.h.val_SI, "W", self.w.val_SI / 1e3
        )

    def v_dependents(self):
        return [self.m, self.p, self.h, self.w]

    def _guess_starting_values(self, units):
        self.h.set_reference_val_SI(4e5)
        self.w.set_reference_val_SI(5)
        self._precalc_guess_values()

    def _precalc_guess_values(self):
        if not self.h.is_var:
            return

        if not self.good_starting_values:
            if self.T.is_set:
                try:
                    self.h.set_reference_val_SI(
                        HAPropsSI("H", "P", self.p.val_SI, "T", self.T.val_SI, "W", self.w.val_SI / 1e3)
                    )
                except ValueError:
                    pass

    def _presolve(self):
        return []

    def _adjust_to_property_limits(self, nw):
        fl = fp.single_fluid(self.fluid_data)
        # if self.h.is_var:
        #     self._adjust_enthalpy(fl)

    @classmethod
    def _result_attributes(cls):
        return ["m", "p", "h", "T", "w", "s", "vol", "v"]

    @classmethod
    def _print_attributes(cls):
        return ["m", "p", "h", "T", "w"]

    def calc_results(self, units):
        self.T.val_SI = self.calc_T()
        self.vol.val_SI = self.calc_vol()
        self.v.val_SI = self.vol.val_SI * self.m.val_SI

        for prop in self._result_attributes():
            param = self.get_attr(prop)
            if not param.is_set:
                param.set_val_from_SI(units)

        self.m.set_val0_from_SI(units)
        self.p.set_val0_from_SI(units)
        self.h.set_val0_from_SI(units)
        self.w.set_val0_from_SI(units)
        self.fluid.val0 = self.fluid.val.copy()

        return True
