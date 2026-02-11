# -*- coding: utf-8
"""Module of class Connection and class Ref.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/connections/humidairconnection.py

SPDX-License-Identifier: MIT
"""

import numpy as np
from CoolProp.CoolProp import HAPropsSI

from tespy.tools import fluid_properties as fp
from tespy.tools.data_containers import FluidComposition as dc_flu
from tespy.tools.data_containers import FluidProperties as dc_prop
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
from tespy.tools.fluid_properties.functions import w_mix_pT_humidair
from tespy.tools.fluid_properties.functions import h_mix_pT
from tespy.tools.fluid_properties.mixtures import _get_fluid_alias, w_mix_fluid_data
from tespy.tools.helpers import seeded_random

from .connection import Connection
from .connection import connection_registry


@connection_registry
class HAConnection(Connection):

    def get_parameters(self):
        return {
            "fluid": dc_flu(
                d=1e-5,
                description="mass fractions of the fluid composition"
            ),
            "m": dc_prop(
                quantity="mass_flow",
                description="mass flow of dry air (system variable)"
            ),
            "mHA": dc_prop(
                quantity="mass_flow",
                description="mass flow of humid air"
            ),
            "mH2O": dc_prop(
                quantity="mass_flow",
                description="mass flow of liquid or solid water not contained in humid air"
            ),
            "p": dc_prop(
                quantity="pressure",
                description="absolute pressure of the fluid (system variable)"
            ),
            "h": dc_prop(
                quantity="enthalpy",
                description="dry air mass specific enthalpy (system variable)"
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
            "r": dc_prop(
                func=self.r_func,
                dependents=self.r_dependents,
                num_eq=1,
                quantity="ratio",
                description="relative humidity"
            ),
            "fluid_balance": dc_simple(
                func=self.fluid_balance_func,
                deriv=self.fluid_balance_deriv,
                _val=False, num_eq_sets=1,
                dependents=self.fluid_balance_dependents,
                description="apply an equation which closes the fluid balance with at least two unknown fluid mass fractions"
            ),
        }

    def _parameter_specification(self, key, value):
        if key == "w" or key == "w0":
            # specification of w is equivalent to specification of fluid
            # composition for humid air
            air = 1 / (1 + value)
            if key == "w":
                self.set_attr(fluid={"air": air, "water": 1 - air})
            else:
                self.set_attr(fluid0={"air": air, "water": 1 - air})
        else:
            super()._parameter_specification(key, value)

    # for HAConnection mixing rule cannot be modified, is always humidair
    def _get_mixing_rule(self):
        return "humidair"

    def _set_mixing_rule(self, value):
        if value is not None and value != self.mixing_rule:
            print(value)
            msg = (
                "You cannot change the mixing rule specification for a "
                f"Connection of type {self.__class__.__name__}"
            )
            raise ValueError(msg)

    mixing_rule = property(_get_mixing_rule, _set_mixing_rule)

    def _guess_starting_values(self, units):
        if self.h.is_var and not self.good_starting_values:
            value = seeded_random(self.label)
            T_rand = 280 + value * (300 - 280)
            h = fp.h_mix_pT(1e5, T_rand, self.fluid_data, self.mixing_rule)
            self.h.set_reference_val_SI(h)
            self._precalc_guess_values()

    def _precalc_guess_values(self):
        if not self.h.is_var:
            return

        if not self.good_starting_values:
            if self.T.is_set:
                try:
                    w = self.calc_w()
                    self.h.set_reference_val_SI(
                        HAPropsSI("H", "P", self.p.val_SI, "T", self.T.val_SI, "W", w)
                    )
                except ValueError:
                    pass

    def _presolve(self):

        air_alias = _get_fluid_alias("air", self.fluid_data)
        water_alias = _get_fluid_alias("water", self.fluid_data)
        if not air_alias:
            msg = "air must be present in fluid composition"
            raise ValueError(msg)

        elif not water_alias:
            msg = "water must be present in fluid composition"
            raise ValueError(msg)

        if len(self.fluid.is_var) > 0:
            return []

        presolved_equations = []

        if self.h.is_var and not self.p.is_var:
            if self.T.is_set:
                self.h.set_reference_val_SI(h_mix_pT(self.p.val_SI, self.T.val_SI, self.fluid_data, self.mixing_rule))
                self.h._potential_var = False
                if "T" in self._equation_set_lookup.values():
                    presolved_equations += ["T"]

        presolved_equations = [
            key for parameter in presolved_equations
            for key, value in self._equation_set_lookup.items()
            if value == parameter
        ]
        return presolved_equations

    def _adjust_to_property_limits(self, nw):

        if self.p.is_var:
            if self.p.val_SI < 100:
                self.p.val_SI = 101
            elif self.p.val_SI > 100e5:
                self.p.val_SI = 99e5

        if self.h.is_var:
            # TODO: check minimum temperature how it matches minimum humidity ratio
            d = self.h._reference_container._d
            hmin = HAPropsSI("H", "T", -50 + 273.15, "P", self.p.val_SI, "R", 1)
            if self.h.val_SI < hmin:
                delta = max(abs(self.h.val_SI * d), d) * 5
                self.h.set_reference_val_SI(hmin + delta)

            else:
                # TODO: where to get reasonable hmax from?!
                hmax = HAPropsSI("H", "T", 300 + 273.15, "P", self.p.val_SI, "R", 0)
                if self.h.val_SI > hmax:
                    delta = max(abs(self.h.val_SI * d), d) * 5
                    self.h.set_reference_val_SI(hmax - delta)

    @classmethod
    def _result_attributes(cls):
        return ["m", "mHA", "mH2O", "p", "h", "T", "w", "s", "vol", "v", "r"]

    @classmethod
    def _print_attributes(cls):
        return ["m", "mHA", "mH2O", "p", "h", "T", "w", "r"]

    def calc_r(self):
        w = self.calc_w()
        try:
            return HAPropsSI("R", "P", self.p.val_SI, "T", self.T.val_SI, "W", w)
        except ValueError as e:
            value = str(e).split("value (")[1].split(")")[0]
            return float(value)

    def r_func(self):
        return self.r.val_SI - self.calc_r()

    def r_dependents(self):
        water_alias = _get_fluid_alias("H2O", self.fluid_data)
        # water alias is already a set
        return {
            "scalars": [self.p, self.h],
            "vectors": [{self.fluid: water_alias}]
        }

    def calc_w(self):
        return w_mix_pT_humidair(self.p.val_SI, self.T.val_SI, self.fluid_data)

    def calc_results(self, units):
        self.T.val_SI = self.calc_T()
        self.vol.val_SI = self.calc_vol()  # Mixture volume per mass of dry air
        self.v.val_SI = self.vol.val_SI * self.m.val_SI  # Mixture volume flow rate
        # handle the water fraction
        self.w.val_SI = self.calc_w()
        # # Convert from kg water/kg dry air to mass fraction of water in humid air
        # x_h2o = self.w.val_SI / (1 + self.w.val_SI)
        # x_air = 1 - x_h2o
        self.mHA.val_SI = self.m.val_SI * (1 + self.w.val_SI)
        w_mixture = w_mix_fluid_data(self.fluid_data)
        # Calculate kg water/kg dry air that is not in the humid air
        delta_w = w_mixture - self.w.val_SI
        self.mH2O.val_SI = self.m.val_SI * delta_w       
        self.r.val_SI = self.calc_r()
        # if self.r.val_SI > 1:
        #     self.r.val_SI = np.nan

        for prop in self._result_attributes():
            param = self.get_attr(prop)
            if not param.is_set:
                param.set_val_from_SI(units)

        self.m.set_val0_from_SI(units)
        self.p.set_val0_from_SI(units)
        self.h.set_val0_from_SI(units)
        self.fluid.val0 = self.fluid.val.copy()

        return True
