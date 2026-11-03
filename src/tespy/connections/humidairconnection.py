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
from tespy.tools.data_containers import ConnectionParameter
from tespy.tools.data_containers import FluidComposition
from tespy.tools.data_containers import OptionalEquation
from tespy.tools.fluid_properties.functions import h_mix_pT
from tespy.tools.fluid_properties.mixtures import _get_fluid_alias
from tespy.tools.fluid_properties.mixtures import w_mix_fluid_data
from tespy.tools.fluid_properties.mixtures import w_mix_pT_humidair
from tespy.tools.fluid_properties.mixtures import w_mix_pTrh_humidair
from tespy.tools.helpers import seeded_random

from .connection import Connection
from .connection import connection_registry


@connection_registry
class HAConnection(Connection):
    r"""
    Connection for humid air with the dry air mass flow as system variable.

    The fluid composition is a mixture of dry air and water, the moisture
    content is specified through the humidity ratio :code:`w` or the
    relative humidity :code:`r`. Mass specific properties refer to the mass
    of dry air.

    Parameters
    ----------

    design : list
        List containing design parameters (stated as String).

    design_path : str
        Path to the individual design case for this connection.

    e_M : float
        Specific mechanical exergy. Quantity: :code:`specific_energy`. Optional
        result - calculated on request only, e.g. by an exergy evaluation.

    E_M : float
        Mechanical exergy flow. Quantity: :code:`power`. Optional result -
        calculated on request only, e.g. by an exergy evaluation.

    e_PH : float
        Specific physical exergy. Quantity: :code:`specific_energy`. Optional
        result - calculated on request only, e.g. by an exergy evaluation.

    E_PH : float
        Physical exergy flow. Quantity: :code:`power`. Optional result -
        calculated on request only, e.g. by an exergy evaluation.

    e_T : float
        Specific thermal exergy. Quantity: :code:`specific_energy`. Optional
        result - calculated on request only, e.g. by an exergy evaluation.

    E_T : float
        Thermal exergy flow. Quantity: :code:`power`. Optional result -
        calculated on request only, e.g. by an exergy evaluation.

    fluid : dict
        Mass fractions of the fluid composition.

    fluid_balance : bool
        Apply an equation which closes the fluid balance with at least two
        unknown fluid mass fractions.
        Equation: :py:meth:`fluid_balance_func <tespy.connections.connection.Connection.fluid_balance_func>`.

    h : float
        Dry air mass specific enthalpy (system variable).

    label : str
        The label of the connection.

    local_design : bool
        Treat this connection in design mode in an offdesign calculation.

    local_offdesign : bool
        Treat this connection in offdesign mode in a design calculation.

    m : float
        Mass flow of dry air (system variable).

    mH2O : float
        Mass flow of liquid or solid water not contained in humid air. Result
        only - calculated after each solve.

    mHA : float
        Mass flow of humid air. Result only - calculated after each solve.

    offdesign : list
        List containing offdesign parameters (stated as String).

    p : float
        Absolute pressure of the fluid (system variable).

    printout : bool
        Include this connection in the network's results printout.

    r : float
        Relative humidity.
        Equation: :py:meth:`r_func <tespy.connections.humidairconnection.HAConnection.r_func>`.

    s : float
        Specific entropy of the fluid. Result only - calculated after each
        solve.

    T : float
        Temperature of the fluid.
        Equation: :py:meth:`T_func <tespy.connections.connection.Connection.T_func>`.

    v : float
        Volumetric flow of the fluid.
        Equation: :py:meth:`v_func <tespy.connections.connection.Connection.v_func>`.

    vol : float
        Specific volume of the fluid. Result only - calculated after each solve.

    w : float
        Mass of water per mass of dry air.
    """

    _state_specifications = ("p", "h", "T", "w")

    def get_parameters(self):
        return {
            "fluid": FluidComposition(
                d=1e-5,
                description="mass fractions of the fluid composition"
            ),
            "m": ConnectionParameter(
                quantity="mass_flow",
                description="mass flow of dry air (system variable)"
            ),
            "mHA": ConnectionParameter(
                calc=self._calc_mHA,
                quantity="mass_flow",
                description="mass flow of humid air"
            ),
            "mH2O": ConnectionParameter(
                calc=self._calc_mH2O,
                quantity="mass_flow",
                description="mass flow of liquid or solid water not contained in humid air"
            ),
            "p": ConnectionParameter(
                quantity="pressure",
                description="absolute pressure of the fluid (system variable)"
            ),
            "h": ConnectionParameter(
                quantity="enthalpy",
                description="dry air mass specific enthalpy (system variable)"
            ),
            "w": ConnectionParameter(
                quantity="ratio",
                description="mass of water per mass of dry air"
            ),
            "T": ConnectionParameter(
                func=self.T_func,
                dependents=self.T_dependents,
                num_eq=1,
                quantity="temperature",
                description="temperature of the fluid"
            ),
            "v": ConnectionParameter(
                func=self.v_func,
                dependents=self.v_dependents,
                num_eq=1,
                quantity="volumetric_flow",
                description="volumetric flow of the fluid"
            ),
            "vol": ConnectionParameter(
                calc=self.calc_vol,
                quantity="specific_volume",
                description="specific volume of the fluid"
            ),
            "s": ConnectionParameter(
                calc=self.calc_s,
                quantity="entropy",
                description="specific entropy of the fluid"
            ),
            "r": ConnectionParameter(
                func=self.r_func,
                dependents=self.r_dependents,
                num_eq=1,
                quantity="ratio",
                description="relative humidity"
            ),
            "fluid_balance": OptionalEquation(
                func=self.fluid_balance_func, num_eq_sets=1,
                dependents=self.fluid_balance_dependents,
                description="apply an equation which closes the fluid balance with at least two unknown fluid mass fractions"
            ),
            **self._exergy_results(),
        }

    def _parameter_specification(self, key, value):
        if key == "w" or key == "w0":
            if value is None:
                self.fluid.is_set = set()
            else:
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
            msg = (
                "You cannot change the mixing rule specification for a "
                f"Connection of type {self.__class__.__name__}"
            )
            raise ValueError(msg)

    mixing_rule = property(_get_mixing_rule, _set_mixing_rule)

    def _temperature_hint(self):
        if self.T.is_set:
            return self.T.val_SI
        if not np.isnan(self.T.val0):
            # holds the SI converted temperature guess from the seed pass
            return self.T.val_SI
        return None

    def _apply_temperature_field(self, field, covered, seeded):
        h_sources = []
        if self not in field or not self.h.is_var:
            return h_sources
        reference = self.h._reference_container
        if reference in covered:
            return h_sources
        try:
            w = self.calc_w()
            h = HAPropsSI("H", "P", self.p.val_SI, "T", field[self], "W", w)
        except ValueError:
            return h_sources
        self.h.set_reference_val_SI(h)
        covered.add(reference)
        h_sources.append(reference)
        return h_sources

    def _guess_starting_values(self, units, covered):
        h_sources = []
        if self.h.is_var and not self.good_starting_values:
            reference = self.h._reference_container
            if reference not in covered:
                value = seeded_random(self.label)
                T_rand = 280 + value * (300 - 280)
                h = fp.h_mix_pT(1e5, T_rand, self.fluid_data, self.mixing_rule)
                self.h.set_reference_val_SI(h)
                covered.add(reference)
                h_sources.append(reference)
            if self._precalc_guess_values():
                covered.add(reference)
                if reference not in h_sources:
                    h_sources.append(reference)
        return h_sources

    def _precalc_guess_values(self):
        if not self.h.is_var or self.good_starting_values:
            return False

        if self.T.is_set:
            try:
                w = self.calc_w()
                self.h.set_reference_val_SI(
                    HAPropsSI("H", "P", self.p.val_SI, "T", self.T.val_SI, "W", w)
                )
                return True
            except ValueError:
                pass
        return False

    def _refine_two_phase_guess(self, seeded=False):
        return False

    def _precalc_guess_values_for_references(self):
        """precalculate starting values for specified temperature
        references
        """
        pass

    def _presolve(self):
        self._presolve_determinations = []

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
                self._presolve_determinations.append(
                    {"property": "h", "via": ["T"], "requires": ["p"]}
                )
                if "T" in self._equation_set_lookup.values():
                    presolved_equations += ["T"]

        presolved_equations = [
            key for parameter in presolved_equations
            for key, value in self._equation_set_lookup.items()
            if value == parameter
        ]
        return presolved_equations

    def _property_bounds(self, prop, nw):
        if prop == "p":
            return 101, 99e5

        elif prop == "h":
            # TODO: check minimum temperature how it matches minimum humidity ratio
            hmin = HAPropsSI("H", "T", -50 + 273.15, "P", self.p.val_SI, "R", 1)
            # TODO: where to get reasonable hmax from?!
            hmax = HAPropsSI("H", "T", 300 + 273.15, "P", self.p.val_SI, "R", 0)
            d = self.h._reference_container._d
            delta = min(
                max(abs(self.h.val_SI * d), d) * 5, (hmax - hmin) / 100
            )
            return hmin + delta, hmax - delta

        return None

    def _finalize_starting_values(self, units, covered, seeded, nw):
        super()._finalize_starting_values(units, covered, seeded, nw)
        if self._project_composition_to_r():
            self._precalc_guess_values()

    def _adjust_to_property_limits(self, nw):
        self._project_composition_to_r()

    def _project_composition_to_r(self):
        """Move a supersaturated variable composition onto the specified
        relative humidity.

        Above saturation the humidity ratio is capped, the relative humidity
        is constant 1 and the equation loses its sensitivity to the
        composition.
        """
        if not self.r.is_set:
            return False

        water_alias = list(_get_fluid_alias("water", self.fluid_data))[0]
        air_alias = list(_get_fluid_alias("air", self.fluid_data))[0]
        if not {water_alias, air_alias} <= self.fluid.is_var:
            return False

        T = self._T_for_evaluation()
        w = self.fluid.val[water_alias] / self.fluid.val[air_alias]
        if w <= w_mix_pTrh_humidair(self.p.val_SI, T, 1.0):
            return False

        w = w_mix_pTrh_humidair(self.p.val_SI, T, self.r.val_SI)
        self.fluid.set_reference_val(water_alias, w / (1 + w))
        self.fluid.set_reference_val(air_alias, 1 / (1 + w))
        return True

    @classmethod
    def _result_attributes(cls):
        return ["m", "mHA", "mH2O", "p", "h", "T", "w", "s", "vol", "v", "r"]

    @classmethod
    def _print_attributes(cls):
        return ["m", "mHA", "mH2O", "p", "h", "T", "w", "r"]

    def _T_for_evaluation(self):
        if self.T.is_set:
            return self.T.val_SI
        return self.calc_T()

    def calc_r(self):
        T = self._T_for_evaluation()
        w = w_mix_pT_humidair(self.p.val_SI, T, self.fluid_data)
        try:
            return HAPropsSI("R", "P", self.p.val_SI, "T", T, "W", w)
        except ValueError as e:
            value = str(e).split("value (")[1].split(")")[0]
            return float(value)

    def r_func(self, **kwargs):
        return self.r.val_SI - self.calc_r()

    def r_dependents(self):
        water_alias = _get_fluid_alias("water", self.fluid_data)
        air_alias = _get_fluid_alias("air", self.fluid_data)
        # the humidity ratio relates the water to the air mass fraction, so
        # the equation depends on both
        return {
            "scalars": [self.p, self.h],
            "vectors": [{self.fluid: water_alias | air_alias}]
        }

    def calc_w(self):
        return w_mix_pT_humidair(self.p.val_SI, self.T.val_SI, self.fluid_data)

    def _calc_mHA(self):
        """Mass flow of humid air from the mass flow of dry air."""
        return self.m.val_SI * (1 + self.w.val_SI)

    def _calc_mH2O(self):
        """Mass flow of liquid or solid water not contained in the humid
        air."""
        w_mixture = w_mix_fluid_data(self.fluid_data)
        return self.m.val_SI * (w_mixture - self.w.val_SI)

    def calc_results(self, units, skip_postprocess):
        self.reset_optional_results()
        self.m.set_val0_from_SI(units)
        self.p.set_val0_from_SI(units)
        self.h.set_val0_from_SI(units)
        self.fluid.val0 = self.fluid.val.copy()
        # the temperature guess is one-shot, see Connection.calc_results
        self.T.val0 = np.nan

        if skip_postprocess:
            return True

        self.T.val_SI = self.calc_T()
        self.w.val_SI = self.calc_w()
        self._calc_result_parameters()
        self.v.val_SI = self.vol.val_SI * self.m.val_SI  # Mixture volume flow rate
        self.r.val_SI = self.calc_r()
        # if self.r.val_SI > 1:
        #     self.r.val_SI = np.nan

        for prop in self._result_attributes():
            param = self.get_attr(prop)
            if not param.is_set:
                param.set_val_from_SI(units)

        return True
