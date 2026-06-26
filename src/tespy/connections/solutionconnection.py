# -*- coding: utf-8
"""Module for the SolutionConnection class.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/connections/solutionconnection.py

SPDX-License-Identifier: MIT
"""
from tespy.tools.data_containers import FluidProperties as dc_prop
from tespy.tools.fluid_properties.functions import h_mix_pT
from tespy.tools.fluid_properties.mixtures import _get_fluid_alias
from tespy.tools.helpers import seeded_random

from .connection import Connection
from .connection import connection_registry


@connection_registry
class SolutionConnection(Connection):
    """Connection for binary LiBr-water absorption-cycle streams.

    Locks the mixing rule to :code:`"libr_water"`, which evaluates
    thermodynamic properties via CoolProp's :code:`INCOMP::LiBr` backend
    (Patek-Klomfar correlations) as a function of temperature and LiBr mass
    fraction.

    The LiBr mass fraction :code:`xi` can be set directly as a convenience
    parameter; it maps internally to
    :code:`fluid={"INCOMP::LiBr": xi, "H2O": 1 - xi}`.

    Notes
    -----
    - Only the liquid-solution side is modelled here.  The refrigerant
      (steam) leaving or entering the absorber/desorber uses a plain
      :class:`.Connection` with :code:`fluid={"H2O": 1}`.
    - Specifying :code:`fluid={"LiBr": xi, "H2O": 1 - xi}` (without the
      :code:`INCOMP::` prefix) is also accepted; the prefix is injected
      automatically.
    """

    def get_parameters(self):
        params = super().get_parameters()
        params["xi"] = dc_prop(
            quantity="ratio",
            description="LiBr mass fraction in solution (convenience parameter)"
        )
        return params

    def _get_mixing_rule(self):
        return "libr_water"

    def _set_mixing_rule(self, value):
        if value is not None and value != self.mixing_rule:
            msg = (
                "You cannot change the mixing rule specification for a "
                f"Connection of type {self.__class__.__name__}."
            )
            raise ValueError(msg)

    mixing_rule = property(_get_mixing_rule, _set_mixing_rule)

    def _parameter_specification(self, key, value):
        if key in ("xi", "xi0"):
            if value is None:
                self.fluid.is_set = set()
            else:
                fluid_spec = {"INCOMP::LiBr": value, "H2O": 1.0 - value}
                if key == "xi":
                    self.set_attr(fluid=fluid_spec)
                else:
                    self.set_attr(fluid0=fluid_spec)
        else:
            super()._parameter_specification(key, value)

    def _fluid_specification(self, key, value):
        if key == "fluid" and isinstance(value, dict):
            translated = {}
            for fluid_name, fraction in value.items():
                if "::" not in fluid_name and fluid_name.lower() == "libr":
                    fluid_name = f"INCOMP::{fluid_name}"
                translated[fluid_name] = fraction
            value = translated
        super()._fluid_specification(key, value)

    def _presolve(self):
        water_alias = _get_fluid_alias("H2O", self.fluid_data)
        if not water_alias:
            msg = (
                f"H2O must be present in the fluid composition of "
                f"SolutionConnection {self.label!r}."
            )
            raise ValueError(msg)

        if len(self.fluid.is_var) > 0:
            return []

        presolved_equations = []
        if self.h.is_var and not self.p.is_var and self.T.is_set:
            self.h.set_reference_val_SI(
                h_mix_pT(self.p.val_SI, self.T.val_SI, self.fluid_data, self.mixing_rule)
            )
            self.h._potential_var = False
            if "T" in self._equation_set_lookup.values():
                presolved_equations += ["T"]

        return [
            key
            for parameter in presolved_equations
            for key, val in self._equation_set_lookup.items()
            if val == parameter
        ]

    def _precalc_guess_values(self):
        if not self.h.is_var:
            return
        if not self.good_starting_values and self.T.is_set:
            try:
                self.h.set_reference_val_SI(
                    h_mix_pT(self.p.val_SI, self.T.val_SI, self.fluid_data, self.mixing_rule)
                )
            except Exception:
                pass

    def _adjust_to_property_limits(self, nw):
        ref = self.fluid._reference_container
        if ref is not None and "LiBr" in self.fluid.is_var:
            xi = self.fluid.val.get("LiBr", 0)
            if xi > 0.74 or xi < 0.01:
                xi_clipped = min(max(xi, 0.01), 0.74)
                ref.val["LiBr"] = xi_clipped
                ref.val["H2O"] = 1.0 - xi_clipped
        super()._adjust_to_property_limits(nw)

    def _guess_starting_values(self, units):
        super()._guess_starting_values(units)
        if self.h.is_var and not self.good_starting_values:
            rand = seeded_random(self.label)
            T_rand = 310 + rand * (420 - 310)
            try:
                h = h_mix_pT(1e5, T_rand, self.fluid_data, self.mixing_rule)
                self.h.set_reference_val_SI(h)
            except Exception:
                pass
            self._precalc_guess_values()

    def calc_xi(self):
        """Return solved LiBr mass fraction from fluid composition."""
        water = _get_fluid_alias("H2O", self.fluid_data)
        if water:
            return 1.0 - self.fluid.val[next(iter(water))]
        return sum(self.fluid.val.values())

    def calc_results(self, units, skip_postprocess):
        if not skip_postprocess:
            self.xi.val_SI = self.calc_xi()
        return super().calc_results(units, skip_postprocess)

    @classmethod
    def _result_attributes(cls):
        return ["m", "p", "h", "T", "v", "s", "vol", "xi"]

    @classmethod
    def _print_attributes(cls):
        return ["m", "p", "h", "T", "xi"]
