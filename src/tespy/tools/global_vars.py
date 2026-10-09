# -*- coding: utf-8

"""Module for global variables used by other modules of the tespy package.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/tools/global_vars.py

SPDX-License-Identifier: MIT
"""
import CoolProp as CP

ERR = 1e-6
# convergence acceptance thresholds of the solver, the block wise and the
# simultaneous solution share these definitions. The scaled tolerances
# apply to residuals divided by their per equation scale (the response of
# the equation to order one relative changes of its variables) and to
# increments relative to the variable magnitudes. The scaled residual
# tolerance must stay well above the relative noise floor of the mixture
# property inversions (tol_rel of newton_with_kwargs)
RESIDUAL_TOLERANCE = ERR ** 0.5
INCREMENT_TOLERANCE = ERR ** 0.25
SCALED_RESIDUAL_TOLERANCE = 1e-7
SCALED_INCREMENT_TOLERANCE = 1e-7
# relative tolerance of the postprocessing limit checks: a value beyond a
# declared bound by less than this fraction of its scale is numerical noise
# of the converged solution and snapped onto the bound instead of reported
LIMIT_RTOL = 1e-6
GAS_CONSTANT_UNI = 8.314462618
GRAVITY = 9.80665

_display = {"mode": "compact"}


def set_display_mode(mode):
    """Set how tespy objects render in print, repl and debugger output.

    Parameters
    ----------
    mode : str
        :code:`"none"` for plain python reprs, :code:`"compact"` for
        constructor-like one-liners (default), :code:`"extensive"` for
        tabular views.
    """
    if mode not in ("none", "compact", "extensive"):
        msg = "The display mode must be 'none', 'compact' or 'extensive'."
        raise ValueError(msg)
    _display["mode"] = mode


def get_display_mode():
    return _display["mode"]


class FluidAliases:
    # this method should be able to handle different wrappers:
    # if not in CoolProp just take whatever the fluid name is

    def __init__(self):
        self.fluids = {}

    def get_fluid(self, fluid):
        if fluid not in self.fluids:
            try:
                self.fluids[fluid] = set(
                    alias.replace(" ", "")
                    for alias in CP.CoolProp.get_aliases(fluid)
                )
            except (RuntimeError, ValueError):  # RuntimeError: CoolProp < 8, ValueError: CoolProp >= 8
                self.fluids[fluid] = set([fluid])

        return self.fluids[fluid]


FLUID_ALIASES = FluidAliases()


class CombustionGases:

    def __init__(self):
        self.fluids={
            "hydrogen": {"hf": 0, "LHV": None},
            "methane": {"hf": -74.6, "LHV": None},
            "ethane": {"hf": -84.0, "LHV": None},
            "propane":{"hf": -103.8, "LHV": None},
            "butane":{"hf": -125.7, "LHV": None},
            "nDodecane":{"hf": -289.4, "LHV": None},
            "Dichloroethane":{"hf": -129.79, "LHV": None},
            "CO":{"hf": -110.5, "LHV": None},
        }

    def add_fluid(self, fluid, hf= None, LHV=None):
        """Add a new fluid to the possible combustion gases.

        Parameters
        ----------
        fluid : str
            name of the fluid. Must be a valid fluid within the fluid property
            backend.
        hf: float
            specific enthalpy of formation at standard conditions in kJ/mol
        LHV: float
            lower heating value of the fuel in J/kg
        """
        self.fluids[fluid] = {"hf": hf, "LHV": LHV}


COMBUSTION_FLUIDS = CombustionGases()
