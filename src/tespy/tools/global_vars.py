# -*- coding: utf-8

"""Module for global variables used by other modules of the tespy package.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/tools/global_vars.py

SPDX-License-Identifier: MIT
"""
import CoolProp as CP

ERR = 1e-6
GAS_CONSTANT_UNI = 8.314462618
GRAVITY = 9.80665


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
