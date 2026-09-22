# -*- coding: utf-8

"""Base class and registry for fluid property wrappers.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/tools/fluid_properties/wrappers/base.py

SPDX-License-Identifier: MIT
"""


def wrapper_registry(type):
    wrapper_registry.items[type.__name__] = type
    return type


wrapper_registry.items = {}


@wrapper_registry
class FluidPropertyWrapper:

    def __init__(self, fluid, back_end=None, **kwargs) -> None:
        """Base class for fluid property wrappers

        Parameters
        ----------
        fluid : str
            Name of the fluid.
        back_end : str, optional
            Name of the back end, by default None
        """
        self.back_end = back_end
        self.fluid = fluid
        self.mixture_type = None

    def _not_implemented(self) -> None:
        raise NotImplementedError(
            f"Method is not implemented for {self.__class__.__name__}."
        )

    def _identify_mixture(self):
        """Parse the fluid name to identify, if and what kind of mixture we are
        working with
        """
        if "[" in self.fluid:
            if "|" not in self.fluid:
                msg = (
                    f"The fluid {self.fluid} requires the specification of "
                    "mass, volume or molar based composition information."
                    "You can do this by appending '|' and 'mass' at the end "
                    "of the fluid string. For example, "
                    "'NAMEOFFLUID[0.5]|mass' to indicate a mass based mixture."
                )
                raise ValueError(msg)

            self.fluid, self.mixture_type = self.fluid.split("|")
            allowed = ["mass", "molar", "volume"]
            if self.mixture_type not in allowed:
                msg = (
                    "For the specification of the composition type you have "
                    f"to select from {', '.join(allowed)}."
                )
                raise ValueError(msg)

        if "&" in self.fluid:
            _fluids_with_fractions = self.fluid.split("&")
        else:
            _fluids_with_fractions = [self.fluid]

        fluid_names = []
        fractions = []
        for fluid in _fluids_with_fractions:
            if "[" in fluid:
                _fluid_name, _fraction = fluid.split("[")
                _fraction = float(_fraction.replace("]", ""))
                fractions += [_fraction]
            else:
                _fluid_name = fluid
            fluid_names += [_fluid_name]

        self.fractions = fractions
        self.fluid = "&".join(fluid_names)

    def isentropic(self, p_1, h_1, p_2):
        self._not_implemented()

    def _is_below_T_critical(self, T):
        self._not_implemented()

    def T_ph(self, p, h):
        self._not_implemented()

    def T_ps(self, p, s):
        self._not_implemented()

    def h_pT(self, p, T):
        self._not_implemented()

    def h_ps(self, p, s):
        self._not_implemented()

    def h_QT(self, Q, T):
        self._not_implemented()

    def h_pQ(self, p, Q):
        self._not_implemented()

    def s_QT(self, Q, T):
        self._not_implemented()

    def T_sat(self, p):
        self._not_implemented()

    def T_dew(self, p):
        return self.T_sat(p)

    def T_bubble(self, p):
        return self.T_sat(p)

    def p_sat(self, T):
        self._not_implemented()

    def p_dew(self, T):
        return self.p_sat(T)

    def p_bubble(self, T):
        return self.p_sat(T)

    def p_sat_TQ(self, T, Q):
        self._not_implemented()

    def Q_ph(self, p, h):
        self._not_implemented()

    def phase_ph(self, p, h):
        self._not_implemented()

    def d_ph(self, p, h):
        self._not_implemented()

    def d_pT(self, p, T):
        self._not_implemented()

    def d_QT(self, Q, T):
        self._not_implemented()

    def viscosity_ph(self, p, h):
        self._not_implemented()

    def viscosity_pT(self, p, T):
        self._not_implemented()

    def conductivity_ph(self, p, h):
        self._not_implemented()

    def conductivity_pT(self, p, T):
        self._not_implemented()

    def s_ph(self, p, h):
        self._not_implemented()

    def s_pT(self, p, T):
        self._not_implemented()
