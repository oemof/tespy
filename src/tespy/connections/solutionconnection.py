# -*- coding: utf-8
"""Module for the SolutionConnection class.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/connections/solutionconnection.py

SPDX-License-Identifier: MIT
"""
from tespy.tools import logger
from tespy.tools.data_containers import ConnectionParameter
from tespy.tools.fluid_properties.helpers import single_fluid
from tespy.tools.fluid_properties.solutions import SolutionComponentWrapper
from tespy.tools.fluid_properties.solutions import get_solution_pair
from tespy.tools.fluid_properties.solutions import solution_pair_names
from tespy.tools.fluid_properties.solutions import solution_pair_of_component
from tespy.tools.global_vars import FLUID_ALIASES

from .connection import Connection
from .connection import connection_registry


@connection_registry
class SolutionConnection(Connection):
    r"""
    Connection for the streams of a working pair in absorption cycles.

    The properties are evaluated by a registered working pair, e.g.
    :code:`"libr-water"` (LiBr-water, Pátek and Klomfar 2006) or
    :code:`"nh3-water"` (NH3-water, IAPWS G4-01). The working pair is the
    mixing rule of the connection. It is taken from

    - the :code:`mixing_rule` specification,
    - the fluid names, e.g. a fluid vector containing LiBr,
    - the other connections of the same fluid branch,

    in this order, and falls back to :code:`"libr-water"` if none of them
    provides it.

    The mass fraction :code:`xi` of the key component of the pair (LiBr for
    LiBr-water, NH3 for NH3-water) can be set directly as a convenience
    parameter; it maps internally to the fluid vector of both components.
    :code:`xi0` is the corresponding starting value.

    Besides the results of a :class:`.Connection`, a solution connection
    reports the mass fraction :code:`xi`, the vapor mass fraction :code:`x`
    of the two-phase solution, the bubble point temperature
    :code:`T_bubble` and the subcooling :code:`td_bubble`. The same
    parameters (and :code:`T_dew`, :code:`td_dew` for pairs with a volatile
    absorbent) can be specified.

    .. note::

        In a closed solution loop the :class:`.CycleCloser` does not close
        the variable composition of its outlet: specify
        :code:`fluid_balance=True` on the connection leaving the cycle
        closer.

    Example
    -------
    >>> from tespy.components import Sink, Source
    >>> from tespy.connections import SolutionConnection
    >>> from tespy.networks import Network
    >>> nw = Network(iterinfo=False)
    >>> c = SolutionConnection(
    ...     Source("source"), "out1", Sink("sink"), "in1",
    ...     mixing_rule="libr-water"
    ... )
    >>> c.set_attr(xi=0.6, p=1e4, m=1, td_bubble=5)
    >>> nw.add_conns(c)
    >>> nw.solve("design")
    >>> round(c.T_bubble.val_SI - c.T.val_SI, 6)
    5.0
    """

    def get_parameters(self):
        params = super().get_parameters()
        # the analytical derivative of the temperature misses the
        # composition, the dependents cover it
        params["T"].deriv = None
        params["xi"] = ConnectionParameter(
            quantity="ratio",
            description="mass fraction of the key component of the working pair",
            calc=self.calc_xi,
        )
        params["x"] = ConnectionParameter(
            func=self.x_func, dependents=self.x_dependents, num_eq=1,
            quantity="quality",
            description="vapor mass fraction of the two-phase solution",
            calc=self.calc_x,
        )
        params["T_bubble"] = ConnectionParameter(
            quantity="temperature",
            description="bubble point temperature of the solution",
            calc=self.calc_T_bubble,
        )
        params["td_bubble"] = ConnectionParameter(
            func=self.td_bubble_func, dependents=self.td_bubble_dependents,
            num_eq=1, quantity="temperature_difference",
            description="subcooling temperature difference to the bubble point",
            calc=self.calc_td_bubble,
        )
        params["T_dew"] = ConnectionParameter(
            quantity="temperature",
            description="dew point temperature of the solution",
            calc=self.calc_T_dew,
        )
        params["td_dew"] = ConnectionParameter(
            func=self.td_dew_func, dependents=self.td_dew_dependents,
            num_eq=1, quantity="temperature_difference",
            description="superheating temperature difference to the dew point",
            calc=self.calc_td_dew,
        )
        return params

    @property
    def pair(self):
        """The working pair of the connection, None if it is unknown."""
        rule = self.mixing_rule
        if rule is None:
            return None
        return get_solution_pair(rule)

    def _get_mixing_rule(self):
        if self._mixing_rule is not None:
            return self._mixing_rule
        fluids = set(self.fluid.val) | self.fluid.is_set | set(self.fluid.val0)
        for fluid in fluids:
            pair = solution_pair_of_component(fluid)
            if pair is not None:
                return pair.name
        return None

    def _set_mixing_rule(self, value):
        if value is None:
            self._mixing_rule = None
            return
        if value not in solution_pair_names():
            if self.mixing_rule is not None:
                value = self.mixing_rule
            elif value == "ideal-cond":
                # no working pair specified anywhere in the fluid branch
                value = "libr-water"
            else:
                available = "', '".join(solution_pair_names())
                msg = (
                    f"The mixing rule of the {self.__class__.__name__} "
                    f"{self.label} must be a working pair, available are "
                    f"'{available}', '{value}' was given."
                )
                logger.error(msg)
                raise ValueError(msg)
        current = self.mixing_rule
        if current is not None and current != value:
            msg = (
                f"The working pair of the {self.__class__.__name__} "
                f"{self.label} is '{current}', it cannot be changed to "
                f"'{value}'."
            )
            logger.error(msg)
            raise ValueError(msg)
        self._mixing_rule = value
        self._apply_pending_composition()

    mixing_rule = property(_get_mixing_rule, _set_mixing_rule)

    def __init__(self, source, outlet_id, target, inlet_id, label=None, **kwargs):
        self._mixing_rule = None
        self._pending_composition = {}
        self._user_composition0 = False
        super().__init__(source, outlet_id, target, inlet_id, label, **kwargs)

    def _parameter_specification(self, key, value):
        if key in ("xi", "xi0"):
            self._pending_composition[key] = value
            if self.mixing_rule is not None:
                self._apply_pending_composition()
        else:
            super()._parameter_specification(key, value)

    def _apply_pending_composition(self):
        pair = self.pair
        if pair is None:
            return
        for key, value in self._pending_composition.items():
            if key == "xi":
                if value is None:
                    self.fluid.is_set -= set(pair.components)
                else:
                    self.set_attr(fluid=self._composition(pair, value))
            elif value is not None:
                self.set_attr(fluid0=self._composition(pair, value))
                self._user_composition0 = True
        self._pending_composition = {}
        # the absorbent of a pair without volatile absorbent has no pure
        # fluid properties
        if not pair.volatile_absorbent:
            self.fluid.engine[pair.absorbent] = SolutionComponentWrapper
        for component in pair.components:
            known = set(self.fluid.val) | self.fluid.is_set | set(self.fluid.val0)
            if not FLUID_ALIASES.get_fluid(component) & known:
                self.fluid.val0[component] = 0.5

    @staticmethod
    def _composition(pair, w):
        other = pair.refrigerant if pair.key == pair.absorbent else pair.absorbent
        return {pair.key: w, other: 1 - w}

    def _fluid_specification(self, key, value):
        super()._fluid_specification(key, value)
        if key == "fluid0":
            self._user_composition0 = True
        if key == "fluid":
            pair = self.pair
            if pair is not None and not pair.volatile_absorbent:
                self.fluid.engine[pair.absorbent] = SolutionComponentWrapper

    def _key_component(self):
        pair = self.pair
        aliases = FLUID_ALIASES.get_fluid(pair.key) & set(self.fluid.val)
        return next(iter(aliases), pair.key)

    def _composition_dependents(self):
        return {self.fluid: {self._key_component()} & self.fluid.is_var}

    def x_dependents(self):
        return {
            "scalars": [[self.p, self.h]],
            "vectors": [self._composition_dependents()],
        }

    def td_bubble_dependents(self):
        return self.x_dependents()

    def td_dew_dependents(self):
        return self.x_dependents()

    def T_dependents(self):
        return self.x_dependents()

    def calc_xi(self):
        """Return the mass fraction of the key component of the pair."""
        return self.pair.key_fraction(self.fluid_data)

    def _guess_starting_values(self, units, covered):
        if not self.good_starting_values and not self._user_composition0:
            self._guess_equilibrium_composition()
        h_sources = super()._guess_starting_values(units, covered)
        if self.h.is_var and not self.good_starting_values:
            reference = self.h._reference_container
            if reference not in covered:
                h = self._h_guess()
                if h is not None:
                    self.h.set_reference_val_SI(h)
                    covered.add(reference)
                    h_sources.append(reference)
        return h_sources

    def _equilibrium_composition_guess(self):
        """Key component mass fraction of a solution leaving a sorption
        component at its bubble point at the guessed state, also through a
        cycle closer, or None."""
        source = self.source
        if source.__class__.__name__ == "CycleCloser":
            inlet = source.inl[0]
            if isinstance(inlet, SolutionConnection):
                return inlet._equilibrium_composition_guess()
            return None
        solution_outlet = getattr(source, "_solution_outlet", None)
        if solution_outlet is None or solution_outlet() is not self:
            return None
        T = self._temperature_hint()
        if T is None:
            return None
        pair = self.pair
        try:
            return pair.w_bubble(self.p.val_SI, T)
        except ValueError:
            return None

    def _guess_equilibrium_composition(self):
        """Start the variable composition at the equilibrium composition of
        the sorption component the solution leaves."""
        if not self.fluid.is_var:
            return
        w = self._equilibrium_composition_guess()
        if w is None:
            return
        composition = self._composition(self.pair, w)
        for fluid in self.fluid.is_var:
            for component, value in composition.items():
                if fluid in FLUID_ALIASES.get_fluid(component):
                    self.fluid.set_reference_val(fluid, value)

    def _h_guess(self):
        """Enthalpy of the subcooled solution at the connection pressure."""
        pair = self.pair
        w = pair.key_fraction(self.fluid_data)
        try:
            T = pair.T_bubble(self.p.val_SI, w) - 5
            return pair.h_pT(self.p.val_SI, T, w)
        except (ValueError, ZeroDivisionError):
            return None

    def _property_bounds(self, prop, nw):
        if prop == "h" and single_fluid(self.fluid_data) is None:
            pair = self.pair
            w = pair.key_fraction(self.fluid_data)
            p = self.p.val_SI
            try:
                lower = pair.h_pT(p, pair.T_min - 20, w)
                upper = pair.h_pT(p, pair.T_max + 20, w)
            except (ValueError, ZeroDivisionError):
                return super()._property_bounds(prop, nw)
            return lower, upper
        return super()._property_bounds(prop, nw)

    def _adjust_to_property_limits(self, nw):
        if self.h.is_var and single_fluid(self.fluid_data) is None and self.it < 30:
            pair = self.pair
            w = pair.key_fraction(self.fluid_data)
            if self.x.is_set:
                self.h.set_reference_val_SI(pair.h_pQ(self.p.val_SI, self.x.val_SI, w))
            elif self.td_bubble.is_set:
                T = pair.T_bubble(self.p.val_SI, w) - self.td_bubble.val_SI
                self.h.set_reference_val_SI(pair.h_pT(self.p.val_SI, T, w))
            return
        super()._adjust_to_property_limits(nw)

    def calc_results(self, units, skip_postprocess):
        converged = super().calc_results(units, skip_postprocess)
        if not skip_postprocess and single_fluid(self.fluid_data) is None:
            try:
                self.phase.val = self.calc_phase()
            except ValueError:
                self.phase.val = "phase not recognized"
        return converged

    def _precalc_guess_values(self):
        assigned = super()._precalc_guess_values()
        if self.h.is_var and self.td_bubble.is_set and not self.good_starting_values:
            pair = self.pair
            w = pair.key_fraction(self.fluid_data)
            try:
                T = pair.T_bubble(self.p.val_SI, w) - self.td_bubble.val_SI
                self.h.set_reference_val_SI(pair.h_pT(self.p.val_SI, T, w))
                assigned = True
            except ValueError:
                pass
        return assigned

    @classmethod
    def _result_attributes(cls):
        return [
            "m", "p", "h", "T", "v", "s", "vol", "x", "xi", "td_bubble",
            "T_bubble", "td_dew", "T_dew"
        ]

    @classmethod
    def _print_attributes(cls):
        return ["m", "p", "h", "T", "xi", "x", "phase"]
