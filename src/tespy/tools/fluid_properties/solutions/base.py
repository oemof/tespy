# -*- coding: utf-8

"""Module for the working pair interface of absorption solutions.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/tools/fluid_properties/solutions/base.py

SPDX-License-Identifier: MIT
"""
from functools import partial

from scipy.optimize import brentq

from tespy.tools.fluid_properties.mixtures import MIXING_RULES
from tespy.tools.fluid_properties.wrappers.base import FluidPropertyWrapper
from tespy.tools.fluid_properties.wrappers.base import wrapper_registry
from tespy.tools.global_vars import FLUID_ALIASES

_SOLUTION_PAIRS = {}


class SolutionPair:
    r"""
    Thermodynamic properties of a binary working pair of absorption cycles.

    A working pair consists of a refrigerant and an absorbent. All state
    functions take the mass fraction :code:`w` of the key component, e.g.
    LiBr in LiBr-water or NH3 in NH3-water, and SI units. A state given by
    pressure, temperature and overall composition is unique for a binary
    pair, also inside the two-phase region.

    Subclasses set the class attributes and implement the state functions.
    Registered pairs are available as mixing rule of the same name, see
    :func:`register_solution_pair`.

    Attributes
    ----------
    name : str
        Name of the pair, also the name of the mixing rule.
    refrigerant : str
        Name of the refrigerant.
    absorbent : str
        Name of the absorbent.
    key : str
        Component the mass fraction :code:`w` refers to.
    volatile_absorbent : bool
        Whether the absorbent is present in the vapor phase.
    molar_mass : dict
        Molar mass of both components in kg/mol.
    T_min, T_max : float
        Validity range of the temperature in K.
    w_max : float
        Upper limit of the key component mass fraction.
    """

    name = None
    refrigerant = None
    absorbent = None
    key = None
    volatile_absorbent = False
    molar_mass = {}
    T_min = None
    T_max = None
    w_max = 1.0

    @property
    def components(self):
        return (self.refrigerant, self.absorbent)

    def key_fraction(self, fluid_data):
        """Return the mass fraction of the key component in *fluid_data*."""
        aliases = FLUID_ALIASES.get_fluid(self.key) & fluid_data.keys()
        if not aliases:
            return 0.0
        return fluid_data[next(iter(aliases))]["mass_fraction"]

    def p_bubble(self, T, w):
        raise NotImplementedError

    def T_bubble(self, p, w):
        raise NotImplementedError

    def w_bubble(self, p, T):
        raise NotImplementedError

    def p_dew(self, T, w):
        raise NotImplementedError

    def T_dew(self, p, w):
        raise NotImplementedError

    def p_TQ(self, T, Q, w):
        """Pressure of the state at temperature *T* with the vapor mass
        fraction *Q*."""
        if Q <= 0:
            return self.p_bubble(T, w)
        if Q >= 1:
            return self.p_dew(T, w)

        def residual(p):
            return self.Q_ph(p, self.h_pT(p, T, w), w) - Q

        return brentq(
            residual, self.p_dew(T, w), self.p_bubble(T, w),
            xtol=1e-8, rtol=1e-12
        )

    def y_eq(self, p, w):
        """Key component mass fraction of the vapor in equilibrium with a
        liquid of composition *w* at its bubble point."""
        raise NotImplementedError

    def w_eq(self, p, y):
        """Key component mass fraction of the liquid in equilibrium with a
        vapor of composition *y* at its dew point."""
        raise NotImplementedError

    def h_pT(self, p, T, w):
        raise NotImplementedError

    def s_pT(self, p, T, w):
        raise NotImplementedError

    def v_pT(self, p, T, w):
        raise NotImplementedError

    def T_ph(self, p, h, w, T0=None):
        raise NotImplementedError

    def T_ps(self, p, s, w, T0=None):
        raise NotImplementedError

    def h_pQ(self, p, Q, w):
        raise NotImplementedError

    def Q_ph(self, p, h, w):
        raise NotImplementedError

    def phase_ph(self, p, h, w):
        raise NotImplementedError

    def s_ph(self, p, h, w, T0=None):
        return self.s_pT(p, self.T_ph(p, h, w, T0=T0), w)

    def v_ph(self, p, h, w, T0=None):
        return self.v_pT(p, self.T_ph(p, h, w, T0=T0), w)


@wrapper_registry
class SolutionComponentWrapper(FluidPropertyWrapper):
    r"""
    Placeholder wrapper for a component that only exists in a solution.

    Components like LiBr have no pure fluid properties. The wrapper carries
    the molar mass and the temperature range of the working pair, all
    properties are evaluated by the pair through its mixing rule.
    """

    def __init__(self, fluid, back_end=None, **kwargs):
        super().__init__(fluid, back_end)
        for pair in _SOLUTION_PAIRS.values():
            component = next((
                c for c in pair.components
                if fluid in FLUID_ALIASES.get_fluid(c)
            ), None)
            if component is not None and not pair.volatile_absorbent:
                break
        else:
            msg = (
                f"The fluid {fluid} is not the absorbent of a registered "
                "working pair without volatile absorbent."
            )
            raise ValueError(msg)
        self._molar_mass = pair.molar_mass[component]
        self._T_min = pair.T_min
        self._T_max = pair.T_max
        self._p_min = 0
        self._p_max = 1e9


def _h_pT(pair, p, T, fluid_data, **kwargs):
    return pair.h_pT(p, T, pair.key_fraction(fluid_data))


def _s_pT(pair, p, T, fluid_data, **kwargs):
    return pair.s_pT(p, T, pair.key_fraction(fluid_data))


def _v_pT(pair, p, T, fluid_data, **kwargs):
    return pair.v_pT(p, T, pair.key_fraction(fluid_data))


def _T_ph(pair, p, h, fluid_data, T0=None, **kwargs):
    return pair.T_ph(p, h, pair.key_fraction(fluid_data), T0=T0)


def _T_ps(pair, p, s, fluid_data, T0=None, **kwargs):
    return pair.T_ps(p, s, pair.key_fraction(fluid_data), T0=T0)


def _phase_ph(pair, p, h, fluid_data, **kwargs):
    return pair.phase_ph(p, h, pair.key_fraction(fluid_data))


def _Q_ph(pair, p, h, fluid_data, **kwargs):
    return pair.Q_ph(p, h, pair.key_fraction(fluid_data))


def _h_pQ(pair, p, Q, fluid_data, **kwargs):
    return pair.h_pQ(p, Q, pair.key_fraction(fluid_data))


def _T_bubble_p(pair, p, fluid_data, **kwargs):
    return pair.T_bubble(p, pair.key_fraction(fluid_data))


def _T_dew_p(pair, p, fluid_data, **kwargs):
    return pair.T_dew(p, pair.key_fraction(fluid_data))


def _p_bubble_T(pair, T, fluid_data, **kwargs):
    return pair.p_bubble(T, pair.key_fraction(fluid_data))


def _p_dew_T(pair, T, fluid_data, **kwargs):
    return pair.p_dew(T, pair.key_fraction(fluid_data))


def _p_TQ(pair, T, Q, fluid_data, **kwargs):
    return pair.p_TQ(T, Q, pair.key_fraction(fluid_data))


def register_solution_pair(pair):
    r"""
    Register a working pair and make it available as mixing rule.

    Parameters
    ----------
    pair : SolutionPair
        Instance of the working pair, its :code:`name` is the name of the
        mixing rule.
    """
    _SOLUTION_PAIRS[pair.name] = pair
    MIXING_RULES.register(
        pair.name,
        h_pT=partial(_h_pT, pair),
        s_pT=partial(_s_pT, pair),
        v_pT=partial(_v_pT, pair),
        T_ph=partial(_T_ph, pair),
        T_ps=partial(_T_ps, pair),
        phase_ph=partial(_phase_ph, pair),
        Q_ph=partial(_Q_ph, pair),
        h_pQ=partial(_h_pQ, pair),
        T_bubble_p=partial(_T_bubble_p, pair),
        T_dew_p=partial(_T_dew_p, pair),
        p_bubble_T=partial(_p_bubble_T, pair),
        p_dew_T=partial(_p_dew_T, pair),
        p_TQ=partial(_p_TQ, pair),
    )


def get_solution_pair(name):
    r"""
    Return the registered working pair *name*.

    Parameters
    ----------
    name : str
        Name of the working pair, e.g. :code:`"libr-water"`.

    Returns
    -------
    SolutionPair
    """
    if name not in _SOLUTION_PAIRS:
        available = "', '".join(sorted(_SOLUTION_PAIRS))
        msg = (
            f"The working pair '{name}' is not available. Available pairs "
            f"are '{available}'."
        )
        raise KeyError(msg)
    return _SOLUTION_PAIRS[name]


def solution_pair_names():
    """Return the names of all registered working pairs."""
    return list(_SOLUTION_PAIRS)


def solution_pairs_of_component(fluid):
    r"""
    Return the working pairs whose key component is *fluid*.

    Parameters
    ----------
    fluid : str
        Fluid name.

    Returns
    -------
    list
        Registered :code:`SolutionPair` instances.
    """
    return [
        pair for pair in _SOLUTION_PAIRS.values()
        if fluid in FLUID_ALIASES.get_fluid(pair.key)
    ]


def solution_pair_of_component(fluid):
    r"""
    Return the working pair whose key component is *fluid*, if it is unique.

    Parameters
    ----------
    fluid : str
        Fluid name.

    Returns
    -------
    SolutionPair, None
        The pair, :code:`None` if no or more than one registered pair has
        *fluid* as key component.
    """
    pairs = solution_pairs_of_component(fluid)
    if len(pairs) == 1:
        return pairs[0]
    return None
