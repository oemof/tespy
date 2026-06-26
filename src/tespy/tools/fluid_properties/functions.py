# -*- coding: utf-8

"""Module for fluid property functions.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/tools/fluid_properties/functions.py

SPDX-License-Identifier: MIT
"""
from CoolProp.CoolProp import HAPropsSI

from tespy.tools.global_vars import FLUID_ALIASES
from tespy.tools.logger import logger

from .helpers import get_number_of_fluids
from .helpers import get_pure_fluid
from .helpers import inverse_temperature_mixture
from .mixtures import MIXING_RULES
from .mixtures import T_mix_ph_libr_water
from .mixtures import phase_mix_ph_libr_water
from .mixtures import w_mix_ph_humidair
from .mixtures import w_mix_ps_humidair


def isentropic(p_1, h_1, p_2, fluid_data, mixing_rule=None, T0=None, T0_out=None):
    r"""
    Calculate specific enthalpy at isentropic outlet conditions.

    For a pure fluid delegates directly to the fluid wrapper. For mixtures
    the isentropic enthalpy is obtained by computing the inlet entropy, then
    inverting :func:`T_mix_ps` at the outlet pressure, and finally evaluating
    :func:`h_mix_pT`.

    Parameters
    ----------
    p_1 : float
        Inlet pressure in Pa.
    h_1 : float
        Inlet specific enthalpy in J/kg.
    p_2 : float
        Outlet pressure in Pa.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Mixing rule identifier; ignored for pure fluids.
    T0 : float, optional
        Temperature starting value in K for the inlet entropy inversion.
    T0_out : float, optional
        Temperature starting value in K for the outlet enthalpy inversion.

    Returns
    -------
    float
        Isentropic specific enthalpy at the outlet in J/kg.

    Notes
    -----
    For mixtures:

    .. math::

        h_2^\text{is} = h\!\left(p_2,\; T\!\left(p_2,\; s(p_1, h_1)\right)\right)
    """
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["wrapper"].isentropic(p_1, h_1, p_2)
    else:
        s_1 = s_mix_ph(p_1, h_1, fluid_data, mixing_rule, T0=T0)
        T_2 = T_mix_ps(p_2, s_1, fluid_data, mixing_rule, T0=T0_out)
        return h_mix_pT(p_2, T_2, fluid_data, mixing_rule)


def _exergy_splitting_in_two_phase(h, s, p, pamb, Tamb, fluid_data):
    r"""
    Check for the special two-phase ambient-state case in exergy calculations.

    If the fluid is a real pure fluid in two-phase state at ambient
    temperature, the thermal exergy is zero and only mechanical exergy is
    returned. Returns :code:`None` for all other cases.

    Parameters
    ----------
    h : float
        Specific enthalpy in J/kg.
    s : float
        Specific entropy in J/(kg K).
    p : float
        Pressure in Pa.
    pamb : float
        Ambient pressure in Pa.
    Tamb : float
        Ambient temperature in K.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.

    Returns
    -------
    tuple or None
        :code:`(e_T, e_M)` in J/kg when the special case applies, otherwise
        :code:`None`.
    """
    pure_fluid = get_pure_fluid(fluid_data)
    if pure_fluid["wrapper"].mixture_type is None:
        # real pure fluid
        if pure_fluid["wrapper"].phase_ph(p, h) == "tp":
            if round(pure_fluid["wrapper"].T_ph(p, h), 6) == round(Tamb, 6):
                h0 = h_mix_pT(pamb, Tamb, fluid_data)
                s0 = s_mix_pT(pamb, Tamb, fluid_data)
                ex_mech = (h - h0) - Tamb * (s - s0)
                return 0.0, ex_mech

    return None


def _physical_exergy_at_min_temperature(h, s, pamb, Tamb, fluid_data, fluid):
    r"""
    Fallback for CoolProp domain errors at ambient conditions.

    Used when the fluid's minimum valid temperature exceeds the ambient
    temperature. Evaluates exergy at :math:`T_\text{min} + 10^{-6}` K and
    returns zero mechanical exergy.

    Parameters
    ----------
    h : float
        Specific enthalpy in J/kg.
    s : float
        Specific entropy in J/(kg K).
    pamb : float
        Ambient pressure in Pa.
    Tamb : float
        Ambient temperature in K.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    fluid : dict
        Single-fluid entry from *fluid_data*.

    Returns
    -------
    tuple
        :code:`(e_ph, 0.0)` - total physical exergy and zero mechanical share,
        both in J/kg.
    """
    Tmin = fluid["wrapper"]._T_min
    logger.warning(
        "Physical exergy at ambient state (p0=%.3f, T0=%.2f) for %s "
        "is outside the FluidWrapper temperature limits --> evaluating at "
        f"the fluid's Tmin {Tmin} K as fallback.",
        pamb, Tamb, fluid["wrapper"].fluid
    )
    # stay marginally above the hard limit to avoid further CoolProp errors
    Tmin_eval = Tmin + 1e-6

    h_min = h_mix_pT(pamb, Tmin_eval, fluid_data)
    s_min = s_mix_pT(pamb, Tmin_eval, fluid_data)

    ex_ph = (h - h_min) - Tamb * (s - s_min)
    return ex_ph, 0.0


def calc_physical_exergy(h, s, p, pamb, Tamb, fluid_data, mixing_rule=None, T0=None):
    r"""
    Calculate specific physical exergy.

    Physical exergy is allocated to a thermal and a mechanical share according
    to :cite:`Morosuk2019`.

    Parameters
    ----------
    h : float
        Specific enthalpy in J/kg.
    s : float
        Specific entropy in J/(kg K).
    p : float
        Pressure in Pa.
    pamb : float
        Ambient pressure in Pa.
    Tamb : float
        Ambient temperature in K.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Mixing rule identifier; ignored for pure fluids.
    T0 : float, optional
        Temperature starting value in K for iterative property inversions.

    Returns
    -------
    tuple
        Specific thermal and mechanical exergy
        (:math:`e^\text{T}`, :math:`e^\text{M}`) in J/kg.

    Notes
    -----
    .. math::

        e^\text{T} = \left( h - h \left( p, T_0 \right) \right) -
        T_0 \cdot \left(s - s\left(p, T_0\right)\right)

        e^\text{M}=\left(h\left(p,T_0\right)-h\left(p_0,T_0\right)\right)
        -T_0\cdot\left(s\left(p, T_0\right)-s\left(p_0,T_0\right)\right)

        e^\text{PH} = e^\text{T} + e^\text{M}
    """
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        if pure_fluid["wrapper"]._T_min > Tamb:
            return _physical_exergy_at_min_temperature(
                h, s, pamb, Tamb, fluid_data, pure_fluid
            )
        else:
            ex = _exergy_splitting_in_two_phase(h, s, p, pamb, Tamb, fluid_data)
            if ex is not None:
                return ex[0], ex[1]

    h_T0_p = h_mix_pT(p, Tamb, fluid_data, mixing_rule)
    s_T0_p = s_mix_pT(p, Tamb, fluid_data, mixing_rule)
    ex_therm = (h - h_T0_p) - Tamb * (s - s_T0_p)
    h0 = h_mix_pT(pamb, Tamb, fluid_data, mixing_rule)
    s0 = s_mix_pT(pamb, Tamb, fluid_data, mixing_rule)
    ex_mech = (h_T0_p - h0) - Tamb * (s_T0_p - s0)
    return ex_therm, ex_mech


def T_mix_ph(p, h, fluid_data, mixing_rule=None, T0=None):
    r"""
    Calculate temperature from pressure and specific enthalpy.

    For a pure fluid delegates to the fluid wrapper. For humid air uses
    :func:`CoolProp.CoolProp.HAPropsSI`. For other mixtures inverts
    :func:`h_mix_pT` iteratively via :func:`.helpers.inverse_temperature_mixture`.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    h : float
        Specific enthalpy in J/kg.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Mixing rule identifier; ignored for pure fluids.
    T0 : float, optional
        Temperature starting value in K for the iterative inversion.

    Returns
    -------
    float
        Temperature in K.
    """
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["wrapper"].T_ph(p, h)
    else:
        if mixing_rule == "humidair":
            w = w_mix_ph_humidair(p, h, fluid_data)
            return HAPropsSI("T", "P", p, "H", h, "W", w)
        elif mixing_rule == "libr_water":
            return T_mix_ph_libr_water(p, h, fluid_data, T0)
        else:
            kwargs = {
                "p": p, "target_value": h, "fluid_data": fluid_data, "T0": T0,
                "f": MIXING_RULES.T_ph(mixing_rule)
            }
            return inverse_temperature_mixture(**kwargs)


def dT_mix_pdh(p, h, fluid_data, mixing_rule=None, T0=None):
    r"""
    Calculate the partial derivative of temperature with respect to enthalpy at constant pressure.

    Uses a central finite difference with step :math:`d = 0.1` J/kg.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    h : float
        Specific enthalpy in J/kg.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Mixing rule identifier; ignored for pure fluids.
    T0 : float, optional
        Temperature starting value in K for the iterative inversion.

    Returns
    -------
    float
        :math:`\partial T / \partial h \big|_p` in K kg/J.

    Notes
    -----
    .. math::

        \frac{\partial T}{\partial h}\bigg|_p
        \approx \frac{T(p,\,h+d) - T(p,\,h-d)}{2\,d}, \quad d = 0.1 \text{ J/kg}
    """
    d = 1e-1
    upper = T_mix_ph(p, h + d, fluid_data, mixing_rule=mixing_rule, T0=T0)
    lower = T_mix_ph(p, h - d, fluid_data, mixing_rule=mixing_rule, T0=upper)
    return (upper - lower) / (2 * d)


def dT_mix_dph(p, h, fluid_data, mixing_rule=None, T0=None):
    r"""
    Calculate the partial derivative of temperature with respect to pressure at constant enthalpy.

    Uses a central finite difference with step :math:`d = 0.1` Pa.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    h : float
        Specific enthalpy in J/kg.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Mixing rule identifier; ignored for pure fluids.
    T0 : float, optional
        Temperature starting value in K for the iterative inversion.

    Returns
    -------
    float
        :math:`\partial T / \partial p \big|_h` in K/Pa.

    Notes
    -----
    .. math::

        \frac{\partial T}{\partial p}\bigg|_h
        \approx \frac{T(p+d,\,h) - T(p-d,\,h)}{2\,d}, \quad d = 0.1 \text{ Pa}
    """
    d = 1e-1
    upper = T_mix_ph(p + d, h, fluid_data, mixing_rule=mixing_rule, T0=T0)
    lower = T_mix_ph(p - d, h, fluid_data, mixing_rule=mixing_rule, T0=upper)
    return (upper - lower) / (2 * d)


def dT_mix_ph_dfluid(p, h, fluid, fluid_data, mixing_rule=None, T0=None):
    r"""
    Calculate the partial derivative of temperature with respect to a fluid mass fraction.

    Perturbs the mass fraction of *fluid* by :math:`d = 10^{-5}` in each
    direction and applies a central finite difference.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    h : float
        Specific enthalpy in J/kg.
    fluid : str
        Name of the fluid whose mass fraction is perturbed.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Mixing rule identifier; ignored for pure fluids.
    T0 : float, optional
        Temperature starting value in K for the iterative inversion.

    Returns
    -------
    float
        :math:`\partial T / \partial x_i \big|_{p,h}` in K.

    Notes
    -----
    .. math::

        \frac{\partial T}{\partial x_i}\bigg|_{p,h}
        \approx \frac{T(p,h,x_i+d) - T(p,h,x_i-d)}{2\,d},
        \quad d = 10^{-5}
    """
    d = 1e-5
    fluid_data[fluid]["mass_fraction"] += d
    upper = T_mix_ph(p, h, fluid_data, mixing_rule=mixing_rule, T0=T0)
    fluid_data[fluid]["mass_fraction"] -= 2 * d
    lower = T_mix_ph(p, h, fluid_data, mixing_rule=mixing_rule, T0=upper)
    fluid_data[fluid]["mass_fraction"] += d
    return (upper - lower) / (2 * d)


def h_mix_pT(p, T, fluid_data, mixing_rule=None):
    r"""
    Calculate specific enthalpy from pressure and temperature.

    For a pure fluid delegates to the fluid wrapper. For mixtures dispatches
    to the mixing rule registered under *mixing_rule*.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Mixing rule identifier; ignored for pure fluids.

    Returns
    -------
    float
        Specific enthalpy in J/kg.
    """
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["wrapper"].h_pT(p, T)
    else:
        return MIXING_RULES.h_pT(mixing_rule)(p, T, fluid_data)


def h_mix_pQ(p, Q, fluid_data, mixing_rule=None):
    r"""
    Calculate specific enthalpy from pressure and vapour quality.

    Valid for single-component :code:`fluid_data`, including mixture-backend
    wrappers handled entirely inside the :class:`.FluidPropertyWrapper`.
    Raises :code:`ValueError` for multi-component :code:`fluid_data`.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    Q : float
        Vapour quality (0 = saturated liquid, 1 = saturated vapour).
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Ignored.

    Returns
    -------
    float
        Specific enthalpy in J/kg.

    Raises
    ------
    ValueError
        If :code:`fluid_data` contains more than one fluid.
    """
    return _single_fluid_wrapper(fluid_data).h_pQ(p, Q)


def dh_mix_dpQ(p, Q, fluid_data, mixing_rule=None):
    r"""
    Calculate the partial derivative of saturation enthalpy with respect to pressure.

    Uses a central finite difference with step :math:`d = 0.1` Pa.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    Q : float
        Vapour quality (0 = saturated liquid, 1 = saturated vapour).
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Ignored.

    Returns
    -------
    float
        :math:`\partial h / \partial p \big|_Q` in J/(kg Pa).

    Notes
    -----
    .. math::

        \frac{\partial h}{\partial p}\bigg|_Q
        \approx \frac{h(p+d,\,Q) - h(p-d,\,Q)}{2\,d}, \quad d = 0.1 \text{ Pa}
    """
    d = 0.1
    upper = h_mix_pQ(p + d, Q, fluid_data)
    lower = h_mix_pQ(p - d, Q, fluid_data)
    return (upper - lower) / (2 * d)


def _single_fluid_wrapper(fluid_data):
    r"""
    Return the fluid wrapper when :code:`fluid_data` contains exactly one entry.

    A single entry may represent either a chemically pure fluid or a
    mixture-backend fluid (e.g. a refrigerant blend handled entirely inside
    the :class:`.wrappers.FluidPropertyWrapper`). Both are treated identically
    by TESPy: no mixing rules are applied and all property calls are
    delegated directly to the wrapper.

    Parameters
    ----------
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.

    Returns
    -------
    FluidPropertyWrapper
        The wrapper of the single fluid entry.

    Raises
    ------
    ValueError
        If :code:`fluid_data` contains more than one fluid.
    """
    if get_number_of_fluids(fluid_data) != 1:
        raise ValueError(
            "This function is not available for multi-component fluid data."
        )
    return get_pure_fluid(fluid_data)["wrapper"]


def Q_mix_ph(p, h, fluid_data, mixing_rule=None):
    r"""
    Calculate vapour quality from pressure and specific enthalpy.

    Valid for single-component :code:`fluid_data`, including mixture-backend
    wrappers handled entirely inside the :class:`.FluidPropertyWrapper`.
    Raises :code:`ValueError` for multi-component :code:`fluid_data`.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    h : float
        Specific enthalpy in J/kg.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Ignored.

    Returns
    -------
    float
        Vapour quality (0 = saturated liquid, 1 = saturated vapour).

    Raises
    ------
    ValueError
        If :code:`fluid_data` contains more than one fluid.
    """
    return _single_fluid_wrapper(fluid_data).Q_ph(p, h)


_MIXING_RULE_PHASE = {
    "ideal": "g",
    "ideal-cond": "g",
    "forced-gas": "g",
    "humidair": "g",
    "incompressible": "l",
}


def phase_mix_ph(p, h, fluid_data, mixing_rule=None):
    r"""
    Return the thermodynamic phase state from pressure and specific enthalpy.

    For single-component :code:`fluid_data` (including mixture-backend
    wrappers) the phase is determined by the underlying
    :class:`.FluidPropertyWrapper`. For multi-component mixtures the phase is
    derived from the :code:`mixing_rule`: gas-type rules (:code:`"ideal"`,
    :code:`"ideal-cond"`, :code:`"forced-gas"`, :code:`"humidair"`) return
    :code:`"g"`; :code:`"incompressible"` returns :code:`"l"`.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    h : float
        Specific enthalpy in J/kg.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Mixing rule identifier; required for multi-component :code:`fluid_data`.

    Returns
    -------
    str
        Phase identifier (e.g. :code:`"tp"` for two-phase, :code:`"g"` for
        gas, :code:`"l"` for liquid).

    Raises
    ------
    ValueError
        If :code:`fluid_data` contains more than one fluid and
        :code:`mixing_rule` is not recognised.
    """
    if get_number_of_fluids(fluid_data) != 1:
        if mixing_rule == "libr_water":
            return phase_mix_ph_libr_water(p, h, fluid_data)
        if mixing_rule not in _MIXING_RULE_PHASE:
            raise ValueError(
                f"Cannot determine phase for multi-component fluid data with "
                f"mixing_rule={mixing_rule!r}. Known rules: "
                + ", ".join(_MIXING_RULE_PHASE)
            )
        return _MIXING_RULE_PHASE[mixing_rule]
    return get_pure_fluid(fluid_data)["wrapper"].phase_ph(p, h)


def p_sat_T(T, fluid_data, mixing_rule=None):
    r"""
    Calculate saturation pressure from temperature.

    Valid for single-component :code:`fluid_data`, including mixture-backend
    wrappers handled entirely inside the :class:`.FluidPropertyWrapper`.
    Raises :code:`ValueError` for multi-component :code:`fluid_data`.

    Parameters
    ----------
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Ignored.

    Returns
    -------
    float
        Saturation pressure in Pa.

    Raises
    ------
    ValueError
        If :code:`fluid_data` contains more than one fluid.
    """
    return _single_fluid_wrapper(fluid_data).p_sat(T)


def p_sat_TQ(T, Q, fluid_data, mixing_rule=None):
    r"""
    Calculate saturation pressure from temperature and vapor quality.

    For pure fluids this is identical to :func:`p_sat_T` for any :code:`Q`.
    For zeotropic mixture backends the pressure corresponds to the two-phase
    state with quality :code:`Q` at temperature :code:`T` (bubble point at
    :code:`Q=0`, dew point at :code:`Q=1`).

    Parameters
    ----------
    T : float
        Temperature in K.
    Q : float
        Vapor quality (0 = bubble point, 1 = dew point).
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Ignored.

    Returns
    -------
    float
        Saturation pressure in Pa.

    Raises
    ------
    ValueError
        If :code:`fluid_data` contains more than one fluid.
    """
    return _single_fluid_wrapper(fluid_data).p_sat_TQ(T, Q)


def T_sat_p(p, fluid_data, mixing_rule=None):
    r"""
    Calculate saturation temperature from pressure.

    Valid for single-component :code:`fluid_data`, including mixture-backend
    wrappers handled entirely inside the :class:`.FluidPropertyWrapper`.
    Raises :code:`ValueError` for multi-component :code:`fluid_data`.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Ignored.

    Returns
    -------
    float
        Saturation temperature in K.

    Raises
    ------
    ValueError
        If :code:`fluid_data` contains more than one fluid.
    """
    return _single_fluid_wrapper(fluid_data).T_sat(p)


def T_dew_p(p, fluid_data, mixing_rule=None):
    r"""
    Calculate dew point temperature from pressure.

    Valid for single-component :code:`fluid_data`, including mixture-backend
    wrappers handled entirely inside the :class:`.FluidPropertyWrapper`.
    Raises :code:`ValueError` for multi-component :code:`fluid_data`.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Ignored.

    Returns
    -------
    float
        Dew point temperature in K.

    Raises
    ------
    ValueError
        If :code:`fluid_data` contains more than one fluid.
    """
    return _single_fluid_wrapper(fluid_data).T_dew(p)


def p_dew_T(T, fluid_data, mixing_rule=None):
    r"""
    Calculate dew point pressure from temperature.

    Valid for single-component :code:`fluid_data`, including mixture-backend
    wrappers handled entirely inside the :class:`.FluidPropertyWrapper`.
    Raises :code:`ValueError` for multi-component :code:`fluid_data`.

    Parameters
    ----------
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Ignored.

    Returns
    -------
    float
        Dew point pressure in Pa.

    Raises
    ------
    ValueError
        If :code:`fluid_data` contains more than one fluid.
    """
    return _single_fluid_wrapper(fluid_data).p_dew(T)


def T_bubble_p(p, fluid_data, mixing_rule=None):
    r"""
    Calculate bubble point temperature from pressure.

    Valid for single-component :code:`fluid_data`, including mixture-backend
    wrappers handled entirely inside the :class:`.FluidPropertyWrapper`.
    Raises :code:`ValueError` for multi-component :code:`fluid_data`.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Ignored.

    Returns
    -------
    float
        Bubble point temperature in K.

    Raises
    ------
    ValueError
        If :code:`fluid_data` contains more than one fluid.
    """
    return _single_fluid_wrapper(fluid_data).T_bubble(p)


def p_bubble_T(T, fluid_data, mixing_rule=None):
    r"""
    Calculate bubble point pressure from temperature.

    Valid for single-component :code:`fluid_data`, including mixture-backend
    wrappers handled entirely inside the :class:`.FluidPropertyWrapper`.
    Raises :code:`ValueError` for multi-component :code:`fluid_data`.

    Parameters
    ----------
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Ignored.

    Returns
    -------
    float
        Bubble point pressure in Pa.

    Raises
    ------
    ValueError
        If :code:`fluid_data` contains more than one fluid.
    """
    return _single_fluid_wrapper(fluid_data).p_bubble(T)


def dT_sat_dp(p, fluid_data, mixing_rule=None):
    r"""
    Calculate the derivative of saturation temperature with respect to pressure.

    Uses a central finite difference with step :math:`d = 0.01` Pa.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Ignored.

    Returns
    -------
    float
        :math:`dT_\text{sat} / dp` in K/Pa.

    Notes
    -----
    .. math::

        \frac{dT_\text{sat}}{dp}
        \approx \frac{T_\text{sat}(p+d) - T_\text{sat}(p-d)}{2\,d},
        \quad d = 0.01 \text{ Pa}
    """
    d = 0.01
    upper = T_sat_p(p + d, fluid_data)
    lower = T_sat_p(p - d, fluid_data)
    return (upper - lower) / (2 * d)


def s_mix_ph(p, h, fluid_data, mixing_rule=None, T0=None):
    r"""
    Calculate specific entropy from pressure and specific enthalpy.

    For a pure fluid delegates to the fluid wrapper. For mixtures first
    determines the temperature via :func:`T_mix_ph`, then evaluates
    :func:`s_mix_pT`.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    h : float
        Specific enthalpy in J/kg.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Mixing rule identifier; ignored for pure fluids.
    T0 : float, optional
        Temperature starting value in K for the iterative inversion.

    Returns
    -------
    float
        Specific entropy in J/(kg K).
    """
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["wrapper"].s_ph(p, h)
    else:
        T = T_mix_ph(p, h, fluid_data, mixing_rule, T0)
        return s_mix_pT(p, T, fluid_data, mixing_rule)


def s_mix_pT(p, T, fluid_data, mixing_rule=None):
    r"""
    Calculate specific entropy from pressure and temperature.

    For a pure fluid delegates to the fluid wrapper. For mixtures dispatches
    to the mixing rule registered under *mixing_rule*.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Mixing rule identifier; ignored for pure fluids.

    Returns
    -------
    float
        Specific entropy in J/(kg K).
    """
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["wrapper"].s_pT(p, T)
    else:
        return MIXING_RULES.s_pT(mixing_rule)(p, T, fluid_data)


def T_mix_ps(p, s, fluid_data, mixing_rule=None, T0=None):
    r"""
    Calculate temperature from pressure and specific entropy.

    For a pure fluid delegates to the fluid wrapper. For humid air uses
    :func:`CoolProp.CoolProp.HAPropsSI`. For other mixtures inverts
    :func:`s_mix_pT` iteratively via :func:`.helpers.inverse_temperature_mixture`.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    s : float
        Specific entropy in J/(kg K).
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Mixing rule identifier; ignored for pure fluids.
    T0 : float, optional
        Temperature starting value in K for the iterative inversion.

    Returns
    -------
    float
        Temperature in K.
    """
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["wrapper"].T_ps(p, s)
    else:
        if mixing_rule == "humidair":
            w = w_mix_ps_humidair(p, s, fluid_data)
            return HAPropsSI("T", "P", p, "S", s, "W", w)
        else:
            kwargs = {
                "p": p, "target_value": s, "fluid_data": fluid_data, "T0": T0,
                "f": MIXING_RULES.T_ps(mixing_rule)
            }
            return inverse_temperature_mixture(**kwargs)


def v_mix_ph(p, h, fluid_data, mixing_rule=None, T0=None):
    r"""
    Calculate specific volume from pressure and specific enthalpy.

    For a pure fluid returns the reciprocal of the wrapper density. For humid
    air uses :func:`CoolProp.CoolProp.HAPropsSI`. For other mixtures first
    determines the temperature via :func:`T_mix_ph`, then evaluates
    :func:`v_mix_pT`.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    h : float
        Specific enthalpy in J/kg.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Mixing rule identifier; ignored for pure fluids.
    T0 : float, optional
        Temperature starting value in K for the iterative inversion.

    Returns
    -------
    float
        Specific volume in m³/kg.
    """
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return 1 / pure_fluid["wrapper"].d_ph(p, h)
    else:
        if mixing_rule == "humidair":
            w = w_mix_ph_humidair(p, h, fluid_data)
            return HAPropsSI("V", "P", p, "H", h, "W", w)
        else:
            T = T_mix_ph(p, h, fluid_data, mixing_rule, T0)
            return v_mix_pT(p, T, fluid_data, mixing_rule)


def dv_mix_dph(p, h, fluid_data, mixing_rule=None, T0=None):
    r"""
    Calculate the partial derivative of specific volume with respect to pressure at constant enthalpy.

    Uses a central finite difference with step :math:`d = 0.1` Pa.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    h : float
        Specific enthalpy in J/kg.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Mixing rule identifier; ignored for pure fluids.
    T0 : float, optional
        Temperature starting value in K for the iterative inversion.

    Returns
    -------
    float
        :math:`\partial v / \partial p \big|_h` in m³/(kg Pa).

    Notes
    -----
    .. math::

        \frac{\partial v}{\partial p}\bigg|_h
        \approx \frac{v(p+d,\,h) - v(p-d,\,h)}{2\,d}, \quad d = 0.1 \text{ Pa}
    """
    d = 1e-1
    upper = v_mix_ph(p + d, h, fluid_data, mixing_rule=mixing_rule, T0=T0)
    lower = v_mix_ph(p - d, h, fluid_data, mixing_rule=mixing_rule, T0=upper)
    return (upper - lower) / (2 * d)


def dv_mix_pdh(p, h, fluid_data, mixing_rule=None, T0=None):
    r"""
    Calculate the partial derivative of specific volume with respect to enthalpy at constant pressure.

    Uses a central finite difference with step :math:`d = 0.1` J/kg.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    h : float
        Specific enthalpy in J/kg.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Mixing rule identifier; ignored for pure fluids.
    T0 : float, optional
        Temperature starting value in K for the iterative inversion.

    Returns
    -------
    float
        :math:`\partial v / \partial h \big|_p` in m³ kg/J.

    Notes
    -----
    .. math::

        \frac{\partial v}{\partial h}\bigg|_p
        \approx \frac{v(p,\,h+d) - v(p,\,h-d)}{2\,d}, \quad d = 0.1 \text{ J/kg}
    """
    d = 1e-1
    upper = v_mix_ph(p, h + d, fluid_data, mixing_rule=mixing_rule, T0=T0)
    lower = v_mix_ph(p, h - d, fluid_data, mixing_rule=mixing_rule, T0=upper)
    return (upper - lower) / (2 * d)


def v_mix_pT(p, T, fluid_data, mixing_rule=None):
    r"""
    Calculate specific volume from pressure and temperature.

    For a pure fluid returns the reciprocal of the wrapper density. For
    mixtures dispatches to the mixing rule registered under *mixing_rule*.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Mixing rule identifier; ignored for pure fluids.

    Returns
    -------
    float
        Specific volume in m³/kg.
    """
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return 1 / pure_fluid["wrapper"].d_pT(p, T)
    else:
        return MIXING_RULES.v_pT(mixing_rule)(p, T, fluid_data)


def viscosity_mix_ph(p, h, fluid_data, mixing_rule=None, T0=None):
    r"""
    Calculate dynamic viscosity from pressure and specific enthalpy.

    For a pure fluid delegates to the fluid wrapper. For humid air uses
    :func:`CoolProp.CoolProp.HAPropsSI`. For other mixtures first determines
    the temperature via :func:`T_mix_ph`, then evaluates
    :func:`viscosity_mix_pT`.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    h : float
        Specific enthalpy in J/kg.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Mixing rule identifier; ignored for pure fluids.
    T0 : float, optional
        Temperature starting value in K for the iterative inversion.

    Returns
    -------
    float
        Dynamic viscosity in Pa s.
    """
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["wrapper"].viscosity_ph(p, h)
    else:
        if mixing_rule == "humidair":
            w = w_mix_ph_humidair(p, h, fluid_data)
            return HAPropsSI("Visc", "P", p, "H", h, "W", w)
        else:
            T = T_mix_ph(p, h, fluid_data, mixing_rule, T0)
            return viscosity_mix_pT(p, T, fluid_data, mixing_rule)


def viscosity_mix_pT(p, T, fluid_data, mixing_rule=None):
    r"""
    Calculate dynamic viscosity from pressure and temperature.

    For a pure fluid delegates to the fluid wrapper. For mixtures dispatches
    to the mixing rule registered under *mixing_rule*.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    T : float
        Temperature in K.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Mixing rule identifier; ignored for pure fluids.

    Returns
    -------
    float
        Dynamic viscosity in Pa s.
    """
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["wrapper"].viscosity_pT(p, T)
    else:
        return MIXING_RULES.viscosity_pT(mixing_rule)(p, T, fluid_data)


def conductivity_mix_ph(p, h, fluid_data, mixing_rule=None, T0=None):
    r"""
    Calculate thermal conductivity from pressure and specific enthalpy.

    Only implemented for single-component :code:`fluid_data`. Raises
    :code:`NotImplementedError` for multi-component :code:`fluid_data`.

    Parameters
    ----------
    p : float
        Pressure in Pa.
    h : float
        Specific enthalpy in J/kg.
    fluid_data : dict
        Fluid property data:
        :code:`{fluid: {"mass_fraction": float, "wrapper": FluidPropertyWrapper}}`.
    mixing_rule : str, optional
        Ignored.
    T0 : float, optional
        Ignored.

    Returns
    -------
    float
        Thermal conductivity in W/(m K).

    Raises
    ------
    NotImplementedError
        If :code:`fluid_data` contains more than one fluid.
    """
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["wrapper"].conductivity_ph(p, h)
    else:
        msg = (
            "Calculation of thermal conductivity is not implemented for "
            "TESPy based mixtures. You are happily invited to contribute it!"
        )
        raise NotImplementedError(msg)
