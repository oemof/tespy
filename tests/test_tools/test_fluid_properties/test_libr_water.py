# -*- coding: utf-8

"""Module for testing the LiBr-water working pair.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_tools/test_fluid_properties/test_libr_water.py

SPDX-License-Identifier: MIT
"""
import numpy as np
from CoolProp.CoolProp import PropsSI
from pytest import approx
from pytest import fixture
from pytest import mark

M_H2O = 0.018015268
M_LIBR = 0.08685

# check values of Patek and Klomfar (2006), Table 7: LiBr mole fraction,
# temperature in K, pressure in Pa, molar density in mol/m3, molar enthalpy in
# J/mol and molar entropy in J/(mol K)
PATEK_KLOMFAR_CHECK_VALUES = [
    (0.05, 300, 3025.1805, 54148.9, 1603.9, 7.79057),
    (0.05, 450, 835097.47, 48984.9, 12189, 36.5288),
    (0.10, 300, 2286.4858, 52985.4, 1445.1, 7.66416),
    (0.10, 450, 647702.12, 48550.2, 11555.2, 34.9553),
    (0.30, 350, 2237.3986, 47826.4, 9072.1, 15.3214),
    (0.40, 450, 43075.149, 45941.8, 21024.4, 33.3788),
]


def _mass_fraction(x):
    return x * M_LIBR / (x * M_LIBR + (1 - x) * M_H2O)


def _molar_mass(x):
    return x * M_LIBR + (1 - x) * M_H2O


@fixture(scope="module")
def pair():
    from tespy.tools.fluid_properties.solutions import get_solution_pair
    return get_solution_pair("libr-water")


@mark.parametrize("x, T, p, rho, h, s", PATEK_KLOMFAR_CHECK_VALUES)
def test_check_values(pair, x, T, p, rho, h, s):
    w = _mass_fraction(x)
    M = _molar_mass(x)
    assert pair.p_bubble(T, w) == approx(p, rel=1e-6)
    p_bubble = pair.p_bubble(T, w)
    assert pair.h_pT(p_bubble, T, w) == approx(h / M, rel=1e-4)
    assert pair.s_pT(p_bubble, T, w) == approx(s / M, rel=1e-4)
    assert pair.v_pT(p_bubble, T, w) == approx(1 / (rho * M), rel=1e-4)


@mark.parametrize("T", [280, 320, 360, 400, 440])
def test_pure_water_limit(pair, T):
    w = 1e-12
    p = pair.p_bubble(T, w)
    assert p == approx(PropsSI("P", "T", T, "Q", 0, "water"), rel=1e-4)
    # the formulation brings its own correlations for saturated liquid
    # water, which deviate slightly from IAPWS-95
    h = PropsSI("H", "T", T, "Q", 0, "water")
    s = PropsSI("S", "T", T, "Q", 0, "water")
    assert pair.h_pT(p * 1.01, T, w) == approx(h, abs=500)
    assert pair.s_pT(p * 1.01, T, w) == approx(s, abs=1.5)


STATES = [(T, w) for T in [290, 330, 370, 410, 450] for w in [0.3, 0.5, 0.65]]


@mark.parametrize("T, w", STATES)
def test_bubble_point_round_trips(pair, T, w):
    p = pair.p_bubble(T, w)
    assert pair.T_bubble(p, w) == approx(T, abs=1e-6)
    assert pair.w_bubble(p, T) == approx(w, abs=1e-8)


@mark.parametrize("T, w", STATES)
def test_bubble_point_increases_with_concentration(pair, T, w):
    assert pair.p_bubble(T, w + 0.01) < pair.p_bubble(T, w)


@mark.parametrize("T, w", STATES)
@mark.parametrize("p_factor", [2.0, 0.5])
def test_temperature_inversion_round_trips(pair, T, w, p_factor):
    p = pair.p_bubble(T, w) * p_factor
    h = pair.h_pT(p, T, w)
    s = pair.s_pT(p, T, w)
    assert pair.T_ph(p, h, w) == approx(T, abs=1e-6)
    assert pair.T_ps(p, s, w) == approx(T, abs=1e-6)
    assert pair.s_ph(p, h, w) == approx(s, abs=1e-6)


@mark.parametrize("T, w", STATES)
def test_flash_balances(pair, T, w):
    p = pair.p_bubble(T, w) * 0.5
    h = pair.h_pT(p, T, w)
    q = pair.Q_ph(p, h, w)
    assert 0 < q < 1
    assert pair.phase_ph(p, h, w) == "tp"
    w_liquid = pair.w_bubble(p, T)
    # the vapor is pure water, all LiBr stays in the liquid
    assert (1 - q) * w_liquid == approx(w, rel=1e-8)
    h_vapor = PropsSI("H", "P", p, "T", T, "water")
    h_liquid = pair.h_pT(p, T, w_liquid)
    assert h == approx((1 - q) * h_liquid + q * h_vapor, rel=1e-8)
    assert pair.h_pQ(p, q, w) == approx(h, rel=1e-8)


@mark.parametrize("T, w", STATES)
def test_subcooled_liquid(pair, T, w):
    p = pair.p_bubble(T, w) * 2
    h = pair.h_pT(p, T, w)
    assert pair.Q_ph(p, h, w) == 0
    assert pair.phase_ph(p, h, w) == "l"


@mark.parametrize("T, w", STATES)
def test_enthalpy_continuous_at_bubble_point(pair, T, w):
    p = pair.p_bubble(T, w)
    h_liquid = pair.h_pT(p, T - 1e-6, w)
    h_two_phase = pair.h_pT(p, T + 1e-6, w)
    assert h_two_phase == approx(h_liquid, abs=1)
    assert pair.h_pQ(p, 0, w) == approx(pair.h_pT(p, T, w), abs=1e-3)


@mark.parametrize("T, w", STATES)
def test_compressed_liquid_pump_work(pair, T, w):
    # incompressible liquid: isentropic compression changes h by v dp
    p = pair.p_bubble(T, w)
    dp = 1e5
    v = pair.v_pT(p, T, w)
    assert pair.s_pT(p + dp, T, w) == approx(pair.s_pT(p, T, w), abs=1e-9)
    assert (
        pair.h_pT(p + dp, T, w) - pair.h_pT(p, T, w)
        == approx(v * dp, rel=1e-6)
    )


def test_heat_of_dilution(pair):
    # diluting a 60 % solution with liquid water to 50 % at 40 degC is
    # exothermic by some tens of kJ per kg of the 60 % solution, this heat of
    # dilution is missing in CoolProp's INCOMP::LiBr
    T = 313.15
    p = 1e5
    m_water = 0.6 / 0.5 - 1
    h_water = PropsSI("H", "P", p, "T", T, "water")
    dh = (
        (1 + m_water) * pair.h_pT(p, T, 0.5)
        - pair.h_pT(p, T, 0.6) - m_water * h_water
    )
    assert -60e3 < dh < -25e3


def test_vapor_is_pure_water(pair):
    assert pair.y_eq(5e3, 0.6) == 0
    assert not pair.volatile_absorbent


@mark.parametrize(
    "p, T, w_reference",
    [
        # CSL-0043 absorber and generator outlet (ASHRAE fits of McNeely)
        (872.575113905957, 313.15, 0.5830767457715),
        (9594.998843682066, 403.15, 0.7477190694506),
    ],
)
def test_equilibrium_composition_against_ashrae(pair, p, T, w_reference):
    # the formulations differ by about half a percentage point, CoolProp's
    # INCOMP::LiBr fit of Patek and Klomfar gives 0.5777 at the first state
    assert pair.w_bubble(p, T) == approx(w_reference, abs=1e-2)


def test_mixing_rule_dispatch(pair):
    from tespy.tools.fluid_properties import Q_mix_ph
    from tespy.tools.fluid_properties import T_mix_ph
    from tespy.tools.fluid_properties import T_mix_ps
    from tespy.tools.fluid_properties import h_mix_pQ
    from tespy.tools.fluid_properties import h_mix_pT
    from tespy.tools.fluid_properties import phase_mix_ph
    from tespy.tools.fluid_properties import s_mix_pT
    from tespy.tools.fluid_properties.functions import T_bubble_p
    from tespy.tools.fluid_properties.functions import p_bubble_T
    from tespy.tools.fluid_properties.solutions import SolutionComponentWrapper
    from tespy.tools.fluid_properties.wrappers import CoolPropWrapper

    w = 0.55
    data = {
        "LiBr": {"wrapper": SolutionComponentWrapper("LiBr"), "mass_fraction": w},
        "H2O": {"wrapper": CoolPropWrapper("H2O"), "mass_fraction": 1 - w},
    }
    rule = "libr-water"
    p, T = 5e3, 340
    h = h_mix_pT(p, T, data, rule)
    s = s_mix_pT(p, T, data, rule)
    assert h == approx(pair.h_pT(p, T, w))
    assert s == approx(pair.s_pT(p, T, w))
    assert T_mix_ph(p, h, data, rule) == approx(T, abs=1e-6)
    assert T_mix_ps(p, s, data, rule) == approx(T, abs=1e-6)
    assert T_bubble_p(p, data, rule) == approx(pair.T_bubble(p, w))
    assert p_bubble_T(T, data, rule) == approx(pair.p_bubble(T, w))
    assert Q_mix_ph(p, h, data, rule) == approx(pair.Q_ph(p, h, w))
    assert phase_mix_ph(p, h, data, rule) == pair.phase_ph(p, h, w)
    assert h_mix_pQ(p, 0, data, rule) == approx(pair.h_pQ(p, 0, w))
    assert np.isfinite(h)


def test_component_wrapper_carries_pair_data(pair):
    from tespy.tools.fluid_properties.solutions import SolutionComponentWrapper

    wrapper = SolutionComponentWrapper("LiBr")
    assert wrapper._molar_mass == approx(M_LIBR)
    assert wrapper._T_min == pair.T_min
    assert wrapper._T_max == pair.T_max


@mark.xfail(
    reason="SolutionConnection on the working pair, added in phase 2",
    strict=True,
)
def test_solution_connection_uses_pair(pair):
    from tespy.components import Sink
    from tespy.components import Source
    from tespy.connections import SolutionConnection
    from tespy.tools.fluid_properties import h_mix_pT

    c = SolutionConnection(
        Source("source"), "out1", Sink("sink"), "in1",
        mixing_rule="libr-water"
    )
    c.set_attr(xi=0.55)
    c._create_fluid_wrapper()
    assert c.mixing_rule == "libr-water"
    assert h_mix_pT(5e3, 340, c.fluid_data, c.mixing_rule) == approx(
        pair.h_pT(5e3, 340, 0.55)
    )
