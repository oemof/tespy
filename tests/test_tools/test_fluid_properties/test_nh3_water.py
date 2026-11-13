# -*- coding: utf-8

"""Module for testing the NH3-water working pair.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_tools/test_fluid_properties/test_nh3_water.py

SPDX-License-Identifier: MIT
"""
import numpy as np
from pytest import approx
from pytest import fixture
from pytest import importorskip
from pytest import mark
from pytest import raises

importorskip("teqp")

# reference: CoolSolve Library CSL-0249 (nh3h2o_absorption_chiller.sol),
# IAPWS G4-01 (Tillner-Roth and Friend 1998) as implemented in CoolSolve;
# per state: pressure in Pa, temperature in degC, NH3 mass fraction, specific
# enthalpy in J/kg, specific entropy in J/(kg K), vapor quality (nan if the
# state is not two-phase)
P_LOW = 2.905681496458e+05
P_HIGH = 1.554537370181e+06
CSL_0249_STATES = {
    1: (P_LOW, 35.0, 4.172235440499e-01, 5.377746443431e+04, 9.334761143742e+02, 0),
    2: (P_HIGH, 3.532424076706e+01, 4.172235440499e-01, 5.627376225467e+04, 9.367142262533e+02, np.nan),
    3: (P_HIGH, 8.491217025590e+01, 4.172235440499e-01, 2.875468603829e+05, 1.631549908719e+03, np.nan),
    4: (P_HIGH, 120.0, 3.067414564222e-01, 4.370375413989e+05, 1.930351202155e+03, 0),
    5: (P_HIGH, 6.072696853694e+01, 3.067414564222e-01, 1.618823071612e+05, 1.172345733674e+03, np.nan),
    6: (P_LOW, 5.706949132980e+01, 3.067414564222e-01, 1.618823071612e+05, 1.177161890492e+03, 1.054540031105e-02),
    7: (P_HIGH, 120.0, 9.034029548463e-01, 1.930343094954e+06, 6.393963327270e+03, 1),
    8: (P_HIGH, 4.612317847674e+01, 0.9995, 1.654229789422e+06, 5.694425556381e+03, 1),
    10: (P_HIGH, 40.0, 0.9995, 5.331945717912e+05, 2.116309358482e+03, 0),
    11: (P_HIGH, 5.839901118200e+00, 0.9995, 3.702955528523e+05, 1.565808892588e+03, np.nan),
    12: (P_LOW, -9.999257775245e+00, 0.9995, 3.702955528523e+05, 1.581201933578e+03, 5.685365084651e-02),
    13: (P_LOW, -8.800141259714e+00, 0.9995, 1.583403098026e+06, 6.190295436686e+03, 0.99),
    14: (P_LOW, 5.475790428811e+01, 0.9995, 1.746302116964e+06, 6.748182629827e+03, np.nan),
}


@fixture(scope="module")
def pair():
    from tespy.tools.fluid_properties.solutions import get_solution_pair
    return get_solution_pair("nh3-water")


def test_pressure_levels(pair):
    # bubble point of the refrigerant at the evaporator and condenser
    assert pair.p_bubble(263.15, 0.9995) == approx(P_LOW, rel=1e-6)
    assert pair.p_bubble(313.15, 0.9995) == approx(P_HIGH, rel=1e-6)
    assert pair.T_bubble(P_LOW, 0.9995) == approx(263.15, abs=1e-5)


@mark.parametrize("state", [1, 4])
def test_equilibrium_liquid_composition(pair, state):
    p, T, w, *_ = CSL_0249_STATES[state]
    assert pair.w_bubble(p, T + 273.15) == approx(w, abs=1e-7)


def test_equilibrium_vapor_composition(pair):
    p, _, w, *_ = CSL_0249_STATES[4]
    assert pair.y_eq(p, w) == approx(CSL_0249_STATES[7][2], abs=1e-7)
    assert pair.volatile_absorbent


def test_dew_point(pair):
    p, T, y, *_ = CSL_0249_STATES[8]
    assert pair.T_dew(p, y) == approx(T + 273.15, abs=1e-5)
    p, T, y, *_ = CSL_0249_STATES[7]
    assert pair.T_dew(p, y) == approx(T + 273.15, abs=1e-5)


@mark.parametrize("state", [1, 4, 5, 10, 11])
def test_liquid_states_from_temperature(pair, state):
    p, T, w, h, s, _ = CSL_0249_STATES[state]
    assert pair.h_pT(p, T + 273.15, w) == approx(h, abs=1)
    assert pair.s_pT(p, T + 273.15, w) == approx(s, abs=1e-2)


@mark.parametrize("state", [7, 8])
def test_saturated_vapor_states(pair, state):
    p, T, w, h, s, _ = CSL_0249_STATES[state]
    assert pair.h_pQ(p, 1, w) == approx(h, abs=1)
    assert pair.T_ph(p, h, w) == approx(T + 273.15, abs=1e-5)
    assert pair.s_ph(p, h, w) == approx(s, abs=1e-2)


@mark.parametrize("state", [2, 3, 6, 12, 13, 14])
def test_states_from_enthalpy(pair, state):
    p, T, w, h, s, q = CSL_0249_STATES[state]
    assert pair.T_ph(p, h, w) == approx(T + 273.15, abs=1e-5)
    assert pair.s_ph(p, h, w) == approx(s, abs=1e-2)
    if np.isfinite(q):
        assert pair.Q_ph(p, h, w) == approx(q, abs=1e-7)


def test_evaporator_outlet_from_quality(pair):
    p, T, w, h, _, q = CSL_0249_STATES[13]
    assert pair.h_pQ(p, q, w) == approx(h, abs=1)


@mark.parametrize("state", [3, 5, 14])
def test_single_phase_states(pair, state):
    p, _, w, h, *_ = CSL_0249_STATES[state]
    expected = "l" if state in (3, 5) else "g"
    assert pair.phase_ph(p, h, w) == expected


@mark.parametrize("state", [6, 12, 13])
def test_two_phase_states(pair, state):
    p, _, w, h, *_ = CSL_0249_STATES[state]
    assert pair.phase_ph(p, h, w) == "tp"


@mark.parametrize("p", [2e5, 5e5, 10e5, 20e5])
@mark.parametrize("w", [0.1, 0.3, 0.5, 0.7, 0.9, 0.99])
def test_flash_robustness(pair, p, w):
    T_bubble = pair.T_bubble(p, w)
    T_dew = pair.T_dew(p, w)
    assert T_dew > T_bubble
    for T in np.linspace(T_bubble - 30, T_dew + 30, 12):
        h = pair.h_pT(p, T, w)
        assert pair.T_ph(p, h, w) == approx(T, abs=1e-5)
        s = pair.s_pT(p, T, w)
        assert pair.T_ps(p, s, w) == approx(T, abs=1e-5)
        q = pair.Q_ph(p, h, w)
        assert 0 <= q <= 1


def test_mixing_rule_dispatch(pair):
    from tespy.tools.fluid_properties import T_mix_ph
    from tespy.tools.fluid_properties import h_mix_pT
    from tespy.tools.fluid_properties.functions import T_dew_p
    from tespy.tools.fluid_properties.wrappers import CoolPropWrapper

    w = 0.4
    data = {
        "NH3": {"wrapper": CoolPropWrapper("NH3"), "mass_fraction": w},
        "H2O": {"wrapper": CoolPropWrapper("H2O"), "mass_fraction": 1 - w},
    }
    p, T = 5e5, 360
    h = h_mix_pT(p, T, data, "nh3-water")
    assert h == approx(pair.h_pT(p, T, w))
    assert T_mix_ph(p, h, data, "nh3-water") == approx(T, abs=1e-6)
    assert T_dew_p(p, data, "nh3-water") == approx(pair.T_dew(p, w))


def test_missing_teqp_names_the_package(monkeypatch):
    import sys

    from tespy.tools.fluid_properties.solutions.nh3_water import AmmoniaWater

    monkeypatch.setitem(sys.modules, "teqp", None)
    with raises(ImportError, match="teqp"):
        AmmoniaWater().p_bubble(300, 0.5)


def _root_with_failures_above(x):
    if x > 2.5:
        raise ValueError
    return x - 2


@mark.parametrize("x0", [0, 2.7, 10])
def test_bracketing_refines_towards_failing_region(x0):
    from tespy.tools.fluid_properties.solutions.nh3_water import _increasing_root

    assert _increasing_root(_root_with_failures_above, x0, 1, 1e-12) == approx(2)


def test_bracketing_fails_without_valid_points():
    from tespy.tools.fluid_properties.solutions.nh3_water import _increasing_root

    def always_failing(x):
        raise ValueError

    with raises(ValueError, match="bracketed"):
        _increasing_root(always_failing, 0, 1, 1e-12)
