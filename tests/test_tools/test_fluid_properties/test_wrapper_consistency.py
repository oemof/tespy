# -*- coding: utf-8

"""Shared consistency tests for all fluid property wrappers.

Every wrapper is registered as a case with its valid test states and
capability flags, all roundtrip and consistency tests run against each
case. To cover a new wrapper, add a case to WRAPPER_CASES.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_tools/test_fluid_properties/test_wrapper_consistency.py

SPDX-License-Identifier: MIT
"""
import importlib.util
from dataclasses import dataclass
from dataclasses import field

import pytest
from CoolProp.CoolProp import get_global_param_string

from tespy.tools.fluid_properties.wrappers import CoolPropWrapper
from tespy.tools.fluid_properties.wrappers import FluidPropertyWrapper
from tespy.tools.fluid_properties.wrappers import wrapper_registry


def _module_missing(name):
    return importlib.util.find_spec(name) is None


skipif_no_refprop = pytest.mark.skipif(
    get_global_param_string("REFPROP_version") == "n/a",
    reason="This test requires REFPROP.",
)
skipif_no_thermopack = pytest.mark.skipif(
    _module_missing("thermopack"), reason="This test requires thermopack."
)
skipif_no_iapws = pytest.mark.skipif(
    _module_missing("iapws"), reason="This test requires iapws."
)
skipif_no_pyromat = pytest.mark.skipif(
    _module_missing("pyromat"), reason="This test requires pyromat."
)


def _thermopack(fluid, back_end):
    from tespy.tools.fluid_properties.wrappers import ThermopackWrapper
    return ThermopackWrapper(fluid, back_end)


def _iapws(back_end):
    from tespy.tools.fluid_properties.wrappers import IAPWSWrapper
    return IAPWSWrapper("H2O", back_end)


def _pyromat(fluid, back_end):
    from tespy.tools.fluid_properties.wrappers import PyromatWrapper
    return PyromatWrapper(fluid, back_end)


@dataclass
class Case:
    """A wrapper under test with valid states and capability flags."""
    factory: callable
    p: float  # pressure for gas state tests (Pa)
    T: float  # gas state temperature at p (K)
    T_sat: float  # temperature on the saturation curve (K)
    two_phase: bool = True
    mixture: bool = False
    rtol: float = 1e-6
    marks: tuple = ()


WRAPPER_CASES = {
    "coolprop-heos-water": Case(
        factory=lambda: CoolPropWrapper("Water"),
        p=1e5, T=500, T_sat=350,
    ),
    "coolprop-heos-propane": Case(
        factory=lambda: CoolPropWrapper("Propane"),
        p=5e5, T=320, T_sat=280,
    ),
    "coolprop-refprop-mixture": Case(
        factory=lambda: CoolPropWrapper(
            "Propane[0.5]&Isobutane[0.5]|molar", "REFPROP"
        ),
        p=5e5, T=330, T_sat=280, mixture=True,
        marks=(skipif_no_refprop,),
    ),
    "iapws-if97": Case(
        factory=lambda: _iapws("IF97"),
        p=1e5, T=500, T_sat=350, rtol=1e-4,
        marks=(skipif_no_iapws,),
    ),
    "iapws-if95": Case(
        factory=lambda: _iapws("IF95"),
        p=1e5, T=500, T_sat=350, rtol=1e-4,
        marks=(skipif_no_iapws,),
    ),
    "pyromat-mp-water": Case(
        factory=lambda: _pyromat("H2O", "mp"),
        p=1e5, T=500, T_sat=350, rtol=1e-4,
        marks=(skipif_no_pyromat,),
    ),
    "pyromat-ig-air": Case(
        factory=lambda: _pyromat("air", "ig"),
        p=1e5, T=500, T_sat=None, two_phase=False, rtol=1e-4,
        marks=(skipif_no_pyromat,),
    ),
    "thermopack-pr-co2": Case(
        factory=lambda: _thermopack("CO2", "PR"),
        p=3e6, T=320, T_sat=280,
        marks=(skipif_no_thermopack,),
    ),
    "thermopack-srk-propane": Case(
        factory=lambda: _thermopack("C3", "SRK"),
        p=5e5, T=320, T_sat=280,
        marks=(skipif_no_thermopack,),
    ),
    "thermopack-pr-mixture": Case(
        factory=lambda: _thermopack("C3[0.6]&IC4[0.4]|molar", "PR"),
        p=5e5, T=330, T_sat=280, mixture=True,
        marks=(skipif_no_thermopack,),
    ),
}


@pytest.fixture(
    params=[
        pytest.param(case, id=key, marks=case.marks)
        for key, case in WRAPPER_CASES.items()
    ],
    scope="module"
)
def case(request):
    c = request.param
    return c, c.factory()


class TestRoundtrips:

    def test_T_ph(self, case):
        c, w = case
        h = w.h_pT(c.p, c.T)
        assert w.T_ph(c.p, h) == pytest.approx(c.T, rel=c.rtol)

    def test_T_ps(self, case):
        c, w = case
        s = w.s_pT(c.p, c.T)
        assert w.T_ps(c.p, s) == pytest.approx(c.T, rel=c.rtol)

    def test_h_ps_s_ph(self, case):
        c, w = case
        h = w.h_pT(c.p, c.T)
        s = w.s_ph(c.p, h)
        assert w.h_ps(c.p, s) == pytest.approx(h, rel=c.rtol)

    def test_isentropic_keeps_entropy(self, case):
        c, w = case
        h_1 = w.h_pT(c.p, c.T)
        s_1 = w.s_ph(c.p, h_1)
        h_2 = w.isentropic(c.p, h_1, c.p * 1.5)
        assert w.s_ph(c.p * 1.5, h_2) == pytest.approx(s_1, rel=c.rtol)

    def test_d_ph_matches_d_pT(self, case):
        c, w = case
        h = w.h_pT(c.p, c.T)
        assert w.d_ph(c.p, h) == pytest.approx(w.d_pT(c.p, c.T), rel=c.rtol)


class TestSaturation:

    @pytest.fixture(autouse=True)
    def _skip_single_phase(self, case):
        if not case[0].two_phase:
            pytest.skip("Wrapper case does not support two-phase states.")

    def test_T_sat_p_sat(self, case):
        c, w = case
        if c.mixture:
            # p_sat and T_sat of the CoolPropWrapper flash at different
            # qualities, for zeotropic mixtures this is not a roundtrip
            pytest.skip("T_sat(p_sat(T)) only holds for pure fluids.")
        assert w.T_sat(w.p_sat(c.T_sat)) == pytest.approx(c.T_sat, rel=c.rtol)

    def test_T_bubble_p_bubble(self, case):
        c, w = case
        assert w.T_bubble(w.p_bubble(c.T_sat)) == pytest.approx(
            c.T_sat, rel=c.rtol
        )

    def test_T_dew_p_dew(self, case):
        c, w = case
        assert w.T_dew(w.p_dew(c.T_sat)) == pytest.approx(c.T_sat, rel=c.rtol)

    def test_glide(self, case):
        c, w = case
        p = w.p_sat(c.T_sat)
        if c.mixture:
            assert w.T_dew(p) > w.T_bubble(p)
        else:
            assert w.T_dew(p) == pytest.approx(w.T_bubble(p), rel=c.rtol)

    def test_h_pQ_monotone(self, case):
        c, w = case
        p = w.p_sat(c.T_sat)
        h = [w.h_pQ(p, Q) for Q in [0, 0.4, 1]]
        assert h[0] < h[1] < h[2]

    def test_Q_ph_h_pQ(self, case):
        c, w = case
        p = w.p_sat(c.T_sat)
        assert w.Q_ph(p, w.h_pQ(p, 0)) == pytest.approx(0, abs=1e-2)
        assert w.Q_ph(p, w.h_pQ(p, 1)) == pytest.approx(1, abs=1e-2)
        Q = w.Q_ph(p, w.h_pQ(p, 0.4))
        if c.mixture:
            # for mixtures h_pQ may interpolate between bubble and dew
            # state and Q_ph is mass based, only rough agreement expected
            assert 0.2 < Q < 0.6
        else:
            assert Q == pytest.approx(0.4, abs=1e-2)

    def test_h_QT_matches_h_pQ(self, case):
        c, w = case
        p = w.p_sat(c.T_sat)
        assert w.h_QT(0, c.T_sat) == pytest.approx(
            w.h_pQ(w.p_bubble(c.T_sat), 0), rel=1e-3
        )
        assert w.h_QT(1, c.T_sat) == pytest.approx(
            w.h_pQ(w.p_dew(c.T_sat), 1), rel=1e-3
        )


class TestReferenceCrossCheck:
    """Coarse comparison against CoolProp HEOS to catch unit conversion
    errors, e.g. molar instead of mass basis."""

    @pytest.fixture(autouse=True)
    def _reference(self, case):
        c, w = case
        if c.mixture or w.__class__ is CoolPropWrapper and w.back_end == "HEOS":
            pytest.skip("Case is its own reference or has no pure fluid.")
        names = {"H2O": "Water", "CO2": "CO2", "C3": "Propane", "air": "Air"}
        self.ref = CoolPropWrapper(names[w.fluid])

    def test_T_sat(self, case):
        c, w = case
        if not c.two_phase:
            pytest.skip("Wrapper case does not support two-phase states.")
        assert w.T_sat(c.p) == pytest.approx(self.ref.T_sat(c.p), abs=2)

    def test_h_pT_relative_to_reference_state(self, case):
        c, w = case
        dh = w.h_pT(c.p, c.T + 50) - w.h_pT(c.p, c.T)
        dh_ref = self.ref.h_pT(c.p, c.T + 50) - self.ref.h_pT(c.p, c.T)
        assert dh == pytest.approx(dh_ref, rel=5e-2)

    def test_d_pT(self, case):
        c, w = case
        assert w.d_pT(c.p, c.T) == pytest.approx(
            self.ref.d_pT(c.p, c.T), rel=5e-2
        )


class TestFluidStringParsing:

    def test_pure_fluid(self):
        w = FluidPropertyWrapper("CO2")
        w._identify_mixture()
        assert w.fluid == "CO2"
        assert w.fractions == []
        assert w.mixture_type is None

    @pytest.mark.parametrize("mixture_type", ["mass", "molar", "volume"])
    def test_mixture(self, mixture_type):
        w = FluidPropertyWrapper(f"C3[0.6]&IC4[0.4]|{mixture_type}")
        w._identify_mixture()
        assert w.fluid == "C3&IC4"
        assert w.fractions == [0.6, 0.4]
        assert w.mixture_type == mixture_type

    def test_missing_composition_type_raises(self):
        w = FluidPropertyWrapper("C3[0.6]&IC4[0.4]")
        with pytest.raises(ValueError):
            w._identify_mixture()

    def test_invalid_composition_type_raises(self):
        w = FluidPropertyWrapper("C3[0.6]&IC4[0.4]|volumetric")
        with pytest.raises(ValueError):
            w._identify_mixture()


class TestWrapperRegistry:

    def test_shipped_wrappers_registered(self):
        for name in [
            "FluidPropertyWrapper", "CoolPropWrapper",
            "IncompressibleFluidWrapper", "IAPWSWrapper", "PyromatWrapper",
            "ThermopackWrapper"
        ]:
            assert name in wrapper_registry.items

    def test_custom_wrapper_registers(self):
        @wrapper_registry
        class MyWrapper(FluidPropertyWrapper):
            pass

        assert wrapper_registry.items["MyWrapper"] is MyWrapper
        del wrapper_registry.items["MyWrapper"]
