# -*- coding: utf-8 -*-

"""Phase-vector checks for the sc->l hot side and l->tp->g cold side.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_sectioned_hx_sc_phases.py

SPDX-License-Identifier: MIT
"""
import numpy as np
from types import SimpleNamespace

from tespy.components.heat_exchangers.sectioned import SectionedHeatExchanger
from tespy.tools.fluid_properties.wrappers import CoolPropWrapper


def _conn(p, h, wrapper, name="Water"):
    fd = {name: {"wrapper": wrapper, "mass_fraction": 1.0}}
    return SimpleNamespace(
        p=SimpleNamespace(val_SI=p),
        h=SimpleNamespace(val_SI=h),
        fluid=SimpleNamespace(val={name: 1.0}),
        fluid_data=fd,
    )


class TestSectionedHXScToLiquidPhases:
    """Hot side crosses T_crit at P > P_crit (sc->l).
    Cold side crosses the full saturation dome at P < P_crit (l->tp->g).

    Water: P_crit = 22.064 MPa, T_crit = 647.1 K.

    Note: sc/l boundary detection is disabled pending a CoolProp fix for
    h_pT(p, T_crit) at p > p_crit (fixed in CoolProp 8). The hot side is
    therefore treated as a single supercritical zone.
    """

    def setup_method(self):
        w = CoolPropWrapper("Water")

        # hot side: 30 MPa (> P_crit), 700 K -> 350 K crosses T_crit
        p_h = 30e6
        self.c_h_in  = _conn(p_h, w.h_pT(p_h, 700.0), w)   # sc
        self.c_h_out = _conn(p_h, w.h_pT(p_h, 350.0), w)   # l

        # cold side: 1 MPa (< P_crit), 300 K -> 600 K crosses saturation dome
        p_c = 1e6
        self.c_c_in  = _conn(p_c, w.h_pT(p_c, 300.0), w)   # l
        self.c_c_out = _conn(p_c, w.h_pT(p_c, 600.0), w)   # g

    def test_hot_zone_phases(self):
        _, zone_phases = SectionedHeatExchanger._get_moving_steps(
            self.c_h_in, self.c_h_out
        )
        # single supercritical zone — no T_crit boundary until CoolProp 8
        assert zone_phases == [3], f"got {zone_phases}"

    def test_cold_zone_phases(self):
        _, zone_phases = SectionedHeatExchanger._get_moving_steps(
            self.c_c_in, self.c_c_out
        )
        # liquid (0) -> two-phase (1) -> gas (2)
        assert zone_phases == [0, 1, 2], f"got {zone_phases}"

    def test_hot_no_interior_step(self):
        steps, _ = SectionedHeatExchanger._get_moving_steps(
            self.c_h_in, self.c_h_out
        )
        assert steps == [0, 1], f"expected no interior boundary, got {steps}"

    def test_cold_step_boundaries_are_ordered(self):
        steps, _ = SectionedHeatExchanger._get_moving_steps(
            self.c_c_in, self.c_c_out
        )
        assert len(steps) == 4
        assert steps[0] < steps[1] < steps[2] < steps[3]

    def test_section_phases_hot(self):
        steps_h, zp_h = SectionedHeatExchanger._get_moving_steps(
            self.c_h_in, self.c_h_out
        )
        steps_c, zp_c = SectionedHeatExchanger._get_moving_steps(
            self.c_c_in, self.c_c_out
        )
        steps_all = np.unique(np.r_[np.linspace(0, 1, 51), steps_h, steps_c])
        phases_h = SectionedHeatExchanger._section_phases(
            steps_all, np.array(steps_h), zp_h
        )
        assert set(phases_h) == {3}, (
            f"hot side should be all supercritical (3), got {sorted(set(phases_h))}"
        )

    def test_section_phases_cold(self):
        steps_h, zp_h = SectionedHeatExchanger._get_moving_steps(
            self.c_h_in, self.c_h_out
        )
        steps_c, zp_c = SectionedHeatExchanger._get_moving_steps(
            self.c_c_in, self.c_c_out
        )
        steps_all = np.unique(np.r_[np.linspace(0, 1, 51), steps_h, steps_c])
        phases_c = SectionedHeatExchanger._section_phases(
            steps_all, np.array(steps_c), zp_c
        )
        assert set(phases_c) == {0, 1, 2}, (
            f"cold side should have liquid (0), two-phase (1), gas (2), got {sorted(set(phases_c))}"
        )

    def test_combined_phase_vector_sequence(self):
        """The merged grid produces three distinct (ph_hot, ph_cold) regions:
        (sc,l) -> (sc,tp) -> (sc,g). Hot side is a single sc zone throughout.
        """
        steps_h, zp_h = SectionedHeatExchanger._get_moving_steps(
            self.c_h_in, self.c_h_out
        )
        steps_c, zp_c = SectionedHeatExchanger._get_moving_steps(
            self.c_c_in, self.c_c_out
        )
        steps_all = np.unique(np.r_[np.linspace(0, 1, 51), steps_h, steps_c])
        phases_h = SectionedHeatExchanger._section_phases(
            steps_all, np.array(steps_h), zp_h
        )
        phases_c = SectionedHeatExchanger._section_phases(
            steps_all, np.array(steps_c), zp_c
        )

        pairs = list(zip(phases_h, phases_c))
        distinct = [pairs[0]]
        for pair in pairs[1:]:
            if pair != distinct[-1]:
                distinct.append(pair)

        assert distinct == [(3, 0), (3, 1), (3, 2)], (
            f"unexpected combined phase sequence: {distinct}"
        )


class TestPhaseVectorOrderingScToLiquid:
    """Complete phase-vector check for sc->l hot side, l->tp->g cold side.

    Moving boundary: sections defined by phase boundaries only (3 sections).
    Sectioned num_sections=4: phase boundaries inserted into uniform grid (6 sections).

    0=liquid, 1=two-phase, 2=gas, 3=supercritical.
    Hot side is a single sc zone (no T_crit boundary until CoolProp 8).
    """

    def setup_method(self):
        w = CoolPropWrapper("Water")
        p_h, p_c = 30e6, 1e6
        self.c_h_in  = _conn(p_h, w.h_pT(p_h, 700.0), w)
        self.c_h_out = _conn(p_h, w.h_pT(p_h, 350.0), w)
        self.c_c_in  = _conn(p_c, w.h_pT(p_c, 300.0), w)
        self.c_c_out = _conn(p_c, w.h_pT(p_c, 600.0), w)

    def test_moving_boundary_phase_vectors(self):
        steps_h, zp_h = SectionedHeatExchanger._get_moving_steps(self.c_h_in, self.c_h_out)
        steps_c, zp_c = SectionedHeatExchanger._get_moving_steps(self.c_c_in, self.c_c_out)
        # MovingBoundary: steps_all = [0, x_liq, x_gas, 1] → 3 sections
        steps_all = np.unique(np.r_[steps_h, steps_c])
        phases_h = SectionedHeatExchanger._section_phases(steps_all, np.array(steps_h), zp_h)
        phases_c = SectionedHeatExchanger._section_phases(steps_all, np.array(steps_c), zp_c)
        assert len(phases_h) == 3
        assert phases_h == [3, 3, 3]
        assert phases_c == [0, 1, 2]

    def test_sectioned_4_sections_phase_vectors(self):
        steps_h, zp_h = SectionedHeatExchanger._get_moving_steps(self.c_h_in, self.c_h_out)
        steps_c, zp_c = SectionedHeatExchanger._get_moving_steps(self.c_c_in, self.c_c_out)
        # Sectioned num_sections=4: uniform [0,.25,.5,.75,1] + x_liq, x_gas → 6 sections
        steps_all = np.unique(np.r_[np.linspace(0, 1, 5), steps_h, steps_c])
        phases_h = SectionedHeatExchanger._section_phases(steps_all, np.array(steps_h), zp_h)
        phases_c = SectionedHeatExchanger._section_phases(steps_all, np.array(steps_c), zp_c)
        assert len(phases_h) == 6
        assert phases_h == [3, 3, 3, 3, 3, 3]
        assert phases_c == [0, 1, 1, 1, 1, 2]
