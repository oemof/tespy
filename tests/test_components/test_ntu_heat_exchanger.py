# -*- coding: utf-8

"""Module for testing the effectiveness-NTU based heat exchanger components.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_ntu_heat_exchangers.py

SPDX-License-Identifier: MIT
"""
import math

import ht
import numpy as np
import pytest
from CoolProp.CoolProp import PropsSI as PSI
from pytest import approx
from pytest import mark

from tespy.components import HeatExchanger
from tespy.components import NTUHeatExchanger
from tespy.components import ParallelFlowHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.components.heat_exchangers.ntu import FLOW_ARRANGEMENTS
from tespy.components.heat_exchangers.ntu import calc_epsilon
from tespy.components.heat_exchangers.ntu import calc_ntu
from tespy.connections import Connection
from tespy.networks import Network
from tespy.tools.fluid_properties.wrappers import IncompressibleFluidWrapper

# Relations cross validated against the ht library, as
# (relation of calc_epsilon, ht subtype, number of shell passes, C_r values).
# C_r=1 is excluded for shell and tube, since ht does not implement the
# analytic limit of that relation.
_HT_PAIRS = [
    ("counterflow", "counterflow", 1, [0.3, 0.7, 1.0]),
    ("parallelflow", "parallel", 1, [0.3, 0.7, 1.0]),
    ("crossflow_both_unmixed", "crossflow approximate", 1, [0.3, 0.7, 1.0]),
    ("crossflow_min_mixed", "crossflow, mixed Cmin", 1, [0.3, 0.7, 1.0]),
    ("crossflow_max_mixed", "crossflow, mixed Cmax", 1, [0.3, 0.7, 1.0]),
    ("shell_and_tube", "S&T", 1, [0.3, 0.7]),
    ("shell_and_tube", "S&T", 2, [0.3, 0.7]),
    ("shell_and_tube", "S&T", 3, [0.3, 0.7]),
]

# ht does not implement a crossflow relation with both fluids mixed, so it
# cannot be cross validated. It is covered by the round trip and the
# reachability tests only.
_HT_UNSUPPORTED = {"crossflow_both_mixed"}

# flow arrangement of the component -> relations of calc_epsilon it resolves
# to; the one side mixed arrangements are assigned by capacity rate at
# runtime, so both relations have to be validated for either of them
_ARRANGEMENT_RELATIONS = {
    "counterflow": {"counterflow"},
    "parallelflow": {"parallelflow"},
    "crossflow_both_unmixed": {"crossflow_both_unmixed"},
    "crossflow_hot_mixed": {"crossflow_min_mixed", "crossflow_max_mixed"},
    "crossflow_cold_mixed": {"crossflow_min_mixed", "crossflow_max_mixed"},
    "crossflow_both_mixed": {"crossflow_both_mixed"},
    "shell_and_tube": {"shell_and_tube"},
}


class TestEpsNTURelations:

    @mark.parametrize("arrangement,num_shell_passes", [
        ("counterflow", 1),
        ("parallelflow", 1),
        ("crossflow_both_unmixed", 1),
        ("crossflow_min_mixed", 1),
        ("crossflow_max_mixed", 1),
        ("crossflow_both_mixed", 1),
        ("shell_and_tube", 1),
        ("shell_and_tube", 2),
        ("shell_and_tube", 3),
    ])
    def test_round_trip(self, arrangement, num_shell_passes):
        for ntu in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0]:
            for C_r in [0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0]:
                eps = calc_epsilon(ntu, C_r, arrangement, num_shell_passes)
                assert 0 < eps <= 1
                ntu_back = calc_ntu(eps, C_r, arrangement, num_shell_passes)
                eps_back = calc_epsilon(
                    ntu_back, C_r, arrangement, num_shell_passes
                )
                assert approx(eps_back, rel=1e-7) == eps

    def test_textbook_values(self):
        # any arrangement collapses to the same relation for C_r=0
        assert approx(calc_epsilon(1.0, 0.0, "counterflow")) == 1 - math.exp(-1)
        assert approx(calc_epsilon(1.0, 0.0, "shell_and_tube")) == 1 - math.exp(-1)
        # balanced counterflow: NTU / (1 + NTU)
        assert approx(calc_epsilon(2.0, 1.0, "counterflow")) == 2 / 3
        # balanced parallelflow: (1 - exp(-2 NTU)) / 2
        assert approx(
            calc_epsilon(2.0, 1.0, "parallelflow")
        ) == (1 - math.exp(-4)) / 2
        # balanced crossflow with both fluids unmixed (approximation)
        assert approx(
            calc_epsilon(1.0, 1.0, "crossflow_both_unmixed")
        ) == 1 - math.exp(math.exp(-1) - 1)
        # for C_r=1 both one side mixed relations coincide
        assert approx(
            calc_epsilon(1.0, 1.0, "crossflow_min_mixed")
        ) == calc_epsilon(1.0, 1.0, "crossflow_max_mixed")

    def test_effectiveness_above_arrangement_maximum(self):
        # parallelflow cannot exceed 1 / (1 + C_r)
        assert np.isnan(calc_ntu(0.9, 1.0, "parallelflow"))
        assert np.isnan(calc_ntu(0.999, 1.0, "crossflow_both_mixed"))
        assert np.isnan(calc_ntu(0.99, 1.0, "shell_and_tube"))
        assert np.isnan(calc_ntu(1.0, 0.0, "counterflow"))

    def test_invalid_arrangement(self):
        with pytest.raises(ValueError, match="flow arrangement"):
            calc_epsilon(1.0, 0.5, "crossflow")
        with pytest.raises(ValueError, match="flow arrangement"):
            calc_ntu(0.5, 0.5, "crossflow")

    def test_cross_validation_against_ht_library(self):
        # C_r=0 for any arrangement corresponds to the ht 'boiler' subtype
        for arrangement, subtype, num_shell_passes, C_r_range in _HT_PAIRS:
            kwargs = (
                {"n_shell_tube": num_shell_passes} if subtype == "S&T" else {}
            )
            for ntu in [0.1, 0.5, 1.0, 2.0, 5.0]:
                for C_r in C_r_range:
                    reference = ht.effectiveness_from_NTU(
                        ntu, C_r, subtype=subtype, **kwargs
                    )
                    assert approx(reference, abs=1e-12) == calc_epsilon(
                        ntu, C_r, arrangement, num_shell_passes
                    )
                    assert approx(ntu, abs=1e-9) == calc_ntu(
                        reference, C_r, arrangement, num_shell_passes
                    )
                reference = ht.effectiveness_from_NTU(ntu, 0.0, subtype="boiler")
                assert approx(reference, abs=1e-12) == calc_epsilon(
                    ntu, 0.0, arrangement, num_shell_passes
                )

    def test_every_flow_arrangement_is_cross_validated(self):
        """Guard against adding an arrangement without cross validation."""
        assert set(_ARRANGEMENT_RELATIONS) == set(FLOW_ARRANGEMENTS)

        validated = {relation for relation, *_ in _HT_PAIRS}
        for arrangement, relations in _ARRANGEMENT_RELATIONS.items():
            if arrangement in _HT_UNSUPPORTED:
                continue
            missing = relations - validated
            assert not missing, (
                f"the flow arrangement '{arrangement}' resolves to the "
                f"relations {sorted(missing)}, which are not cross validated "
                "against ht - add them to _HT_PAIRS, or to _HT_UNSUPPORTED if "
                "ht does not implement them"
            )

    def test_crossflow_approximation_error(self):
        """Bound the error of the crossflow both unmixed approximation.

        The implemented relation is the explicit approximation of the exact
        double series solution, which ht provides as its :code:`'crossflow'`
        subtype. The largest deviation over the range below is 3.33 % at
        NTU=0.5 and C_r=1.
        """
        worst = 0.0
        for ntu in [0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0]:
            for C_r in [0.1, 0.25, 0.5, 0.75, 1.0]:
                approximation = calc_epsilon(ntu, C_r, "crossflow_both_unmixed")
                exact = ht.effectiveness_from_NTU(ntu, C_r, subtype="crossflow")
                worst = max(worst, abs(approximation - exact) / exact)

        assert worst < 0.035


class TestNTUHeatExchanger:

    def setup_method(self):
        self.nw = Network(iterinfo=False)
        self.nw.units.set_defaults(**{
            "pressure": "bar", "pressure_difference": "bar",
            "temperature": "degC"
        })

    def setup_network(self, instance):
        self.c1 = Connection(Source('inlet 1'), 'out1', instance, 'in1')
        self.c2 = Connection(instance, 'out1', Sink('outlet 1'), 'in1')
        self.c3 = Connection(Source('inlet 2'), 'out1', instance, 'in2')
        self.c4 = Connection(instance, 'out2', Sink('outlet 2'), 'in1')
        self.nw.add_conns(self.c1, self.c2, self.c3, self.c4)

    def setup_reference_network(self, reference):
        nw_ref = Network(iterinfo=False)
        nw_ref.units.set_defaults(**{
            "pressure": "bar", "pressure_difference": "bar",
            "temperature": "degC"
        })
        c1 = Connection(Source("source 1"), "out1", reference, "in1")
        c2 = Connection(reference, "out1", Sink("sink 1"), "in1")
        c3 = Connection(Source("source 2"), "out1", reference, "in2")
        c4 = Connection(reference, "out2", Sink("sink 2"), "in1")
        nw_ref.add_conns(c1, c2, c3, c4)
        return nw_ref, c1, c2, c3, c4

    def test_counterflow_matches_lmtd(self):
        instance = NTUHeatExchanger("heat exchanger")
        self.setup_network(instance)

        self.c1.set_attr(fluid={"air": 1}, m=1, T=85, p=1)
        self.c3.set_attr(fluid={"water": 1}, m=0.5, T=25, p=3)
        instance.set_attr(
            flow_arrangement="counterflow", pr1=1, pr2=1, UA=1500
        )

        self.nw.solve("design")
        self.nw.assert_convergence()

        reference = HeatExchanger("reference heat exchanger")
        nw_ref, c1, c2, c3, c4 = self.setup_reference_network(reference)
        c1.set_attr(fluid={"air": 1}, m=1, T=85, p=1)
        c3.set_attr(fluid={"water": 1}, m=0.5, T=25, p=3)
        reference.set_attr(pr1=1, pr2=1, UA=1500)
        nw_ref.solve("design")
        nw_ref.assert_convergence()

        # the models are equivalent for constant capacity rates, real fluid
        # cp variation causes a small deviation
        assert approx(self.c2.T.val_SI, abs=0.15) == c2.T.val_SI
        assert approx(self.c4.T.val_SI, abs=0.15) == c4.T.val_SI

    def test_parallelflow_matches_lmtd(self):
        instance = NTUHeatExchanger("heat exchanger")
        self.setup_network(instance)

        self.c1.set_attr(fluid={"air": 1}, m=1, T=85, p=1)
        self.c3.set_attr(fluid={"water": 1}, m=0.5, T=25, p=3)
        instance.set_attr(flow_arrangement="parallelflow", pr1=1, pr2=1, UA=1500)

        self.nw.solve("design")
        self.nw.assert_convergence()

        reference = ParallelFlowHeatExchanger("reference heat exchanger")
        nw_ref, c1, c2, c3, c4 = self.setup_reference_network(reference)
        c1.set_attr(fluid={"air": 1}, m=1, T=85, p=1)
        c3.set_attr(fluid={"water": 1}, m=0.5, T=25, p=3)
        reference.set_attr(pr1=1, pr2=1, UA=1500)
        nw_ref.solve("design")
        nw_ref.assert_convergence()

        assert approx(self.c2.T.val_SI, abs=0.15) == c2.T.val_SI
        assert approx(self.c4.T.val_SI, abs=0.15) == c4.T.val_SI

    def test_crossflow_offdesign(self):
        instance = NTUHeatExchanger("heat exchanger")
        self.setup_network(instance)

        self.c1.set_attr(fluid={"air": 1}, m=1, T=85, p=1)
        self.c2.set_attr(T=45, design=["T"])
        self.c3.set_attr(fluid={"water": 1}, m=3, T=25, p=3)
        instance.set_attr(
            flow_arrangement="crossflow_both_unmixed", pr1=1, pr2=1,
            offdesign=["UA"]
        )

        self.nw.solve("design")
        self.nw.assert_convergence()
        design_state = self.nw.save(as_dict=True)

        # air is the minimum capacity rate side, so the design effectiveness
        # is the air temperature ratio (85 - 45) / (85 - 25)
        eps_design = instance.eps.val
        assert approx(eps_design) == 2 / 3
        assert approx(eps_design) == calc_epsilon(
            instance.NTU.val, instance.C_r.val, "crossflow_both_unmixed"
        )

        self.c1.set_attr(m=0.75)
        self.nw.solve("offdesign", design_path=design_state)
        self.nw.assert_convergence()

        assert approx(instance.UA.val_SI) == instance.UA.design
        assert approx(instance.eps.val) == calc_epsilon(
            instance.NTU.val, instance.C_r.val, "crossflow_both_unmixed"
        )
        # the reduced air flow lowers the minimum capacity rate at fixed UA,
        # so the number of transfer units and the effectiveness increase
        assert instance.NTU.val > instance.NTU.design
        assert instance.eps.val > eps_design

    def test_condensing_hot_side(self):
        instance = NTUHeatExchanger("heat exchanger")
        self.setup_network(instance)

        self.c1.set_attr(fluid={"water": 1}, T=50, x=1)
        self.c2.set_attr(x=0)
        self.c3.set_attr(fluid={"water": 1}, m=30, T=20, p=2)
        instance.set_attr(
            flow_arrangement="shell_and_tube", dp1=0, dp2=0, Q=-2e6
        )

        self.nw.solve("design")
        self.nw.assert_convergence()

        # isothermal hot side: C_r = 0, all arrangements collapse to the
        # same relation
        C_cold = self.c3.m.val_SI * (
            (self.c4.h.val_SI - self.c3.h.val_SI)
            / (self.c4.T.val_SI - self.c3.T.val_SI)
        )
        eps = 2e6 / (C_cold * (self.c1.T.val_SI - self.c3.T.val_SI))
        NTU = -math.log(1 - eps)

        assert approx(instance.C_r.val) == 0
        assert approx(instance.eps.val) == eps
        assert approx(instance.NTU.val) == NTU
        assert approx(instance.UA.val_SI) == NTU * C_cold

    def test_coolsolve_exchangers1(self):
        # counterflow oil cooler with given UA from the CoolSolve examples.
        instance = NTUHeatExchanger("heat exchanger")
        self.setup_network(instance)

        UA = 12.5 * 400
        c_oil = 2100
        c_w = PSI("C", "T", 293.15, "P", 1e5, "water")
        C_oil = 2 * c_oil
        C_w = 0.48 * c_w
        omega = C_w / C_oil
        NTU = UA / C_w
        e = math.exp(-NTU * (1 - omega))
        epsilon = (1 - e) / (1 - omega * e)
        Q_ref = epsilon * C_w * 80

        oil_data = {
            "temperature_data": np.array([273.15, 473.15]),
            "heat_capacity_data": np.array([c_oil, c_oil], dtype=float),
            "density_data": np.array([850.0, 850.0]),
            "viscosity_data": np.array([5e-3, 5e-3])
        }
        self.c1.set_attr(
            fluid={"oil": 1},
            fluid_engines={"oil": IncompressibleFluidWrapper},
            fluid_wrapper_kwargs={"oil": oil_data},
            m=2, T=100, p=1
        )
        self.c3.set_attr(fluid={"water": 1}, m=0.48, T=20, p=1)
        instance.set_attr(flow_arrangement="counterflow", dp1=0, dp2=0, UA=UA)

        self.nw.solve("design")
        self.nw.assert_convergence()

        # deviation from the water heat capacity: the reference uses the inlet
        # value, the component the effective mean value
        assert approx(-instance.Q.val, rel=2e-3) == Q_ref
        assert approx(self.c2.T.val, abs=0.1) == 100 - Q_ref / C_oil
        assert approx(self.c4.T.val, abs=0.1) == 20 + Q_ref / C_w

    def test_coolsolve_exchangers2(self):
        # crossflow water heater solving for UA from the CoolSolve examples.
        instance = NTUHeatExchanger("heat exchanger")
        self.setup_network(instance)

        c_gaz = 1000
        c_w = PSI("C", "T", 273.15 + 80, "P", 5e5, "water")
        Q_ref = c_w * 90
        M_gaz = Q_ref / (c_gaz * 200)
        C_gaz = M_gaz * c_gaz
        omega = C_gaz / c_w
        epsilon = Q_ref / (C_gaz * 265)
        NTU_ref = calc_ntu(epsilon, omega, "crossflow_both_unmixed")
        UA_ref = NTU_ref * C_gaz

        gas_data = {
            "temperature_data": np.array([273.15, 673.15]),
            "heat_capacity_data": np.array([c_gaz, c_gaz], dtype=float),
            "density_data": np.array([0.9, 0.9]),
            "viscosity_data": np.array([3e-5, 3e-5])
        }
        self.c1.set_attr(
            fluid={"gas": 1},
            fluid_engines={"gas": IncompressibleFluidWrapper},
            fluid_wrapper_kwargs={"gas": gas_data},
            T=300, p=1
        )
        self.c2.set_attr(T=100)
        self.c3.set_attr(fluid={"water": 1}, m=1, T=35, p=5)
        self.c4.set_attr(T=125)
        instance.set_attr(
            flow_arrangement="crossflow_both_unmixed", dp1=0, dp2=0
        )

        self.nw.solve("design")
        self.nw.assert_convergence()

        assert approx(-instance.Q.val, rel=2e-3) == Q_ref
        assert approx(self.c1.m.val, rel=2e-3) == M_gaz
        assert approx(instance.UA.val, rel=5e-3) == UA_ref

    def test_coolsolve_exchangers3(self):
        # shell and tube steam condenser from the CoolSolve examples.
        instance = NTUHeatExchanger("heat exchanger")
        self.setup_network(instance)

        Q_ref = 2000e6
        T_w_ex = 20
        # the published model evaluates cp at the mean water temperature,
        # which itself depends on the outlet temperature
        for _ in range(5):
            cp_w = PSI("C", "T", 273.15 + (20 + T_w_ex) / 2, "P", 1e5, "water")
            T_w_ex = 20 + Q_ref / (30000 * cp_w)
        C_w = 30000 * cp_w
        epsilon = Q_ref / (C_w * 30)
        NTU_ref = -math.log(1 - epsilon)

        self.c1.set_attr(fluid={"water": 1}, T=50, x=1)
        self.c2.set_attr(x=0)
        self.c3.set_attr(fluid={"water": 1}, m=30000, T=20, p=1)
        instance.set_attr(
            flow_arrangement="shell_and_tube", dp1=0, dp2=0, Q=-Q_ref
        )

        self.nw.solve("design")
        self.nw.assert_convergence()

        assert approx(self.c4.T.val, abs=0.1) == T_w_ex
        assert approx(instance.eps.val, rel=5e-3) == epsilon
        assert approx(instance.UA.val, rel=1e-2) == NTU_ref * C_w
        assert instance.C_r.val == 0

    def test_labothappy_reference(self):
        # water to water counterflow plate heat exchanger cross-validated
        # against the LaboThapPy HexeNTU component, the UA value stems from its
        # plate geometry model (Gnielinski correlation, plate conduction and
        # fouling)
        instance = NTUHeatExchanger("heat exchanger")
        self.setup_network(instance)

        self.c1.set_attr(fluid={"water": 1}, m=0.014, T=90, p=4)
        self.c3.set_attr(fluid={"water": 1}, m=0.08, T=12, p=4)
        instance.set_attr(
            flow_arrangement="counterflow", dp1=0, dp2=0,
            UA=132.069586
        )

        self.nw.solve("design")
        self.nw.assert_convergence()

        # LaboThapPy evaluates the heat capacities at the inlet states, the
        # component uses effective mean values
        assert approx(-instance.Q.val, rel=2e-3) == 3961.899
        assert approx(self.c2.T.val_SI, abs=0.1) == 295.54
        assert approx(self.c4.T.val_SI, abs=0.1) == 296.98

    def test_unspecified_arrangement(self):
        """The flow arrangement has no default and has to be specified."""
        instance = NTUHeatExchanger("heat exchanger")
        self.setup_network(instance)

        self.c1.set_attr(fluid={"air": 1}, m=1, T=85, p=1)
        self.c3.set_attr(fluid={"water": 1}, m=3, T=25, p=3)
        instance.set_attr(pr1=1, pr2=1, UA=1500)

        assert not instance.flow_arrangement.is_set
        with pytest.raises(ValueError, match="not specified"):
            self.nw.solve("design", init_only=True)

        # specifying it resolves the error
        instance.set_attr(flow_arrangement="counterflow")
        self.nw.solve("design")
        self.nw.assert_convergence()

    def test_invalid_arrangement(self):
        instance = NTUHeatExchanger("heat exchanger")
        self.setup_network(instance)

        self.c1.set_attr(fluid={"air": 1}, m=1, T=85, p=1)
        self.c2.set_attr(T=45)
        self.c3.set_attr(fluid={"water": 1}, m=3, T=25, p=3)
        instance.set_attr(flow_arrangement="crossflow", pr1=1, pr2=1)

        with pytest.raises(ValueError, match="flow arrangement"):
            self.nw.solve("design", init_only=True)
