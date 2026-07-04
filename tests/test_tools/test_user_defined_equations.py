# -*- coding: utf-8

"""Module for testing UserDefinedEquation.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tests/test_user_defined_equations.py

SPDX-License-Identifier: MIT
"""

import pytest
from pytest import approx

from tespy.components import Compressor
from tespy.components import CycleCloser
from tespy.components import SimpleHeatExchanger
from tespy.components import Valve
from tespy.connections import Connection
from tespy.networks import Network
from tespy.tools import UserDefinedEquation
from tespy.tools.fluid_properties.functions import T_dew_p


def _build_refrigeration_cycle():
    """Return a solved refrigeration cycle and its key connections.

    Topology: evaporator -> superheater -> compressor -> condenser -> valve -> cc -> evaporator

    The evaporator and superheater are both SimpleHeatExchanger components so
    that the superheat can be controlled independently of the saturation state.
    """
    nw = Network()
    nw.iterinfo = False
    nw.units.set_defaults(
        temperature="°C", pressure="bar", pressure_difference="bar",
        power="kW", heat="kW"
    )

    evaporator = SimpleHeatExchanger("evaporator")
    superheater = SimpleHeatExchanger("superheater")
    compressor = Compressor("compressor")
    condenser = SimpleHeatExchanger("condenser")
    valve = Valve("valve")
    cc = CycleCloser("cc")

    c1 = Connection(evaporator, "out1", superheater, "in1", label="c1")
    c1b = Connection(superheater, "out1", compressor, "in1", label="c1b")
    c2 = Connection(compressor, "out1", condenser, "in1", label="c2")
    c3 = Connection(condenser, "out1", valve, "in1", label="c3")
    c4 = Connection(valve, "out1", cc, "in1", label="c4")
    c5 = Connection(cc, "out1", evaporator, "in1", label="c5")

    nw.add_conns(c1, c1b, c2, c3, c4, c5)

    c1.set_attr(fluid={"R600": 1}, T_dew=10, td_dew=2)
    c1b.set_attr(td_dew=12)
    c3.set_attr(T_bubble=60, td_bubble=5)

    evaporator.set_attr(dp=0)
    superheater.set_attr(dp=1)
    condenser.set_attr(dp=0, Q=-100)
    compressor.set_attr(eta_s=0.8)

    nw.solve("design")

    return nw, c1, c1b


def _make_superheat_ude(conns, superheat):
    """Build a UDE that fixes the superheating referenced to the dew temperature
    at the evaporator outlet pressure (c1), not the superheater outlet pressure."""
    def func(ude):
        c1, c1b = ude.conns
        return (
            c1b.calc_T() - T_dew_p(c1.p.val_SI, c1.fluid_data)
            - ude.params["superheat"]
        )

    def dependents(ude):
        c1, c1b = ude.conns
        return [c1.p, c1b.p, c1b.h]

    return UserDefinedEquation(
        "superheat_ude",
        func,
        dependents,
        conns=list(conns),
        params={"superheat": superheat},
    )


class TestUserDefinedEquationIsSet:

    def setup_method(self):
        nw, c1, c1b = _build_refrigeration_cycle()
        self.nw = nw
        self.c1 = c1
        self.c1b = c1b

        # Replace the simple td_dew spec on c1b with our UDE
        c1b.set_attr(td_dew=None)
        self.ude = _make_superheat_ude([c1, c1b], superheat=10)
        nw.add_ude(self.ude)
        nw.solve("design")
        # Reference result: superheat should be ~10 K above dew temperature at c1 pressure
        self.T_c1b_with_ude = c1b.T.val

    def test_default_is_set_true(self):
        ude = _make_superheat_ude([self.c1, self.c1b], superheat=10)
        assert ude.is_set is True

    def test_is_set_rejects_non_bool(self):
        with pytest.raises(TypeError):
            self.ude.is_set = "yes"

    def test_deactivate_removes_equation_from_preprocess(self):
        self.ude.is_set = False
        self.ude._preprocess(0)
        assert self.ude.num_eq == 0
        assert self.ude._equation_set_lookup == {}

    def test_activate_includes_equation_in_preprocess(self):
        self.ude.is_set = True
        self.ude._preprocess(0)
        assert self.ude.num_eq == 1
        assert 0 in self.ude._equation_set_lookup

    def test_deactivate_changes_solve_result(self):
        """Deactivating the UDE frees a degree of freedom; adding a direct
        td_dew spec should then give a different temperature at c1b."""
        self.ude.is_set = False
        # With UDE inactive we need another constraint - use a different superheat
        self.c1b.set_attr(td_dew=6)
        self.nw.solve("design")

        T_without_ude = self.c1b.T.val
        assert T_without_ude != approx(self.T_c1b_with_ude, rel=1e-3)

    def test_reactivate_restores_solve_result(self):
        """Reactivating the UDE (after removing the competing constraint)
        should reproduce the original result."""
        # Deactivate and solve with a different spec
        self.ude.is_set = False
        self.c1b.set_attr(td_dew=6)
        self.nw.solve("design")

        # Restore UDE, drop the competing constraint
        self.c1b.set_attr(td_dew=None)
        self.ude.is_set = True
        self.nw.solve("design")

        assert self.c1b.T.val == approx(self.T_c1b_with_ude, rel=1e-3)

    def test_ude_stays_in_network_after_deactivation(self):
        self.ude.is_set = False
        assert "superheat_ude" in self.nw.user_defined_eq

    def test_toggle_without_network_resolve(self):
        """is_set can be toggled on the object without calling solve."""
        self.ude.is_set = False
        assert self.ude.is_set is False
        self.ude.is_set = True
        assert self.ude.is_set is True
