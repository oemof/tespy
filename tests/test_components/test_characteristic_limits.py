# -*- coding: utf-8

"""Tests of the range check of characteristics after convergence.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_characteristic_limits.py

SPDX-License-Identifier: MIT
"""
import logging

import pytest

from tespy.components import Motor
from tespy.components import PowerSink
from tespy.components import PowerSource
from tespy.components import Sink
from tespy.components import Source
from tespy.components import TurboCompressor
from tespy.components import Valve
from tespy.connections import Connection
from tespy.connections import PowerConnection
from tespy.networks import Network
from tespy.tools.characteristics import CharLine
from tespy.tools.characteristics import CharMap
from tespy.tools.characteristics import load_default_char as ldc


def _domain_messages(caplog):
    return [m for m in caplog.messages if "Operating point outside" in m]


class TestMotorEfficiencyLine:

    def setup_method(self):
        self.nw = Network(iterinfo=False)
        self.motor = Motor("motor")
        self.e1 = PowerConnection(PowerSource("grid"), "power", self.motor, "power_in")
        self.e2 = PowerConnection(self.motor, "power_out", PowerSink("load"), "power")
        self.nw.add_conns(self.e1, self.e2)
        self.motor.set_attr(eta=0.95, design=["eta"], offdesign=["eta_char"])
        self.e2.set_attr(E=100e3)

    def _design(self, extrapolate):
        self.motor.set_attr(eta_char=CharLine(
            x=[0.5, 1.0, 1.2], y=[0.9, 1.0, 1.02], extrapolate=extrapolate
        ))
        self.nw.solve("design")
        assert self.nw.status == 0
        assert self.motor.characteristics_in_use() == []
        return self.nw.save()

    def test_inside_range(self, caplog):
        design = self._design(extrapolate=False)
        self.e2.set_attr(E=80e3)
        with caplog.at_level(logging.WARNING, logger="TESPyLogger"):
            self.nw.solve("offdesign", design_path=design)
        assert self.nw.status == 0
        assert self.motor.characteristics_in_use() == ["eta_char"]
        assert self.motor.get_char_inputs("eta_char") == pytest.approx((0.8,))
        assert _domain_messages(caplog) == []

    def test_outside_range_evaluates_at_bound(self, caplog):
        design = self._design(extrapolate=False)
        self.e2.set_attr(E=30e3)
        with caplog.at_level(logging.WARNING, logger="TESPyLogger"):
            self.nw.solve("offdesign", design_path=design)
        assert self.nw.status == 1
        assert self.motor.eta.val == pytest.approx(0.95 * 0.9)
        messages = _domain_messages(caplog)
        assert len(messages) == 1
        assert "eta_char: x = 0.3 below minimum value (0.5)" in messages[0]
        assert "evaluated at that bound" in messages[0]

    def test_outside_range_extrapolates(self, caplog):
        design = self._design(extrapolate=True)
        self.e2.set_attr(E=30e3)
        with caplog.at_level(logging.WARNING, logger="TESPyLogger"):
            self.nw.solve("offdesign", design_path=design)
        assert self.nw.status == 0
        assert self.motor.eta.val == pytest.approx(0.95 * 0.86)
        messages = _domain_messages(caplog)
        assert len(messages) == 1
        assert "extrapolated" in messages[0]


class TestTurboCompressorMap:

    def setup_method(self):
        self.nw = Network(iterinfo=False)
        self.nw.units.set_defaults(pressure="bar", temperature="degC")
        self.comp = TurboCompressor("compressor")
        self.c1 = Connection(Source("source"), "out1", self.comp, "in1")
        self.c2 = Connection(self.comp, "out1", Sink("sink"), "in1")
        self.nw.add_conns(self.c1, self.c2)
        self.c1.set_attr(fluid={"N2": 1}, m=1, p=1, T=5)
        self.c2.set_attr(p=6)
        self.comp.set_attr(eta_s=0.8)
        self.nw.solve("design")
        assert self.nw.status == 0
        self.design = self.nw.save()
        self.c2.set_attr(p=None)
        self.comp.set_attr(
            char_map_pr={"char_func": ldc("TurboCompressor", "char_map_pr", "DEFAULT", CharMap), "is_set": True},
            char_map_eta_s={"char_func": ldc("TurboCompressor", "char_map_eta_s", "DEFAULT", CharMap), "is_set": True},
            eta_s=None, igva=0
        )

    def test_design_point_inside_map(self, caplog):
        with caplog.at_level(logging.WARNING, logger="TESPyLogger"):
            self.nw.solve("offdesign", design_path=self.design)
        assert self.nw.status == 0
        assert self.comp.characteristics_in_use() == ["char_map_eta_s", "char_map_pr"]
        assert self.comp.get_char_inputs("char_map_pr") == pytest.approx((1.0, 1.0))
        assert _domain_messages(caplog) == []

    def test_outside_speedlines(self, caplog):
        self.c1.set_attr(T=-30)
        with caplog.at_level(logging.WARNING, logger="TESPyLogger"):
            self.nw.solve("offdesign", design_path=self.design)
        assert self.nw.converged
        assert self.nw.status == 1
        messages = _domain_messages(caplog)
        assert any("char_map_pr: x =" in m and "above maximum value" in m for m in messages)
        assert any("char_map_eta_s: x =" in m for m in messages)

    def test_igva_scales_mass_flow_axis(self):
        self.nw.solve("offdesign", design_path=self.design)
        self.comp.set_attr(igva=20)
        self.nw.solve("offdesign", design_path=self.design)
        x, y_unscaled = self.comp.get_char_inputs("char_map_pr")
        _, y = self.comp._char_map_inputs()
        assert y_unscaled == pytest.approx(y / 0.8)


def test_valve_Kv_char_inputs_are_the_opening():
    nw = Network(iterinfo=False)
    nw.units.set_defaults(pressure="bar", temperature="degC")
    valve = Valve("valve")
    c1 = Connection(Source("source"), "out1", valve, "in1")
    c2 = Connection(valve, "out1", Sink("sink"), "in1")
    nw.add_conns(c1, c2)
    c1.set_attr(fluid={"water": 1}, p=5, T=20, m=10)
    c2.set_attr(p=4)
    valve.set_attr(
        Kv_char={"char_func": CharLine(x=[0, 0.5, 1], y=[0, 100, 200]), "is_set": True},
        opening="var"
    )
    nw.solve("design")
    assert nw.status == 0
    assert valve.characteristics_in_use() == ["Kv_char"]
    assert valve.get_char_inputs("Kv_char") == (valve.opening.val_SI,)
