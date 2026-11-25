# -*- coding: utf-8

"""Module for testing the value bounds the solver applies to the variables.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_value_bounds.py

SPDX-License-Identifier: MIT
"""
import logging

import pytest
from CoolProp.CoolProp import HAPropsSI

from tespy.components import Compressor
from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Turbine
from tespy.connections import Connection
from tespy.connections import HAConnection
from tespy.networks import Network
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.fluid_properties.mixtures import w_mix_pT_humidair

HYDROGEN_RICH = {"H2": 0.9, "N2": 0.1}


def _turbomachine_network(
        fluid, T_in, p_in, p_out, machine=Turbine, h_range=None, **kwargs
    ):
    """The outlet enthalpy is solved by the component's efficiency equation,
    the outlet connection does not own an equation itself."""
    nw = Network(iterinfo=False)
    if h_range is not None:
        nw.h_range = h_range
    cp = machine("machine")
    c1 = Connection(Source("source"), "out1", cp, "in1", label="c1")
    c2 = Connection(cp, "out1", Sink("sink"), "in1", label="c2")
    nw.add_conns(c1, c2)
    c1.set_attr(fluid=fluid, T=T_in, p=p_in, m=1, **kwargs)
    c2.set_attr(p=p_out)
    cp.set_attr(eta_s=0.9)
    return nw


def _solve_and_collect_log(nw, caplog, **kwargs):
    with caplog.at_level(logging.DEBUG, logger="TESPyLogger"):
        nw.solve("design", **kwargs)
    return [msg for msg in caplog.messages if "did not converge" in msg]


class TestIterationGate:

    def test_mixture_bound_released_in_block_solve(self, caplog):
        """The outlet enthalpy of 9.3e6 J/kg is outside the range applied in
        the first iterations, the block must still converge on its own."""
        nw = _turbomachine_network(
            HYDROGEN_RICH, 1273, 10e5, 1e5, h_range=[1e3, 7e6]
        )
        failed_blocks = _solve_and_collect_log(nw, caplog)
        assert nw.status == 0
        assert failed_blocks == []
        assert nw.get_conn("c2").h.val_SI == pytest.approx(9.3336e6, rel=1e-4)


class TestMixtureEnthalpyBounds:

    @pytest.mark.parametrize("block_solve", [True, False])
    @pytest.mark.parametrize("machine, T_in, p_in, p_out, h_out", [
        (Turbine, 1273, 10e5, 1e5, 9.3336e6),
        (Compressor, 400, 1e5, 10e5, 1.03144e7),
    ])
    def test_enthalpy_above_former_default_range(
            self, caplog, block_solve, machine, T_in, p_in, p_out, h_out
        ):
        nw = _turbomachine_network(
            HYDROGEN_RICH, T_in, p_in, p_out, machine=machine
        )
        failed_blocks = _solve_and_collect_log(
            nw, caplog, block_solve=block_solve
        )
        assert nw.status == 0
        assert failed_blocks == []
        assert nw.get_conn("c2").h.val_SI == pytest.approx(h_out, rel=1e-4)

    @pytest.mark.parametrize("block_solve", [True, False])
    def test_enthalpy_below_former_default_range(self, caplog, block_solve):
        """The reference state of the pyromat ideal gases is the enthalpy of
        formation, the enthalpy of cold air is negative."""
        pytest.importorskip("pyromat")
        from tespy.tools.fluid_properties.wrappers import PyromatWrapper

        nw = _turbomachine_network(
            {"ig::N2": 0.79, "ig::O2": 0.21}, 400, 10e5, 1e5,
            fluid_engines={"N2": PyromatWrapper, "O2": PyromatWrapper}
        )
        failed_blocks = _solve_and_collect_log(
            nw, caplog, block_solve=block_solve
        )
        assert nw.status == 0
        assert failed_blocks == []
        assert nw.get_conn("c2").h.val_SI < 0

    @pytest.mark.parametrize("fluid, mixing_rule", [
        ({"water": 0.01, "air": 0.99}, "ideal-cond"),
        ({"N2": 0.72, "CO2": 0.15, "water": 0.08, "O2": 0.05}, "ideal-cond"),
        ({"R134a": 0.5, "air": 0.5}, "ideal"),
    ])
    def test_temperature_can_be_recovered_at_bounds(self, fluid, mixing_rule):
        nw = _turbomachine_network(
            fluid, 400, 10e5, 1e5, mixing_rule=mixing_rule
        )
        nw.solve("design", init_only=True)
        c2 = nw.get_conn("c2")
        p = c2.p.val_SI
        lower, upper = c2._property_bounds("h", nw)
        T_lower = T_mix_ph(p, lower, c2.fluid_data, c2.mixing_rule)
        T_upper = T_mix_ph(p, upper, c2.fluid_data, c2.mixing_rule)
        assert T_lower < T_upper
        assert lower == pytest.approx(
            h_mix_pT(p, T_lower, c2.fluid_data, c2.mixing_rule)
        )
        assert upper == pytest.approx(
            h_mix_pT(p, T_upper, c2.fluid_data, c2.mixing_rule)
        )

    def test_upper_bound_beyond_temperature_range_of_fuel(self):
        """A trace of methane must not keep the flue gas below the upper
        temperature limit of methane."""
        fluid = {
            "N2": 0.72, "CO2": 0.15, "water": 0.079, "O2": 0.05, "CH4": 0.001
        }
        nw = _turbomachine_network(
            fluid, 400, 10e5, 1e5, mixing_rule="ideal-cond"
        )
        nw.solve("design", init_only=True)
        c2 = nw.get_conn("c2")
        _, upper = c2._property_bounds("h", nw)
        assert c2.fluid.wrapper["CH4"]._T_max < 1000
        assert upper > h_mix_pT(
            c2.p.val_SI, 1773.15, c2.fluid_data, c2.mixing_rule
        )

    def test_bounds_independent_of_previous_pressure(self):
        nw = _turbomachine_network(
            {"N2": 0.72, "CO2": 0.15, "water": 0.08, "O2": 0.05},
            400, 10e5, 1e5, mixing_rule="ideal-cond"
        )
        nw.solve("design", init_only=True)
        c2 = nw.get_conn("c2")
        reference = c2._property_bounds("h", nw)
        c2.p.set_reference_val_SI(200)
        c2._property_bounds("h", nw)
        c2.p.set_reference_val_SI(1e5)
        assert c2._property_bounds("h", nw) == reference

    def test_user_range_tightens_bounds(self):
        nw = _turbomachine_network(
            HYDROGEN_RICH, 1273, 10e5, 1e5, h_range=[2e6, 7e6]
        )
        nw.solve("design", init_only=True)
        assert nw.get_conn("c2")._property_bounds("h", nw) == (2e6, 7e6)

    def test_user_range_does_not_widen_bounds(self):
        nw = _turbomachine_network(
            HYDROGEN_RICH, 1273, 10e5, 1e5, h_range=[-1e9, 1e9]
        )
        nw.solve("design", init_only=True)
        c2 = nw.get_conn("c2")
        lower, upper = c2._property_bounds("h", nw)
        assert -1e9 < lower < upper < 1e9


class TestEnthalpyRangeOfNetwork:

    def test_not_set_by_default(self):
        nw = Network()
        assert nw.h_range is None
        assert "h_range" not in nw._serialize()

    def test_user_range_exported_and_loaded(self):
        nw = _turbomachine_network(
            HYDROGEN_RICH, 800, 10e5, 1e5, h_range=[2e6, 8e6]
        )
        nw.solve("design")
        loaded = Network.from_dict(nw.export())
        assert list(loaded.h_range_SI) == pytest.approx([2e6, 8e6])

    def test_former_default_in_file_is_not_a_specification(self):
        nw = _turbomachine_network(
            HYDROGEN_RICH, 800, 10e5, 1e5, h_range=[1e3, 7e6]
        )
        nw.solve("design")
        loaded = Network.from_dict(nw.export())
        assert loaded.h_range is None


def _humid_air_heater(T_in, T_out, p=101325, **humidity):
    """The outlet enthalpy is a variable, the heat is specified to reach
    the outlet temperature without a change of the humidity ratio."""
    nw = Network(iterinfo=False)
    heater = SimpleHeatExchanger("heater")
    c1 = HAConnection(Source("source"), "out1", heater, "in1", label="c1")
    c2 = HAConnection(heater, "out1", Sink("sink"), "in1", label="c2")
    nw.add_conns(c1, c2)
    if "r" in humidity:
        w = HAPropsSI("W", "T", T_in, "P", p, "R", humidity["r"])
        c1.set_attr(
            r=humidity["r"], fluid_balance=True,
            fluid0={"air": 0.99, "water": 0.01}
        )
    else:
        w = humidity["w"]
        c1.set_attr(w=w)
    c1.set_attr(p=p, T=T_in, m=1)
    # water above saturation is not part of the humid air enthalpy
    try:
        w_out = min(w, HAPropsSI("W", "T", T_out, "P", p, "R", 1))
    except ValueError:
        # no saturated state above the boiling temperature
        w_out = w
    heater.set_attr(pr=1, Q=(
        HAPropsSI("H", "T", T_out, "P", p, "W", w_out)
        - HAPropsSI("H", "T", T_in, "P", p, "W", w)
    ))
    return nw


class TestHumidAirEnthalpyBounds:

    @pytest.mark.parametrize("T_in, T_out, humidity", [
        (293.15, 308.15, {"r": 0.5}),
        (318.15, 338.15, {"r": 1}),
        # below the former lower limit of saturated air at -50 °C
        (243.15, 213.15, {"w": 1e-5}),
        (243.15, 153.15, {"w": 1e-5}),
        # above the former upper limit of dry air at 300 °C
        (333.15, 353.15, {"r": 0.9}),
        (323.15, 453.15, {"w": 0.05}),
        (333.15, 523.15, {"w": 0.1}),
        (293.15, 603.15, {"w": 0.005}),
    ])
    def test_heating_and_cooling(self, T_in, T_out, humidity):
        nw = _humid_air_heater(T_in, T_out, **humidity)
        nw.solve("design")
        assert nw.status == 0
        assert nw.get_conn("c2").T.val_SI == pytest.approx(T_out, abs=1e-3)

    def test_bounds_use_humidity_of_stream(self):
        nw = _humid_air_heater(333.15, 523.15, w=0.1)
        nw.solve("design", init_only=True)
        c2 = nw.get_conn("c2")
        lower, upper = c2._property_bounds("h", nw)
        p = c2.p.val_SI
        assert lower < HAPropsSI("H", "T", 135, "P", p, "W", 0)
        assert upper > HAPropsSI("H", "T", 615, "P", p, "W", 0.1)

    def test_bounds_independent_of_previous_pressure(self):
        nw = _humid_air_heater(333.15, 523.15, w=0.1)
        nw.solve("design", init_only=True)
        c2 = nw.get_conn("c2")
        reference = c2._property_bounds("h", nw)
        c2.p.set_reference_val_SI(90e5)
        c2._property_bounds("h", nw)
        c2.p.set_reference_val_SI(101325)
        assert c2._property_bounds("h", nw) == reference


def test_humidity_ratio_without_saturation_limit():
    """Above the boiling temperature of water there is no saturation limit
    for the humidity ratio."""
    fluid_data = {
        "air": {"mass_fraction": 1 / 1.05},
        "water": {"mass_fraction": 0.05 / 1.05},
    }
    assert w_mix_pT_humidair(101325, 473.15, fluid_data) == pytest.approx(0.05)
