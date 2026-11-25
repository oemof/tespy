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
from tespy.tools.fluid_properties import h_mix_pT

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

    def test_bounds_follow_temperature_range_of_mixture(self):
        nw = _turbomachine_network(
            {"water": 0.01, "air": 0.99}, 400, 10e5, 1e5,
            mixing_rule="ideal-cond"
        )
        nw.solve("design", init_only=True)
        c2 = nw.get_conn("c2")
        lower, upper = c2._property_bounds("h", nw)
        T_low, T_high = c2._mixture_T_limits
        # the triple point of the water does not limit the humid air
        assert T_low < 273.15
        assert T_high > 1273.15
        assert lower < upper
        assert lower == pytest.approx(
            h_mix_pT(c2.p.val_SI, T_low, c2.fluid_data, c2.mixing_rule)
        )
        assert upper == pytest.approx(
            h_mix_pT(c2.p.val_SI, T_high, c2.fluid_data, c2.mixing_rule)
        )

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
