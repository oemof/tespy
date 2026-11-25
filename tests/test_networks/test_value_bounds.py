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

