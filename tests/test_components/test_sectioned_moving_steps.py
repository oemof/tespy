# -*- coding: utf-8

"""Tests of the phase boundary steps of the sectioned heat exchangers with
an end state on the saturation line.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_sectioned_moving_steps.py

SPDX-License-Identifier: MIT
"""
from CoolProp.CoolProp import PropsSI as PSI
from pytest import approx
from pytest import mark

from tespy.components import SectionedHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.networks import Network


def _condenser():
    nw = Network(iterinfo=False)
    nw.units.set_defaults(temperature="degC", pressure="bar")
    hx = SectionedHeatExchanger("condenser")
    c1 = Connection(Source("hot in"), "out1", hx, "in1")
    c2 = Connection(hx, "out1", Sink("hot out"), "in1")
    c3 = Connection(Source("cold in"), "out1", hx, "in2")
    c4 = Connection(hx, "out2", Sink("cold out"), "in1")
    nw.add_conns(c1, c2, c3, c4)
    c1.set_attr(fluid={"NH3": 1}, m=1, td_dew=20, T=60)
    c2.set_attr(td_bubble=5)
    c3.set_attr(fluid={"water": 1}, p=1, T=30)
    c4.set_attr(T=40)
    hx.set_attr(dp1=0.5, dp2=0.1)
    nw.solve("design")
    nw.assert_convergence()
    return nw, hx, c1, c2


@mark.parametrize(
    "saturated_inlet, saturated_outlet, zone_phases",
    [
        # hot outlet on the bubble line: no liquid zone
        (False, True, [1, 2]),
        # hot inlet on the dew line: no gas zone
        (True, False, [0, 1]),
        # both ends on the saturation line: two-phase only
        (True, True, [1]),
    ]
)
def test_moving_steps_without_zero_width_zones(
        saturated_inlet, saturated_outlet, zone_phases
    ):
    """With pressure drop, an end state which is gas or liquid but on the
    saturation line within the tolerance does not have a gas or liquid
    zone: its step must be left out together with the zone."""
    nw, hx, c1, c2 = _condenser()
    assert hx._moving_steps(0)[1] == [0, 1, 2]

    # 1 J/kg next to the saturation line the end states are not two-phase,
    # the steps are identified in the liquid to gas branch
    if saturated_inlet:
        h_dew = PSI("H", "P", c1.p.val_SI, "Q", 1, "NH3")
        c1.set_attr(td_dew=None, T=None, p=c1.p.val, h=h_dew + 1)
    if saturated_outlet:
        h_bubble = PSI("H", "P", c2.p.val_SI, "Q", 0, "NH3")
        c2.set_attr(td_bubble=None, h=h_bubble - 1)

    nw.solve("design")
    nw.assert_convergence()

    assert c1.calc_phase() == "g"
    assert c2.calc_phase() == "l"

    steps, phases = hx._moving_steps(0)
    assert phases == zone_phases
    assert len(steps) == len(phases) + 1
    assert steps[0] == 0 and steps[-1] == 1
    assert all(0 < step < 1 for step in steps[1:-1])

    # every section is assigned to one of the zones
    assert set(hx.phase_hot_per_section.val_SI) == set(zone_phases)
    assert approx(sum(hx.Q_per_section.val_SI)) == -hx.Q.val_SI
