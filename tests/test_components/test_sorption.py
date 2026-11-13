# -*- coding: utf-8

"""Module for testing the sorption components.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_sorption.py

SPDX-License-Identifier: MIT
"""
from CoolProp.CoolProp import PropsSI
from pytest import approx
from pytest import mark

from tespy.components import Absorber
from tespy.components import CooledAbsorber
from tespy.components import Desorber
from tespy.components import HeatedDesorber
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.connections import SolutionConnection
from tespy.networks import Network

P_EV = PropsSI("P", "T", 278.15, "Q", 0, "water")
P_CD = PropsSI("P", "T", 318.15, "Q", 0, "water")


def _network():
    nw = Network()
    nw.iterinfo = False
    return nw


def _xi(c):
    return c.fluid.val["LiBr"]


def _sorption_heat(inlets, outlets):
    return (
        sum(c.m.val_SI * c.h.val_SI for c in outlets)
        - sum(c.m.val_SI * c.h.val_SI for c in inlets)
    )


def _assert_sorption_balances(solution_in, solution_out, vapor):
    if vapor.source is solution_in.target:
        inlets, outlets = [solution_in], [solution_out, vapor]
    else:
        inlets, outlets = [solution_in, vapor], [solution_out]
    m_in = sum(c.m.val_SI for c in inlets)
    m_out = sum(c.m.val_SI for c in outlets)
    assert m_in == approx(m_out, rel=1e-8)
    assert (
        solution_in.m.val_SI * _xi(solution_in)
        == approx(solution_out.m.val_SI * _xi(solution_out), rel=1e-8)
    )
    p = {c.p.val_SI for c in [solution_in, solution_out, vapor]}
    assert max(p) == approx(min(p), rel=1e-10)
    assert sum(solution_out.fluid.val.values()) == approx(1, abs=1e-10)
    return _sorption_heat(inlets, outlets)


def absorber_network():
    nw = _network()
    absorber = Absorber("absorber")
    poor = SolutionConnection(
        Source("poor solution"), "out1", absorber, "in1", label="poor"
    )
    vapor = Connection(Source("vapor"), "out1", absorber, "in2", label="vapor")
    rich = SolutionConnection(
        absorber, "out1", Sink("rich solution"), "in1", label="rich"
    )
    nw.add_conns(poor, vapor, rich)
    return nw, absorber, poor, vapor, rich


def solve_absorber_temperature_specification():
    nw, absorber, poor, vapor, rich = absorber_network()
    poor.set_attr(xi=0.62, T=320)
    vapor.set_attr(fluid={"H2O": 1}, p=P_EV, m=1, x=1)
    rich.set_attr(T=306.15)
    nw.solve("design")
    return nw, absorber, poor, vapor, rich


def test_absorber_temperature_specification():
    nw, absorber, poor, vapor, rich = solve_absorber_temperature_specification()
    nw.assert_convergence()
    Q = _assert_sorption_balances(poor, rich, vapor)
    assert absorber.Q.val_SI == approx(Q, rel=1e-6)
    assert _xi(rich) < _xi(poor)
    assert poor.m.val_SI == approx(
        vapor.m.val_SI * _xi(rich) / (_xi(poor) - _xi(rich)), rel=1e-8
    )
    assert absorber.Q.val_SI < 0


def test_absorber_enthalpy_specification():
    _, _, _, _, reference = solve_absorber_temperature_specification()

    nw, absorber, poor, vapor, rich = absorber_network()
    poor.set_attr(xi=0.62, T=320)
    vapor.set_attr(fluid={"H2O": 1}, p=P_EV, m=1, x=1)
    rich.set_attr(h=reference.h.val_SI)
    nw.solve("design")
    nw.assert_convergence()
    assert rich.T.val_SI == approx(306.15, abs=1e-4)
    assert _xi(rich) == approx(_xi(reference), abs=1e-8)


def test_absorber_poor_solution_mass_flow_specification():
    nw, absorber, poor, vapor, rich = absorber_network()
    poor.set_attr(xi=0.62, T=320, m=10)
    vapor.set_attr(fluid={"H2O": 1}, p=P_EV, x=1)
    rich.set_attr(T=306.15)
    nw.solve("design")
    nw.assert_convergence()
    _assert_sorption_balances(poor, rich, vapor)
    assert vapor.m.val_SI == approx(10 * (_xi(poor) / _xi(rich) - 1), rel=1e-8)


def test_absorber_rich_solution_composition_specification():
    nw, absorber, poor, vapor, rich = absorber_network()
    poor.set_attr(xi=0.62, T=320)
    vapor.set_attr(fluid={"H2O": 1}, m=1, x=1)
    rich.set_attr(T=306.15, fluid={"LiBr": 0.5, "H2O": 0.5})
    nw.solve("design")
    nw.assert_convergence()
    _assert_sorption_balances(poor, rich, vapor)
    # the pressure follows from the equilibrium at the solution outlet
    assert P_EV / 3 < vapor.p.val_SI < P_EV * 3


def test_absorber_heat_specification_recovers_mass_flow():
    nw, absorber, poor, vapor, rich = solve_absorber_temperature_specification()
    Q = absorber.Q.val_SI
    m_poor = poor.m.val_SI
    vapor.set_attr(m=None)
    absorber.set_attr(Q=Q * 1.5)
    nw.solve("design")
    nw.assert_convergence()
    assert vapor.m.val_SI == approx(1.5, rel=1e-6)
    assert poor.m.val_SI == approx(m_poor * 1.5, rel=1e-6)


def desorber_network(vapor_state="inlet-equilibrium"):
    nw = _network()
    desorber = Desorber("desorber", vapor_state=vapor_state)
    rich = SolutionConnection(
        Source("rich solution"), "out1", desorber, "in1", label="rich"
    )
    poor = SolutionConnection(
        desorber, "out1", Sink("poor solution"), "in1", label="poor"
    )
    vapor = Connection(desorber, "out2", Sink("vapor"), "in1", label="vapor")
    nw.add_conns(rich, poor, vapor)
    return nw, desorber, rich, poor, vapor


def solve_desorber_temperature_specification():
    nw, desorber, rich, poor, vapor = desorber_network()
    rich.set_attr(xi=0.55, p=P_CD, m=10, T=315)
    poor.set_attr(T=360)
    vapor.set_attr(fluid={"H2O": 1})
    nw.solve("design")
    return nw, desorber, rich, poor, vapor


def test_desorber_temperature_specification():
    nw, desorber, rich, poor, vapor = solve_desorber_temperature_specification()
    nw.assert_convergence()
    Q = _assert_sorption_balances(rich, poor, vapor)
    assert desorber.Q.val_SI == approx(Q, rel=1e-6)
    assert _xi(poor) > _xi(rich)
    assert desorber.Q.val_SI > 0
    # the vapor leaves in equilibrium with the entering solution
    assert vapor.T.val_SI == approx(rich.T_bubble.val_SI, abs=1e-6)
    assert vapor.T.val_SI > PropsSI("T", "P", P_CD, "Q", 1, "water")


def test_desorber_vapor_in_equilibrium_with_outlet():
    nw, desorber, rich, poor, vapor = desorber_network("outlet-equilibrium")
    rich.set_attr(xi=0.55, p=P_CD, m=10, T=315)
    poor.set_attr(T=360)
    vapor.set_attr(fluid={"H2O": 1})
    nw.solve("design")
    nw.assert_convergence()
    assert vapor.T.val_SI == approx(360, abs=1e-6)


def test_desorber_free_vapor_state():
    nw, desorber, rich, poor, vapor = desorber_network("free")
    rich.set_attr(xi=0.55, p=P_CD, m=10, T=315)
    poor.set_attr(T=360)
    vapor.set_attr(fluid={"H2O": 1}, x=1)
    nw.solve("design")
    nw.assert_convergence()
    Q = _assert_sorption_balances(rich, poor, vapor)
    assert desorber.Q.val_SI == approx(Q, rel=1e-6)
    assert vapor.T.val_SI == approx(
        PropsSI("T", "P", P_CD, "Q", 1, "water"), abs=1e-6
    )


def test_desorber_approach_to_equilibrium():
    nw, desorber, rich, poor, vapor = desorber_network()
    rich.set_attr(xi=0.55, p=P_CD, m=10, T=315)
    poor.set_attr(T=360)
    vapor.set_attr(fluid={"H2O": 1})
    desorber.set_attr(td_eq=2)
    nw.solve("design")
    nw.assert_convergence()
    assert poor.td_bubble.val_SI == approx(2, abs=1e-6)
    assert desorber.td_eq.val_SI == approx(2, abs=1e-6)


def test_desorber_enthalpy_specification():
    _, _, ref_rich, ref_poor, _ = solve_desorber_temperature_specification()

    nw, desorber, rich, poor, vapor = desorber_network()
    rich.set_attr(xi=0.55, p=P_CD, m=10, h=ref_rich.h.val_SI)
    poor.set_attr(h=ref_poor.h.val_SI)
    vapor.set_attr(fluid={"H2O": 1})
    nw.solve("design")
    nw.assert_convergence()
    assert poor.T.val_SI == approx(360, abs=1e-4)
    assert _xi(poor) == approx(_xi(ref_poor), abs=1e-8)


def test_desorber_poor_solution_mass_flow_specification():
    nw, desorber, rich, poor, vapor = solve_desorber_temperature_specification()
    poor.set_attr(T=None, m=9.5)
    nw.solve("design")
    nw.assert_convergence()
    _assert_sorption_balances(rich, poor, vapor)
    assert vapor.m.val_SI == approx(0.5, rel=1e-8)
    assert _xi(poor) == approx(0.55 * 10 / 9.5, rel=1e-8)


def test_desorber_poor_solution_composition_specification():
    nw, desorber, rich, poor, vapor = solve_desorber_temperature_specification()
    rich.set_attr(p=None)
    poor.set_attr(fluid={"LiBr": 0.6, "H2O": 0.4})
    nw.solve("design")
    nw.assert_convergence()
    _assert_sorption_balances(rich, poor, vapor)
    assert poor.m.val_SI == approx(10 * 0.55 / 0.6, rel=1e-8)


def cooled_absorber_network():
    nw = _network()
    absorber = CooledAbsorber("absorber")
    poor = SolutionConnection(
        Source("poor solution"), "out1", absorber, "in1", label="poor"
    )
    cool_in = Connection(
        Source("cooling water in"), "out1", absorber, "in2", label="cool_in"
    )
    vapor = Connection(Source("vapor"), "out1", absorber, "in3", label="vapor")
    rich = SolutionConnection(
        absorber, "out1", Sink("rich solution"), "in1", label="rich"
    )
    cool_out = Connection(
        absorber, "out2", Sink("cooling water out"), "in1", label="cool_out"
    )
    nw.add_conns(poor, cool_in, vapor, rich, cool_out)
    absorber.set_attr(dp2=0)
    poor.set_attr(xi=0.62, T=320, m=10)
    vapor.set_attr(fluid={"H2O": 1}, p=P_EV, m=1, x=1)
    cool_in.set_attr(fluid={"H2O": 1}, p=2e5, T=293.15)
    cool_out.set_attr(T=298.15)
    return nw, absorber, poor, vapor, rich, cool_in, cool_out


def test_cooled_absorber_specification_changes():
    nw, absorber, poor, vapor, rich, cool_in, cool_out = (
        cooled_absorber_network()
    )
    nw.solve("design")
    nw.assert_convergence()
    Q_sorption = _assert_sorption_balances(poor, rich, vapor)
    Q_cooling = _sorption_heat([cool_in], [cool_out])
    assert absorber.Q.val_SI == approx(Q_sorption, rel=1e-6)
    assert absorber.Q.val_SI == approx(-Q_cooling, rel=1e-6)
    ttd_l = absorber.ttd_l.val_SI
    m_vapor = vapor.m.val_SI

    vapor.set_attr(m=None)
    absorber.set_attr(ttd_l=ttd_l)
    nw.solve("design")
    nw.assert_convergence()
    assert vapor.m.val_SI == approx(m_vapor, rel=1e-5)
    Q = absorber.Q.val_SI

    absorber.set_attr(ttd_l=None, Q=Q)
    nw.solve("design")
    nw.assert_convergence()
    assert vapor.m.val_SI == approx(m_vapor, rel=1e-5)
    assert absorber.ttd_l.val_SI == approx(ttd_l, abs=1e-4)


def heated_desorber_network():
    nw = _network()
    desorber = HeatedDesorber("desorber")
    rich = SolutionConnection(
        Source("rich solution"), "out1", desorber, "in1", label="rich"
    )
    heat_in = Connection(
        Source("heating water in"), "out1", desorber, "in2", label="heat_in"
    )
    poor = SolutionConnection(
        desorber, "out1", Sink("poor solution"), "in1", label="poor"
    )
    heat_out = Connection(
        desorber, "out2", Sink("heating water out"), "in1", label="heat_out"
    )
    vapor = Connection(desorber, "out3", Sink("vapor"), "in1", label="vapor")
    nw.add_conns(rich, heat_in, poor, heat_out, vapor)
    desorber.set_attr(dp2=0)
    rich.set_attr(xi=0.55, p=P_CD, m=10, T=315)
    poor.set_attr(T=365)
    vapor.set_attr(fluid={"H2O": 1})
    heat_in.set_attr(fluid={"H2O": 1}, p=5e5, T=385)
    heat_out.set_attr(T=375)
    return nw, desorber, rich, poor, vapor, heat_in, heat_out


def test_heated_desorber_specification_changes():
    nw, desorber, rich, poor, vapor, heat_in, heat_out = (
        heated_desorber_network()
    )
    nw.solve("design")
    nw.assert_convergence()
    Q_sorption = _assert_sorption_balances(rich, poor, vapor)
    Q_heating = _sorption_heat([heat_in], [heat_out])
    assert desorber.Q.val_SI == approx(Q_heating, rel=1e-6)
    assert desorber.Q.val_SI == approx(-Q_sorption, rel=1e-6)
    assert heat_in.m.val_SI > 0
    assert desorber.ttd_u.val_SI == approx(20, abs=1e-6)
    m_heat = heat_in.m.val_SI

    poor.set_attr(T=None)
    desorber.set_attr(ttd_u=20)
    nw.solve("design")
    nw.assert_convergence()
    assert heat_in.m.val_SI == approx(m_heat, rel=1e-5)
    Q = desorber.Q.val_SI

    desorber.set_attr(ttd_u=None, Q=Q)
    nw.solve("design")
    nw.assert_convergence()
    assert poor.T.val_SI == approx(365, abs=1e-4)


@mark.parametrize(
    "solve",
    [
        solve_absorber_temperature_specification,
        solve_desorber_temperature_specification,
    ],
)
def test_solution_outlet_at_bubble_point(solve):
    _, component, *_ = solve()
    outlet = component.outl[0]
    assert outlet.td_bubble.val_SI == approx(0, abs=1e-6)
    assert outlet.T_bubble.val_SI == approx(outlet.T.val_SI, abs=1e-6)
