# -*- coding: utf-8

"""Module for testing components of type equilibrium separator.
This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_equilibrium_separator.py
SPDX-License-Identifier: MIT
"""
from pytest import approx
from pytest import fixture

from tespy.components import EquilibriumSeparator
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.networks import Network


@fixture()
def separator_network_setup():
    nw = Network(iterinfo=False)
    nw.units.set_defaults(**{
        "pressure": "bar", "pressure_difference": "bar",
        "temperature": "degC"
    })
    so = Source("inflow")
    sil = Sink("liquid outflow")
    sig = Sink("gas outflow")
    es = EquilibriumSeparator("separator")
    c1 = Connection(so, "out1", es, "in1", label="1")
    c2 = Connection(es, "out1", sil, "in1", label="2")
    c3 = Connection(es, "out2", sig, "in1", label="3")
    nw.add_conns(c1, c2, c3)
    yield nw


def energy_balance_residual(nw):
    c1 = nw.get_conn("1")
    c2 = nw.get_conn("2")
    c3 = nw.get_conn("3")
    return (
        c1.m.val_SI * c1.h.val_SI
        - c2.m.val_SI * c2.h.val_SI
        - c3.m.val_SI * c3.h.val_SI
    )


def test_liquid_inlet(separator_network_setup):
    nw = separator_network_setup
    c1, c2, c3 = nw.get_conn(["1", "2", "3"])
    c1.set_attr(fluid={"water": 1}, m=10, p=1.5, td_bubble=10)
    nw.solve("design")
    nw.assert_convergence()

    assert c3.m.val_SI == approx(0, abs=1e-9)
    assert c2.m.val_SI == approx(c1.m.val_SI)
    assert c2.h.val_SI == approx(c1.h.val_SI)
    assert c3.calc_Q() == approx(1)
    assert energy_balance_residual(nw) == approx(0, abs=1e-6)


def test_gas_inlet(separator_network_setup):
    nw = separator_network_setup
    c1, c2, c3 = nw.get_conn(["1", "2", "3"])
    c1.set_attr(fluid={"water": 1}, m=10, p=1.5, td_dew=10)
    nw.solve("design")
    nw.assert_convergence()

    assert c2.m.val_SI == approx(0, abs=1e-9)
    assert c3.m.val_SI == approx(c1.m.val_SI)
    assert c3.h.val_SI == approx(c1.h.val_SI)
    assert c2.calc_Q() == approx(0, abs=1e-9)
    assert energy_balance_residual(nw) == approx(0, abs=1e-6)


def test_two_phase_inlet_splits_by_quality(separator_network_setup):
    nw = separator_network_setup
    c1, c2, c3 = nw.get_conn(["1", "2", "3"])
    c1.set_attr(fluid={"water": 1}, m=10, p=1.5, x=0.6)
    nw.solve("design")
    nw.assert_convergence()

    assert c3.m.val_SI == approx(0.6 * c1.m.val_SI)
    assert c2.m.val_SI == approx(0.4 * c1.m.val_SI)
    assert c2.calc_Q() == approx(0, abs=1e-9)
    assert c3.calc_Q() == approx(1)
    assert energy_balance_residual(nw) == approx(0, abs=1e-6)


def test_saturated_liquid_inlet(separator_network_setup):
    nw = separator_network_setup
    c1, c2, c3 = nw.get_conn(["1", "2", "3"])
    c1.set_attr(fluid={"water": 1}, m=10, p=1.5, x=0)
    nw.solve("design")
    nw.assert_convergence()

    assert c3.m.val_SI == approx(0, abs=1e-9)
    assert c2.h.val_SI == approx(c1.h.val_SI)


def test_saturated_gas_inlet(separator_network_setup):
    nw = separator_network_setup
    c1, c2, c3 = nw.get_conn(["1", "2", "3"])
    c1.set_attr(fluid={"water": 1}, m=10, p=1.5, x=1)
    nw.solve("design")
    nw.assert_convergence()

    assert c2.m.val_SI == approx(0, abs=1e-9)
    assert c3.h.val_SI == approx(c1.h.val_SI)


def test_outlet_mass_flow_determines_inlet_state(separator_network_setup):
    nw = separator_network_setup
    c1, c2, c3 = nw.get_conn(["1", "2", "3"])
    c1.set_attr(fluid={"water": 1}, m=10, T=111.35)
    c3.set_attr(m=9.5)
    nw.solve("design")
    nw.assert_convergence()

    assert c1.calc_Q() == approx(0.95)
    assert c1.calc_T_sat() == approx(c1.T.val_SI)


def test_supercritical_inlet_does_not_converge(separator_network_setup):
    nw = separator_network_setup
    c1 = nw.get_conn("1")
    c1.set_attr(fluid={"water": 1}, m=10, p=250, T=400)
    nw.solve("design")
    assert not nw.converged


def test_phase_change_on_resolve(separator_network_setup):
    nw = separator_network_setup
    c1, c2, c3 = nw.get_conn(["1", "2", "3"])
    c1.set_attr(fluid={"water": 1}, m=10, p=1.5, x=0.6)
    nw.solve("design")
    nw.assert_convergence()

    c1.set_attr(x=None, td_bubble=10)
    nw.solve("design")
    nw.assert_convergence()
    assert c3.m.val_SI == approx(0, abs=1e-9)
    assert c2.h.val_SI == approx(c1.h.val_SI)

    c1.set_attr(td_bubble=None, td_dew=10)
    nw.solve("design")
    nw.assert_convergence()
    assert c2.m.val_SI == approx(0, abs=1e-9)
    assert c3.h.val_SI == approx(c1.h.val_SI)

    c1.set_attr(td_dew=None, x=0.3)
    nw.solve("design")
    nw.assert_convergence()
    assert c3.m.val_SI == approx(0.3 * c1.m.val_SI)
