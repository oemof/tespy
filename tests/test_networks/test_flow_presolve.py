# -*- coding: utf-8

"""Module for testing the flow variable presolve of the starting values.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_flow_presolve.py

SPDX-License-Identifier: MIT
"""
import math

from CoolProp.CoolProp import PropsSI
from pytest import approx

from tespy.components import HeatExchanger
from tespy.components import Merge
from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Splitter
from tespy.connections import Connection
from tespy.networks import Network
from tespy.tools.helpers import UserDefinedEquation


def _network():
    nw = Network(iterinfo=False)
    nw.units.set_defaults(
        temperature="degC", pressure="bar", pressure_difference="bar",
        temperature_difference="delta_degC"
    )
    return nw


def _single_r134a_connection(**specs):
    nw = _network()
    c = Connection(Source("source"), "out1", Sink("sink"), "in1", label="c1")
    nw.add_conns(c)
    c.set_attr(fluid={"R134a": 1}, m=0.1, **specs)
    nw.solve("design", init_only=True)
    return c


def test_superheating_specification_seeds_enthalpy_off_the_dome():
    # a superheated port must not start on the saturation line, the
    # specified distance to the dew point is the better starting value
    c = _single_r134a_connection(p0=9.5, td_dew=50)
    T_dew = PropsSI("T", "P", 9.5e5, "Q", 1, "R134a")
    expected = PropsSI("H", "P", 9.5e5, "T", T_dew + 50, "R134a")
    assert c.h.val_SI == approx(expected, rel=1e-3)


def test_subcooling_specification_seeds_enthalpy_off_the_dome():
    c = _single_r134a_connection(p0=9.5, td_bubble=16)
    T_bubble = PropsSI("T", "P", 9.5e5, "Q", 0, "R134a")
    expected = PropsSI("H", "P", 9.5e5, "T", T_bubble - 16, "R134a")
    assert c.h.val_SI == approx(expected, rel=1e-3)


def test_negative_dew_point_difference_seeds_subcooled_state():
    c = _single_r134a_connection(p0=9.5, td_dew=-10)
    T_dew = PropsSI("T", "P", 9.5e5, "Q", 1, "R134a")
    expected = PropsSI("H", "P", 9.5e5, "T", T_dew - 10, "R134a")
    assert c.h.val_SI == approx(expected, rel=1e-3)


def test_zero_dew_point_difference_seeds_next_to_the_dome():
    c = _single_r134a_connection(p0=9.5, td_dew=0)
    h_dew = PropsSI("H", "P", 9.5e5, "Q", 1, "R134a")
    assert c.h.val_SI == approx(h_dew + 1e3, rel=1e-6)


def test_gas_state_specification_seeds_next_to_the_dome():
    c = _single_r134a_connection(p=9.5, state="g")
    h_dew = PropsSI("H", "P", 9.5e5, "Q", 1, "R134a")
    assert c.h.val_SI == approx(h_dew + 1e3, rel=1e-6)


def test_seeded_enthalpy_on_the_specified_side_is_kept():
    # a user value in the superheated region is kept even though it does not
    # match the specified 50 K exactly
    c = _single_r134a_connection(p0=9.5, h0=430e3, td_dew=50)
    assert c.h.val_SI == approx(430e3, rel=1e-9)


def test_seeded_enthalpy_on_the_wrong_side_is_replaced_by_the_target():
    c = _single_r134a_connection(p0=9.5, h0=250e3, td_dew=50)
    T_dew = PropsSI("T", "P", 9.5e5, "Q", 1, "R134a")
    expected = PropsSI("H", "P", 9.5e5, "T", T_dew + 50, "R134a")
    assert c.h.val_SI == approx(expected, rel=1e-3)


def _bypass_coil_network():
    # air split into a contacted branch cooled on a coil and a bypass, the
    # contacted fraction follows from a user defined equation on the
    # temperatures: the split is invisible to the structural presolve and the
    # mixing energy balance is inconsistent with the guessed enthalpies.
    # Bypass interpretation of the contact effectiveness of the wet regime of
    # the CoolSolve example cooling_coil.eescode
    # (github.com/CoolProp/CoolSolve, examples), here with dry air
    nw = _network()
    splitter = Splitter("splitter", num_out=2)
    coil = SimpleHeatExchanger("coil")
    merge = Merge("merge", num_in=2)
    c1 = Connection(Source("inlet"), "out1", splitter, "in1", label="c1")
    c2 = Connection(splitter, "out1", coil, "in1", label="c2")
    c3 = Connection(coil, "out1", merge, "in1", label="c3")
    c4 = Connection(splitter, "out2", merge, "in2", label="c4")
    c5 = Connection(merge, "out1", Sink("outlet"), "in1", label="c5")
    nw.add_conns(c1, c2, c3, c4, c5)
    c1.set_attr(fluid={"air": 1}, m=1, T=30, p=1)
    c3.set_attr(T=5)

    def contact_func(ude):
        c_in, c_branch, c_contact = ude.conns
        epsilon = (
            (c_in.calc_T() - ude.params["T_target"])
            / (c_in.calc_T() - c_contact.calc_T())
        )
        return c_branch.m.val_SI - epsilon * c_in.m.val_SI

    def contact_dependents(ude):
        c_in, c_branch, c_contact = ude.conns
        return [c_in.m, c_branch.m, c_in.p, c_in.h, c_contact.p, c_contact.h]

    nw.add_ude(UserDefinedEquation(
        "contact effectiveness", contact_func, contact_dependents,
        conns=[c1, c2, c3], params={"T_target": 293.15}
    ))
    return nw, c1, c2, c4


def test_flow_presolve_keeps_branch_flows_feasible():
    nw, c1, c2, c4 = _bypass_coil_network()
    nw.solve("design", init_only=True)
    assert 0 < c2.m.val_SI < c1.m.val_SI
    assert 0 < c4.m.val_SI < c1.m.val_SI


def test_bypass_split_converges_without_starting_values():
    nw, c1, c2, c4 = _bypass_coil_network()
    nw.solve("design")
    nw.assert_convergence()
    assert c2.m.val_SI + c4.m.val_SI == approx(c1.m.val_SI)
    assert c2.m.val_SI == approx(0.4, rel=1e-6)


def test_zone_split_with_property_dependent_closures_converges():
    # refrigerant condensing in two zones in series, the air split in
    # parallel with the zone heat flows closed by user defined equations
    # on the states: only the mass balances determine the split in the
    # presolve, the energy balances of the zones carry flat enthalpy
    # differences at the start. Reduced two zone variant of the three zone
    # condenser below, with a constant UA per kg of air instead of the
    # published heat transfer correlations
    nw = _network()
    tp = HeatExchanger("two phase zone")
    sc = HeatExchanger("subcooling zone")
    splitter = Splitter("air split", num_out=2)
    merge = Merge("air mix", num_in=2)
    r1 = Connection(Source("refrigerant in"), "out1", tp, "in1", label="r1")
    r2 = Connection(tp, "out1", sc, "in1", label="r2")
    r3 = Connection(sc, "out1", Sink("refrigerant out"), "in1", label="r3")
    a1 = Connection(Source("air in"), "out1", splitter, "in1", label="a1")
    a2 = Connection(splitter, "out1", tp, "in2", label="a2")
    a3 = Connection(splitter, "out2", sc, "in2", label="a3")
    a4 = Connection(tp, "out2", merge, "in1", label="a4")
    a5 = Connection(sc, "out2", merge, "in2", label="a5")
    a6 = Connection(merge, "out1", Sink("air out"), "in1", label="a6")
    nw.add_conns(r1, r2, r3, a1, a2, a3, a4, a5, a6)
    # the condensing pressure is the unknown closed by the air split
    r1.set_attr(fluid={"R134a": 1}, m=0.1, x=1)
    r2.set_attr(x=0)
    r3.set_attr(td_bubble=5)
    a1.set_attr(fluid={"air": 1}, m=1.6, T=20, p=1.01325)
    tp.set_attr(dp1=0, dp2=0)
    sc.set_attr(dp1=0)

    def zone_func(ude):
        # effectiveness-NTU with the air heat capacity at the mean of the
        # zone inlet and outlet temperature
        r_in, r_out, a_in, a_out = ude.conns
        T_mean = 0.5 * (a_in.calc_T() + a_out.calc_T())
        cp_air = PropsSI("C", "T", T_mean, "P", a_in.p.val_SI, "air")
        C_air = a_in.m.val_SI * cp_air
        NTU = ude.params["UA_per_kg"] * a_in.m.val_SI / C_air
        epsilon = 1 - math.exp(-NTU)
        dT = r_in.calc_T() - a_in.calc_T()
        return (
            r_in.m.val_SI * (r_in.h.val_SI - r_out.h.val_SI)
            - epsilon * C_air * dT
        )

    def zone_dependents(ude):
        r_in, r_out, a_in, a_out = ude.conns
        return [
            r_in.m, r_in.p, r_in.h, r_out.p, r_out.h,
            a_in.m, a_in.p, a_in.h, a_out.p, a_out.h
        ]

    for label, conns in (("tp", [r1, r2, a2, a4]), ("sc", [r2, r3, a3, a5])):
        nw.add_ude(UserDefinedEquation(
            f"zone {label}", zone_func, zone_dependents, conns=conns,
            params={"UA_per_kg": 2500}
        ))
    nw.solve("design")
    nw.assert_convergence()
    assert a2.m.val_SI + a3.m.val_SI == approx(a1.m.val_SI)
    assert a2.m.val_SI > a3.m.val_SI > 0
    assert 6e5 < r1.p.val_SI < 9e5


def test_three_zone_condenser_converges_without_starting_values():
    # desuperheating, two phase and subcooling zone in series on the
    # refrigerant, the air split in parallel by the zone area fractions:
    # the condensing pressure follows from the closure of the split. The
    # zone closures are user defined counterflow eps-NTU relations with the
    # heat capacities at the zone mean temperatures, so the presolve sees
    # only the mass balances and the energy balances of the zones.
    # Parameters of the CoolSolve example condenser_3zones.eescode
    # (github.com/CoolProp/CoolSolve, examples) after Cuevas C. and Winandy
    # E. 2002, Simplified 3 zones modelling of a car air-conditioning
    # condenser, IIR Zero Leakage - Minimum Charge Conference, Stockholm, with
    # the resistance coefficients rounded
    nw = _network()
    zones = {
        "sh": HeatExchanger("desuperheating zone"),
        "tp": HeatExchanger("two phase zone"),
        "sc": HeatExchanger("subcooling zone"),
    }
    splitter = Splitter("air split", num_out=3)
    merge = Merge("air mix", num_in=3)
    r1 = Connection(Source("refrigerant in"), "out1", zones["sh"], "in1")
    r2 = Connection(zones["sh"], "out1", zones["tp"], "in1")
    r3 = Connection(zones["tp"], "out1", zones["sc"], "in1")
    r4 = Connection(zones["sc"], "out1", Sink("refrigerant out"), "in1")
    a1 = Connection(Source("air in"), "out1", splitter, "in1")
    a_in = [
        Connection(splitter, f"out{i + 1}", he, "in2")
        for i, he in enumerate(zones.values())
    ]
    a_out = [
        Connection(he, "out2", merge, f"in{i + 1}")
        for i, he in enumerate(zones.values())
    ]
    a6 = Connection(merge, "out1", Sink("air out"), "in1")
    nw.add_conns(r1, r2, r3, r4, a1, a6, *a_in, *a_out)

    r1.set_attr(fluid={"R134a": 1}, m=0.1, td_dew=50)
    r2.set_attr(x=1)
    r3.set_attr(x=0)
    r4.set_attr(td_bubble=16)
    a1.set_attr(fluid={"air": 1}, m=1.6, T=20, p=1.01325)
    zones["sh"].set_attr(dp1=0, dp2=0)
    zones["tp"].set_attr(dp1=0)
    zones["sc"].set_attr(dp1=0)
    resistance = {"sh": 3.4e-5, "tp": 7.4e-6, "sc": 7.8e-5}
    R_air = 2e-4 * 1.6 ** -0.5

    def zone_func(ude):
        r_in, r_out, c_in, c_out = ude.conns
        T_air = 0.5 * (c_in.calc_T() + c_out.calc_T())
        cp_air = PropsSI("C", "T", T_air, "P", c_in.p.val_SI, "air")
        C_air = c_in.m.val_SI * cp_air
        alpha = c_in.m.val_SI / 1.6
        UA = alpha / (R_air + ude.params["R"] * r_in.m.val_SI ** -0.8)
        if ude.params["two_phase"]:
            T_hot = r_in.calc_T_sat()
            C_min = C_air
            epsilon = 1 - math.exp(-UA / C_min)
        else:
            T_hot = r_in.calc_T()
            p_mean = 0.5 * (r_in.p.val_SI + r_out.p.val_SI)
            T_mean = 0.5 * (T_hot + r_out.calc_T())
            C_ref = r_in.m.val_SI * PropsSI("C", "T", T_mean, "P", p_mean, "R134a")
            C_min = min(C_ref, C_air)
            omega = C_min / max(C_ref, C_air)
            e = math.exp(-UA / C_min * (1 - omega))
            epsilon = (1 - e) / (1 - omega * e)
        return (
            r_in.m.val_SI * (r_in.h.val_SI - r_out.h.val_SI)
            - epsilon * C_min * (T_hot - c_in.calc_T())
        )

    def zone_dependents(ude):
        r_in, r_out, c_in, c_out = ude.conns
        return [
            r_in.m, r_in.p, r_in.h, r_out.p, r_out.h,
            c_in.m, c_in.p, c_in.h, c_out.p, c_out.h
        ]

    for (label, R), r_in, r_out, c_in, c_out in zip(
            resistance.items(), (r1, r2, r3), (r2, r3, r4), a_in, a_out
        ):
        nw.add_ude(UserDefinedEquation(
            f"zone {label}", zone_func, zone_dependents,
            conns=[r_in, r_out, c_in, c_out],
            params={"R": R, "two_phase": label == "tp"}
        ))

    nw.solve("design", max_iter=100)
    nw.assert_convergence()
    assert sum(c.m.val_SI for c in a_in) == approx(1.6)
    assert all(c.m.val_SI > 0 for c in a_in)
    assert r1.p.val_SI == approx(9.49e5, rel=2e-2)
