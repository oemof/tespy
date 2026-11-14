# -*- coding: utf-8

"""Module for testing single-effect absorption chiller models.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_models/test_absorption_chillers.py

SPDX-License-Identifier: MIT
"""
import importlib.util

from CoolProp.CoolProp import PropsSI
from pytest import approx
from pytest import mark
from pytest import param

from tespy.components import Absorber
from tespy.components import CycleCloser
from tespy.components import Desorber
from tespy.components import HeatExchanger
from tespy.components import Pump
from tespy.components import SimpleHeatExchanger
from tespy.components import Valve
from tespy.connections import Connection
from tespy.connections import SolutionConnection
from tespy.networks import Network

T_EV = 278.15
T_ABS = 313.15
T_CD = 318.15
T_GEN = 403.15
TTD_SHX = 10
ETA_PUMP = 0.7


def libr_water_reference():
    """Sequential solution of the single-effect cycle with the pair functions.

    Same cycle as CoolSolve Library CSL-0043 apart from the solution heat
    exchanger (terminal temperature difference instead of an enthalpy
    effectiveness) and the solution pump (isentropic efficiency instead of no
    pump work).
    """
    from tespy.tools.fluid_properties.solutions import get_solution_pair

    pair = get_solution_pair("libr-water")
    p_ev = PropsSI("P", "T", T_EV, "Q", 0, "water")
    p_cd = PropsSI("P", "T", T_CD, "Q", 0, "water")
    w_weak = pair.w_bubble(p_ev, T_ABS)
    w_strong = pair.w_bubble(p_cd, T_GEN)

    m_refrigerant = 1
    m_strong = m_refrigerant * w_weak / (w_strong - w_weak)
    m_weak = m_strong + m_refrigerant

    h_absorber_out = pair.h_pT(p_ev, T_ABS, w_weak)
    h_pump_out_s = pair.h_pT(p_cd, T_ABS, w_weak)
    h_pump_out = h_absorber_out + (h_pump_out_s - h_absorber_out) / ETA_PUMP
    T_pump_out = pair.T_ph(p_cd, h_pump_out, w_weak)

    h_generator_out = pair.h_pT(p_cd, T_GEN, w_strong)
    h_shx_hot_out = pair.h_pT(p_cd, T_pump_out + TTD_SHX, w_strong)
    Q_shx = m_strong * (h_generator_out - h_shx_hot_out)
    h_generator_in = h_pump_out + Q_shx / m_weak

    h_vapor = PropsSI("H", "P", p_cd, "T", T_GEN, "water")
    h_condenser_out = PropsSI("H", "P", p_cd, "Q", 0, "water")
    h_evaporator_out = PropsSI("H", "P", p_ev, "Q", 1, "water")

    Q_gen = (
        m_refrigerant * h_vapor + m_strong * h_generator_out
        - m_weak * h_generator_in
    )
    Q_ev = m_refrigerant * (h_evaporator_out - h_condenser_out)
    Q_cd = m_refrigerant * (h_condenser_out - h_vapor)
    Q_abs = (
        m_weak * h_absorber_out - m_strong * h_shx_hot_out
        - m_refrigerant * h_evaporator_out
    )
    P_pump = m_weak * (h_pump_out - h_absorber_out)
    return {
        "p_ev": p_ev, "p_cd": p_cd, "w_weak": w_weak, "w_strong": w_strong,
        "m_weak": m_weak, "m_strong": m_strong, "Q_gen": Q_gen, "Q_ev": Q_ev,
        "Q_cd": Q_cd, "Q_abs": Q_abs, "Q_shx": Q_shx, "P_pump": P_pump,
        "h_generator_in": h_generator_in, "COP": Q_ev / Q_gen,
    }


def libr_water_chiller(**solve_kwargs):
    nw = Network()
    nw.iterinfo = False

    generator = Desorber("generator", vapor_state="outlet-equilibrium")
    condenser = SimpleHeatExchanger("condenser")
    refrigerant_valve = Valve("refrigerant valve")
    evaporator = SimpleHeatExchanger("evaporator")
    absorber = Absorber("absorber")
    closer = CycleCloser("solution cycle closer")
    pump = Pump("solution pump")
    shx = HeatExchanger("solution heat exchanger")
    solution_valve = Valve("solution valve")

    s0 = SolutionConnection(
        absorber, "out1", closer, "in1", label="s0", mixing_rule="libr-water"
    )
    s1 = SolutionConnection(closer, "out1", pump, "in1", label="s1")
    s2 = SolutionConnection(pump, "out1", shx, "in2", label="s2")
    s3 = SolutionConnection(shx, "out2", generator, "in1", label="s3")
    s4 = SolutionConnection(generator, "out1", shx, "in1", label="s4")
    s5 = SolutionConnection(shx, "out1", solution_valve, "in1", label="s5")
    s6 = SolutionConnection(solution_valve, "out1", absorber, "in1", label="s6")
    r1 = Connection(generator, "out2", condenser, "in1", label="r1")
    r2 = Connection(condenser, "out1", refrigerant_valve, "in1", label="r2")
    r3 = Connection(refrigerant_valve, "out1", evaporator, "in1", label="r3")
    r4 = Connection(evaporator, "out1", absorber, "in2", label="r4")
    nw.add_conns(s0, s1, s2, s3, s4, s5, s6, r1, r2, r3, r4)

    condenser.set_attr(dp=0)
    evaporator.set_attr(dp=0)
    pump.set_attr(eta_s=ETA_PUMP)
    shx.set_attr(dp1=0, dp2=0, ttd_l=TTD_SHX)

    r1.set_attr(fluid={"H2O": 1})
    r2.set_attr(x=0, T=T_CD)
    r4.set_attr(x=1, T=T_EV, m=1)
    s0.set_attr(T=T_ABS)
    # the cycle closer does not close the variable composition
    s1.set_attr(fluid_balance=True)
    s4.set_attr(T=T_GEN)

    nw.solve("design", **solve_kwargs)
    return nw, {
        "generator": generator, "condenser": condenser,
        "evaporator": evaporator, "absorber": absorber, "pump": pump,
        "shx": shx, "s0": s0, "s3": s3, "s4": s4, "r1": r1,
    }


def test_libr_water_chiller_against_sequential_solution():
    reference = libr_water_reference()
    nw, parts = libr_water_chiller()
    nw.assert_convergence()

    assert parts["s0"].p.val_SI == approx(reference["p_ev"], rel=1e-6)
    assert parts["s4"].p.val_SI == approx(reference["p_cd"], rel=1e-6)
    assert parts["s0"].xi.val_SI == approx(reference["w_weak"], abs=1e-7)
    assert parts["s4"].xi.val_SI == approx(reference["w_strong"], abs=1e-7)
    assert parts["s0"].m.val_SI == approx(reference["m_weak"], rel=1e-5)
    assert parts["s4"].m.val_SI == approx(reference["m_strong"], rel=1e-5)
    assert parts["s3"].h.val_SI == approx(
        reference["h_generator_in"], rel=1e-5
    )
    assert parts["r1"].T.val_SI == approx(T_GEN, abs=1e-5)

    assert parts["generator"].Q.val_SI == approx(reference["Q_gen"], rel=1e-5)
    assert parts["absorber"].Q.val_SI == approx(reference["Q_abs"], rel=1e-5)
    assert parts["condenser"].Q.val_SI == approx(reference["Q_cd"], rel=1e-5)
    assert parts["evaporator"].Q.val_SI == approx(reference["Q_ev"], rel=1e-5)
    assert parts["shx"].Q.val_SI == approx(-reference["Q_shx"], rel=1e-5)
    assert parts["pump"].P.val_SI == approx(reference["P_pump"], rel=1e-5)

    balance = (
        parts["generator"].Q.val_SI + parts["evaporator"].Q.val_SI
        + parts["pump"].P.val_SI + parts["absorber"].Q.val_SI
        + parts["condenser"].Q.val_SI
    )
    assert balance == approx(0, abs=1e-6 * parts["generator"].Q.val_SI)
    COP = parts["evaporator"].Q.val_SI / parts["generator"].Q.val_SI
    assert COP == approx(reference["COP"], rel=1e-5)
    # plausibility against CSL-0043 (ASHRAE fits, enthalpy effectiveness 0.75
    # of the solution heat exchanger, no pump work): COP 0.7466
    assert 0.65 < COP < 0.85


# reference: CoolSolve Library CSL-0249 (nh3h2o_absorption_chiller.sol)
CSL_0249 = {
    "P_low": 2.905681496458e+05, "P_high": 1.554537370181e+06,
    "x_strong": 4.172235440499e-01, "x_weak": 3.067414564222e-01,
    "x_generator_vapor": 9.034029548463e-01,
    "m_strong": 5.168811742633e-02, "m_weak": 4.344482519739e-02,
    "m_generator_vapor": 9.570939525832e-03, "m_reflux": 1.327647296888e-03,
    "m_refrigerant": 8.243292228944e-03,
    "Q_gen": 2.201922902116e+04, "Q_rect": 4.258665747408e+03,
    "Q_cd": 9.241020897867e+03, "Q_abs": 1.864857131076e+04,
    "Q_SHX": 1.195407105360e+04, "W_pump": 1.290289348697e+02,
    "COP": 4.541485076697e-01, "f_circ": 6.270324524568e+00,
    "T_3": 8.491217025590e+01, "T_6": 5.706949132980e+01,
    "T_8": 4.612317847674e+01, "T_14": 5.475790428811e+01,
}


def nh3_water_chiller(**solve_kwargs):
    from tespy.components import Merge
    from tespy.components import Rectifier
    from tespy.tools import UserDefinedEquation
    from tespy.tools.fluid_properties.solutions import get_solution_pair

    pair = get_solution_pair("nh3-water")
    x_ref = 0.9995
    p_low = pair.p_bubble(263.15, x_ref)

    nw = Network()
    nw.iterinfo = False

    absorber = Absorber("absorber")
    closer = CycleCloser("solution cycle closer")
    pump = Pump("solution pump")
    shx = HeatExchanger("solution heat exchanger")
    merge = Merge("reflux merge")
    generator = Desorber("generator", vapor_state="outlet-equilibrium")
    rectifier = Rectifier("rectifier")
    solution_valve = Valve("solution valve")
    condenser = SimpleHeatExchanger("condenser")
    rhx = HeatExchanger("refrigerant heat exchanger")
    refrigerant_valve = Valve("refrigerant valve")
    evaporator = SimpleHeatExchanger("evaporator")

    def sc(source, outlet, target, inlet, label, **kwargs):
        return SolutionConnection(
            source, outlet, target, inlet, label=label, **kwargs
        )

    c0 = sc(absorber, "out1", closer, "in1", "0", mixing_rule="nh3-water")
    c1 = sc(closer, "out1", pump, "in1", "1")
    c2 = sc(pump, "out1", shx, "in2", "2")
    c3 = sc(shx, "out2", merge, "in1", "3")
    c3m = sc(merge, "out1", generator, "in1", "3m")
    c4 = sc(generator, "out1", shx, "in1", "4")
    c5 = sc(shx, "out1", solution_valve, "in1", "5")
    c6 = sc(solution_valve, "out1", absorber, "in1", "6")
    c7 = sc(generator, "out2", rectifier, "in1", "7")
    c8 = sc(rectifier, "out1", condenser, "in1", "8")
    c9 = sc(rectifier, "out2", merge, "in2", "9")
    c10 = sc(condenser, "out1", rhx, "in1", "10")
    c11 = sc(rhx, "out1", refrigerant_valve, "in1", "11")
    c12 = sc(refrigerant_valve, "out1", evaporator, "in1", "12")
    c13 = sc(evaporator, "out1", rhx, "in2", "13")
    c14 = sc(rhx, "out2", absorber, "in2", "14")
    nw.add_conns(
        c0, c1, c2, c3, c3m, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14
    )

    for component in [condenser, evaporator]:
        component.set_attr(dp=0)
    for component in [shx, rhx]:
        component.set_attr(dp1=0, dp2=0)
    pump.set_attr(eta_s=0.6)
    evaporator.set_attr(Q=1e4)

    c0.set_attr(T=308.15, p=p_low)
    c1.set_attr(fluid_balance=True)
    c4.set_attr(T=393.15)
    c8.set_attr(xi=x_ref)
    c10.set_attr(x=0, T=313.15)
    c13.set_attr(x=0.99)

    def shx_effectiveness(ude):
        hot_in, cold_in, hot_out = ude.conns
        return hot_out.calc_T() - (
            hot_in.calc_T() - 0.7 * (hot_in.calc_T() - cold_in.calc_T())
        )

    def shx_dependents(ude):
        return {
            "scalars": [[v for c in ude.conns for v in (c.p, c.h)]],
            "vectors": [{c.fluid: {"NH3"} & c.fluid.is_var for c in ude.conns}],
        }

    nw.add_ude(
        UserDefinedEquation(
            "solution heat exchanger effectiveness", shx_effectiveness,
            dependents=shx_dependents, conns=[c4, c2, c5]
        ),
        UserDefinedEquation(
            "refrigerant heat exchanger effectiveness", shx_effectiveness,
            dependents=shx_dependents, conns=[c10, c13, c11]
        ),
    )
    nw.solve("design", **solve_kwargs)
    return nw, {
        "absorber": absorber, "pump": pump, "shx": shx,
        "generator": generator, "rectifier": rectifier,
        "condenser": condenser, "evaporator": evaporator, "rhx": rhx,
        "c1": c1, "c3": c3, "c4": c4, "c6": c6, "c7": c7, "c8": c8,
        "c9": c9, "c10": c10, "c14": c14,
    }


@mark.skipif(
    importlib.util.find_spec("teqp") is None, reason="teqp is not installed"
)
def test_nh3_water_chiller_against_coolsolve():
    nw, parts = nh3_water_chiller()
    nw.assert_convergence()
    ref = CSL_0249

    assert parts["c1"].p.val_SI == approx(ref["P_low"], rel=1e-5)
    assert parts["c4"].p.val_SI == approx(ref["P_high"], rel=1e-5)
    assert parts["c1"].xi.val_SI == approx(ref["x_strong"], abs=1e-5)
    assert parts["c4"].xi.val_SI == approx(ref["x_weak"], abs=1e-5)
    assert parts["c7"].xi.val_SI == approx(ref["x_generator_vapor"], abs=1e-5)
    assert parts["c9"].xi.val_SI == approx(ref["x_weak"], abs=1e-5)

    assert parts["c1"].m.val_SI == approx(ref["m_strong"], rel=1e-3)
    assert parts["c4"].m.val_SI == approx(ref["m_weak"], rel=1e-3)
    assert parts["c7"].m.val_SI == approx(ref["m_generator_vapor"], rel=1e-3)
    assert parts["c9"].m.val_SI == approx(ref["m_reflux"], rel=1e-3)
    assert parts["c8"].m.val_SI == approx(ref["m_refrigerant"], rel=1e-3)

    assert parts["c3"].T.val_SI - 273.15 == approx(ref["T_3"], abs=1e-2)
    assert parts["c6"].T.val_SI - 273.15 == approx(ref["T_6"], abs=1e-2)
    assert parts["c8"].T.val_SI - 273.15 == approx(ref["T_8"], abs=1e-2)
    assert parts["c14"].T.val_SI - 273.15 == approx(ref["T_14"], abs=1e-2)

    # generator and rectifier heat follow the CoolSolve split: the reflux
    # enters the generator
    assert parts["generator"].Q.val_SI == approx(ref["Q_gen"], rel=1e-3)
    assert parts["rectifier"].Q.val_SI == approx(-ref["Q_rect"], rel=1e-3)
    assert parts["condenser"].Q.val_SI == approx(-ref["Q_cd"], rel=1e-3)
    assert parts["absorber"].Q.val_SI == approx(-ref["Q_abs"], rel=1e-3)
    assert parts["shx"].Q.val_SI == approx(-ref["Q_SHX"], rel=1e-3)
    assert parts["pump"].P.val_SI == approx(ref["W_pump"], rel=1e-3)
    COP = parts["evaporator"].Q.val_SI / parts["generator"].Q.val_SI
    assert COP == approx(ref["COP"], rel=1e-3)
    f_circ = parts["c1"].m.val_SI / parts["c8"].m.val_SI
    assert f_circ == approx(ref["f_circ"], rel=1e-3)

    # the refrigerant heat exchanger of the reference is infeasible: the
    # effectiveness of 0.7 refers to the liquid side with the larger heat
    # capacity rate, the cold vapor leaves at 54.8 degC above the 40 degC of
    # the entering liquid; TESPy reproduces the state and flags it
    assert parts["c14"].T.val_SI > parts["c10"].T.val_SI
    assert parts["rhx"].eff_cold.val_SI > 1
    assert nw.status == 1


@mark.parametrize("builder", [
    libr_water_chiller,
    param(nh3_water_chiller, marks=mark.skipif(
        importlib.util.find_spec("teqp") is None,
        reason="teqp is not installed"
    )),
])
def test_chiller_block_decomposition(builder):
    nw, _ = builder(block_solve=False)
    reference = {
        c.label: (c.m.val_SI, c.p.val_SI, c.h.val_SI)
        for c in nw.conns["object"]
    }
    assert nw.get_structural_analysis() == []

    # the declared dependents of the sorption equations are exact
    nw.presolve("design")
    check = nw.get_dependents_check()
    assert check["missing"] == {}
    assert check["unused_declared"] == {}
    assert check["not_probed"] == []

    nw, _ = builder(block_solve_fluid_composition=True)
    assert nw.status in (0, 1)
    for c in nw.conns["object"]:
        assert (c.m.val_SI, c.p.val_SI, c.h.val_SI) == approx(
            reference[c.label], rel=1e-6
        )
