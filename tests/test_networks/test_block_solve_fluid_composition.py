# -*- coding: utf-8

"""Module for testing the structural analysis and the block-wise solution of
networks with variable fluid composition.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_block_solve_fluid_composition.py

SPDX-License-Identifier: MIT
"""
from pytest import approx
from pytest import fixture
from pytest import mark
from pytest import raises

from tespy.components import Compressor
from tespy.components import DiabaticCombustionChamber
from tespy.components import Drum
from tespy.components import HeatExchanger
from tespy.components import Merge
from tespy.components import MovingBoundaryHeatExchanger
from tespy.components import Pump
from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Splitter
from tespy.components import Turbine
from tespy.connections import Connection
from tespy.connections import HAConnection
from tespy.networks import Network
from tespy.solver.problem import Problem
from tespy.tools.data_containers import ComponentParameter as dc_cp
from tespy.tools.helpers import TESPyNetworkError
from tespy.tools.helpers import UserDefinedEquation

AIR = {"Ar": 0.0129, "N2": 0.7553, "CO2": 0.0004, "O2": 0.2314}
FUEL = {"CO2": 0.04, "CH4": 0.96}


def _network():
    nw = Network(iterinfo=False)
    nw.units.set_defaults(
        temperature="degC", pressure="bar", pressure_difference="bar"
    )
    return nw


def gas_turbine():
    """Open cycle gas turbine, the flue gas composition follows from the
    combustion and is a constant for everything downstream."""
    nw = _network()
    air = Source("air source")
    fuel = Source("fuel source")
    cp = Compressor("compressor")
    cc = DiabaticCombustionChamber("combustion chamber")
    tu = Turbine("turbine")
    cooler = SimpleHeatExchanger("flue gas cooler")
    stack = Sink("stack")

    c1 = Connection(air, "out1", cp, "in1", label="1")
    c2 = Connection(cp, "out1", cc, "in1", label="2")
    c3 = Connection(cc, "out1", tu, "in1", label="3")
    c4 = Connection(tu, "out1", cooler, "in1", label="4")
    c5 = Connection(cooler, "out1", stack, "in1", label="5")
    c6 = Connection(fuel, "out1", cc, "in2", label="6")
    nw.add_conns(c1, c2, c3, c4, c5, c6)

    cp.set_attr(eta_s=0.85, pr=15)
    cc.set_attr(pr=0.97, eta=0.98)
    tu.set_attr(eta_s=0.9)
    cooler.set_attr(pr=0.98)
    c1.set_attr(p=1, T=20, fluid=AIR)
    c3.set_attr(T=1200)
    c5.set_attr(p=1, T=150)
    c6.set_attr(m=1, p=16, T=20, fluid=FUEL)
    return nw


def gas_turbine_hrsg():
    """Gas turbine with a heat recovery steam generator of economizer, drum
    with evaporator circulation and superheater: the flue gas side has a
    variable composition, the water steam side a fixed one."""
    nw = _network()
    air = Source("air source")
    fuel = Source("fuel source")
    cp = Compressor("compressor")
    cc = DiabaticCombustionChamber("combustion chamber")
    tu = Turbine("turbine")
    sh = HeatExchanger("superheater")
    ev = HeatExchanger("evaporator")
    eco = HeatExchanger("economizer")
    stack = Sink("stack")
    feed = Source("feed water")
    pump = Pump("feed pump")
    drum = Drum("drum")
    st = Turbine("steam turbine")
    exhaust = Sink("steam exhaust")

    c1 = Connection(air, "out1", cp, "in1", label="1")
    c2 = Connection(cp, "out1", cc, "in1", label="2")
    c3 = Connection(cc, "out1", tu, "in1", label="3")
    c4 = Connection(tu, "out1", sh, "in1", label="4")
    c5 = Connection(sh, "out1", ev, "in1", label="5")
    c6 = Connection(ev, "out1", eco, "in1", label="6")
    c7 = Connection(eco, "out1", stack, "in1", label="7")
    c8 = Connection(fuel, "out1", cc, "in2", label="8")
    nw.add_conns(c1, c2, c3, c4, c5, c6, c7, c8)

    w1 = Connection(feed, "out1", pump, "in1", label="w1")
    w2 = Connection(pump, "out1", eco, "in2", label="w2")
    w3 = Connection(eco, "out2", drum, "in1", label="w3")
    w4 = Connection(drum, "out1", ev, "in2", label="w4")
    w5 = Connection(ev, "out2", drum, "in2", label="w5")
    w6 = Connection(drum, "out2", sh, "in2", label="w6")
    w7 = Connection(sh, "out2", st, "in1", label="w7")
    w8 = Connection(st, "out1", exhaust, "in1", label="w8")
    nw.add_conns(w1, w2, w3, w4, w5, w6, w7, w8)

    cp.set_attr(eta_s=0.85, pr=15)
    cc.set_attr(pr=0.97, eta=0.98)
    tu.set_attr(eta_s=0.9)
    sh.set_attr(pr1=0.99, pr2=0.97)
    ev.set_attr(pr1=0.99, ttd_l=10)
    eco.set_attr(pr1=0.99, pr2=0.98)
    pump.set_attr(eta_s=0.75)
    st.set_attr(eta_s=0.85)
    c1.set_attr(p=1, T=20, fluid=AIR)
    c3.set_attr(T=1200)
    c7.set_attr(p=1)
    c8.set_attr(m=1, p=16, T=20, fluid=FUEL)
    w1.set_attr(p=0.1, x=0, fluid={"water": 1})
    w2.set_attr(p=45)
    w3.set_attr(td_bubble=5)
    w5.set_attr(x=0.5)
    w7.set_attr(T=450)
    w8.set_attr(p=0.1)
    return nw


def gas_turbine_oxygen():
    """The oxygen fraction of the flue gas determines the air mass flow."""
    nw = gas_turbine()
    nw.get_conn("3").set_attr(T=None, fluid={"O2": 0.14})
    return nw


def merge_chain(cooler_class=None):
    """Gas streams of different composition are mixed in two stages with a
    compression and a cooling in between."""
    nw = _network()
    s1 = Source("nitrogen")
    s2 = Source("oxygen")
    s3 = Source("carbon dioxide")
    m1 = Merge("merge 1", num_in=2)
    m2 = Merge("merge 2", num_in=2)
    cp = Compressor("compressor")
    if cooler_class is None:
        cooler = SimpleHeatExchanger("cooler")
    else:
        cooler = cooler_class("cooler")
    tu = Turbine("expander")
    sink = Sink("sink")

    c1 = Connection(s1, "out1", m1, "in1", label="1")
    c2 = Connection(s2, "out1", m1, "in2", label="2")
    c3 = Connection(m1, "out1", cp, "in1", label="3")
    c4 = Connection(cp, "out1", cooler, "in1", label="4")
    c5 = Connection(cooler, "out1", m2, "in1", label="5")
    c6 = Connection(s3, "out1", m2, "in2", label="6")
    c7 = Connection(m2, "out1", tu, "in1", label="7")
    c8 = Connection(tu, "out1", sink, "in1", label="8")
    nw.add_conns(c1, c2, c3, c4, c5, c6, c7, c8)

    c1.set_attr(m=3, p=1, T=20, fluid={"N2": 1})
    c2.set_attr(m=1, T=40, fluid={"O2": 1})
    c6.set_attr(m=0.5, T=200, fluid={"CO2": 1})
    cp.set_attr(eta_s=0.8, pr=8)
    cooler.set_attr(pr=0.98)
    if cooler_class is None:
        c5.set_attr(T=80)
    else:
        cooler.set_attr(T_out=80 + 273.15)
    tu.set_attr(eta_s=0.85)
    c8.set_attr(p=1)
    return nw


def recirculation():
    """Flue gas recirculation, the composition is part of a recycle and
    shares one block with the fluid property equations."""
    nw = _network()
    air = Source("air source")
    fuel = Source("fuel source")
    mix = Merge("recirculation merge", num_in=2)
    cp = Compressor("compressor")
    cc = DiabaticCombustionChamber("combustion chamber")
    tu = Turbine("turbine")
    cooler = SimpleHeatExchanger("flue gas cooler")
    split = Splitter("recirculation splitter", num_out=2)
    stack = Sink("stack")

    c0 = Connection(air, "out1", mix, "in1", label="0")
    c1 = Connection(mix, "out1", cp, "in1", label="1")
    c2 = Connection(cp, "out1", cc, "in1", label="2")
    c3 = Connection(cc, "out1", tu, "in1", label="3")
    c4 = Connection(tu, "out1", cooler, "in1", label="4")
    c5 = Connection(cooler, "out1", split, "in1", label="5")
    c6 = Connection(split, "out1", stack, "in1", label="6")
    c7 = Connection(split, "out2", mix, "in2", label="7")
    c8 = Connection(fuel, "out1", cc, "in2", label="8")
    nw.add_conns(c0, c1, c2, c3, c4, c5, c6, c7, c8)

    cp.set_attr(eta_s=0.85, pr=15)
    cc.set_attr(pr=0.97, eta=0.98)
    tu.set_attr(eta_s=0.9)
    cooler.set_attr(pr=0.98)
    c0.set_attr(p=1, T=20, fluid=AIR)
    c3.set_attr(T=1200)
    c5.set_attr(T=60)
    c7.set_attr(m=10)
    c8.set_attr(m=1, p=16, T=20, fluid=FUEL)
    return nw


def _humid_air_cooler():
    nw = _network()
    air_in = Source("air source")
    air_out = Sink("air sink")
    ref_in = Source("refrigerant source")
    ref_out = Sink("refrigerant sink")
    hex = MovingBoundaryHeatExchanger("air cooler")

    a1 = Connection(ref_in, "out1", hex, "in2", label="a1")
    a2 = Connection(hex, "out2", ref_out, "in1", label="a2")
    c1 = HAConnection(air_in, "out1", hex, "in1", label="c1")
    c2 = HAConnection(hex, "out1", air_out, "in1", label="c2")
    nw.add_conns(a1, a2, c1, c2)

    a1.set_attr(fluid={"NH3": 1}, m=2, x=0.3)
    a2.set_attr(td_dew=5, T_dew=-20)
    c2.set_attr(T=5)
    hex.set_attr(dp1=0, dp2=0)
    return nw


def humid_air_condensing():
    """Relative humidity and temperature at the inlet determine the
    composition, the air leaves the cooler saturated."""
    nw = _humid_air_cooler()
    nw.get_conn("c1").set_attr(
        p=1, T=15, r=0.9, fluid0={"air": 0.9, "water": 0.1},
        fluid_balance=True
    )
    return nw


def humid_air_unsaturated():
    """As :code:`humid_air_condensing`, but the air leaves above its dew
    point, so the outlet temperature depends on the composition."""
    nw = humid_air_condensing()
    nw.get_conn("c2").set_attr(T=14)
    return nw


def humid_air_mass_flow():
    """Relative humidity and mass flow at the inlet: the inlet temperature,
    and with it the composition, couples to the heat exchanger."""
    nw = _humid_air_cooler()
    nw.get_conn("c1").set_attr(
        p=1, m=90, r=0.9, fluid0={"air": 0.9, "water": 0.1},
        fluid_balance=True
    )
    return nw


class LegacyCooler(SimpleHeatExchanger):
    """Custom component as written without any knowledge of the composition:
    the outlet temperature equation evaluates a fluid property and declares
    pressure and enthalpy only."""

    def get_parameters(self):
        parameters = super().get_parameters()
        parameters["T_out"] = dc_cp(
            num_eq_sets=1,
            func=self.T_out_func,
            dependents=self.T_out_dependents,
            description="outlet temperature in K"
        )
        return parameters

    def T_out_func(self):
        return self.outl[0].calc_T() - self.T_out.val

    def T_out_dependents(self):
        return [self.outl[0].p, self.outl[0].h]


MODELS = {
    "gas_turbine": gas_turbine,
    "gas_turbine_hrsg": gas_turbine_hrsg,
    "gas_turbine_oxygen": gas_turbine_oxygen,
    "merge_chain": merge_chain,
    "humid_air_condensing": humid_air_condensing,
    "humid_air_unsaturated": humid_air_unsaturated,
    "humid_air_mass_flow": humid_air_mass_flow,
}


def _assert_same_results(nw, reference, rel=1e-5):
    for c_ref in reference.conns["object"]:
        if not hasattr(c_ref, "fluid"):
            continue
        c = nw.get_conn(c_ref.label)
        for prop in ("m", "p", "h"):
            assert c.get_attr(prop).val_SI == approx(
                c_ref.get_attr(prop).val_SI, rel=rel
            ), f"{c_ref.label} ({prop})"
        for fluid, value in c_ref.fluid.val.items():
            assert c.fluid.val[fluid] == approx(value, rel=1e-6, abs=1e-9), (
                f"{c_ref.label} ({fluid})"
            )


@mark.parametrize("build", MODELS.values(), ids=MODELS.keys())
def test_block_solve_matches_simultaneous_solve(build):
    reference = build()
    reference.solve("design", block_solve=False)
    reference.assert_convergence()

    nw = build()
    nw.solve("design")
    nw.assert_convergence()
    _assert_same_results(nw, reference)


def test_custom_component_without_composition_dependents():
    nw = merge_chain(cooler_class=LegacyCooler)
    nw.solve("design")
    nw.assert_convergence()
    assert nw.get_conn("5").T.val == approx(80)

    reference = merge_chain()
    reference.solve("design")
    _assert_same_results(nw, reference)


def test_structural_analysis_underdetermined_gas_turbine():
    nw = gas_turbine()
    nw.get_conn("3").set_attr(T=None)
    with raises(TESPyNetworkError):
        nw.solve("design")
    assert nw.status == 11

    analysis = nw.get_structural_analysis()
    assert [part["kind"] for part in analysis] == ["underdetermined"]


def test_structural_analysis_overdetermined_gas_turbine():
    nw = gas_turbine()
    nw.get_conn("3").set_attr(fluid={"O2": 0.14})
    with raises(TESPyNetworkError):
        nw.solve("design")
    assert nw.status == 12

    analysis = nw.get_structural_analysis()
    assert [part["kind"] for part in analysis] == ["overdetermined"]


@fixture
def fallbacks(monkeypatch):
    """Simultaneous solutions of the full system during the test, the
    escalation of the block-wise solution."""
    calls = []
    solve_simultaneous = Problem._solve_simultaneous

    def record(problem, *args, **kwargs):
        calls.append(problem)
        return solve_simultaneous(problem, *args, **kwargs)

    monkeypatch.setattr(Problem, "_solve_simultaneous", record)
    return calls


def _assert_solved_block_wise(nw, fallbacks):
    assert fallbacks == []
    assert all(block.status == 0 for block in nw.problem.decompose().blocks)


@mark.parametrize("build", MODELS.values(), ids=MODELS.keys())
def test_block_solve_fluid_composition_matches_simultaneous_solve(build, fallbacks):
    reference = build()
    reference.solve("design", block_solve=False)
    reference.assert_convergence()
    fallbacks.clear()

    nw = build()
    nw.solve("design", block_solve_fluid_composition=True)
    nw.assert_convergence()
    _assert_same_results(nw, reference)
    _assert_solved_block_wise(nw, fallbacks)


def test_recirculation_block_solve_fluid_composition(fallbacks):
    nw = recirculation()
    nw.solve("design", block_solve_fluid_composition=True)
    nw.assert_convergence()
    _assert_solved_block_wise(nw, fallbacks)

    air, mixed, recirculated = (nw.get_conn(label) for label in "017")
    stack, fuel = nw.get_conn("6"), nw.get_conn("8")
    assert stack.m.val_SI == approx(air.m.val_SI + fuel.m.val_SI)
    assert nw.get_conn("3").T.val == approx(1200)
    assert nw.get_conn("5").T.val == approx(60)
    # the recycle closes: the compressor inlet is the mixture of the fresh
    # air with the flue gas, which differs from the air in its composition
    assert recirculated.fluid.val["CO2"] > 10 * air.fluid.val["CO2"]
    for fluid in mixed.fluid.val:
        assert mixed.m.val_SI * mixed.fluid.val[fluid] == approx(
            air.m.val_SI * air.fluid.val.get(fluid, 0)
            + recirculated.m.val_SI * recirculated.fluid.val[fluid],
            rel=1e-6, abs=1e-9
        )


def test_offdesign_block_solve_fluid_composition(tmp_path, fallbacks):
    results = []
    for kwargs in ({"block_solve_fluid_composition": True}, {"block_solve": False}):
        nw = gas_turbine_hrsg()
        nw.solve("design", **kwargs)
        nw.assert_convergence()
        path = f"{tmp_path}/design_{len(results)}.json"
        nw.save(path)

        for label in ("superheater", "evaporator", "economizer"):
            nw.get_comp(label).set_attr(offdesign=["UA"])
        nw.get_comp("evaporator").set_attr(design=["ttd_l"])
        nw.get_conn("w3").set_attr(design=["td_bubble"])
        nw.get_conn("w7").set_attr(design=["T"])
        nw.get_conn("8").set_attr(m=0.9)
        if not results:
            fallbacks.clear()
        nw.solve("offdesign", design_path=path, **kwargs)
        nw.assert_convergence()
        if not results:
            _assert_solved_block_wise(nw, fallbacks)
        results.append(nw)

    nw, reference = results
    _assert_same_results(nw, reference)
    assert nw.get_conn("w7").T.val != approx(450)


def test_pause_and_continue_block_solve_fluid_composition():
    reference = gas_turbine()
    reference.solve("design", block_solve=False)

    nw = gas_turbine()
    nw.solve(
        "design", max_iter=2, pause_on_block_failure=True,
        block_solve_fluid_composition=True
    )
    assert nw.status == 20
    nw.solve_continue(max_iter=50)
    nw.assert_convergence()
    assert nw.problem._composition_derivatives
    _assert_same_results(nw, reference)


def test_inferred_composition_dependence():
    nw = gas_turbine_hrsg()
    nw.solve("design", block_solve_fluid_composition=True)
    problem = nw.problem

    composition_columns = {
        col for col, data in problem.variables_dict.items()
        if data["variable"] == "fluid"
    }
    inferred = {
        problem._equation_lookup[eq_num][0]
        for eq_num in problem._inferred_incidence
    }
    # the property equations of the flue gas side read the composition,
    # the water steam side, the pure mass and fluid balances and the
    # energy balances in enthalpy do not
    assert {"turbine", "evaporator", "3"} <= inferred
    assert {
        "feed pump", "steam turbine", "drum", "superheater", "economizer"
    }.isdisjoint(inferred)
    for eq_num, columns in problem._inferred_incidence.items():
        assert set(columns) <= composition_columns
        assert set(columns).isdisjoint(problem._declared_incidence[eq_num])
        assert set(columns) <= set(problem.get_incidence()[eq_num])

    # the flow presolve works on the declared incidence
    for eq_num in problem._inferred_incidence:
        assert composition_columns.isdisjoint(
            problem._declared_incidence[eq_num]
        ) or problem._equation_lookup[eq_num][0] == "combustion chamber"

    # the water steam side does not depend on the composition, only the
    # flows determined by the energy balances of the hrsg follow it
    blocks = problem.decompose().blocks
    composition_blocks = [
        block for block in blocks
        if composition_columns & set(block.variables)
    ]
    assert len(composition_blocks) == 1
    for block in blocks:
        labels = {label for label, _ in block.equation_labels}
        if labels & {"feed pump", "steam turbine", "drum"}:
            assert composition_columns.isdisjoint(block.variables)
        if labels & {"superheater", "evaporator", "economizer"}:
            assert block.id >= composition_blocks[0].id


def test_no_inference_without_variable_composition():
    nw = Network(iterinfo=False)
    so = Source("source")
    cp = Compressor("compressor")
    si = Sink("sink")
    c1 = Connection(so, "out1", cp, "in1", label="1")
    c2 = Connection(cp, "out1", si, "in1", label="2")
    nw.add_conns(c1, c2)
    c1.set_attr(fluid={"N2": 0.7, "O2": 0.3}, m=1, p=1e5, T=293.15)
    cp.set_attr(pr=3, eta_s=0.8)
    nw.solve("design", block_solve_fluid_composition=True)
    nw.assert_convergence()
    problem = nw.problem
    assert not problem._variable_composition()
    assert problem._inferred_incidence == {}
    assert problem._declared_incidence is problem.get_incidence()


def test_custom_component_block_solve_fluid_composition():
    nw = merge_chain(cooler_class=LegacyCooler)
    nw.solve("design", block_solve_fluid_composition=True)
    nw.assert_convergence()
    assert nw.get_conn("5").T.val == approx(80)

    inferred = [
        nw.problem._equation_lookup[eq_num]
        for eq_num in nw.problem._inferred_incidence
    ]
    assert ("cooler", ("T_out", 0)) in inferred


class InletDeclaredCooler(LegacyCooler):
    """Wrong declaration: the outlet temperature equation lists the inlet
    state instead of the outlet state."""

    def T_out_dependents(self):
        return [self.inl[0].p, self.inl[0].h]


@mark.parametrize(
    "build",
    list(MODELS.values()) + [recirculation],
    ids=list(MODELS.keys()) + ["recirculation"]
)
def test_dependents_check_models(build):
    nw = build()
    nw.solve("design", block_solve_fluid_composition=True)
    nw.presolve("design")
    check = nw.get_dependents_check()
    assert check["missing"] == {}
    assert check["not_probed"] == []
    nw.print_dependents_check()


def test_dependents_check_detects_missing_dependents():
    nw = merge_chain(cooler_class=LegacyCooler)
    nw.solve("design", block_solve_fluid_composition=True)
    reference = nw.get_conn("5").h.val_SI

    nw = merge_chain(cooler_class=InletDeclaredCooler)
    nw.get_conn("5").set_attr(h0=reference)
    nw.presolve("design")
    check = nw.get_dependents_check()
    # inlet and outlet share the composition of their branch, so the inlet
    # state brings in the mass fractions by inference, only the outlet
    # enthalpy is missing
    missing = check["missing"][("cooler", ("T_out", 0))]
    assert len(missing) == 1
    assert missing[0].startswith("h")


def test_dependents_check_requires_loaded_variables():
    nw = gas_turbine()
    with raises(TESPyNetworkError):
        nw.get_dependents_check()
    nw.solve("design")
    with raises(TESPyNetworkError):
        nw.get_dependents_check()


class ExplicitCooler(LegacyCooler):
    """The outlet temperature equation lists the mass fractions of its
    outlet itself and switches the inference off."""

    def get_parameters(self):
        parameters = super().get_parameters()
        parameters["T_out"].infer_composition = False
        return parameters

    def T_out_dependents(self):
        outlet = self.outl[0]
        return {
            "scalars": [[outlet.p, outlet.h]],
            "vectors": [{outlet.fluid: outlet.fluid.is_var}]
        }


class WronglyFlaggedCooler(LegacyCooler):
    """The inference is switched off without listing the mass fractions."""

    def get_parameters(self):
        parameters = super().get_parameters()
        parameters["T_out"].infer_composition = False
        return parameters


def _inferred_labels(nw):
    return {
        nw.problem._equation_lookup[eq_num]
        for eq_num in nw.problem._inferred_incidence
    }


def test_core_energy_balances_not_inferred():
    nw = merge_chain()
    nw.solve("design", block_solve_fluid_composition=True)
    inferred = _inferred_labels(nw)
    assert ("merge 1", ("energy_balance_constraints", 0)) not in inferred
    assert ("merge 2", ("energy_balance_constraints", 0)) not in inferred
    assert ("compressor", ("eta_s", 0)) in inferred


def test_combustion_energy_balance_declares_composition():
    nw = gas_turbine()
    nw.solve("design")
    problem = nw.problem
    assert ("combustion chamber", ("eta", 0)) not in _inferred_labels(nw)

    eq_num = next(
        eq_num for eq_num, label in problem._equation_lookup.items()
        if label == ("combustion chamber", ("eta", 0))
    )
    composition_columns = {
        col for col, data in problem.variables_dict.items()
        if data["variable"] == "fluid"
    }
    assert composition_columns
    assert composition_columns <= set(problem._declared_incidence[eq_num])
    for col in composition_columns:
        assert problem.jacobian[eq_num, col] != 0


def test_custom_component_explicit_composition_dependents(fallbacks):
    nw = merge_chain(cooler_class=ExplicitCooler)
    nw.solve("design", block_solve_fluid_composition=True)
    nw.assert_convergence()
    assert nw.get_conn("5").T.val == approx(80)
    assert ("cooler", ("T_out", 0)) not in _inferred_labels(nw)
    _assert_solved_block_wise(nw, fallbacks)

    nw.presolve("design")
    assert nw.get_dependents_check()["missing"] == {}


def test_custom_component_wrongly_flagged(fallbacks):
    nw = merge_chain(cooler_class=WronglyFlaggedCooler)
    nw.solve("design", block_solve_fluid_composition=True)
    nw.assert_convergence()
    assert nw.get_conn("5").T.val == approx(80)

    blocks = nw.problem.decompose().blocks
    wrong_order = bool(fallbacks)

    nw.presolve("design")
    missing = nw.get_dependents_check()["missing"]
    assert ("cooler", ("T_out", 0)) in missing

    # nothing orders the equation relative to the composition it reads. If
    # it is solved first, the residual check after the block phase detects
    # the wrong block order and the simultaneous solution recovers
    equation_block = next(
        block.id for block in blocks
        if ("cooler", ("T_out", 0)) in block.equation_labels
    )
    composition_block = next(
        block.id for block in blocks
        if set(missing[("cooler", ("T_out", 0))]) & set(block.variable_labels)
    )
    assert equation_block != composition_block
    assert wrong_order == (equation_block < composition_block)


def test_skip_fallbacks_at_failed_verification(fallbacks):
    reference = merge_chain(cooler_class=WronglyFlaggedCooler)
    reference.solve("design", block_solve_fluid_composition=True)
    reference.assert_convergence()
    # the model only is a case for the option, if the residual check fails
    assert fallbacks
    fallbacks.clear()

    nw = merge_chain(cooler_class=WronglyFlaggedCooler)
    nw.solve("design", block_solve_fluid_composition=True, skip_fallbacks=True)
    assert nw.status == 2
    assert fallbacks == []
    assert all(block.status == 0 for block in nw.problem.decompose().blocks)


def _outlet_temperature_ude(infer_composition, deriv=None):
    nw = merge_chain()
    c5 = nw.get_conn("5")
    c5.set_attr(T=None)

    def func(ude):
        return ude.conns[0].calc_T() - 80 - 273.15

    if infer_composition:
        def dependents(ude):
            return [ude.conns[0].p, ude.conns[0].h]
    else:
        def dependents(ude):
            c = ude.conns[0]
            return {
                "scalars": [[c.p, c.h]],
                "vectors": [{c.fluid: c.fluid.is_var}]
            }

    ude = UserDefinedEquation(
        "outlet temperature", func, dependents, deriv=deriv, conns=[c5],
        infer_composition=infer_composition
    )
    nw.add_ude(ude)
    return nw


@mark.parametrize("infer_composition", [True, False])
def test_user_defined_equation_composition(infer_composition, fallbacks):
    nw = _outlet_temperature_ude(infer_composition)
    nw.solve("design", block_solve_fluid_composition=True)
    nw.assert_convergence()
    assert nw.get_conn("5").T.val == approx(80)
    _assert_solved_block_wise(nw, fallbacks)

    inferred = ("outlet temperature", ("equation", 0)) in _inferred_labels(nw)
    assert inferred == infer_composition

    nw.presolve("design")
    assert nw.get_dependents_check()["missing"] == {}


def test_user_defined_equation_derivatives_with_inference(fallbacks):
    calls = []

    def deriv(increment_filter, k, dependents=None, ude=None):
        calls.append(k)
        c = ude.conns[0]
        ude.partial_derivative(c.p)
        ude.partial_derivative(c.h)

    nw = _outlet_temperature_ude(True, deriv=deriv)
    nw.solve("design", block_solve_fluid_composition=True)
    nw.assert_convergence()
    assert calls
    assert nw.get_conn("5").T.val == approx(80)
    _assert_solved_block_wise(nw, fallbacks)

    assert ("outlet temperature", ("equation", 0)) in _inferred_labels(nw)
