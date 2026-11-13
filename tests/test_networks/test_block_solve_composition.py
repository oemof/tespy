# -*- coding: utf-8

"""Module for testing the structural analysis and the block-wise solution of
networks with variable fluid composition.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_block_solve_composition.py

SPDX-License-Identifier: MIT
"""
from pytest import approx
from pytest import mark
from pytest import raises

from tespy.components import Compressor
from tespy.components import DiabaticCombustionChamber
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
from tespy.tools.data_containers import ComponentParameter as dc_cp
from tespy.tools.helpers import TESPyNetworkError

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
    """Gas turbine with a heat recovery steam generator: the flue gas side
    has a variable composition, the water steam side a fixed one."""
    nw = _network()
    air = Source("air source")
    fuel = Source("fuel source")
    cp = Compressor("compressor")
    cc = DiabaticCombustionChamber("combustion chamber")
    tu = Turbine("turbine")
    hrsg = HeatExchanger("hrsg")
    stack = Sink("stack")
    feed = Source("feed water")
    pump = Pump("feed pump")
    st = Turbine("steam turbine")
    exhaust = Sink("steam exhaust")

    c1 = Connection(air, "out1", cp, "in1", label="1")
    c2 = Connection(cp, "out1", cc, "in1", label="2")
    c3 = Connection(cc, "out1", tu, "in1", label="3")
    c4 = Connection(tu, "out1", hrsg, "in1", label="4")
    c5 = Connection(hrsg, "out1", stack, "in1", label="5")
    c6 = Connection(fuel, "out1", cc, "in2", label="6")
    w1 = Connection(feed, "out1", pump, "in1", label="w1")
    w2 = Connection(pump, "out1", hrsg, "in2", label="w2")
    w3 = Connection(hrsg, "out2", st, "in1", label="w3")
    w4 = Connection(st, "out1", exhaust, "in1", label="w4")
    nw.add_conns(c1, c2, c3, c4, c5, c6, w1, w2, w3, w4)

    cp.set_attr(eta_s=0.85, pr=15)
    cc.set_attr(pr=0.97, eta=0.98)
    tu.set_attr(eta_s=0.9)
    hrsg.set_attr(pr1=0.98, pr2=0.95)
    pump.set_attr(eta_s=0.75)
    st.set_attr(eta_s=0.85)
    c1.set_attr(p=1, T=20, fluid=AIR)
    c3.set_attr(T=1200)
    c5.set_attr(p=1, T=150)
    c6.set_attr(m=1, p=16, T=20, fluid=FUEL)
    w1.set_attr(p=0.1, x=0, fluid={"water": 1})
    w2.set_attr(p=40)
    w3.set_attr(T=450)
    w4.set_attr(p=0.1)
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


@mark.xfail(
    reason=(
        "the composition derivatives of the fluid property equations are "
        "missing in the jacobian, the recycle does not converge"
    ),
    strict=True
)
def test_recirculation_simultaneous_solve():
    nw = recirculation()
    nw.solve("design", block_solve=False)
    nw.assert_convergence()


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


def _fallback_messages(caplog):
    return [
        record.getMessage() for record in caplog.records
        if "simultaneous" in record.getMessage()
        or "did not converge" in record.getMessage()
    ]


@mark.parametrize("build", MODELS.values(), ids=MODELS.keys())
def test_block_solve_composition_matches_simultaneous_solve(build, caplog):
    reference = build()
    reference.solve("design", block_solve=False)
    reference.assert_convergence()

    nw = build()
    with caplog.at_level("INFO", logger="TESPyLogger"):
        nw.solve("design", block_solve_composition=True)
    nw.assert_convergence()
    _assert_same_results(nw, reference)
    assert _fallback_messages(caplog) == []


def test_recirculation_block_solve_composition():
    nw = recirculation()
    nw.solve("design", block_solve_composition=True)
    nw.assert_convergence()

    stack, recirculated = nw.get_conn("6"), nw.get_conn("7")
    air, fuel = nw.get_conn("0"), nw.get_conn("8")
    assert stack.m.val_SI == approx(air.m.val_SI + fuel.m.val_SI)
    for fluid, value in stack.fluid.val.items():
        assert recirculated.fluid.val[fluid] == approx(value)


def test_inferred_composition_dependence():
    nw = gas_turbine_hrsg()
    nw.solve("design", block_solve_composition=True)
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
    # the water steam side and the pure mass and fluid balances do not
    assert {"turbine", "hrsg"} <= inferred
    assert {"feed pump", "steam turbine"}.isdisjoint(inferred)
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
    # water mass flow from the energy balance of the hrsg follows it
    blocks = problem.decompose().blocks
    composition_blocks = [
        block for block in blocks
        if composition_columns & set(block.variables)
    ]
    assert len(composition_blocks) == 1
    for block in blocks:
        labels = {label for label, _ in block.equation_labels}
        if labels & {"feed pump", "steam turbine", "hrsg"}:
            assert composition_columns.isdisjoint(block.variables)
        if "hrsg" in labels:
            assert block.id > composition_blocks[0].id


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
    nw.solve("design", block_solve_composition=True)
    nw.assert_convergence()
    problem = nw.problem
    assert not problem._variable_composition()
    assert problem._inferred_incidence == {}
    assert problem._declared_incidence is problem.get_incidence()


def test_custom_component_block_solve_composition():
    nw = merge_chain(cooler_class=LegacyCooler)
    nw.solve("design", block_solve_composition=True)
    nw.assert_convergence()
    assert nw.get_conn("5").T.val == approx(80)

    inferred = [
        nw.problem._equation_lookup[eq_num]
        for eq_num in nw.problem._inferred_incidence
    ]
    assert ("cooler", ("T_out", 0)) in inferred
