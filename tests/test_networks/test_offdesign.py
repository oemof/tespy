# -*- coding: utf-8

"""Module for testing local_offdesign specification.
This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_local_offdesign.py
SPDX-License-Identifier: MIT
"""
import numpy as np
from pytest import approx
from pytest import raises

from tespy.components import Compressor
from tespy.components import Merge
from tespy.components import Motor
from tespy.components import PowerSink
from tespy.components import PowerSource
from tespy.components import Pump
from tespy.components import SectionedHeatExchanger
from tespy.components import Sink
from tespy.components import SolarCollector
from tespy.components import Source
from tespy.components import Splitter
from tespy.components import Valve
from tespy.connections import Connection
from tespy.connections import PowerConnection
from tespy.networks import Network
from tespy.tools import CharLine
from tespy.tools import load_default_char
from tespy.tools.helpers import TESPyNetworkError


def _build_network_and_designs():
    """
    Build a compressor -> heat-exchanger network and save two design cases.

    Design 1: refrigerant mass flow = 1 kg/s
    Design 2: refrigerant mass flow = 0.8 kg/s

    Returns a dict with the network, all components and connections, the two
    save-file paths, and the parsed JSON dicts for both design cases.
    """
    nw = Network()
    nw.units.set_defaults(
        pressure="bar", temperature="degC", pressure_difference="kPa"
    )

    so1, so2 = Source("source 1"), Source("source 2")
    si1, si2 = Sink("sink 1"), Sink("sink 2")
    compressor = Compressor("compressor")
    heatex = SectionedHeatExchanger("heat exchanger")

    c0 = Connection(so1, "out1", compressor, "in1", label="c0")
    c1 = Connection(compressor, "out1", heatex, "in1", label="c1")
    c2 = Connection(heatex, "out1", si1, "in1", label="c2")
    c3 = Connection(so2, "out1", heatex, "in2", label="c3")
    c4 = Connection(heatex, "out2", si2, "in1", label="c4")

    nw.add_conns(c0, c1, c2, c3, c4)

    compressor.set_attr(eta_s=0.8, design=["eta_s"], offdesign=["eta_s_char"])
    c0.set_attr(x=1, T_dew=20)
    c1.set_attr(fluid={"R290": 1}, T_dew=50, m=1)
    c2.set_attr(x=0)
    c3.set_attr(fluid={"water": 1}, p=2, T=30)
    kA_char = load_default_char(
        "HeatExchanger", "kA_char1", "DEFAULT", CharLine
    )
    heatex.set_attr(
        dp1=2, dp2=10, td_pinch=5, design=["td_pinch"], offdesign=["UA_char"],
        kA_char1=kA_char, kA_char2=kA_char
    )

    # design case 1: m = 1 kg/s
    nw.solve("design")
    design1 = nw.save(as_dict=True)

    # design case 2: m = 0.8 kg/s
    c1.set_attr(m=0.8)
    c1.set_attr(T_dew=52)
    nw.solve("design")
    design2 = nw.save(as_dict=True)

    return dict(
        nw=nw,
        compressor=compressor, heatex=heatex,
        c0=c0, c1=c1, c2=c2, c3=c3, c4=c4,
        design1=design1, design2=design2,
    )


def _build_isolated_design():
    """
    Build a standalone heat-exchanger network (no compressor) and save its
    design. Labels match those used in the bigger network built by
    _build_network_and_designs.
    """
    nw = Network()
    nw.units.set_defaults(
        pressure="bar",
        temperature="degC",
        pressure_difference="kPa"
    )

    so1 = Source("source 1")
    so2 = Source("source 2")
    si1 = Sink("sink 1")
    si2 = Sink("sink 2")

    heatex = SectionedHeatExchanger("heat exchanger")

    c1 = Connection(so1, "out1", heatex, "in1", label="c1")
    c2 = Connection(heatex, "out1", si1, "in1", label="c2")
    c3 = Connection(so2, "out1", heatex, "in2", label="c3")
    c4 = Connection(heatex, "out2", si2, "in1", label="c4")

    nw.add_conns(c1, c2, c3, c4)

    c1.set_attr(fluid={"R290": 1}, T_dew=50, td_dew=25, m=1)
    c2.set_attr(x=0)
    c3.set_attr(fluid={"water": 1}, p=2, T=30)
    heatex.set_attr(
        dp1=2, dp2=10, td_pinch=5, design=["td_pinch"], offdesign=["UA_char"]
    )

    nw.solve("design")
    design = nw.save(as_dict=True)

    return dict(design=design)


def _build_isolated_design_with_different_labels():
    """
    Like _build_isolated_design but with different component and connection
    labels. The heat exchanger type is the same, so the single-type fallback
    should resolve the component, and port-based matching should resolve the
    connections.
    """
    nw = Network()
    nw.units.set_defaults(
        pressure="bar",
        temperature="degC",
        pressure_difference="kPa"
    )

    so1 = Source("source A")
    so2 = Source("source B")
    si1 = Sink("sink A")
    si2 = Sink("sink B")

    heatex = SectionedHeatExchanger("heat exchanger with different label")

    c1 = Connection(so1, "out1", heatex, "in1", label="b1")
    c2 = Connection(heatex, "out1", si1, "in1", label="b2")
    c3 = Connection(so2, "out1", heatex, "in2", label="b3")
    c4 = Connection(heatex, "out2", si2, "in1", label="b4")

    nw.add_conns(c1, c2, c3, c4)

    c1.set_attr(fluid={"R290": 1}, T_dew=50, td_dew=25, m=1.1)
    c2.set_attr(x=0)
    c3.set_attr(fluid={"water": 1}, p=2, T=30)
    heatex.set_attr(
        dp1=2, dp2=10, td_pinch=5, design=["td_pinch"], offdesign=["UA_char"]
    )

    nw.solve("design")
    design = nw.save(as_dict=True)

    return dict(design=design)


def test_individual_design_path_offdesign():
    """
    Component design values: heat exchanger uses its local design path,
    compressor stays on the global design path.
    """
    o = _build_network_and_designs()
    nw, compressor, heatex = o["nw"], o["compressor"], o["heatex"]
    c1 = o["c1"]
    design1, design2 = o["design1"], o["design2"]

    nw.solve("offdesign", design_path=design1, init_only=True)

    assert heatex.UA.design == design1["Component"]["SectionedHeatExchanger"]["heat exchanger"]["UA"]
    assert compressor.eta_s.design == design1["Component"]["Compressor"]["compressor"]["eta_s"]

    assert heatex._conn_design(c1, "m") == c1.m.design
    assert heatex._conn_design(c1, "T") == c1.T.design

    assert compressor._conn_design(c1, "m") == c1.m.design
    assert compressor._conn_design(c1, "T") == c1.T.design

    heatex.set_attr(design_path=design2)
    nw.solve("offdesign", design_path=design1, init_only=True)

    assert heatex.UA.design == design2["Component"]["SectionedHeatExchanger"]["heat exchanger"]["UA"]
    assert compressor.eta_s.design == design1["Component"]["Compressor"]["compressor"]["eta_s"]

    assert approx(heatex._conn_design(c1, "m")) != compressor._conn_design(c1, "m")
    assert approx(heatex._conn_design(c1, "T")) != compressor._conn_design(c1, "T")
    assert heatex._conn_design(c1, "m") != c1.m.design
    assert heatex._conn_design(c1, "T") != c1.T.design


def test_local_offdesign_adjacent_components_different_paths():
    """
    Two adjacent components (compressor -> heat exchanger, sharing connection
    c1) can each have local_offdesign=True pointing to *different* design
    paths. Every component must independently use its own local design
    reference for the shared connection; they must not interfere with each
    other.
    """
    o = _build_network_and_designs()
    nw, compressor, heatex = o["nw"], o["compressor"], o["heatex"]
    design1, design2 = o["design1"], o["design2"]
    c1 = o["c1"]

    compressor.set_attr(design_path=design1, local_offdesign=True)
    heatex.set_attr(design_path=design2, local_offdesign=True)
    nw.solve("design")
    nw.assert_convergence()

    # Each component loads its own design reference
    assert compressor.eta_s.design == design1["Component"]["Compressor"]["compressor"]["eta_s"]
    assert heatex.UA.design == design2["Component"]["SectionedHeatExchanger"]["heat exchanger"]["UA"]

    # The shared connection c1 gets independent local references per component
    local_c1_m_compressor = compressor._conn_design(c1, "m")
    local_c1_m_heatex = heatex._conn_design(c1, "m")

    expected_m_design1 = design1["Connection"]["Connection"]["c1"]["m"]
    expected_m_design2 = design2["Connection"]["Connection"]["c1"]["m"]

    assert local_c1_m_compressor == approx(expected_m_design1)
    assert local_c1_m_heatex == approx(expected_m_design2)

    # The two local references for c1 must differ from each other
    assert local_c1_m_compressor != approx(local_c1_m_heatex)

    # The connection's own .design attribute must be nan because we are in
    # design mode
    assert np.isnan(c1.m.design)


def test_unset_local_offdesign_clears_data():
    """
    After a run with local_offdesign=True a run with local_offdesign=False
    must not be able to reference values from the old design_path
    """
    o = _build_network_and_designs()
    nw, heatex = o["nw"], o["heatex"]
    design2 = o["design2"]
    c1 = o["c1"]

    heatex.set_attr(design_path=design2, local_offdesign=True)
    nw.solve("design", init_only=True)

    expected_m = design2["Connection"]["Connection"]["c1"]["m"]
    assert approx(heatex._conn_design(c1, "m")) == expected_m

    heatex.set_attr(local_offdesign=False)
    nw.solve("design", init_only=True)

    assert np.isnan(heatex._conn_design(c1, "m"))


def test_offdesign_connection_design_stays_global():
    """
    After an offdesign solve the adjacent connections' own .design
    attributes must still reflect the global design path. The local design
    reference lives exclusively in _local_connection_design_state on the
    component; the connection objects themselves are not modified.
    """
    o = _build_network_and_designs()
    nw, heatex = o["nw"], o["heatex"]
    c1 = o["c1"]
    design1, design2 = o["design1"], o["design2"]

    heatex.set_attr(design_path=design2)
    nw.solve("offdesign", design_path=design1)

    # c1 is adjacent to heatex, but its .design comes from the global design1.
    # Mass flow is in kg/s both as SI and as saved unit, so no conversion needed.
    expected_c1_m = design1["Connection"]["Connection"]["c1"]["m"]   # 1 kg/s
    assert c1.m.design == approx(expected_c1_m)
    assert c1.m.design != approx(heatex._conn_design(c1, "m"))


def test_offdesign_non_local_component_uses_global(tmp_path):
    """
    A component without local_offdesign=True must have an empty
    _local_connection_design_state and _conn_design must return the
    connection's own (global) .design value unchanged.
    """
    o = _build_network_and_designs()
    nw, compressor, heatex = o["nw"], o["compressor"], o["heatex"]
    c0, c1 = o["c0"], o["c1"]
    design1, design2 = o["design1"], o["design2"]

    heatex.set_attr(design_path=design2)
    nw.solve("offdesign", design_path=design1)

    # Compressor has no local_offdesign -> empty local state
    assert compressor._local_connection_design_state == {}

    # _conn_design must fall through to the connection's own .design
    assert compressor._conn_design(c0, "m") == c0.m.design
    assert compressor._conn_design(c1, "m") == c1.m.design


def test_local_design_path_reverts_to_global_when_unset():
    """
    After unsetting a component's design_path the global is used again
    """
    o = _build_network_and_designs()
    nw, heatex = o["nw"], o["heatex"]
    design1, design2 = o["design1"], o["design2"]
    c1 = o["c1"]

    # Enable local offdesign -> UA from design2
    heatex.set_attr(design_path=design2)
    nw.solve("offdesign", design_path=design1)
    assert approx(heatex.UA.design) == design2["Component"]["SectionedHeatExchanger"]["heat exchanger"]["UA"]
    assert approx(heatex._conn_design(c1, "m")) != c1.m.design

    # Unset the component-level design_path, _conn_design() must immediately
    # use the global .design because it checks if design_path is available
    # before consulting _local_connection_design_state.
    heatex.set_attr(design_path=None)
    assert approx(heatex._conn_design(c1, "m")) == c1.m.design
    # component specific properties keep their value until next model
    # initialization
    assert approx(heatex.UA.design) == design2["Component"]["SectionedHeatExchanger"]["heat exchanger"]["UA"]

    # UA.design reverts to the global value once _load_offdesign_state() is
    # triggered.  Clearing design_path sets new_design=True on the component,
    # which is enough to trigger a reload on the next solve.
    nw.solve("offdesign", design_path=design1)
    assert approx(heatex.UA.design) == design1["Component"]["SectionedHeatExchanger"]["heat exchanger"]["UA"]


def test_isolated_design_matching_labels():
    """
    Load a local_offdesign component's design state from an isolated design
    file where all labels match. The UA design value and adjacent connection
    design values must come from the isolated design, not the global path.
    """
    o = _build_network_and_designs()
    iso = _build_isolated_design()

    nw, heatex = o["nw"], o["heatex"]
    design1 = o["design1"]

    # The isolated design has the same UA value as design1 (same conditions),
    # so we confirm the component and connections were loaded by checking the
    # _local_connection_design_state is populated.
    heatex.set_attr(design_path=iso["design"])
    nw.solve("offdesign", design_path=design1)
    nw.assert_convergence()

    # Component design value (UA) loaded from the isolated design
    assert heatex.UA.design == approx(
        iso["design"]["Component"]["SectionedHeatExchanger"]["heat exchanger"]["UA"]
    )

    # All adjacent connections must be resolved in the local state
    for label in ("c1", "c2", "c3", "c4"):
        assert label in heatex._local_connection_design_state, (
            f"connection '{label}' missing from _local_connection_design_state"
        )

    # c1 mass flow in the isolated design is 1 kg/s (same as global design1)
    local_c1_m = heatex._local_connection_design_state["c1"]["m"]
    expected = iso["design"]["Connection"]["Connection"]["c1"]["m"]
    assert local_c1_m == approx(expected)


def test_isolated_design_different_labels():
    """
    Load a local_offdesign component's design state from an isolated design
    file where the component label and connection labels differ from the main
    network. The single-type fallback should resolve the component, and
    port-based matching should resolve the connections.
    """
    o = _build_network_and_designs()
    iso = _build_isolated_design_with_different_labels()

    nw, heatex = o["nw"], o["heatex"]
    design1 = o["design1"]

    heatex.set_attr(design_path=iso["design"], local_offdesign=True)
    nw.solve("offdesign", design_path=design1)
    nw.assert_convergence()

    # Component design value loaded via single-type fallback
    isolated_comp_label = "heat exchanger with different label"
    assert heatex.UA.design == approx(
        iso["design"]["Component"]["SectionedHeatExchanger"][isolated_comp_label]["UA"]
    )

    # All adjacent connections resolved via port matching
    for label in ("c1", "c2", "c3", "c4"):
        assert label in heatex._local_connection_design_state, (
            f"connection '{label}' missing from _local_connection_design_state "
            "after port-based matching"
        )

    # c1 connects to heatex via in1; in the isolated design b1 also connects
    # to the (only) heat exchanger via in1 -> must be matched
    local_c1_m = heatex._local_connection_design_state["c1"]["m"]
    # b1 in the isolated design has m=1 kg/s
    assert local_c1_m == approx(
        iso["design"]["Connection"]["Connection"]["b1"]["m"]
    )


def test_UA_char_char_expr_offdesign_design_reference():
    """
    Verify that the UA modification uses different references for their
    mass flows in offdesign with individual design_path and global design path
    """
    o = _build_network_and_designs()
    nw, heatex = o["nw"], o["heatex"]
    c1 = o["c1"]
    design1, design2 = o["design1"], o["design2"]

    c1.set_attr(m=0.9)

    heatex.set_attr(design_path=design2)
    nw.solve("offdesign", design_path=design1)
    nw.assert_convergence()

    f1 = heatex.get_char_expr("m", **heatex.kA_char1.char_params)
    f2 = heatex.get_char_expr("m", **heatex.kA_char2.char_params)

    fUA1 = heatex.kA_char1.char_func.evaluate(f1)
    fUA2 = heatex.kA_char2.char_func.evaluate(f2)

    fUA_a = 2 / (1 / fUA1 + 1 / fUA2)

    assert approx(fUA_a * heatex.UA.design) == heatex.UA.val_SI

    # change to different specification: no design path of component
    heatex.set_attr(design_path=None)
    nw.solve("offdesign", design_path=design1)
    nw.assert_convergence()

    f1 = heatex.get_char_expr("m", **heatex.kA_char1.char_params)
    f2 = heatex.get_char_expr("m", **heatex.kA_char2.char_params)

    fUA1 = heatex.kA_char1.char_func.evaluate(f1)
    fUA2 = heatex.kA_char2.char_func.evaluate(f2)

    fUA_b = 2 / (1 / fUA1 + 1 / fUA2)

    assert approx(fUA_b * heatex.UA.design) == heatex.UA.val_SI
    assert approx(fUA_a) != fUA_b


def test_design_without_topology_info_mismatched_labels_raises():
    """
    An old-style design file (no source/target/source_id/target_id topology
    columns) must work fine when connection labels match, but must raise a
    KeyError with a clear message when labels differ and topology-based
    fallback is therefore unavailable.
    """
    o = _build_network_and_designs()
    iso = _build_isolated_design_with_different_labels()

    nw, heatex = o["nw"], o["heatex"]
    design1 = o["design1"]

    # Strip topology columns to simulate a pre-0.9.15 export
    data = iso["design"]
    for conns in data["Connection"].values():
        for conn_data in conns.values():
            for col in ("source", "target", "source_id", "target_id"):
                conn_data.pop(col, None)

    heatex.set_attr(design_path=data)
    with raises(KeyError):
        nw.solve("offdesign", design_path=design1, init_only=True)


def _build_power_network_and_designs():
    """
    Build a power-only model (PowerSource -> Motor/inverter -> Motor -> PowerSink)
    and save two design cases.

    Design 1: inverter eta=0.95, inlet power E=100 kW
    Design 2: inverter eta=0.90, inlet power E=80 kW

    Returns a dict with the network, all components and connections, and the
    dicts for both design cases.
    """
    nw = Network()
    nw.units.set_defaults(
        pressure="bar", temperature="degC", pressure_difference="kPa",
        power="kW"
    )

    so = PowerSource("source")
    si = PowerSink("sink")
    inverter = Motor("inverter")
    motor = Motor("motor")

    e1 = PowerConnection(so, "power", inverter, "power_in", label="e0")
    e2 = PowerConnection(inverter, "power_out", motor, "power_in", label="e1")
    e3 = PowerConnection(motor, "power_out", si, "power", label="e2")

    nw.add_conns(e1, e2, e3)

    inverter.set_attr(eta=0.95, design=["eta"], offdesign=["eta_char"])
    motor.set_attr(eta=0.95)

    # design case 1: eta=0.95, E=100 kW
    e1.set_attr(E=100)
    nw.solve("design")
    design1 = nw.save(as_dict=True)

    # design case 2: eta=0.90, E=80 kW
    inverter.set_attr(eta=0.90)
    e1.set_attr(E=80)
    nw.solve("design")
    design2 = nw.save(as_dict=True)

    return dict(
        nw=nw,
        inverter=inverter, motor=motor,
        e1=e1, e2=e2, e3=e3,
        design1=design1, design2=design2,
    )


def test_power_individual_design_path_offdesign(tmp_path):
    """
    Check if with a shared power connection, the reference to its design
    value is different from two different components, when they use different
    design paths
    """
    o = _build_power_network_and_designs()
    nw, inverter, motor = o["nw"], o["inverter"], o["motor"]
    e1, e2 = o["e1"], o["e2"]
    design1, design2 = o["design1"], o["design2"]

    # --- Global design1 as baseline ---
    nw.solve("offdesign", design_path=design1, init_only=True)

    assert inverter.eta.design == approx(
        design1["Component"]["Motor"]["inverter"]["eta"]
    )
    assert motor.eta.design == approx(
        design1["Component"]["Motor"]["motor"]["eta"]
    )
    assert approx(inverter._conn_design(e2, "E")) == motor._conn_design(e2, "E")

    # --- Individual design2 for inverter and its inlet connection ---
    inverter.set_attr(design_path=design2)
    nw.solve("offdesign", design_path=design1, init_only=True)

    assert inverter.eta.design == approx(
        design2["Component"]["Motor"]["inverter"]["eta"]
    )
    # motor has no individual path -> still global design1
    assert motor.eta.design == approx(
        design1["Component"]["Motor"]["motor"]["eta"]
    )
    assert approx(inverter._conn_design(e2, "E")) != motor._conn_design(e2, "E")


def test_eta_char_char_expr_design_mode_local_offdesign():
    """
    When the inverter runs with local_offdesign=True and design_path=design2 in
    a design-mode solve, its eta_char_func must evaluate the characteristic
    expression using the design2 power reference. If E changes, the expression
    needs to change.
    """
    o = _build_power_network_and_designs()
    nw, inverter = o["nw"], o["inverter"]
    e1 = o["e1"]
    design2 = o["design2"]

    # design2 conditions still active: inverter.eta=0.90, e1.E=80 kW
    inverter.set_attr(local_offdesign=True, design_path=design2)
    nw.solve("design")
    nw.assert_convergence()

    # eta.design must come from design2
    assert inverter.eta.design == approx(0.90)

    # eta_char expr = E_current / E_design; at the design2 operating point = 1.0
    expr = e1.E.val_SI / inverter._conn_design(e1, "E")
    assert expr == approx(1.0)

    e1.set_attr(E=100)

    nw.solve("design")
    nw.assert_convergence()
    # eta_char expr = E_current / E_design = 1.25
    expr = e1.E.val_SI / inverter._conn_design(e1, "E")
    assert expr == approx(100 / 80)


def test_eta_char_char_expr_offdesign_design_reference():
    """
    Verify that eta_char_func uses the correct design power reference in
    offdesign mode depending on whether the inlet connection has an individual
    design_path.
    """
    o = _build_power_network_and_designs()
    nw, inverter = o["nw"], o["inverter"]
    e1 = o["e1"]
    design1, design2 = o["design1"], o["design2"]

    e1.set_attr(E=90)

    # individual design_path=design2
    inverter.set_attr(design_path=design2)
    nw.solve("offdesign", design_path=design1)
    nw.assert_convergence()

    assert approx(inverter._conn_design(e1, "E")) != e1.E.design
    assert inverter.eta.design == approx(0.90)
    expr_a = e1.E.val_SI / inverter._conn_design(e1, "E")

    # no individual design_path (global design1)
    inverter.set_attr(design_path=None)
    nw.solve("offdesign", design_path=design1)
    nw.assert_convergence()

    assert approx(inverter._conn_design(e1, "E")) == e1.E.design
    assert inverter.eta.design == approx(0.95)
    expr_b = e1.E.val_SI / inverter._conn_design(e1, "E")

    # The two design references and characteristic expressions must differ
    assert expr_a != approx(expr_b)


class TestNetworkIndividualOffdesign:

    def setup_Network_individual_offdesign(self):
        """Set up network for individual offdesign tests."""
        self.nw = Network()
        self.nw.units.set_defaults(**{
            "pressure": "bar", "pressure_difference": "bar",
            "temperature": "degC", "volumetric_flow": "m3/s"
        })

        so = Source('source')
        sp = Splitter('splitter', num_out=2)
        self.pump1 = Pump('pump 1')
        self.sc1 = SolarCollector('collector field 1')
        v1 = Valve('valve1')
        self.pump2 = Pump('pump 2')
        self.sc2 = SolarCollector('collector field 2')
        v2 = Valve('valve2')
        me = Merge('merge', num_in=2)
        si = Sink('sink')

        self.pump1.set_attr(
            eta_s=0.8, design=['eta_s'], offdesign=['eta_s_char']
        )
        self.pump2.set_attr(
            eta_s=0.8, design=['eta_s'], offdesign=['eta_s_char']
        )
        self.sc1.set_attr(
            pr=0.95, lkf_lin=3.33, lkf_quad=0.011, A=1252, E=700,
            Tamb=20, eta_opt=0.92, design=['pr'], offdesign=['zeta_d4']
        )
        self.sc2.set_attr(
            pr=0.95, lkf_lin=3.5, lkf_quad=0.011, A=700, E=800,
            Tamb=20, eta_opt=0.92, design=['pr'], offdesign=['zeta_d4']
        )

        fl = {'H2O': 1}
        inlet = Connection(so, 'out1', sp, 'in1', T=50, p=3, fluid=fl)
        outlet = Connection(me, 'out1', si, 'in1', p=3)

        self.sp_p1 = Connection(sp, 'out1', self.pump1, 'in1')
        self.p1_sc1 = Connection(self.pump1, 'out1', self.sc1, 'in1')
        self.sc1_v1 = Connection(self.sc1, 'out1', v1, 'in1', p=3.1, T=90)
        v1_me = Connection(v1, 'out1', me, 'in1')

        self.sp_p2 = Connection(sp, 'out2', self.pump2, 'in1')
        self.p2_sc2 = Connection(self.pump2, 'out1', self.sc2, 'in1')
        self.sc2_v2 = Connection(self.sc2, 'out1', v2, 'in1', p=3.1, m=0.1)
        v2_me = Connection(v2, 'out1', me, 'in2')

        self.nw.add_conns(
            inlet, outlet, self.sp_p1, self.p1_sc1, self.sc1_v1,
            v1_me, self.sp_p2, self.p2_sc2, self.sc2_v2, v2_me
        )

    def test_individual_design_path_on_connections_and_components(self):
        """Test individual design path specification."""
        self.setup_Network_individual_offdesign()
        self.nw.solve('design')
        self.nw.assert_convergence()
        self.sc2_v2.set_attr(m=0)
        self.nw.solve('design')
        self.nw.assert_convergence()
        design1 = self.nw.save(as_dict=True)
        v1_design = self.sc1_v1.v.val_SI
        zeta_sc1_design = self.sc1.zeta.val

        self.sc2_v2.set_attr(T=95, state='l', m=None)
        self.sc1_v1.set_attr(m=0.001, T=None)
        self.nw.solve('design')
        self.nw.assert_convergence()
        design2 = self.nw.save(as_dict=True)
        v2_design = self.sc2_v2.v.val_SI
        zeta_sc2_design = self.sc2.zeta.val

        self.sc1_v1.set_attr(m=None)
        self.sc1_v1.set_attr(design=['T'], offdesign=['v'], state='l')
        self.sc2_v2.set_attr(design=['T'], offdesign=['v'], state='l')

        self.sc2.set_attr(design_path=design2)
        self.pump2.set_attr(design_path=design2)
        self.sp_p2.set_attr(design_path=design2)
        self.p2_sc2.set_attr(design_path=design2)
        self.sc2_v2.set_attr(design_path=design2)
        self.nw.solve('offdesign', design_path=design1)
        self.nw.assert_convergence()

        self.sc1.set_attr(E=500)
        self.sc2.set_attr(E=950)

        self.nw.solve('offdesign', design_path=design1)
        self.nw.assert_convergence()
        self.sc2_v2.set_attr(design_path=None)

        # volumetric flow comparison
        msg = f"Design path was set to None, is {self.sc2_v2.design_path}."
        assert self.sc2_v2.design_path is None, msg

        # volumetric flow comparison
        msg = (
            f"Value of volumetric flow must be {v1_design}, is "
            f"{self.sc1_v1.v.val_SI}."
        )
        assert round(v1_design, 5) == round(self.sc1_v1.v.val_SI, 5), msg

        msg = (
            f"Value of volumetric flow must be {v2_design}, is "
            f"{self.sc2_v2.v.val_SI}."
        )
        assert round(v2_design, 5) == round(self.sc2_v2.v.val_SI, 5), msg

        # zeta value of solar collector comparison
        msg = (
            f"Value of zeta must be {zeta_sc1_design}, is {self.sc1.zeta.val}."
        )
        assert round(zeta_sc1_design, 0) == round(self.sc1.zeta.val, 0), msg

        msg = (
            f"Value of zeta must be {zeta_sc2_design}, is {self.sc2.zeta.val}."
        )
        assert round(zeta_sc2_design, 0) == round(self.sc2.zeta.val, 0), msg

    def test_local_offdesign_on_connections_and_components(self):
        """Test local offdesign feature."""
        self.setup_Network_individual_offdesign()
        self.nw.solve('design')
        self.nw.assert_convergence()
        self.sc2_v2.set_attr(m=0)
        self.nw.solve('design')
        self.nw.assert_convergence()
        design1 = self.nw.save(as_dict=True)

        self.sc1_v1.set_attr(design=['T'], offdesign=['v'], state='l')
        self.sc2_v2.set_attr(design=['T'], offdesign=['v'], state='l')

        self.sc1.set_attr(local_offdesign=True, design_path=design1)
        self.pump1.set_attr(local_offdesign=True, design_path=design1)
        self.sp_p1.set_attr(local_offdesign=True, design_path=design1)
        self.p1_sc1.set_attr(local_offdesign=True, design_path=design1)
        self.sc1_v1.set_attr(local_offdesign=True, design_path=design1)
        self.sc1.set_attr(E=500)

        self.sc2_v2.set_attr(T=95, m=None)
        self.nw.solve('design')
        self.nw.assert_convergence()

        # connections and components on side 1 must have switched to offdesign

        msg = (
            'Solar collector outlet temperature must be different from design '
            f'value {round(self.sc1_v1.T.design - 273.15, 1)}, is '
            f'{round(self.sc1_v1.T.val, 1)}.'
        )
        assert self.sc1_v1.T.design > self.sc1_v1.T.val, msg

        msg = "Parameter eta_s_char must be set for pump one."
        assert self.pump1.eta_s_char.is_set, msg

        msg = (
            "Parameter v must be set for connection from solar collector1 to "
            "pump1."
        )
        assert self.sc1_v1.v.is_set, msg

    def test_missing_design_path_local_offdesign_on_connections(self):
        """Test missing design path on connections in local offdesign mode."""
        self.setup_Network_individual_offdesign()
        self.nw.solve('design')
        self.nw.assert_convergence()
        self.sc2_v2.set_attr(m=0)
        self.nw.solve('design')
        self.nw.assert_convergence()
        design_state = self.nw.save(as_dict=True)

        self.sc1_v1.set_attr(design=['T'], offdesign=['v'], state='l')
        self.sc2_v2.set_attr(design=['T'], offdesign=['v'], state='l')

        self.sc1.set_attr(local_offdesign=True, design_path=design_state)
        self.pump1.set_attr(local_offdesign=True, design_path=design_state)
        self.sp_p1.set_attr(local_offdesign=True, design_path=design_state)
        self.p1_sc1.set_attr(local_offdesign=True, design_path=design_state)
        self.sc1_v1.set_attr(local_offdesign=True)
        self.sc1.set_attr(E=500)

        self.sc2_v2.set_attr(T=95, m=None)
        try:
            self.nw.solve('design', init_only=True)
        except TESPyNetworkError:
            pass
