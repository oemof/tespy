# -*- coding: utf-8

"""Tests of the files written by tespy 0.11.4 for loading in tespy 0.12.

Version 0.12 removes the deprecated parameter names. 0.11.4 writes exports,
states and specification files with the current names only, so that a model
saved with 0.11.4 loads in 0.12.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_forward_compatible_files.py

SPDX-License-Identifier: MIT
"""
import warnings
from copy import deepcopy

import pytest

from tespy.components import HeatExchanger
from tespy.components import PolynomialCompressor
from tespy.components import Pump
from tespy.components import SectionedHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Valve
from tespy.components.displacementmachinery.polynomial_compressor import swept_volume_from_displacement
from tespy.connections import Connection
from tespy.networks import Network
from tespy.networks.network import FORMAT_VERSION
from tespy.tools.characteristics import CharLine
from tespy.tools.characteristics import CharMap

DEPRECATED_NAMES = {
    "kA", "kA_char", "kA_char1", "kA_char2", "kA_group", "kA_char_group",
    "zeta", "zeta1", "zeta2", "td_log", "rpm", "eta_vol_group_rpm",
    "eta_vol_poly_group_rpm", "re_exp_r", "re_exp_sf", "refrigerant_index",
    "UA_cecchinato",
}


def _two_sided(cls):
    nw = Network(iterinfo=False)
    nw.units.set_defaults(
        pressure="bar", pressure_difference="bar", temperature="degC"
    )
    hx = cls("hx")
    c1 = Connection(Source("hot in"), "out1", hx, "in1", label="1")
    c2 = Connection(hx, "out1", Sink("hot out"), "in1", label="2")
    c3 = Connection(Source("cold in"), "out1", hx, "in2", label="3")
    c4 = Connection(hx, "out2", Sink("cold out"), "in1", label="4")
    nw.add_conns(c1, c2, c3, c4)
    c1.set_attr(fluid={"water": 1}, p=5, T=120, m=1)
    c2.set_attr(T=60)
    c3.set_attr(fluid={"water": 1}, p=2, T=20, m=2)
    return nw, hx


def test_format_version_is_2():
    assert FORMAT_VERSION == 2
    nw, hx = _two_sided(HeatExchanger)
    hx.set_attr(pr1=1, pr2=1)
    nw.solve("design")
    for data in (nw.export(), nw.save(as_dict=True), nw.save_specifications()):
        assert data["meta"]["format_version"] == 2


def test_deprecated_names_are_not_written():
    nw, hx = _two_sided(HeatExchanger)
    hx.set_attr(pr1=1, pr2=1, design=["pr1"], offdesign=["zeta1_d4"])
    nw.solve("design")
    exported = nw.export()["Component"]["HeatExchanger"]["hx"]
    assert not DEPRECATED_NAMES & set(exported)
    specs = nw.save_specifications()["Component"]["hx"]
    assert not DEPRECATED_NAMES & set(specs)


def test_specification_of_deprecated_name_moves_to_successor():
    """Files of earlier versions may hold the specification under the
    deprecated name, e.g. kA, and the name in the design lists."""
    nw, hx = _two_sided(HeatExchanger)
    hx.set_attr(pr1=1, pr2=1)
    nw.solve("design")
    data = nw.export()
    legacy = data["Component"]["HeatExchanger"]["hx"]
    legacy["kA"] = {**legacy.pop("UA"), "val": 1200.0, "is_set": True}
    legacy["zeta1"] = {"val": 100.0, "is_set": True}
    legacy["offdesign"] = ["zeta1", "kA_char"]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        loaded = Network.from_dict(deepcopy(data))

    exported = loaded.export()["Component"]["HeatExchanger"]["hx"]
    assert not DEPRECATED_NAMES & set(exported)
    assert exported["UA"]["is_set"]
    assert exported["UA"]["val"] == pytest.approx(1200)
    assert exported["zeta1_d4"]["is_set"]
    assert exported["offdesign"] == ["zeta1_d4", "UA_char"]


def test_rpm_is_written_as_frequency():
    nw = Network(iterinfo=False)
    nw.units.set_defaults(
        pressure="bar", pressure_difference="bar", temperature="degC"
    )
    compressor = PolynomialCompressor("compressor")
    c1 = Connection(Source("source"), "out1", compressor, "in1", label="1")
    c2 = Connection(compressor, "out1", Sink("sink"), "in1", label="2")
    nw.add_conns(c1, c2)
    c1.set_attr(fluid={"R134a": 1}, x=1, T=5)
    c2.set_attr(p=6, T=150)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        compressor.set_attr(
            eta_vol=0.8, rpm=1500, dissipation_ratio=0,
            reference_state={
                "T_sh": 10, "T_sc": 0, "frequency_poly": 50,
                "displacement": 214, "frequency_displacement": 20
            },
            design=["rpm"]
        )
        nw.solve("design")
    nw.assert_convergence()

    exported = nw.export()["Component"]["PolynomialCompressor"]["compressor"]
    assert not DEPRECATED_NAMES & set(exported)
    assert exported["frequency"]["is_set"]
    assert exported["frequency"]["val"] == pytest.approx(25)
    assert exported["frequency"]["unit"] == "1/s"
    assert exported["design"] == ["frequency"]

    state = nw.save(as_dict=True)["Component"]["PolynomialCompressor"]
    assert state["compressor"]["frequency"] == pytest.approx(25)
    assert state["compressor"]["frequency_unit"] == "1/s"

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        reloaded = Network.from_dict(nw.export())
        reloaded.solve("design")
    reloaded.assert_convergence()
    swept_volume = swept_volume_from_displacement(214, 20)
    assert reloaded.get_conn("1").v.val_SI == pytest.approx(
        swept_volume * 25 * 0.8, rel=1e-3
    )


@pytest.mark.parametrize("refrigerant_index", [0, 1])
def test_legacy_cecchinato_is_written_in_hot_cold_convention(refrigerant_index):
    nw, hx = _two_sided(SectionedHeatExchanger)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        hx.set_attr(
            dp1=0, dp2=0, area_ratio=4, alpha_ratio=2, re_exp_r=0.8,
            re_exp_sf=0.5, refrigerant_index=refrigerant_index,
            offdesign=["UA_cecchinato"]
        )
    exported = nw.export()["Component"]["SectionedHeatExchanger"]["hx"]
    assert not DEPRECATED_NAMES & set(exported)
    assert exported["offdesign"] == ["UA_cecchinato_hc"]
    if refrigerant_index == 0:
        expected = {"re_exp_hot": 0.8, "re_exp_cold": 0.5,
                    "area_ratio": 4, "alpha_ratio": 2}
    else:
        expected = {"re_exp_hot": 0.5, "re_exp_cold": 0.8,
                    "area_ratio": 0.25, "alpha_ratio": 0.5}
    for key, value in expected.items():
        assert exported[key]["is_set"]
        assert exported[key]["val"] == pytest.approx(value)


def test_dimensioned_characteristics_declare_si_units():
    nw = Network(iterinfo=False)
    nw.units.set_defaults(
        pressure="bar", pressure_difference="bar", temperature="degC"
    )
    pump = Pump("pump")
    valve = Valve("valve")
    c1 = Connection(Source("source"), "out1", pump, "in1", label="1")
    c2 = Connection(pump, "out1", valve, "in1", label="2")
    c3 = Connection(valve, "out1", Sink("sink"), "in1", label="3")
    nw.add_conns(c1, c2, c3)
    pump.set_attr(
        flow_char={
            "char_func": CharLine([0, 1e-3, 2e-3], [3e5, 2e5, 1e5]),
            "is_set": False
        },
        head_flow_map={"char_func": CharMap(), "is_set": False}
    )
    valve.set_attr(dp_char={
        "char_func": CharLine([0, 1, 2], [0, 1e4, 2e4]), "is_set": False
    })
    exported = nw.export()["Component"]
    flow_char = exported["Pump"]["pump"]["flow_char"]["char_func"]
    assert (flow_char["x_unit"], flow_char["y_unit"]) == ("m3/s", "Pa")
    head_map = exported["Pump"]["pump"]["head_flow_map"]["char_func"]
    assert head_map["z_unit"] == "m"
    dp_char = exported["Valve"]["valve"]["dp_char"]["char_func"]
    assert (dp_char["x_unit"], dp_char["y_unit"]) == ("kg/s", "Pa")

    reloaded = Network.from_dict(nw.export())
    assert reloaded.get_comp("pump").flow_char.char_func.evaluate(1e-3) == (
        pytest.approx(2e5)
    )
