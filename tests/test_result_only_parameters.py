"""Reject output quantities used as model specifications."""

import math

import pint
import pytest

from tespy.components import Condenser
from tespy.components import CycleCloser
from tespy.components import DiabaticCombustionChamber
from tespy.components import HeatExchanger
from tespy.components import NTUHeatExchanger
from tespy.components import Pipe
from tespy.components import PolynomialCompressor
from tespy.components import PolynomialCompressorWithCooling
from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.connections import HAConnection
from tespy.connections import Ref
from tespy.networks import Network

RESULT_PARAMETERS = [
    (HeatExchanger, "lmtd", "delta_degC"),
    (Condenser, "lmtd", "delta_degC"),
    (SimpleHeatExchanger, "lmtd", "delta_degC"),
    (Pipe, "lmtd", "delta_degC"),
    (NTUHeatExchanger, "C_r", ""),
    (CycleCloser, "mass_deviation", "kg/s"),
    (CycleCloser, "fluid_deviation", ""),
    (DiabaticCombustionChamber, "Qloss", "W"),
    (PolynomialCompressor, "Q_diss", "W"),
    (PolynomialCompressorWithCooling, "Q_diss", "W"),
    (PolynomialCompressorWithCooling, "td_minimal", "delta_degC"),
    (Connection, "s", "J/kg/K"),
    (Connection, "vol", "m^3/kg"),
    (HAConnection, "s", "J/kg/K"),
    (HAConnection, "vol", "m^3/kg"),
    (HAConnection, "mHA", "kg/s"),
    (HAConnection, "mH2O", "kg/s"),
]


def make_instance(cls, **kwargs):
    if issubclass(cls, Connection):
        return cls(Source("source"), "out1", Sink("sink"), "in1", **kwargs)
    return cls("heat exchanger", **kwargs)


@pytest.mark.parametrize("cls,parameter,unit", RESULT_PARAMETERS)
@pytest.mark.parametrize("value_kind", ["number", "dict", "quantity", "metadata", "var"])
@pytest.mark.parametrize("on_creation", [False, True])
def test_result_only_specification_rejected(cls, parameter, unit, value_kind, on_creation):
    if value_kind == "number":
        value = 10
    elif value_kind == "var":
        value = "var"
    elif value_kind == "quantity":
        value = pint.UnitRegistry().Quantity(10, unit)
    else:
        value = {"val": 10, "unit": unit, "is_set": True}
        if value_kind == "metadata":
            value["result_only"] = False

    if on_creation:
        with pytest.raises(TypeError, match="result-only"):
            make_instance(cls, **{parameter: value})
    else:
        instance = make_instance(cls)
        container = instance.get_attr(parameter)
        container.set_attr(val=2, unit=unit)
        before = container._serialize()
        with pytest.raises(TypeError, match="result-only"):
            instance.set_attr(**{parameter: value})
        assert container._serialize() == before
        assert container.result_only


@pytest.mark.parametrize("cls,parameter,unit", RESULT_PARAMETERS)
def test_result_only_scalar_can_be_unset(cls, parameter, unit):
    instance = make_instance(cls)
    instance.set_attr(**{parameter: None})
    assert not instance.get_attr(parameter).is_set


@pytest.mark.parametrize("parameter", ["s", "vol"])
def test_result_only_reference_still_rejected(parameter):
    conn = make_instance(Connection)
    with pytest.raises(NotImplementedError, match=f"Referencing {parameter}"):
        conn.set_attr(**{parameter: Ref(conn, 1, 0)})
    assert not conn.get_attr(parameter).is_set


@pytest.mark.parametrize("cls,parameter,unit", [
    item for item in RESULT_PARAMETERS if issubclass(item[0], Connection)
])
def test_result_only_starting_value_rejected(cls, parameter, unit):
    conn = make_instance(cls)
    with pytest.raises(TypeError, match="result-only"):
        conn.set_attr(**{f"{parameter}0": 10})
    assert math.isnan(conn.get_attr(parameter).val0)
    conn.set_attr(**{f"{parameter}0": None})


@pytest.mark.parametrize("cls", [PolynomialCompressor, PolynomialCompressorWithCooling])
@pytest.mark.parametrize("parameter,value", [
    ("P", 1000), ("eta_s", 0.7), ("eta_vol", 0.8),
    ("frequency", 50), ("frequency", "var"), ("rpm", 1500), ("rpm", "var"),
])
def test_calculated_group_input_still_accepted(cls, parameter, value):
    instance = make_instance(cls, **{parameter: value})
    container = instance.get_attr(parameter)
    assert container.is_set
    if value == "var":
        assert container.is_var
    else:
        assert container.val == value


@pytest.mark.parametrize("cls,old,new", [
    (HeatExchanger, "kA", "UA"),
    (HeatExchanger, "zeta1", "zeta1_d4"),
    (SimpleHeatExchanger, "zeta", "zeta_d4"),
])
def test_deprecated_input_alias_still_accepted(cls, old, new):
    instance = make_instance(cls)
    with pytest.warns(FutureWarning):
        instance.set_attr(**{old: 10})
    assert instance.get_attr(new).is_set
    assert instance.get_attr(new).val == 10


def test_result_only_outputs_survive_solve_and_network_roundtrip():
    hx = SimpleHeatExchanger("heat exchanger", Q=-10000, pr=1, Tamb=280)
    inlet = Connection(Source("source"), "out1", hx, "in1", label="inlet")
    outlet = Connection(hx, "out1", Sink("sink"), "in1", label="outlet")
    inlet.set_attr(m=1, p=1e5, T=300, fluid={"water": 1})
    network = Network(iterinfo=False)
    network.add_conns(inlet, outlet)
    network.solve("design")
    network.assert_convergence()

    expected = [hx.lmtd.val_SI, outlet.s.val_SI, outlet.vol.val_SI]
    assert all(math.isfinite(value) and value > 0 for value in expected)
    restored = Network.from_dict(network.export())
    restored.solve("design")
    restored.assert_convergence()
    hx = restored.get_comp("heat exchanger")
    outlet = restored.get_conn("outlet")
    for container, value in zip([hx.lmtd, outlet.s, outlet.vol], expected):
        assert not container.is_set
        assert container.val_SI == pytest.approx(value)
