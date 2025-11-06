import inspect
import sys

import pytest
from pytest import mark

from tespy.components import Subsystem
from tespy.components.component import component_registry
from tespy.tools.data_containers import ComponentProperties as dc_cp

ALL_COMPONENT_CLASSES = [
    obj for _, obj in inspect.getmembers(sys.modules["tespy.components"])
    # exclude the Subsystem component as it is just a wrapper
    if inspect.isclass(obj) and obj is not Subsystem
]

@mark.parametrize("obj", ALL_COMPONENT_CLASSES)
def test_all_classes_in_registry(obj):
    msg = (
        f"The class {obj.__name__} is missing in the "
        "component_registry"
    )
    assert obj in component_registry.items.values(), msg


QUANTITY_EXEMPTIONS = {
    "CycleCloser": {"fluid_deviation"},
    "CombustionEngine": {"zeta1", "zeta2"},
    "PolynomialCompressor": {"rpm", "re_exp_r", "re_exp_sf"},
    "PolynomialCompressorWithCooling": {"rpm", "re_exp_r", "re_exp_sf"},
    "HeatExchanger": {"zeta1", "zeta2"},
    "Condenser": {"zeta1", "zeta2"},
    "Desuperheater": {"zeta1", "zeta2"},
    "SectionedHeatExchanger": {"zeta1", "zeta2", "re_exp_r", "re_exp_sf"},
    "MovingBoundaryHeatExchanger": {"zeta1", "zeta2", "re_exp_r", "re_exp_sf"},
    "SimpleHeatExchanger": {"zeta"},
    "ParabolicTrough": {"zeta", "c_1", "c_2", "iam_1", "iam_2"},
    "ParallelFlowHeatExchanger": {"zeta1", "zeta2"},
    "SolarCollector": {"zeta", "lkf_lin", "lkf_quad"},
    "Pipe": {"zeta"},
    "Valve": {"zeta"},
    "FuelCell": {"zeta"},
    "WaterElectrolyzer": {"zeta"}
}

def properties_of(instance):
    return [
        prop
        for prop, container in instance.get_parameters().items()
        if isinstance(container, dc_cp)
    ]

def pytest_generate_tests(metafunc):
    if "cls_name" in metafunc.fixturenames and "prop" in metafunc.fixturenames:
        params = []
        for name, cls in component_registry.items.items():
            instance = cls("")
            for prop in properties_of(instance):
                params.append(pytest.param(name, prop, id=f"{cls.__name__}::{prop}"))
        metafunc.parametrize("cls_name,prop", params)

def test_property_value_not_none(cls_name, prop):

    instance = component_registry.items[cls_name]("")
    value = instance.get_attr(prop)

    condition = (
        prop in QUANTITY_EXEMPTIONS.get(cls_name, set())
        or value.quantity is not None
    )

    assert condition, f"Quantity for {prop} of {cls_name} must not be None"
