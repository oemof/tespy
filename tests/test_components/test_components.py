import inspect
import sys

import pytest
from pytest import mark

from tespy.components import Subsystem
from tespy.components.component import component_registry
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import SimpleDataContainer as dc_simple


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
    "SimpleHeatExchanger": {"zeta", "ks_HW"},
    "ParabolicTrough": {"zeta", "c_1", "c_2", "iam_1", "iam_2", "ks_HW"},
    "ParallelFlowHeatExchanger": {"zeta1", "zeta2"},
    "SolarCollector": {"zeta", "lkf_lin", "lkf_quad", "ks_HW"},
    "Pipe": {"zeta", "ks_HW"},
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


def properties_with_eq_of(instance):
    return [
        prop
        for prop, container in instance.get_parameters().items()
        if (
            not isinstance(container, dc_simple)
            and container.func is not None or container.structure_matrix is not None
        )
    ]


def generate_class_property_params(property_fn):
    params = []
    for name, cls in component_registry.items.items():
        instance = cls("")
        for prop in property_fn(instance):
            params.append(pytest.param(name, prop, id=f"{cls.__name__}::{prop}"))
    return params


def pytest_generate_tests(metafunc):
    if {"cls_name", "prop"} <= set(metafunc.fixturenames):
        if metafunc.function.__name__ == "test_property_value_not_none":
            metafunc.parametrize(
                "cls_name,prop",
                generate_class_property_params(properties_of)
            )

        elif metafunc.function.__name__ == "test_num_equations_with_func_or_structure_matrix":
            metafunc.parametrize(
                "cls_name,prop",
                generate_class_property_params(properties_with_eq_of)
            )


def test_property_value_not_none(cls_name, prop):

    instance = component_registry.items[cls_name]("")
    value = instance.get_attr(prop)

    condition = (
        prop in QUANTITY_EXEMPTIONS.get(cls_name, set())
        or value.quantity is not None
    )

    assert condition, f"Quantity for {prop} of {cls_name} must not be None"


def test_num_equations_with_func_or_structure_matrix(cls_name, prop):

    instance = component_registry.items[cls_name]("")
    value = instance.get_attr(prop)

    condition = value.num_eq_sets > 0

    assert condition, f"The parameter {prop} of {cls_name} lacks `num_eq_sets` specification"
