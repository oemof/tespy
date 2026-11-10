# -*- coding: utf-8

"""Tests spanning all connection classes of the registry.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_connections/test_registry.py

SPDX-License-Identifier: MIT
"""
import inspect
import sys

import pytest
from pytest import mark

from tespy.components import PowerSink
from tespy.components import PowerSource
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.connections import HeatConnection
from tespy.connections import PowerConnection
from tespy.connections import Ref
from tespy.connections.connection import ConnectionBase
from tespy.connections.connection import connection_registry
from tespy.connections.humidairconnection import HAConnection
from tespy.tools.data_containers import FluidProperties as dc_prop

ALL_CONNECTION_CLASSES = [
    obj for _, obj in inspect.getmembers(sys.modules["tespy.connections"])
    # exclude the Subsystem component as it is just a wrapper
    if inspect.isclass(obj)
    and obj not in {ConnectionBase, Ref}
]

@mark.parametrize("obj", ALL_CONNECTION_CLASSES)
def test_all_classes_in_registry(obj):
    msg = (
        f"The class {obj.__name__} is missing in the "
        "connection_registry"
    )
    assert obj in connection_registry.items.values(), msg


def make_connection(cls):
    if cls == Connection or cls == HAConnection:
        return cls(Source(""), "out1", Sink(""), "in1")
    elif cls == PowerConnection:
        return cls(PowerSource(""), "power", PowerSink(""), "power")
    elif cls == HeatConnection:
        from tespy.components import HeatSink
        from tespy.components import HeatSource
        return cls(HeatSource(""), "heat", HeatSink(""), "heat")
    else:
        raise NotImplementedError(
            f"The connection class {cls} is not implemented in testing"
        )


QUANTITY_EXEMPTIONS = {}


def properties_of(instance):
    return [
        prop
        for prop, container in instance.get_parameters().items()
        if isinstance(container, dc_prop)
    ]


def pytest_generate_tests(metafunc):
    if "cls_name" in metafunc.fixturenames and "prop" in metafunc.fixturenames:
        params = []
        for name, cls in connection_registry.items.items():
            instance = make_connection(cls)
            for prop in properties_of(instance):
                params.append(pytest.param(name, prop, id=f"{cls.__name__}::{prop}"))
        metafunc.parametrize("cls_name,prop", params)


def test_property_value_not_none(cls_name, prop):

    instance = make_connection(connection_registry.items[cls_name])
    value = instance.get_attr(prop)

    condition = (
        prop in QUANTITY_EXEMPTIONS.get(cls_name, set())
        or value.quantity is not None
    )

    assert condition, f"Quantity for {prop} of {cls_name} must not be None"
