from tespy.components.component import component_registry
from tespy.components import Subsystem
import sys, inspect
from pytest import mark
import pytest


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
