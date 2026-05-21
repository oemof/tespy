# -*- coding: utf-8

"""Utility for generating a JSON schema description of TESPy components
and connections, listing each parameter with its data-container type,
physical quantity and description.
"""

import json
import warnings

from tespy.tools.data_containers import ComponentCharacteristicMaps
from tespy.tools.data_containers import ComponentCharacteristics
from tespy.tools.data_containers import ComponentMandatoryConstraints
from tespy.tools.data_containers import ComponentProperties
from tespy.tools.data_containers import FluidComposition
from tespy.tools.data_containers import FluidProperties
from tespy.tools.data_containers import GroupedComponentCharacteristics
from tespy.tools.data_containers import GroupedComponentProperties
from tespy.tools.data_containers import ReferencedFluidProperties
from tespy.tools.data_containers import SimpleDataContainer

# Human-readable names for each data-container class
_DC_TYPE_NAMES = {
    ComponentProperties: "ComponentProperty",
    ComponentCharacteristics: "ComponentCharacteristic",
    ComponentCharacteristicMaps: "ComponentCharacteristicMap",
    GroupedComponentProperties: "GroupedComponentProperties",
    GroupedComponentCharacteristics: "GroupedComponentCharacteristics",
    ComponentMandatoryConstraints: "ComponentMandatoryConstraint",
    FluidProperties: "FluidProperties",
    FluidComposition: "FluidComposition",
    ReferencedFluidProperties: "ReferencedFluidProperties",
    SimpleDataContainer: "SimpleDataContainer",
}


def _dc_type_name(dc_instance):
    """Return the human-readable data-container type name."""
    for cls, name in _DC_TYPE_NAMES.items():
        if type(dc_instance) is cls:
            return name
    # Fallback: use the class name directly
    return type(dc_instance).__name__


def _serialize_parameter(name, dc_instance):
    """Build a dict describing one parameter from its data-container."""
    entry = {
        "name": name,
        "container_type": _dc_type_name(dc_instance),
        "description": getattr(dc_instance, "description", None),
    }

    # Physical quantity (present on ComponentProperties / FluidProperties
    # and subclasses, absent on SimpleDataContainer etc.)
    quantity = getattr(dc_instance, "quantity", None)
    if quantity is not None:
        entry["quantity"] = quantity

    # For grouped containers, list the member parameter names
    if isinstance(dc_instance, (GroupedComponentProperties, GroupedComponentCharacteristics)):
        entry["elements"] = list(getattr(dc_instance, "elements", []))

    return entry


def _instantiate_component(cls):
    """
    Return a minimally initialised instance of a component class.

    Node-like classes declare get_parameters() as a @staticmethod and do not
    need an instance.  For all other classes we run the base Component.__init__
    with a harmless dummy label, which is the safest path because some
    get_parameters() overrides inspect self.inl / self.outl counts.
    """
    import inspect
    if isinstance(inspect.getattr_static(cls, "get_parameters"), staticmethod):
        return None  # no instance needed

    from tespy.components.component import Component
    instance = cls.__new__(cls)
    Component.__init__(instance, "_schema_probe_")
    return instance


def _get_parameters_for_class(cls, instance=None):
    """Return the get_parameters() dict for *cls*, reusing *instance* if given."""
    import inspect
    if isinstance(inspect.getattr_static(cls, "get_parameters"), staticmethod):
        return cls.get_parameters()
    if instance is None:
        instance = _instantiate_component(cls)
    return instance.get_parameters()



def generate_component_schema(
    exclude_base_classes=("Component",),
    as_json=True,
    indent=2,
):
    """
    Return a serialized description of all registered TESPy component classes.

    For every class in the component registry (except those listed in
    *exclude_base_classes*) the function collects every parameter returned by
    ``get_parameters()`` and records:

    * the parameter name
    * the data-container type (e.g. ``ComponentProperty``, ``ComponentCharacteristic``)
    * the physical quantity (e.g. ``"pressure"``, ``"efficiency"``) when available
    * the human-readable description

    Parameters
    ----------
    exclude_base_classes : tuple[str], optional
        Names of registry entries to skip.  Defaults to ``("Component",)``.
    as_json : bool, optional
        When *True* (default) return a JSON string; when *False* return the
        raw nested dict.
    indent : int, optional
        JSON indentation level (only used when *as_json* is True).

    Returns
    -------
    str | dict
        JSON string or dict depending on *as_json*.
    """
    warnings.warn(
        "tespy.tools.schema is not yet stable and may change without notice in future releases.",
        FutureWarning,
        stacklevel=2,
    )
    from tespy.components.component import component_registry

    schema = {}

    for class_name, cls in component_registry.items.items():
        if class_name in exclude_base_classes:
            continue

        try:
            instance = _instantiate_component(cls)
            parameters = _get_parameters_for_class(cls, instance)
        except Exception as exc:
            schema[class_name] = {"error": str(exc)}
            continue

        params_list = [
            _serialize_parameter(param_name, dc)
            for param_name, dc in parameters.items()
        ]

        # e.g. "tespy.components.turbomachinery.turbine" → "turbomachinery"
        module_parts = cls.__module__.split(".")
        submodule = module_parts[2] if len(module_parts) > 2 else cls.__module__

        ports = cls.port_schema()

        schema[class_name] = {
            "module": cls.__module__,
            "submodule": submodule,
            "inlets": ports["inlets"],
            "outlets": ports["outlets"],
            "powerinlets": ports["powerinlets"],
            "poweroutlets": ports["poweroutlets"],
            "heatinlets": ports["heatinlets"],
            "heatoutlets": ports["heatoutlets"],
            "parameters": params_list,
        }

    if as_json:
        return json.dumps(schema, indent=indent)
    return schema


def generate_connection_schema(
    exclude_base_classes=(),
    as_json=True,
    indent=2,
):
    """
    Return a serialized description of all registered TESPy connection classes.

    Parameters
    ----------
    exclude_base_classes : tuple[str], optional
        Names of registry entries to skip.  Defaults to
        ``("Connection", "PowerConnection")``.
    as_json : bool, optional
        When *True* (default) return a JSON string; when *False* return the
        raw nested dict.
    indent : int, optional
        JSON indentation level (only used when *as_json* is True).

    Returns
    -------
    str | dict
        JSON string or dict depending on *as_json*.
    """
    warnings.warn(
        "tespy.tools.schema is not yet stable and may change without notice in future releases.",
        FutureWarning,
        stacklevel=2,
    )
    from tespy.connections.connection import connection_registry

    schema = {}

    for class_name, cls in connection_registry.items.items():
        if class_name in exclude_base_classes:
            continue

        # Connection classes need two component arguments; use a lightweight
        # approach: call get_parameters() on a bare (uninitialized) instance
        # via a staticmethod check, or fall back to calling the unbound method
        # with a sentinel self.
        try:
            import inspect as _inspect
            static_gp = _inspect.getattr_static(cls, "get_parameters", None)
            if isinstance(static_gp, staticmethod):
                parameters = cls.get_parameters()
            else:
                # Connection.get_parameters() only accesses self for method
                # lookup; create a hollow proxy to avoid full __init__
                instance = cls.__new__(cls)
                parameters = cls.get_parameters(instance)
        except Exception as exc:
            schema[class_name] = {"error": str(exc)}
            continue

        params_list = [
            _serialize_parameter(param_name, dc)
            for param_name, dc in parameters.items()
        ]

        schema[class_name] = {
            "module": cls.__module__,
            "parameters": params_list,
        }

    if as_json:
        return json.dumps(schema, indent=indent)
    return schema
