# -*- coding: utf-8

"""Tests for tespy.tools.schema - component and connection schema generation.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tests/test_tools/test_schema.py

SPDX-License-Identifier: MIT
"""
import json

import pytest

from tespy.components.component import component_registry
from tespy.connections.connection import connection_registry
from tespy.tools.schema import _dc_type_name
from tespy.tools.schema import _serialize_parameter
from tespy.tools.schema import generate_component_schema
from tespy.tools.schema import generate_connection_schema

# ---------------------------------------------------------------------------
# Fixtures - generate schemas once per session for speed
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def component_schema():
    return generate_component_schema(as_json=False)


@pytest.fixture(scope="session")
def connection_schema():
    return generate_connection_schema(as_json=False)


# ---------------------------------------------------------------------------
# Return-type tests
# ---------------------------------------------------------------------------

class TestReturnTypes:
    def test_component_schema_returns_json_string_by_default(self):
        result = generate_component_schema()
        assert isinstance(result, str)
        parsed = json.loads(result)
        assert isinstance(parsed, dict)

    def test_component_schema_returns_dict_when_requested(self, component_schema):
        assert isinstance(component_schema, dict)

    def test_connection_schema_returns_json_string_by_default(self):
        result = generate_connection_schema()
        assert isinstance(result, str)
        parsed = json.loads(result)
        assert isinstance(parsed, dict)

    def test_connection_schema_returns_dict_when_requested(self, connection_schema):
        assert isinstance(connection_schema, dict)


# ---------------------------------------------------------------------------
# Component schema - presence and structure
# ---------------------------------------------------------------------------

class TestComponentSchemaStructure:
    def test_base_component_excluded_by_default(self, component_schema):
        assert "Component" not in component_schema

    def test_known_components_present(self, component_schema):
        expected = set(component_registry.items.keys()) - {"Component"}
        for name in expected:
            assert name in component_schema, f"{name} missing from component schema"

    def test_no_component_has_error(self, component_schema):
        errors = {k: v for k, v in component_schema.items() if "error" in v}
        assert errors == {}, f"Schema generation errors: {errors}"

    def test_component_entry_has_required_keys(self, component_schema):
        required = {
            "module", "submodule", "inlets", "outlets",
            "powerinlets", "poweroutlets",
            "heatinlets", "heatoutlets",
            "parameters"
        }
        for name, entry in component_schema.items():
            assert required <= entry.keys(), (
                f"{name} is missing keys: {required - entry.keys()}"
            )

    def test_parameters_is_a_list(self, component_schema):
        for name, entry in component_schema.items():
            assert isinstance(entry["parameters"], list), (
                f"{name}.parameters is not a list"
            )

    def test_each_parameter_has_required_keys(self, component_schema):
        required = {"name", "container_type", "description"}
        for comp_name, entry in component_schema.items():
            for param in entry["parameters"]:
                assert required <= param.keys(), (
                    f"{comp_name}.{param.get('name')} missing keys: "
                    f"{required - param.keys()}"
                )

    def test_submodule_extracted_correctly(self, component_schema):
        assert component_schema["Turbine"]["submodule"] == "turbomachinery"
        assert component_schema["HeatExchanger"]["submodule"] == "heat_exchangers"
        assert component_schema["Source"]["submodule"] == "basics"
        assert component_schema["Merge"]["submodule"] == "nodes"
        assert component_schema["Motor"]["submodule"] == "power"

    def test_exclusion_respected(self):
        schema = generate_component_schema(
            exclude_base_classes=("Turbine", "Compressor"), as_json=False
        )
        assert "Turbine" not in schema
        assert "Compressor" not in schema
        assert "Pump" in schema


# ---------------------------------------------------------------------------
# Component schema - specific parameter details
# ---------------------------------------------------------------------------

class TestComponentParameters:
    def _get_param(self, component_schema, component, param_name):
        params = {p["name"]: p for p in component_schema[component]["parameters"]}
        assert param_name in params, (
            f"Parameter '{param_name}' not found in {component}. "
            f"Available: {list(params)}"
        )
        return params[param_name]

    def test_turbine_eta_s_is_component_property(self, component_schema):
        p = self._get_param(component_schema, "Turbine", "eta_s")
        assert p["container_type"] == "ComponentProperty"
        assert p["quantity"] == "efficiency"

    def test_polynomial_compressor_has_grouped_properties(self, component_schema):
        p = self._get_param(component_schema, "PolynomialCompressor", "eta_s_group")
        assert p["container_type"] == "GroupedComponentProperties"
        assert "elements" in p
        assert isinstance(p["elements"], list)
        assert len(p["elements"]) > 0

    def test_turbine_pr_is_component_property(self, component_schema):
        p = self._get_param(component_schema, "Turbine", "pr")
        assert p["container_type"] == "ComponentProperty"

    def test_motor_eta_is_efficiency(self, component_schema):
        p = self._get_param(component_schema, "Motor", "eta")
        assert p["container_type"] == "ComponentProperty"
        assert p["quantity"] == "efficiency"

    def test_motor_eta_char_is_characteristic(self, component_schema):
        p = self._get_param(component_schema, "Motor", "eta_char")
        assert p["container_type"] == "ComponentCharacteristic"

    def test_subsystem_interface_num_inter_is_simple(self, component_schema):
        p = self._get_param(component_schema, "SubsystemInterface", "num_inter")
        assert p["container_type"] == "SimpleDataContainer"


# ---------------------------------------------------------------------------
# Port schema
# ---------------------------------------------------------------------------

class TestPortSchema:
    def test_fixed_ports_have_ports_list(self, component_schema):
        inlets = component_schema["CycleCloser"]["inlets"]
        assert inlets["type"] == "fixed"
        assert "ports" in inlets
        assert inlets["ports"] == ["in1"]

    def test_variable_ports_have_parameter_and_pattern(self, component_schema):
        inlets = component_schema["Merge"]["inlets"]
        assert inlets["type"] == "variable"
        assert inlets["parameter"] == "num_in"
        assert "pattern" in inlets
        assert "min" in inlets

    def test_subsystem_interface_variable_in_and_out(self, component_schema):
        entry = component_schema["SubsystemInterface"]
        assert entry["inlets"]["type"] == "variable"
        assert entry["outlets"]["type"] == "variable"
        assert entry["inlets"]["parameter"] == "num_inter"
        assert entry["outlets"]["parameter"] == "num_inter"

    def test_splitter_has_variable_outlets(self, component_schema):
        outlets = component_schema["Splitter"]["outlets"]
        assert outlets["type"] == "variable"
        assert outlets["parameter"] == "num_out"

    def test_motor_has_power_ports_only(self, component_schema):
        entry = component_schema["Motor"]
        assert entry["inlets"]["ports"] == []
        assert entry["outlets"]["ports"] == []
        assert entry["powerinlets"]["type"] == "fixed"
        assert entry["powerinlets"]["ports"] == ["power_in"]
        assert entry["poweroutlets"]["type"] == "fixed"
        assert entry["poweroutlets"]["ports"] == ["power_out"]
        assert entry["heatinlets"]["ports"] == []
        assert entry["heatoutlets"]["ports"] == []

    def test_source_has_no_inlets(self, component_schema):
        inlets = component_schema["Source"]["inlets"]
        assert inlets["type"] == "fixed"
        assert inlets["ports"] == []

    def test_sink_has_no_outlets(self, component_schema):
        outlets = component_schema["Sink"]["outlets"]
        assert outlets["type"] == "fixed"
        assert outlets["ports"] == []

    def test_simple_heat_exchanger_heat_and_power_ports(self, component_schema):
        entry = component_schema["SimpleHeatExchanger"]
        assert entry["powerinlets"]["ports"] == ["heat"]
        assert entry["poweroutlets"]["ports"] == ["heat"]
        assert entry["heatinlets"]["ports"] == ["heat"]
        assert entry["heatoutlets"]["ports"] == ["heat"]

    def test_heat_sink_has_heat_inlet(self, component_schema):
        entry = component_schema["HeatSink"]
        assert entry["heatinlets"]["type"] == "fixed"
        assert entry["heatinlets"]["ports"] == ["heat"]
        assert entry["powerinlets"]["ports"] == []

    def test_heat_source_has_heat_outlet(self, component_schema):
        entry = component_schema["HeatSource"]
        assert entry["heatoutlets"]["type"] == "fixed"
        assert entry["heatoutlets"]["ports"] == ["heat"]
        assert entry["poweroutlets"]["ports"] == []

    def test_power_sink_has_power_inlet(self, component_schema):
        entry = component_schema["PowerSink"]
        assert entry["powerinlets"]["ports"] == ["power"]
        assert entry["heatinlets"]["ports"] == []

    def test_power_source_has_power_outlet(self, component_schema):
        entry = component_schema["PowerSource"]
        assert entry["poweroutlets"]["ports"] == ["power"]
        assert entry["heatoutlets"]["ports"] == []

    def test_heat_bus_has_variable_heat_ports(self, component_schema):
        entry = component_schema["HeatBus"]
        assert entry["heatinlets"]["type"] == "variable"
        assert entry["heatinlets"]["pattern"] == "heat_in{n}"
        assert entry["heatoutlets"]["type"] == "variable"
        assert entry["heatoutlets"]["pattern"] == "heat_out{n}"
        assert entry["powerinlets"]["ports"] == []
        assert entry["poweroutlets"]["ports"] == []

    def test_turbine_has_no_power_ports_by_default(self, component_schema):
        entry = component_schema["Turbine"]
        assert entry["powerinlets"]["type"] == "fixed"
        assert entry["powerinlets"]["ports"] == []
        assert entry["poweroutlets"]["type"] == "fixed"


# ---------------------------------------------------------------------------
# Connection schema
# ---------------------------------------------------------------------------

class TestConnectionSchema:
    def test_connection_types_present(self, connection_schema):
        expected = set(connection_registry.items.keys())
        for name in expected:
            assert name in connection_schema, f"{name} missing from component schema"

    def test_connection_entry_has_required_keys(self, connection_schema):
        for name, entry in connection_schema.items():
            assert "module" in entry
            assert "parameters" in entry

    def test_no_connection_has_error(self, connection_schema):
        errors = {k: v for k, v in connection_schema.items() if "error" in v}
        assert errors == {}, f"Schema generation errors: {errors}"

    def test_connection_has_fluid_parameters(self, connection_schema):
        params = {p["name"] for p in connection_schema["Connection"]["parameters"]}
        assert {"m", "p", "h", "T", "x", "fluid"} <= params

    def test_power_connection_has_E_parameter(self, connection_schema):
        params = {p["name"]: p for p in connection_schema["PowerConnection"]["parameters"]}
        assert "E" in params
        assert params["E"]["container_type"] == "FluidProperties"

    def test_fluid_composition_container_type(self, connection_schema):
        params = {p["name"]: p for p in connection_schema["Connection"]["parameters"]}
        assert params["fluid"]["container_type"] == "FluidComposition"

    def test_specific_connection_excluded(self):
        exclude = "Connection"
        expected = set(connection_registry.items.keys()) - {exclude}
        schema = generate_connection_schema(
            exclude_base_classes=(exclude,), as_json=False
        )
        for name in expected:
            assert name in schema, f"{name} missing from component schema"
        assert exclude not in schema


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

class TestDcTypeName:
    def test_known_container_types(self):
        from tespy.tools.data_containers import ComponentCharacteristicMaps
        from tespy.tools.data_containers import ComponentCharacteristics
        from tespy.tools.data_containers import ComponentProperties
        from tespy.tools.data_containers import FluidComposition
        from tespy.tools.data_containers import FluidProperties
        from tespy.tools.data_containers import GroupedComponentCharacteristics
        from tespy.tools.data_containers import GroupedComponentProperties
        from tespy.tools.data_containers import SimpleDataContainer

        cases = [
            # TODO: Name of container?!
            (ComponentProperties(), "ComponentProperty"),
            (ComponentCharacteristics(), "ComponentCharacteristic"),
            (ComponentCharacteristicMaps(), "ComponentCharacteristicMap"),
            (GroupedComponentProperties(), "GroupedComponentProperties"),
            (GroupedComponentCharacteristics(), "GroupedComponentCharacteristics"),
            (FluidProperties(), "FluidProperties"),
            (FluidComposition(), "FluidComposition"),
            (SimpleDataContainer(), "SimpleDataContainer"),
        ]
        for instance, expected in cases:
            assert _dc_type_name(instance) == expected, (
                f"Expected {expected} for {type(instance).__name__}"
            )

    def test_unknown_type_falls_back_to_class_name(self):
        class MyContainer:
            pass
        assert _dc_type_name(MyContainer()) == "MyContainer"


class TestSerializeParameter:
    def test_quantity_included_when_present(self):
        from tespy.tools.data_containers import ComponentProperties as dc_cp
        dc = dc_cp(quantity="pressure")
        result = _serialize_parameter("pr", dc)
        assert result["quantity"] == "pressure"
        assert result["name"] == "pr"
        assert result["container_type"] == "ComponentProperty"

    def test_quantity_omitted_when_absent(self):
        from tespy.tools.data_containers import SimpleDataContainer as dc_simple
        dc = dc_simple()
        result = _serialize_parameter("flag", dc)
        assert "quantity" not in result

    def test_grouped_includes_elements(self):
        from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
        dc = dc_gcp(elements=["a", "b", "c"])
        result = _serialize_parameter("group", dc)
        assert result["elements"] == ["a", "b", "c"]

    def test_non_grouped_has_no_elements(self):
        from tespy.tools.data_containers import ComponentProperties as dc_cp
        dc = dc_cp()
        result = _serialize_parameter("pr", dc)
        assert "elements" not in result
