# -*- coding: utf-8

"""Tests of the data container roles.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_tools/test_data_containers.py

SPDX-License-Identifier: MIT
"""
import warnings

import numpy as np
import pytest

from tespy.components import HeatExchanger
from tespy.components import Pump
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Turbine
from tespy.connections import Connection
from tespy.networks import Network
from tespy.tools import data_containers as dc
from tespy.tools.characteristics import CharLine
from tespy.tools.characteristics import CharMap
from tespy.tools.characteristics import EN12900Polynomial
from tespy.tools.characteristics import Polynomial
from tespy.tools.characteristics import characteristic_from_dict
from tespy.tools.units import Units


class TestConfigParameter:

    def test_dtype_is_validated(self):
        param = dc.ConfigParameter(dtype="int")
        param.accept(3)
        assert param.val == 3 and param.is_set
        with pytest.raises(TypeError):
            param.accept("three")

    def test_untyped_accepts_anything(self):
        param = dc.ConfigParameter()
        param.accept({"a": 1})
        assert param.val == {"a": 1}
        param.accept(None)
        assert not param.is_set

    def test_has_no_equation(self):
        assert not hasattr(dc.ConfigParameter(), "num_eq_sets")


class TestOptionalEquation:

    def test_switch_accepts_bool_only(self):
        switch = dc.OptionalEquation(num_eq_sets=1)
        switch.accept(True)
        assert switch.is_set and switch.val is True
        with pytest.raises(TypeError):
            switch.accept(1)
        switch.accept(None)
        assert not switch.is_set

    def test_val_aliases_is_set(self):
        group = dc.ParameterGroup(elements=["a", "b"])
        group.set_attr(val=True)
        assert group.is_set
        assert group._serialize() == {"is_set": True}

    def test_turbine_cone_is_optional_equation(self):
        turbine = Turbine("turbine")
        assert isinstance(turbine.cone, dc.OptionalEquation)
        assert not isinstance(turbine.cone, dc.ParameterGroup)
        turbine.set_attr(cone=True)
        assert turbine.cone.is_set

    def test_heat_exchanger_characteristic_equation_is_optional(self):
        hx = HeatExchanger("hx")
        assert isinstance(hx.UA_char, dc.OptionalEquation)
        assert not isinstance(hx.UA_char, dc.ParameterGroup)
        assert hx.UA_char.elements == ["UA_char1", "UA_char2"]
        assert isinstance(hx.kA_char, dc.OptionalEquation)
        from tespy.components import Pipe
        assert isinstance(Pipe("pipe").darcy_group, dc.ParameterGroup)


class TestArrayParameter:

    def test_result_array_rejects_input(self):
        array = dc.ArrayParameter(quantity="heat")
        assert array.result_only
        with pytest.raises(TypeError):
            array.accept([1, 2, 3])

    def test_input_array_converts_units(self):
        array = dc.ArrayParameter(quantity="pressure", is_result=False)
        assert not array.result_only
        array.accept([1, 2])
        units = Units()
        units.set_defaults(pressure="bar")
        array.set_SI_from_val(units)
        np.testing.assert_allclose(array.val_SI, [1e5, 2e5])
        assert array._serialize(units)["unit"] == "bar"


class TestScalarParameters:

    def test_component_parameter_var_switch(self):
        param = dc.ComponentParameter(_allows_var=True)
        param.accept("var")
        assert param.is_var and param.is_set
        param.accept(3)
        assert not param.is_var and param.val == 3
        with pytest.raises(ValueError):
            dc.ComponentParameter().accept("var")

    def test_starting_value(self):
        param = dc.ConnectionParameter(quantity="pressure")
        param.val0 = 5
        assert param.val0 == 5
        assert np.isnan(dc.ComponentParameter().val0)

    def test_config_parameter_accepts_empty_markers(self):
        param = dc.ConfigParameter(dtype="str")
        param.set_attr(val=np.nan, is_set=False)
        param.set_attr(val=None)
        assert not param.is_set

    def test_result_only_scalar_refuses_values(self):
        hx = HeatExchanger("hx")
        assert hx.lmtd.result_only
        with pytest.raises(TypeError):
            hx.set_attr(lmtd=5)
        hx.set_attr(lmtd=None)
        assert not hx.UA.result_only
        assert not hx.Q.result_only

    def test_group_element_is_input(self):
        from tespy.components import PolynomialCompressor
        compressor = PolynomialCompressor("compressor")
        assert compressor.P.is_result
        assert not compressor.P.result_only
        assert compressor.Q_diss.result_only

    def test_power_connection_reference(self):
        from tespy.components import PowerSink
        from tespy.components import PowerSource
        from tespy.connections import PowerConnection
        from tespy.connections import Ref
        source = PowerSource("source")
        sink = PowerSink("sink")
        e1 = PowerConnection(source, "power", sink, "power", label="e1")
        e2 = PowerConnection(source, "power", sink, "power", label="e2")
        e2.set_attr(E=Ref(e1, 0.5, 0))
        assert e2.E.reference.is_set
        assert not e2.E.is_set
        assert "E_ref" in e2.parameters
        e2.set_attr(E=5)
        assert not e2.E.reference.is_set
        assert e2.E.is_set

    def test_siblings_share_base(self):
        assert issubclass(dc.ComponentParameter, dc.ScalarParameter)
        assert issubclass(dc.ConnectionParameter, dc.ScalarParameter)
        assert not issubclass(dc.ComponentParameter, dc.ConnectionParameter)


class TestCharacteristicData:

    def test_data_is_set_when_provided_equation_is_not(self):
        data = dc.CharacteristicData()
        data.accept(CharLine())
        assert data.is_set
        equation = dc.CharacteristicEquation()
        equation.accept(CharLine())
        assert not equation.is_set
        equation.accept({"char_func": CharLine(), "is_set": True})
        assert equation.is_set

    def test_signature_is_checked(self):
        data = dc.CharacteristicData()
        data.accept(CharLine())
        with pytest.raises(ValueError):
            data.accept(CharMap())
        with pytest.raises(ValueError):
            dc.CharacteristicData(input_quantities=[None, None]).accept(
                CharLine()
            )

    def test_any_representation_with_matching_signature(self):
        data = dc.CharacteristicData(input_quantities=[None, None])
        data.accept(EN12900Polynomial([1] + [0] * 9))
        assert data.evaluate(3, 4) == 1
        data.accept(CharMap())
        assert data.evaluate(0.5, 0.5) == 1

    def test_dimensioned_axes_follow_network_units(self):
        units = Units()
        units.set_defaults(
            volumetric_flow="m3/h", pressure_difference="bar"
        )
        data = dc.CharacteristicEquation(
            param="v", char_params={"type": "abs"},
            output_quantity="pressure_difference"
        )
        data.accept(CharLine(x=[0, 3600], y=[2, 1]))
        data._bind_units(units, "pump: flow_char")
        with pytest.warns(UserWarning, match="pump: flow_char"):
            assert data.evaluate(0.5) == pytest.approx(1.5e5)
        assert data.input_bounds(0) == pytest.approx((0, 1))

    def test_si_network_does_not_warn(self):
        data = dc.CharacteristicEquation(
            param="v", char_params={"type": "abs"},
            output_quantity="pressure_difference"
        )
        data.accept(CharLine(x=[0, 1], y=[2e5, 1e5]))
        data._bind_units(Units())
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            assert data.evaluate(0.5) == pytest.approx(1.5e5)

    def test_default_data_are_si(self):
        units = Units()
        units.set_defaults(volumetric_flow="m3/h")
        data = dc.CharacteristicEquation(param="v", char_params={"type": "abs"})
        data.set_default(CharLine(x=[0, 1], y=[0, 1]))
        data._bind_units(units)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            assert data.evaluate(0.5) == pytest.approx(0.5)

    def test_axis_units_of_the_data_convert(self):
        units = Units()
        data = dc.CharacteristicEquation(
            param="v", char_params={"type": "abs"},
            output_quantity="pressure_difference"
        )
        data.accept(CharLine(x=[0, 3600], y=[2, 1], x_unit="m3/h", y_unit="bar"))
        data._bind_units(units)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            assert data.evaluate(0.5) == pytest.approx(1.5e5)
        assert data.input_bounds(0) == pytest.approx((0, 1))
        exported = data._serialize(units)["char_func"]
        assert exported["x_unit"] == "m3/h"
        assert exported["y_unit"] == "bar"

    def test_carried_units_are_read(self):
        units = Units()
        units.set_defaults(temperature="K")
        data = dc.CharacteristicData(
            input_quantities=["temperature", "temperature"]
        )
        data.accept(EN12900Polynomial([0, 1, 0, 0, 0, 0, 0, 0, 0, 0]))
        data._bind_units(units)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            assert data.evaluate(283.15, 300) == pytest.approx(10)
        kelvin = EN12900Polynomial([0, 1, 0, 0, 0, 0, 0, 0, 0, 0])
        kelvin.input_units = ["K", "K"]
        data.accept(kelvin)
        assert data.evaluate(283.15, 300) == pytest.approx(283.15)

    def test_quantity_arrays_carry_units(self):
        units = Units()
        units.set_defaults(volumetric_flow="m3/h", pressure_difference="bar")
        Q_ = units.ureg.Quantity
        line = CharLine(x=Q_([0, 1], "m3/s"), y=Q_([2e5, 1e5], "Pa"))
        assert line.x_unit == str(Q_(1, "m3/s").units)
        assert line.input_units == [line.x_unit]
        data = dc.CharacteristicEquation(
            param="v", char_params={"type": "abs"},
            output_quantity="pressure_difference"
        )
        data.accept(line)
        data._bind_units(units)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            assert data.evaluate(0.5) == pytest.approx(1.5e5)
        restored = characteristic_from_dict(line._serialize())
        assert restored.x_unit == line.x_unit
        assert restored.y_unit == line.y_unit
        char_map = CharMap(x_unit="Hz", y_unit="m3/h", z_unit="m")
        assert char_map.input_units == ["Hz", "m3/h"]
        assert char_map.output_unit == "m"

    def test_relative_axes_are_dimensionless(self):
        data = dc.CharacteristicEquation(param="v")
        assert data._resolved_input_quantities() == [None]
        data = dc.CharacteristicEquation(param="v", char_params={"type": "abs"})
        assert data._resolved_input_quantities() == ["volumetric_flow"]

    def test_serialization_roundtrip(self):
        for char in [CharLine(), CharMap(), EN12900Polynomial([1] * 10),
                     Polynomial([1, 2], exponents=[(0,), (1,)])]:
            restored = characteristic_from_dict(char._serialize())
            assert type(restored) is type(char)
            assert restored._serialize() == char._serialize()

    def test_legacy_serialization_without_type(self):
        assert isinstance(characteristic_from_dict({"x": [0, 1], "y": [1, 1]}), CharLine)
        assert isinstance(
            characteristic_from_dict({"x": [0, 1], "y": [[0, 1]] * 2, "z": [[0, 1]] * 2}),
            CharMap
        )


class TestPolynomials:

    def test_en12900_fixed_form(self):
        poly = EN12900Polynomial(range(10))
        t_e, t_c = 2.0, 3.0
        expected = (
            0 + 1 * t_e + 2 * t_c + 3 * t_e**2 + 4 * t_e * t_c + 5 * t_c**2
            + 6 * t_e**3 + 7 * t_e**2 * t_c + 8 * t_e * t_c**2 + 9 * t_c**3
        )
        assert poly.evaluate(t_e, t_c) == pytest.approx(expected)
        with pytest.raises(ValueError):
            EN12900Polynomial([1, 2, 3])

    def test_general_polynomial(self):
        poly = Polynomial([1, 2, 3], exponents=[(0, 0), (1, 0), (0, 2)])
        assert poly.num_inputs == 2
        assert poly.evaluate(2, 3) == pytest.approx(1 + 4 + 27)
        with pytest.raises(ValueError):
            poly.evaluate(1)
        assert poly.input_bounds(0) is None


class TestPumpFlowCharacteristicUnits:

    def _network(self, **defaults):
        nw = Network(iterinfo=False)
        nw.units.set_defaults(
            pressure="bar", temperature="degC", **defaults
        )
        source = Source("source")
        sink = Sink("sink")
        pump = Pump("pump")
        inlet = Connection(source, "out1", pump, "in1")
        outlet = Connection(pump, "out1", sink, "in1")
        nw.add_conns(inlet, outlet)
        inlet.set_attr(fluid={"water": 1}, p=1, T=20, v=1.5)
        pump.set_attr(eta_s=0.8)
        return nw, pump, inlet, outlet

    def test_flow_char_follows_network_units(self):
        nw, pump, inlet, outlet = self._network(
            volumetric_flow="m3/h", pressure_difference="bar"
        )
        inlet.set_attr(v=1.5 * 3600)
        char = CharLine(
            x=np.array([0, 0.4, 0.8, 1.2, 1.6, 2]) * 3600,
            y=np.array([15, 14, 12, 9, 5, 0])
        )
        pump.set_attr(flow_char={"char_func": char, "is_set": True})
        with pytest.warns(UserWarning, match="pump: flow_char"):
            nw.solve("design")
        nw.assert_convergence()
        assert outlet.p.val - inlet.p.val == pytest.approx(6, abs=1e-3)

    def test_flow_char_with_declared_units(self):
        nw, pump, inlet, outlet = self._network(
            volumetric_flow="m3/h", pressure_difference="bar"
        )
        inlet.set_attr(v=1.5 * 3600)
        char = CharLine(
            x=np.array([0, 0.4, 0.8, 1.2, 1.6, 2]) * 3600,
            y=np.array([15, 14, 12, 9, 5, 0]), x_unit="m3/h", y_unit="bar"
        )
        pump.set_attr(flow_char={"char_func": char, "is_set": True})
        nw.solve("design")
        nw.assert_convergence()
        assert outlet.p.val - inlet.p.val == pytest.approx(6, abs=1e-3)
        exported = pump.flow_char._serialize(nw.units)["char_func"]
        assert exported["x_unit"] == "m3/h"


class TestDeprecatedNames:

    def test_old_names_warn_and_build_new_containers(self):
        with pytest.warns(FutureWarning):
            param = dc.ComponentProperties(quantity="pressure")
        assert isinstance(param, dc.ComponentParameter)
        with pytest.warns(FutureWarning):
            switch = dc.GroupedComponentCharacteristics()
        assert isinstance(switch, dc.OptionalEquation)
        assert not isinstance(switch, dc.ParameterGroup)
        with pytest.warns(FutureWarning):
            char_map = dc.ComponentCharacteristicMaps()
        assert char_map.num_inputs == 2

    def test_old_names_recognise_new_instances(self):
        hx = HeatExchanger("hx")
        assert isinstance(hx.UA, dc.ComponentProperties)
        assert isinstance(hx.UA_char1, dc.ComponentCharacteristics)
        assert not isinstance(hx.UA_char1, dc.ComponentCharacteristicMaps)
        assert isinstance(hx.Q_sections, dc.ComponentArrayProperties)
        assert not isinstance(hx.UA, dc.FluidProperties)
        turbine = Turbine("turbine")
        assert isinstance(turbine.cone, dc.GroupedComponentCharacteristics)
        assert not isinstance(turbine.cone, dc.GroupedComponentProperties)

    def test_result_only_is_available(self):
        hx = HeatExchanger("hx")
        assert hx.Q_sections.result_only
        assert not hx.UA.result_only
        assert not hx.UA_char1.result_only
