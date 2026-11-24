import json
import math

import pytest

from tespy.plotting.description import FORMAT_VERSION
from tespy.plotting.description import Annotation
from tespy.plotting.description import Axis
from tespy.plotting.description import PlotDescription
from tespy.plotting.description import Series


def _description():
    return PlotDescription(
        kind="qt_diagram",
        title="condenser",
        x_axis=Axis("heat transferred", "kW"),
        y_axis=Axis("temperature", "°C"),
        series=[
            Series(
                "hot side", x=[0.0, 10.0], y=[60.0, 80.0], role="hot",
                fields={"T": [60.0, 80.0], "dT": [5.0, float("nan")]},
                units={"T": "°C", "dT": "K"},
                references=[{"component": "condenser"}] * 2,
                reference={"component": "condenser", "side": "hot"}
            )
        ],
        annotations=[Annotation("segment", x=[0, 0], y=[55, 60], text="pinch 5 K")],
        info=["heat pump"],
    )


def test_axis_title_with_and_without_unit():
    assert Axis("temperature", "°C").title == "temperature in °C"
    assert Axis("COP").title == "COP"


def test_axis_scale_is_checked():
    with pytest.raises(ValueError, match="log"):
        Axis("p", "bar", scale="logarithmic")


def test_series_kind_and_role_are_checked():
    with pytest.raises(ValueError, match="line"):
        Series("a", x=[0], y=[0], kind="pie")
    with pytest.raises(ValueError, match="hot"):
        Series("a", x=[0], y=[0], role="lukewarm")


def test_series_fields_need_one_value_per_point():
    with pytest.raises(ValueError, match="T"):
        Series("a", x=[0, 1], y=[0, 1], fields={"T": [1.0]})


def test_dict_has_format_meta():
    data = _description().to_dict()
    assert data["meta"] == {"format": "tespy-plot", "version": FORMAT_VERSION}
    assert data["kind"] == "qt_diagram"


def test_json_round_trip():
    description = _description()
    restored = PlotDescription.from_json(description.to_json())
    assert restored.series[0].references == description.series[0].references
    assert restored.annotations[0].text == "pinch 5 K"
    assert restored.x_axis.title == "heat transferred in kW"


def test_json_has_no_nan_literals():
    text = _description().to_json()
    assert "NaN" not in text
    data = json.loads(text)
    assert data["series"][0]["fields"]["dT"][1] is None


def test_from_dict_rejects_other_formats():
    data = _description().to_dict()
    data["meta"]["format"] = "something-else"
    with pytest.raises(ValueError, match="tespy-plot"):
        PlotDescription.from_dict(data)


def test_from_dict_rejects_newer_versions():
    data = _description().to_dict()
    data["meta"]["version"] = FORMAT_VERSION + 1
    with pytest.raises(ValueError, match="version"):
        PlotDescription.from_dict(data)


def test_numpy_values_are_exported_as_floats():
    import numpy as np
    description = PlotDescription(
        kind="sweep", x_axis=Axis("x"), y_axis=Axis("y"),
        series=[Series("a", x=np.array([1, 2]), y=np.array([np.inf, 2.0]))]
    )
    data = json.loads(description.to_json())
    assert data["series"][0]["x"] == [1.0, 2.0]
    assert data["series"][0]["y"][0] is None
    assert not math.isnan(description.series[0].y[1])
