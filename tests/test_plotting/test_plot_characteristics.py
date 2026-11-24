import numpy as np
import pytest

from tespy.components import Compressor
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.networks import Network
from tespy.plotting import characteristic_data
from tespy.plotting import plot_characteristic
from tespy.plotting import plot_characteristic_plotly
from tespy.tools.characteristics import CharLine
from tespy.tools.characteristics import CharMap
from tespy.tools.characteristics import EN12900Polynomial
from tespy.tools.characteristics import Polynomial


def _by_role(description, role):
    return [s for s in description.series if s.role == role]


@pytest.fixture(scope="module")
def compressor_offdesign(tmp_path_factory):
    nw = Network(iterinfo=False)
    source, sink = Source("source"), Sink("sink")
    compressor = Compressor("compressor")
    a = Connection(source, "out1", compressor, "in1", label="a")
    b = Connection(compressor, "out1", sink, "in1", label="b")
    nw.add_conns(a, b)
    a.set_attr(fluid={"air": 1}, m=1, p=1e5, T=293.15)
    compressor.set_attr(
        pr=3, eta_s=0.8, design=["eta_s"], offdesign=["eta_s_char"]
    )
    nw.solve("design")
    path = str(tmp_path_factory.mktemp("design") / "design.json")
    nw.save(path)
    a.set_attr(m=0.7)
    nw.solve("offdesign", design_path=path)
    assert nw.status == 0
    return compressor


def test_char_line_support_points():
    line = CharLine([0, 1, 2], [1, 0.9, 0.7], x_unit="m3/h")
    description = characteristic_data(line)
    support = _by_role(description, "category")[0]
    assert list(support.x) == [0, 1, 2]
    assert list(support.y) == [1, 0.9, 0.7]
    assert description.x_axis.unit == "m3/h"


def test_char_line_valid_range_is_shaded():
    description = characteristic_data(CharLine([0, 1, 2], [1, 0.9, 0.7]))
    regions = [a for a in description.annotations if a.kind == "region"]
    assert len(regions) == 2
    assert all(a.role == "invalid" for a in regions)
    assert regions[0].x[1] == 0
    assert regions[1].x[0] == 2


def test_char_line_operating_points_inside_and_outside():
    line = CharLine([0, 1, 2], [1, 0.9, 0.7])
    description = characteristic_data(line, operating_points=[(1.5,), (3,)])
    inside = _by_role(description, "cycle")[0]
    outside = _by_role(description, "invalid")[0]
    assert inside.x == [1.5]
    assert inside.y == [pytest.approx(0.8)]
    assert outside.x == [3]
    assert description.x_axis.range[1] > 3


def test_char_line_extrapolation_is_dashed():
    line = CharLine([0, 1, 2], [1, 0.9, 0.7], extrapolate=True)
    description = characteristic_data(line, operating_points=[(3,)])
    dashed = [s for s in description.series if s.dash == "dash"]
    assert len(dashed) == 1
    assert dashed[0].x[-1] == 3
    assert dashed[0].y[-1] == pytest.approx(line.evaluate(3))


def test_char_map_one_curve_per_speed_line():
    char_map = CharMap(
        x=[0.9, 1.0, 1.1],
        y=[[0.5, 1.0], [0.6, 1.1], [0.7, 1.2]],
        z=[[0.8, 0.85], [0.82, 0.88], [0.84, 0.9]],
    )
    description = characteristic_data(char_map, operating_points=[(1.0, 0.8)])
    curves = _by_role(description, "scale")
    assert len(curves) == 3
    assert list(curves[0].x) == [0.5, 1.0]
    point = _by_role(description, "cycle")[0]
    assert point.y[0] == pytest.approx(char_map.evaluate(1.0, 0.8))


def test_polynomial_needs_ranges():
    with pytest.raises(ValueError, match="ranges"):
        characteristic_data(Polynomial([1, 0.5], exponents=[(0,), (1,)]))


def test_polynomial_with_one_input():
    polynomial = Polynomial([1, 0.5], exponents=[(0,), (1,)])
    description = characteristic_data(polynomial, ranges=[(0, 4)], points=5)
    line = description.series[0]
    np.testing.assert_allclose(line.x, [0, 1, 2, 3, 4])
    np.testing.assert_allclose(line.y, [1, 1.5, 2, 2.5, 3])


def test_en12900_polynomial_contour():
    polynomial = EN12900Polynomial(
        [1000, 20, -5, 0, 0, 0, 0, 0, 0, 0], value_unit="W"
    )
    description = characteristic_data(
        polynomial, ranges=[(-20, 10), (30, 60)], points=7
    )
    contour = description.series[0]
    assert contour.kind == "contour"
    assert np.shape(contour.z) == (7, 7)
    assert contour.z[0][0] == pytest.approx(polynomial.evaluate(-20, 30))
    assert description.x_axis.label == "T_evap"
    assert description.x_axis.unit == "°C"
    assert description.z_axis.unit == "W"


def test_polynomial_with_three_inputs_raises():
    polynomial = Polynomial([1], exponents=[(1, 1, 1)])
    with pytest.raises(ValueError, match="two inputs"):
        characteristic_data(polynomial, ranges=[(0, 1)] * 3)


def test_component_characteristic_marks_the_operating_point(compressor_offdesign):
    data = compressor_offdesign.eta_s_char
    description = characteristic_data(compressor_offdesign, "eta_s_char")
    point = _by_role(description, "cycle")[0]
    inputs = data._inputs_to_data(compressor_offdesign.get_char_inputs("eta_s_char"))
    assert point.x[0] == pytest.approx(inputs[0])
    assert point.y[0] == pytest.approx(data.char_func.evaluate(*inputs))
    assert point.references[0] == {
        "component": "compressor", "characteristic": "eta_s_char"
    }
    assert "compressor" in description.title


def test_component_needs_the_characteristic_name(compressor_offdesign):
    with pytest.raises(ValueError, match="eta_s_char"):
        characteristic_data(compressor_offdesign)


def test_matplotlib_plot():
    ax = plot_characteristic(CharLine([0, 1, 2], [1, 0.9, 0.7]))
    assert len(ax.lines) >= 1


def test_plotly_plot_hover():
    figure = plot_characteristic_plotly(
        CharLine([0, 1, 2], [1, 0.9, 0.7]), operating_points=[(1.5,)]
    )
    point = next(t for t in figure.data if t.meta["role"] == "cycle")
    assert "x = " in point.hovertemplate


@pytest.mark.parametrize("char", [
    CharLine([0, 1], [1, 1]),
    CharMap(x=[0, 1], y=[[0, 1], [0, 1]], z=[[1, 1], [1, 1]]),
])
def test_old_plot_methods_are_deprecated(char, tmp_path):
    with pytest.warns(FutureWarning, match="0.14"):
        char.plot(str(tmp_path / "char.png"), "title", "x", "y")


def test_polynomial_contour_is_drawn_by_both_renderers():
    polynomial = EN12900Polynomial([1000, 20, -5, 0, 0, 0, 0, 0, 0, 0])
    kwargs = {"ranges": [(-20, 10), (30, 60)], "operating_points": [(-5, 45)]}
    ax = plot_characteristic(polynomial, **kwargs)
    assert len(ax.figure.axes) == 2
    figure = plot_characteristic_plotly(polynomial, **kwargs)
    assert figure.data[0].type == "contour"
