import numpy as np
import pytest
from matplotlib.colors import to_hex

from tespy.plotting import plot_QT_diagram
from tespy.plotting import plot_QT_diagram_plotly
from tespy.plotting import qt_diagram_data
from tespy.plotting.style import get_style


def _series(description, role, kind):
    return [s for s in description.series if s.role == role and s.kind == kind]


def test_condenser_curves_span_the_inlet_and_outlet_states(heat_pump):
    condenser = heat_pump.get_comp("condenser")
    description = qt_diagram_data(condenser)
    assert description.kind == "qt_diagram"
    hot = _series(description, "hot", "line")[0]
    cold = _series(description, "cold", "line")[0]
    assert max(hot.y) == pytest.approx(heat_pump.get_conn("2").T.val)
    assert min(hot.y) == pytest.approx(heat_pump.get_conn("3").T.val)
    assert min(cold.y) == pytest.approx(heat_pump.get_conn("w1").T.val)
    assert max(cold.y) == pytest.approx(heat_pump.get_conn("w2").T.val)
    assert max(hot.x) == pytest.approx(abs(condenser.Q.val))


def test_units_are_network_units(heat_pump):
    description = qt_diagram_data(heat_pump.get_comp("condenser"))
    assert description.x_axis.unit == "kW"
    assert description.y_axis.unit == "°C"


def test_moving_boundary_sections_are_marked(heat_pump):
    condenser = heat_pump.get_comp("condenser")
    markers = _series(qt_diagram_data(condenser), "hot", "markers")[0]
    np.testing.assert_allclose(markers.x, condenser.Q_sections.val, atol=1e-9)


def test_plain_heat_exchanger_shows_the_phase_change(heat_pump):
    evaporator = heat_pump.get_comp("evaporator")
    assert len(evaporator.Q_sections.val) == 2
    markers = _series(qt_diagram_data(evaporator), "cold", "markers")[0]
    assert len(markers.x) == 3
    assert any(abs(T) < 1e-6 for T in markers.y)


def test_pinch_is_annotated(heat_pump):
    condenser = heat_pump.get_comp("condenser")
    description = qt_diagram_data(condenser)
    pinch = next(a for a in description.annotations if a.kind == "segment")
    assert pinch.y[1] - pinch.y[0] == pytest.approx(condenser.td_pinch.val, abs=1e-6)
    assert "pinch" in pinch.text


def test_hover_fields_and_references(heat_pump):
    hot = _series(qt_diagram_data(heat_pump.get_comp("condenser")), "hot", "line")[0]
    assert {"Q", "T hot", "T cold", "ΔT"} <= set(hot.fields)
    assert hot.reference == {"component": "condenser", "side": "hot"}


def test_composite_of_two_heat_exchangers(heat_pump):
    condenser = heat_pump.get_comp("condenser")
    evaporator = heat_pump.get_comp("evaporator")
    description = qt_diagram_data([condenser, evaporator])
    hot = _series(description, "hot", "line")
    assert len(hot) == 2
    total = abs(condenser.Q.val) + abs(evaporator.Q.val)
    assert max(hot[1].x) == pytest.approx(total)
    boundary = next(a for a in description.annotations if a.kind == "vline")
    assert boundary.x[0] == pytest.approx(abs(condenser.Q.val))


def test_other_components_raise(heat_pump):
    with pytest.raises(TypeError, match="HeatExchanger"):
        qt_diagram_data(heat_pump.get_comp("compressor"))


def test_matplotlib_plot(heat_pump):
    ax = plot_QT_diagram(heat_pump.get_comp("condenser"))
    colors = {to_hex(line.get_color()) for line in ax.lines}
    assert get_style().hot in colors
    assert get_style().cold in colors


def test_plotly_plot(heat_pump):
    figure = plot_QT_diagram_plotly(heat_pump.get_comp("condenser"))
    hot = next(t for t in figure.data if t.meta["role"] == "hot")
    assert "ΔT = " in hot.hovertemplate
    assert hot.customdata[0][-1] == {"component": "condenser", "side": "hot"}
