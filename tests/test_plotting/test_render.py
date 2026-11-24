import matplotlib
import pytest

matplotlib.use("Agg")

from matplotlib import pyplot as plt  # noqa: E402
from matplotlib.collections import QuadMesh  # noqa: E402
from matplotlib.colors import to_hex  # noqa: E402
from matplotlib.offsetbox import OffsetImage  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from tespy.plotting import settings  # noqa: E402
from tespy.plotting.description import Annotation  # noqa: E402
from tespy.plotting.description import Axis  # noqa: E402
from tespy.plotting.description import PlotDescription  # noqa: E402
from tespy.plotting.description import Series  # noqa: E402
from tespy.plotting.render import draw_matplotlib  # noqa: E402
from tespy.plotting.render import draw_plotly  # noqa: E402
from tespy.plotting.style import get_style  # noqa: E402


@pytest.fixture(autouse=True)
def _clean():
    settings.reset_defaults()
    yield
    settings.reset_defaults()
    plt.close("all")


def _lines():
    return PlotDescription(
        kind="qt_diagram",
        title="condenser",
        x_axis=Axis("heat transferred", "kW"),
        y_axis=Axis("temperature", "°C"),
        series=[
            Series(
                "hot side", x=[0.0, 10.0], y=[60.0, 80.0], role="hot",
                kind="line+markers", fields={"T": [60.0, 80.0]},
                units={"T": "°C"}, references=[{"connection": "a"}, {"connection": "b"}]
            ),
            Series("cold side", x=[0.0, 10.0], y=[50.0, 55.0], role="cold"),
        ],
        annotations=[
            Annotation("segment", x=[0, 0], y=[50, 60], text="pinch 10 K"),
            Annotation("region", x=[8, 10], y=[None, None], role="invalid"),
        ],
        info=["heat pump"],
    )


def _bars():
    return PlotDescription(
        kind="component_values",
        x_axis=Axis("power", "kW"),
        y_axis=Axis("component", categories=["compressor", "pump", "fan"]),
        series=[
            Series("turbomachinery", x=[10.0, 1.0, 0.5], y=["compressor", "pump", "fan"], kind="bar", slot=0),
        ],
    )


def _heatmap():
    return PlotDescription(
        kind="operating_map",
        x_axis=Axis("T cond", "°C"),
        y_axis=Axis("T evap", "°C"),
        z_axis=Axis("COP"),
        series=[Series("COP", x=[30, 40, 50], y=[0, 10], z=[[4, 3, 2.5], [5, 4, 3]], kind="heatmap")],
    )


def test_matplotlib_lines_use_style_roles():
    style = get_style("tespy")
    ax = draw_matplotlib(_lines())
    lines = {line.get_label(): line for line in ax.lines}
    assert to_hex(lines["hot side"].get_color()) == style.hot
    assert to_hex(lines["cold side"].get_color()) == style.cold
    assert ax.get_xlabel() == "heat transferred in kW"
    assert ax.get_ylabel() == "temperature in °C"
    assert ax.get_title(loc="left") == "condenser"
    assert len(ax.get_legend().get_texts()) == 2


def test_matplotlib_draws_into_given_axes():
    _, ax = plt.subplots()
    assert draw_matplotlib(_lines(), ax=ax) is ax


def test_matplotlib_annotations():
    style = get_style("tespy")
    ax = draw_matplotlib(_lines())
    assert "pinch 10 K" in [t.get_text() for t in ax.texts]
    regions = [p for p in ax.patches if isinstance(p, Rectangle)]
    assert any(to_hex(p.get_facecolor()) == style.invalid for p in regions)


def test_matplotlib_logo_follows_switch():
    assert len(draw_matplotlib(_lines()).findobj(OffsetImage)) == 1
    assert len(draw_matplotlib(_lines(), logo=False).findobj(OffsetImage)) == 0
    settings.set_defaults(logo=False)
    assert len(draw_matplotlib(_lines()).findobj(OffsetImage)) == 0


def test_matplotlib_does_not_change_global_settings():
    before = dict(matplotlib.rcParams)
    draw_matplotlib(_lines(), style="tespy-dark")
    assert matplotlib.rcParams["axes.facecolor"] == before["axes.facecolor"]


def test_matplotlib_dark_style_colors_the_axes():
    style = get_style("tespy-dark")
    ax = draw_matplotlib(_lines(), style=style)
    assert to_hex(ax.get_facecolor()) == style.surface


def test_matplotlib_log_axis():
    description = _lines()
    description.y_axis.scale = "log"
    assert draw_matplotlib(description).get_yscale() == "log"


def test_matplotlib_bars_follow_category_order():
    ax = draw_matplotlib(_bars())
    labels = [t.get_text() for t in ax.get_yticklabels()]
    assert labels == ["compressor", "pump", "fan"]
    assert len(ax.patches) >= 3


def test_matplotlib_heatmap_with_colorbar():
    ax = draw_matplotlib(_heatmap())
    assert any(isinstance(c, QuadMesh) for c in ax.collections)
    colorbar_axes = [a for a in ax.figure.axes if a is not ax]
    assert colorbar_axes[0].get_ylabel() == "COP"


def test_plotly_lines_use_style_roles():
    style = get_style("tespy")
    figure = draw_plotly(_lines())
    traces = {t.name: t for t in figure.data if t.name}
    assert traces["hot side"].line.color == style.hot
    assert traces["cold side"].line.color == style.cold
    assert figure.layout.xaxis.title.text == "heat transferred in kW"
    assert figure.layout.paper_bgcolor == style.surface


def test_plotly_hover_shows_fields_with_units():
    trace = next(t for t in draw_plotly(_lines()).data if t.name == "hot side")
    assert "T = " in trace.hovertemplate
    assert "°C" in trace.hovertemplate
    assert trace.customdata[0][-1] == {"connection": "a"}


def test_plotly_figure_carries_the_description():
    figure = draw_plotly(_lines())
    assert figure.layout.meta["tespy"]["kind"] == "qt_diagram"


def test_plotly_logo_follows_switch():
    assert len(draw_plotly(_lines()).layout.images) == 1
    assert len(draw_plotly(_lines(), logo=False).layout.images) == 0


def test_plotly_info_lines_are_shown():
    texts = [a.text for a in draw_plotly(_lines()).layout.annotations]
    assert any("heat pump" in str(t) for t in texts)


def test_plotly_annotations():
    figure = draw_plotly(_lines())
    assert any("pinch 10 K" in str(a.text) for a in figure.layout.annotations)
    assert any(s.fillcolor == get_style("tespy").invalid for s in figure.layout.shapes)


def test_plotly_log_axis():
    description = _lines()
    description.y_axis.scale = "log"
    assert draw_plotly(description).layout.yaxis.type == "log"


def test_plotly_bars():
    trace = draw_plotly(_bars()).data[0]
    assert trace.type == "bar"
    assert trace.orientation == "h"


def test_plotly_heatmap_uses_style_colorscale():
    style = get_style("tespy")
    trace = draw_plotly(_heatmap(), style=style).data[0]
    assert trace.type == "heatmap"
    assert [c for _, c in trace.colorscale] == list(style.sequential)
