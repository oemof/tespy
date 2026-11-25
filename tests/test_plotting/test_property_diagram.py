import importlib.metadata
from unittest.mock import MagicMock

import pytest
from conftest import build_heat_pump
from matplotlib import pyplot as plt
from matplotlib.offsetbox import OffsetImage

from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.networks import Network
from tespy.plotting import clear_diagram_cache
from tespy.plotting import plot_property_diagram
from tespy.plotting import plot_property_diagram_plotly
from tespy.plotting import property_diagram_data
from tespy.plotting.description import PlotDescription
from tespy.plotting.render import draw_matplotlib
from tespy.plotting.render import draw_plotly
from tespy.tools.fluid_properties.wrappers import IAPWSWrapper

REFRIGERANT = ["1", "2", "3", "4", "0"]


def _states(description):
    return next(s for s in description.series if s.role == "cycle")


def _one_pipe(fluid, **kwargs):
    nw = Network(iterinfo=False)
    source, sink = Source("source"), Sink("sink")
    pipe = SimpleHeatExchanger("pipe")
    a = Connection(source, "out1", pipe, "in1", label="a")
    b = Connection(pipe, "out1", sink, "in1", label="b")
    nw.add_conns(a, b)
    a.set_attr(fluid=fluid, m=1, p=1e5, T=300, **kwargs)
    pipe.set_attr(pr=1, Q=1e4)
    nw.solve("design")
    return nw


def test_data_of_the_refrigerant_branch(heat_pump):
    description = property_diagram_data(heat_pump, "1")
    assert description.kind == "property_diagram"
    assert description.extra["fluid"] == "HEOS::R290"
    states = _states(description)
    assert sorted(r["connection"] for r in states.references) == sorted(REFRIGERANT)


def test_state_fields_are_in_network_units(heat_pump):
    states = _states(property_diagram_data(heat_pump, "3"))
    index = [r["connection"] for r in states.references].index("3")
    assert states.fields["T"][index] == pytest.approx(heat_pump.get_conn("3").T.val)
    assert states.fields["T"][index] == pytest.approx(47, abs=1e-3)
    assert states.units["T"] == "°C"
    assert states.units["p"] == "bar"


def test_log_ph_axes(heat_pump):
    description = property_diagram_data(heat_pump, "1")
    assert description.x_axis.quantity == "h"
    assert description.y_axis.quantity == "p"
    assert description.y_axis.scale == "log"
    assert description.x_axis.unit == "kJ/kg"
    states = _states(description)
    assert description.y_axis.range[0] < min(states.y)
    assert description.y_axis.range[1] > max(states.y)


def test_Ts_axes(heat_pump):
    description = property_diagram_data(heat_pump, "1", diagram_type="Ts")
    assert description.x_axis.quantity == "s"
    assert description.y_axis.quantity == "T"
    assert description.y_axis.scale == "linear"


def test_paths_follow_the_components(heat_pump):
    paths = property_diagram_data(heat_pump, "1").extra["paths"]
    assert ["1", "2"] in paths
    assert ["2", "3"] in paths
    assert ["4", "0"] in paths
    assert not any("w1" in pair for pair in paths)


def test_water_branch_is_separate(heat_pump):
    description = property_diagram_data(heat_pump, "w1")
    assert description.extra["fluid"] == "HEOS::water"
    assert {r["connection"] for r in _states(description).references} == {"w1", "w2"}


def test_unknown_diagram_type_raises(heat_pump):
    with pytest.raises(ValueError, match="Ts"):
        property_diagram_data(heat_pump, "1", diagram_type="xy")


def test_unknown_connection_raises(heat_pump):
    with pytest.raises(KeyError, match="nonexistent"):
        property_diagram_data(heat_pump, "nonexistent")


def test_mixture_raises():
    nw = _one_pipe({"N2": 0.79, "O2": 0.21})
    with pytest.raises(ValueError, match="mixture"):
        property_diagram_data(nw, "a")


def test_other_property_back_end_raises():
    nw = _one_pipe({"H2O": 1}, fluid_engines={"H2O": IAPWSWrapper})
    with pytest.raises(ValueError, match="CoolProp"):
        property_diagram_data(nw, "a")


def test_coolprop_mixture_is_passed_with_its_basis():
    nw = _one_pipe({"R32[0.7]&R125[0.3]|mass": 1})
    description = property_diagram_data(nw, "a")
    assert description.extra["fluid"] == "HEOS::R32[0.7]&R125[0.3]"
    assert description.extra["basis"] == "mass"


def test_matplotlib_diagram(heat_pump):
    ax = plot_property_diagram(heat_pump, "1")
    assert "kJ/kg" in ax.get_xlabel()
    assert ax.get_yscale() == "log"
    assert len(ax.findobj(OffsetImage)) == 1


def test_matplotlib_diagram_without_logo(heat_pump):
    ax = plot_property_diagram(heat_pump, "1", logo=False)
    assert len(ax.findobj(OffsetImage)) == 0


def test_matplotlib_diagram_into_given_axes(heat_pump):
    _, ax = plt.subplots()
    assert plot_property_diagram(heat_pump, "1", "Ts", ax=ax) is ax


def test_diagram_hook_receives_the_coolplot_diagram(heat_pump):
    from coolplot import PropertyDiagram
    hook = MagicMock()
    plot_property_diagram(heat_pump, "1", diagram_hook=hook)
    hook.assert_called_once()
    assert isinstance(hook.call_args.args[0], PropertyDiagram)


def test_plotly_diagram_hover_and_references(heat_pump):
    figure = plot_property_diagram_plotly(heat_pump, "1")
    states = [t for t in figure.data if t.customdata is not None]
    assert len(states) == 1
    assert len(states[0].customdata) == len(REFRIGERANT)
    assert "T = " in states[0].hovertemplate
    references = [row[-1]["connection"] for row in states[0].customdata]
    assert sorted(references) == sorted(REFRIGERANT)
    assert figure.layout.meta["tespy"]["kind"] == "property_diagram"


def test_plotly_diagram_information_block(heat_pump):
    figure = plot_property_diagram_plotly(heat_pump, "1", info=["heat pump"])
    texts = " ".join(str(a.text) for a in figure.layout.annotations)
    assert f"TESPy {importlib.metadata.version('tespy')}" in texts
    assert "heat pump" in texts
    assert "Lemmon" in texts
    assert len(figure.layout.images) == 2


def test_plotly_diagram_without_logo(heat_pump):
    figure = plot_property_diagram_plotly(heat_pump, "1", logo=False)
    assert len(figure.layout.images) == 0


def test_exported_description_can_be_drawn_again(heat_pump):
    text = property_diagram_data(heat_pump, "1").to_json()
    description = PlotDescription.from_json(text)
    assert draw_matplotlib(description).get_yscale() == "log"
    assert len(draw_plotly(description).data) > 1


def test_region_contains_the_states_with_a_margin(heat_pump):
    T_min, T_max, p_min, p_max = property_diagram_data(heat_pump, "1").extra["tp_limits"]
    connections = [heat_pump.get_conn(label) for label in REFRIGERANT]
    T = [c.T.val_SI for c in connections]
    p = [c.p.val_SI for c in connections]
    assert T_min <= min(T) - 40 and T_max >= max(T) + 40
    assert p_min <= min(p) / 2 and p_max >= max(p) * 2
    # on a coarse grid, so that similar cycles share the region
    assert T_min % 50 == 0 and T_max % 50 == 0


def test_similar_cycles_share_the_region(heat_pump):
    other = build_heat_pump(T_evap=3, T_cond=52)
    assert (
        property_diagram_data(other, "1").extra["tp_limits"]
        == property_diagram_data(heat_pump, "1").extra["tp_limits"]
    )


def test_region_can_be_given(heat_pump):
    limits = [200.0, 450.0, 1e5, 60e5]
    description = property_diagram_data(heat_pump, "1", tp_limits=limits)
    assert description.extra["tp_limits"] == limits
    diagrams = []
    plot_property_diagram(heat_pump, "1", tp_limits=limits, diagram_hook=diagrams.append)
    assert diagrams[0].tp_limits.T_max == 450.0


def test_isolines_are_cached_between_plots(heat_pump):
    clear_diagram_cache()
    diagrams = []
    plot_property_diagram(heat_pump, "1", diagram_hook=diagrams.append)
    cache = diagrams[0].cache
    computed = cache.misses
    assert computed > 0
    plot_property_diagram(heat_pump, "1", style="tespy-dark", diagram_hook=diagrams.append)
    plot_property_diagram_plotly(heat_pump, "1", diagram_hook=diagrams.append)
    assert all(d.cache is cache for d in diagrams)
    isolines = [i for i in diagrams[0].scene.ids() if i.startswith("iso/") and i.count("/") == 2]
    # only the cycle paths may be computed again, no isoline
    assert cache.misses - computed < len(isolines)
    assert cache.hits >= 2 * len(isolines)
    clear_diagram_cache()
    plot_property_diagram(heat_pump, "1", diagram_hook=diagrams.append)
    assert diagrams[-1].cache is not cache


@pytest.mark.parametrize("diagram_type", ["ph", "Ts"])
def test_at_most_eight_lines_per_family(heat_pump, diagram_type):
    diagrams = []
    plot_property_diagram(heat_pump, "1", diagram_type, diagram_hook=diagrams.append)
    families = {}
    for item_id in diagrams[0].scene.ids():
        parts = item_id.split("/")
        if parts[0] == "iso" and len(parts) == 3:
            families.setdefault(parts[1], set()).add(parts[2])
    assert families
    assert all(len(values) <= 8 for values in families.values())
    assert len(families["T" if diagram_type == "ph" else "p"]) >= 5


def test_specific_volume_field_is_not_the_volumetric_flow(heat_pump):
    states = _states(property_diagram_data(heat_pump, "1"))
    index = [r["connection"] for r in states.references].index("2")
    connection = heat_pump.get_conn("2")
    assert states.fields["v"][index] == pytest.approx(connection.vol.val)
    assert states.fields["v"][index] != pytest.approx(connection.v.val)


def test_every_path_is_listed_once(heat_pump):
    paths = [tuple(pair) for pair in property_diagram_data(heat_pump, "1").extra["paths"]]
    assert len(paths) == len(set(paths))


@pytest.mark.parametrize("diagram_type", ["ph", "Ts"])
def test_automatic_badge_location(heat_pump, diagram_type):
    ax = plot_property_diagram(heat_pump, "1", diagram_type, logo_location="auto")
    assert len(ax.findobj(OffsetImage)) == 1
    figure = plot_property_diagram_plotly(heat_pump, "1", diagram_type, logo_location="auto")
    assert len(figure.layout.images) == 2


@pytest.mark.parametrize("location", ["lower right", "upper left", "auto"])
def test_plotly_logo_is_inside_the_badge_box(heat_pump, location):
    figure = plot_property_diagram_plotly(heat_pump, "1", logo_location=location)
    box, image = figure.layout.shapes[0], figure.layout.images[-1]
    assert image.yanchor == "top" and image.xanchor == "left"
    assert box.x0 < image.x and image.x + image.sizex < box.x1
    assert box.y0 < image.y - image.sizey and image.y < box.y1
    # the lines of the badge stay inside the box, below the logo
    lines = [
        a for a in figure.layout.annotations
        if a.xref == "paper" and ("Lemmon" in str(a.text) or "TESPy" in str(a.text))
    ]
    assert len(lines) >= 2
    assert all(box.y0 < a.y < image.y - image.sizey for a in lines)


def test_plotly_logo_is_not_covered_by_the_badge_box(heat_pump):
    """Plotly draws shapes above images: the fill of the box would lie on
    the logo, so the background is an image below the logo instead."""
    figure = plot_property_diagram_plotly(heat_pump, "1")
    box = figure.layout.shapes[0]
    assert box.fillcolor in (None, "rgba(0,0,0,0)")
    background, logo = figure.layout.images
    assert background.sizing == "stretch"
    assert background.x == pytest.approx(box.x0) and background.y == pytest.approx(box.y1)
    assert background.sizex == pytest.approx(box.x1 - box.x0)
    assert background.sizey == pytest.approx(box.y1 - box.y0)
    assert logo.opacity in (None, 1)
