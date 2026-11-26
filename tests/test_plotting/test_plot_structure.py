import pytest
from conftest import build_heat_pump
from matplotlib.colors import to_hex
from matplotlib.patches import FancyBboxPatch

from tespy.plotting import incidence_matrix_data
from tespy.plotting import plot_incidence_matrix
from tespy.plotting import plot_incidence_matrix_plotly
from tespy.plotting.structure import COUPLED_BLOCK
from tespy.plotting.structure import EARLIER_BLOCK
from tespy.plotting.structure import EMPTY
from tespy.plotting.structure import ENTRY
from tespy.plotting.style import get_style


@pytest.fixture(scope="module")
def coupled():
    """Heat pump with a pinch specification, which couples three equations."""
    nw = build_heat_pump()
    nw.get_conn("3").set_attr(T_bubble=None)
    nw.get_comp("condenser").set_attr(td_pinch=5)
    nw.solve("design")
    assert nw.status == 0
    return nw


def _entries(description):
    """Cells with a dependency and their kind."""
    z = description.series[0].z
    return {
        (i, j): z[i][j] for i in range(len(z)) for j in range(len(z[i]))
        if z[i][j] != EMPTY
    }


def test_incidence_matrix_has_one_entry_per_dependency(coupled):
    incidence = coupled.problem.get_incidence()
    description = incidence_matrix_data(coupled)
    assert description.kind == "incidence_matrix"
    z = description.series[0].z
    assert len(z) == len(incidence)
    assert len(z[0]) == len({v for deps in incidence.values() for v in deps})
    entries = _entries(description)
    assert len(entries) == sum(len(deps) for deps in incidence.values())
    assert set(entries.values()) == {ENTRY}


def test_rows_and_columns_are_labeled(coupled):
    description = incidence_matrix_data(coupled)
    equations, variables = description.y_axis.tick_labels, description.x_axis.tick_labels
    assert len(equations) == len(description.series[0].z)
    assert len(set(equations)) == len(equations)
    assert len(set(variables)) == len(variables)
    assert any(label.startswith("condenser.") for label in equations)
    # the first equation is at the top
    assert description.y_axis.range[0] > description.y_axis.range[1]


def test_block_order_is_lower_triangular(coupled):
    description = incidence_matrix_data(coupled, block_order=True)
    blocks = description.extra["blocks"]
    assert [b["size"] for b in blocks] == [3, 1, 1]
    block_of_row = {
        row: b for b in blocks
        for row in range(b["start"], b["start"] + b["size"])
    }
    entries = _entries(description)
    for (row, column), kind in entries.items():
        block = block_of_row[row]
        assert column < block["start"] + block["size"]
        if column < block["start"]:
            assert kind == EARLIER_BLOCK
        elif block["size"] > 1:
            assert kind == COUPLED_BLOCK
        else:
            assert kind == ENTRY
    assert {EARLIER_BLOCK, COUPLED_BLOCK, ENTRY} == set(entries.values())


def test_matplotlib_plot(coupled):
    style = get_style("tespy")
    description = incidence_matrix_data(coupled, block_order=True)
    ax = plot_incidence_matrix(coupled, block_order=True)
    cells = [p for p in ax.patches if isinstance(p, FancyBboxPatch) and p.get_fill()]
    z = description.series[0].z
    assert len(cells) == len(z) * len(z[0])
    colors = {to_hex(p.get_facecolor()) for p in cells}
    assert colors == {style.grid, style.sequential[1], style.color(0), style.color(1)}
    # the colors mark kinds of entries: no color bar, no frame
    assert len(ax.figure.axes) == 1
    assert not any(spine.get_visible() for spine in ax.spines.values())
    assert ax.yaxis_inverted()
    labels = [t.get_text() for t in ax.get_yticklabels()]
    assert labels == description.y_axis.tick_labels


def test_plotly_cells_are_rounded_shapes_in_the_style_colors(coupled):
    style = get_style("tespy")
    description = incidence_matrix_data(coupled, block_order=True)
    figure = plot_incidence_matrix_plotly(coupled, block_order=True)
    shapes = [
        shape for shape in figure.layout.shapes
        if shape.type == "path" and shape.line.dash != "dash"
    ]
    # one shape per kind of cell, with one rounded square per cell
    assert {shape.fillcolor for shape in shapes} == {
        style.grid, style.sequential[1], style.color(0), style.color(1)
    }
    z = description.series[0].z
    assert sum(shape.path.count("M") for shape in shapes) == len(z) * len(z[0])
    assert all("Q" in shape.path and shape.xref == "x" and shape.yref == "y" for shape in shapes)
    assert all(shape.line.width == 0 for shape in shapes)
    # the heat map only carries the hover text
    trace = figure.data[0]
    assert trace.type == "heatmap" and trace.opacity == 0 and trace.showscale is False


def test_plotly_first_equation_is_at_the_top(coupled):
    """Numbered axes with a reversed y axis: the diagonal runs from the
    upper left to the lower right, also after autoscale."""
    description = incidence_matrix_data(coupled, block_order=True)
    figure = plot_incidence_matrix_plotly(coupled, block_order=True)
    trace, layout = figure.data[0], figure.layout
    n = len(description.y_axis.tick_labels)
    assert list(trace.y) == list(range(n)) and list(trace.x) == list(range(n))
    assert layout.yaxis.type == "linear" and layout.xaxis.type == "linear"
    assert layout.yaxis.autorange == "reversed"
    assert layout.xaxis.autorange in (None, True)
    assert layout.yaxis.range is None
    assert list(layout.yaxis.ticktext) == description.y_axis.tick_labels
    assert list(layout.yaxis.tickvals) == list(range(n))
    assert list(layout.xaxis.ticktext) == description.x_axis.tick_labels
    # the button that restores the first view is not needed next to autoscale
    assert "resetScale2d" in layout.modebar.remove


def test_plotly_hover_names_equation_and_variable(coupled):
    description = incidence_matrix_data(coupled, block_order=True)
    trace = plot_incidence_matrix_plotly(coupled, block_order=True).data[0]
    assert "%{text}" in trace.hovertemplate
    assert description.y_axis.tick_labels[0] in trace.text[0][0]
    assert description.x_axis.tick_labels[1] in trace.text[0][1]


def test_plotly_cells_are_square(coupled):
    """The axes are locked to each other, whatever room the labels take."""
    layout = plot_incidence_matrix_plotly(coupled, block_order=True).layout
    assert layout.yaxis.scaleanchor == "x" and layout.yaxis.scaleratio == 1
    assert layout.yaxis.constrain == "domain" and layout.xaxis.constrain == "domain"


def test_coupled_blocks_have_a_boundary(coupled):
    description = incidence_matrix_data(coupled, block_order=True)
    boxes = [a for a in description.annotations if a.kind == "box"]
    assert len(boxes) == 1
    assert boxes[0].x == [-0.5, 2.5] and boxes[0].y == [-0.5, 2.5]
    assert not incidence_matrix_data(coupled).annotations


def test_matplotlib_block_boundary_is_a_dashed_rounded_outline(coupled):
    ax = plot_incidence_matrix(coupled, block_order=True)
    outlines = [
        p for p in ax.patches
        if isinstance(p, FancyBboxPatch) and not p.get_fill()
    ]
    assert len(outlines) == 1
    assert outlines[0].get_linestyle() == "--"
    assert outlines[0].get_x() == pytest.approx(-0.5)
    assert outlines[0].get_width() == pytest.approx(3.0)


def test_plotly_block_boundary_is_a_dashed_rounded_outline(coupled):
    figure = plot_incidence_matrix_plotly(coupled, block_order=True)
    outlines = [s for s in figure.layout.shapes if s.line.dash == "dash"]
    assert len(outlines) == 1
    assert outlines[0].type == "path" and "Q" in outlines[0].path
    assert outlines[0].fillcolor in (None, "rgba(0,0,0,0)")


def test_plotly_logo_is_anchored_to_the_matrix(coupled):
    """The cells are square, so the matrix can be smaller than the plot
    area: the logo and its lines follow the matrix, not the plot area."""
    figure = plot_incidence_matrix_plotly(coupled, block_order=True)
    image = figure.layout.images[0]
    assert (image.xref, image.yref) == ("x domain", "y domain")
    assert image.x > 1 and image.y == pytest.approx(1.0)
    info = next(a for a in figure.layout.annotations if "TESPy" in str(a.text))
    assert (info.xref, info.yref) == ("x domain", "y domain")
    assert info.x == pytest.approx(image.x) and info.y < image.y
