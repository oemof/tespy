# -*- coding: utf-8

"""Draw plot descriptions with matplotlib or Plotly.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/plotting/render.py

SPDX-License-Identifier: MIT
"""
import struct

from tespy.plotting.branding import LOGO_WIDTH_PT
from tespy.plotting.branding import add_logo
from tespy.plotting.branding import info_lines
from tespy.plotting.branding import logo_data_uri
from tespy.plotting.branding import logo_path
from tespy.plotting.settings import resolve
from tespy.plotting.style import DASHES
from tespy.plotting.style import MARKERS
from tespy.plotting.style import get_style

#: Size of new Plotly figures in pixels.
PLOTLY_SIZE = (760, 520)
PLOTLY_MARGIN = {"l": 70, "r": 30, "t": 50, "b": 60}
PX_PER_PT = 4 / 3


def series_color(series, style):
    """Color of a series from its role."""
    if series.role == "category":
        return style.color(series.slot)
    if series.role == "reference":
        return style.muted
    if series.role == "scale":
        from matplotlib.colors import to_hex
        return to_hex(style.colormap()(series.color_value))
    return getattr(style, series.role)


def _annotation_color(annotation, style):
    if annotation.role == "invalid":
        return style.invalid
    return style.ink_secondary


def _legend_labels(description):
    """Legend label per series, one entry per group."""
    seen = set()
    labels = []
    for series in description.series:
        key = series.group or series.name
        if not series.show_legend or not series.name or key in seen:
            labels.append(None)
            continue
        seen.add(key)
        labels.append(key)
    return labels


def draw_matplotlib(description, ax=None, style=None, logo=None, logo_location=None, **kwargs):
    """Draw a plot description with matplotlib.

    Parameters
    ----------
    description : tespy.plotting.description.PlotDescription
        The plot.

    ax : matplotlib.axes.Axes, optional
        Axes to draw into, a new figure if not given.

    style : str, Style, optional
        Style of the plot, the default style if not given.

    logo : bool, optional
        Show the logo and the information block.

    logo_location : str, optional
        Corner of the logo.

    kwargs
        Passed on to the drawing of special plot kinds, e.g.
        :code:`diagram_hook` for property diagrams.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if description.kind == "property_diagram":
        from tespy.plotting.cycles import draw_property_diagram_matplotlib
        return draw_property_diagram_matplotlib(
            description, ax=ax, style=style, logo=logo,
            logo_location=logo_location, **kwargs
        )

    from matplotlib import pyplot as plt

    style = get_style(style)
    logo, location = resolve(logo, logo_location)
    with plt.rc_context(style.rc_params()):
        if ax is None:
            _, ax = plt.subplots(figsize=(7.5, 5), layout="constrained")
        style_matplotlib_axes(ax, style)
        labels = _legend_labels(description)
        for series, label in zip(description.series, labels):
            _draw_series_matplotlib(ax, series, label, style, description)
        _set_axes_matplotlib(ax, description)
        for annotation in description.annotations:
            _draw_annotation_matplotlib(ax, annotation, style)
        if description.title:
            ax.set_title(description.title)
        if sum(label is not None for label in labels) > 1:
            ax.legend()
        if logo:
            add_logo(ax, style, location, description.info)
    return ax


def style_matplotlib_axes(ax, style):
    """Apply the style to axes that were created outside the style."""
    ax.set_facecolor(style.surface)
    ax.figure.set_facecolor(style.surface)
    for name, spine in ax.spines.items():
        spine.set_color(style.frame)
        spine.set_visible(name in ("left", "bottom") or style.frame == style.ink)
    ax.tick_params(colors=style.frame, labelcolor=style.ink_secondary)
    ax.xaxis.label.set_color(style.ink)
    ax.yaxis.label.set_color(style.ink)
    ax.title.set_color(style.ink)
    ax.grid(True, color=style.grid, linewidth=0.5)
    ax.set_axisbelow(True)


def _draw_series_matplotlib(ax, series, label, style, description):
    color = series_color(series, style)
    label = label or "_nolegend_"
    if series.kind == "bar":
        categories = description.y_axis.categories or list(series.y)
        positions = [categories.index(c) for c in series.y]
        ax.barh(
            positions, series.x, color=color, label=label, height=0.62,
            hatch=style.hatch(series.slot) or None, edgecolor=style.surface,
            linewidth=0
        )
        return
    if series.kind in ("heatmap", "contour"):
        cmap = style.colormap()
        if series.kind == "heatmap":
            mappable = ax.pcolormesh(
                series.x, series.y, series.z, cmap=cmap, shading="nearest",
                edgecolors=style.surface, linewidth=1.0
            )
        else:
            mappable = ax.contourf(series.x, series.y, series.z, levels=12, cmap=cmap)
            lines = ax.contour(
                series.x, series.y, series.z, levels=mappable.levels,
                colors=style.surface, linewidths=0.5
            )
            ax.clabel(lines, fontsize=style.font_size * 0.7, colors=style.ink)
        ax.grid(False)
        colorbar = ax.figure.colorbar(mappable, ax=ax, pad=0.02)
        colorbar.outline.set_visible(False)
        if description.z_axis is not None:
            colorbar.set_label(description.z_axis.title)
        return
    dash = DASHES[series.dash or style.dash(series.slot)]
    kwargs = {"color": color, "label": label, "linestyle": dash}
    if series.kind in ("markers", "line+markers"):
        kwargs.update({
            "marker": MARKERS[style.marker(series.slot)],
            "markeredgecolor": style.surface, "markeredgewidth": 0.8
        })
    if series.kind == "markers":
        kwargs["linestyle"] = "none"
    ax.plot(series.x, series.y, zorder=3, **kwargs)


def _set_axes_matplotlib(ax, description):
    for axis, set_label, set_scale, set_limits in (
            (description.x_axis, ax.set_xlabel, ax.set_xscale, ax.set_xlim),
            (description.y_axis, ax.set_ylabel, ax.set_yscale, ax.set_ylim)):
        set_label(axis.title)
        set_scale(axis.scale)
        if axis.range is not None:
            set_limits(*axis.range)
    categories = description.y_axis.categories
    if categories:
        ax.set_yticks(range(len(categories)), categories)
        ax.set_ylim(len(categories) - 0.5, -0.5)
        ax.grid(False, axis="y")


def _draw_annotation_matplotlib(ax, annotation, style):
    from matplotlib.patches import Rectangle
    from matplotlib.transforms import blended_transform_factory

    color = _annotation_color(annotation, style)
    small = style.font_size * 0.85
    x, y = annotation.x, annotation.y
    if annotation.kind == "segment":
        ax.plot(x, y, color=color, linewidth=0.8, zorder=4, label="_nolegend_")
        ax.annotate(
            annotation.text, ((x[0] + x[1]) / 2, (y[0] + y[1]) / 2),
            xytext=(6, 0), textcoords="offset points", va="center",
            color=color, fontsize=small
        )
    elif annotation.kind == "region":
        x_open, y_open = x[0] is None, y[0] is None
        transform = blended_transform_factory(
            ax.transAxes if x_open else ax.transData,
            ax.transAxes if y_open else ax.transData
        )
        x0, x1 = (0, 1) if x_open else x
        y0, y1 = (0, 1) if y_open else y
        ax.add_patch(Rectangle(
            (x0, y0), x1 - x0, y1 - y0, transform=transform, facecolor=color,
            alpha=0.12, linewidth=0, zorder=1
        ))
    elif annotation.kind == "vline":
        ax.axvline(x[0], color=color, linewidth=0.8, linestyle="--", zorder=2)
        ax.text(
            x[0], 0.99, f" {annotation.text}", va="top", ha="left",
            color=color, fontsize=small,
            transform=blended_transform_factory(ax.transData, ax.transAxes)
        )
    else:
        ax.annotate(
            annotation.text, (x[0], y[0]), xytext=(6, 6),
            textcoords="offset points", color=color, fontsize=small
        )


def _png_aspect(path):
    with open(path, "rb") as f:
        width, height = struct.unpack(">II", f.read(24)[16:24])
    return width / height


def plotly_layout(style, width=None, height=None):
    """Plotly layout settings of a style."""
    width, height = width or PLOTLY_SIZE[0], height or PLOTLY_SIZE[1]
    font_px = style.font_size * PX_PER_PT
    axis = {
        "gridcolor": style.grid, "linecolor": style.frame, "showline": True,
        "zeroline": False, "tickfont": {"color": style.ink_secondary},
        "title": {"font": {"color": style.ink}}, "ticks": "outside",
        "tickcolor": style.frame
    }
    return {
        "width": width, "height": height, "margin": dict(PLOTLY_MARGIN),
        "paper_bgcolor": style.surface, "plot_bgcolor": style.surface,
        "font": {
            "family": ", ".join(style.font_family), "color": style.ink,
            "size": font_px
        },
        "colorway": list(style.categorical),
        "xaxis": dict(axis), "yaxis": dict(axis),
        "hoverlabel": {
            "bgcolor": style.surface, "bordercolor": style.frame,
            "font": {"color": style.ink}
        },
        "legend": {
            "bgcolor": "rgba(0,0,0,0)", "font": {"color": style.ink_secondary}
        },
        "hovermode": "closest",
    }


def draw_plotly(description, style=None, logo=None, logo_location=None,
                width=None, height=None, **kwargs):
    """Draw a plot description with Plotly.

    Every point shows its values when hovered, the references of the points
    to network objects are the last entry of their :code:`customdata`. The
    whole description is stored in :code:`figure.layout.meta["tespy"]`.

    Parameters
    ----------
    description : tespy.plotting.description.PlotDescription
        The plot.

    style : str, Style, optional
        Style of the plot, the default style if not given.

    logo : bool, optional
        Show the logo and the information block.

    logo_location : str, optional
        Corner of the logo.

    width, height : int, optional
        Size of the figure in pixels.

    kwargs
        Passed on to the drawing of special plot kinds.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    if description.kind == "property_diagram":
        from tespy.plotting.cycles import draw_property_diagram_plotly
        return draw_property_diagram_plotly(
            description, style=style, logo=logo, logo_location=logo_location,
            width=width, height=height, **kwargs
        )

    import plotly.graph_objects as go

    style = get_style(style)
    logo, location = resolve(logo, logo_location)
    figure = go.Figure(layout=plotly_layout(style, width, height))
    labels = _legend_labels(description)
    for series, label in zip(description.series, labels):
        figure.add_trace(_plotly_trace(series, label, style, description))
    _set_axes_plotly(figure, description)
    for annotation in description.annotations:
        _draw_annotation_plotly(figure, annotation, style)
    if description.title:
        figure.update_layout(title={"text": f"<b>{description.title}</b>", "x": 0.02})
    figure.update_layout(
        showlegend=sum(label is not None for label in labels) > 1,
        meta={"tespy": description.to_dict()}
    )
    if logo:
        add_plotly_logo(figure, style, location, info_lines(description.info))
    return figure


def _hovertemplate(series, description):
    lines = [f"<b>{series.name}</b>"] if series.name else []
    if series.fields:
        for i, key in enumerate(series.fields):
            unit = series.units.get(key) or ""
            lines.append(f"{key} = %{{customdata[{i}]:.4g}} {unit}".rstrip())
    elif series.kind == "bar":
        lines.append(f"%{{y}}: %{{x:.4g}} {description.x_axis.unit or ''}".rstrip())
    elif series.kind in ("heatmap", "contour"):
        z_title = description.z_axis.title if description.z_axis else "value"
        lines += [
            f"{description.x_axis.title} = %{{x:.4g}}",
            f"{description.y_axis.title} = %{{y:.4g}}",
            f"{z_title} = %{{z:.4g}}",
        ]
    else:
        lines += [
            f"{description.x_axis.title} = %{{x:.4g}}",
            f"{description.y_axis.title} = %{{y:.4g}}",
        ]
    return "<br>".join(lines) + "<extra></extra>"


def _customdata(series):
    if not series.fields and series.references is None:
        return None
    columns = [list(values) for values in series.fields.values()]
    references = series.references or [series.reference] * len(series.x)
    return [
        [column[i] for column in columns] + [references[i]]
        for i in range(len(series.x))
    ]


def _plotly_trace(series, label, style, description):
    import plotly.graph_objects as go

    color = series_color(series, style)
    common = {
        "name": label or series.name, "showlegend": label is not None,
        "legendgroup": series.group or series.name,
        "meta": {"reference": series.reference, "role": series.role},
        "hovertemplate": _hovertemplate(series, description),
    }
    if series.kind == "bar":
        return go.Bar(
            x=series.x, y=series.y, orientation="h",
            marker={"color": color, "line": {"width": 0}},
            customdata=_customdata(series), **common
        )
    if series.kind in ("heatmap", "contour"):
        trace = go.Heatmap if series.kind == "heatmap" else go.Contour
        colorbar = {"outlinewidth": 0}
        if description.z_axis is not None:
            colorbar["title"] = {"text": description.z_axis.title}
        extra = {"xgap": 2, "ygap": 2} if series.kind == "heatmap" else {}
        return trace(
            x=series.x, y=series.y, z=series.z,
            colorscale=style.colorscale(), colorbar=colorbar, **common, **extra
        )
    mode = {"line": "lines", "markers": "markers", "line+markers": "lines+markers"}
    return go.Scatter(
        x=series.x, y=series.y, mode=mode[series.kind],
        line={"color": color, "width": 2.5, "dash": series.dash or style.dash(series.slot)},
        marker={
            "color": color, "size": 9, "symbol": style.marker(series.slot),
            "line": {"color": style.surface, "width": 1.5}
        },
        customdata=_customdata(series), **common
    )


def _set_axes_plotly(figure, description):
    import math

    for name, axis in (("xaxis", description.x_axis), ("yaxis", description.y_axis)):
        settings = {"title": {"text": axis.title}, "type": axis.scale}
        if axis.range is not None:
            settings["range"] = (
                [math.log10(v) for v in axis.range] if axis.scale == "log"
                else list(axis.range)
            )
        if axis.categories:
            settings.update({
                "type": "category", "categoryorder": "array",
                "categoryarray": list(axis.categories), "autorange": "reversed"
            })
        figure.update_layout({name: settings})


def _draw_annotation_plotly(figure, annotation, style):
    color = _annotation_color(annotation, style)
    font = {"color": color, "size": style.font_size * PX_PER_PT * 0.85}
    x, y = annotation.x, annotation.y
    if annotation.kind == "segment":
        figure.add_shape(
            type="line", x0=x[0], x1=x[1], y0=y[0], y1=y[1],
            line={"color": color, "width": 1}
        )
        figure.add_annotation(
            x=(x[0] + x[1]) / 2, y=(y[0] + y[1]) / 2, text=annotation.text,
            showarrow=False, xanchor="left", xshift=6, font=font
        )
    elif annotation.kind == "region":
        x_open, y_open = x[0] is None, y[0] is None
        figure.add_shape(
            type="rect", xref="paper" if x_open else "x",
            yref="paper" if y_open else "y",
            x0=0 if x_open else x[0], x1=1 if x_open else x[1],
            y0=0 if y_open else y[0], y1=1 if y_open else y[1],
            fillcolor=color, opacity=0.12, line={"width": 0}, layer="below"
        )
    elif annotation.kind == "vline":
        figure.add_shape(
            type="line", xref="x", yref="paper", x0=x[0], x1=x[0], y0=0,
            y1=1, line={"color": color, "width": 1, "dash": "dash"}
        )
        figure.add_annotation(
            x=x[0], y=1, xref="x", yref="paper", text=annotation.text,
            showarrow=False, xanchor="left", yanchor="top", font=font
        )
    else:
        figure.add_annotation(
            x=x[0], y=y[0], text=annotation.text, showarrow=False,
            xanchor="left", yanchor="bottom", xshift=6, yshift=6, font=font
        )


def add_plotly_logo(figure, style, location, lines, box=None):
    """Add the logo and the information lines to a Plotly figure.

    Parameters
    ----------
    figure : plotly.graph_objects.Figure
        The figure.

    style : Style
        Decides the logo variant and the text colors.

    location : str
        Corner of the plot area.

    lines : list
        Lines of text, drawn next to the logo.

    box : tuple, optional
        Paper coordinates :code:`(x0, y0, x1, y1)` of a box the logo goes
        next to instead of drawing the lines, e.g. the badge of CoolPlot.
    """
    width = figure.layout.width or PLOTLY_SIZE[0]
    height = figure.layout.height or PLOTLY_SIZE[1]
    margin = figure.layout.margin
    plot_width = width - (margin.l or 0) - (margin.r or 0)
    plot_height = height - (margin.t or 0) - (margin.b or 0)
    logo_px = LOGO_WIDTH_PT * PX_PER_PT
    size_x = logo_px / plot_width
    size_y = logo_px / _png_aspect(logo_path(style.logo, "png")) / plot_height
    upper = location.startswith("upper")
    right = location.endswith("right")
    gap = 4 / plot_height

    if box is None:
        font_px = style.font_size * PX_PER_PT * 0.75
        text_height = len(lines) * font_px * 1.3 / plot_height
        x = 0.99 if right else 0.01
        text_y = 0.99 - size_y - gap if upper else 0.01
        figure.add_annotation(
            x=x, y=text_y, xref="paper", yref="paper", showarrow=False,
            text="<br>".join(lines), align="right" if right else "left",
            xanchor="right" if right else "left",
            yanchor="top" if upper else "bottom",
            font={"color": style.ink_secondary, "size": font_px}
        )
        image_y = 0.99 if upper else 0.01 + text_height + gap
    else:
        x0, y0, x1, y1 = box
        x = x1 if right else x0
        image_y = y0 - gap if upper else y1 + gap
    figure.add_layout_image(
        source=logo_data_uri(style.logo), xref="paper", yref="paper", x=x,
        y=image_y, sizex=size_x, sizey=size_y, sizing="contain",
        xanchor="right" if right else "left",
        yanchor="top" if upper else "bottom", layer="above"
    )
