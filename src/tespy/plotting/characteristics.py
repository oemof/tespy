# -*- coding: utf-8

"""Plots of characteristic lines, maps and polynomials.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/plotting/characteristics.py

SPDX-License-Identifier: MIT
"""
import numpy as np

from tespy.components.component import Component
from tespy.plotting.description import Annotation
from tespy.plotting.description import Axis
from tespy.plotting.description import PlotDescription
from tespy.plotting.description import Series
from tespy.plotting.description import display_unit
from tespy.plotting.render import draw_matplotlib
from tespy.plotting.render import draw_plotly
from tespy.tools.characteristics import CharLine
from tespy.tools.characteristics import CharMap


def _padded(low, high, padding=0.08):
    delta = (high - low) or abs(high) * 0.1 or 1.0
    return [low - padding * delta, high + padding * delta]


def _source(characteristic, key):
    """Characteristic function, axis names, axis units, unit registry,
    operating point and references of a characteristic or a component
    characteristic."""
    if not isinstance(characteristic, Component):
        from tespy.tools.units import Units

        char = characteristic
        names = list(char.axes) if char.axes is not None else (
            [f"x{i + 1}" for i in range(char.num_inputs)] + ["y"]
        )
        units = [getattr(char, f"{name}_unit", None) for name in names]
        return char, names, units, Units(), [], None, ""

    component = characteristic
    if key is None:
        options = ", ".join(component.characteristics_in_use()) or "e.g. eta_s_char"
        msg = (
            f"Name the characteristic of {component.label} to plot, "
            f"{options}."
        )
        raise ValueError(msg)
    data = component.get_attr(key)
    char = data.char_func
    names = list(data.axes) if data.axes is not None else (
        [f"x{i + 1}" for i in range(char.num_inputs)] + ["y"]
    )
    _, output = data._axis_units()
    units = [data.input_unit(i) for i in range(len(names) - 1)]
    units.append(None if output is None else output[1])
    inputs = component.get_char_inputs(key)
    points = [] if inputs is None else [tuple(data._inputs_to_data(inputs))]
    reference = {"component": component.label, "characteristic": key}
    return char, names, units, data._units, points, reference, f"{component.label}: {key}"


def _point_series(char, points, names, units, reference, x_index, y_of):
    """Operating points inside and outside the range of the data."""
    series = []
    groups = {"cycle": [], "invalid": []}
    for point in points:
        role = "invalid" if char.domain_violations(*point) else "cycle"
        groups[role].append(point)
    for role, chosen in groups.items():
        if not chosen:
            continue
        outputs = [char.evaluate(*p) for p in chosen]
        fields = {
            name: [p[i] for p in chosen] for i, name in enumerate(names[:-1])
        }
        fields[names[-1]] = outputs
        name = "operating point" if role == "cycle" else "outside the data"
        series.append(Series(
            name, x=[p[x_index] for p in chosen], y=y_of(chosen, outputs),
            kind="markers", role=role, slot=1, fields=fields,
            units={n: u for n, u in zip(names, units) if u},
            references=[reference] * len(chosen) if reference else None,
        ))
    return series


def characteristic_data(characteristic, key=None, operating_points=None,
                        ranges=None, points=50, title=None, info=None):
    """Collect the data of a plot of a characteristic.

    Parameters
    ----------
    characteristic : CharLine, CharMap, Polynomial, Component
        The characteristic, or a component of a solved network together
        with the name of its characteristic. The operating point of the
        component is marked.

    key : str, optional
        Name of the characteristic of the component, e.g.
        :code:`"eta_s_char"`.

    operating_points : list, optional
        Further operating points, one tuple of inputs per point in the
        units of the data.

    ranges : list, optional
        Lower and upper value per input for polynomials, which have no
        range of their own.

    points : int
        Number of points per input of polynomials.

    title : str, optional
        Title of the plot.

    info : list, optional
        Further lines of the information block.

    Returns
    -------
    tespy.plotting.description.PlotDescription
    """
    char, names, units, registry, own_points, reference, default_title = (
        _source(characteristic, key)
    )
    all_points = own_points + [tuple(p) for p in operating_points or []]
    labels = [display_unit(registry, u) if u else None for u in units]

    if isinstance(char, CharLine):
        series, annotations, axes = _char_line(char, all_points, names, labels, reference)
    elif isinstance(char, CharMap):
        series, annotations, axes = _char_map(char, all_points, names, labels, reference)
    else:
        series, annotations, axes = _polynomial(
            char, all_points, names, labels, reference, ranges, points
        )
    x_axis, y_axis, z_axis = axes
    return PlotDescription(
        kind="characteristic", x_axis=x_axis, y_axis=y_axis, z_axis=z_axis,
        series=series, annotations=annotations,
        title=default_title if title is None else title,
        info=list(info or []),
    )


def _char_line(char, points, names, labels, reference):
    x_values = [char.x[0], char.x[-1]] + [p[0] for p in points]
    x_range = _padded(min(x_values), max(x_values))
    series = [Series(
        "data", x=char.x, y=char.y, kind="line+markers",
        fields={names[0]: char.x, names[1]: char.y},
        units={n: u for n, u in zip(names, labels) if u},
    )]
    for point in points:
        x = point[0]
        if char.extrapolate and not char.x[0] <= x <= char.x[-1]:
            start = char.x[0] if x < char.x[0] else char.x[-1]
            series.append(Series(
                "extrapolated", x=[start, x],
                y=[char.evaluate(start), char.evaluate(x)], role="reference",
                dash="dash", show_legend=False,
            ))
    series += _point_series(
        char, points, names, labels, reference, 0, lambda _, outputs: outputs
    )
    annotations = [
        Annotation("region", x=[x_range[0], char.x[0]], y=[None, None], role="invalid"),
        Annotation("region", x=[char.x[-1], x_range[1]], y=[None, None], role="invalid"),
    ]
    axes = (
        Axis(names[0], labels[0], range=x_range, quantity=names[0]),
        Axis(names[1], labels[1], quantity=names[1]),
        None,
    )
    return series, annotations, axes


def _char_map(char, points, names, labels, reference):
    series = []
    n = len(char.x)
    for i, x in enumerate(char.x):
        series.append(Series(
            f"{names[0]} = {x:.4g}", x=char.y[i], y=char.z[i], kind="line+markers",
            role="scale", color_value=0.35 + 0.65 * i / max(n - 1, 1),
            fields={names[1]: char.y[i], names[2]: char.z[i]},
            units={k: u for k, u in zip(names[1:], labels[1:]) if u},
        ))
    for point in points:
        yarr, zarr = char.evaluate_x(point[0])
        series.append(Series(
            f"{names[0]} = {point[0]:.4g}", x=yarr, y=zarr, role="reference",
            dash="dash", show_legend=False,
        ))
    series += _point_series(
        char, points, names, labels, reference, 1, lambda _, outputs: outputs
    )
    axes = (
        Axis(names[1], labels[1], quantity=names[1]),
        Axis(names[2], labels[2], quantity=names[2]),
        None,
    )
    return series, [], axes


def _polynomial(char, points, names, labels, reference, ranges, n):
    if char.num_inputs > 2:
        msg = (
            f"Plots are available for polynomials with one or two inputs, "
            f"this one has {char.num_inputs}."
        )
        raise ValueError(msg)
    if ranges is None or len(ranges) != char.num_inputs:
        msg = (
            "A polynomial has no range of its own, pass ranges with a lower "
            "and upper value per input, e.g. ranges=[(-20, 10), (30, 60)]."
        )
        raise ValueError(msg)
    grids = [np.linspace(low, high, n) for low, high in ranges]
    units = {k: u for k, u in zip(names, labels) if u}
    if char.num_inputs == 1:
        x = grids[0]
        series = [Series(
            "polynomial", x=x, y=[char.evaluate(v) for v in x],
            fields={names[0]: x, names[1]: [char.evaluate(v) for v in x]},
            units=units,
        )]
        series += _point_series(
            char, points, names, labels, reference, 0, lambda _, outputs: outputs
        )
        axes = (
            Axis(names[0], labels[0], quantity=names[0]),
            Axis(names[1], labels[1], quantity=names[1]),
            None,
        )
        return series, [], axes

    x, y = grids
    z = [[char.evaluate(xi, yi) for xi in x] for yi in y]
    series = [Series("polynomial", x=x, y=y, z=z, kind="contour", role="scale")]
    series += _point_series(
        char, points, names, labels, reference, 0,
        lambda chosen, _: [p[1] for p in chosen]
    )
    axes = (
        Axis(names[0], labels[0], quantity=names[0]),
        Axis(names[1], labels[1], quantity=names[1]),
        Axis(names[2], labels[2], quantity=names[2]),
    )
    return series, [], axes


def plot_characteristic(characteristic, key=None, ax=None, style=None, logo=None,
                        logo_location=None, operating_points=None, ranges=None,
                        points=50, title=None, info=None):
    """Plot a characteristic with matplotlib.

    See :py:func:`characteristic_data` for the data and
    :py:func:`tespy.plotting.render.draw_matplotlib` for the drawing
    parameters.

    Returns
    -------
    matplotlib.axes.Axes
    """
    description = characteristic_data(
        characteristic, key, operating_points, ranges, points, title, info
    )
    return draw_matplotlib(
        description, ax=ax, style=style, logo=logo, logo_location=logo_location
    )


def plot_characteristic_plotly(characteristic, key=None, style=None, logo=None,
                               logo_location=None, operating_points=None,
                               ranges=None, points=50, title=None, info=None,
                               width=None, height=None):
    """Plot a characteristic with Plotly.

    See :py:func:`characteristic_data` for the data and
    :py:func:`tespy.plotting.render.draw_plotly` for the drawing parameters.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    description = characteristic_data(
        characteristic, key, operating_points, ranges, points, title, info
    )
    return draw_plotly(
        description, style=style, logo=logo, logo_location=logo_location,
        width=width, height=height
    )
