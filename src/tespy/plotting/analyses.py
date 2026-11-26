# -*- coding: utf-8

"""Plots of results: parameter sweeps and operating maps.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/plotting/analyses.py

SPDX-License-Identifier: MIT
"""
from tespy.plotting.description import Axis
from tespy.plotting.description import PlotDescription
from tespy.plotting.description import Series
from tespy.plotting.render import draw_matplotlib
from tespy.plotting.render import draw_plotly


def _columns(table):
    """Columns of a pandas DataFrame or a dict of lists as lists."""
    return {key: list(values) for key, values in dict(table).items()}


def sweep_data(table, x, y, group=None, units=None, title=None, info=None):
    """Collect the results of a parameter sweep as lines.

    Parameters
    ----------
    table : pandas.DataFrame, dict
        Results, one row per solve.

    x, y : str
        Columns on the horizontal and vertical axis.

    group : str, optional
        Column that splits the rows into lines, e.g. the working fluid.

    units : dict, optional
        Unit per column.

    title : str, optional
        Title of the plot.

    info : list, optional
        Further lines of the information block.

    Returns
    -------
    tespy.plotting.description.PlotDescription
    """
    columns = _columns(table)
    units = units or {}
    groups = [None] if group is None else list(dict.fromkeys(columns[group]))
    series = []
    for slot, value in enumerate(groups):
        rows = [
            i for i in range(len(columns[x]))
            if group is None or columns[group][i] == value
        ]
        series.append(Series(
            str(y if value is None else value),
            x=[columns[x][i] for i in rows], y=[columns[y][i] for i in rows],
            kind="line+markers", slot=slot,
            fields={x: [columns[x][i] for i in rows], y: [columns[y][i] for i in rows]},
            units={k: units[k] for k in (x, y) if k in units},
            reference=None if value is None else {group: value},
        ))
    return PlotDescription(
        kind="sweep",
        x_axis=Axis(x, units.get(x), quantity=x),
        y_axis=Axis(y, units.get(y), quantity=y),
        series=series, title=title or "", info=list(info or []),
    )


def operating_map_data(table, x, y, z, units=None, title=None, info=None):
    """Collect results on a grid of two parameters as a heat map.

    Parameters
    ----------
    table : pandas.DataFrame, dict
        Results, one row per solve.

    x, y : str
        Columns of the two parameters.

    z : str
        Column of the result, missing combinations stay empty.

    units : dict, optional
        Unit per column.

    title : str, optional
        Title of the plot.

    info : list, optional
        Further lines of the information block.

    Returns
    -------
    tespy.plotting.description.PlotDescription
    """
    columns = _columns(table)
    units = units or {}
    x_values = sorted(set(columns[x]))
    y_values = sorted(set(columns[y]))
    grid = [[float("nan")] * len(x_values) for _ in y_values]
    for xi, yi, zi in zip(columns[x], columns[y], columns[z]):
        grid[y_values.index(yi)][x_values.index(xi)] = zi
    return PlotDescription(
        kind="operating_map",
        x_axis=Axis(x, units.get(x), quantity=x),
        y_axis=Axis(y, units.get(y), quantity=y),
        z_axis=Axis(z, units.get(z), quantity=z),
        series=[Series(z, x=x_values, y=y_values, z=grid, kind="heatmap", role="scale")],
        title=title or "", info=list(info or []),
    )


def plot_sweep(table, x, y, group=None, units=None, ax=None, style=None,
               logo=None, logo_location=None, title=None, info=None):
    """Plot the results of a parameter sweep with matplotlib, see
    :py:func:`sweep_data`.

    Returns
    -------
    matplotlib.axes.Axes
    """
    description = sweep_data(table, x, y, group, units, title, info)
    return draw_matplotlib(
        description, ax=ax, style=style, logo=logo, logo_location=logo_location
    )


def plot_sweep_plotly(table, x, y, group=None, units=None, style=None,
                      logo=None, logo_location=None, title=None, info=None,
                      width=None, height=None):
    """Plot the results of a parameter sweep with Plotly, see
    :py:func:`sweep_data`.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    description = sweep_data(table, x, y, group, units, title, info)
    return draw_plotly(
        description, style=style, logo=logo, logo_location=logo_location,
        width=width, height=height
    )


def plot_operating_map(table, x, y, z, units=None, ax=None, style=None,
                       logo=None, logo_location=None, title=None, info=None):
    """Plot results on a grid of two parameters with matplotlib, see
    :py:func:`operating_map_data`.

    Returns
    -------
    matplotlib.axes.Axes
    """
    description = operating_map_data(table, x, y, z, units, title, info)
    return draw_matplotlib(
        description, ax=ax, style=style, logo=logo, logo_location=logo_location
    )


def plot_operating_map_plotly(table, x, y, z, units=None, style=None, logo=None,
                              logo_location=None, title=None, info=None,
                              width=None, height=None):
    """Plot results on a grid of two parameters with Plotly, see
    :py:func:`operating_map_data`.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    description = operating_map_data(table, x, y, z, units, title, info)
    return draw_plotly(
        description, style=style, logo=logo, logo_location=logo_location,
        width=width, height=height
    )
