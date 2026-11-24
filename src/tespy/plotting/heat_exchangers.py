# -*- coding: utf-8

"""Q-T diagrams of heat exchangers.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/plotting/heat_exchangers.py

SPDX-License-Identifier: MIT
"""
import numpy as np

from tespy.components.heat_exchangers.base import HeatExchanger
from tespy.components.heat_exchangers.sectioned import SectionedHeatExchanger
from tespy.plotting.description import Annotation
from tespy.plotting.description import Axis
from tespy.plotting.description import PlotDescription
from tespy.plotting.description import Series
from tespy.plotting.description import display_unit
from tespy.plotting.render import draw_matplotlib
from tespy.plotting.render import draw_plotly


def _steps(component, points):
    """Section boundaries of the component and phase boundaries of both
    sides as fractions of the heat transferred, and a dense set of steps
    including them."""
    boundaries = [np.asarray(component._assign_steps(), dtype=float)]
    for side in range(2):
        steps, _ = SectionedHeatExchanger._get_moving_steps(
            component.inl[side], component.outl[side]
        )
        boundaries.append(np.asarray(steps, dtype=float))
    boundaries = np.unique(np.concatenate(boundaries))
    dense = np.unique(np.r_[np.linspace(0, 1, points), boundaries])
    return boundaries, dense


def qt_diagram_data(components, points=50, title=None, info=None):
    """Collect the data of a Q-T diagram of one or more heat exchangers.

    Several heat exchangers are drawn next to each other in the given order,
    e.g. the economizer, evaporator and superheater of a steam generator.

    Parameters
    ----------
    components : HeatExchanger, list
        Heat exchanger or list of heat exchangers of a solved network.

    points : int
        Number of points of each curve, the section and phase boundaries are
        added.

    title : str, optional
        Title of the plot.

    info : list, optional
        Further lines of the information block.

    Returns
    -------
    tespy.plotting.description.PlotDescription
    """
    if not isinstance(components, (list, tuple)):
        components = [components]
    for component in components:
        if not isinstance(component, HeatExchanger):
            msg = (
                f"{component.label} is a {type(component).__name__}, Q-T "
                "diagrams are available for components of the class "
                "HeatExchanger and its subclasses."
            )
            raise TypeError(msg)

    units = components[0].Q_sections._units
    heat_unit = units.default["heat"]
    T_unit = units.default["temperature"]
    dT_unit = units.default["temperature_difference"]

    def heat(values):
        return units.convert(np.asarray(values), "W", heat_unit)

    def temperature(values):
        return units.convert(np.asarray(values), "K", T_unit)

    def difference(values):
        return units.convert(np.asarray(values), "delta_degC", dT_unit)

    labels = {
        "Q": display_unit(units, heat_unit),
        "T": display_unit(units, T_unit),
        "dT": display_unit(units, dT_unit),
    }
    series, annotations = [], []
    offset = 0.0
    pinch = None
    for component in components:
        boundaries, dense = _steps(component, points)
        for steps, kind in ((dense, "line"), (boundaries, "markers")):
            T_hot, T_cold = component._get_T_at_steps(steps)
            Q = offset + heat(component._get_Q_cumsum_steps(steps))
            dT = difference(T_hot - T_cold)
            fields = {
                "Q": Q, "T hot": temperature(T_hot),
                "T cold": temperature(T_cold), "ΔT": dT
            }
            field_units = {
                "Q": labels["Q"], "T hot": labels["T"],
                "T cold": labels["T"], "ΔT": labels["dT"]
            }
            for side, T, slot in (("hot", T_hot, 0), ("cold", T_cold, 1)):
                series.append(Series(
                    f"{component.label} {side} side", x=Q, y=temperature(T),
                    kind=kind, role=side, slot=slot, fields=fields,
                    units=field_units, group=f"{side} side",
                    reference={"component": component.label, "side": side},
                    show_legend=kind == "line",
                ))
            if kind == "line":
                i = int(np.argmin(dT))
                if pinch is None or dT[i] < pinch[0]:
                    pinch = (dT[i], Q[i], temperature(T_cold[i]), temperature(T_hot[i]))
        if offset > 0:
            annotations.append(Annotation(
                "vline", x=[offset], y=[None], text=component.label
            ))
        offset = series[-1].x[-1]

    dT, Q, T_cold, T_hot = pinch
    annotations.append(Annotation(
        "segment", x=[Q, Q], y=[T_cold, T_hot],
        text=f"pinch {dT:.1f} {labels['dT']}"
    ))
    default_title = "Q-T diagram of " + ", ".join(c.label for c in components)
    return PlotDescription(
        kind="qt_diagram",
        x_axis=Axis("heat transferred", labels["Q"], quantity="Q"),
        y_axis=Axis("temperature", labels["T"], quantity="T"),
        series=series, annotations=annotations,
        title=default_title if title is None else title,
        info=list(info or []),
    )


def plot_QT_diagram(components, ax=None, style=None, logo=None, logo_location=None,
                    title=None, info=None, points=50):
    """Plot the Q-T diagram of one or more heat exchangers with matplotlib.

    See :py:func:`qt_diagram_data` for the data and
    :py:func:`tespy.plotting.render.draw_matplotlib` for the drawing
    parameters.

    Returns
    -------
    matplotlib.axes.Axes
    """
    description = qt_diagram_data(components, points, title, info)
    return draw_matplotlib(
        description, ax=ax, style=style, logo=logo, logo_location=logo_location
    )


def plot_QT_diagram_plotly(components, style=None, logo=None, logo_location=None,
                           title=None, info=None, points=50, width=None,
                           height=None):
    """Plot the Q-T diagram of one or more heat exchangers with Plotly.

    See :py:func:`qt_diagram_data` for the data and
    :py:func:`tespy.plotting.render.draw_plotly` for the drawing parameters.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    description = qt_diagram_data(components, points, title, info)
    return draw_plotly(
        description, style=style, logo=logo, logo_location=logo_location,
        width=width, height=height
    )
