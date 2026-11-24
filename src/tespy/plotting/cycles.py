# -*- coding: utf-8

"""Property diagrams of the fluid in a part of a network, drawn by CoolPlot.

The diagrams are available for fluids of the CoolProp back end. The state
points are the results of TESPy, CoolPlot draws the isolines, the saturation
dome and the paths between the states.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/plotting/cycles.py

SPDX-License-Identifier: MIT
"""
import math
import warnings

from tespy.plotting.branding import LOGO_WIDTH_PT
from tespy.plotting.branding import info_lines
from tespy.plotting.branding import logo_path
from tespy.plotting.description import Axis
from tespy.plotting.description import PlotDescription
from tespy.plotting.description import Series
from tespy.plotting.description import display_unit
from tespy.plotting.settings import resolve
from tespy.plotting.style import get_style
from tespy.tools.fluid_properties.helpers import single_fluid
from tespy.tools.fluid_properties.wrappers.coolprop import CoolPropWrapper

#: Diagram types, vertical axis first, and their (x, y) quantities.
DIAGRAM_TYPES = {
    "ph": ("h", "p"),
    "Ts": ("s", "T"),
    "hs": ("s", "h"),
    "ps": ("s", "p"),
    "prho": ("rho", "p"),
    "Trho": ("rho", "T"),
    "pT": ("T", "p"),
}

# label and network quantity of the state values
_QUANTITIES = {
    "T": ("temperature", "temperature"),
    "p": ("pressure", "pressure"),
    "h": ("specific enthalpy", "enthalpy"),
    "s": ("specific entropy", "entropy"),
    "v": ("specific volume", "specific_volume"),
    "m": ("mass flow", "mass_flow"),
}
_DENSITY_UNIT = "kg/m3"
_BASIS = {"mass": "mass", "molar": "mole"}
_CYCLE = "process"


def _branch(nw, connection_label):
    if connection_label not in nw.conns.index:
        msg = f"There is no connection with the label {connection_label}."
        raise KeyError(msg)
    connection = nw.get_conn(connection_label)
    for branch in nw.fluid_wrapper_branches.values():
        if connection in branch["connections"]:
            return branch["connections"], branch["components"]


def _fluid(connection):
    fluid = single_fluid(connection.fluid_data)
    if fluid is None:
        msg = (
            f"The fluid of connection {connection.label} is a mixture of "
            f"{', '.join(connection.fluid.val)}. Property diagrams are "
            "available for single CoolProp fluids only."
        )
        raise ValueError(msg)
    wrapper = connection.fluid.wrapper[fluid]
    if not isinstance(wrapper, CoolPropWrapper):
        msg = (
            f"The fluid {fluid} of connection {connection.label} uses "
            f"{type(wrapper).__name__}. Property diagrams are available for "
            "fluids of the CoolProp back end only."
        )
        raise ValueError(msg)
    components = wrapper.fluid.split("&")
    basis = None
    if wrapper.mixture_type is not None:
        if wrapper.mixture_type not in _BASIS:
            msg = (
                f"The fluid {fluid} has a {wrapper.mixture_type} based "
                "composition. Property diagrams need a mass or molar based "
                "composition."
            )
            raise ValueError(msg)
        basis = _BASIS[wrapper.mixture_type]
        components = [
            f"{name}[{fraction}]"
            for name, fraction in zip(components, wrapper.fractions)
        ]
    return f"{wrapper.back_end}::{'&'.join(components)}", basis


def _paths(connections, components):
    """Pairs of connection labels a component connects by a process."""
    paths = []
    for component in components:
        pairs = []
        if component.num_i == component.num_o:
            pairs = list(zip(component.inl, component.outl))
        elif component.num_i == 1:
            pairs = [(component.inl[0], o) for o in component.outl]
        elif component.num_o == 1:
            pairs = [(i, component.outl[0]) for i in component.inl]
        paths += [
            [a.label, b.label] for a, b in pairs
            if a in connections and b in connections
        ]
    return paths


def _state_value(connection, quantity):
    if quantity == "rho":
        return 1 / connection.vol.val_SI
    return connection.get_attr(quantity).val


def _axis_range(values, scale, padding=0.3):
    low, high = min(values), max(values)
    if scale == "log":
        low, high = math.log10(low), math.log10(high)
    delta = (high - low) or abs(high) * 0.1 or 1.0
    low, high = low - padding * delta, high + padding * delta
    if scale == "log":
        return [10 ** low, 10 ** high]
    return [low, high]


def property_diagram_data(nw, connection_label, diagram_type="ph", title=None, info=None):
    """Collect the data of a property diagram of a part of a network.

    The part is the fluid branch of the connection: all connections with the
    same fluid that are linked by components.

    Parameters
    ----------
    nw : tespy.networks.network.Network
        Solved network.

    connection_label : str
        Label of any connection of the branch.

    diagram_type : str
        Type of diagram, vertical axis first, see :py:data:`DIAGRAM_TYPES`.

    title : str, optional
        Title of the plot.

    info : list, optional
        Further lines of the information block.

    Returns
    -------
    tespy.plotting.description.PlotDescription
    """
    if diagram_type not in DIAGRAM_TYPES:
        msg = (
            f"The diagram type {diagram_type} is unknown, use one of "
            f"{', '.join(DIAGRAM_TYPES)}."
        )
        raise ValueError(msg)
    connections, components = _branch(nw, connection_label)
    connections = sorted(connections, key=lambda c: list(nw.conns.index).index(c.label))
    fluid, basis = _fluid(nw.get_conn(connection_label))
    units = nw.units.default

    fields = {
        key: [c.get_attr(key).val for c in connections] for key in _QUANTITIES
    }
    field_units = {
        key: display_unit(nw.units, units[quantity])
        for key, (_, quantity) in _QUANTITIES.items()
    }
    axes = []
    for quantity in DIAGRAM_TYPES[diagram_type]:
        scale = "log" if quantity == "p" else "linear"
        values = [_state_value(c, quantity) for c in connections]
        if quantity == "rho":
            label, unit = "density", _DENSITY_UNIT
        else:
            label, unit = _QUANTITIES[quantity][0], field_units[quantity]
        axes.append((Axis(
            label, unit, scale=scale, range=_axis_range(values, scale),
            quantity=quantity
        ), values))
    (x_axis, x), (y_axis, y) = axes

    states = Series(
        fluid, x=x, y=y, kind="markers", role="cycle", fields=fields,
        units=field_units,
        references=[{"connection": c.label} for c in connections],
    )
    return PlotDescription(
        kind="property_diagram", x_axis=x_axis, y_axis=y_axis,
        series=[states], title=title or "", info=list(info or []),
        extra={
            "fluid": fluid, "basis": basis, "diagram_type": diagram_type,
            "paths": _paths(connections, components),
            "states_SI": {
                c.label: {"p": c.p.val_SI, "h": c.h.val_SI} for c in connections
            },
            "units": {
                quantity: units[quantity]
                for quantity in [q for _, q in _QUANTITIES.values()]
                + ["temperature_difference"]
            },
        },
    )


def _unit_system(units):
    """CoolPlot unit system of the network units."""
    from coolplot.core.units import Unit
    from coolplot.core.units import get_unit_system

    from tespy.tools.units import Units

    converter = Units()
    system = get_unit_system("SI")
    for key, quantity in (("T", "temperature"), ("p", "pressure"),
                          ("h", "enthalpy"), ("s", "entropy"),
                          ("u", "enthalpy"), ("v", "specific_volume")):
        unit = units[quantity]
        si = converter.default[quantity]
        offset = converter.convert(0, si, unit)
        scale = converter.convert(1, si, unit) - offset
        difference = units["temperature_difference"] if key == "T" else unit
        system = system.with_unit(key, Unit(
            display_unit(converter, unit), scale, offset,
            display_unit(converter, difference)
        ))
    return system


def _theme(style):
    """CoolPlot theme of a style."""
    from dataclasses import replace

    from coolplot.style.theme import get_theme
    from coolplot.style.values import Font
    from coolplot.style.values import LineStyle
    from coolplot.style.values import MarkerStyle
    from coolplot.style.values import TextStyle

    def font(size, **kwargs):
        return Font(family=tuple(style.font_family), size_pt=size, **kwargs)

    base = get_theme("dark" if style.logo == "dark" else "default")
    dashes = {"solid": "solid", "dash": "dashed", "dashdot": "dashdot", "dot": "dotted"}
    if style.isolines is None:
        categorical, slot_dashes = (style.muted,), tuple(dashes[d] for d in style.dashes)
        isolines = None
    else:
        categorical, slot_dashes = style.categorical[:4], ("solid",)
        isolines = dict(style.isolines)
        isolines.setdefault("v", isolines.get("rho", style.muted))
    surface = style.surface
    return base.with_(
        name=style.name,
        categorical=categorical,
        isoline_colors=isolines,
        slot_dashes=slot_dashes,
        saturation=LineStyle(color=style.ink, width_pt=1.5, casing_color=surface, casing_width_pt=1.0),
        quality=LineStyle(color=style.muted, width_pt=0.5),
        isoline_label=TextStyle(
            font=font(7.5), color=style.ink_secondary, halo_color=surface,
            halo_width_pt=2.0
        ),
        processes=(LineStyle(
            color=style.cycle, width_pt=2.0, casing_color=surface,
            casing_width_pt=1.5
        ),),
        state_points=(MarkerStyle(
            shape="o", size_pt=6.0, face_color=style.cycle, edge_color=surface,
            edge_width_pt=1.0
        ),),
        state_number=TextStyle(
            font=font(8.0, weight="bold"), color=style.ink, halo_color=surface,
            halo_width_pt=2.0
        ),
        axes=replace(
            base.axes, background=surface, figure_background=surface,
            frame_color=style.frame, tick_color=style.frame,
            tick_label=TextStyle(font=font(8.0), color=style.ink_secondary),
            axis_label=TextStyle(font=font(9.5), color=style.ink),
            title=TextStyle(font=font(11.0, weight="bold"), color=style.ink),
            grid_major=LineStyle(color=style.grid, width_pt=0.5),
            legend=replace(base.axes.legend, visible=False),
        ),
    )


def _diagram(description, style, branding, diagram_hook):
    """CoolPlot diagram of a property diagram description."""
    from coolplot import PropertyDiagram
    from coolplot.spec.fluid import FluidSpec

    extra = description.extra
    fluid = FluidSpec.parse(extra["fluid"], basis=extra["basis"])
    diagram = PropertyDiagram(
        fluid, extra["diagram_type"], units=_unit_system(extra["units"]),
        theme=_theme(style)
    )
    diagram.settings = diagram.settings.with_(branding=branding)
    diagram.set_default_lines(labels=True, rounding=True)
    diagram.set_view(
        x_limits=description.x_axis.range, y_limits=description.y_axis.range
    )
    labels = [r["connection"] for r in description.series[0].references]
    numbers = {label: i + 1 for i, label in enumerate(labels)}
    diagram.add_cycle(
        _CYCLE,
        {numbers[label]: extra["states_SI"][label] for label in labels},
        loops=[[numbers[a], numbers[b]] for a, b in extra["paths"]],
        units="SI", numbers=False, labels={n: label for label, n in numbers.items()},
    )
    if description.title:
        diagram.title = description.title
    if diagram_hook is not None:
        diagram_hook(diagram)
    return diagram


def _branding(style, logo, location, lines, picture):
    from coolplot.spec.settings import BrandingSettings

    if not logo:
        return BrandingSettings(references="none", versions=False)
    pictures = (logo_path(style.logo, "svg"), logo_path(style.logo, "png"))
    return BrandingSettings(
        logo=pictures if picture else (), logo_width_pt=LOGO_WIDTH_PT,
        references="short", versions=True, lines=tuple(lines),
        location=location,
    )


def draw_property_diagram_matplotlib(description, ax=None, style=None, logo=None,
                                     logo_location=None, diagram_hook=None):
    """Draw a property diagram description with matplotlib, see
    :py:func:`tespy.plotting.render.draw_matplotlib`.

    Parameters
    ----------
    diagram_hook : callable, optional
        Called with the :code:`coolplot.PropertyDiagram` before it is drawn,
        e.g. to change isolines or limits.
    """
    from coolplot.render.mpl import MatplotlibRenderer
    from matplotlib import pyplot as plt

    style = get_style(style)
    logo, location = resolve(logo, logo_location)
    branding = _branding(style, logo, location, info_lines(description.info), True)
    # CoolPlot keeps the isoline labels clear of its own, centered title
    params = {**style.rc_params(), "axes.titlelocation": "center"}
    with plt.rc_context(params), warnings.catch_warnings():
        # CoolPlot reports isolines outside the view, which is expected
        warnings.filterwarnings("ignore", "Skipping the line", UserWarning)
        diagram = _diagram(description, style, branding, diagram_hook)
        if ax is None:
            _, ax = plt.subplots(figsize=(8, 6), layout="constrained")
        diagram.draw(MatplotlibRenderer(ax=ax))
    return ax


def draw_property_diagram_plotly(description, style=None, logo=None, logo_location=None,
                                 width=None, height=None, diagram_hook=None):
    """Draw a property diagram description with Plotly, see
    :py:func:`tespy.plotting.render.draw_plotly`.

    Parameters
    ----------
    diagram_hook : callable, optional
        Called with the :code:`coolplot.PropertyDiagram` before it is drawn.
    """
    from coolplot.render.plotly import PlotlyRenderer

    from tespy.plotting.render import _customdata
    from tespy.plotting.render import _hovertemplate
    from tespy.plotting.render import add_plotly_logo

    style = get_style(style)
    logo, location = resolve(logo, logo_location)
    # CoolPlot draws the text of the badge in Plotly, TESPy adds the logo
    branding = _branding(style, logo, location, info_lines(description.info), False)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "Skipping the line", UserWarning)
        diagram = _diagram(description, style, branding, diagram_hook)
        size = {"width": width or 800, "height": height or 600}
        figure = diagram.draw(PlotlyRenderer(**size)).figure

    states = description.series[0]
    for trace in figure.data:
        if (trace.meta or {}).get("item") == f"cycle/{_CYCLE}/states":
            trace.customdata = _customdata(states)
            trace.hovertemplate = _hovertemplate(states, description)
            trace.meta = {**trace.meta, "role": "cycle"}
        else:
            trace.hoverinfo = "skip"
    figure.update_layout(
        meta={"tespy": description.to_dict()},
        hoverlabel={
            "bgcolor": style.surface, "bordercolor": style.frame,
            "font": {"color": style.ink}
        },
    )
    if logo and figure.layout.shapes:
        box = figure.layout.shapes[0]
        add_plotly_logo(
            figure, style, location, [], box=(box.x0, box.y0, box.x1, box.y1)
        )
    return figure


def plot_property_diagram(nw, connection_label, diagram_type="ph", ax=None,
                          style=None, logo=None, logo_location=None,
                          title=None, info=None, diagram_hook=None):
    """Plot a property diagram of a part of a network with matplotlib.

    Parameters
    ----------
    nw : tespy.networks.network.Network
        Solved network.

    connection_label : str
        Label of any connection of the fluid branch to plot.

    diagram_type : str
        Type of diagram, vertical axis first, see :py:data:`DIAGRAM_TYPES`.

    ax : matplotlib.axes.Axes, optional
        Axes to draw into, a new figure if not given.

    style : str, Style, optional
        Style of the plot.

    logo : bool, optional
        Show the logo and the information block.

    logo_location : str, optional
        Corner of the logo and the information block.

    title : str, optional
        Title of the plot.

    info : list, optional
        Further lines of the information block.

    diagram_hook : callable, optional
        Called with the :code:`coolplot.PropertyDiagram` before it is drawn,
        e.g. to change isolines or limits.

    Returns
    -------
    matplotlib.axes.Axes
    """
    description = property_diagram_data(nw, connection_label, diagram_type, title, info)
    return draw_property_diagram_matplotlib(
        description, ax=ax, style=style, logo=logo,
        logo_location=logo_location, diagram_hook=diagram_hook
    )


def plot_property_diagram_plotly(nw, connection_label, diagram_type="ph",
                                 style=None, logo=None, logo_location=None,
                                 title=None, info=None, width=None,
                                 height=None, diagram_hook=None):
    """Plot a property diagram of a part of a network with Plotly.

    Hovering a state shows its values and the label of its connection. See
    :py:func:`plot_property_diagram` for the parameters.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    description = property_diagram_data(nw, connection_label, diagram_type, title, info)
    return draw_property_diagram_plotly(
        description, style=style, logo=logo, logo_location=logo_location,
        width=width, height=height, diagram_hook=diagram_hook
    )
