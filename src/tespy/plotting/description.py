# -*- coding: utf-8

"""Plots as plain data.

A :py:class:`PlotDescription` holds everything a plot shows: axes with
units, series with the values behind every point, annotations and the lines
of the information block. Every point can carry a reference to the network
object it belongs to (connection or component label). The matplotlib and
the Plotly renderer draw the same description, and the description can be
exported to JSON.

The format is provisional, its version is :py:data:`FORMAT_VERSION`.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/plotting/description.py

SPDX-License-Identifier: MIT
"""
import json
import math
from dataclasses import asdict
from dataclasses import dataclass
from dataclasses import field

FORMAT = "tespy-plot"
FORMAT_VERSION = 0

SERIES_KINDS = (
    "line", "markers", "line+markers", "bar", "heatmap", "contour", "cells"
)
SERIES_ROLES = (
    "category", "hot", "cold", "cycle", "reference", "invalid", "scale"
)
ANNOTATION_KINDS = ("segment", "region", "box", "vline", "text")
AXIS_SCALES = ("linear", "log")


def display_unit(units, unit):
    """Unit as shown in labels, e.g. °C for degC.

    Units with a division keep the spelling of the network, because pint
    reorders their parts (kJ/kg/K would read kJ/K/kg).

    Parameters
    ----------
    units : tespy.tools.units.Units
        Unit registry of the network.

    unit : str
        Unit string.
    """
    pretty = f"{units.parsed_unit(unit):~P}"
    return unit if "/" in pretty else pretty


def _plain(value):
    """Python floats and lists instead of numpy types, None for NaN and
    infinite values."""
    if isinstance(value, dict):
        return {k: _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)) or hasattr(value, "tolist"):
        value = value.tolist() if hasattr(value, "tolist") else value
        if not isinstance(value, list):
            return _plain(value)
        return [_plain(v) for v in value]
    if isinstance(value, bool) or value is None or isinstance(value, str):
        return value
    if isinstance(value, (int, float)):
        value = float(value)
        return value if math.isfinite(value) else None
    return value


@dataclass
class Axis:
    """An axis of a plot.

    Parameters
    ----------
    label : str
        Name of the quantity.

    unit : str, optional
        Unit of the values.

    scale : str
        :code:`"linear"` or :code:`"log"`.

    range : list, optional
        Lower and upper limit, chosen by the renderer if not given.

    categories : list, optional
        Category names in display order for bar plots.

    quantity : str, optional
        Short name of the quantity, e.g. :code:`"T"`.

    tick_labels : list, optional
        Labels of the positions 0, 1, 2, ... of the axis, e.g. the
        equations of a matrix plot.
    """

    label: str
    unit: str = None
    scale: str = "linear"
    range: list = None
    categories: list = None
    quantity: str = None
    tick_labels: list = None

    def __post_init__(self):
        if self.scale not in AXIS_SCALES:
            msg = f"The axis scale is one of {', '.join(AXIS_SCALES)}, not {self.scale}."
            raise ValueError(msg)

    @property
    def title(self):
        if self.unit:
            return f"{self.label} in {self.unit}"
        return self.label


@dataclass
class Series:
    """A series of a plot.

    Parameters
    ----------
    name : str
        Name in the legend and the hover box.

    x, y : list
        Values of the points, category names on the category axis of bar
        plots.

    kind : str
        How the series is drawn, one of :py:data:`SERIES_KINDS`.

    role : str
        What the series stands for, decides its color, one of
        :py:data:`SERIES_ROLES`.

    slot : int
        Categorical slot for the role :code:`"category"`, also selects the
        marker and dash.

    z : list, optional
        Values on the grid of :code:`x` and :code:`y` (one row per y value)
        for heat maps and contour plots. For the kind :code:`"cells"` the
        kind of every cell of a matrix, see
        :py:data:`tespy.plotting.render.cell_colors`.

    fields : dict, optional
        Further values per point, shown when hovering, e.g. :code:`{"T":
        [...]}`.

    units : dict, optional
        Units of the fields.

    references : list, optional
        Per point the network object it belongs to, e.g.
        :code:`{"connection": "4"}`.

    reference : dict, optional
        Network object of the whole series.

    dash : str, optional
        Dash of the line, the dash of the slot if not given.

    group : str, optional
        Series of a group share one legend entry.

    color_value : float, optional
        Position on the sequential colormap for the role :code:`"scale"`,
        between 0 and 1.

    show_legend : bool
        List the series in the legend.
    """

    name: str
    x: list
    y: list
    kind: str = "line"
    role: str = "category"
    slot: int = 0
    z: list = None
    fields: dict = field(default_factory=dict)
    units: dict = field(default_factory=dict)
    references: list = None
    reference: dict = None
    dash: str = None
    group: str = None
    color_value: float = None
    show_legend: bool = True

    def __post_init__(self):
        if self.kind not in SERIES_KINDS:
            msg = f"The series kind is one of {', '.join(SERIES_KINDS)}, not {self.kind}."
            raise ValueError(msg)
        if self.role not in SERIES_ROLES:
            msg = f"The series role is one of {', '.join(SERIES_ROLES)}, not {self.role}."
            raise ValueError(msg)
        for key, values in self.fields.items():
            if len(values) != len(self.x):
                msg = (
                    f"The field {key} of the series {self.name} needs one "
                    f"value per point ({len(self.x)})."
                )
                raise ValueError(msg)


@dataclass
class Annotation:
    """A mark that is not data: a pinch, a boundary or an invalid range.

    Parameters
    ----------
    kind : str
        :code:`"segment"` (a line between two points), :code:`"region"` (a
        filled rectangle, :code:`None` for an open side), :code:`"box"` (a
        dashed rounded outline), :code:`"vline"` (a vertical line) or
        :code:`"text"`.

    x, y : list
        Coordinates.

    text : str
        Label of the mark.

    role : str
        Decides the color, see :py:data:`SERIES_ROLES`.
    """

    kind: str
    x: list
    y: list
    text: str = ""
    role: str = "reference"

    def __post_init__(self):
        if self.kind not in ANNOTATION_KINDS:
            msg = f"The annotation kind is one of {', '.join(ANNOTATION_KINDS)}, not {self.kind}."
            raise ValueError(msg)


@dataclass
class PlotDescription:
    """Everything a plot shows, independent of the plotting library.

    Parameters
    ----------
    kind : str
        Kind of the plot, e.g. :code:`"qt_diagram"`.

    x_axis, y_axis : Axis
        Axes of the plot.

    series : list
        :py:class:`Series` of the plot.

    annotations : list
        :py:class:`Annotation` of the plot.

    title : str
        Title of the plot.

    z_axis : Axis, optional
        Color axis of heat maps and contour plots.

    info : list
        Further lines of the information block.

    extra : dict
        Data specific to the kind of plot, e.g. the fluid of a property
        diagram.
    """

    kind: str
    x_axis: Axis
    y_axis: Axis
    series: list
    annotations: list = field(default_factory=list)
    title: str = ""
    z_axis: Axis = None
    info: list = field(default_factory=list)
    extra: dict = field(default_factory=dict)

    def to_dict(self):
        """The description as plain data with a format header."""
        data = _plain(asdict(self))
        return {"meta": {"format": FORMAT, "version": FORMAT_VERSION}, **data}

    def to_json(self, **kwargs):
        """The description as JSON text, keywords go to :code:`json.dumps`."""
        return json.dumps(self.to_dict(), allow_nan=False, **kwargs)

    @classmethod
    def from_dict(cls, data):
        """Create a description from :py:meth:`to_dict` output."""
        data = dict(data)
        meta = data.pop("meta", {})
        if meta.get("format") != FORMAT:
            msg = f"The data are not a {FORMAT} description."
            raise ValueError(msg)
        if meta.get("version", 0) > FORMAT_VERSION:
            msg = (
                f"The description has format version {meta['version']}, this "
                f"version of TESPy reads up to version {FORMAT_VERSION}."
            )
            raise ValueError(msg)
        for key in ("x_axis", "y_axis", "z_axis"):
            if data.get(key) is not None:
                data[key] = Axis(**data[key])
        data["series"] = [Series(**s) for s in data["series"]]
        data["annotations"] = [Annotation(**a) for a in data.get("annotations", [])]
        return cls(**data)

    @classmethod
    def from_json(cls, text):
        """Create a description from :py:meth:`to_json` output."""
        return cls.from_dict(json.loads(text))
