# -*- coding: utf-8

"""Styles for all plots of TESPy.

A :py:class:`Style` holds colors and line styles by their role in a plot
(surface, ink, hot side, cold side, categories, ...). The matplotlib
settings, the colormaps, the Plotly layout and the CoolPlot theme of a plot
are all derived from the same style, so that plots of different kinds look
alike.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/plotting/style.py

SPDX-License-Identifier: MIT
"""
import dataclasses
import tomllib
from dataclasses import dataclass
from dataclasses import field

#: Marker names of a style and their matplotlib symbols, the names are the
#: ones Plotly uses.
MARKERS = {
    "circle": "o",
    "square": "s",
    "diamond": "D",
    "triangle-up": "^",
    "triangle-down": "v",
    "cross": "P",
    "x": "X",
    "pentagon": "p",
}

#: Dash names of a style and their matplotlib line styles, the names are the
#: ones Plotly uses.
DASHES = {
    "solid": "-",
    "dash": "--",
    "dashdot": "-.",
    "dot": ":",
}

_COLOR_ROLES = (
    "surface", "ink", "ink_secondary", "muted", "grid", "frame", "hot",
    "cold", "cycle", "invalid"
)
_COLOR_LISTS = ("categorical", "sequential", "diverging")


def _to_hex(color):
    if isinstance(color, str) and color.startswith("#"):
        digits = color[1:].lower()
        if len(digits) == 3:
            digits = "".join(d * 2 for d in digits)
        if len(digits) == 6 and all(d in "0123456789abcdef" for d in digits):
            return f"#{digits}"
    from matplotlib.colors import to_hex
    return to_hex(color)


@dataclass(frozen=True)
class Style:
    r"""
    Colors, markers, line styles and fonts of a plot, by their role.

    Parameters
    ----------
    name : str
        Name the style is registered under.

    surface, ink, ink_secondary, muted, grid, frame : str
        Background, main text, secondary text (ticks, legend), recessive
        marks, grid lines and the axes frame.

    categorical : tuple
        Colors of categories (fluids, components, series) in a fixed order.

    hot, cold : str
        Hot and cold side of heat exchangers.

    cycle : str
        Process lines and state points in property diagrams.

    invalid : str
        Marks and areas outside the valid range of data.

    sequential : tuple
        Colors of the sequential colormap, from the surface away.

    diverging : tuple
        Colors of the diverging colormap with a neutral midpoint.

    isolines : dict
        Isoline color per quantity (:code:`"T"`, :code:`"p"`, :code:`"h"`,
        :code:`"s"`, :code:`"rho"`, :code:`"u"`) in property diagrams,
        :code:`None` to use the categorical colors with dashes.

    markers : tuple
        Marker per category slot, see :py:data:`MARKERS`.

    dashes : tuple
        Dash per category slot, see :py:data:`DASHES`.

    hatches : tuple
        matplotlib hatch per category slot for filled areas.

    font_family : tuple
        Fonts in order of preference.

    font_size : float
        Base font size in points.

    logo : str
        Variant of the TESPy logo, :code:`"light"` or :code:`"dark"`.

    Example
    -------
    >>> from tespy.plotting.style import get_style
    >>> company = get_style("tespy").derive("company", cycle="#00395B")
    >>> company.cycle
    '#00395b'
    >>> company.surface == get_style("tespy").surface
    True
    """

    name: str
    surface: str
    ink: str
    ink_secondary: str
    muted: str
    grid: str
    frame: str
    categorical: tuple
    hot: str
    cold: str
    cycle: str
    invalid: str
    sequential: tuple
    diverging: tuple
    isolines: dict = field(default=None, hash=False)
    markers: tuple = (
        "circle", "square", "diamond", "triangle-up", "triangle-down", "cross"
    )
    dashes: tuple = ("solid",)
    hatches: tuple = ("",)
    font_family: tuple = ("DejaVu Sans", "Arial", "Helvetica", "sans-serif")
    font_size: float = 9.0
    logo: str = "light"

    def __post_init__(self):
        for role in _COLOR_ROLES:
            object.__setattr__(self, role, _to_hex(getattr(self, role)))
        for role in _COLOR_LISTS:
            colors = tuple(_to_hex(c) for c in getattr(self, role))
            object.__setattr__(self, role, colors)
        if not self.categorical:
            msg = f"Style {self.name}: give at least one categorical color."
            raise ValueError(msg)
        for role in ("sequential", "diverging"):
            if len(getattr(self, role)) < 2:
                msg = f"Style {self.name}: the {role} colormap needs two colors."
                raise ValueError(msg)
        if self.isolines is not None:
            isolines = {k: _to_hex(v) for k, v in dict(self.isolines).items()}
            object.__setattr__(self, "isolines", isolines)
        for role, options in (("markers", MARKERS), ("dashes", DASHES)):
            values = tuple(getattr(self, role))
            unknown = [v for v in values if v not in options]
            if not values or unknown:
                msg = (
                    f"Style {self.name}: unknown {role[:-1]} {unknown}, use "
                    f"{', '.join(options)}."
                )
                raise ValueError(msg)
            object.__setattr__(self, role, values)
        for role in ("hatches", "font_family"):
            object.__setattr__(self, role, tuple(getattr(self, role)) or ("",))
        if self.logo not in ("light", "dark"):
            msg = f"Style {self.name}: logo is 'light' or 'dark'."
            raise ValueError(msg)

    def derive(self, name, **changes):
        """Copy of the style under a new name with some roles changed."""
        return dataclasses.replace(self, name=name, **changes)

    def color(self, slot):
        """Categorical color of a slot, the colors repeat after the last."""
        return self.categorical[slot % len(self.categorical)]

    def marker(self, slot):
        """Marker name of a slot."""
        return self.markers[slot % len(self.markers)]

    def dash(self, slot):
        """Dash name of a slot."""
        return self.dashes[slot % len(self.dashes)]

    def hatch(self, slot):
        """matplotlib hatch of a slot."""
        return self.hatches[slot % len(self.hatches)]

    def to_dict(self):
        """The style as plain data, e.g. for a TOML or JSON file."""
        data = dataclasses.asdict(self)
        for key, value in data.items():
            if isinstance(value, tuple):
                data[key] = list(value)
        return data

    @classmethod
    def from_dict(cls, data):
        """Create a style from plain data.

        With the key :code:`based_on` the data only need the roles that
        differ from that registered style.
        """
        data = dict(data)
        based_on = data.pop("based_on", None)
        if based_on is not None:
            name = data.pop("name", based_on)
            return get_style(based_on).derive(name, **data)
        return cls(**data)

    @classmethod
    def from_toml(cls, path):
        """Create a style from a TOML file, see :py:meth:`from_dict`."""
        with open(path, "rb") as f:
            return cls.from_dict(tomllib.load(f))

    def colormap(self, kind="sequential"):
        """matplotlib colormap of the style.

        Parameters
        ----------
        kind : str
            :code:`"sequential"` or :code:`"diverging"`.
        """
        from matplotlib.colors import LinearSegmentedColormap

        colors = self._colormap_colors(kind)
        name = self.name if kind == "sequential" else f"{self.name}_{kind}"
        return LinearSegmentedColormap.from_list(name, colors)

    def colorscale(self, kind="sequential"):
        """Plotly colorscale of the style, see :py:meth:`colormap`."""
        colors = self._colormap_colors(kind)
        n = len(colors) - 1
        return [[i / n, color] for i, color in enumerate(colors)]

    def _colormap_colors(self, kind):
        if kind not in ("sequential", "diverging"):
            msg = f"The colormap kind is 'sequential' or 'diverging', not {kind}."
            raise ValueError(msg)
        return getattr(self, kind)

    def rc_params(self):
        """matplotlib settings of the style."""
        return {
            "font.family": list(self.font_family),
            "font.size": self.font_size,
            "axes.titlesize": self.font_size * 1.15,
            "axes.titleweight": "bold",
            "axes.titlelocation": "left",
            "axes.labelsize": self.font_size,
            "axes.labelcolor": self.ink,
            "axes.edgecolor": self.frame,
            "axes.linewidth": 0.8,
            "axes.facecolor": self.surface,
            "axes.grid": True,
            "axes.axisbelow": True,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": self.surface,
            "savefig.facecolor": self.surface,
            "text.color": self.ink,
            "xtick.color": self.frame,
            "ytick.color": self.frame,
            "xtick.labelcolor": self.ink_secondary,
            "ytick.labelcolor": self.ink_secondary,
            "xtick.labelsize": self.font_size * 0.9,
            "ytick.labelsize": self.font_size * 0.9,
            "grid.color": self.grid,
            "grid.linewidth": 0.5,
            "lines.linewidth": 2.0,
            "lines.markersize": 5.5,
            "legend.frameon": False,
            "legend.fontsize": self.font_size * 0.9,
            "legend.labelcolor": self.ink_secondary,
            "hatch.color": self.surface,
        }


TESPY = Style(
    name="tespy",
    surface="#ffffff",
    ink="#14212b",
    ink_secondary="#4d5a64",
    muted="#8a949c",
    grid="#e4e8eb",
    frame="#b9c1c7",
    categorical=(
        "#1b6298", "#e0662b", "#1f9d8d", "#b98900", "#8a4fa3", "#c8375a",
        "#3b4fc4", "#5a8a00"
    ),
    hot="#d1492e",
    cold="#1b6298",
    cycle="#1f567d",
    invalid="#d03b3b",
    sequential=("#eaf1f7", "#a9c6dd", "#5f94bd", "#1f567d", "#0e2f47"),
    diverging=("#1b6298", "#8fb3d1", "#ecebe8", "#e6a189", "#b8321f"),
    isolines={
        "T": "#d1492e", "p": "#8a4fa3", "h": "#b98900", "s": "#1b6298",
        "rho": "#1f9d8d", "u": "#c8375a"
    },
)

TESPY_DARK = Style(
    name="tespy-dark",
    surface="#161b1f",
    ink="#f1f4f6",
    ink_secondary="#c2cad0",
    muted="#87929a",
    grid="#272e34",
    frame="#465058",
    categorical=(
        "#3f8fd0", "#dd6a35", "#1fa797", "#b88a0c", "#9f6fd0", "#e0506e",
        "#5f78e8", "#6aa31a"
    ),
    hot="#e2674a",
    cold="#3f8fd0",
    cycle="#6fb3e8",
    invalid="#e66767",
    sequential=("#1a2a36", "#1f567d", "#5f94bd", "#a9c6dd", "#eaf1f7"),
    diverging=("#3f8fd0", "#2a4d6b", "#3a3c3e", "#6e3b2c", "#e2674a"),
    isolines={
        "T": "#e2674a", "p": "#9f6fd0", "h": "#b88a0c", "s": "#3f8fd0",
        "rho": "#1fa797", "u": "#e0506e"
    },
    logo="dark",
)

TESPY_PAPER = Style(
    name="tespy-paper",
    surface="#ffffff",
    ink="#000000",
    ink_secondary="#333333",
    muted="#7a7a7a",
    grid="#d9d9d9",
    frame="#000000",
    categorical=("#000000", "#1f567d", "#6b6b6b", "#9a9a9a"),
    hot="#000000",
    cold="#1f567d",
    cycle="#1f567d",
    invalid="#7a7a7a",
    sequential=("#f4f4f4", "#bdbdbd", "#737373", "#252525"),
    diverging=("#1f567d", "#8fa9bf", "#f1f1f1", "#a3a3a3", "#252525"),
    isolines=None,
    dashes=("solid", "dash", "dashdot", "dot"),
    hatches=("", "////", "....", "xxxx"),
    font_family=(
        "DejaVu Serif", "Times New Roman", "Liberation Serif", "serif"
    ),
)

_REGISTRY = {s.name: s for s in (TESPY, TESPY_PAPER, TESPY_DARK)}


def register_style(style, overwrite=False):
    """Register a style under its name.

    Parameters
    ----------
    style : Style
        The style.

    overwrite : bool
        Replace a registered style of the same name.
    """
    if style.name in _REGISTRY and not overwrite:
        msg = (
            f"A style named {style.name} exists, pass overwrite=True to "
            "replace it."
        )
        raise ValueError(msg)
    _REGISTRY[style.name] = style


def available_styles():
    """Names of the registered styles."""
    return list(_REGISTRY)


def get_style(style=None):
    """Return a style.

    Parameters
    ----------
    style : str, Style, optional
        Name of a registered style or a style, the default style if
        :code:`None` (see :py:func:`tespy.plotting.settings.set_defaults`).

    Returns
    -------
    Style
    """
    if isinstance(style, Style):
        return style
    if style is None:
        from tespy.plotting.settings import get_defaults
        style = get_defaults()["style"]
    try:
        return _REGISTRY[style]
    except KeyError:
        msg = (
            f"There is no style named {style}, the registered styles are "
            f"{', '.join(_REGISTRY)}."
        )
        raise KeyError(msg) from None
