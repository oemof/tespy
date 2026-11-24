# -*- coding: utf-8

"""TESPy logo and information block of plots.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/plotting/branding.py

SPDX-License-Identifier: MIT
"""
import base64
import importlib.metadata
import os

import tespy
from tespy.plotting.settings import resolve
from tespy.plotting.style import get_style

#: Width of the logo in points.
LOGO_WIDTH_PT = 60.0


def logo_path(variant="light", fmt="png"):
    """Path of a logo file shipped with TESPy.

    Parameters
    ----------
    variant : str
        :code:`"light"` for light and :code:`"dark"` for dark surfaces.

    fmt : str
        :code:`"png"` or :code:`"svg"`.
    """
    if variant not in ("light", "dark"):
        msg = f"The logo variant is 'light' or 'dark', not {variant}."
        raise ValueError(msg)
    if fmt not in ("png", "svg"):
        msg = f"The logo format is 'png' or 'svg', not {fmt}."
        raise ValueError(msg)
    return os.path.join(tespy.__datapath__, "logos", f"logo_tespy_{variant}.{fmt}")


def logo_data_uri(variant="light"):
    """The SVG logo as data URI, for HTML output."""
    with open(logo_path(variant, "svg"), "rb") as f:
        encoded = base64.b64encode(f.read()).decode()
    return f"data:image/svg+xml;base64,{encoded}"


def logo_aspect(variant="light"):
    """Width over height of the logo."""
    from matplotlib.image import imread
    height, width = imread(logo_path(variant, "png")).shape[:2]
    return width / height


def info_lines(lines=None):
    """Lines of the information block: the TESPy version and user lines.

    Parameters
    ----------
    lines : str, list, optional
        Further lines.
    """
    if isinstance(lines, str):
        lines = [lines]
    return [f"TESPy {importlib.metadata.version('tespy')}"] + list(lines or [])


def add_logo(ax, style=None, location=None, lines=None):
    """Add the logo and the information block to matplotlib axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw into.

    style : str, Style, optional
        Style of the plot, decides the logo variant and text colors.

    location : str, optional
        Corner of the axes, the default location if not given.

    lines : list, optional
        Further lines of the information block, below the TESPy version.

    Returns
    -------
    matplotlib.offsetbox.AnchoredOffsetbox
    """
    from matplotlib.image import imread
    from matplotlib.offsetbox import AnchoredOffsetbox
    from matplotlib.offsetbox import OffsetImage
    from matplotlib.offsetbox import TextArea
    from matplotlib.offsetbox import VPacker

    style = get_style(style)
    _, location = resolve(True, location)
    image = imread(logo_path(style.logo, "png"))
    zoom = LOGO_WIDTH_PT / image.shape[1] * 72 / ax.figure.dpi
    children = [OffsetImage(image, zoom=zoom, interpolation="hanning")]
    text = "\n".join(info_lines(lines))
    children.append(TextArea(text, textprops={
        "color": style.ink_secondary, "fontsize": style.font_size * 0.65
    }))
    box = AnchoredOffsetbox(
        loc=location, frameon=False, pad=0.3, borderpad=0.4,
        child=VPacker(children=children, align="left", pad=0, sep=2)
    )
    box.set_zorder(5)
    ax.add_artist(box)
    return box
