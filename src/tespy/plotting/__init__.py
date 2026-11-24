# -*- coding: utf-8

"""Plots of TESPy networks, components and results.

Each plot kind has a function that collects the data of the plot
(:code:`*_data`), a function that draws it with matplotlib (:code:`plot_*`)
and one that draws it with Plotly (:code:`plot_*_plotly`). matplotlib,
Plotly and CoolPlot are imported when a plot is drawn, not on import.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/plotting/__init__.py

SPDX-License-Identifier: MIT
"""
from tespy.plotting.cycles import plot_property_diagram  # noqa: F401
from tespy.plotting.cycles import plot_property_diagram_plotly  # noqa: F401
from tespy.plotting.cycles import property_diagram_data  # noqa: F401
from tespy.plotting.render import draw_matplotlib  # noqa: F401
from tespy.plotting.render import draw_plotly  # noqa: F401
from tespy.plotting.settings import get_defaults  # noqa: F401
from tespy.plotting.settings import reset_defaults  # noqa: F401
from tespy.plotting.settings import set_defaults  # noqa: F401
from tespy.plotting.style import Style  # noqa: F401
from tespy.plotting.style import available_styles  # noqa: F401
from tespy.plotting.style import get_style  # noqa: F401
from tespy.plotting.style import register_style  # noqa: F401
