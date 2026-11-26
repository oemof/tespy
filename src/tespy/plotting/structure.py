# -*- coding: utf-8

"""Plots of the structure of the system of equations: the incidence matrix
and its block lower triangular form.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/plotting/structure.py

SPDX-License-Identifier: MIT
"""
from tespy.plotting.description import Annotation
from tespy.plotting.description import Axis
from tespy.plotting.description import PlotDescription
from tespy.plotting.description import Series
from tespy.plotting.render import draw_matplotlib
from tespy.plotting.render import draw_plotly

#: Kinds of cells: no dependency, a dependency on a variable of a block
#: solved before, a dependency, and a dependency inside a block of several
#: equations that are solved together.
EMPTY = 0
EARLIER_BLOCK = 1
ENTRY = 2
COUPLED_BLOCK = 3
# room around the matrix, so that the outline of a block at its edge is not
# clipped
_EDGE = 0.1


def incidence_matrix_data(nw, block_order=False, title=None, info=None):
    """Collect the incidence matrix of the system of equations of a network:
    which equation depends on which variable.

    The matrix holds the equations and variables that remain after the
    presolving of the network.

    Parameters
    ----------
    nw : tespy.networks.network.Network
        Network with a built problem, e.g. after solving.

    block_order : bool
        Sort equations and variables by the blocks of the decomposition in
        solve order. The matrix is then block lower triangular: a block
        only depends on its own variables and on those of earlier blocks.
        Blocks of several equations that are solved together are
        highlighted and outlined, dependencies on earlier blocks are
        recessive.

    title : str, optional
        Title of the plot.

    info : list, optional
        Further lines of the information block.

    Returns
    -------
    tespy.plotting.description.PlotDescription
    """
    problem = nw.problem
    incidence = problem.get_incidence()
    equations = problem.get_equations()
    blocks = []
    if block_order:
        decomposition = problem.decompose()
        rows = [eq for block in decomposition.blocks for eq in block.equations]
        columns = [v for block in decomposition.blocks for v in block.variables]
        start = 0
        for block in decomposition.blocks:
            blocks.append({
                "id": block.id, "kind": block.kind, "start": start,
                "size": len(block.equations)
            })
            start += len(block.equations)
    else:
        rows = sorted(incidence)
        columns = sorted({v for deps in incidence.values() for v in deps})

    block_of_row = {
        row: block for block in blocks
        for row in range(block["start"], block["start"] + block["size"])
    }
    position = {v: j for j, v in enumerate(columns)}
    z = [[EMPTY] * len(columns) for _ in rows]
    for i, eq in enumerate(rows):
        block = block_of_row.get(i)
        for v in incidence[eq]:
            j = position[v]
            if block is None:
                z[i][j] = ENTRY
            elif j < block["start"]:
                z[i][j] = EARLIER_BLOCK
            else:
                z[i][j] = COUPLED_BLOCK if block["size"] > 1 else ENTRY

    row_labels = [
        f"{equations[eq][0]}.{problem.format_eq_name(equations[eq][1])}"
        for eq in rows
    ]
    column_labels = [problem.format_var_label(v) for v in columns]
    # the outline runs in the gaps between the cells
    annotations = [
        Annotation(
            "box",
            x=[b["start"] - 0.5, b["start"] + b["size"] - 0.5],
            y=[b["start"] - 0.5, b["start"] + b["size"] - 0.5],
        )
        for b in blocks if b["size"] > 1
    ]
    if title is None:
        title = (
            f"Block lower triangular form, {len(blocks)} blocks" if block_order
            else f"Incidence matrix, {len(rows)} equations"
        )
    return PlotDescription(
        kind="incidence_matrix",
        x_axis=Axis(
            "variable", range=[-0.5 - _EDGE, len(columns) - 0.5 + _EDGE],
            tick_labels=column_labels
        ),
        y_axis=Axis(
            "equation", range=[len(rows) - 0.5 + _EDGE, -0.5 - _EDGE],
            tick_labels=row_labels
        ),
        series=[Series(
            "dependency", x=list(range(len(columns))), y=list(range(len(rows))),
            z=z, kind="cells", show_legend=False
        )],
        annotations=annotations, title=title, info=list(info or []),
        extra={"blocks": blocks, "block_order": block_order},
    )


def plot_incidence_matrix(nw, block_order=False, ax=None, style=None, logo=None,
                          logo_location=None, title=None, info=None):
    """Plot the incidence matrix of a network or its block lower triangular
    form with matplotlib, see :py:func:`incidence_matrix_data`.

    Returns
    -------
    matplotlib.axes.Axes
    """
    description = incidence_matrix_data(nw, block_order, title, info)
    return draw_matplotlib(
        description, ax=ax, style=style, logo=logo, logo_location=logo_location
    )


def plot_incidence_matrix_plotly(nw, block_order=False, style=None, logo=None,
                                 logo_location=None, title=None, info=None,
                                 width=None, height=None):
    """Plot the incidence matrix of a network or its block lower triangular
    form with Plotly, see :py:func:`incidence_matrix_data`. Hovering an
    entry shows its equation and variable.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    description = incidence_matrix_data(nw, block_order, title, info)
    return draw_plotly(
        description, style=style, logo=logo, logo_location=logo_location,
        width=width, height=height
    )
