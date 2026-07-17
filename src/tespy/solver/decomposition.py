# -*- coding: utf-8

"""Module for the structural decomposition of the equation system.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/solver/decomposition.py

SPDX-License-Identifier: MIT
"""
from dataclasses import dataclass
from dataclasses import field

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.csgraph import maximum_bipartite_matching


@dataclass
class Block:
    r"""A block of the block lower triangular form of the equation system.

    - :code:`"scalar"`: a single equation in a single variable
    - :code:`"square"`: a group of mutually coupled equations and variables
    - :code:`"underdetermined"`: variables not covered by a maximum matching
      of the incidence together with the equations connected to them
    - :code:`"overdetermined"`: equations not covered by a maximum matching
      of the incidence together with the variables involved in them
    """
    id: int
    kind: str
    equations: list
    variables: list
    equation_labels: list = field(default_factory=list)
    variable_labels: list = field(default_factory=list)
    status: int = None
    residual_history: list = field(default_factory=list)


@dataclass
class Decomposition:
    r"""Ordered block structure of the equation system.

    The square and scalar blocks are ordered such that every block only
    depends on variables of preceding blocks. :code:`precedence` maps a block
    id to the ids of the blocks it depends on.
    """
    blocks: list
    precedence: dict

    @property
    def defective_blocks(self):
        return [
            block for block in self.blocks
            if block.kind in ("underdetermined", "overdetermined")
        ]


def dulmage_mendelsohn(incidence, num_variables):
    r"""Compute the Dulmage-Mendelsohn decomposition of an incidence.

    Parameters
    ----------
    incidence : dict
        Lookup of equation number to the list of variable numbers the
        equation depends on.
    num_variables : int
        Total number of variables.

    Returns
    -------
    Decomposition
        The over- and underdetermined parts (if any) and the square part
        split into its strongly connected blocks in precedence order.
    """
    num_equations = max(incidence, default=-1) + 1
    if num_equations == 0 and num_variables == 0:
        return Decomposition(blocks=[], precedence={})
    if num_equations == 0:
        block = Block(
            id=0, kind="underdetermined", equations=[],
            variables=list(range(num_variables))
        )
        return Decomposition(blocks=[block], precedence={0: []})
    if num_variables == 0:
        block = Block(
            id=0, kind="overdetermined",
            equations=list(range(num_equations)), variables=[]
        )
        return Decomposition(blocks=[block], precedence={0: []})

    eq_vars = {eq: set(dependents) for eq, dependents in incidence.items()}
    var_eqs = {}
    for eq, dependents in eq_vars.items():
        for col in dependents:
            var_eqs.setdefault(col, set()).add(eq)

    rows = [eq for eq, dependents in eq_vars.items() for _ in dependents]
    cols = [col for dependents in eq_vars.values() for col in dependents]
    matrix = csr_matrix(
        (np.ones(len(rows)), (rows, cols)),
        shape=(num_equations, num_variables)
    )
    matching = maximum_bipartite_matching(matrix, perm_type="column")
    col_match = np.full(num_variables, -1)
    for eq, col in enumerate(matching):
        if col >= 0:
            col_match[col] = eq

    # overdetermined part: alternating paths starting from unmatched equations
    over_eqs = {eq for eq in range(num_equations) if matching[eq] < 0}
    over_vars = set()
    stack = list(over_eqs)
    while stack:
        eq = stack.pop()
        for col in eq_vars.get(eq, []):
            if col in over_vars:
                continue
            over_vars.add(col)
            matched_eq = col_match[col]
            if matched_eq >= 0 and matched_eq not in over_eqs:
                over_eqs.add(matched_eq)
                stack.append(matched_eq)

    # underdetermined part: alternating paths starting from unmatched variables
    under_vars = {col for col in range(num_variables) if col_match[col] < 0}
    under_eqs = set()
    stack = list(under_vars)
    while stack:
        col = stack.pop()
        for eq in var_eqs.get(col, []):
            if eq in under_eqs:
                continue
            under_eqs.add(eq)
            matched_col = matching[eq]
            if matched_col >= 0 and matched_col not in under_vars:
                under_vars.add(matched_col)
                stack.append(matched_col)

    # square part: nodes are the matched equation/variable pairs, an edge
    # from the node owning a variable to every other node whose equation
    # depends on that variable defines the solve precedence
    square_eqs = sorted(
        eq for eq in range(num_equations)
        if eq not in over_eqs and eq not in under_eqs
    )
    node_of_eq = {eq: node for node, eq in enumerate(square_eqs)}
    num_nodes = len(square_eqs)

    edges = set()
    for eq in square_eqs:
        for col in eq_vars.get(eq, []):
            owner = col_match[col]
            if owner != eq and owner in node_of_eq:
                edges.add((node_of_eq[owner], node_of_eq[eq]))

    blocks = []
    precedence = {}
    if num_nodes > 0:
        if edges:
            edge_rows, edge_cols = zip(*edges)
        else:
            edge_rows, edge_cols = (), ()
        graph = csr_matrix(
            (np.ones(len(edges)), (edge_rows, edge_cols)),
            shape=(num_nodes, num_nodes)
        )
        _, labels = connected_components(
            graph, directed=True, connection="strong"
        )

        label_edges = {
            (labels[tail], labels[head])
            for tail, head in edges
            if labels[tail] != labels[head]
        }
        prerequisites = {}
        for tail, head in label_edges:
            prerequisites.setdefault(head, set()).add(tail)

        # topological order of the block condensation
        remaining = set(labels)
        order = []
        while remaining:
            free = sorted(
                label for label in remaining
                if not (prerequisites.get(label, set()) & remaining)
            )
            order += free
            remaining -= set(free)

        label_to_id = {label: block_id for block_id, label in enumerate(order)}
        for label in order:
            block_id = label_to_id[label]
            equations = sorted(
                int(square_eqs[node]) for node in range(num_nodes)
                if labels[node] == label
            )
            variables = sorted(int(matching[eq]) for eq in equations)
            blocks.append(Block(
                id=block_id,
                kind="scalar" if len(equations) == 1 else "square",
                equations=equations,
                variables=variables
            ))
            precedence[block_id] = sorted(
                label_to_id[tail] for tail, head in label_edges
                if head == label
            )

    if under_vars:
        block_id = len(blocks)
        blocks.append(Block(
            id=block_id,
            kind="underdetermined",
            equations=sorted(int(eq) for eq in under_eqs),
            variables=sorted(int(col) for col in under_vars)
        ))
        precedence[block_id] = []
    if over_eqs:
        block_id = len(blocks)
        blocks.append(Block(
            id=block_id,
            kind="overdetermined",
            equations=sorted(int(eq) for eq in over_eqs),
            variables=sorted(int(col) for col in over_vars)
        ))
        precedence[block_id] = []

    return Decomposition(blocks=blocks, precedence=precedence)
