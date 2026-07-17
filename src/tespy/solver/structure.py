# -*- coding: utf-8

"""Module for the structure graph of the equation system.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/solver/structure.py

SPDX-License-Identifier: MIT
"""
from tespy.tools import helpers as hlp


class StructureGraph:
    r"""
    Persistent graph representation of the equation system structure.

    Built from the structure matrix assembled during preprocessing, the graph
    holds three classes of structural information:

    - :code:`affine_adjacency`: variable pairs linked by a structure matrix
      row with exactly two entries, forming the affine elimination graph used
      by the presolver
    - :code:`linear_rows`: structure matrix rows with more or fewer than two
      entries (e.g. mass balances of merges and splitters), which cannot be
      used for pairwise elimination
    - :code:`nonlinear_incidence`: the incidence of the iterated solver
      equations, attached after the solver preparation via
      :code:`attach_solver_incidence`
    """

    def __init__(self, structure_matrix, rhs, variable_lookup,
                 equation_set_lookup, equation_set_origin):
        self._variable_lookup = variable_lookup
        self._equation_set_lookup = equation_set_lookup
        self._equation_set_origin = equation_set_origin

        self.affine_adjacency = {}
        self.edge_eq_idx = {}
        self.edge_offsets = {}
        self.edges_with_factors = []
        self.linear_rows = {}
        self.nonlinear_incidence = {}

        self._build(structure_matrix, rhs)

    def _build(self, sparse_matrix, rhs):
        # Extract edges and offsets from rows with two non-zero entries
        rows = {k[0] for k in sparse_matrix}
        # sorting needs to be applied to always have same orientation on edges
        # otherwise duplicate edges are not found if one is just in reverse
        rows_with_cols = {
            row: sorted([k[1] for k in sparse_matrix if k[0] == row])
            for row in rows
        }
        for row, cols in rows_with_cols.items():
            if len(cols) != 2:
                self.linear_rows[row] = cols
                continue

            non_zero_values = (
                sparse_matrix[(row, cols[0])], sparse_matrix[(row, cols[1])]
            )
            col1, col2 = cols
            val1, val2 = non_zero_values
            factor = -val1 / val2
            offset = rhs[row] / val2
            self.edges_with_factors.append((col1, col2, factor))
            self.edge_offsets[(col1, col2)] = offset
            if (col1, col2) in self.edge_eq_idx:
                variables = self._variable_labels([col1, col2])
                equations = self._equation_labels(
                    [self.edge_eq_idx[(col1, col2)], row]
                )
                var_str = ", ".join(f"{lbl} ({prop})" for lbl, prop in variables)
                eq_str = ", ".join(f"{lbl}.{eq}" for lbl, eq in equations)
                msg = (
                    "Two variables are directly linked by two equations. "
                    "This overdetermines the problem.\n"
                    f"  Variables:  {var_str}\n"
                    f"  Equations:  {eq_str}"
                )
                raise hlp.TESPyNetworkError(msg)

            self.edge_eq_idx[(col1, col2)] = row

        for col1, col2, factor in self.edges_with_factors:
            if col1 not in self.affine_adjacency:
                self.affine_adjacency[col1] = []
            if col2 not in self.affine_adjacency:
                self.affine_adjacency[col2] = []

            # Add edge with factor and reverse edge with reciprocal value
            self.affine_adjacency[col1].append((col2, factor))
            self.affine_adjacency[col2].append((col1, 1 / factor))

    def _variable_labels(self, number_list) -> list:
        return [
            (v["object"].label, v["property"])
            for k, v in self._variable_lookup.items()
            if k in number_list
        ]

    def _equation_labels(self, number_list) -> list:
        return [self._equation_set_lookup[num] for num in number_list]

    def find_cycle(self):
        """Find a cycle in the affine elimination graph.

        Returns
        -------
        set | None
            Set of the variable numbers forming a cycle or :code:`None` in
            case the graph is acyclic.
        """
        graph = {
            k: [x[0] for x in v] for k, v in self.affine_adjacency.items()
        }
        visited = set()
        parent = {}

        def dfs(node, prev):
            visited.add(node)
            for neighbor in graph.get(node, []):
                if neighbor not in visited:
                    parent[neighbor] = node
                    result = dfs(neighbor, node)
                    if result:
                        return result
                elif neighbor != prev:
                    # Cycle found, reconstruct it
                    cycle = [neighbor, node]
                    while cycle[-1] != neighbor:
                        cycle.append(parent[cycle[-1]])
                    cycle.reverse()
                    return cycle
            return None

        for node in graph:
            if node not in visited:
                parent[node] = None
                cycle = dfs(node, None)
                if cycle:
                    return set(cycle)

        return None

    def raise_error_if_cycle(self, cycle):
        edge_list = [
            e[:2] for e in self.edges_with_factors
            if e[0] in cycle or e[1] in cycle
        ]
        cycling_eqs = [
            v for k, v in self.edge_eq_idx.items() if k in edge_list
        ]
        variable_names = self._variable_labels(cycle)
        equations = self._equation_labels(cycling_eqs)
        var_str = ", ".join(f"{lbl} ({prop})" for lbl, prop in variable_names)
        eq_str = ", ".join(f"{lbl}.{eq}" for lbl, eq in equations)
        msg = (
            "A circular dependency has been detected. This overdetermines the problem.\n"
            f"  Variables:  {var_str}\n"
            f"  Equations:  {eq_str}"
        )
        raise hlp.TESPyNetworkError(msg)

    def find_linear_dependent_variables(self):
        """Find the sets of affinely linked variables.

        Returns
        -------
        list
            List of dictionaries per set of linked variables with the
            variables, the reference variable, the factors and offsets
            relative to the reference and the structural equations consumed
            by the elimination.
        """
        if len(self.affine_adjacency) == 0 and len(self.linear_rows) == 0:
            return []

        cycle = self.find_cycle()
        if cycle is not None:
            self.raise_error_if_cycle(cycle)

        adjacency_list = self.affine_adjacency
        eq_idx = self.edge_eq_idx
        rhs_offsets = self.edge_offsets

        # Find connected components and compute factors/offsets
        visited = set()
        variables_factors_offsets = []

        def dfs_component(node, current_factor, current_offset):
            """DFS to calculate factors and offsets relative to the reference
            variable."""
            stack = [(node, current_factor, current_offset)]
            factors = {node: current_factor}
            offsets = {node: current_offset}
            equation_indices = {}

            while stack:
                curr_node, curr_factor, curr_offset = stack.pop()
                visited.add(curr_node)

                for neighbor, edge_factor in adjacency_list.get(curr_node, []):
                    if neighbor not in factors:  # Process unvisited neighbor
                        # Calculate edge offset
                        idx_f = neighbor, curr_node
                        idx_b = curr_node, neighbor
                        edge_offset = (
                            rhs_offsets.get(idx_b, 0.0)
                            or -rhs_offsets.get(idx_f, 0.0)
                        )
                        # Determine which equation to use
                        if idx_f in rhs_offsets:
                            equation_indices[idx_f] = eq_idx[idx_f]
                        else:
                            equation_indices[idx_b] = eq_idx[idx_b]

                        # Compute new factor and offset
                        new_factor = curr_factor * edge_factor
                        new_offset = curr_offset * edge_factor + edge_offset

                        # Store and continue traversal
                        factors[neighbor] = new_factor
                        offsets[neighbor] = new_offset
                        stack.append((neighbor, new_factor, new_offset))

            return factors, offsets, equation_indices

        # Process each connected component
        for node in adjacency_list:
            if node not in visited:
                reference = node
                factors, offsets, equation_indices = dfs_component(
                    reference, 1.0, 0.0
                )

                variables_factors_offsets.append({
                    'variables': list(factors.keys()),
                    'reference': reference,
                    'factors': factors,
                    'offsets': offsets,
                    'equation_indices': equation_indices
                })

        return variables_factors_offsets

    def mass_flow_branches(self):
        """Find the branches of equal mass flow given by the topology.

        Only mass flow equalities originating from mandatory component
        constraints are followed. Mass flow edges imposed by specifications
        (e.g. referenced mass flow) do not merge physically separate
        branches. Connections without any topological mass flow link form a
        branch of their own.

        Returns
        -------
        list
            List of sets of structural variable numbers of the mass flow
            variables per branch.
        """
        m_cols = {
            col for col, data in self._variable_lookup.items()
            if data["property"] == "m"
        }
        adjacency = {col: [] for col in m_cols}
        for (col1, col2), row in self.edge_eq_idx.items():
            if (
                    col1 in m_cols and col2 in m_cols
                    and self._equation_set_origin[row] == "topology"
                ):
                adjacency[col1].append(col2)
                adjacency[col2].append(col1)

        branches = []
        visited = set()
        for col in m_cols:
            if col in visited:
                continue
            branch = set()
            stack = [col]
            while stack:
                current = stack.pop()
                if current in branch:
                    continue
                branch.add(current)
                stack += adjacency[current]
            visited |= branch
            branches.append(branch)

        return branches

    def attach_solver_incidence(self, incidence_matrix, variables_dict):
        """Attach the incidence of the iterated solver equations.

        The solver incidence maps equation numbers to solver variable columns.
        It is translated back to the structural variable space here, so all
        three edge classes of the graph refer to the same variable numbering.
        """
        self.nonlinear_incidence = {
            eq: sorted({
                sm_col
                for j_col in dependents
                for sm_col in variables_dict[j_col]["_represents"]
            })
            for eq, dependents in incidence_matrix.items()
        }

    def structural_key(self):
        """Get a hashable summary of the problem structure.

        The key is identical for two solves with unchanged topology and
        specification locations and changes when a specification moves,
        appears or vanishes. Specified values do not enter the key, so it can
        be used to cache structure-derived results (e.g. a decomposition)
        across parameter sweeps.
        """
        specified = []
        for col, data in self._variable_lookup.items():
            container = data["object"].get_attr(data["property"])
            is_set = container.is_set
            if isinstance(is_set, set):
                specified.append((col, frozenset(is_set)))
            elif is_set:
                specified.append((col, True))

        return (
            frozenset(self.edge_eq_idx.items()),
            frozenset(
                (row, tuple(cols)) for row, cols in self.linear_rows.items()
            ),
            frozenset(specified),
            frozenset(
                (eq, tuple(cols))
                for eq, cols in self.nonlinear_incidence.items()
            ),
        )
