# -*- coding: utf-8

"""Module for the solver problem representation.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/solver/problem.py

SPDX-License-Identifier: MIT
"""
import math
from time import time

import numpy as np
from numpy.linalg import norm

from tespy.connections.connection import ConnectionBase
from tespy.solver.decomposition import dulmage_mendelsohn
from tespy.solver.strategies import BlockDriver
from tespy.solver.strategies import search_reducing_step
from tespy.solver.structure import StructureGraph
from tespy.tools import helpers as hlp
from tespy.tools import logger
from tespy.tools.data_containers import ScalarVariable as dc_scavar
from tespy.tools.data_containers import VectorVariable as dc_vecvar
from tespy.tools.global_vars import ERR
from tespy.tools.global_vars import SCALED_INCREMENT_TOLERANCE
from tespy.tools.global_vars import SCALED_RESIDUAL_TOLERANCE

# weight of the phase expectation priors in the starting value least
# squares: low enough that any information transported by the guess edges
# dominates, only holding the absolute level where nothing else reaches
PRIOR_WEIGHT = 0.1

# Only require cupy if Cuda shall be used
try:
    import cupy as cu
except ModuleNotFoundError:
    cu = None


class Problem:
    r"""
    Solver representation of a prepared network.

    The problem holds the variable space, the equation lookups and the state
    of the newton algorithm (residual vector, jacobian and increment). It is
    created by the network in the preprocessing stage of every :code:`solve`
    call and remains accessible through :code:`Network.problem` afterwards
    for inspection.
    """

    def __init__(self, network):
        self.network = network

        self.variable_counter = 0
        self.variables_dict = {}

        self.num_comp_eq = 0
        self.num_conn_eq = 0
        self.num_ude_eq = 0

        self._equation_lookup = {}
        self._equation_obj_lookup = {}
        self._incidence_matrix = {}
        self._incidence_matrix_dense = None

        self.status = 99
        self.iter = 0
        self.residual = None
        self.jacobian = None
        self.increment = None
        self.residual_scales = None
        self.variable_weights = None
        self.lin_dep = False
        self.singularity_msg = ""

        self.structure_graph = None
        self._decomposition = None
        self._paused_driver = None

        self.num_original_variables = 0
        self.num_original_equations = 0

    def decompose(self):
        """Get the Dulmage-Mendelsohn decomposition of the solver incidence.

        The result is cached for the lifetime of the problem instance.

        Returns
        -------
        tespy.solver.decomposition.Decomposition
            Ordered block structure of the equation system.
        """
        if self._decomposition is None:
            # variable mass fractions of one fluid vector are coupled through
            # their normalization outside of the equation system, so they
            # must stay within one block. Block-wise solving falls back to
            # the simultaneous solution for variable compositions anyway -
            # this grouping keeps the structural diagnostics (determination
            # check, defective block reports, get_blocks) honest for such
            # networks
            coupled = {}
            for j_col, data in self.variables_dict.items():
                if data["variable"] == "fluid":
                    coupled.setdefault(id(data["obj"]), set()).add(j_col)
            groups = [group for group in coupled.values() if len(group) > 1]
            self._decomposition = dulmage_mendelsohn(
                self._incidence_matrix, self.variable_counter,
                coupled_variables=groups
            )
            for block in self._decomposition.blocks:
                block.equation_labels = [
                    self._equation_lookup[eq] for eq in block.equations
                ]
                block.variable_labels = [
                    self.format_var_label(col) for col in block.variables
                ]
        return self._decomposition

    def get_structural_analysis(self):
        """Get the over- and under-determined parts of the problem.

        Returns
        -------
        list
            One dictionary per defective part of the maximum matching of
            the incidence with its kind (:code:`"overdetermined"` or
            :code:`"underdetermined"`), the equations as tuples of object
            label and equation name and the short variable labels. An empty
            list means the problem is structurally sound.
        """
        return [
            {
                "kind": block.kind,
                "equations": list(block.equation_labels),
                "variables": list(block.variable_labels),
            }
            for block in self.decompose().defective_blocks
        ]

    def _structural_summary(self):
        parts = []
        for block in self.decompose().defective_blocks:
            num_eq = len(block.equations)
            num_var = len(block.variables)
            eq_word = "equation" if num_eq == 1 else "equations"
            var_word = "variable" if num_var == 1 else "variables"
            if block.kind == "overdetermined":
                parts.append(
                    f"an over-determined part ({num_eq} {eq_word} competing "
                    f"for {num_var} {var_word})"
                )
            else:
                parts.append(
                    f"an under-determined part ({num_var} {var_word} "
                    f"constrained by {num_eq} {eq_word})"
                )
        if not parts:
            return ""
        return (
            "\nStructural analysis - the problem contains "
            f"{' and '.join(parts)}."
        )

    def build(self):
        """Construct the problem from the prepared network.

        - Assemble the structure matrix and the structure graph and reduce the
          variable space by the affine linear dependencies.
        - Presolve fluid vectors and fluid properties iteratively.
        - Assign the solver variable space and collect the equations and their
          incidence.
        """
        self._create_structure_matrix()

        self.num_original_variables = len(self._variable_lookup)
        self.num_original_equations = len(self._equation_set_lookup)
        num_affine = len(self._presolved_equations)
        msg = (
            f"Original problem: {self.num_original_variables} variables and "
            f"{self.num_original_equations} equations, of which {num_affine} "
            "equations were consumed by the affine variable elimination."
        )
        logger.debug(msg)

        self._presolve()

        msg = (
            "Fluid property presolving consumed "
            f"{len(self._presolved_equations) - num_affine} further "
            "equations."
        )
        logger.debug(msg)

        self._prepare_for_solver()

        num_remaining_eq = (
            self.num_comp_eq + self.num_conn_eq + self.num_ude_eq
        )
        msg = (
            f"Reduced problem for the solver: {self.variable_counter} "
            f"variables and {num_remaining_eq} equations."
        )
        logger.debug(msg)

    def _create_structure_matrix(self):
        """Assemble the structural equations and eliminate affine groups.

        Enumerates all potential variables, collects the structural (symbolic
        linear) equation rows of every connection, component and user defined
        equation, and derives the affine variable groups from the resulting
        graph. Every group receives one shared reference container the members
        attach to with their :code:`factor` and :code:`offset`. The equations
        consumed by this elimination are recorded in
        :code:`self._presolved_equations`.
        """
        self._structure_matrix = {}
        self._rhs = {}
        self._variable_lookup = {}
        self._object_to_variable_lookup = {}
        self._equation_set_lookup = {}
        self._equation_set_origin = {}
        self._presolved_equations = []
        self._reference_container_lookup = {}

        num_vars = self._prepare_variables()

        self._reassign_ude_objects()

        sum_eq = 0
        sum_eq = self._preprocess_network_parts(self.network.conns["object"], sum_eq)
        sum_eq = self._preprocess_network_parts(self.network.comps["object"], sum_eq)
        sum_eq = self._preprocess_network_parts(self.network.user_defined_eq.values(), sum_eq)

        self.structure_graph = StructureGraph(
            self._structure_matrix, self._rhs,
            self._variable_lookup, self._equation_set_lookup,
            self._equation_set_origin
        )
        _linear_dependencies = (
            self.structure_graph.find_linear_dependent_variables()
        )
        _linear_dependent_variables = [
            var for linear_dependents in _linear_dependencies
            for var in linear_dependents["variables"]
        ]
        # variables without any affine partner become singleton groups
        # referencing themselves: every variable owns a reference
        # container this way and no consumer has to distinguish grouped
        # from ungrouped variables
        _missing_variables = [
            {
                "variables": [var],
                "reference": var,
                "factors": {var: 1.0},
                "offsets": {var: 0.0},
                "equation_indices": {}
            }
            for var in set(range(num_vars)) - set(_linear_dependent_variables)
        ]
        self._variable_dependencies = _missing_variables + _linear_dependencies

        for linear_dependents in self._variable_dependencies:
            reference_variable = self._variable_lookup[
                linear_dependents["reference"]
            ]
            reference = reference_variable["object"].get_attr(
                reference_variable["property"]
            )
            d = reference.d
            if hasattr(reference, "min_val"):
                min_val = reference.min_val
                max_val = reference.max_val
            else:
                min_val = None
                max_val = None

            if reference_variable["property"] != "fluid":
                DataContainer = dc_scavar
                reference_container = DataContainer(
                    _is_var=True,
                    _d=d,
                    min_val=min_val,
                    max_val=max_val,
                )
            else:
                DataContainer = dc_vecvar
                reference_container = DataContainer(
                    _d=d
                )

            self._reference_container_lookup[
                linear_dependents['reference']
            ] = reference_container

            for variable in linear_dependents["variables"]:
                variable_data = self._variable_lookup[variable]
                if variable_data["property"] != reference_variable["property"]:
                    msg = (
                        "There is a direct linear dependency between two "
                        "variables of different properties. This is unexpected "
                        "and might not work as intended: "
                        f"{reference_variable['object'].label}: "
                        f"{reference_variable['property']} and "
                        f"{variable_data['object'].label}: "
                        f"{variable_data['property']}."
                    )
                    logger.warning(msg)

                container = variable_data["object"].get_attr(variable_data["property"])
                container._reference_container = reference_container
                container._factor = linear_dependents["factors"][variable]
                container._offset = linear_dependents["offsets"][variable]

        # the user specifications were written to the local containers
        # before the reference containers existed: push them through the
        # affine relation now, so a value specified on any member of a
        # group reaches the shared reference
        for conn in self.network.conns["object"]:
            for prop, container in conn.get_variables().items():
                if conn.get_attr(prop).is_set:
                    conn.get_attr(prop).set_reference_val_SI(conn.get_attr(prop)._val_SI)

        # collect all presolved equations
        self._presolved_equations = [
            indices
            for dependents in self._variable_dependencies
            for indices in dependents["equation_indices"].values()
        ]

    def _prepare_variables(self):
        """Assign a structural column to every potential variable.

        Scalar connection properties, one column per fluid composition
        vector (the fractions are only split into individual solver
        variables later, after the vector groups are resolved) and the
        component variables. The :code:`_potential_var` flag marks what
        is not directly specified.

        Returns
        -------
        int
            Number of structural variables.
        """
        num_vars = 0
        for conn in self.network.conns["object"]:
            for prop, container in conn.get_variables().items():
                # flag potential variables
                container._potential_var = not container.is_set
                container.sm_col = num_vars
                num_vars += 1

                self._variable_lookup[container.sm_col] = {
                    "object": conn, "property": prop
                }
                if conn not in self._object_to_variable_lookup:
                    self._object_to_variable_lookup[conn] = {}
                self._object_to_variable_lookup[conn].update(
                    {prop: container.sm_col}
                )

            if hasattr(conn, "fluid"):
                # fluid is handled separately
                container = conn.fluid
                container.sm_col = num_vars
                num_vars += 1
                self._variable_lookup[container.sm_col] = {
                    "object": conn, "property": "fluid"
                }
                if conn not in self._object_to_variable_lookup:
                    self._object_to_variable_lookup[conn] = {}
                self._object_to_variable_lookup[conn].update(
                    {"fluid": container.sm_col}
                )

        for comp in self.network.comps["object"]:
            for prop, container in comp.get_variables().items():
                container.sm_col = num_vars
                num_vars += 1

                self._variable_lookup[container.sm_col] = {
                    "object": comp, "property": prop
                }
                if comp not in self._object_to_variable_lookup:
                    self._object_to_variable_lookup[comp] = {}
                self._object_to_variable_lookup[comp].update(
                    {prop: container.sm_col}
                )
        return num_vars

    def _reassign_ude_objects(self):
        """Re-point user defined equations to this network's objects.

        The user may have created the equation with objects of a
        different network instance, e.g. after loading a network from
        file: the labels identify the objects, the instances are looked
        up here.
        """
        for ude in self.network.user_defined_eq.values():
            ude.conns = [self.network.get_conn(c.label) for c in ude.conns]
            ude.comps = [self.network.get_comp(c.label) for c in ude.comps]

    def _preprocess_network_parts(self, parts, eq_counter):
        """Collect the structural equation rows of a group of objects.

        Every object writes its rows against a global row numbering
        starting at :code:`eq_counter`; the rows, right hand sides and
        the equation set lookups are merged into the problem level
        dictionaries and every equation set is classified by its origin.

        Parameters
        ----------
        parts : iterable
            Connection, component or user defined equation objects.

        eq_counter : int
            Global structural row number to continue from.

        Returns
        -------
        int
            Structural row number after the group's equations.
        """
        for obj in parts:
            obj._preprocess(eq_counter)
            self._structure_matrix.update(obj._structure_matrix)
            self._rhs.update(obj._rhs)
            eq_map = {
                eq_num: (obj.label, eq_name)
                for eq_num, eq_name in obj._equation_set_lookup.items()
            }
            self._equation_set_lookup.update(eq_map)
            # mandatory (and bypass) constraints exist independently of the
            # parametrization, all other equations are specification imposed
            self._equation_set_origin.update({
                eq_num: (
                    "topology"
                    if eq_name in getattr(obj, "constraints", {})
                    else "specification"
                )
                for eq_num, eq_name in obj._equation_set_lookup.items()
            })
            eq_counter += obj.num_eq

        return eq_counter

    def _presolve(self):
        """Determine variables from specifications in fixpoint rounds.

        After the fluid vectors are resolved and the user specifications
        are imposed on the affine groups, the connections check whether
        their specifications determine further variables (e.g. a
        temperature at known pressure fixes the enthalpy). The rounds
        are event driven: a group of variables that became known only
        triggers a recheck of the connections actually holding one of
        its members, until nothing new is determined.
        """
        self._presolve_fluid_vectors()

        # impose the user specifications on the affine groups first, so the
        # first round of connection presolving already sees them as known
        self._presolve_linear_dependents()

        # only scalar references are tracked: the fluid vectors were
        # fully resolved above and their is_var holds a set of fraction
        # names, not a boolean
        known_references = set()
        for linear_dependents in self._variable_dependencies:
            reference = linear_dependents["reference"]
            if self._variable_lookup[reference]["property"] == "fluid":
                continue
            container = self._reference_container_lookup[reference]
            if not container.is_var:
                known_references.add(id(container))

        # a determination on one connection can only enable further
        # determinations on the connections holding a variable of a group
        # that became known, so only those have to be checked again
        connections = set(self.network.conns['object'])
        to_check = self.network.conns['object'].tolist()
        rounds = 0
        while to_check:
            rounds += 1
            for c in to_check:
                self._presolved_equations += c._presolve()
            self._presolve_linear_dependents()

            newly_known = []
            for linear_dependents in self._variable_dependencies:
                reference = linear_dependents["reference"]
                if self._variable_lookup[reference]["property"] == "fluid":
                    continue
                container = self._reference_container_lookup[reference]
                if (
                        not container.is_var
                        and id(container) not in known_references
                    ):
                    known_references.add(id(container))
                    newly_known.append(linear_dependents)

            candidates = {}
            for group in newly_known:
                for variable in group["variables"]:
                    obj = self._variable_lookup[variable]["object"]
                    if obj in connections:
                        candidates[obj] = None
            to_check = list(candidates)

        msg = f"Fluid property presolving finished after {rounds} rounds."
        logger.debug(msg)

    def _presolve_fluid_vectors(self):
        """Resolve the fluid composition of every fluid vector group.

        Per affine group of fluid vectors the specified mass fractions
        of all members are combined: fully determined compositions are
        imposed, otherwise the unspecified fractions become solver
        variables and are registered in the variable space (one column
        per variable fraction).

        Limitations: factors and offsets between the vectors of a group
        are ignored, and a branch of constant composition consisting of
        a single connection is not detected as constant.
        """
        for linear_dependents in self._variable_dependencies:
            reference = linear_dependents["reference"]

            if self._variable_lookup[reference]["property"] != "fluid":
                continue

            all_connections = [
                self._variable_lookup[var]["object"]
                for var in linear_dependents["variables"]
            ]
            reference_container = self._reference_container_lookup[reference]
            reference_conn = all_connections[0]

            fluid_specs = [f for c in all_connections for f in c.fluid.is_set]
            fluid0 = {
                f: value for c in all_connections
                for f, value in c.fluid.val0.items()
            }
            if len(fluid_specs) == 0:

                if len(reference_conn._potential_fluids) > 1:
                    reference_container.is_var = {
                        f for f in reference_conn._potential_fluids
                    }
                    reference_container.val = {
                        f: 1 / len(reference_container.is_var)
                        for f in reference_container.is_var
                    }
                    # load up specification of starting values if any are
                    # available
                    reference_container.val.update(fluid0)
                else:
                    reference_container.val[
                        list(reference_conn._potential_fluids)[0]
                    ] = 1

            elif len(fluid_specs) != len(set(fluid_specs)):
                msg = (
                    "The mass fraction of a single fluid has been been "
                    "specified more than once in the following linear branch "
                    "of connections: "
                    f"{', '.join([c.label for c in all_connections])}."
                )
                raise hlp.TESPyNetworkError(msg)
            else:
                fixed_fractions = {
                    f: c.fluid._val[f]
                    for c in all_connections
                    for f in fluid_specs
                    if f in c.fluid._is_set
                }
                mass_fraction_sum = sum(fixed_fractions.values())
                if mass_fraction_sum > 1 + ERR:
                    msg = (
                        "The mass fraction of fluids within a linear branch "
                        "of connections cannot exceed 1: "
                        f"{', '.join([c.label for c in all_connections])}."
                    )
                    raise ValueError(msg)
                elif mass_fraction_sum < 1 - ERR:
                    # set the fluids with specified mass fraction
                    # remaining fluids are variable, create wrappers for them
                    all_fluids = reference_conn._potential_fluids
                    num_remaining_fluids = len(all_fluids) - len(fixed_fractions)
                    if num_remaining_fluids == 1:
                        missing_fluid = list(
                            set(all_fluids) - set(fixed_fractions.keys())
                        )[0]
                        fixed_fractions[missing_fluid] = 1 - mass_fraction_sum
                        variable = set()
                    else:
                        missing_fluids = (
                            set(all_fluids) - set(fixed_fractions.keys())
                        )
                        variable = {f for f in missing_fluids}

                else:
                    # fluid mass fraction is 100 %, all other fluids are 0 %
                    all_fluids = reference_container.val.keys()
                    remaining_fluids = (
                        reference_container.val.keys() - fixed_fractions.keys()
                    )
                    for f in remaining_fluids:
                        fixed_fractions[f] = 0

                    variable = set()

                reference_container.val.update(fixed_fractions)
                reference_container.is_set = {f for f in fixed_fractions}
                reference_container.is_var = variable
                # the remaining share is split evenly as starting value,
                # user provided starting fractions override the split.
                # NOTE: that override has caused convergence issues in
                # some models, e.g. the basic gas turbine tutorial
                num_var = len(variable)
                for f in variable:
                    reference_container.val[f] = (1 - mass_fraction_sum) / num_var
                    if f in fluid0:
                        reference_container.val[f] = fluid0[f]

            for fluid in reference_container.is_var:
                reference_container._J_col[fluid] = self.variable_counter
                self.variables_dict[self.variable_counter] = {
                    "obj": reference_container,
                    "variable": "fluid",
                    "fluid": fluid,
                    "_represents": linear_dependents["variables"]
                }
                self.variable_counter += 1

    def _presolve_linear_dependents(self):
        """Impose user specifications on their affine groups.

        A group with exactly one specified member is known: its
        reference stops being a variable. More than one specification
        inside a group over-determines the group and raises, naming the
        involved variables.
        """
        for linear_dependents in self._variable_dependencies:
            reference = linear_dependents["reference"]

            if self._variable_lookup[reference]["property"] == "fluid":
                continue

            all_containers = [
                self._variable_lookup[var]["object"].get_attr(
                    self._variable_lookup[var]["property"]
                ) for var in linear_dependents["variables"]
            ]

            number_specifications = sum(
                [not c._potential_var for c in all_containers]
            )
            if number_specifications > 1:
                variables_properties = [
                    f"{self._variable_lookup[var]['object'].label} "
                    f"({self._variable_lookup[var]['property']})"
                    for var in linear_dependents["variables"]
                ]
                var_str = ", ".join(variables_properties)
                msg = (
                    "You specified more than one variable within a set of "
                    "linearly dependent variables.\n"
                    f"  Variables:  {var_str}"
                )
                raise hlp.TESPyNetworkError(msg)
            elif number_specifications == 1:
                reference_data = self._variable_lookup[reference]
                reference_container = reference_data["object"].get_attr(
                    reference_data["property"]
                )._reference_container
                reference_container.is_var = False

    def _prepare_for_solver(self):
        """Assign the solver variable space and collect the equations.

        Every reference that is still a variable after presolving
        receives a solver column, the network parts contribute their
        iterated equations and incidence, and the solver incidence is
        attached to the structure graph.
        """
        for variable in self._variable_dependencies:
            reference = self._variable_lookup[variable["reference"]]
            represents = variable["variables"]
            if reference["property"] != "fluid":
                self._assign_variable_space(reference, represents)

        eq_counter = 0

        _eq_counter = self._prepare_network_parts(self.network.comps["object"], eq_counter)
        self.num_comp_eq = _eq_counter - eq_counter
        eq_counter = _eq_counter

        _eq_counter = self._prepare_network_parts(self.network.conns["object"], eq_counter)
        self.num_conn_eq = _eq_counter - eq_counter
        eq_counter = _eq_counter

        _eq_counter = self._prepare_network_parts(self.network.user_defined_eq.values(), eq_counter)
        self.num_ude_eq = _eq_counter - eq_counter
        eq_counter = _eq_counter

        self.structure_graph.attach_solver_incidence(
            self._incidence_matrix, self.variables_dict
        )

    def _prepare_network_parts(self, parts, eq_counter):
        """Collect the iterated equations and incidence of a group of
        objects.

        Every object numbers its remaining equations against the global
        counter, skipping the equation sets consumed by the affine
        elimination and the presolving
        (:code:`self._presolved_equations`). The incidence rows are
        assembled from the declared scalar and vector dependents.

        Parameters
        ----------
        parts : iterable
            Connection, component or user defined equation objects.

        eq_counter : int
            Global solver equation number to continue from.

        Returns
        -------
        int
            Solver equation number after the group's equations.
        """
        for obj in parts:
            eq_counter = obj._prepare_for_solver(self._presolved_equations, eq_counter)
            eq_map = {
                eq_num: (obj.label, eq_name)
                for eq_num, eq_name in obj._equation_lookup.items()
            }
            self._equation_lookup.update(eq_map)
            self._equation_obj_lookup.update(
                {eq_num: obj for eq_num in obj._equation_lookup}
            )

            dependents_map = {
                eq_num: [dependent.J_col for dependent in dependents]
                for eq_num, dependents in obj._equation_scalar_dependents_lookup.items()
            }
            self._incidence_matrix.update(dependents_map)

            dependents_map = {
                eq_num: [
                    dependent.J_col[key]
                    for dependent, keys in dependents.items()
                    for key in keys
                ]
                for eq_num, dependents in obj._equation_vector_dependents_lookup.items()
            }
            for eq_num, dependents in dependents_map.items():
                if eq_num in self._incidence_matrix:
                    self._incidence_matrix[eq_num] += dependents

                else:
                    self._incidence_matrix[eq_num] = dependents

        return eq_counter

    def _assign_variable_space(self, reference, represents):
        """Register a variable reference container in the solver space.

        The container receives its jacobian column and the entry records
        which structural variables the column represents.

        Parameters
        ----------
        reference : dict
            Entry of :code:`self._variable_lookup` for the group's
            reference variable, holding the owning object and the
            property name.

        represents : list
            Structural variable numbers the solver column represents.
        """
        container = reference["object"].get_attr(reference["property"])._reference_container
        if container.is_var:
            container.J_col = self.variable_counter
            self.variables_dict[self.variable_counter] = {
                "obj": container,
                "variable": reference["property"],
                "fluid": None,
                "_represents": represents
            }
            self.variable_counter += 1

    def unload_variables(self):
        """Detach all containers from their reference containers.

        Detaching writes the final value (through factor and offset)
        back into every member container, so the results remain after
        the shared references are released.
        """
        for dependents in self._variable_dependencies:
            for variable_num in dependents["variables"]:
                variable_dict = self._variable_lookup[variable_num]
                variable = variable_dict["object"].get_attr(variable_dict["property"])
                variable.detach()

    def starting_temperature_field(self):
        """Reconcile a temperature per connection for the starting values.

        The components provide approximate temperature relations between
        their ports, including the terminal temperature differences of
        heat exchangers - so known temperature levels reach coupled
        circuits the enthalpy relations cannot connect. The known
        temperatures - specifications, values derived from saturation
        specifications and converted fixed states - anchor a least squares
        problem per connected component of that graph; the temperature
        hints of the component state expectations act as soft priors.
        Returns a dict mapping connections to reconciled temperatures,
        including the known ones.
        """
        edges = []
        neighbors = {}
        for cp in self.network.comps["object"]:
            for c_from, c_to, offset, weight in cp._initial_temperature_edges():
                edges.append((c_from, c_to, offset, weight))
                neighbors.setdefault(c_from, []).append(c_to)
                neighbors.setdefault(c_to, []).append(c_from)

        fixed = {}
        priors = {}
        for c in self.network.conns["object"]:
            hint = c._temperature_hint()
            if hint is not None:
                fixed[c] = hint
                continue
            state = (
                c._declared_state() if hasattr(c, "_declared_state") else None
            )
            if state is not None and state.get("T") is not None:
                priors[c] = state["T"]

        field = dict(fixed)
        visited = set()
        for start in neighbors:
            if start in visited or start in fixed:
                continue
            component = []
            anchored = False
            stack = [start]
            visited.add(start)
            while stack:
                node = stack.pop()
                component.append(node)
                if node in priors:
                    anchored = True
                for neighbor in neighbors[node]:
                    if neighbor in fixed:
                        anchored = True
                    elif neighbor not in visited:
                        visited.add(neighbor)
                        stack.append(neighbor)

            if not anchored:
                continue

            index = {c: i for i, c in enumerate(component)}
            rows = []
            rhs = []
            for c_from, c_to, offset, weight in edges:
                from_in = c_from in index
                to_in = c_to in index
                if not from_in and not to_in:
                    continue
                row = np.zeros(len(component))
                if from_in and to_in:
                    row[index[c_to]] = 1.0
                    row[index[c_from]] = -1.0
                    residual = offset
                elif to_in:
                    if c_from not in fixed:
                        continue
                    row[index[c_to]] = 1.0
                    residual = offset + fixed[c_from]
                else:
                    if c_to not in fixed:
                        continue
                    row[index[c_from]] = -1.0
                    residual = offset - fixed[c_to]
                rows.append(row * weight)
                rhs.append(residual * weight)

            for c, value in priors.items():
                if c not in index:
                    continue
                row = np.zeros(len(component))
                row[index[c]] = PRIOR_WEIGHT
                rows.append(row)
                rhs.append(value * PRIOR_WEIGHT)

            solution = np.linalg.lstsq(
                np.array(rows), np.array(rhs), rcond=None
            )[0]
            for c, value in zip(component, solution):
                if not np.isnan(value):
                    field[c] = float(value)
        return field

    def propagate_starting_values(self, seeded, props):
        """Propagate starting values through approximate affine relations.

        The components provide rough affine guesses between their inlet and
        outlet properties, mapped into the space of the reference variables
        through the exact affine relations of the variable groups. The
        known values - specified, presolved or one of the passed source
        containers - anchor a least squares problem over every connected
        component of the guess graph, minimizing the sum of the squared
        edge residuals. Along a chain this reproduces the composed guesses,
        between several known values it interpolates, and around a cycle
        the closure residual of the guesses distributes over the edges
        instead of accumulating at an arbitrary meeting point. Returns the
        set of reference containers that received a value; unreached
        variables are left to the local component anchors and generic
        values.
        """
        properties = {}
        known = []
        for sm_col, container in self._reference_container_lookup.items():
            if type(container) == dc_vecvar:
                continue
            prop = self._variable_lookup[sm_col]["property"]
            properties[container] = prop
            if prop in props and not container.is_var:
                known.append(container)

        edges = []
        neighbors = {}
        for cp in self.network.comps["object"]:
            for edge in cp._initial_affine_edges():
                c_from, c_to, factor, offset = edge[:4]
                weight = edge[4] if len(edge) > 4 else 1.0
                ref_from = c_from._reference_container
                ref_to = c_to._reference_container
                if ref_from is None or ref_to is None or ref_from is ref_to:
                    continue
                if (
                        properties.get(ref_from) not in props
                        or properties.get(ref_to) not in props
                    ):
                    continue
                # with member = f * ref + o on both sides, the approximate
                # relation to = factor * from + offset maps to reference space
                a = factor * c_from._factor / c_to._factor
                b = (
                    factor * c_from._offset + offset - c_to._offset
                ) / c_to._factor
                edges.append((ref_from, ref_to, a, b, weight))
                neighbors.setdefault(ref_from, []).append(ref_to)
                neighbors.setdefault(ref_to, []).append(ref_from)

        fixed = {}
        for ref in dict.fromkeys(seeded + known):
            if properties.get(ref) in props and not np.isnan(ref.val_SI):
                fixed[ref] = ref.val_SI

        # the phase expectations of the components act as soft priors: they
        # hold the absolute level where no fixed value reaches, but any
        # real information transported by the edges dominates them
        priors = {}
        if "h" in props:
            for c in self.network.conns["object"]:
                if not hasattr(c, "_state_prior"):
                    continue
                ref = c.h._reference_container
                if (
                        ref is None or not getattr(ref, "is_var", False)
                        or ref in fixed or ref in priors
                    ):
                    continue
                value = c._state_prior()
                if value is not None:
                    priors[ref] = value
                    neighbors.setdefault(ref, [])

        assigned = set()
        visited = set()
        for start in neighbors:
            if start in visited or start in fixed:
                continue
            # collect the connected component of unknowns around this node,
            # the known values act as its boundary
            component = []
            anchored = False
            stack = [start]
            visited.add(start)
            while stack:
                node = stack.pop()
                component.append(node)
                if node in priors:
                    anchored = True
                for neighbor in neighbors[node]:
                    if neighbor in fixed:
                        anchored = True
                    elif neighbor not in visited:
                        visited.add(neighbor)
                        stack.append(neighbor)

            if not anchored:
                continue

            index = {ref: i for i, ref in enumerate(component)}
            rows = []
            rhs = []
            for ref_from, ref_to, a, b, weight in edges:
                from_in = ref_from in index
                to_in = ref_to in index
                if not from_in and not to_in:
                    continue
                row = np.zeros(len(component))
                if from_in and to_in:
                    row[index[ref_to]] = 1.0
                    row[index[ref_from]] = -a
                    residual = b
                elif to_in:
                    if ref_from not in fixed:
                        continue
                    row[index[ref_to]] = 1.0
                    residual = b + a * fixed[ref_from]
                else:
                    if ref_to not in fixed:
                        continue
                    row[index[ref_from]] = -a
                    residual = b - fixed[ref_to]
                rows.append(row * weight)
                rhs.append(residual * weight)

            for ref, value in priors.items():
                if ref not in index:
                    continue
                row = np.zeros(len(component))
                row[index[ref]] = PRIOR_WEIGHT
                rows.append(row)
                rhs.append(value * PRIOR_WEIGHT)

            solution = np.linalg.lstsq(
                np.array(rows), np.array(rhs), rcond=None
            )[0]
            if any(
                np.isnan(value)
                or (properties[ref] in ("p", "m") and value <= 0)
                for ref, value in zip(component, solution)
            ):
                continue
            for ref, value in zip(component, solution):
                ref._val_SI = float(value)
                assigned.add(ref)
        return assigned

    def presolve_flow_variables(self, seeded):
        """Solve the flow variables on the completed property field.

        With pressures and enthalpies frozen at their starting values,
        the equations determining the flows - mass balances, energy
        balances, heat and power specifications, energy connector and
        bus balances - are linear in the mass and energy flow
        variables. A single newton step restricted to those variables
        therefore solves them exactly with respect to the guessed
        property field, replacing the random mass flow starting values
        by consistent levels and branch ratios. Seeded values stay
        untouched. Returns the number of updated variables.
        """
        columns = [
            key for key, data in self.variables_dict.items()
            if data["variable"] in ("m", "E") and data["obj"] not in seeded
        ]
        if not any(
                self.variables_dict[col]["variable"] == "m"
                for col in columns
            ):
            return 0

        column_set = set(columns)
        # the frozen property field only reconciles pressures and
        # enthalpies - equations depending on a guessed fluid composition
        # (combustion stoichiometry and energy balances) would inject
        # arbitrary flow information and are left to the solver
        fluid_columns = set()
        for data in self.variables_dict.values():
            if data["variable"] == "fluid":
                fluid_columns.update(data["obj"].J_col.values())
        rows = [
            row for row, cols in self._incidence_matrix.items()
            if not column_set.isdisjoint(cols)
            and fluid_columns.isdisjoint(cols)
        ]
        if not rows:
            return 0

        n = self.variable_counter
        self.residual = np.zeros(n)
        self.jacobian = np.zeros((n, n))
        self.increment = np.ones(n)
        self.residual_scales = np.ones(n)
        self.variable_weights = np.ones(n)
        self.increment_filter = np.absolute(self.increment) < ERR ** 2
        try:
            self._solve_equations()
        except Exception as e:
            # best effort pass: any failure of the evaluation on the
            # guessed values - including badly parametrized equations,
            # which raise their descriptive errors in the actual solve -
            # just skips the flow presolve
            msg = (
                "Skipping the flow variable presolve, the equation "
                f"evaluation on the starting values failed: {e}"
            )
            logger.debug(msg)
            return 0
        finally:
            for obj in (
                    self.network.comps["object"].tolist()
                    + self.network.conns["object"].tolist()
                ):
                obj.it = 0

        rows = [
            row for row in rows
            if np.isfinite(self.residual[row])
            and np.isfinite(self.jacobian[row, columns]).all()
        ]
        if not rows:
            return 0

        matrix = self.jacobian[np.ix_(rows, columns)]
        rhs = -self.residual[rows]
        # equilibrate the rows: energy balances live on a scale of watts,
        # mass balances on kilograms per second - without normalization
        # the least squares compromise would sacrifice the exact balances
        # to the approximate energy rows
        norms = np.linalg.norm(matrix, axis=1)
        keep = norms > 0
        if not keep.any():
            return 0
        matrix = matrix[keep] / norms[keep, None]
        rhs = rhs[keep] / norms[keep]

        delta = np.linalg.lstsq(matrix, rhs, rcond=None)[0]
        if not np.isfinite(delta).all():
            return 0        # deliberately stay short of the exact solution of the linearized
        # system: some models have degenerate jacobians exactly on the
        # solution manifold of the flow equations (recirculation with a
        # specified split ratio), the same way starting enthalpies are kept
        # off the saturation lines
        delta *= 0.9

        for col, increment in zip(columns, delta):
            data = self.variables_dict[col]
            value = data["obj"]._val_SI + float(increment)
            # a previously positive mass flow crossing zero means the
            # retained rows underdetermine this flow (e.g. a bare total
            # mass balance distributing its correction arbitrarily) - it
            # keeps a fraction of its previous value instead
            if (
                    data["variable"] == "m" and data["obj"]._val_SI > 0
                    and value <= 0
                ):
                value = data["obj"]._val_SI * 0.1
            data["obj"]._val_SI = value
        return len(columns)

    def check_determination(self):
        r"""Check, if the number of supplied parameters is sufficient."""
        _hint = (
            "\nInspect the problem with\n"
            "  - nw.print_structural_analysis() for the over- or "
            "under-determined part of the model,\n"
            "  - nw.print_variables() and nw.print_equations() for the "
            "variables and equations that are present,\n"
            "  - nw.print_equations_with_dependents() for the variables "
            "each equation depends on,\n"
            "  - nw.print_incidence_matrix() for a compact overview."
        )
        n = self.num_comp_eq + self.num_conn_eq + self.num_ude_eq
        _reduction = (
            f"\nThe original problem holds {self.num_original_variables} "
            f"variables and {self.num_original_equations} equations. "
            f"Presolving consumed {len(self._presolved_equations)} equations "
            f"and reduced the problem to {self.variable_counter} variables "
            f"and {n} equations."
        )
        if n > self.variable_counter:
            msg = (
                f"You have provided too many parameters: {self.variable_counter} "
                f"required, {n} supplied. Aborting calculation!{_reduction}"
                f"{self._structural_summary()}{_hint}"
            )
            logger.error(msg)
            self.status = 12
            raise hlp.TESPyNetworkError(msg)
        elif n < self.variable_counter:
            msg = (
                f"You have not provided enough parameters: {self.variable_counter} "
                f"required, {n} supplied. Aborting calculation!{_reduction}"
                f"{self._structural_summary()}{_hint}"
            )
            logger.error(msg)
            self.status = 11
            raise hlp.TESPyNetworkError(msg)

    def no_progress_message(self):
        return (
            "The solver does not seem to make any progress, aborting "
            "calculation. Scaled residual value is "
            "{:.2e}{}".format(
                self.scaled_residual().max(), self._worst_equation_hint()
            ) +
            "\nPossible reasons include:\n"
            " - fluid properties moving outside the valid range of the "
            "property database (consider adjusting p_range or h_range),\n"
            " - an impossible constraint that can never be satisfied \n"
            " - bad starting values causing the Newton solver to diverge.\n"
            "Use nw.print_residuals() to identify which equations have "
            "the largest residuals."
        )

    @property
    def paused(self):
        """Whether the solution process is paused at a failed block."""
        return self._paused_driver is not None

    def _release_pause(self):
        if self._paused_driver is not None:
            self._paused_driver._unfreeze()
            self._paused_driver = None

    def solve_loop(self, max_iter, min_iter, use_cuda, robust_relax,
                   oscillation_damping, iterinfo, print_results=True,
                   block_solve=False, pause_on_block_failure=None):
        r"""Loop of the newton algorithm."""
        self.max_iter = max_iter
        self.min_iter = min_iter
        self.use_cuda = use_cuda
        self.robust_relax = robust_relax
        self.oscillation_damping = oscillation_damping
        self._auto_damping = False

        if self._paused_driver is not None:
            driver = self._paused_driver
            if pause_on_block_failure is not None:
                driver.pause_on_block_failure = pause_on_block_failure
            self.start_time = time()
            driver.resume(iterinfo, print_results)
            self.end_time = time()
            return

        if pause_on_block_failure and not block_solve:
            msg = (
                "Pausing on a failed block is not possible for the "
                "simultaneous solution, the option has no effect."
            )
            logger.info(msg)

        n = self.variable_counter
        self._incidence_matrix_dense = np.zeros((n, n))
        for row, cols in self._incidence_matrix.items():
            self._incidence_matrix_dense[row, cols] = 1

        # parameter definitions
        self.residual_history = np.array([])
        self.residual = np.zeros([self.variable_counter])
        self.increment = np.ones([self.variable_counter])
        self.jacobian = np.zeros((self.variable_counter, self.variable_counter))
        self.residual_scales = np.ones([self.variable_counter])
        self.variable_weights = np.ones([self.variable_counter])
        self._prev_residual = None

        self.start_time = time()

        if block_solve:
            BlockDriver(
                self, pause_on_block_failure=bool(pause_on_block_failure)
            ).solve(iterinfo, print_results)
        else:
            self._solve_simultaneous(iterinfo, print_results)

        self.end_time = time()

    def _solve_simultaneous(self, iterinfo, print_results=True):
        if iterinfo:
            self._print_iterinfo_head(print_results)

        self._auto_damping = False
        for self.iter in range(self.max_iter):
            self.increment_filter = np.absolute(self.increment) < ERR ** 2
            modified = self._solve_iteration()
            self.residual_history = np.append(
                self.residual_history, self.scaled_residual().max()
            )
            if (
                    self.iter >= 3 and not self._auto_damping
                    and self._detect_residual_oscillation(
                        self.residual_history[-4:]
                    )
                ):
                self._auto_damping = True
                msg = (
                    "The residual norm alternates between two values, "
                    "dampening the oscillating increments."
                )
                logger.debug(msg)
            if iterinfo:
                self._print_iterinfo_body(print_results)

            if self.lin_dep:
                self.status = 3
                break

            elif self.iter > 40:
                if (
                    all(
                        self.residual_history[(self.iter - 3):] >= self.residual_history[-3] * 0.95
                    ) and self.residual_history[-1] >= self.residual_history[-2] * 0.95
                ):
                    self.status = 2
                    break

            # accept once two consecutive residuals and the increment are
            # small in an iteration without interference of the value
            # heuristics
            if (
                    self.iter >= self.min_iter - 1
                    and not modified
                    and (
                        self.residual_history[-2:]
                        < SCALED_RESIDUAL_TOLERANCE
                    ).all()
                    # the increment should also be small, but it does not need
                    # to be that small
                    and (
                        abs(self.increment) / self.variable_weights
                        < SCALED_INCREMENT_TOLERANCE
                    ).all()
                ):
                self.status = 0
                break

        self.end_time = time()

        if iterinfo:
            self._print_iterinfo_tail(print_results)

        if self.iter == self.max_iter - 1:
            msg = (
                f"Reached maximum iteration count ({self.max_iter}), "
                "calculation stopped. Scaled residual value is "
                "{:.2e}{}. ".format(
                    self.scaled_residual().max(),
                    self._worst_equation_hint()
                ) +
                "\nPossible reasons include:\n"
                " - fluid properties moving outside the valid range of the "
                "property database (consider adjusting p_range or h_range),\n"
                " - an impossible constraint that can never be satisfied \n"
                " - bad starting values causing the Newton solver to diverge.\n"
                "Use nw.print_residuals() to identify which equations have "
                "the largest residuals."
            )
            logger.warning(msg)
            self.status = 2

    def _print_iterinfo_head(self, print_results=True):
        """Print head of convergence progress."""
        # Start with defining the format here
        self.iterinfo_fmt = ' {iter:5s} | {residual:10s} | {progress:10s} '
        self.iterinfo_fmt += '| {massflow:10s} | {pressure:10s} | {enthalpy:10s} '
        self.iterinfo_fmt += '| {fluid:10s} | {energy:10s} | {component:10s} '
        # Use the format to create the first logging entry
        msg = self.iterinfo_fmt.format(
            iter='iter',
            residual='residual',
            progress='progress',
            massflow='massflow',
            pressure='pressure',
            enthalpy='enthalpy',
            fluid='fluid',
            energy='energy',
            component='component'
        )
        logger.progress(0, msg)
        msg2 = '-' * 7 + '+------------' * 8

        logger.progress(0, msg2)
        if print_results and not logger.console_logging_enabled():
            print('\n' + msg + '\n' + msg2)

    def _print_iterinfo_body(self, print_results=True):
        """Print convergence progress."""
        m = [k for k, v in self.variables_dict.items() if v["variable"] == "m"]
        p = [k for k, v in self.variables_dict.items() if v["variable"] == "p"]
        h = [k for k, v in self.variables_dict.items() if v["variable"] == "h"]
        fl = [k for k, v in self.variables_dict.items() if v["variable"] == "fluid"]
        e = [k for k, v in self.variables_dict.items() if v["variable"] == "E"]
        cp = [k for k in self.variables_dict if k not in m + p + h + fl + e]

        iter_str = str(self.iter + 1)
        residual_norm = self.scaled_residual().max()
        residual = 'NaN'
        progress = 'NaN'
        massflow = 'NaN'
        pressure = 'NaN'
        enthalpy = 'NaN'
        fluid = 'NaN'
        energy = 'NaN'
        component = 'NaN'

        progress_val = -1

        if not np.isnan(residual_norm):
            residual = '{:.2e}'.format(residual_norm)

            if not self.lin_dep:
                massflow = '{:.2e}'.format(norm(self.increment[m]))
                pressure = '{:.2e}'.format(norm(self.increment[p]))
                enthalpy = '{:.2e}'.format(norm(self.increment[h]))
                fluid = '{:.2e}'.format(norm(self.increment[fl]))
                energy  = '{:.2e}'.format(norm(self.increment[e]))
                component  = '{:.2e}'.format(norm(self.increment[cp]))

            # 0 % at a scaled residual of 1 (variables off by their own
            # magnitude), 100 % at the acceptance tolerance
            if residual_norm > np.finfo(float).eps * 100:
                progress_scaled = (
                    math.log(max(residual_norm, SCALED_RESIDUAL_TOLERANCE))
                    / math.log(SCALED_RESIDUAL_TOLERANCE)
                )
                progress_val = int(max(0, min(1, progress_scaled)) * 100)
            else:
                progress_val = 100

            progress = '{:d} %'.format(progress_val)

        msg = self.iterinfo_fmt.format(
            iter=iter_str,
            residual=residual,
            progress=progress,
            massflow=massflow,
            pressure=pressure,
            enthalpy=enthalpy,
            fluid=fluid,
            energy=energy,
            component=component
        )
        logger.progress(progress_val, msg)
        if print_results and not logger.console_logging_enabled():
            print(msg)

    def _print_iterinfo_tail(self, print_results=True):
        """Print tail of convergence progress."""
        num_iter = self.iter + 1
        clc_time = self.end_time - self.start_time
        num_ips = num_iter / clc_time if clc_time > 1e-10 else np.inf
        msg = '-' * 7 + '+------------' * 7
        logger.progress(100, msg)
        msg = (
            "Total iterations: {0:d}, Calculation time: {1:.2f} s, "
            "Iterations per second: {2:.2f}"
        ).format(num_iter, clc_time, num_ips)
        logger.debug(msg)
        if print_results and not logger.console_logging_enabled():
            print(msg)
        return

    def _fill_jacobian_surrogates(self, equations=None, variables=None):
        """Restore invertibility for all-zero rows and find better steps.

        For each row that is entirely zero but expected to have non-zero
        entries (per the incidence matrix), inserts 1 in the expected positions
        so the Jacobian can be inverted for all other variables.
        Subsequently searches value of associated variable(s) to find the
        increment for the affected variable(s) that reduces that equation's
        residual.

        Returns a dict {col: step} of increment overrides to apply after the
        inversion.
        """
        overrides = {}
        if equations is None:
            zero_rows = self._check_all_zero_rows(self.jacobian)
            columns = None
        else:
            # in the restricted case only the block's columns matter: both
            # the zero row check and the surrogate placeholders must ignore
            # other columns, otherwise a stale placeholder of a non block
            # column masks the all-zero row in subsequent iterations
            zero_rows = [
                row for row in equations
                if not self.jacobian[row, variables].any()
            ]
            columns = set(variables)
        for row in zero_rows:
            for col in self._incidence_matrix.get(row, []):
                if columns is not None and col not in columns:
                    continue
                if self.jacobian[row, col] == 0.0:
                    self.jacobian[row, col] = 1.0
                    step = search_reducing_step(self, row, col)
                    if step is not None:
                        overrides[col] = step
        return overrides

    def _invert_jacobian(self, equations=None, variables=None):
        """Compute Newton step, storing result in self.increment. Sets self.lin_dep."""
        self.lin_dep = False
        self.increment = self.residual * 0
        if len(self.variables_dict) == 0:
            return

        self._check_residual_and_jacobian_for_nan(equations)

        overrides = self._fill_jacobian_surrogates(equations, variables)

        if equations is None:
            jacobian = self.jacobian
            residual = self.residual
        else:
            jacobian = self.jacobian[np.ix_(equations, variables)]
            residual = self.residual[equations]

        try:
            if self.use_cuda and equations is None:
                self.increment = cu.asnumpy(cu.dot(
                    cu.linalg.inv(cu.asarray(self.jacobian)),
                    -cu.asarray(self.residual)
                ))
            else:
                row_scales = np.abs(jacobian).max(axis=1)
                row_scales[row_scales == 0] = 1.0
                J_eq = jacobian / row_scales[:, None]

                col_scales = np.maximum(np.abs(J_eq).max(axis=0), 1e-10)
                J_sc = J_eq / col_scales[None, :]
                r_eq = -residual / row_scales

                increment = np.linalg.solve(J_sc, r_eq) / col_scales
                if equations is None:
                    self.increment = increment
                else:
                    self.increment[variables] = increment
        except np.linalg.LinAlgError:
            self.lin_dep = True
            return

        for col, step in overrides.items():
            self.increment[col] = step

    def _check_residual_and_jacobian_for_nan(self, equations=None):
        """Raise an informative error if a NaN or inf entry is found.

        A single NaN in the residual or Jacobian typically contaminates the
        whole solution vector once the linear system is solved, so the bad
        value has to be located here, before the solve, to be able to point
        at the equation that produced it. Otherwise it propagates silently
        into the next iteration's variable values and only surfaces several
        iterations later as an unrelated and confusing error.
        """
        if equations is None:
            nan_rows = set(np.where(~np.isfinite(self.residual))[0])
            nan_rows |= set(np.where(~np.isfinite(self.jacobian).any(axis=1))[0])
        else:
            nan_rows = {
                row for row in equations
                if not np.isfinite(self.residual[row])
                or not np.isfinite(self.jacobian[row]).any()
            }
        if len(nan_rows) == 0:
            return

        eq_str = ", ".join(
            f"{lbl}.{self.format_eq_name(eq_name)}"
            for lbl, eq_name in self._get_equations_by_number(nan_rows).values()
        )
        msg = (
            "The residual or one of the partial derivatives of the "
            f"following equation(s) is NaN or inf: {eq_str}. Please check "
            "the inputs to and the implementation of these equations, e.g. "
            "a UserDefinedEquation."
        )
        logger.error(msg)
        raise hlp.TESPyNetworkError(msg)

    def _diagnose_singularity(self):
        """Build singularity_msg after a failed matrix solve."""
        if self.iter == 0:
            summary = self._structural_summary()
            if summary != "":
                self.singularity_msg = (
                    "Detected singularity in Jacobian matrix. This "
                    "singularity is caused by the structure of your problem "
                    "and NOT a numerical issue. Double check your setup."
                    f"{summary}\nUse nw.print_structural_analysis() for the "
                    "affected equations and variables.\n"
                )
                self._find_linear_dependencies(self.jacobian)
                return

        if self.iter == 0 and np.linalg.matrix_rank(self._incidence_matrix_dense) < self._incidence_matrix_dense.shape[0]:
            self.singularity_msg = (
                "Detected singularity in Jacobian matrix. This singularity "
                "is most likely caused by the parametrization of your "
                "problem and NOT a numerical issue. Double check your "
                "setup.\n"
            )
            self._find_linear_dependencies(self.jacobian)
            return

        expected_entries = self._incidence_matrix_dense.astype(bool)
        actual_entries = self.jacobian.astype(bool)
        rows, cols = np.where(expected_entries != actual_entries)

        missing_entries = []
        for row, col in zip(rows, cols):
            lbl, eq_name = self._equation_lookup[row]
            eq_str = f"{lbl}.{self.format_eq_name(eq_name)}"
            var_str = self.format_var_label(col)
            missing_entries += [f"{eq_str}: {var_str}"]

        entries_str = ", ".join(missing_entries)
        self.singularity_msg = (
            "Found singularity in Jacobian matrix, calculation aborted! "
            "The setup of your problem seems to be solvable. It failed "
            "due to partial derivatives in the Jacobian being zero where "
            "a non-zero was expected, or vice versa. This usually lies in "
            "starting value selection or bad convergence.\n"
            "  The following equation/variable pairs may have an "
            f"unexpected zero/non-zero partial derivative:  {entries_str}\n"
        )
        self._find_linear_dependencies(self.jacobian)

    def _find_linear_dependencies(self, matrix):
        all_zero_cols = self._check_all_zero_columns(matrix)
        all_zero_rows = self._check_all_zero_rows(matrix)
        if len(all_zero_cols) + len(all_zero_rows) == 0:
            eq_indices = self._cauchy_schwarz_inequality(matrix)
            eq_str = ", ".join(
                f"{lbl}.{self.format_eq_name(eq_name)}"
                for lbl, eq_name in self._get_equations_by_number(eq_indices).values()
            )
            self.singularity_msg += (
                "The following equations form a linear dependency:\n"
                f"  {eq_str}\n"
            )
        else:
            if len(all_zero_cols) > 0:
                var_str = ", ".join(self.format_var_label(i) for i in all_zero_cols)
                self.singularity_msg += (
                    "The following variables are not associated with any equation:\n"
                    f"  {var_str}\n"
                )
            if len(all_zero_rows) > 0:
                eq_str = ", ".join(
                    f"{lbl}.{self.format_eq_name(eq_name)}"
                    for lbl, eq_name in self._get_equations_by_number(all_zero_rows).values()
                )
                self.singularity_msg += (
                    "The following equations do not depend on any variable:\n"
                    f"  {eq_str}\n"
                )

    def _check_all_zero_columns(self, matrix):
        return np.where((matrix == 0).all(axis=0))[0]

    def _check_all_zero_rows(self, matrix):
        return np.where((matrix == 0).all(axis=1))[0]

    def _cauchy_schwarz_inequality(self, matrix):
        n = matrix.shape[0]
        dependent_equations = []
        for i in range(n):
            for j in range(n):
                if i != j:
                    inner_product = np.inner(
                        matrix[i,:],
                        matrix[j,:]
                    )
                    norm_i = np.linalg.norm(matrix[i,:])
                    norm_j = np.linalg.norm(matrix[j,:])

                    if np.abs(inner_product - norm_j * norm_i) < 1e-5:
                        dependent_equations += [i]
        return list(set(dependent_equations))

    def _detect_residual_oscillation(self, history):
        """Detect a two-cycle of the residual norm.

        Newton bouncing between two states produces an alternating
        residual norm: the last two values repeat the two before them
        while differing from each other. Such a cycle is stable up to
        floating point noise and does not resolve on its own within any
        reasonable iteration budget.
        """
        if len(history) < 4:
            return False
        a, b, c, d = history[-1], history[-2], history[-3], history[-4]
        if a == 0 or b == 0:
            return False
        return (
            abs(a - c) < 1e-3 * abs(a)
            and abs(b - d) < 1e-3 * abs(b)
            and abs(a - b) > 1e-1 * max(abs(a), abs(b))
        )

    def _dampen_oscillating_increments(self):
        """Halve increments for variables whose residuals changed sign since the last iteration.

        A sign change means the Newton step overshot the zero of that equation.
        Halving the dependent columns' increments keeps the bracket from growing
        and makes subsequent iterations behave like bisection for those variables.
        """
        if self._prev_residual is None:
            return
        sign_flips = np.where(self._prev_residual * self.residual < 0)[0]
        for row in sign_flips:
            cols = np.nonzero(self.jacobian[row, :])[0]
            self.increment[cols] *= 0.5

    def _update_variables(self, variables=None):
        # cast dtype to float from numpy float64
        # this is necessary to keep the doctests running and note make them
        # look ugly all over the place
        # I have yet to come up with a better idea, or vectorize all operations
        # which requires major changes in tespy
        increment = [float(val) for val in self.increment]
        # the J_cols here point to actual variables, no need to call to
        # get_J_col yet
        relax = 1
        if self.robust_relax:
            relax = 0.05 + 0.95 * min(1, self.iter / (0.25 * self.max_iter))

        for key, data in self.variables_dict.items():
            if variables is not None and key not in variables:
                continue
            if data["variable"] in ["m", "h", "E"]:
                container = data["obj"]
                container._val_SI += increment[container.J_col] * relax
            elif data["variable"] in ["p"]:
                container = data["obj"]
                p_relax = max(
                    1, -2 * increment[container.J_col] / container.val_SI
                )
                container._val_SI += increment[container.J_col] / p_relax
            elif data["variable"] == "fluid":
                container = data["obj"]
                inc = increment[container.J_col[data["fluid"]]]
                value = container.val[data["fluid"]]
                if value > ERR:
                    f_relax = max(1, -2 * inc / value)
                else:
                    f_relax = 1

                container.val[data["fluid"]] += inc / f_relax

                if container.val[data["fluid"]] < ERR:
                    container.val[data["fluid"]] = 0
                elif container.val[data["fluid"]] > 1 - ERR:
                    container.val[data["fluid"]] = 1
            else:
                # add increment
                data["obj"]._val_SI += increment[data["obj"].J_col] * relax

                # keep value within specified value range
                if data["obj"].val_SI < data["obj"].min_val:
                    data["obj"].val_SI = data["obj"].min_val
                elif data["obj"].val_SI > data["obj"].max_val:
                    data["obj"].val_SI = data["obj"].max_val

    def _variable_values(self):
        values = np.empty(self.variable_counter)
        for key, data in self.variables_dict.items():
            values[key] = hlp.get_variable_value(data)
        return values

    def _adapt_to_variable_bounds(self, connections=None, components=None):
        """Apply value bounds and convergence heuristics.

        Returns whether any variable value was modified. As long as the
        heuristics interfere, residual and increment do not describe the
        state the next iteration will see, so convergence must not be
        accepted.
        """
        if connections is None:
            connections = self.network.conns['object']
        if components is None:
            components = self.network.comps['object']

        values_before = self._variable_values()

        # this could be in a different place, its kind of in between
        # network and connection
        if self.iter < 10:
            for data in self.variables_dict.values():
                if type(data["obj"]) == dc_vecvar:
                    total_mass_fractions = sum(data["obj"].val.values())
                    # renormalizing an already normalized vector would
                    # rewrite the fractions in their last bits and flag a
                    # modification, vetoing the convergence acceptance
                    if abs(total_mass_fractions - 1) < 1e-12:
                        continue
                    for fluid in data["obj"].is_var:
                        data["obj"]._val[fluid] /= total_mass_fractions

        # gate the heuristics on the scaled increment: an absolute norm
        # over mixed units opens the gate for meaningless steps, e.g.
        # sub-watt energy flow increments on a megawatt scale plant at a
        # warm started solution
        increment_moving = (
            np.abs(self.increment) / self.variable_weights
        ).max() > 1e-6

        if increment_moving:
            self._apply_property_bounds(connections)
            for c in connections:
                c._adjust_to_property_limits(self.network)

            for cp in components:
                cp._adjust_to_property_limits()

        # second check based on component heuristics
        # - for first three iterations
        # - only if the increment is sufficiently large
        # - only in design case
        if (
                self.iter < 3
                and increment_moving
                and self.network.mode == "design"
            ):
            for cp in components:
                cp.convergence_check()

            self._apply_property_bounds(connections)
            for c in connections:
                c._adjust_to_property_limits(self.network)
        # detect if any heuristic modifications were applied. Changes
        # below the solver's resolution do not count: e.g. the mass
        # fraction renormalization and the converged solution can disagree
        # by less than the acceptance tolerance indefinitely - such a
        # rewrite must not veto the convergence acceptance
        if self.variable_counter == 0:
            return False
        relative_change = (
            np.abs(self._variable_values() - values_before)
            / self.variable_weights
        )
        return bool(relative_change.max() > SCALED_RESIDUAL_TOLERANCE)

    def _apply_property_bounds(self, connections):
        """Clamp variables to the tightest bounds of their represented ones.

        Every connection provides the bounds of its own m, p and h in its
        value space. The bounds are mapped through the affine relation into
        the space of the reference variable and intersected, so a variable
        representing several linearly dependent ones is adjusted exactly
        once instead of every connection clamping the shared reference in
        turn.
        """
        connections = set(connections)
        # the enthalpy bounds depend on pressure, so pressures are adjusted
        # first
        for props in (("m", "p"), ("h",)):
            for key, data in self.variables_dict.items():
                if data["variable"] not in props or not data["obj"].is_var:
                    continue

                lower, upper = -np.inf, np.inf
                labels = []
                for var in data["_represents"]:
                    member = self._variable_lookup[var]
                    obj = member["object"]
                    if obj not in connections:
                        continue

                    bounds = obj._property_bounds(
                        member["property"], self.network
                    )
                    if bounds is None:
                        continue

                    labels.append(obj.label)
                    container = obj.get_attr(member["property"])
                    f = container._factor
                    lo, hi = (
                        (b - container._offset) / f if b is not None else None
                        for b in bounds
                    )
                    if f < 0:
                        lo, hi = hi, lo
                    if lo is not None and lo > lower:
                        lower = lo
                    if hi is not None and hi < upper:
                        upper = hi

                if lower > upper:
                    msg = (
                        "The value bounds of variable "
                        f"{self.format_var_label(key)} mapped from its "
                        f"represented variables ({', '.join(labels)}) are "
                        "incompatible, skipping the value adjustment."
                    )
                    logger.warning(msg)
                    continue

                if data["obj"]._val_SI < lower:
                    data["obj"]._val_SI = lower
                elif data["obj"]._val_SI > upper:
                    data["obj"]._val_SI = upper
                else:
                    continue
                msg = (
                    f"Value of variable {self.format_var_label(key)} out of "
                    f"bounds, adjusted to {data['obj']._val_SI}."
                )
                logger.debug(msg)

    def _solve_iteration(self):
        r"""
        Control iteration step of the newton algorithm.

        - Calculate the residual value for each equation
        - Calculate the jacobian matrix
        - Calculate new values for variables
        - Restrict fluid properties to value ranges
        - Check component parameters for consistency
        """
        self._solve_equations()
        self._update_residual_scales()
        self._invert_jacobian()

        if self.lin_dep:
            self._diagnose_singularity()
            return True

        if self.oscillation_damping or self._auto_damping:
            self._dampen_oscillating_increments()

        self._update_variables()
        modified = self._adapt_to_variable_bounds()
        self._prev_residual = self.residual.copy()
        return modified

    def _update_residual_scales(self):
        """Per equation scales from the current jacobian.

        The scale of an equation is its response to order one relative
        changes of the variables it depends on, so the scaled residual
        :code:`|r| / s` reads as the relative change of the variables
        required to close the equation. This makes the convergence
        acceptance independent of the physical magnitude of the
        equations - a heat exchanger energy balance on a megawatt scale
        plant would otherwise have to close to fractions of a watt,
        below the noise of the fluid property evaluations chained into
        its residual.
        """
        values = np.abs(self._variable_values())
        self.variable_weights = np.maximum(values, 1.0)
        self.residual_scales = np.maximum(
            np.abs(self.jacobian) @ self.variable_weights, 1.0
        )

    def scaled_residual(self, rows=None):
        """Residuals divided by their per equation scales."""
        if rows is None:
            return np.abs(self.residual) / self.residual_scales
        return np.abs(self.residual[rows]) / self.residual_scales[rows]

    def _solve_equations(self, objects=None, residual_only=False):
        r"""
        Calculate the residual and derivatives of all equations.
        """
        if objects is None:
            to_solve = (
                self.network.comps["object"].tolist()
                + self.network.conns["object"].tolist()
                + list(self.network.user_defined_eq.values())
            )
        else:
            to_solve = objects

        if residual_only:
            for obj in to_solve:
                hlp.solve_residuals(obj)
                if len(obj.residual) > 0:
                    rows = list(obj.residual.keys())
                    self.residual[rows] = list(obj.residual.values())
                obj.it += 1
            return

        for obj in to_solve:
            hlp.solve(obj, self.increment_filter)
            if len(obj.jacobian) > 0:
                rows = list(obj.residual.keys())
                data = list(obj.residual.values())
                self.residual[rows] = data

                rows = [k[0] for k in obj.jacobian]
                columns = [k[1] for k in obj.jacobian]
                data = list(obj.jacobian.values())
                self.jacobian[rows, columns] = data

            obj.it += 1

    def format_var_label(self, v_idx):
        v_data = self.variables_dict[v_idx]
        v_type = v_data["variable"]
        if v_type == "fluid" and v_data["fluid"] is not None:
            return f"{v_data['fluid']}{v_idx}"
        return f"{v_type}{v_idx}"

    @staticmethod
    def format_eq_name(eq_name):
        if isinstance(eq_name, tuple):
            name, sub_idx = eq_name
            return f"{name}{{{sub_idx}}}" if sub_idx > 0 else name
        return eq_name

    def _get_equations_by_number(self, number_list) -> dict:
        """Get the actual equations after presolving the problem by equation
        number

        Returns
        -------
        dict
            Lookup with equation number as index and tuple of label and
            parameter defining the equation. In case one parameter defines
            multiple equations, the same equation is repeated.
        """
        return {
            k: v for k, v in self._equation_lookup.items() if k in number_list
        }

    def get_sorted_residual_index(self) -> list[int]:
        """Get the equation indices sorted by scaled residual.

        Returns
        -------
        list[int]
            List of equation numbers, the index values.
        """
        return list(np.argsort(self.scaled_residual())[::-1])

    def _worst_equation_hint(self):
        """Name the equation with the largest scaled residual."""
        if self.residual is None or len(self.residual) == 0:
            return ""
        worst = int(np.argmax(self.scaled_residual()))
        equations = self._get_equations_by_number([worst])
        for obj, name in equations.values():
            label = obj.label if hasattr(obj, "label") else str(obj)
            return f" ({label}: {self.format_eq_name(name)})"
        return ""

    def get_linear_dependent_variables(self) -> list:
        """Get the sets of linearly dependent variables by label."""
        variable_list = []
        for dependents in self._variable_dependencies:
            variables = [
                self._variable_lookup[v] for v in dependents["variables"]
            ]
            variable_list += [
                [(v["object"].label, v["property"]) for v in variables]
            ]
        return variable_list

    def get_presolved_equations(self) -> list:
        """Get the equations consumed by presolving by label."""
        return [
            v for k, v in self._equation_set_lookup.items()
            if k in self._presolved_equations
        ]

    def get_variables_before_presolve(self) -> list:
        """Get the original variables by label."""
        return [
            (v["object"].label, v["property"])
            for v in self._variable_lookup.values()
        ]

    def get_presolved_variables(self) -> list:
        """Get the variables solved by presolving by label."""
        represented_variables = []
        for v in self.variables_dict.values():
            represented_variables += v["_represents"]
        if len(self.variables_dict) == 0 and len(self._presolved_equations) == 0:
            return []
        return [
            (v["object"].label, v["property"])
            for key, v in self._variable_lookup.items()
            if key not in represented_variables
        ]

    def get_variables(self) -> dict:
        """Get the variables with the original variables they represent."""
        return {
            (key, data["variable"]):
            [
                (
                    self._variable_lookup[v]["object"].label,
                    self._variable_lookup[v]["property"]
                ) for v in data["_represents"]
            ]
            for key, data in self.variables_dict.items()
        }

    def get_variable_values(self, block=None) -> dict:
        """Get the variables with the current values of the original
        variables they represent, optionally restricted to the variables
        of the block with the given id."""
        if block is not None:
            blocks = {b.id: b for b in self.decompose().blocks}
            if block not in blocks:
                msg = f"There is no block with id {block}."
                raise KeyError(msg)
            columns = set(blocks[block].variables)

        result = {}
        for key, data in self.variables_dict.items():
            if block is not None and key not in columns:
                continue
            represents = []
            for v in data["_represents"]:
                lookup = self._variable_lookup[v]
                container = lookup["object"].get_attr(lookup["property"])
                if data["variable"] == "fluid":
                    value = container.val[data["fluid"]]
                else:
                    value = container.val_SI
                represents.append(
                    (lookup["object"].label, lookup["property"], value)
                )
            result[(key, data["variable"])] = represents
        return result

    def get_block_jacobian(self, block) -> dict:
        """Get the linear system of a failed block at its state of failure.

        The jacobian and residual restricted to the block's equations and
        variables are recorded as evaluated in the failing iteration of the
        last solution attempt (surrogate placeholders for all-zero rows
        appear as their fill-in value of 1). They are only available for
        blocks whose solution failed. The result also holds the expected
        entries declared by the incidence matrix, so evaluated entries can
        be compared against the declared dependencies, e.g. to find
        derivatives that are zero although the equation declares the
        dependency.
        """
        blocks = {b.id: b for b in self.decompose().blocks}
        if block not in blocks:
            msg = f"There is no block with id {block}."
            raise KeyError(msg)

        data = blocks[block]
        if data.jacobian is None:
            msg = (
                f"There is no recorded jacobian for block {block}: it is "
                "only recorded in case the solution of the block fails."
            )
            raise hlp.TESPyNetworkError(msg)

        expected = np.zeros(data.jacobian.shape, dtype=bool)
        for i, eq in enumerate(data.equations):
            declared = set(self._incidence_matrix.get(eq, []))
            for j, col in enumerate(data.variables):
                expected[i, j] = col in declared

        return {
            "equations": [
                f"{lbl}.{self.format_eq_name(name)}"
                for lbl, name in data.equation_labels
            ],
            "variables": list(data.variable_labels),
            "jacobian": data.jacobian.copy(),
            "residual": data.residual.copy(),
            "expected": expected,
            "values": list(data.variable_values),
        }

    def get_block_states(self, block, at="current") -> list:
        """Get the states of the connections a block touches.

        The properties of every connection involved in the block's
        equations and variables. What is reported comes from the
        connection classes, e.g. mass flow, pressure, enthalpy, temperature
        and phase for fluid connections. With :code:`at="current"` the
        states are evaluated from the current variable values - while the
        solution process is paused at the block this is its entry state,
        including any modification made in between. :code:`at="failure"`
        returns the states as recorded in the failing iteration of the
        block's last solution attempt, only available for blocks whose
        solution failed.
        """
        blocks = {b.id: b for b in self.decompose().blocks}
        if block not in blocks:
            msg = f"There is no block with id {block}."
            raise KeyError(msg)

        if at == "current":
            return self._connection_states(blocks[block])
        elif at == "failure":
            data = blocks[block]
            if data.connection_states is None:
                msg = (
                    f"There are no recorded states for block {block}: they "
                    "are only recorded in case the solution of the block "
                    "fails."
                )
                raise hlp.TESPyNetworkError(msg)
            return [dict(state) for state in data.connection_states]
        else:
            msg = "The at parameter must either be 'current' or 'failure'."
            raise ValueError(msg)

    def _connection_states(self, block):
        """Collect the state of every connection the block touches."""
        columns = set(block.variables)
        objects = [self._equation_obj_lookup[eq] for eq in block.equations]
        for col in block.variables:
            for sm_col in self.variables_dict[col]["_represents"]:
                objects.append(self._variable_lookup[sm_col]["object"])

        states = []
        seen = set()
        for obj in objects:
            if isinstance(obj, ConnectionBase):
                connections = [obj]
            else:
                connections = (
                    getattr(obj, "inl", []) + getattr(obj, "outl", [])
                )
            for c in connections:
                if not isinstance(c, ConnectionBase) or c.label in seen:
                    continue
                seen.add(c.label)
                state = {"label": c.label, "block_variables": []}
                for prop, value, container in c._debug_state():
                    state[prop] = value
                    if container is None:
                        continue
                    try:
                        if container.is_var and container.J_col in columns:
                            state["block_variables"].append(prop)
                    except ValueError:
                        pass
                states.append(state)
        return states

    def set_variable_value(self, obj, prop, value):
        """Set the SI value of a variable through any of the linearly
        dependent variables it represents.

        Raises a :code:`KeyError` in case the object does not hold a
        variable of the given name and a :code:`TESPyNetworkError` in case
        the property is not part of the variable space.
        """
        lookup = self._object_to_variable_lookup
        if obj not in lookup or prop not in lookup[obj]:
            msg = f"The object {obj.label} does not have a variable {prop}."
            raise KeyError(msg)

        if self._paused_driver is not None:
            block = self._paused_driver.paused_block
            allowed = {
                sm_col for col in block.variables
                for sm_col in self.variables_dict[col]["_represents"]
            }
            if lookup[obj][prop] not in allowed:
                var_str = ", ".join(block.variable_labels)
                msg = (
                    f"The solution process is paused at block {block.id}: "
                    "only the variables of that block can be modified "
                    f"({var_str})."
                )
                raise hlp.TESPyNetworkError(msg)

        container = obj.get_attr(prop)
        if prop == "fluid":
            msg = (
                "Setting fluid mass fractions through set_variable_value is "
                "not supported."
            )
            raise hlp.TESPyNetworkError(msg)
        if not container.is_var:
            msg = (
                f"{obj.label} ({prop}) is not part of the variable space, "
                "its value is fixed by specification or presolving."
            )
            raise hlp.TESPyNetworkError(msg)

        container.set_reference_val_SI(value)

    def get_linear_dependents_by_object(self, obj, prop) -> list:
        """Get the linearly dependent variables of a variable by label.

        Raises a :code:`KeyError` in case the object does not hold a
        variable of the given name.
        """
        if obj not in self._object_to_variable_lookup:
            msg = f"The object {obj.label} does not have any variables."
            raise KeyError(msg)

        if prop not in self._object_to_variable_lookup[obj]:
            msg = f"The object {obj.label} does not have a variable {prop}."
            raise KeyError(msg)

        variable_idx = self._object_to_variable_lookup[obj][prop]
        for dependents in self._variable_dependencies:
            if variable_idx in dependents["variables"]:
                variables = [
                    self._variable_lookup[v] for v in dependents["variables"]
                ]
                return [
                    (v["object"].label, v["property"]) for v in variables
                ]
        msg = (
            f"Variable index {variable_idx} not found in any dependency "
            "group."
        )
        raise KeyError(msg)

    def _get_variables_by_number(self, number_list) -> dict:
        return {
            (key, data["variable"]):
            [
                (
                    self._variable_lookup[v]["object"].label,
                    self._variable_lookup[v]["property"]
                ) for v in data["_represents"]
            ]
            for key, data in self.variables_dict.items()
            if key in number_list
        }

    def get_equations(self) -> dict:
        """Get the lookup of solver equation numbers to labels."""
        return self._equation_lookup

    def get_incidence(self) -> dict:
        """Get the lookup of solver equation numbers to the variable
        numbers each equation depends on."""
        return self._incidence_matrix

    def get_equations_with_dependents(self) -> dict:
        """Get the equations with the variables they depend on."""
        dependencies = {}
        for eq_idx, dependents in self._incidence_matrix.items():
            dependencies.update({
                self._equation_lookup[eq_idx]:
                list(self._get_variables_by_number(dependents).keys())
            })
        return dependencies

    def structure_variables(self) -> list:
        """Get the serializable description of the structural variables."""
        affine_map = {}
        for dependents in self._variable_dependencies:
            for var in dependents["variables"]:
                affine_map[var] = (
                    dependents["reference"],
                    dependents["factors"][var],
                    dependents["offsets"][var]
                )

        variables = []
        for col, data in sorted(self._variable_lookup.items()):
            obj = data["object"]
            prop = data["property"]
            container = obj.get_attr(prop)
            reference, factor, offset = affine_map.get(col, (col, 1.0, 0.0))

            if prop == "fluid":
                if len(container.is_var) > 0:
                    state = "variable"
                    solver_index = {
                        fluid: int(j_col)
                        for fluid, j_col in container.J_col.items()
                    }
                elif len(container._is_set) > 0:
                    state = "specified"
                    solver_index = None
                else:
                    state = "presolved"
                    solver_index = None
            elif container.is_set:
                state = "specified"
                solver_index = None
            elif container.is_var:
                state = "variable"
                solver_index = int(container.J_col)
            else:
                state = "presolved"
                solver_index = None

            variables.append({
                "id": int(col),
                "object": obj.label,
                "object_type": (
                    "connection" if isinstance(obj, ConnectionBase)
                    else "component"
                ),
                "class_name": type(obj).__name__,
                "property": prop,
                "quantity": getattr(container, "quantity", None),
                "state": state,
                "reference": int(reference),
                "factor": float(factor),
                "offset": float(offset),
                "solver_index": solver_index,
            })
        return variables

    def structure_equations(self) -> list:
        """Get the serializable description of the structural equations."""
        graph = self.structure_graph
        edge_by_row = {row: edge for edge, row in graph.edge_eq_idx.items()}
        solver_equations = {}
        for eq_num, (label, eq_name) in self._equation_lookup.items():
            name = eq_name[0] if isinstance(eq_name, tuple) else eq_name
            solver_equations.setdefault((label, name), []).append(int(eq_num))

        equations = []
        for row, (label, name) in sorted(self._equation_set_lookup.items()):
            if row in edge_by_row:
                kind = "affine"
                related = [int(col) for col in edge_by_row[row]]
            elif row in graph.linear_rows:
                kind = "linear"
                related = [int(col) for col in graph.linear_rows[row]]
            else:
                kind = "nonlinear"
                related = []

            state = (
                "consumed" if row in self._presolved_equations else "active"
            )
            indices = (
                solver_equations.get((label, name), [])
                if state == "active" else []
            )
            if kind == "nonlinear" and indices:
                related = sorted({
                    int(col)
                    for eq_num in indices
                    for col in graph.nonlinear_incidence.get(eq_num, [])
                })

            equations.append({
                "id": int(row),
                "object": label,
                "name": name,
                "kind": kind,
                "origin": self._equation_set_origin[row],
                "state": state,
                "variables": related,
                "solver_indices": indices,
            })
        return equations
