# -*- coding: utf-8

"""Module for OptimizationProblem class.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/tools/optimization.py

SPDX-License-Identifier: MIT
"""
import warnings

import numpy as np
import pandas as pd
from pymoo.core.problem import ElementwiseProblem

from tespy.tools.helpers import merge_dicts
from tespy.tools.logger import logger


class OptimizationProblem(ElementwiseProblem):
    r"""
    The OptimizationProblem handles the optimization.

    - Set up the optimization problems by specifying constraints, upper and
      lower bounds for the decision variables and selection of the objective
      function.
    - Run the optimization, see
      :py:meth:`tespy.tools.optimization.OptimizationProblem.run`.
    - Provide the optimization results DataFrame in the
      :code:`.individuals` attribute of the :code:`OptimizationProblem` class.

    Parameters
    ----------
    model : custom class
        Object of some class, which provides all the methods required by the
        optimization suite, see the Example section for a downloadable
        template of the implementation.

    variables : dict
        Dictionary containing the decision variables and their respective
        bounds.

    constraints : dict
        Dictionary containing the constraints for the model.

    objective : list
        Name of the objective(s). :code:`objective` is passed to the
        :code:`get_objectives` method of your tespy model instance.

    minimize : list
        Minimize (:code:`True`) or maximize an objective (:code:`False`). Must
        be passed as list in same order as objective. E.g.
        :code:`minimize=[True, False]` if the first objective should be
        minimized and the second one maximized.

    Example
    -------
    For an example please check out
    :ref:`this section <tutorial_optimization_label>` in the docs.
    """

    def __init__(self, model, variables={}, constraints={}, objective=[], minimize=None, kpi={}):
        self.model = model
        default_variables = {"Connections": {}, "Components": {}}
        default_constraints = {
            "lower limits": {"Connections": {}, "Components": {}},
            "upper limits": {"Connections": {}, "Components": {}}
        }
        default_kpi = {"Connections": {}, "Components": {}}
        # merge the passed values into the default dictionary structure
        self.variables = merge_dicts(variables, default_variables)
        self.constraints = merge_dicts(constraints, default_constraints)
        self.kpi = merge_dicts(kpi, default_kpi)

        self._build_objective(objective, minimize)
        self._build_variables()
        self._build_kpi()

        self.input_dict = self.variables.copy()

        self.nic = 0
        self.constraint_list = []
        self._build_constraints("upper")
        self._build_constraints("lower")

        self.log = []

        super().__init__(
            n_var=len(self.variable_list),
            n_obj=len(self.objective_list),
            n_ieq_constr=len(self.constraint_list),
            n_eq_constr=0,
            xl=self._bounds[0],
            xu=self._bounds[1]
        )

    def _build_objective(self, objective, minimize) -> None:
        if not isinstance(objective, list):
            msg = "The objective(s) must be passed as a list."
            raise TypeError(msg)

        self.objective_list = objective
        self.nobj = len(self.objective_list)

        if minimize is None:
            self.minimize = [True for _ in self.objective_list]
        elif len(minimize) != self.nobj:
            msg = (
                "If you supply the minimize argument the number of values in "
                "the list must be identical to the number of objectives."
            )
            raise ValueError(msg)
        else:
            self.minimize = minimize

    def _build_variables(self) -> None:
        self.variable_list = []
        self._bounds = [[], []]
        for obj, data in self.variables.items():
            for label, params in data.items():
                if obj in ["Connections", "Components"]:
                    for param in params:
                        self._bounds[0] += [
                            self.variables[obj][label][param]['min']
                        ]
                        self._bounds[1] += [
                            self.variables[obj][label][param]['max']
                        ]
                        self.variable_list += [obj + '-' + label + '-' + param]
                else:
                    self._bounds[0] += [self.variables[obj][label]['min']]
                    self._bounds[1] += [self.variables[obj][label]['max']]
                    self.variable_list += [obj + '-' + label]

    def _build_kpi(self) -> None:
        self.kpi_list = []
        for obj, data in self.kpi.items():
            for label, params in data.items():
                if obj in ["Connections", "Components"]:
                    for param in params:
                        self.kpi_list += [obj + '-' + label + '-' + param]
                else:
                    self.kpi_list += [obj + '-' + label]

    def _build_constraints(self, border: str) -> None:
        for obj, data in self.constraints[f"{border} limits"].items():
            for label, constraints in data.items():
                for param, constraint in constraints.items():
                    self.nic += 1
                    if isinstance(constraint, str):
                        right_side = "-".join(self.constraints[constraint])
                    else:
                        right_side = str(constraint)

                    direction = ">=" if border == "lower" else "<="
                    self.constraint_list += [
                        f"{obj}-{label}-{param}{direction}{right_side}"
                    ]

    def _evaluate_constraints(self, border: str):
        evaluation = []
        for obj, data in self.constraints[f"{border} limits"].items():
            for label, constraints in data.items():
                for param, constraint in constraints.items():
                    if isinstance(constraint, str):
                        # this is an internal reference to another attribute in
                        # the model
                        c = (
                            self.model.get_param(
                                *self.constraints[constraint]
                            ) - self.model.get_param(obj, label, param)
                        )
                    else:
                        # this is the constraint as an actual numerical vaue
                        c = (
                            constraint -
                            self.model.get_param(obj, label, param)
                        )
                    if border == "lower":
                        evaluation += [c]
                    else:
                        evaluation += [-c]
        return evaluation

    def _evaluate(self, x: np.ndarray, out: dict, *args, **kwargs) -> None:
        i = 0
        for obj, data in self.variables.items():
            for label, params in data.items():
                if obj in ["Connections", "Components"]:
                    for param in params:
                        self.input_dict[obj][label][param] = x[i]
                        i += 1
                else:
                    self.input_dict[obj][label] = x[i]
                    i += 1

        self.model.solve_model(**self.input_dict)
        fitness = self.model.get_objectives(self.objective_list)

        kpi = []
        for obj, data in self.kpi.items():
            for label, params in data.items():
                if obj in ["Connections", "Components"]:
                    for param in params:
                        kpi += [self.model.get_param(obj, label, param)]
                else:
                    kpi += [self.model.get_param(obj, label, None)]

        # negate the fitness function evaluation for minimize = False
        # parenthesis around the -1 are required!
        out["F"] = [
            (-1) ** (sense + 1) * f for f, sense in zip(fitness, self.minimize)
        ]

        cu = self._evaluate_constraints("upper")
        cl = self._evaluate_constraints("lower")

        out["G"] = cu + cl

        log_entry = {
            **{self.variable_list[i]: val for i, val in enumerate(x)},
            **{self.objective_list[i]: val for i, val in enumerate(fitness)},
            **{self.constraint_list[i]: val for i, val in enumerate(cu + cl)},
            **{self.kpi_list[i]: val for i, val in enumerate(kpi)}
        }
        self.log.append(log_entry)

    def fitness(self, x):
        """Evaluate the fitness function of an individual.

        Parameters
        ----------
        x : list
            List of the decision variables' values of the current individual.

        Returns
        -------
        fitness : list
            A list containing the fitness function evaluation as well as the
            evaluation of the upper and lower constraints.
        """
        i = 0
        for obj, data in self.variables.items():
            for label, params in data.items():
                if obj in ["Connections", "Components"]:
                    for param in params:
                        self.input_dict[obj][label][param] = x[i]
                        i += 1
                else:
                    self.input_dict[obj][label] = x[i]
                    i += 1

        self.model.solve_model(**self.input_dict)
        fitness = self.model.get_objectives(self.objective_list)

        # negate the fitness function evaluation for minimize = False
        # parenthesis around the -1 are required!
        fitness = [
            (-1) ** (sense + 1) * f for f, sense in zip(fitness, self.minimize)
        ]

        cu = self.collect_constraints("upper")
        cl = self.collect_constraints("lower")

        return fitness + cu + cl

    def get_nobj(self):
        """Return number of objectives."""
        return self.nobj

    # inequality constraints (equality constraints not required)
    def get_nic(self):
        """Return number of inequality constraints."""
        return self.nic

    def get_bounds(self):
        """Return bounds of decision variables."""
        return self._bounds

    def _process_generation_data(self, evo, pop):
        """Process the data of the individuals within one evolution.

        Parameters
        ----------
        evo : int
            Evolution number.

        pop : pygmo.population
            PyGMO population object.
        """
        for individual, (x, obj) in enumerate(zip(pop.get_x(), pop.get_f())):
            self.individuals.loc[(evo, individual), self.variable_list] = x
            self.individuals.loc[
                (evo, individual),
                self.objective_list + self.constraint_list
            ] = obj

    def run(self, algo, pop, num_ind, num_evo):
        """Run the optimization algorithm.

        Parameters
        ----------
        algo : pygmo.core.algorithm
            PyGMO optimization algorithm.

        pop : pygmo.core.population
            PyGMO population.

        num_ind : int
            Number of individuals.

        num_evo : int
            Number of evolutions.
        """
        msg = (
            "The pygmo API is deprecated and will be removed in the next "
            "major release. Please use the pymoo API instead in the future."
        )
        logger.warning(msg)
        warnings.warn(msg, FutureWarning)

        self.individuals = pd.DataFrame(index=range(num_evo * num_ind))

        self.individuals["evo"] = [
            evo for evo in range(num_evo) for _ in range(num_ind)
        ]
        self.individuals["ind"] = [
            ind for _ in range(num_evo) for ind in range(num_ind)
        ]

        self.individuals.set_index(["evo", "ind"], inplace=True)

        evo = 0
        for evo in range(num_evo - 1):
            self._process_generation_data(evo, pop)
            pop = algo.evolve(pop)

        if num_evo > 1:
            evo += 1

        self._process_generation_data(evo, pop)
        for obj, sense in zip(self.objective_list, self.minimize):
            self.individuals[obj] = (
                self.individuals[obj] * (-1) ** (sense + 1)
            )
        return pop
