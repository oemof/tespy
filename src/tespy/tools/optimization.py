# -*- coding: utf-8

"""Module for OptimizationProblem class.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/tools/optimization.py

SPDX-License-Identifier: MIT
"""
import warnings

import numpy as np

try:
    # this is to make this import of tespy possible without the optional
    # dependency pymoo available
    from pymoo.core.problem import ElementwiseProblem

except ModuleNotFoundError:
    class ElementwiseProblem:

        def __init__(self, **kwargs):
            pass

from tespy.tools.logger import logger


class _NestedModelAdapter:
    """Wraps a nested-API model to present the flat get_parameter/solve_model interface."""

    def __init__(self, model, param_mapping):
        self._model = model
        self._mapping = param_mapping

    def solve_model(self, **flat_kwargs):
        nested = {}
        for flat_name, value in flat_kwargs.items():
            obj, label, param = self._mapping[flat_name]
            if param is not None:
                nested.setdefault(obj, {}).setdefault(label, {})[param] = value
            else:
                nested.setdefault(obj, {})[label] = value
        self._model.solve_model(**nested)

    def get_parameter(self, name):
        if name in self._mapping:
            obj, label, param = self._mapping[name]
            return self._model.get_param(obj, label, param)
        return self._model.get_objective(name)

    def get_objectives(self, objective_list):
        return self._model.get_objectives(objective_list)

    def penalize(self, fitness, c):
        return self._model.penalize(fitness, c)


def _is_nested(variables, constraints, kpi):
    if variables and not any(
        isinstance(v, dict) and ("min" in v or "max" in v)
        for v in variables.values()
    ):
        return True
    if "lower limits" in constraints or "upper limits" in constraints:
        return True
    if isinstance(kpi, dict):
        return True
    return False


def _translate_nested(variables, constraints, kpi):
    """Translate legacy nested dicts to flat equivalents.

    Returns (flat_vars, flat_constraints, flat_kpi, param_mapping).
    param_mapping: flat_name -> (obj, label, param)  - param is None for Customs.
    """
    param_mapping = {}

    flat_vars = {}
    for obj, data in variables.items():
        if obj in ("Connections", "Components"):
            for label, params in data.items():
                for param, bounds in params.items():
                    flat_name = f"{obj}-{label}-{param}"
                    flat_vars[flat_name] = bounds
                    param_mapping[flat_name] = (obj, label, param)
        else:
            for label, bounds in data.items():
                flat_name = f"{obj}-{label}"
                flat_vars[flat_name] = bounds
                param_mapping[flat_name] = (obj, label, None)

    ref_map = {}
    for key, value in constraints.items():
        if key not in ("lower limits", "upper limits") and isinstance(value, list):
            obj, label, param = value
            ref_flat = f"{obj}-{label}-{param}"
            ref_map[key] = ref_flat
            param_mapping[ref_flat] = (obj, label, param)

    flat_constraints = {}
    for border, direction in [("lower", "min"), ("upper", "max")]:
        limit_key = f"{border} limits"
        if limit_key not in constraints:
            continue
        for obj, data in constraints[limit_key].items():
            if obj in ("Connections", "Components"):
                for label, params in data.items():
                    for param, val in params.items():
                        flat_name = f"{obj}-{label}-{param}"
                        param_mapping.setdefault(flat_name, (obj, label, param))
                        bound = ref_map[val] if isinstance(val, str) and val in ref_map else val
                        flat_constraints.setdefault(flat_name, {})[direction] = bound
            else:
                for label, val in data.items():
                    flat_name = f"{obj}-{label}"
                    param_mapping.setdefault(flat_name, (obj, label, None))
                    bound = ref_map[val] if isinstance(val, str) and val in ref_map else val
                    flat_constraints.setdefault(flat_name, {})[direction] = bound

    flat_kpi = []
    if isinstance(kpi, dict):
        for obj, data in kpi.items():
            if obj in ("Connections", "Components"):
                for label, params in data.items():
                    for param in params:
                        flat_name = f"{obj}-{label}-{param}"
                        flat_kpi.append(flat_name)
                        param_mapping.setdefault(flat_name, (obj, label, param))
            else:
                for label in data:
                    flat_name = f"{obj}-{label}"
                    flat_kpi.append(flat_name)
                    param_mapping.setdefault(flat_name, (obj, label, None))
    else:
        flat_kpi = list(kpi)

    return flat_vars, flat_constraints, flat_kpi, param_mapping


class OptimizationProblem(ElementwiseProblem):
    r"""
    The OptimizationProblem handles the optimization.

    Parameters
    ----------
    model : tespy.models.ModelTemplate
        Model instance providing :code:`set_parameters`, :code:`get_parameter`,
        :code:`get_objectives` and :code:`solve_model`.

    variables : dict
        Flat dictionary of decision variables and their bounds, e.g.

        .. code-block:: python

            {
                "extraction pressure 1": {"min": 1, "max": 40},
                "extraction pressure 2": {"min": 1, "max": 40},
            }

    constraints : dict
        Flat dictionary of parameter constraints. Values are numeric bounds or
        a parameter name string for cross-parameter references, e.g.

        .. code-block:: python

            {
                "extraction pressure 1": {"min": "extraction pressure 2"},
                "some temperature": {"min": 20, "max": 200},
            }

    objective : list
        Names of the objective parameters as defined in the model's
        :code:`_parameter_lookup`.

    minimize : list
        :code:`True` to minimize, :code:`False` to maximize. One entry per
        objective, in the same order.

    kpi : list
        Parameter names to log at each evaluation in addition to the
        objectives, e.g. :code:`["hpt power", "hpt pressure ratio"]`.

    penalty_instead_of_constraints : bool
        If :code:`True`, constraints are passed to :code:`model.penalize`
        instead of being enforced directly. Default :code:`False`.

    Example
    -------
    For an example please check out
    :ref:`this section <tutorial_optimization_label>` in the docs.
    """

    def __init__(
        self,
        model,
        variables=None,
        constraints=None,
        objective=None,
        minimize=None,
        kpi=None,
        penalty_instead_of_constraints=False,
    ):
        if variables is None:
            variables = {}
        if constraints is None:
            constraints = {}
        if objective is None:
            objective = []
        if kpi is None:
            kpi = []

        if _is_nested(variables, constraints, kpi):
            warnings.warn(
                "Passing nested dictionaries to OptimizationProblem is "
                "deprecated and will be removed in a future release. Use flat "
                "parameter name dictionaries instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            variables, constraints, kpi, param_mapping = _translate_nested(variables, constraints, kpi)
            model = _NestedModelAdapter(model, param_mapping)

        self.model = model
        self.variables = variables
        self.constraints = constraints
        self.kpi = kpi

        self._build_objective(objective, minimize)
        self._build_variables()
        self._build_kpi()

        self.nic = 0
        self.constraint_list = []
        self._build_constraints()

        self.penalty_instead_of_constraints = penalty_instead_of_constraints
        n_ieq_constr = 0 if self.penalty_instead_of_constraints else len(self.constraint_list)

        self.log = []

        super().__init__(
            n_var=len(self.variable_list),
            n_obj=len(self.objective_list),
            n_ieq_constr=n_ieq_constr,
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
        for param_name, bounds in self.variables.items():
            self._bounds[0].append(bounds["min"])
            self._bounds[1].append(bounds["max"])
            self.variable_list.append(param_name)

    def _build_kpi(self) -> None:
        self.kpi_list = list(self.kpi)

    def _build_constraints(self) -> None:
        for param_name, bounds in self.constraints.items():
            if "min" in bounds:
                right = bounds["min"] if isinstance(bounds["min"], str) else str(bounds["min"])
                self.constraint_list.append(f"{param_name}>={right}")
                self.nic += 1
            if "max" in bounds:
                right = bounds["max"] if isinstance(bounds["max"], str) else str(bounds["max"])
                self.constraint_list.append(f"{param_name}<={right}")
                self.nic += 1

    def _evaluate_constraints(self) -> list:
        evaluation = []
        for param_name, bounds in self.constraints.items():
            actual = self.model.get_parameter(param_name)
            if "min" in bounds:
                val = bounds["min"]
                limit = self.model.get_parameter(val) if isinstance(val, str) else val
                evaluation.append(limit - actual)
            if "max" in bounds:
                val = bounds["max"]
                limit = self.model.get_parameter(val) if isinstance(val, str) else val
                evaluation.append(actual - limit)
        return evaluation

    def _evaluate(self, x: np.ndarray, out: dict, *_args, **_kwargs) -> None:
        self.model.solve_model(**dict(zip(self.variable_list, x)))
        fitness = self.model.get_objectives(self.objective_list)
        kpi = [self.model.get_parameter(k) for k in self.kpi_list]

        _fitness = [(-1) ** (sense + 1) * f for f, sense in zip(fitness, self.minimize)]
        c = self._evaluate_constraints()

        if self.penalty_instead_of_constraints:
            _fitness = self.model.penalize(_fitness, c)
        else:
            out["G"] = c

        out["F"] = [value if not np.isnan(value) else 1e21 for value in _fitness]

        log_entry = {
            **{self.variable_list[i]: val for i, val in enumerate(x)},
            **{self.objective_list[i]: val for i, val in enumerate(fitness)},
            **{self.constraint_list[i]: val for i, val in enumerate(c)},
            **{self.kpi_list[i]: val for i, val in enumerate(kpi)},
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

        cu = self._evaluate_constraints("upper")
        cl = self._evaluate_constraints("lower")

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
            "The optimization using pygmo is deprecated and will be  removed "
            "in the next major release, please use pymoo in the future. To "
            "get pymoo you can install tespy with extra dependency 'opt': pip "
            "install tespy[opt]."
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
