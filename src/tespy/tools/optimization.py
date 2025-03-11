try:
    import pygmo as pg
except ImportError:
    pg = None

import pandas as pd

from tespy.tools.helpers import merge_dicts


class OptimizationProblem:
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

    Note
    ----
    For the required structure of the input dictionaries see the example in
    below.

    Installation of pygmo via pip is not available for Windows and OSX users
    currently. Please use conda instead or refer to their
    `documentation <https://esa.github.io/pygmo2/>`_.

    Example
    -------
    For an example please go to the tutorials section of TESPy's online
    documentation.
    """

    def __init__(self, model, variables={}, constraints={}, objective=[]):
        if pg is None:
            msg = (
                "For this function of TESPy pygmo has to be installed. Either"
                " use pip (Linux users only) or conda to install the latest"
                " pygmo version."
            )
            raise ImportError(msg)

        self.model = model
        default_variables = {"Connections": {}, "Components": {}}
        default_constraints = {
            "lower limits": {"Connections": {}, "Components": {}},
            "upper limits": {"Connections": {}, "Components": {}}
        }
        # merge the passed values into the default dictionary structure
        self.variables = merge_dicts(variables, default_variables)
        self.constraints = merge_dicts(constraints, default_constraints)
        self.variable_list = []
        self.constraint_list = []

        if not isinstance(objective, list):
            msg = "The objective(s) must be passed as a list."
            raise TypeError(msg)

        self.objective_list = objective
        self.nobj = len(self.objective_list)

        self.bounds = [[], []]
        for obj, data in self.variables.items():
            for label, params in data.items():
                if obj in ["Connections", "Components"]:
                    for param in params:
                        self.bounds[0] += [
                            self.variables[obj][label][param]['min']
                        ]
                        self.bounds[1] += [
                            self.variables[obj][label][param]['max']
                        ]
                        self.variable_list += [obj + '-' + label + '-' + param]
                else:
                    self.bounds[0] += [self.variables[obj][label]['min']]
                    self.bounds[1] += [self.variables[obj][label]['max']]
                    self.variable_list += [obj + '-' + label]

        self.input_dict = self.variables.copy()

        self.nic = 0
        self.collect_constraints("upper", build=True)
        self.collect_constraints("lower", build=True)

    def collect_constraints(self, border, build=False):
        """Collect the constraints

        Parameters
        ----------
        border : str
            "upper" or "lower", determine which constraints to collect.
        build : bool, optional
            If True, the constraints are evaluated and returned, by default
            False

        Returns
        -------
        tuple
            Return the upper and lower constraints evaluation lists.
        """
        evaluation = []
        for obj, data in self.constraints[f'{border} limits'].items():
            for label, constraints in data.items():
                for param, constraint in constraints.items():
                    # to build the equations
                    if build:
                        self.nic += 1
                        if isinstance(constraint, str):
                            right_side = '-'.join(self.constraints[constraint])
                        else:
                            right_side = str(constraint)

                        direction = '>=' if border == 'lower' else '<='
                        self.constraint_list += [
                            obj + '-' + label + '-' + param + direction +
                            right_side
                        ]
                    # to get the constraints evaluation
                    else:
                        if isinstance(constraint, str):
                            c = (
                                self.model.get_param(
                                    *self.constraints[constraint]
                                ) - self.model.get_param(obj, label, param)
                            )
                        else:
                            c = (
                                constraint -
                                self.model.get_param(obj, label, param)
                            )
                        if border == 'lower':
                            evaluation += [c]
                        else:
                            evaluation += [-c]

        if build:
            return None
        else:
            return evaluation

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
        f1 = self.model.get_objectives(self.objective_list)

        cu = self.collect_constraints("upper")
        cl = self.collect_constraints("lower")

        return f1 + cu + cl

    def get_nobj(self):
        """Return number of objectives."""
        return self.nobj

    # inequality constraints (equality constraints not required)
    def get_nic(self):
        """Return number of inequality constraints."""
        return self.nic

    def get_bounds(self):
        """Return bounds of decision variables."""
        return self.bounds

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
        return pop
