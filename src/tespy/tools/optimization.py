try:
    import pygmo as pg
except ImportError:
    msg = (
        "For this function of TESPy pygmo has to be installed. Either use "
        "pip (Linux users only) or conda to install the latest pygmo "
        "version. It is also possible to install tespy using the [opt] "
        "option: pip install tespy[opt]."
    )
    raise ImportError(msg)

import pandas as pd
from collections import OrderedDict


def _nested_OrderedDict(dictionary):
    """Create a nested OrderedDict from a nested dict.

    Parameters
    ----------
    dictionary : dict
        Nested dict.

    Returns
    -------
    dictionary : collections.OrderedDict
        Nested OrderedDict.
    """
    dictionary = OrderedDict(dictionary)
    for key, value in dictionary.items():
        if isinstance(value, dict):
            dictionary[key] = _nested_OrderedDict(value)

    return dictionary


class OptimizationProblem:
    r"""
    The OptimizationProblem handles the optimization.

    - Set up the optimization problems by specifying constraints, upper and
      lower bounds for the decision variables and selection of the objective
      function.
    - Run the optimization, see
      :py:meth:`tespy.tools.optimization.OptimizationProblem.run`.
    - Provide the optimization results DataFrame
      :code:`OptimizationProblem.individuals`.

    Parameters
    ----------
    model : custom class
        Object of some class, which provides all the methods required by the
        optimization suite, see ...

    variables : dict
        Dictionary containing the decision variables and their respective
        bounds.

    constraints : dict
        Dictionary containing the constraints for the model.

    objective : str
        Name of the objective. :code:`objective` is passed to the
        :code:`get_objective` method of your tespy model instance.

    Note
    ----
    Installation of pygmo via pip is not available for Windows and OSX users
    currently. Please use conda instead or refer to their
    `documentation <https://esa.github.io/pygmo2/>`_

    Example
    -------
    Add some example code here, maybe refer to some repositories.
    """

    def __init__(self, model, variables, constraints, objective):
        self.model = model
        # use OrderedDicts to have consistent order of variables,
        # constraints (and objectives in the future)
        self.variables = _nested_OrderedDict(variables)
        self.constraints = _nested_OrderedDict(constraints)
        self.objective = objective
        self.variable_list = []
        self.constraint_list = []

        self.objective_list = [objective]
        self.nobj = len(self.objective_list)

        self.bounds = [[], []]
        for obj, data in self.variables.items():
            for label, params in data.items():
                for param in params:
                    self.bounds[0] += [self.variables[obj][label][param]['min']]
                    self.bounds[1] += [self.variables[obj][label][param]['max']]
                    self.variable_list += [obj + '-' + label + '-' + param]

        self.input_dict = self.variables.copy()

        self.nic = 0
        for obj, data in self.constraints['upper limits'].items():
            for label, constraints in data.items():
                for param, constraint in constraints.items():
                    self.nic += 1
                    self.constraint_list += [
                        obj + '-' + label + '-' + param + ' <= ' +
                        str(constraint)
                    ]

        for obj, data in self.constraints['lower limits'].items():
            for label, constraints in data.items():
                for param, constraint in constraints.items():
                    self.nic += 1
                    self.constraint_list += [
                        obj + '-' + label + '-' + param + ' >= ' +
                        str(constraint)
                    ]

    def fitness(self, x):
        """Fitness function."""
        i = 0
        for obj, data in self.variables.items():
            for label, params in data.items():
                for param in params:
                    self.input_dict[obj][label][param] = x[i]
                    i += 1

        self.model.solve_model(**self.input_dict)
        f1 = [self.model.get_objective(self.objective)]

        cu = [
            self.model.get_param(obj, label, parameter) - constraint
            for obj, data in self.constraints['upper limits'].items()
            for label, constraints in data.items()
            for parameter, constraint in constraints.items()
        ]
        cl = [
            constraint - self.model.get_param(obj, label, parameter)
            for obj, data in self.constraints['lower limits'].items()
            for label, constraints in data.items()
            for parameter, constraint in constraints.items()
        ]

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

    # throw individuals, params_list, objectives_list, constraint_list into optimization class!
    def _process_generation_data(self, gen, pop):

        individual = 0
        for x in pop.get_x():
            self.individuals.loc[(gen, individual), self.variable_list] = x
            individual += 1

        individual = 0
        for objective in pop.get_f():
            self.individuals.loc[(gen, individual), self.objective_list + self.constraint_list] = objective
            individual += 1

        self.individuals['valid'] = (
            self.individuals[self.constraint_list] < 0
        ).all(axis='columns')

    def run(self, algo, num_ind, num_gen):

        self.individuals = pd.DataFrame(
            index=range(num_gen * num_ind)
        )

        self.individuals["gen"] = [
            gen for gen in range(num_gen) for ind in range(num_ind)
        ]
        self.individuals["ind"] = [
            ind for gen in range(num_gen) for ind in range(num_ind)
        ]

        self.individuals.set_index(["gen", "ind"], inplace=True)

        algo = pg.algorithm(algo)
        pop = pg.population(pg.problem(self), size=num_ind)

        gen = 0
        for gen in range(num_gen - 1):
            self._process_generation_data(gen, pop)
            print()
            print('Evolution: {}'.format(gen))
            for i in range(len(self.objective_list)):
                print(self.objective_list[i] + ': {}'.format(round(pop.champion_f[i], 4)))
            for i in range(len(self.variable_list)):
                print(self.variable_list[i] + ': {}'.format(round(pop.champion_x[i], 4)))
            pop = algo.evolve(pop)

        gen += 1
        self._process_generation_data(gen, pop)

        print()

        print('Final evolution: {}'.format(gen))
        for i in range(len(self.objective_list)):
            print(self.objective_list[i] + ': {}'.format(round(pop.champion_f[i], 4)))
        for i in range(len(self.variable_list)):
            print(self.variable_list[i] + ': {}'.format(round(pop.champion_x[i], 4)))

        print()
