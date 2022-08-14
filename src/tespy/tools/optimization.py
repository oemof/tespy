try:
    import pygmo as pg
except ImportError:
    pg = None

import pandas as pd

from tespy.tools.helpers import merge_dicts
from tespy.tools.helpers import nested_OrderedDict


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

    objective : str
        Name of the objective. :code:`objective` is passed to the
        :code:`get_objective` method of your tespy model instance.

    Note
    ----
    For the required structure of the input dictionaries see the example in
    below.

    Installation of pygmo via pip is not available for Windows and OSX users
    currently. Please use conda instead or refer to their
    `documentation <https://esa.github.io/pygmo2/>`_.

    Example
    -------
    This example shows the optimization of the thermal efficiency of the
    `SamplePlant` with respect to the pressure value at the intermediate
    extration of the turbine. You can find the example code on the GitHub page
    of TESPy:
    :download:`example implementation </../tutorial/optimization_example.py>`.

    To use the API, you need to define a class holding a TESPy network. You
    create an instance of your plant class, i.e. :code:`plant = SamplePlant()`.
    Then create an instance of the class
    :py:class:`tespy.tools.optimization.OptimizationProblem` and pass

    - the plant instance,
    - the variables,
    - the constraints and
    - the objective function name.

    For the optimization problem in this example, it can be formulated as
    unconstrained problem by defining the lower and the upper limits for the
    variable values, the constraints parameter can be left out. The objective
    function of your plant (:code:`get_objective`), should return the
    evaluation of the objective function. You can define multiple objective
    functions, which can be accessed by the name of the objective. In the
    example code only the thermal efficiency is defined, therefore the
    :code:`objective` keyword does not need to be defined. The keywod is mainly
    of use, if you want to quickly change the evaluation.

    The only variable in this example is the extraction pressure at the
    turbine. The upper limit is 50 bar and the lower limit 0.4 bar. Of course,
    it is possible to use multiple variables and component parameters as
    variables as well. Just provide them in the same structure as in this
    example.

    .. note::

        Please note, that the sense of optimization is always minimization,
        therefore you need to define your objective functions in the
        appropriate way.

    After selection of an appropriate algorithm (differential evolution is a
    good fit for this application) we can start the optimization run. For more
    information on algorithms available in the PyGMO framework and their
    individual specifications please refer to the respective section in their
    online documentation:
    `list of algorithms <https://esa.github.io/pagmo2/overview.html#list-of-algorithms>`__.
    Specify the number of individuals, the number of generations and call the
    :py:meth:`tespy.tools.optimization.OptimizationProblem.run` method of your
    :code:`OptimizationProblem` instance passing the algorithm and the number
    of individials and generations.

    In our sample run, we found an optimal value for the extraction pressure of
    about 4.45 bar for a thermal efficiency of 38.7 %. The results for every
    individual in each generation are stored in the :code:`individuals`
    attribute of the :code:`OptimizationProblem`.
    """

    def __init__(self, model, variables={}, constraints={}, objective="objective"):
        if pg is None:
            msg = (
                "For this function of TESPy pygmo has to be installed. Either use "
                "pip (Linux users only) or conda to install the latest pygmo "
                "version."
            )
            raise ImportError(msg)

        self.model = model
        default_variables = {"Connections": {}, "Components": {}}
        default_constraints = {
            "lower limits": {"Connections": {}, "Components": {}},
            "upper limits": { "Connections": {}, "Components": {}}
        }
        # merge the passed values into the default dictionary structure
        variables = merge_dicts(variables, default_variables)
        constraints = merge_dicts(constraints, default_constraints)

        # pygmo creates a vector for the variables and constraints, which has
        # to be in consistent order. Therefore use OrderedDicts instead of
        # dictionaries
        self.variables = nested_OrderedDict(variables)
        self.constraints = nested_OrderedDict(constraints)
        self.objective = objective
        self.variable_list = []
        self.constraint_list = []

        self.objective_list = [objective]
        self.nobj = len(self.objective_list)

        self.bounds = [[], []]
        for obj, data in self.variables.items():
            for label, params in data.items():
                if obj in ["Connections", "Components"]:
                    for param in params:
                        self.bounds[0] += [self.variables[obj][label][param]['min']]
                        self.bounds[1] += [self.variables[obj][label][param]['max']]
                        self.variable_list += [obj + '-' + label + '-' + param]
                else:
                    self.bounds[0] += [self.variables[obj][label]['min']]
                    self.bounds[1] += [self.variables[obj][label]['max']]
                    self.variable_list += [obj + '-' + label]

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

    def _process_generation_data(self, gen, pop):
        """Process the data of the individuals within one generation.

        Parameters
        ----------
        gen : int
            Generation number.

        pop : pygmo.population
            PyGMO population object.
        """
        individual = 0
        for x in pop.get_x():
            self.individuals.loc[(gen, individual), self.variable_list] = x
            individual += 1

        individual = 0
        for objective in pop.get_f():
            self.individuals.loc[
                (gen, individual),
                self.objective_list + self.constraint_list
            ] = objective
            individual += 1

        self.individuals['valid'] = (
            self.individuals[self.constraint_list] < 0
        ).all(axis='columns')

    def run(self, algo, num_ind, num_gen):
        """Run the optimization algorithm.

        Parameters
        ----------
        algo : pygmo.core
            PyGMO optimization algorithm.

        num_ind : int
            Number of individuals.

        num_gen : int
            Number of generations.
        """

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

        # replace prints with logging
        gen = 0
        for gen in range(num_gen - 1):
            self._process_generation_data(gen, pop)

            print('Evolution: {}'.format(gen))
            for i in range(len(self.objective_list)):
                print(
                    self.objective_list[i] + ': {}'.format(
                        round(pop.champion_f[i], 4)
                    )
                )
            for i in range(len(self.variable_list)):
                print(
                    self.variable_list[i] + ': {}'.format(
                        round(pop.champion_x[i], 4)
                    )
                )
            pop = algo.evolve(pop)

        gen += 1
        self._process_generation_data(gen, pop)

        print('Final evolution: {}'.format(gen))
        for i in range(len(self.objective_list)):
            print(
                self.objective_list[i] + ': {}'.format(
                    round(pop.champion_f[i], 4)
                )
            )
        for i in range(len(self.variable_list)):
            print(
                self.variable_list[i] + ': {}'.format(
                    round(pop.champion_x[i], 4)
                )
            )
