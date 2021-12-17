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

import numpy as np
import pandas as pd
from collections import OrderedDict

from tespy.components.basics.cycle_closer import CycleCloser
from tespy.components.heat_exchangers.condenser import Condenser
from tespy.components.heat_exchangers.heat_exchanger_simple import HeatExchangerSimple
from tespy.components.nodes.merge import Merge
from tespy.components.nodes.splitter import Splitter
from tespy.components.piping.valve import Valve
from tespy.components.turbomachinery.pump import Pump
from tespy.components.turbomachinery.turbine import Turbine
from tespy.connections.bus import Bus
from tespy.connections.connection import Connection
from tespy.networks.network import Network


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
    - Provide the optimization results DataFrame in the
      :code:`.individuals` attribute of the :code:`OptimizationProblem` class.

    Parameters
    ----------
    model : custom class
        Object of some class, which provides all the methods required by the
        optimization suite, see the
        :py:class:`tespy.tools.optimization.SamplePlant` for a template.

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
    `documentation <https://esa.github.io/pygmo2/>`_.

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
        algo : pygmo.algorithm
            PyGMO optimization algorithm.

        num_ind : int
            Number of individuals.

        num_gen : int
            Number of generations.
        """
        print(type(algo))

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
                print(self.objective_list[i] + ': {}'.format(round(pop.champion_f[i], 4)))
            for i in range(len(self.variable_list)):
                print(self.variable_list[i] + ': {}'.format(round(pop.champion_x[i], 4)))
            pop = algo.evolve(pop)

        gen += 1
        self._process_generation_data(gen, pop)

        print('Final evolution: {}'.format(gen))
        for i in range(len(self.objective_list)):
            print(self.objective_list[i] + ': {}'.format(round(pop.champion_f[i], 4)))
        for i in range(len(self.variable_list)):
            print(self.variable_list[i] + ': {}'.format(round(pop.champion_x[i], 4)))


class SamplePlant:
    """Class template for TESPy model usage in optimization module."""
    def __init__(self):

        self.nw = Network(fluids=['water'])
        self.nw.set_attr(p_unit="bar", T_unit="C", h_unit="kJ / kg", iterinfo=False)

        # main cycle components cycle closer
        steam_generator = HeatExchangerSimple("steam generator")
        close_cycle = CycleCloser("cycle closer")

        turbine_hp = Turbine("turbine high pressure")
        turbine_lp = Turbine("turbine low pressure")
        extraction = Splitter("steam extraction splitter", num_out=2)
        preheater = Condenser("feed water preheater")
        valve = Valve("preheater condensate valve")
        waste_steam_merge = Merge("waste steam merge")

        condenser = HeatExchangerSimple("main condenser")
        feed_pump = Pump("feed water pump")

        # Connections

        # main cycle
        c0 = Connection(steam_generator, "out1", close_cycle, "in1", label="0")
        c1 = Connection(close_cycle, "out1", turbine_hp, "in1", label="1")
        c2 = Connection(turbine_hp, "out1", extraction, "in1", label="2")
        c3 = Connection(extraction, "out1", turbine_lp, "in1", label="3")
        c4 = Connection(turbine_lp, "out1", waste_steam_merge, "in1", label="4")
        c5 = Connection(waste_steam_merge, "out1", condenser, "in1", label="5")
        c6 = Connection(condenser, "out1", feed_pump, "in1", label="6")
        c7 = Connection(feed_pump, "out1", preheater, "in2", label="7")
        c8 = Connection(preheater, "out2", steam_generator, "in1", label="8")

        # steam extraction
        c11 = Connection(extraction, "out2", preheater, "in1", label="11")
        c12 = Connection(preheater, "out1", valve, "in1", label="12")
        c13 = Connection(valve, "out1", waste_steam_merge, "in2", label="13")

        self.nw.add_conns(c0, c1, c2, c3, c4, c5, c6, c7, c8, c11, c12, c13)

        # component specifications
        steam_generator.set_attr(pr=0.92)
        turbine_hp.set_attr(eta_s=0.9)
        turbine_lp.set_attr(eta_s=0.9)
        condenser.set_attr(pr=1)
        feed_pump.set_attr(eta_s=0.75)
        preheater.set_attr(ttd_u=5, pr1=1, pr2=0.98)

        # connection specifications
        c1.set_attr(fluid={'water': 1}, p=100, T=600, m=10)
        # pressure at connection 2 will be the parameter to optimize
        c2.set_attr(p=10)
        c6.set_attr(x=0, T=30)

        power_bus = Bus('power output')
        power_bus.add_comps(
            {'comp': turbine_hp, 'char': 0.97},
            {'comp': turbine_lp, 'char': 0.97},
            {'comp': feed_pump, 'char': 0.97, 'base': 'bus'}
        )
        heat_bus = Bus('heat input')
        heat_bus.add_comps({'comp': steam_generator})
        self.nw.add_busses(power_bus, heat_bus)

        self.nw.solve("design")
        self.stable = "_stable"
        self.nw.save(self.stable)

    def get_param(self, obj, label, parameter):
        """Get the value of a parameter in the network's unit system.

        Parameters
        ----------
        obj : str
            Object to get parameter for (Components/Connections).

        label : str
            Label of the object in the TESPy model.

        parameter : str
            Name of the parameter of the object.

        Returns
        -------
        value : float
            Value of the parameter.
        """
        if obj == 'Components':
            return self.nw.get_comp(label).get_attr(parameter).val
        elif obj == 'Connections':
            return self.nw.get_conn(label).get_attr(parameter).val

    def set_params(self, **kwargs):

        if "Connections" in kwargs:
            for c, params in kwargs["Connections"].items():
                self.nw.get_conn(c).set_attr(**params)

        if "Components" in kwargs:
            for c, params in kwargs["Components"].items():
                self.nw.get_comp(c).set_attr(**params)

    def solve_model(self, **kwargs):

        self.set_params(**kwargs)

        self.solved = False
        try:
            self.nw.solve("design")
            if self.nw.res[-1] >= 1e-3 or self.nw.lin_dep:
                self.nw.solve("design", init_only=True, init_path=self.stable)
            else:
                # might need more checks here!
                if any(self.nw.result['Condenser']['Q'] > 0):
                    self.solved = False
                else:
                    self.solved = True
        except:
            self.nw.lin_dep = True
            self.nw.solve("design", init_only=True, init_path=self.stable)

    def get_objective(self):
        if self.solved:
            return -(
                self.nw.busses['power output'].P.val /
                self.nw.busses['heat bus'].P.val
            )
        else:
            return np.nan
