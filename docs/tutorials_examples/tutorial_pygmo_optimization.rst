Cycle optimization example using PyGMO
---------------------------------------

.. contents::
    :depth: 1
    :local:
    :backlinks: top
    

Task
^^^^

Designing a power plant meets multiple different tasks, such as finding the 
optimal fresh steam temperature and pressure to reduce exhaust steam water 
content, or the optimization of extraction pressures to maximize cycle 
efficiency and many more. 
In case of a rather simple power plant topologies the task of finding optimized 
values for e.g. extraction pressures is still managable without any optimization 
tool. As the topology becomes more complexe and boundary conditions come into play 
the usage of additional tools is recommended. 
The following example is intended to show the usage of PyGMO in combination 
with TESPy to maximize the cycle efficiency of a power plant with two extractions.

The source code can be found at the `tespy_examples repository
<https://github.com/oemof/oemof-examples/tree/master/oemof_examples/tespy/clausius_rankine>`_.  


What is PyGMO?
^^^^^^^^^^^^^^

PyGMO (Python Parallel Global Multiobjective Optimizer) is a library that provides 
a large number of evolutionary optimization algroithms. PyGMO can be used to 
solve constrained, unconstrained, single objective and multi objective problems.


Evolutionary Algorithms
+++++++++++++++++++++++

Evolutionary Algorithms are optimization algorithms inspired by biological evolution. 
In a given Population the Algorithm uses the so called fitness function to determine 
the quality of the solutions to each individual (set of decision variables) problem. 
The best possible solution of the population is called champion. Via Mutation, 
Recombination and Selection your Population evolves to find better Solutions. 

Evolutionary Algorithms will never find an exact solution to your problem. 
They are only approximations to the real optimum.


Install PyGMO
+++++++++++++

We recommend installation with conda, which is included with the Anaconda distribution of Python. 
You can install pygmo by opening the anaconda prompt and executing the following::

    conda install -c conda-forge pygmo
    
The installation of all available optimization solvers is included in this process.


Creating your TESPy-Model
^^^^^^^^^^^^^^^^^^^^^^^^^

It is necessary to use object oriented programing in PyGMO. Therefore we creat 
a creat a class PowerPlant which contains our TESPy-Model and a function to return 
the cycle efficiency:

.. code-block:: python

    from tespy.networks import network
    from tespy.components import (
    turbine, splitter, merge, condenser, pump, sink, source,
    heat_exchanger_simple, desuperheater, cycle_closer
    )
    from tespy.connections import connection, bus
    
    class PowerPlant():
        def __init__(self):
            Your TESPy_Model
        
        def calculate_efficiency(self,x):
            # set extraction pressure
            self.nw.connections['extraction1'].set_attr(p=x[0])
            self.nw.connections['extraction2'].set_attr(p=x[1])
            
            self.nw.solve('design')
            self.nw.save('extraction')
                    
            return self.nw.busses['power'].P.val/self.nw.busses['heat'].P.val
        
In calculate_efficiency(self, x) the variable x is a list containing your 
decision variables. This function returns the cycle efficiency for a specific 
set of decision variables.


Creating your PyGMO-Model
^^^^^^^^^^^^^^^^^^^^^^^^^

The optimization in PyGMO starts be defining the problem at hand. You can set 
the number of objectives your problem has in get_nobj(). The number of constraints 
is set in get_nec() (equality constraints) and get_nic() (inequality constraints). 
In get_bounds() you set the bounds of your decision variables. Finally, you define 
your fitness function and constraints in fitness(self, x):

.. code-block:: python

    import pygmo as pg
    
    class optimization_problem():
        def fitness(self, x):
            f1 = 1/self.model.calculate_efficiency(x)
            ci1 = -x[0]+x[1]
            return [f1, ci1]
    
        def get_nobj(self):
            """Return number of objectives."""
            return 1
    
        # equality constraints
        def get_nec(self):
            return 0
    
        # inequality constraints
        def get_nic(self):
            return 1
    
        # # integer dimension
        # def get_nix(self):
        #     return 0
    
        def get_bounds(self):
            """Return bounds of decision variables."""
            return ([1,1], [40,40])
    
        def get_name(self):
            """Return function name."""
            return ""
            
By default PyGMO minimizes the fitness function. Therefore we set the fitness 
function f1 to the reciprocal of the cycle efficiency. We set one inequality 
constraint so that the pressure of the first extraction has to be bigger than 
the second one:

.. math::

    p_{e,1} > p_{e,2}

In PyGMO your inequality constraint has to be in form <0:

.. math::
    - p_{e,1} + p_{e,2} < 0


We expect that the extraction pressure won't be more than 40 bar and not less 
1 bar. Therefore we set the bounds of our decision variables:

.. math::

    1 bar < p_{e,1} < 40 bar
    1 bar < p_{e,2} < 40 bar


Run PyGMO-Optimization
^^^^^^^^^^^^^^^^^^^^^^

The following code shows how to run the PyGMO optimization:

.. code-block:: python

    optimize = optimization_problem()
    optimize.model = PowerPlant()
    prob = pg.problem(optimize)
    
    pop = pg.population(prob, size=20)
    algo = pg.algorithm(pg.nlopt())
    
    for i in range(15):
        print(1/pop.champion_f[0]*100, pop.champion_x)
        p = [pop.champion_x[0], pop.champion_x[1]]
        pop = algo.evolve(pop)


    print()
    print('Efficiency: {} %'.format(round(100/pop.champion_f[0],4)))
    print('Extraction 1: {} bar'.format(round(p[0],4)))
    print('Extraction 2: {} bar'.format(round(p[1],4)))


With optimize you tell PyGMO which problem you want to optimize. In the class 
optimization_problem() we defined our problem be setting fitness function 
and inequality constraint. With optimize.model we set the model we want to optimize. 
In our case we want to optimize the extraction pressures in our PowerPlant(). 
Finally, our problem is set in prob = pg.problem(optimize).

With pop we define the the size of each population in our problem, algo is used 
to set the algorithm. A list of available algorithms can be found in
<https://esa.github.io/pygmo2/overview.html#list-of-algorithms>`_. The choice 
of your algorithm depends on the type of problem. Have you set equality or 
inequality constraints? Do you perform a single- or multi-objective optimization?
 
In a for-loop we evolve and print the champion of our last population.
After 15 generations we print our final champion.
