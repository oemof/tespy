.. _tespy_tutorial_pygmo_optimization_label:

Thermal Power Plant Efficiency Optimization
-------------------------------------------

Task
^^^^
Designing a power plant meets multiple different tasks, such as finding the
optimal fresh steam temperature and pressure to reduce exhaust steam water
content, or the optimization of extraction pressures to maximize cycle
efficiency and many more.

In case of a rather simple power plant topologies the task of finding
optimized values for e.g. extraction pressures is still manageable without any
optimization tool. As the topology becomes more complex and boundary
conditions come into play the usage of additional tools is recommended. The
following tutorial is intended to show the usage of PyGMO in combination with
TESPy to **maximize the cycle efficiency of a power plant with two**
**extractions.**

You can download the code here:
:download:`optimization_example.py </../tutorial/advanced/optimization_example.py>`

.. figure:: /_static/images/tutorials/pygmo_optimization/flowsheet.svg
    :align: center
    :alt: Topology of the power plant
    :figclass: only-light

    Figure: Topology of the power plant

.. figure:: /_static/images/tutorials/pygmo_optimization/flowsheet_darkmode.svg
    :align: center
    :alt: Topology of the power plant
    :figclass: only-dark

    Figure: Topology of the power plant

What is PyGMO?
^^^^^^^^^^^^^^

PyGMO (Python Parallel Global Multiobjective Optimizer, :cite:`Biscani2020`)
is a library that provides numerous evolutionary optimization algorithms. PyGMO
can be used to solve constrained, unconstrained, single objective and multi
objective problems.

Evolutionary Algorithms
+++++++++++++++++++++++

Evolutionary Algorithms (EA) are optimization algorithms inspired by biological
evolution. In a given population the algorithm uses the so-called fitness
function to determine the quality of the solutions to each individual (set of
decision variables) problem. The best possible solution of the population is
called champion. Via mutation, recombination and selection your population
evolves to find better solutions.

EA will never find an exact solution to your problem. They can only give an
approximation for the real optimum.

Install PyGMO
+++++++++++++

.. tab-set::

   .. tab-item:: conda

        With the conda package manager PyGMO is available for Linux, OSX and
        Windows thanks to the infrastructure of
        `conda-forge <https://conda-forge.org/>`_:

        .. code-block:: bash

            conda install -c conda-forge pygmo


        Windows user can perform an installation from source as an alternative
        to conda. For further information on this process we recommend the
        `PyGMO installation <https://esa.github.io/pygmo2/install.html#installation-from-source>`__
        accordingly.


   .. tab-item:: pip (Linux only!)

        On Linux you also have the option to use the
        `pip <https://pip.pypa.io/en/stable/>`_ package installer:

        .. code-block:: bash

            pip install pygmo

Creating your TESPy-Model
^^^^^^^^^^^^^^^^^^^^^^^^^
To use the API, you need to define a class holding a TESPy network. The
initialization of the class should build the network and run an initial
simulation. Furthermore, you have to define methods

- to get component or connection parameters of the plant :code:`get_param`,
- to run a new simulation for every new input from PyGMO :code:`solve_model`
  and
- to return the objective values :code:`get_objectives`.

First, we set up the class with the TESPy network.

.. dropdown:: Display source code for the PowerPlant class

    .. literalinclude:: /../tutorial/advanced/optimization_example.py
        :language: python
        :start-after: [sec_1]
        :end-before: [sec_2]


Next, we add the methods :code:`get_param`, :code:`solve_model` and
:code:`get_objectives`. On top of that, we add a setter working similarly as the
getter. The objective is to maximize thermal efficiency as defined in the
equation below. The :code:`get_objectives` method calls a :code:`get_objective`
method and collects all objectives values. This is useful if you are
implementing pareto or multi-objective problems.

.. math::

    \eta_\mathrm{th}=\frac{|\sum P|}{\dot{Q}_{sg}}

.. attention::

    The sense of optimization is minimization by default. We can change the
    sense for each of the objectives we pass to the :code:`OptimizationProblem`
    class by passing a list with :code:`True` or :code:`False` for each
    inidivual objective **in identical order as the objectives** using the
    :code:`minimize` argument. In this example we only use a single objective,
    so there is not too much, that can go wrong.

We also have to make sure, only the results of physically feasible solutions
are returned. In case we have infeasible solutions, we can simply return
:code:`np.nan`. An infeasible solution is obtained in case the power of a
turbine is positive, the power of a pump is negative, or the heat exchanged
in any of the preheaters is positive. We also check, if the calculation does
converge.

.. dropdown:: Display source code for the class methods

    .. literalinclude:: /../tutorial/advanced/optimization_example.py
        :language: python
        :start-after: [sec_2]
        :end-before: [sec_3]

After this, we import the
:py:class:`tespy.tools.optimization.OptimizationProblem` class and create an
instance of our self defined class, which we pass to an instance of the
OptimizationProblem class. We also have to pass

- the variables to optimize,
- the constraints to consider,
- the objective function name (you could define multiple in the
  :code:`get_objective` method if you wanted) and
- the optimization sense for each objective, here we pass :code:`[False]` since
  we want to maximize efficiency.

We set one inequality constraint, namely that the pressure of the first
extraction has to be higher than the pressure at the second one:

.. math::

    p_{e,1} > p_{e,2}

To do this, we can set a lower limit for the pressure at connection 2 and
reference the pressure at connection 4 as seen in the code:

.. literalinclude:: /../tutorial/advanced/optimization_example.py
    :language: python
    :start-after: [sec_3]
    :end-before: [sec_4]

Before we can run the optimization, we only need to select an appropriate
algorithm. After that we can start the optimization run. For more information
on algorithms available in the PyGMO framework and their individual
specifications please refer to the respective section in their online
documentation:
`list of algorithms <https://esa.github.io/pagmo2/overview.html#list-of-algorithms>`__.
Create an initial population and then specify the number of individuals, the
number of generations and call the
:py:meth:`tespy.tools.optimization.OptimizationProblem.run` method of your
:code:`OptimizationProblem` instance passing the algorithm, the population and
the number of individuals and generations.

Run PyGMO-Optimization
^^^^^^^^^^^^^^^^^^^^^^
The following code then simply runs the PyGMO optimization.

.. literalinclude:: /../tutorial/advanced/optimization_example.py
    :language: python
    :start-after: [sec_4]
    :end-before: [sec_5]

In our run, we got:

.. code:: bash

    Efficiency: 44.82 %
    Extraction 1: 26.462 bar
    Extraction 2: 2.820 bar

.. figure:: /_static/images/tutorials/pygmo_optimization/pygmo_optimization.svg
    :align: center
    :alt: Scatter plot for all individuals during the optimization
    :figclass: only-light

    Figure: Scatter plot for all individuals during the optimization

.. figure:: /_static/images/tutorials/pygmo_optimization/pygmo_optimization_darkmode.svg
    :align: center
    :alt: Scatter plot for all individuals during the optimization
    :figclass: only-dark

    Figure: Scatter plot for all individuals during the optimization

Finally, you can access the individuals in each of the generations, and you
can have a look at you population. For more info on the population API please
visit the pygmo documentation.

.. literalinclude:: /../tutorial/advanced/optimization_example.py
    :language: python
    :start-after: [sec_5]
    :end-before: [sec_6]
