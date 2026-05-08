.. _tutorial_optimization_label:

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
following tutorial is intended to show the usage of pymoo in combination with
TESPy to **maximize the cycle efficiency of a power plant with two**
**extractions.**

You can download the code here:
:download:`optimization_example.py </../tutorial/advanced/optimization_example.py>`

.. figure:: /_static/images/tutorials/optimization/flowsheet.svg
    :align: center
    :alt: Topology of the power plant
    :figclass: only-light

    Figure: Topology of the power plant

.. figure:: /_static/images/tutorials/optimization/flowsheet_darkmode.svg
    :align: center
    :alt: Topology of the power plant
    :figclass: only-dark

    Figure: Topology of the power plant

What is pymoo?
^^^^^^^^^^^^^^

`pymoo <https://pymoo.org/>`__ ( Multi-objective Optimization in Python,
:cite:`pymoo`) is a library that provides numerous optimization algorithms.
pymoo can be used to solve constrained, unconstrained, single objective and
multi objective problems. In this example we will use an evolutionary algorithm
to find the optimal pressure values.

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

Install pymoo
+++++++++++++

Installation of pymoo is straight-forward: Just install tespy with the extra
dependencies "opt", e.g. with pip

.. code-block:: bash

    pip install tespy[opt]

of uv

.. code-block:: bash

    uv add tespy --extra opt

Creating your TESPy-Model
^^^^^^^^^^^^^^^^^^^^^^^^^

We use the :py:class:`ModelTemplate <tespy.models.template.ModelTemplate>` base
class, which provides parameter access, solving and optimization infrastructure
automatically. Inheriting from it requires implementing three methods:

- :code:`_create_network` - assembles the TESPy network and produces a stable
  initial solution. Always call :code:`super()._create_network()` first.
- :code:`_parameter_lookup` - returns a :code:`dict` that maps human-readable
  parameter names to their location in the network. These names are used
  directly by the optimizer, :code:`get_parameter` and :code:`set_parameters`.
- :code:`solve_model` - routes calls to the appropriate solve logic (e.g.
  :code:`solve_model_design` or custom logic).

The :code:`_create_network` method below builds the power plant with two
extraction turbines, sets boundary conditions and runs the initial
design-point solve.

.. dropdown:: Display source code for the SamplePlant class

    .. literalinclude:: /../tutorial/advanced/optimization_example.py
        :language: python
        :start-after: [sec_1]
        :end-before: [sec_2]

Next, :code:`_parameter_lookup` registers lookup between model internals and
keyword parameter specifications for inputs and outputs. Connection attributes
using the :code:`["Connections", label, attr]` form, Component attributes
follow the :code:`["Components", label, attr]` form. The thermal efficiency is
registered as a read-only derived quantity via :code:`{"get": callable}`:

.. math::

    \eta_\text{th}=\frac{|\sum P|}{\dot{Q}_{sg}}

The feasibility check inside :code:`calc_efficiency` returns :code:`np.nan`
for non-converged or physically infeasible solutions, which the optimizer
treats as a penalty automatically.

.. dropdown:: Display source code for :code:`_parameter_lookup`, :code:`calc_efficiency` and :code:`solve_model`

    .. literalinclude:: /../tutorial/advanced/optimization_example.py
        :language: python
        :start-after: [sec_2]
        :end-before: [sec_3]

After instantiating the model, all registered parameters are immediately
accessible via :code:`get_parameter`:

.. literalinclude:: /../tutorial/advanced/optimization_example.py
    :language: python
    :start-after: [sec_3]
    :end-before: [sec_4]

.. attention::

    The sense of optimization is minimization by default. We can change the
    sense for each of the objectives we pass to
    :py:meth:`optimize <tespy.models.template.ModelTemplate.optimize>`
    by passing a list with :code:`True` or :code:`False` for each individual
    objective **in identical order as the objectives** using the
    :code:`minimize_flags` argument.

We set one inequality constraint - the first extraction pressure must exceed
the second:

.. math::

    p_{e,1} > p_{e,2}

This is expressed inside the :code:`constraints` dictionary by referencing the
name of the first parameter as the lower bound of the second. The :code:`kpi`
argument selects additional model outputs to include in the result log. Here
we track the high-pressure turbine power and pressure ratio.

Run pymoo-Optimization
^^^^^^^^^^^^^^^^^^^^^^

:py:meth:`model.optimize <tespy.models.template.ModelTemplate.optimize>` wraps
:py:class:`OptimizationProblem <tespy.tools.optimization.OptimizationProblem>`
and a pymoo algorithm into a single call and returns a :code:`DataFrame` of all
evaluated individuals. For more information on available algorithms see the
`list of algorithms <https://pymoo.org/algorithms/index.html>`__.

.. literalinclude:: /../tutorial/advanced/optimization_example.py
    :language: python
    :start-after: [sec_5]
    :end-before: [sec_6]

In our run, we got the following optimal solution:

.. code:: bash

    Efficiency: 44.82 %
    Extraction 1: 25.753 bar
    Extraction 2: 2.685 bar

.. figure:: /_static/images/tutorials/optimization/optimization_result.svg
    :align: center
    :alt: Scatter plot for all individuals during the optimization
    :figclass: only-light

    Figure: Scatter plot for all individuals during the optimization

.. figure:: /_static/images/tutorials/optimization/optimization_result_darkmode.svg
    :align: center
    :alt: Scatter plot for all individuals during the optimization
    :figclass: only-dark

    Figure: Scatter plot for all individuals during the optimization

You can access the data generated by the optimization from the returned
:code:`DataFrame`, which contains every evaluated individual. The plotting code
in the download shows how to filter for feasible solutions and highlight the
optimum.

.. literalinclude:: /../tutorial/advanced/optimization_example.py
    :language: python
    :start-after: [sec_6]
    :end-before: [sec_7]
