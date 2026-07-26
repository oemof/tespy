.. _modules_solver_label:

Solver
======
This section provides in depth information on how the solver works in TESPy.
The solver transforms the plant and the boundary conditions into a mathematical
problem and finds its solution. In that, several steps are taken:

- Reducing the variable space and problem complexity
- Selecting starting values
- Decomposing the system of equations into subproblems
- Solving the mathematical problem

Furthermore, we show how you can interact with the process to debug your models.
If you are looking for the modeling API have a look at the
:ref:`network section <modules_networks_label>` instead.

.. tip::

    An exhaustive walkthrough at the concrete examples of a simple heat pump
    and an ORC power cycle are available.

.. _solver_pipeline_label:

From model to solution
----------------------
A TESPy model is represented by a nonlinear system of equations. The system
variables are a subset of

* mass flow, pressure, enthalpy and fluid mass fractions of a fluid in every
  :code:`Connection` type instance,
* energy flow of every :code:`PowerConnection` and :code:`HeatConnection` and
* component variables (e.g. a variable compressor frequency) if declared.

Whenever you call :code:`Network.solve`, the following steps are taken:

1. **Check and preprocess**: the topology is checked for consistency,
   components are initialised and for offdesign calculations the
   design/offdesign switch is performed.
2. **Construct and reduce the problem**: all equations are collected, variables
   coupled by linear relations (e.g. fixed `dp` of a heat exchanger) are merged
   and variables directly determined by the specifications is presolved.
3. **Generate starting values** for the remaining variables.
4. **Decompose** the reduced system into blocks that can be solved in parallel
   and/or one after another.
5. **Solve** the blocks in order or the full system simultaneously with the
   Newton-Raphson method.
6. **Postprocess**: all derived properties and component results are
   calculated and checked against their limits.

.. figure placeholder: solve process flowchart (light/dark), planned

The :py:class:`~tespy.networks.network.Network` instance orchestrates the
modeling steps (topology, units, design/offdesign switch, starting values)
and hands the mathematical work to a :py:class:`~tespy.solver.problem.Problem`
instance. That instance is created freshly for every solve and remains
accessible afterwards through the read-only :code:`Network.problem` property.
The connections and components declare their linear equations, provide
residuals and derivatives of their nonlinear equations together with the
list of variables each equation depends on, and contribute starting value
information. The :code:`Problem` assembles all of this, runs the variable space
reduction, builds the incidence and the block decomposition and executes the
solution strategies.

.. figure placeholder: class interaction diagram Network/Problem/
   Connection/Component/StructureGraph, planned

.. _solver_reduction_label:

Problem construction and reduction
----------------------------------
To reduce the size of the system of equations, the variable space is simplified
before solving. For every primary variable a value directly specified by the
user removes that variable from the variable space. On top of that, three
structural steps are performed:

- merging variables coupled by linear relations,
- simplifying the fluid vectors and
- presolving pressure and enthalpy from state specifications if the fluid
  composition is already fixed.

Many component equations couple two variables in an exactly linear way, for
example

- equality of mass flow or fluid composition at inlet and outlet,
- equality of pressure at the inlet and each outlet of a splitter,
- a constant ratio of inlet and outlet pressure through a specified pressure
  ratio,
- the linear relation imposed by a :code:`Ref` specification.

These relations form a graph. Every connected group of variables collapses
into a single variable the solver iterates. In a simple Clausius Rankine cycle,
for example, a single mass flow variable remains for the whole cycle.

The fluid vectors are reduced next: user specified mass fractions are fixed,
and if only a single fluid of a vector remains unknown, its mass fraction is
the difference of the specified ones to 1. Only undetermined fractions stay in
the variable space.

Finally, state specifications are presolved wherever the fluid composition is
fixed: a temperature at known pressure determines the enthalpy, a vapor
quality at known temperature determines pressure and enthalpy, and every
newly determined value can trigger further determinations. This runs in
rounds until nothing new can be derived.

The reduction can be inspected before solving:

.. code-block:: python

    nw.solve("design", init_only=True)
    nw.print_variables_before_presolve()
    nw.print_presolved_variables()
    nw.print_presolved_equations()
    nw.print_variables()
    nw.print_equations()
    nw.print_equations_with_dependents()
    nw.print_incidence_matrix()

- :py:meth:`~tespy.networks.network.Network.print_presolved_variables` prints
  a table of all variables already solved during preprocessing, identified by
  their object label and property, e.g. :code:`c1 (p)`.
- :py:meth:`~tespy.networks.network.Network.print_presolved_equations` prints
  a table of the equations used to resolve those variables, identified by
  object label and equation name, e.g. :code:`c1.T`.
- :py:meth:`~tespy.networks.network.Network.print_variables` prints a table
  of the remaining variables passed to the solver. Each has a short label for
  unique identification (e.g. :code:`h1`) and the original variables it
  represents (indicated by the original label of the :code:`Connection` and
  the respective variable, e.g. :code:`c1.h`).
- :py:meth:`~tespy.networks.network.Network.print_equations` prints a table
  of the remaining equations passed to the solver, identified by object label
  and equation name, e.g. :code:`compressor.eta_s`.
- :py:meth:`~tespy.networks.network.Network.print_equations_with_dependents`
  extends the equations table with the variables each equation depends on.
- :py:meth:`~tespy.networks.network.Network.print_incidence_matrix` prints
  the same table in a matrix form, i.e. the incidence matrix. The equations are
  the rows and variables are the columns, using :code:`x` for a dependency and
  :code:`-` for none.

The underlying data can also be retrieved programmatically via the
corresponding :code:`get_*` methods.

.. _solver_starting_values_label:

Starting values
---------------
The numerical solving requires initial guesses for all variables. **High**
**quality initial values are crucial for convergence speed and stability**,
bad starting values might lead to instability and diverging calculation can
be the result. The initial guesses for the variables are taken over or
generated from the following sources:

1. **User values and previous solutions**: starting values specified through
   the :code:`val0` attributes, the result of a previous solve (unless
   :code:`init_previous=False`) and values loaded from an :code:`init_path`
   are imposed first and are never overwritten be the subsequent steps.
2. **Propagation**: Starting from specifications, presolved values and existing
   initial guesses are taken as starting points from which approximate
   component relations (e.g. a generic pressure ratio across a compressor) fill
   neighboring connection information. This is applied for mass flow and
   pressure only in a first step, enthalpies follow later.
3. **Temperature reconciliation**: some components provide approximate
   temperature relations between their ports, e.g. heat exchangers or
   turbomachinery. Each of those relations carries a weight and then a least
   squares pass anchored at known temperature levels reconciles one temperature
   per connection, carrying known levels across heat exchangers into circuits
   without any own anchor.
4. **Phase aware enthalpy and pressure anchors**: components declare the
   expected phase at their ports (a compressor inlet is gas, a condenser
   outlet is liquid, a drum is on the saturation line). Where the phase is
   known, the reconciled temperature is converted into saturation pressure
   anchors and phase consistent enthalpies.
5. **Flow presolve**: with all pressures and enthalpies guessed, the
   equations determining mass and energy flows, e.g. balances, heat and power
   specifications, are linear in those variables, so one restricted least
   squares step assigns consistent mass flow and energy flow levels.
6. **Generic fallback** values fill whatever remains, and user specified
   states (vapor quality, saturation offsets) are imposed authoritatively at
   the end.

.. figure placeholder: starting value pipeline diagram, planned

To initialise a calculation from a saved state, provide it via the
:code:`init_path` argument with either a file path to a json file or a
:code:`dict` returned by :code:`nw.save(as_dict=True)`. TESPy searches
through the connection entries by label and takes over the starting values
for the system variables where a label matches.

.. note::

    The files do not need to contain all connections of your network. You can
    build your network step by step and initialise the existing parts of your
    network from the :code:`init_path`.

.. _solver_decomposition_label:

Decomposition into blocks
-------------------------
Before solving, TESPy analyses which equation can determine which variable
and in which order the equations can be processed. The input of this
analysis is the incidence matrix, which states which variables each equation
depends on. Three steps turn the incidence into a solution plan (the method is
known as Dulmage-Mendelsohn decomposition, see :cite:`Parker2023` for its
application to the analysis of process models):

1. A **maximum bipartite matching** assigns every equation one variable it is
   responsible for determining. If that assignment is impossible, the
   unmatched equations and variables delimit the over- or under-determined
   part of the model, which
   :py:meth:`~tespy.networks.network.Network.print_structural_analysis`
   reports for debugging.
2. Equations whose assignments depend on each other in a cycle can not be
   ordered. In that way they form a **block** that must be solved
   simultaneously.
3. The blocks are sorted in a way that every block only needs results of blocks
   solved before it.

The result of this is the incidence matrix in block lower triangular form. The
blocks are solved in order and each block is a small subproblem with all its
external inputs already known. The decomposition depends only on the structure
of your model, not on any values. You can inspect the structure in the
following way:

.. code-block:: python

    nw.print_blocks()
    nw.print_incidence_matrix(block_order=True)

- :py:meth:`~tespy.networks.network.Network.print_blocks` lists the blocks in
  solve order with their kind, prerequisites, equations, variables and solve
  status.
- :py:meth:`~tespy.networks.network.Network.print_incidence_matrix` with
  :code:`block_order=True` sorts equations and variables in block solve
  order, marking every equation's entries within its own block with
  :code:`X` and dependencies on variables of preceding blocks with
  :code:`x`, so the block lower triangular form becomes visible.

The :ref:`advanced debugging tutorial <tutorial_advanced_debugging_label>`
walks through the usage of these methods on a concrete model.

.. _networks_solving_label:

Solving
-------
By default the blocks are solved one after another in their precedence order. A
block consisting of a single equation is solved with a scalar newton algorithm,
which falls back to bracketing the root with Brent's method when the residual
oscillates or changes sign between iterations. Coupled blocks are solved with
the Newton-Raphson method restricted to the block's equations and variables.
Once a block is solved, its variables are fixed for all subsequent blocks.

If a failure occurs in the solving process of a block, the following steps are
taken:

- The failed block first restores its variables' values at start of solving and
  then it is solved with all remaining blocks as one coupled system, keeping
  the progress of the already solved blocks.
- If this fails as well, the solver restarts with the simultaneous solution of
  the full system from its initial state, which reproduces the behavior of a
  plain simultaneous Newton solver exactly.

.. tip::

    If a block fails to solve, it is also possible to pause the solver in that
    exact state to try to investigate the reason for the failure. More
    information is available in the section on
    :ref:`interaction with the solver <solver_interaction_label>`.

After a successful block phase the full residual vector is verified and a
simultaneous polish runs in case couplings outside of the declared incidence
left a remaining residual. The block phase can therefore never degrade
robustness compared to solving everything at once.

.. note::

    Currently, networks with variable fluid composition are always solved
    simultaneously. The simultaneous solution can also be enforced with
    :code:`nw.solve(mode, block_solve=False)`.

Newton-Raphson method
^^^^^^^^^^^^^^^^^^^^^
The Newton-Raphson method requires the calculation of residual values for the
equations and of the partial derivatives to all system variables (Jacobian
matrix). In the next step the matrix is inverted and multiplied with the
residual vector to calculate the increment for the system variables. This
process is repeated until every equation's result in the system is "correct".
All equations are of the same structure:

.. math::

    0 = \text{expression}

calculate the residuals

.. math::

    f(\vec{x}_i)

Jacobian matrix J

.. math::

    J(\vec{x})=\left(\begin{array}{cccc}
    \frac{\partial f_1}{\partial x_1} & \frac{\partial f_1}{\partial x_2} &
    \cdots & \frac{\partial f_1}{\partial x_n} \\
    \frac{\partial f_2}{\partial x_1} & \frac{\partial f_2}{\partial x_2} &
    \cdots & \frac{\partial f_2}{\partial x_n} \\
    \vdots & \vdots & \ddots & \vdots \\
    \frac{\partial f_n}{\partial x_1} & \frac{\partial f_n}{\partial x_2} &
    \cdots & \frac{\partial f_n}{\partial x_n}
    \end{array}\right)

derive the increment

.. math::

    \vec{x}_{i+1}=\vec{x}_i-J(\vec{x}_i)^{-1}\cdot f(\vec{x}_i)

while

.. math::

    ||f(\vec{x}_i)|| > \epsilon

For the convergence check each residual is divided by its per equation scale,
which is the response of the equation to relative changes of the variables it
depends on. The scaled residual therefore reads as the relative change of the
variables required to close the equation, which makes the acceptance criterion
independent of the physical magnitude of the plant. E.g. an energy balance of a
megawatt scale heat exchanger is judged by the same relative measure as a
kilowatt scale one. The solution is accepted once the largest scaled residual
and the increments relative to the variable magnitudes fall below the internal
tolerances, and only in an iteration in which the convergence heuristics did
not modify any value.

.. _module_convergence_label:

Convergence stability
^^^^^^^^^^^^^^^^^^^^^
One of the main downsides of the Newton-Raphson method is that the initial
step width is very large and that it does not know physical boundaries, for
example mass fractions smaller than 0 and larger than 1 or negative pressure.
Also, the large step width can adjust enthalpy or pressure to quantities that
are not covered by the fluid property databases. This would cause an inability
e.g. to calculate a temperature from pressure and enthalpy in the next
iteration of the algorithm. To improve convergence stability, the solver
limits the step width and checks the values against boundaries after every
update of the system variables.

**Step limiting**

* Increments on pressure and on fluid mass fractions are relaxed adaptively:
  a single step can at most halve the value. A negative pressure can
  therefore never be the result of a single large step.
* Fluid mass fractions smaller than 0 and larger than 1 are cut off and the
  fractions of a vector are renormalized in the first iterations.
* With :code:`nw.solve(mode, robust_relax=True)` all increments are
  additionally relaxed globally, starting at 5 % of the full step and ramping
  up to the full step over the first quarter of the iterations. This can help
  with very sensitive models at the cost of extra iterations.

**Value boundaries**

* Fluid properties of pure fluids are checked against the available ranges
  of the property database and readjusted if necessary. Specifications
  restrict the ranges further: a connection with a two-phase specification
  (vapor quality, saturation temperature offsets or the state keyword) keeps
  its pressure below the critical pressure, so saturation properties remain
  accessible. The same applies to ports a component places on the saturation
  line, e.g. the outlets of a drum. A specified state keeps the enthalpy on
  the respective side of the two-phase region in the early iterations.
* For mixtures the user specified value ranges (:code:`p_range, h_range`)
  apply in the first iterations of a calculation without good starting
  values. If a temperature is specified, the enthalpy range is further
  restricted to the temperature range covered by all fluids of the mixture.
* Mass flows are kept inside :code:`m_range` on all connections, and
  component variables are kept inside their :code:`min_val` and
  :code:`max_val` bounds.

A variable representing a group of linearly dependent variables collects the
boundaries of all of its represented variables: each boundary is transformed
through the linear relation and the tightest one applies.

**Component heuristics**

In the first iterations of a design calculation the fluid properties of
the connections are checked against the components they are connecting.
For example, the pressure at the outlet of a turbine should be lower than
the pressure at the inlet. If there are violations, the corresponding
variables are manipulated. If you want to look up, what exactly the
convergence check for a specific component does, look out for the
:code:`convergence_check` methods in the
:py:mod:`tespy.components module <tespy.components>`.

These manipulations only engage while the iteration is still taking
meaningful steps, a converged state is never modified and any modification
postpones the acceptance of a solution to the next iteration. During
block-wise solving all checks apply per block, restricted to the block's
connections and components. Additionally, when the residual of an equation
oscillates between two values the increments of the variables that equation
depends on are dampened automatically. The damping can also be activated from
the first iteration on with :code:`nw.solve(mode, oscillation_damping=True)`.

Status codes
^^^^^^^^^^^^
To check if the solver successfully found a solution for your model you can
check the :code:`.status` attribute of the Network class after calling the
:code:`solve` method. It will be

- 0 in case a solution was successfully found
- 1 in case a solution was found, but some parameters violate physical
  limits
- 2 in case no convergence was achieved after completion of the iterations
- 3 in case a linear dependency in the Jacobian matrix is found
- 11 in case the number of specified parameters is too small for the given
  problem
- 12 in case the number of specified parameters is too large for the given
  problem
- 20 in case the solution process is paused at a block that did not
  converge (only with :code:`pause_on_block_failure=True`)
- 99 in case the simulation crashed due to any other reason

The :code:`solve` does not exit with an exception in case the status is
0, 1, 2, 3 or 99. If you want to raise an error in your script, you can call
the :code:`Network.assert_convergence()` method. It will raise an
:code:`AssertionError` if the simulation did not find a converged solution.

.. _solver_interaction_label:

Interacting with the solver
---------------------------
You can also run the solution process in two stages:
:py:meth:`~tespy.networks.network.Network.presolve` prepares and presolves
the problem and stops before solving, so you can inspect the problem and
modify variable values.
:py:meth:`~tespy.networks.network.Network.solve_continue` then runs the
solver and the postprocessing. :code:`solve` is simply the combination of
the two stages. Between the stages,
:py:meth:`~tespy.networks.network.Network.get_variable_values` and
:py:meth:`~tespy.networks.network.Network.print_variable_values` list every
variable with all linearly dependent variables it represents and their
individual SI values, and
:py:meth:`~tespy.networks.network.Network.set_variable_value` imposes an SI
value through any of the represented variables.

.. code-block:: python

    nw.presolve("design")
    nw.print_variable_values()
    nw.set_variable_value("c1", "h", 2.8e6)
    nw.solve_continue()

The full solver information of the most recent :code:`solve` call, e.g. the
variable space, the equation lookups and the residual vector or the jacobian,
is held by a :py:class:`~tespy.solver.problem.Problem` instance, which is
accessible through the read-only :code:`Network.problem` property, e.g.
:code:`nw.problem.residual` or :code:`nw.problem.variables_dict`.

Pausing on a failed block
^^^^^^^^^^^^^^^^^^^^^^^^^
For non-converging models the block-wise solution process can pause at the
first block that fails instead of escalating to the coupled solution stages:

.. code-block:: python

    nw.solve("design", pause_on_block_failure=True)
    # status 20: paused, the log names the failed block, the cause of the
    # failure, its equations and its variables with their current values
    nw.print_variable_values(block=7)
    nw.print_block_jacobian(block=7)
    nw.print_block_states(block=7)
    nw.set_variable_value("c1", "h", 2.8e6)
    nw.solve_continue()

The already solved blocks stay fixed and the failed block's variables are
restored to their values at the block's entry and only the variables of the
failed block can be modified while paused. :code:`solve_continue` retries the
block with the modified values and continues the solution process, pausing
again in case the same or another block fails. Passing
:code:`pause_on_block_failure=False` there hands a still failing block over to
the standard escalation stages instead.

.. _networks_debugging_label:

Debugging
---------
In this section we show you how you can debug your models. Problems surface
in the same order as the stages of the solve process: the topology check, the
problem reduction, the determination check and the solution process itself.
Each stage has its own inspection tools. Two worked tutorials complement this
section: the :ref:`debugging tutorial <tutorial_debugging_label>` covers
specification errors and information extraction, the
:ref:`advanced debugging tutorial <tutorial_advanced_debugging_label>`
covers the solution process.

**Topology**

First, make sure your network topology is set up correctly, TESPy will prompt
an error, if it is not, and provide you with information, which components
are missing connections. Usually, this is the case, when you forgot to add
the connections to the network.

**Problem reduction**

In the first part of the presolving phase, the variable space reduction is
performed. TESPy will prompt errors, in case the parameter specifications in
context of the topology lead to an infeasibility in any of the variables.
This can be, for example

- a circular linear dependency between a set of variables. Typically, the
  mass flow can be over-determined by not including a :code:`CycleCloser`
  component in a circular network. For example, if you are modeling a cycle,
  e.g. the Clausius Rankine cycle, you need to make a cut in the cycle using
  the :code:`CycleCloser` or a :code:`Sink` and a :code:`Source` not to
  over-determine the system. Have a look in the
  :ref:`tutorial section <basics_label>` to understand why this is important
  and how it can be implemented.
- two parallel flows starting and ending in a common point (e.g. from a
  :code:`Splitter` to a :code:`Merge`) and both having linear specifications
  for the change of pressure from the start to the end. Then the common
  inflow and the common outflow pressure would be connected linearly through
  two different ways, which cannot be solved. One of both must be a result.
  Note: the same is of course true for a nonlinear dependency of pressure
  change, but this cannot be detected by the presolving.
- the values of two variables (or more) are directly specified in a set of
  linearly dependent variables. This does not need to be direct
  specification, it can also be indirect, through specifying temperature and
  vapor mass fraction in one location and specifying pressure in a different
  location while the specified pressure is linearly dependent to the
  pressure at the location with specified temperature and vapor mass
  fraction. In this case, the combination of temperature and vapor mass
  fraction determines the saturation pressure and therefore we end up with
  two pressure values fixed.

The error names the variables involved. To see, which variables form a
linearly dependent group and which specifications were consumed by the
reduction, use the
:ref:`inspection methods of the reduction <solver_reduction_label>`, e.g.
:code:`print_presolved_variables` and :code:`print_variables`.

**Determination check**

After the reduction is complete, a check will be carried out, if you specified
a sufficient number of parameters, meaning the exact number matching the number
of equations imposed to the problem. TESPy will prompt an error, if you did not
provide enough or if you provide too many parameters for your calculation. The
structural analysis helps locating the defects:
:py:meth:`~tespy.networks.network.Network.print_structural_analysis` names
the specific variables and equations of the under- or over-determined part
of the model and remains available after the error was raised. Since the
solver does not run in this case, solve with :code:`init_only=True` to use
the other inspection methods while you iterate on the specifications.

A specification error can also hide behind a matching parameter count: if
one part of the model is over-determined while another is under-determined,
the counts balance but no assignment of equations to variables exists.
Block-wise solving is not possible then. In case the structural defect stems
from an incomplete dependency declaration of a custom equation rather than from
the specifications, the simultaneous solution can still be attempted with
:code:`nw.solve(mode, block_solve=False)`.

.. note::

    Always keep in mind, that the system has to find a value for mass flow,
    pressure, enthalpy and the fluid mass fractions. Try to build up your
    network step by step. This way, you can easily determine, which parameters
    are still to be specified.

See the :ref:`debugging tutorial <tutorial_debugging_label>` for a worked
example on specification errors.

**Failures during solving**

With the correct number of specifications the solution process can still
fail, which the :code:`status` attribute reports:

- Status 2 - the iterations completed without convergence. The residual
  development can be inspected with :code:`nw.print_residuals()`. The reason
  often lies in variables driven to the boundaries of their validity or Newton
  producing nearly linear-dependent Jacobian matrices, which induce extremely
  small increments.
- Status 3 - a linear dependency in the Jacobian matrix. Either the values
  of your specifications over-determine one variable while missing out on
  another or the starting values are bad. Typically, this happens in
  combination with equations that require the **calculation of a temperature**,
  e.g. specifying a temperature at some point of the network with unknown
  pressure, or terminal temperature differences at heat exchangers. The
  singularity diagnosis names the suspect equations and variables in the log.
- Status 99 - the calculation crashed, most often because a fluid property
  cannot be evaluated. If this happens in the very first iterations, check
  whether you have specified all values in the correct unit.

To narrow such failures down, the following tools are available:

- The block-wise solution process localizes the failure: the log names the
  first block that did not converge together with its equations and
  variables.
- Passing :code:`pause_on_block_failure=True` freezes the solver in exactly
  that state, so you can inspect the block's residuals, Jacobian and fluid
  states, as shown in the section on
  :ref:`interacting with the solver <solver_interaction_label>`.
- The starting values can be inspected before any iteration is run:
  :code:`nw.presolve(mode)` followed by :code:`nw.print_variable_values()`
  shows the generated initial guesses. The most common root cause is a
  starting point in the wrong region of the fluid - liquid, two-phase or
  vapor - on connections involved in temperature equations. Correct such
  values through the :code:`val0` attributes or with
  :code:`set_variable_value` and continue with :code:`solve_continue`.

The :ref:`advanced debugging tutorial <tutorial_advanced_debugging_label>`
walks through this workflow on a concrete model.

Calculation speed improvement
-----------------------------
For improvement of calculation speed, the calculation of specific derivatives
is skipped, if the change of the corresponding variable was below a
threshold of :code:`1e-12` in the iteration before.

As a user you can take two more measures to improve calculation speed:
Specify primary variables whenever possible/reasonable. This will not only
reduce the variable space but also remove the necessity to calculate partial
derivatives towards them.

Finally, you can skip the postprocessing step.

- It will skip calculating all temperatures, volumetric flows, entropies,
  etc. of connections. Only the variables (mass flow, pressure, enthalpy and
  fluid composition will be available).
- It will skip calculation and evaluation of all component parameters and
  component parameter bounds. Errors in the results are not reported!

You can make use of this if you only want to extract certain results and if
you execute multiple simulations. In this case, you have to calculate the
respective results manually from the variables yourself and you have to make
sure that physical limitations are not violated.
