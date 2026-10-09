.. _tutorial_starting_values_label:

How to Generate Stable Starting Values
======================================
Applying numerical algorithms and methods, the starting value of a variable is
the value used for the first iteration. TESPy generates starting values for all
variables automatically: known values propagate through the components,
temperature levels are estimated and phase consistent enthalpies and saturation
pressures are derived from them. The full process is described in the
:ref:`starting values section <solver_starting_values_label>` of the solver
documentation. With this machinery most models converge without any help from
the user.

With more complex models it can still happen, that the simulation does not
converge from the automatically generated values. In that case, the strategy is
to solve the model with a simpler but robust set of specifications first, for
example directly imposing the saturation temperature levels of a heat pump
cycle instead of terminal temperature differences. Since the solver starts from
the previous solution, the actual specifications can then be imposed in a
second solve: primary variables that hold a solved value do not need to be
guessed anymore.

Here we provide a short tutorial for you to better understand, how this process
could look like at the example of a subcritical heat pump with different
working fluids.

.. note::

    If the heat pump operates in trans- or supercritical range, some
    modifications have to be made on this setup. We plan to include respective
    examples here in the future.

You can download the full code of this example here:
:download:`starting_values.py </../tutorial/advanced/starting_values.py>`

Topology of the heat pump
-------------------------
Following the first tutorial a slightly different topology for a heat pump
with internal heat exchangers is considered instead of dumping the heat to the
ambient. You can see the plant topology in the figure below.

.. figure:: /_static/images/tutorials/heat_pump_starting_values/flowsheet.svg
    :align: center
    :alt: Topology of heat pump with internal heat exchanger
    :figclass: only-light

    Figure: Topology of heat pump with internal heat exchanger

.. figure:: /_static/images/tutorials/heat_pump_starting_values/flowsheet_darkmode.svg
    :align: center
    :alt: Topology of heat pump with internal heat exchanger
    :figclass: only-dark

    Figure: Topology of heat pump with internal heat exchanger

The system consists of a consumer system, a valve, an evaporator system, a
compressor and additionally an internal heat exchanger. In order to simulate
this heat pump, the TESPy model has to be built up. First, the network has to
be initialized, and the refrigerants used have to be specified. This example
shows how to make the heat pump model work with a variety of working fluids
with water on both the heat source and heat sink side of the system.

Running into errors
-------------------
As always, we start by importing the necessary TESPy classes.

.. literalinclude:: /../tutorial/advanced/starting_values.py
    :language: python
    :start-after: [sec_1]
    :end-before: [sec_2]

Then, we can build the network by defining components and connections. The
working fluid will be set with the variable `wf`, `"R290"` is used in the
first setup. This way, we will be able to change the working fluid in a flexible
way.

.. dropdown:: Click to expand to code section

    .. literalinclude:: /../tutorial/advanced/starting_values.py
        :language: python
        :start-after: [sec_2]
        :end-before: [sec_3]

After setting up the topology, the system's parameters should be set in the
following way:

- Heat sink temperature levels (`T` at 23 and 24)
- Heat source temperature levels (`T` at 11 and 13)
- Degree of overheating after the internal heat exchanger (`td_dew` at 2)
- Pinch point temperature difference at the evaporator (`ttd_l`) to derive
  evaporation pressure
- Temperature difference at the condenser (`ttd_u`) to derive condensation
  pressure
- Saturated gaseous state of the working fluid (`x=1`) after leaving the
  evaporator
- Efficiencies of pumps and the compressor (`eta_s`)
- Pressure losses in all heat exchangers (`pr1`, `pr2`, `pr`)
- Consumer heat demand (`Q`)

.. dropdown:: Click to expand to code section

    .. literalinclude:: /../tutorial/advanced/starting_values.py
        :language: python
        :start-after: [sec_3]
        :end-before: [sec_3_end]

The system should be well defined with the parameter settings, however the
solver does not find a solution from the automatically generated starting
values. It isolates the part of the problem it fails on and reports it:

.. error::

    .. code-block:: text

        Block 5 did not converge, solving the remaining 2 blocks simultaneously.
          Cause: no acceptance within the iteration budget of 50 iterations,
          the last scaled residual is 1.86e-02
          Equations: Compressor.eta_s,
          Heat Sink Condenser.energy_balance_constraints,
          Heat Source Evaporator.ttd_l,
          Internal Heat Exchanger.energy_balance_constraints, 1.x, 2.td_dew
          Variables: h0, h1, h4, m6, h7, p11
        The remaining system did not converge either, restarting with the
        simultaneous solution of the full system from its initial state.
        The solver does not seem to make any progress, aborting calculation.
        Scaled residual value is 1.86e-02 (1: x)
        Possible reasons include:
         - fluid properties moving outside the valid range of the property
           database (consider adjusting p_range or h_range),
         - an impossible constraint that can never be satisfied
         - bad starting values causing the Newton solver to diverge.
        Use nw.print_residuals() to identify which equations have the largest
        residuals.

The messages already narrow the failure down: the equation system was
:ref:`decomposed into blocks <solver_decomposition_label>` and every block
solved fine except one. The group of six coupled equations determining the
evaporation side of the cycle failed. Neither solving the remaining blocks in
one coupled system nor restarting with the full system finds the solution, so
the calculation ends with :code:`nw.status == 2`.

.. tip::

    To investigate such a failure interactively, the solver can pause at
    the failing block instead of escalating:

    .. code-block:: python

        nw.solve("design", pause_on_block_failure=True)
        nw.print_block_states(block=5, at="failure")
        nw.print_block_jacobian(block=5)

    The state tables show the fluid states the newton algorithm diverged to,
    the variable values of the failed block can be modified and the solution
    process continued with :code:`nw.solve_continue()`. The section on
    :ref:`interacting with the solver <solver_interaction_label>` describes the
    workflow in detail.

Fixing the errors
-----------------

All equations of the failed block are tied to the pressure and enthalpy levels
of the evaporation side of the cycle, so this is where better starting values
are required. To provide them, it is recommended to fix the saturation levels
of the cycle directly instead of the temperature differences in a first
calculation. In this example, the fixed points can be identified with the help
of the logph diagram which you can see in the figure below.

.. figure:: /_static/images/tutorials/heat_pump_starting_values/logph.svg
    :align: center
    :alt: Logph diagram of propane
    :figclass: only-light

    Figure: Logph diagram of propane

.. figure:: /_static/images/tutorials/heat_pump_starting_values/logph_darkmode.svg
    :align: center
    :alt: Logph diagram of propane
    :figclass: only-dark

    Figure: Logph diagram of propane

A rough estimation of the evaporation and condensation temperature can be
obtained from the temperature levels of the heat source and the heat sink.
Evaporation takes place a few Kelvin below the heat source backflow
temperature, condensation a few Kelvin above the consumer feed flow
temperature. These estimates can be imposed directly through the dew line
temperature :code:`T_dew` at the evaporator outlet and the bubble line
temperature :code:`T_bubble` at the condenser outlet. Each of them fixes the
pressure of its saturation state.

The terminal temperature differences are unset and the saturation temperatures
are set instead.

.. literalinclude:: /../tutorial/advanced/starting_values.py
    :language: python
    :start-after: [sec_4]
    :end-before: [sec_5]

The model was solved successfully and has stored the starting values for any
follow-up. Therefore, we can undo our recent changes and restart the
simulation. For example, the COP is then calculated.

.. literalinclude:: /../tutorial/advanced/starting_values.py
    :language: python
    :start-after: [sec_5]
    :end-before: [sec_6]

Expand fix to any working fluids
--------------------------------
Finally, using this strategy, it is possible to build a generic function,
building a network, that works with a variety of working fluids.

.. dropdown:: Click to expand to code section

    .. literalinclude:: /../tutorial/advanced/starting_values.py
        :language: python
        :start-after: [sec_6]
        :end-before: [sec_7]

We can run that function for different working fluids and plot the results:

.. literalinclude:: /../tutorial/advanced/starting_values.py
    :language: python
    :start-after: [sec_7]
    :end-before: [sec_8]

.. figure:: /_static/images/tutorials/heat_pump_starting_values/COP_by_wf.svg
    :align: center
    :alt: Analysis of the COP using different working fluids
    :figclass: only-light

    Figure: Analysis of the COP using different working fluids

.. figure:: /_static/images/tutorials/heat_pump_starting_values/COP_by_wf_darkmode.svg
    :align: center
    :alt: Analysis of the COP using different working fluids
    :figclass: only-dark

    Figure: Analysis of the COP using different working fluids
