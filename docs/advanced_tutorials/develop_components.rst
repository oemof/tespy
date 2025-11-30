.. _develop_components_tutorial_label:

How to develop a new component
==============================

What should I do if I need a component for a specific simulation that is not
available in TESPy? One solution would be to request support on GitHub or to
raise the issue at the next online or community meeting. But what good is that
if I need the new component right now?

I will just implement it myself! The question becomes how to do this, since I
don't know the code structure. These instructions on **how to develop a
new component** are intended to simplify the transition from user to an
awesome developer.

Implementation of a polynomial compressor with cooling
------------------------------------------------------
The following basic principles should be in your mind when implementing a
component:

* Visualize your task with a flowsheet, cycle diagram and/or the equations you
  want to implement
* Adding the new equations and variables
* Expand the parameter definitions and derived calculations
* Perform a lot of tests
* Write docstrings for the new component
* Implement tests to verify correctness

.. figure:: /_static/images/tutorials/compressor_with_cooling/drawing_pcwc.svg
    :align: center
    :alt: Blackboard drawing of compressor with cooling
    :figclass: only-light

    Figure: Blackboard drawing of compressor with cooling

.. figure:: /_static/images/tutorials/compressor_with_cooling/drawing_pcwc_darkmode.svg
    :align: center
    :alt: Blackboard drawing of compressor with cooling
    :figclass: only-dark

    Figure: Blackboard drawing of compressor with cooling

1. Implement the inputs and outputs
-----------------------------------

First, we create a new class named :code:`PolynomialCompressorWithCooling`.
Then, we have to implement the inputs and outputs. Based on our blackboard
drawing we need two inputs as well as two outputs.

.. note::

    The base component class has pairwise mass flow and fluid composition
    balance:

    * in1 matches with out1
    * in2 matches with out2

    This is automatically expanded with every new pair of ports.

.. attention::
    If you want to add ports with other names or non-paired ports, this may
    break.

.. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
    :language: python
    :start-after: [sec_3]
    :end-before: [sec_4]

After creating the inputs and outputs, a simple model should be built up to
test the code. As usual, we will first create a network with the correct unit
definitions. Then we will compile the components, including the new polynomial
compressor with cooling and its connections. Once these have been added to the
network, we will parameterise the system and solve it in the design mode. The
correctness of the process can be confirmed by checking the mass flow and
fluid composition in results. Since there is no connection between the working
fluid ports of the compressor and the cooling fluid yet, we have to provide
pressure, temperature and a mass flow.

.. dropdown:: Display source code for testing the model

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_1]
        :end-before: [sec_2]

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_9]
        :end-before: [sec_10]

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_10]
        :end-before: [sec_11]

2. Add mandatory constraints
----------------------------

The next step is to add the mandatory constraints. To do this, a method is
created that adds an additional constraint to the dictionary of mandatory
constraints.

.. code-block:: python

    from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc

.. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
    :language: python
    :start-after: [sec_4]
    :end-before: [sec_5]

In this case, it retrieves the basic constraints of the upper class. A new
constraint is added to ensure that the energy balance of the cooling system is
met. Finally, the complete set of constraints is returned.

3. Define new equations
-----------------------

Now we define the methods that are connected to the constraint: The function
returning the residual value of the equation and the list of variables the
equation depends on. The :code:`cooling_energy_balance_func()` method describes
the actual energy balance equation. It calculates the heat dissipated by the
working fluid side in the compressor through the value of the
:code:`dissipation_ratio`. On the cold side this heat should be added to the
cooling fluid, but not in its entirety, only an usable share of it. For now
we can hardcode that usable share in the equation with the variable
:code:`eta_recovery`.

.. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
    :language: python
    :start-after: [sec_5]
    :end-before: [sec_6]

As in the first step, it is recommended to test the new code in between.

.. dropdown:: Display source code for testing the model

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_1]
        :end-before: [sec_2]

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_9]
        :end-before: [sec_10]

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_11]
        :end-before: [sec_12]

    .. error::

        .. code-block:: bash

            You have provided too many parameters: 0 required, 1 supplied. Aborting calculation!

            TESPyNetworkError                                 Traceback (most recent call last):

            Cell In[6], line 8
             5 b2.set_attr(T=25, p=1)
             6 compressor.set_attr(dissipation_ratio=0.1)
             8 nw.solve("design")

            File ~/gitprojects/tespy/src/tespy/networks/network.py:2486, in Network.solve(self, mode, init_path,
            design_path, max_iter, min_iter, init_only, init_previous, use_cuda, print_results, robust_relax)
            2483 msg = 'Starting solver.'
            2484 logger.info(msg)
            2486 self.solve_determination()
            2488 try:
            2489     self.solve_loop(print_results=print_results)

            File ~/gitprojects/tespy/src/tespy/networks/network.py:2603, in Network.solve_determination(self)
            2601     logger.error(msg)
            2602     self.status = 12
            2603     raise hlp.TESPyNetworkError(msg)
            2604 elif n  self.variable_counter:
            2605     msg = (
            2606         f"You have not provided enough parameters: {self.variable_counter} "
            2607         f"required, {n} supplied. Aborting calculation!"
            2608     )

            TESPyNetworkError: You have provided too many parameters: 0
            required, 1 supplied. Aborting calculation!

    With no changes in our original specifications, the model results in an
    error. Since we now have an additional equation, this means that we have to
    specify one parameter less than before, e.g., the cooling mass flow.

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_12]
        :end-before: [sec_13]

4. Expand definitions of the parameters and add derived calculations
--------------------------------------------------------------------

Next, we can define parameters for the :code:`PolynomialCompressorWithCooling`.
We do this similar to the mandatory constraints by calling the
:code:`get_parameters` method, updating the dictionary and returning it. We
can start with the definition for the efficiency of the heat recovery
:code:`eta_recovery`.

For this, we use a data container :code:`dc_cp()`, which creates an object of
the :code:`ComponentProperties` class. These objects describe how TESPy should
handle a specific physical parameter, e.g., temperature, pressure loss,
efficiency, etc. Accordingly, TESPy uses these objects to automate unit
conversion, validation, equation integration and documentation.

.. code-block:: python

    from tespy.tools.data_containers import ComponentProperties as dc_cp
    from tespy.tools.helpers import TESPyComponentError


        def get_parameters(self):
            params = super().get_parameters()
            params["eta_recovery"] = dc_cp()
            return params

Along with the introduction of this parameter, we also update the equation.

.. code-block:: python

        def cooling_energy_balance_func(self):
            residual = (
                self.inl[1].m.val_SI * (self.outl[1].h.val_SI - self.inl[1].h.val_SI)
                + self.inl[0].m.val_SI * (
                    self.outl[0].h.val_SI
                    - self.outl[0].h.val_SI / (1 - self.dissipation_ratio.val_SI)
                    + self.inl[0].h.val_SI * (
                        self.dissipation_ratio.val_SI / (1 - self.dissipation_ratio.val_SI)
                    )
                ) * self.eta_recovery.val_SI
            )
            return residual

And, we can make the specification of :code:`eta_recovery` mandatory if we
want. This can be done by overriding the default :code:`_preprocess` method
like this:

.. code-block:: python

        def _preprocess(self, row_idx):
            if not self.eta_recovery.is_set:
                msg = (
                    f"The component {self.label} of type {self.__class__.__name__}"
                    "requires you to specify the share of heat recovery "
                    "eta_recovery."
                )
                raise TESPyComponentError(msg)

            return super()._preprocess(row_idx)

Once a again, it is recommended to test the code.

.. dropdown:: Display source code for testing the model

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_1]
        :end-before: [sec_2]

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_9]
        :end-before: [sec_10]

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_13]
        :end-before: [sec_14]

    We can also check if changing boundary conditions works and if the results
    seem reasonable:

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_14]
        :end-before: [sec_15]

    .. code-block:: python

        from tespy.tools.fluid_properties import T_mix_ph

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_15]
        :end-before: [sec_16]

The tests are going well. But with very small mass flows temperature at
cooling output can be higher than the temperature on the gas side, as there is
no limit implemented.

Further expansion of parameter definition
+++++++++++++++++++++++++++++++++++++++++

The next step is to define the parameter for the minimum temperature difference
:code:`td_minimal` between the compressor and the cooling medium. In this case,
the attribute :code:`min_val=0` means that this value must not be negative. If
it is, a warning is put out in the postprocessing automatically.

.. code-block:: python

    from tespy.tools.data_containers import ComponentProperties as dc_cp
    from tespy.tools.helpers import TESPyComponentError

        def get_parameters(self):
            params = super().get_parameters()
            params["eta_recovery"] = dc_cp()
            params["td_minimal"] = dc_cp(
                min_val=0
            )
            return params

In addition, :code:`calc_parameters()` is used to calculate the derived
properties after the simulation. In this case, the internal maximum temperature
in the compressor (:code:`T_max_compressor_internal`) and the minimum
temperature difference between the compressor and the cooling fluid
(:code:`td_minimal`) are calculated using the outlet enthalpy of the
compressor.

.. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
    :language: python
    :start-after: [sec_7]
    :end-before: [sec_8]

Further tests are being carried out at this point to examine the further
expansion of the parameter definition.

.. dropdown:: Display source code for testing the model

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_1]
        :end-before: [sec_2]

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_9]
        :end-before: [sec_10]

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_16]
        :end-before: [sec_17]

Final expansion of parameter definition
+++++++++++++++++++++++++++++++++++++++

To take the energy and pressure balance system of the component into account,
the parameter definition is extended one last time. For this purpose, the
pressure loss in the cooling circuit :code:`dp_cooling` is described. As in
the previous step, the value has not to be negative due to the
:code:`min_val=0` attribute. Instead of a :code:`func` we can use the
:code:`structure matrix` for this parameter. The :code:`structure_matrix` is
used in the preprocessing to make a linear connection between inlet and outlet
pressure in a way, that these two variables can then be mapped into a single
one in the numerical solution finding. The :code:`func_params` attribute
specifies the assignment for the internal calculation. Finally, the
:code:`quantity` specified indicates the physical unit.

.. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
    :language: python
    :start-after: [sec_6]
    :end-before: [sec_7]

Subsequent to the simulation, the calculation of the derived thermodynamic
variables is also required for the pressure loss :code:`dp_cooling` in the
cooling flow.

.. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
    :language: python
    :start-after: [sec_8]
    :end-before: [sec_9]

Our extension is also being tested here.

.. dropdown:: Display source code for testing the model

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_1]
        :end-before: [sec_2]

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_9]
        :end-before: [sec_10]

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_17]
        :end-before: [sec_18]

After checking that everything is correct, it's time to pat ourselves on the
back, because we have implemented a :code:`PolynomialCompressorWithCooling` in
TESPy.

5. Further tasks
----------------

Once implementation is complete, the hard work begins. Docstrings make your
code understandable, while tests ensure its reliability. If you want to
contribute your new component to tespy, then you can have a look at the
:ref:`developer guide <developing_label>`. More information on component
customization and implementation can also be found in
:ref:`this section <custom_components_label>`.
