.. _tespy_development_new_comp_label:

How to contribute a new component
---------------------------------

What should I do if I need a component for a specific simulation that is not 
available in TESPy? One solution would be to request support on GitHub or to 
raise the issue at the next online or community meeting. But what good is that 
if I need the new component right now?

I will just implement it myself! The question becomes how to do this, since I 
don't know the code structure. These instructions on **how to contribute a 
new component** are intended to simplify the transition from user to an 
awesome developer.

Implementation of a polynomial compressor with cooling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The following basic principles should be in your mind when implementing a 
component:

* Visualize your task with a circuit diagram, equations and/or basic principle
* Adding the new equations and variables
* Expand the parameter definitions and derived calculations
* Write doc strings for the new component
* Implement tests for checking

.. tip::

    New component inherits from existing "parent" component

.. figure:: /_static/images/tutorials/compressor_with_cooling/drawing.svg
    :align: center
    :alt: Blackboard drawing of compressor with cooling
    :figclass: only-light

    Figure: Blackboard drawing of compressor with cooling

1. Implementing inputs and outputs
++++++++++++++++++++++++++++++++++

First, we have to create a new class and name it
::code::`PolynomialCompressorWithCooling`. Then, we have to implement the 
inputs and outputs. Based on our blackboard drawing we need two inputs as well 
as two ouputs.

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
    :start-after: [sec_6]
    :end-before: [sec_7]

After creating the inputs and outputs, a simple model should be built up to 
test the code. As usual, we will first create a network with the correct unit 
definitions. Then we will compile the components, including the new polynomial 
compressor with cooling and its connections. Once these have been added to the 
network, we will parameterise the system and solve it in the design mode. The 
correctness of the process can be confirmed by checking the mass flow and 
fluid composition in results.

.. dropdown:: Display source code for testing the model

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_1]
        :end-before: [sec_2]

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_11]
        :end-before: [sec_12]

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_12]
        :end-before: [sec_13]

2. Get mandatory constraints
++++++++++++++++++++++++++++

The next step is to get the mandatory constraints. To do this, a method is 
created that adds an additional constraint to the dictionary of mandatory 
constraints.

.. code-block:: python

    from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc

.. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
    :language: python
    :start-after: [sec_7]
    :end-before: [sec_8]

In this case, it retrieves the basic constraints of the upper class. A new 
constraint is added to ensure that the energy balance of the cooling system is 
met. Finally, the complete set of constraints is returned.

3. Define new equations
+++++++++++++++++++++++

After getting the mandatory constraint, the cooling energy balance is to be 
defined. First, the :code:`_preprocess()` method have to be defined. This 
performs validation before initiating the actual preprocessing of the 
component. It checks whether a required parameter (:code:`eta_recovery`) is 
set. If not, it throws an explanatory error (:code:`TESPyComponentError`). If 
so, it calls the default preprocessing of the base class.

.. note::

    During the first run, it may be helpful to ensure that the component 
    parameters are initialized. In this case, it could look like this:

    .. code-block:: python

        def _preprocess(self, row_idx):
            if not self.dissipation_ratio.is_set:
                self.dissipation_ratio.is_set = True
                self.dissipation_ratio.val = 0
                self.dissipation_ratio.val_SI = 0
            return super()._preprocess(row_idx)


.. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
    :language: python
    :start-after: [sec_8]
    :end-before: [sec_9]

In contrast, the cooling energy balance function and its dependencies are to 
be defined. The :code:`cooling_energy_balance_func()` method  describes the 
actual energy balance equation. It calculates the difference between the 
enthalpies of the incoming and outgoing mass flows, taking into account the 
energy loss (:code:`dissipation_ratio`) and a usable portion 
(:code:`usable_share`).

For dependencies, :code:`energy_balance_dependents()` determines which 
variables (mass flows and enthalpies of the inlets and outlets) appear in this 
balance equation and thus depend on the solver.

As in the first step, it is recommended to test the new code in between.

.. dropdown:: Display source code for testing the model

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_1]
        :end-before: [sec_2]

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_11]
        :end-before: [sec_12]

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_13]
        :end-before: [sec_14]

    .. error::

        .. code-block:: python

            You have provided too many parameters: 0 required, 1 supplied. Aborting calculation!
            ------------------------------------------------------------------------------------

            TESPyNetworkError                                 Traceback (most recent call last):

            Cell In[6], line 8
                  5 b2.set_attr(T=25, p=1)
                  6 compressor.set_attr(dissipation_ratio=0.1)
            ----> 8 nw.solve("design")

            File ~/gitprojects/tespy/src/tespy/networks/network.py:2486, in Network.solve(self, mode, init_path,
            design_path, max_iter, min_iter, init_only, init_previous, use_cuda, print_results, robust_relax)
            2483 msg = 'Starting solver.'
            2484 logger.info(msg)
            -> 2486 self.solve_determination()
            2488 try:
            2489     self.solve_loop(print_results=print_results)

            File ~/gitprojects/tespy/src/tespy/networks/network.py:2603, in Network.solve_determination(self)
            2601     logger.error(msg)
            2602     self.status = 12
            -> 2603     raise hlp.TESPyNetworkError(msg)
            2604 elif n < self.variable_counter:
            2605     msg = (
            2606         f"You have not provided enough parameters: {self.variable_counter} "
            2607         f"required, {n} supplied. Aborting calculation!"
            2608     )

            TESPyNetworkError: You have provided too many parameters: 0 
            required, 1 supplied. Aborting calculation!

    The test performed results in an error. Since we now have an additional 
    equation, this means that we have to define one less parameter, e.g., the 
    cooling mass flow.

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_14]
        :end-before: [sec_15]


4. Expand definitions of the parameters and add derived calculations
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Next, the new dimensions of the **polynomial compressor with cooling** are to 
be defined. To achieve this, a parameter definition for the efficiency of the 
heat recovery :code:`eta_recovery` and has to be determined.

.. note::

    In the next steps, the parameter definition is expanded.

.. code-block:: python

    from tespy.tools.data_containers import ComponentProperties as dc_cp
    from tespy.tools.helpers import TESPyComponentError

        def get_parameters(self):
            params = super().get_parameters()
            params["eta_recovery"] = dc_cp()
            return params


Once a again, it is recommended to test the code.

.. dropdown:: Display source code for testing the model

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
cooling output can be higher than temperature at gas output, as there is no 
heat exchanger balance equation implemented.

Further expansion of parameter definition:
Describe the expansion for :code:`td_minimal` and describe what :code:`min_val` is.

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
thermodynamic variables after the simulation. In this case, the internal 
maximum temperature in the compressor (:code:`T_max_compressor_internal`) and 
the minimum temperature difference between the compressor and the cooling 
medium (:code:`td_minimal`) are calculated using the outlet enthalpy of the 
compressor.

.. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
    :language: python
    :start-after: [sec_10]
    :end-before: [sec_11]

.. dropdown:: Display source code for testing the model

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_18]
        :end-before: [sec_19]

Last expansion of parameter definition:
Describe the expansion for :code:`dp_cooling` and which adjustable or 
calculated parameters the component has.

.. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
    :language: python
    :start-after: [sec_9]
    :end-before: [sec_10]

.. dropdown:: Display source code for testing the model

    .. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
        :language: python
        :start-after: [sec_19]
        :end-before: [sec_20]

5. Further tasks
++++++++++++++++
* doc strings
* Tests
    - After that: tests with assert_convergence()
    - Additional tests for different input parameter (edge) cases