.. _tespy_development_new_comp_label:

How to contribute a new component
---------------------------------

What should I do if I need a component for a specific simulation that is not 
available in TESPy? One solution would be to request support on GitHub or to 
raise the issue at the next online or community meeting. But what good is that 
if I need the new component right now?

I will just implement it myself! The question then becomes how to do this, 
since I don't know the code structure. These instructions on 'how to 
contribute a new component' are intended to simplify the transition from user 
to a awesome developer.

Implementation of a polynomial compressor with cooling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The following basic principles should be followed when implementing a 
component:

* Visualize your task with a circuit diagram, equations and/or basic principle
* Adding the new variables
* Adding the new equations
* New component inherits from existing "parent" component

.. figure:: /_static/images/tutorials/compressor_with_cooling/drawing.svg
    :align: center
    :alt: Blackboard drawing of compressor with cooling
    :figclass: only-light

    Figure: Blackboard drawing of compressor with cooling

1. Implementing inputs and outputs and easy model for testing
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

First, we have to create a new class and name it
'PolynomialCompressorWithCooling'. Then, we have to implement the inputs 
and outputs. Based on our design we need two inputs as well as two ouputs.

.. note::

    The base component class has pairwise mass flow and fluid composition 
    balance:

    * in1 matches with out1

    * in2 matches with out2

    This is automatically expanded with every new pair of ports.

.. error::
    If you want to add ports with other names or non-paired ports, this may
    break. 

.. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
    :language: python
    :start-after: [sec_1]
    :end-before: [sec_2]

After creating the inputs and outputs, a simple model should be built up to 
test the code. As usual, we will first create a network with the correct unit 
definitions. Then we will compile the components, including the new polynomial 
compressor with cooling and its connections. Once these have been added to the 
network, we will parameterise the system and solve it in the design mode. The 
correctness of the process can be confirmed by checking the mass flow.

.. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
    :language: python
    :start-after: [sec_2]
    :end-before: [sec_3]

2. Get mandatory constraints
++++++++++++++++++++++++++++

The next step is to get the mandatory constraints. To do this, a method is 
created that adds an additional constraint to the dictionary of mandatory 
constraints.

.. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
    :language: python
    :start-after: [sec_3]
    :end-before: [sec_4]

In this case, it retrieves the basic constraints of the upper class. A new 
constraint is added to ensure that the energy balance of the cooling system is 
met. Finally, the complete set of constraints is returned.

3. Define new equations
+++++++++++++++++++++++

After getting the mandatory constraint, the cooling energy balance is to be 
defined. In a first step, the :code:`_preprocess()` method ensures that the 
:code:`dissipation_ratio` parameter (the proportion of dissipated energy) is 
initialized during the first run.

.. note::

    However, it was rejected later because we decided that it would be better 
    to set it explicitly than to have it happen in the background.

.. literalinclude:: /../tutorial/advanced/compressor_with_cooling.py
    :language: python
    :start-after: [sec_4]
    :end-before: [sec_5]

In contrast, the cooling energy balance function and its dependencies are to 
be defined. The :code:`cooling_energy_balance_func()` method  describes the 
actual energy balance equation: It calculates the difference between the 
enthalpies of the incoming and outgoing flows, taking into account the energy 
loss (:code:`dissipation_ratio`) and a usable portion (:code:`usable_share`).

For dependencies, energy_balance_dependents()` determines which variables
(mass flows and enthalpies of the inlets and outlets) appear in this balance 
equation and thus depend on the solver.

4. Test the new equations
+++++++++++++++++++++++++
..  --> set one less parameter

5. Data container
+++++++++++++++++
.. BUT: At very small massflows, T_cooling_out can be higher than T_gas_out, as 
.. there is no heat exchanger balance equation implemented

6. Structure matrics
++++++++++++++++++++

7. Further tasks
++++++++++++++++
* doc strings
* Tests
    - After that: tests with assert_convergence()
    - Additional tests for different input parameter (edge) cases