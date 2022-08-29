.. _tespy_tutorial_label:

~~~~~~~~~~~~~~~~~~
Advanced Tutorials
~~~~~~~~~~~~~~~~~~
We provide more advanced tutorials for you to better understand how to work
with more complex systems in TESPy.

At the example of different heat pump topologies, you will learn to

- create a more complex model *step by step* and get the idea of designing a
  plant and calculating the offdesign behavior.
- set up a code structure, which allows you to generate stable starting values
  flexibly, helping you to make faster analyses.
- use the inbuilt exergy analysis method in a simple geothermal heat pump
  setting.

Furthermore, we introduce the coupling of TESPy with pygmo in order to create
an optimization problem, which optimizes thermal efficiency of a clausius
rankine power plant.

- CGAM!
- Larger dh system

.. toctree::
    :maxdepth: 1
    :glob:
    :hidden:

    tutorials/heat_pump_steps.rst
    tutorials/starting_values.rst
    tutorials/heat_pump_exergy.rst
    tutorials/pygmo_optimization.rst
    tutorials/combustion_chamber.rst
