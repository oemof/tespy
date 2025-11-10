.. _advanced_tutorials_label:

~~~~~~~~~~~~~~~~~~
Advanced Tutorials
~~~~~~~~~~~~~~~~~~
We provide more advanced tutorials for you to better understand how to work
with more complex systems in TESPy.

The first tutorial highlights the different heat exchanger type classes
available in the framework and how they differentiate to help you selecting the
appropriate one.

At the example of different heat pump topologies, you will learn to

- create a more complex model *step by step* and get the idea of designing a
  plant and calculating the offdesign behavior.
- set up a code structure, which allows you to generate stable starting values
  flexibly, helping you to build reliable setups faster.
- use the inbuilt exergy analysis method in a simple geothermal heat pump
  setting.

Furthermore, we introduce the coupling of TESPy with pygmo in order to create
an optimization problem, which optimizes thermal efficiency of a Clausius
Rankine power plant. Also, there is a tutorial on the implementation of an air
source heat pump in an energy system dispatch optimization problem using
`oemof-solph <https://github.com/oemof/oemof-solph>`__. For that tutorial an air
source heat pump is implemented in various details of modeling complexity and
the results are transferred to the energy system optimization problem.

If you have any questions, ideas for other tutorials or feedback, please reach
out to us. We are looking forward to hearing from you!

.. grid:: 2
    :gutter: 1

    .. grid-item-card:: Overview on available classes for heat exchangers
        :link: tutorial_heat_exchanger
        :link-type: ref

        .. image:: /_static/images/tutorials/heat_exchangers/SectionedHeatExchanger_vs_HeatExchanger.svg
            :class: only-light

        .. image:: /_static/images/tutorials/heat_exchangers/SectionedHeatExchanger_vs_HeatExchanger_darkmode.svg
            :class: only-dark

    .. grid-item-card::  Build complex systems step by step
        :link: tutorial_heat_pump_label
        :link-type: ref

        .. image:: /_static/images/tutorials/heat_pump_stepwise/flowsheet.svg
            :class: only-light

        .. image:: /_static/images/tutorials/heat_pump_stepwise/flowsheet_darkmode.svg
            :class: only-dark

    .. grid-item-card::  Generate stable starting values
        :link: tutorial_starting_values_label
        :link-type: ref

        .. image:: /_static/images/tutorials/heat_pump_starting_values/COP_by_wf.svg
            :class: only-light

        .. image:: /_static/images/tutorials/heat_pump_starting_values/COP_by_wf_darkmode.svg
            :class: only-dark


    .. grid-item-card::  Exergy analysis of a heat pump
        :link: tutorial_heat_pump_exergy_label
        :link-type: ref

        .. image:: /_static/images/tutorials/heat_pump_exergy/diagram_E_D.svg
            :class: only-light

        .. image:: /_static/images/tutorials/heat_pump_exergy/diagram_E_D_darkmode.svg
            :class: only-dark

.. toctree::
    :maxdepth: 1
    :glob:
    :hidden:

    advanced_tutorials/heat_exchangers.rst
    advanced_tutorials/heat_pump_steps.rst
    advanced_tutorials/starting_values.rst
    advanced_tutorials/heat_pump_exergy.rst
