.. _tespy_tutorial_hrsg_label:

Heat Recovery Steam Generator for Process Steam Supply
------------------------------------------------------

We provide the full script presented in this tutorial here:
:download:`stepwise.py </../tutorial/advanced/hrsg.py>`

.. figure:: /_static/images/tutorials/hrsg_steam_supply/flowsheet_steps.svg
    :align: center
    :alt: Topology of the heat recovery steam generator
    :figclass: only-light

    Figure: Topology of the heat recovery steam generator

.. figure:: /_static/images/tutorials/hrsg_steam_supply/flowsheet_steps_darkmode.svg
    :align: center
    :alt: Topology of the heat recovery steam generator
    :figclass: only-dark

    Figure: Topology of the heat recovery steam generator

Overview
^^^^^^^^
In this tutorial we will build a combined heat and power plant utilizing flue
gas (e.g. from a gas turbine) to generate process steam and electricity. The
figure above shows the structure of the system and the 4 steps we take to build
the system.

...

- Table on boundary conditions
- Table on assumptions

....

The challenge here lies in the complexity of the system. With a well-structured
approach the problem can be tackled.

Set up the Network
^^^^^^^^^^^^^^^^^^

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_1]
    :end-before: [sec_2]

Modeling the heat pump: Consumer system
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: /_static/images/tutorials/heat_pump_stepwise/flowsheet_p1.svg
    :align: center
    :alt: First part of the system
    :figclass: only-light

    Figure: First part of the system

.. figure:: /_static/images/tutorials/heat_pump_stepwise/flowsheet_p1_darkmode.svg
    :align: center
    :alt: First part of the system
    :figclass: only-dark

    Figure: First part of the system

Components
++++++++++

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_2]
    :end-before: [sec_3]

Connections
+++++++++++

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_3]
    :end-before: [sec_4]

Parametrization & solving
+++++++++++++++++++++++++

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_5]
    :end-before: [sec_6]


More to follow
