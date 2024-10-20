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

It is a simplified version of an actual plant with some unique features:

- The economizer bypass (connection 28) actually represents a performance
  downgrade for the boiler. The technical reason is to allow it to produce less
  than 30 % of the nominal steam mass flow. With this measure the CHP plant
  can produce almost any steam flow between 40 t/h and 210 t/h, which is a
  unusual wide range.
- The real economizer has eight stages. To simplify, we only model there here.
- Interestingly, the boiler is operated with what could be called a "negative
  approach" (the temperature coming out of the economizer was higher than that
  of the HP drum). This issue is managed by increasing the economizer pressure.
  Although this solution is inefficient from an energy standpoint, it was
  chosen for operational safety, likely to reduce the risk of economizer
  steaming.
- The sweetwater condenser is a solution for HP steam desuperheating.

...

- Table on boundary conditions
- Table on assumptions

....

The challenge here lies in the complexity of the system. With a well-structured
approach the problem can be tackled.

Steam supply and turbine section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: /_static/images/tutorials/hrsg_steam_supply/flowsheet_p1.svg
    :align: center
    :alt: First part of the system
    :figclass: only-light

    Figure: First part of the system

.. figure:: /_static/images/tutorials/hrsg_steam_supply/flowsheet_p1_darkmode.svg
    :align: center
    :alt: First part of the system
    :figclass: only-dark

    Figure: First part of the system

In the first part we will model the steam supply part and the turbine. We start
by importing the :code:`Network` class and creating an instance of it.

.. literalinclude:: /../tutorial/advanced/hrsg.py
    :language: python
    :start-after: [sec_1]
    :end-before: [sec_2]

Components and Connections
++++++++++++++++++++++++++

For the components, we need to import the respective classes as shown in this
part of the flowsheet and create the respective instances.

.. literalinclude:: /../tutorial/advanced/hrsg.py
    :language: python
    :start-after: [sec_2]
    :end-before: [sec_3]

The same applies for the connections. The extraction steam turbine is connected
to the live steam source. At a later stage, this connection will be replaced by
one that directly links the high pressure superheater 2. Similarly, the
connection 15 connects to a sink, which will be replaced by the low pressure
feed pump later.

.. literalinclude:: /../tutorial/advanced/hrsg.py
    :language: python
    :start-after: [sec_3]
    :end-before: [sec_4]

Parametrization & solving
+++++++++++++++++++++++++

For the live steam at the turbine inlet we have to assume a steam mass flow. In
the completed model it will be governed from the flue gas mass flow. On top of
that, we need temperature and pressure information. For the steam supply a
pressure of 13 bar and 10 °C of superheating is required. There also is a
specific mass flow to be provided, which however cannot be directly as this
specification would cause convergence issues. Instead, massflow 3 is fixed
initially. The outlet of the condenser is saturated liquid at a temperature of
35 °C. The condensate pump pumps the pressure back to ambient pressure. It is
assumed that 75 % of the process steam leaving the plant returns in the form of
condensate. The remaining mass flow is replaced by the makeup water.

.. literalinclude:: /../tutorial/advanced/hrsg.py
    :language: python
    :start-after: [sec_4]
    :end-before: [sec_5]

For the components we can specify the isentropic efficiencies of the turbines
and the pumps as well as the pressure ratio of the condensing turbine inlet
valve as well as the condenser pressure loss.

.. literalinclude:: /../tutorial/advanced/hrsg.py
    :language: python
    :start-after: [sec_5]
    :end-before: [sec_6]


More to follow
