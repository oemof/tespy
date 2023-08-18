.. _tespy_basics_district_heating_label:

District Heating Network
========================

.. figure:: /_static/images/basics/district_heating.svg
    :align: center
    :alt: Topology of the district heating network
    :figclass: only-light

    Figure: Topology of the district heating network

.. figure:: /_static/images/basics/district_heating_darkmode.svg
    :align: center
    :alt: Topology of the district heating network
    :figclass: only-dark

    Figure: Topology of the district heating network

The model used in this example is shown the figure. It consists of a central
heating plant and a consumer, represented by a heat exchanger with a control
valve. A much more complex district heating system is included in the
advanced tutorials section.

Download the full script here:
:download:`district_heating.py </../tutorial/basics/district_heating.py>`

Setting up the System
^^^^^^^^^^^^^^^^^^^^^
For this model we have to import the :code:`Network` and :code:`Connection`
classes as well as the respective components. After setting up the network we
can create the components, connect them to the network (as shown in the other)
examples. As a fluid, we will use the incompressibles back-end of CoolProp,
since we only need liquid water. The incompressible back-end has much higher
access speed while preserving high accuracy.


.. tip::

    For more information on the fluid properties in TESPy,
    :ref:`check out this page <tespy_fluid_properties_label>`.

.. dropdown:: Click to expand to code section

    .. literalinclude:: /../tutorial/basics/district_heating.py
        :language: python
        :start-after: [sec_1]
        :end-before: [sec_2]

In the first step, we assume we have a specific heat demand of the consumer
and constant pressure and thermal losses in the pipes. Furthermore, the pump
produces a constant pressure at the feed part of the system. With the control
valve in place the pressure of the return part of the system is then decoupled
from that value. Therefore, we need to set a pressure value at the sink as
well, which should be equal to the pressure at the pump's inlet. The pressure
drop in the valve will then be the residual pressure drop between the feed and
the return part of the system. Lastly, we fix the feed flow and the return
flow (at connection 4) temperature values.

.. literalinclude:: /../tutorial/basics/district_heating.py
    :language: python
    :start-after: [sec_2]
    :end-before: [sec_3]

Design Pipe Dimensions
^^^^^^^^^^^^^^^^^^^^^^
In the second step we will design the pipe's dimensions. There are two tasks
for this:

- Calculate the necessary pipe diameter given a target pressure loss as well
  as length and pipe roughness.
- Calculate the necessary insulation of the pipe based on assumptions
  regarding the heat loss at a given ambient temperature value.

For the first step, we set lengths and roughness of the pipe and the diameter
to :code:`"var"`, indicating the diameter of the pipe should be a variable
value in the calculation.

.. literalinclude:: /../tutorial/basics/district_heating.py
    :language: python
    :start-after: [sec_3]
    :end-before: [sec_4]

In the second step we can fix the diameter to its resulting value and
therefore unset the desired pressure loss first. Then, we set the ambient
temperature of the pipes (we assume the temperature of the ambient is not
affected by the heat loss of the pipe). With the given heat loss, the
:code:`kA` value can be calculated. It is the area independent heat transfer
coefficient.

.. literalinclude:: /../tutorial/basics/district_heating.py
    :language: python
    :start-after: [sec_4]
    :end-before: [sec_5]

.. note::

    In the results you can see, that the pipes' pressure losses are still at
    the desired value after remove the pressure ration specification and using
    the calculated value of the diameter instead.

Changing Operation Conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Next, we want to investigate what happens, in case the

- ambient temperature changes.
- heat load varies.
- overall temperature level in the heating system is reduced.

To do that, we will use similar setups as show in the Rankine cycle
introduction. The :code:`KA` value of both pipes is assumed to be fixed, the
efficiency of the pump and pressure losses in consumer and heat source are
constant as well.

.. dropdown:: Click to expand to code section

    .. literalinclude:: /../tutorial/basics/district_heating.py
        :language: python
        :start-after: [sec_5]
        :end-before: [sec_6]

.. figure:: /_static/images/basics/district_heating_partload.svg
    :align: center
    :alt: Performance of the district heating system at changing operating conditions
    :figclass: only-light

    Figure: Performance of the district heating system at changing operating conditions.

.. figure:: /_static/images/basics/district_heating_partload_darkmode.svg
    :align: center
    :alt: Performance of the district heating system at changing operating conditions
    :figclass: only-dark

    Figure: Performance of the district heating system at changing operating conditions.

.. note::

    The efficiency value is defined as ratio of the heat delivered to the
    consumer to the heat production in the central heating plant.

    .. math::

        \eta = \frac{\dot{Q}_\text{consumer}}{\dot{Q}_\text{production}}
