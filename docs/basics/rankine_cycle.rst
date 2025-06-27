.. _tespy_basics_rankine_cycle_label:

Rankine Cycle
=============

.. figure:: /_static/images/basics/rankine_cycle.svg
    :align: center
    :alt: Topology of the rankine cycle
    :figclass: only-light

    Figure: Topology of the rankine cycle

.. figure:: /_static/images/basics/rankine_cycle_darkmode.svg
    :align: center
    :alt: Topology of the rankine cycle
    :figclass: only-dark

    Figure: Topology of the rankine cycle

Download the full script here:
:download:`rankine.py </../tutorial/basics/rankine.py>`

Setting up the Cycle
^^^^^^^^^^^^^^^^^^^^
We will model the cycle including the cooling water of the condenser. For this
start with the :code:`Network` set up we already know.

.. literalinclude:: /../tutorial/basics/rankine.py
    :language: python
    :start-after: [sec_1]
    :end-before: [sec_2]

Following, we create the components and connect them. The :code:`Condenser` has
a hot side inlet and outlet as well as a cold side inlet and outlet. The hot
side is indicated by using the index 1 for the inlet and outlet :code:`in1` and
:code:`out1`, the cold side uses the index 2 (:code:`in2` and :code:`out2`).

Again, for the closed thermodynamic cycle we have to insert a cycle closer. The
cooling water inlet and the cooling water outlet of the condenser are directly
connected to a :code:`Source` and a :code:`Sink` respectively.

.. literalinclude:: /../tutorial/basics/rankine.py
    :language: python
    :start-after: [sec_2]
    :end-before: [sec_3]

For the parameters, we predefine the pressure losses in the heat exchangers.
For the condenser, the hot side pressure losses are neglected :code:`pr1=1`,
for the cooling water side we assume pressure loss of 2 % :code:`pr2=0.98`. The
steam generator inflicts a pressure loss of 10 %.

The turbine and feed pump will have the isentropic efficiency specified. For
the connection parameters, the fluid has to be defined in both the main cycle
and the cooling water system. Furthermore, the live steam temperature, pressure
and mass flow are set. Lastly, we set the condensation pressure level and the
feed and return flow temperature of the cooling water as well as its feed
pressure.

.. literalinclude:: /../tutorial/basics/rankine.py
    :language: python
    :start-after: [sec_3]
    :end-before: [sec_4]

After running the simulation, for example, we can observe the temperature
differences at the condenser. Instead of directly setting a pressure value for
condensation, we could also set the upper terminal temperature difference
:code:`ttd_u` instead. It is defined as the condensation temperature to cooling
water return flow temperature.

.. tip::

    You will find the documentation of each equation of the components in the
    respective seciton of the API documentation. For example, the condenser
    :py:class:`tespy.components.heat_exchangers.condenser.Condenser`.

.. literalinclude:: /../tutorial/basics/rankine.py
    :language: python
    :start-after: [sec_4]
    :end-before: [sec_5]

After rerunning, we will see that the condensation temperature and pressure
are both automatically calculated by the specified terminal temperature value.

Generating T-s Diagram
^^^^^^^^^^^^^^^^^^^^^^
To visualize the Rankine cycle, we generate a temperature (T) versus entropy (s)
diagram using the fluprodia (Fluid Property Diagram) package.

.. dropdown:: Click to expand to code section

    .. literalinclude:: /../tutorial/basics/rankine.py
        :language: python
        :start-after: [sec_5]
        :end-before: [sec_6]

The steps involved in generating the T-s diagram are as follows:

- Import the Package: Import fluprodia and create an object by passing the
  alias of the fluid.
- Specify the Unit System: Set the unit system for all fluid properties.
- Specify Custom Isolines: Define custom isolines for the diagram.
- Calculate and draw isolines: Calculate and draw the background isolines.
- Calculate and draw process points and change of state
- Save and Export the Diagram: Save and export the completed T-s diagram.

.. figure:: /_static/images/basics/rankine_ts_diagram.svg
    :align: center
    :alt: T-s Diagram of Rankine Cycle
    :figclass: only-light

    Figure: T-s Diagram of Rankine Cycle

.. figure:: /_static/images/basics/rankine_ts_diagram_darkmode.svg
    :align: center
    :alt: T-s Diagram of Rankine Cycle
    :figclass: only-dark

    Figure: T-s Diagram of Rankine Cycle

Besides visualization, this feature is also useful for analysis purposes.
For example, if the T-s diagram forms a closed loop, validating the accuracy of
the model and that the operating fluid completes a successful Rankine Cycle. By
applying fluprodia, we can create and customize different types of diagrams for
all pure and pseudo-pure fluids available in CoolProp. For more information on
fluprodia, we refer users to the
`fluprodia documentation <https://fluprodia.readthedocs.io/en/latest/>`__.

Assess Electrical Power
^^^^^^^^^^^^^^^^^^^^^^^
To assess the electrical power output we want to consider the power generated
by the turbine as well as the power required to drive the feed pump. It is
possible to include both of the component's power values in a single electrical
:code:`PowerConnection`. We can do this by importing the respective component
classes, i.e.

- :py:class:`tespy.components.power.generator.Generator`
- :py:class:`tespy.components.power.motor.Motor`
- :py:class:`tespy.components.power.bus.PowerBus`
- :py:class:`tespy.components.power.sink.PowerSink`

as well as the :py:class:`tespy.connections.powerconnection.PowerConnection`
class. Then we create a :code:`PowerBus` for the electrical distribution, which
takes electricity fromt he :code:`Generator` and distributes it to the pump's
:code:`Motor` and the grid represented by the :code:`PowerSink`.

.. literalinclude:: /../tutorial/basics/rankine.py
    :language: python
    :start-after: [sec_6]
    :end-before: [sec_7]

The results for the :code:`PowerConnections` are printed separately. You can
also set the total desired power production of the system, for example
replacing the mass flow specification at connection 1:

.. literalinclude:: /../tutorial/basics/rankine.py
    :language: python
    :start-after: [sec_7]
    :end-before: [sec_8]

Analyze Efficiency and power generation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In this section, we will analyze the power production and the efficiency
of the cycle, given constant steam mass flow and with varying values for the

- live steam pressure,
- live steam temperature and
- cooling water temperature level.

To do that, we are using a very similar setup as has been used in the
:ref:`heat pump tutorial <tespy_basics_heat_pump_label>`. For the feed water
temperature level we want to set the change in temperature at the condenser
to a constant value. Also, we have to unset the power generation specification
again and use a constant mass flow instead. With :code:`iterinfo=False` we
can disable the printout of the convergence history.

.. dropdown:: Click to expand to code section

    .. literalinclude:: /../tutorial/basics/rankine.py
        :language: python
        :start-after: [sec_8]
        :end-before: [sec_9]

.. figure:: /_static/images/basics/rankine_parametric.svg
    :align: center
    :alt: Parametric analysis of the efficiency and power output
    :figclass: only-light

    Figure: Parametric analysis of the efficiency and power output

.. figure:: /_static/images/basics/rankine_parametric_darkmode.svg
    :align: center
    :alt: Parametric analysis of the efficiency and power output
    :figclass: only-dark

    Figure: Parametric analysis of the efficiency and power output

Part load Simulation
^^^^^^^^^^^^^^^^^^^^
In the part load simulation part, we are starting with a specific design of the
plant and calculate the part load performance with some assumptions on the
component's individual behavior. The table below summarizes the assumptions,
which we will keep as simple as possible at this moment. For more insights
have a look at the step by step
:ref:`heat pump tutorial <tespy_tutorial_heat_pump_label>` or at the
:ref:`Network documentation <tespy_modules_networks_label>`.

+-----------+---------------------------+-------------------------------+
| Component | Assumptions               | Settings                      |
+===========+===========================+===============================+
| Turbine   | cone law applies          | unset inlet pressure and      |
|           |                           | apply cone law                |
+-----------+---------------------------+-------------------------------+
| Condenser | constant heat transfer    | unset terminal temperature    |
|           | coefficient               | difference and set heat       |
|           |                           | transfer coefficient          |
+-----------+---------------------------+-------------------------------+
| Cooling   | constant volumetric flow  | unset return temperature      |
| water     |                           | value and set volumetric flow |
+-----------+---------------------------+-------------------------------+

With these specifications, the following physics are applied to the model:

- Due to the constant volumetric flow of water, the temperature of the cooling
  water returning from the condenser will react to the total heat transferred
  in the condensation: Increased heat transfer means increasing temperature,
  decreased heat transfer means decreased temperature.
- The constant heat transfer coefficient of the condenser will calculate the
  condensation temperature (and therefore pressure) based on the temperature
  regime in the cooling water side:

  - Increase in temperature for the cooling water leads to increased
    condensation temperature (at constant heat transfer).
  - Increase in heat transfer means increase in necessary temperature
    difference at the condenser (at constant cooling water inlet temperature).

- The cone law is a mathematical model to predict the pressure at the turbine's
  inlet based on the deviation from the design conditions. Generally,
  increased mass flow leads to higher inlet pressure (at constant inlet
  temperature and constant outlet pressure). However, this equation is more
  complex, since a lot more parameters are involved compared to the other
  equations applied.

In order to apply these specifications, we can use the :code:`design` and
:code:`offdesign` keywords. The keyword :code:`design` unsets the specified
data in the list in an offdesign calculation. The keyword :code:`offdesign`
automatically sets the respective parameter for the offdesign calculation. In
case the specification refers to a value, the value is taken from the design
mode calculation. In the example, the main condenser's kA value is calculated
in the design simulation and its value will be kept constant through the
offdesign simulations.

.. literalinclude:: /../tutorial/basics/rankine.py
    :language: python
    :start-after: [sec_9]
    :end-before: [sec_10]

We have to save the design state of the network and run the :code:`solve`
method with the :code:`design_path` specified.

.. literalinclude:: /../tutorial/basics/rankine.py
    :language: python
    :start-after: [sec_10]
    :end-before: [sec_11]

Finally, we can alter the mass flow from its design value of 20 kg/s to only
50 % of its value. In this example, we calculate the efficiency and plot it.

.. literalinclude:: /../tutorial/basics/rankine.py
    :language: python
    :start-after: [sec_11]
    :end-before: [sec_12]

.. figure:: /_static/images/basics/rankine_partload.svg
    :align: center
    :alt: Part load electric efficiency of the Rankine cycle
    :figclass: only-light

    Figure: Part load electric efficiency of the Rankine cycle

.. figure:: /_static/images/basics/rankine_partload_darkmode.svg
    :align: center
    :alt: Part load electric efficiency of the Rankine cycle
    :figclass: only-dark

    Figure: Part load electric efficiency of the Rankine cycle
