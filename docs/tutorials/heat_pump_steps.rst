.. _tespy_tutorial_heat_pump_label:

Build a Heat Pump Stepwise
--------------------------

We provide the full script presented in this tutorial here:
:download:`stepwise.py </../tutorial/advanced/stepwise.py>`

.. figure:: /_static/images/tutorials/heat_pump_stepwise/flowsheet.svg
    :align: center
    :alt: Topology of the heat pump system
    :figclass: only-light

    Figure: Topology of the heat pump system

.. figure:: /_static/images/tutorials/heat_pump_stepwise/flowsheet_darkmode.svg
    :align: center
    :alt: Topology of the heat pump system
    :figclass: only-dark

    Figure: Topology of the heat pump system

Task
^^^^
This tutorial introduces you in how to model a heat pump in TESPy. You can see
the plants topology in the figure.

The main purpose of the heat pump is to deliver heat e.g. for the consumers of
a heating system. Thus, the heat pump's parameters will be set in a way, which
supports this target. Generally, if systems are getting more complex, it is
highly recommended to set up your plant in incremental steps. This tutorial
divides the plant in three sections: The consumer part, the valve and the
evaporator and the compressor as last element. Each new section will be
appended to the existing ones.

The system will be built up in a way, that independent of what working fluid
we use, we will be able to generate stable starting values. After achieving
that, the final parameters are specified instead of the initial values.

Set up a Network
^^^^^^^^^^^^^^^^
First, we have to create an instance of the
:py:class:`tespy.networks.network.Network` class. The network is the main
container of the model and will be required in all following sections. First,
it is necessary to specify a list of the fluids used in the plant. In this
example we will work with water (H\ :sub:`2`\O) and ammonia (NH\ :sub:`3`\).
Water is used for the cold side of the heat exchanger, for the consumer and
for the hot side of the environmental temperature. Ammonia is used as
refrigerant within the heat pump circuit. If you don’t specify the unit
system, the variables are set to SI-Units. We also keep the working fluid a
variable to make reusing the script with a different working fluid easy.

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_1]
    :end-before: [sec_2]

We will use °C, bar and kJ/kg as units for temperature and enthalpy.

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
We will start with the consumer as the plant will be designed to deliver a
specific heat flow. From the figure above you can determine the components of
the consumer system: condenser, pump and the consumer (
:py:class:`tespy.components.heat_exchangers.simple.HeatExchangerSimple`
). Additionally we need a source and a sink for the consumer and the heat pump
circuit respectively. We will import all necessary components already in the
first step, so the imports will not need further adjustment.

.. tip::

    We label the sink for the refrigerant :code:`"valve"`, as for our next
    calculation the valve will be attached there instead of the sink. In this
    way, the fluid properties can be initialized at the interface-connection,
    too.

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_2]
    :end-before: [sec_3]

Connections
+++++++++++
In the next steps we will connect the components in order to form the network.
Every connection requires the source, the source id, the target and the target
id as arguments: the source is the component from which the connection
originates, the source id is the outlet id of that component. This applies
analogously to the target. To find all inlet and outlet ids of a component look
up the class documentation of the respective component. An overview of the
components available and the class documentations is provided in the
:ref:`TESPy modules overview <tespy_modules_components_label>`. The
:py:class:`tespy.connections.connection.Ref` class is used specify fluid
property values by referencing fluid properties of different connections. It
is used in a later step.

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_3]
    :end-before: [sec_4]

.. seealso::

    As we learned in the :ref:`basic introduction <tespy_basics_intro_label>`
    instead of just connecting the consumer's outlet to the pump's inlet, we
    must make use of the CycleCloser. Closing a cycle without further
    adjustments will always result in a linear dependency in the fluid and the
    mass flow equations. The
    :py:class:`tespy.components.basics.cycle_closer.CycleCloser` component
    makes sure, the fluid properties pressure and enthalpy are identical at
    the inlet and the outlet. The component will prompt a warning, if the mass
    flow or the fluid composition at its outlet are different to those at its
    inlet. A different solution to this problem, is adding a merge and a
    splitter at some point of your network and connect the second inlet/outlet
    to a source/sink. This causes residual mass flow and residual fluids to
    emerge/drain there.

Parametrization
+++++++++++++++
For the condenser we set pressure ratios on hot and cold side. The consumer
will have pressure losses, too. Further we set the isentropic efficiency for
the pump. The most important parameter is the consumer's heat demand since it
decides the overall mass flow in the systems.

.. tip::

    In this tutorial we will first build the system with parameters that
    ensure stable starting values for a simulation, which in the end will be
    switched to reasonable values for the individual parts of the system. For
    example, instead of the evaporation pressure we will use the terminal
    temperature difference at the condenser instead.

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_4]
    :end-before: [sec_5]

In order to calculate this network further parametrization is necessary, as
e.g. the fluids are not determined yet: At the hot inlet of the condenser we
define the temperature, pressure and the fluid vector. A good guess for
pressure can be obtained from CoolProp's PropsSI function. We know that the
condensation temperature must be higher than the consumer's feed flow
temperature. Therefore we can set the pressure to a slightly higher value of
that temperature's corresponding condensation pressure.

The same needs to be done for the consumer cycle. We suggest to set
the parameters at the pump's inlet. On top, we assume that the consumer
requires a constant inlet temperature. The :code:`CycleCloser` automatically
makes sure, that the fluid's state at the consumer's outlet is the same as at
the pump's inlet.

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_5]
    :end-before: [sec_6]

Solve
+++++
After creating the system, we want to solve our network. Until we have not
set up the full system we will run design case calculations.

.. note::

    In TESPy there are two different types of calculations: design point and
    offdesign calculation.

    Generally, the design calculation is used for designing your system in the
    way you want it to look like. This means, that you might want to specify a
    design point isentropic efficiency, pressure loss or terminal temperature
    difference. After you have designed your system, you are able to make
    offdesign calculations with TESPy. The offdesign calculation is used to
    predict the system's behavior at different points of operation. For this
    case, this might be different ambient temperature, different feed flow
    temperature, or partial load. Add the end of this tutorial, you will learn
    how to run the offdesign calculation.

.. seealso::

    For general information on the solving process in TESPy and available
    parameters check the corresponding section in the
    :ref:`TESPy modules introduction <tespy_modules_networks_label>`.

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_6]
    :end-before: [sec_7]


Valve and evaporator system
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Next we will add the valve and the evaporator system to our existing network.
The figure below indicates the sections we will append in this step. This part
contains of a valve followed by an evaporator with a drum (separating
saturated liquid from saturated gas) and a superheater.

.. figure:: /_static/images/tutorials/heat_pump_stepwise/flowsheet_p2.svg
    :align: center
    :alt: Second part of the system
    :figclass: only-light

    Figure: Second part of the system

.. figure:: /_static/images/tutorials/heat_pump_stepwise/flowsheet_p2_darkmode.svg
    :align: center
    :alt: Second part of the system
    :figclass: only-dark

    Figure: Second part of the system

Components
++++++++++
First, we need to import the new components, which are
:py:class:`tespy.components.nodes.drum.Drum`,
:py:class:`tespy.components.heat_exchangers.base.HeatExchanger` and
:py:class:`tespy.components.piping.valve.Valve`. We will add these components
to the script.

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_7]
    :end-before: [sec_8]

Connections
+++++++++++
Since the old connection :code:`1` lead to a sink, we have to replace this
connection in the network. We can do that by using the method
:code:`del_conns` passing :code:`c1`. After that, we can create the new
connections and add them to the network as we did before.

The valve connects to the drum at the inlet :code:`'in1'`. The drum's outlet
:code:`'out1'` is saturated liquid and connects to the evaporator's cold side
inlet :code:`'in2'`. The inlet reconnects to the drum's inlet :code:`'in2'`.
The superheater's cold side connects to the drum's outlet :code:`'out2'`.
On the ambient side we simply connect the source to the superheater to the
evaporator and finally to the ambient sink. This will add the following
connections to the model:

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_8]
    :end-before: [sec_9]

.. attention::

    The drum is special component, it has an inbuilt CycleCloser, therefore
    although we are technically forming a cycle at the drum's outlet 1 to its
    inlet 2, we do not need to include a CycleCloser here.

Parametrization
+++++++++++++++
Previous parametrization stays untouched. Regarding the evaporator, we specify
pressure ratios on hot side as well as the evaporation pressure, for which we
can obtain a good initial guess based on the ambient temperature level using
CoolProp. From this specification the pinch point layout will be a result,
similar as in waste heat steam generators. The pressure ratio of the cold side
*MUST NOT* be specified in this setup as the drum imposes pressure equality
for all inlets and outlets.

The superheater will also use the pressure ratios on hot and cold side.
Further we set a value for the enthalpy at the working fluid side outlet. This
determines the degree of overheating and is again based on a good guess.

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_9]
    :end-before: [sec_10]

Next step is the connection parametrization: The pressure in the drum and the
enthalpy of the wet steam reentering the drum need to be determined. For the
enthalpy we can specify the vapor mass fraction :code:`x` determining the
degree of evaporation. On the hot side inlet of the superheater we define the
temperature, pressure and the fluid. At last we have to fully determine the
state of the incoming fluid at the superheater's hot side.

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_10]
    :end-before: [sec_11]

Solve
+++++
We can again run a simulation after adding these parts. This step is not
required, but in larger, more complex networks, it is recommended to achieve
better convergence.

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_11]
    :end-before: [sec_12]

Compressor system
^^^^^^^^^^^^^^^^^
To complete the heat pump, we will add the compressor system to our existing
network. This requires to change the connections 0, 6 and 17. The connection 6
has to be changed to include the compressor. After the last compressor stage,
connection 0 has to redefined, since we need to include the CycleCloser of the
working fluid's cycle. The connection 17 has to be connected to the heat
exchanger for intermittent cooling as well as the bypass.

.. figure:: /_static/images/tutorials/heat_pump_stepwise/flowsheet.svg
    :align: center
    :alt: Topology of the heat pump system
    :figclass: only-light

    Figure: Topology of the heat pump system

.. figure:: /_static/images/tutorials/heat_pump_stepwise/flowsheet_darkmode.svg
    :align: center
    :alt: Topology of the heat pump system
    :figclass: only-dark

    Figure: Topology of the heat pump system

Components
++++++++++
This part contains two compressors with intermittent cooling between them. The
cold side of the heat exchanger will be connected to a pump upstream and to
the superheater downstream. The bypass is used to give the system flexibility
in the temperature levels between the heat exchangers. We will also replace
the source for the refrigerant of :code:`c0` at the condenser with another
CycleCloser to make sure the fluid properties after the second compressor are
identical to the fluid properties at the condenser's inlet.

.. tip::

    The intermittent cooling extracts heat from the cycle. As this heat is
    however used to increase the evaporation pressure of the working fluid due
    to the higher temperature level of the heat source, the reduction is very
    limited. We use a two stage compressor, because in a single stage
    compression, the outlet temperature of the refrigerant might violate
    technical boundary conditions of the real-world component.

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_12]
    :end-before: [sec_13]

Connections
+++++++++++
We remove connections 0, 6 and 13 from the network, define the new connections
and add them again.

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_13]
    :end-before: [sec_14]

Parametrization
+++++++++++++++
For the first compressor we set the pressure ratio to the square root of the
full pressure ration between condensation and evaporation. In the first step,
we do not impose isentropic efficiency, because the respective equations are
quite sensitive to good starting value. We will set these values after the
full system has been calculated. The pump's isentropic efficiency value is not
as critical, therefore we set this value. The intermittent cooling imposes
pressure losses on both sides.

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_14]
    :end-before: [sec_15]

Regarding the connections we set enthalpy values for all working fluid side
connections. After the superheater and intermittent cooling the value will be
near saturation (enthalpy value of connection c5), after the compressors it
will be higher.

For the ambient side, we set temperature, pressure and fluid
at connection 11. On top of that, we can specify the temperature of the
ambient water after leaving the intermittent cooler.

With readding of connection 0 we have to set the fluid and the pressure again,
but not the temperature value, because this value will be a result of the
condensation pressure and the given enthalpy at the compressor's outlet.

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_15]
    :end-before: [sec_16]

Solve and Set Final System Parameters
+++++++++++++++++++++++++++++++++++++
Now we solve again. After that, we can exchange our good guesses with actual
useful parameters:

The condensation and evaporation pressure levels will be replaced by terminal
temperature values of the condenser and the evaporator respectively. The lower
terminal temperature value of the evaporator :code:`ttd_l` defines the pinch
point. The upper terminal temperature value :code:`ttd_u` of the condenser
defines the condensation pressure.

The degree of superheating in the superheater will be determined by the upper
terminal temperature instead of the enthalpy value at connection 6. The outlet
enthalpies after both compressors are replaced by the isentropic efficiency
values. Finally, the enthalpy after the intermittent cooling is replaced by
the temperature difference to the boiling point. With this we can ensure, the
working fluid does not start to condensate at the intermittent cooler.

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_16]
    :end-before: [sec_17]

Calculate Partload Performance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
After setting up the full system, we want to predict partload operation at
different values for the consumer's heat demand. Some of the values utilized
in the previous setup will change, if a component is not operated at its
design point. This is individual to every system, so the designer has to
answer the question: Which parameters are design point parameters and how does
the component perform at a different operation point.

.. tip::

    To make the necessary changes to the model, we can specify a design and an
    offdesign attribute, both lists containing component or connection
    parameter names. All parameters specified in the design attribute of a
    component or connection, will be unset in a offdesign calculation, all
    parameters specified in the offdesign attribute of a component or
    connection will be set for the offdesign calculation. The value for these
    parameters is the value derived from the design-calculation.

The changes we want to apply can be summarized as follows:

- All heat exchangers should be calculated based on their heat transfer
  coefficient with a characteristic for correction of that value depending
  on the change of mass flow (:code:`kA_char`). Therefore terminal temperature
  value specifications need to be added to the design parameters. Also, the
  temperature at connection 14 cannot be specified anymore, since it will be a
  result of the intermittent cooler's characteristics.
- Pumps and compressors will have a characteristic function for their
  isentropic efficiency instead of a constant value (:code:`eta_s_char`).
- Pressure drops in components will be a result of the changing mass flow
  through that component given the diameter in the design. The pressure ratio
  will therefore be replaced by :code:`zeta` for all heat exchangers. The zeta
  value is a geometry independent value.

On top of that, for the evaporator the characteristic function of the heat
transfer coefficient should follow different data than the default
characteristic. The name of that line is 'EVAPORATING FLUID' for the cold
side. The default line 'DEFAULT' will be kept for the hot side. These lines
are available in the :ref:`tespy.data <tespy_data_label>` module.

.. attention::

    If you run the offdesign simulation without any changes in the
    specification values, the results must be identical as in the respective
    design case! If they are not, it is likely, something went wrong.

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_17]
    :end-before: [sec_18]

.. seealso::

    If you want to learn more about handling characteristic functions you
    should have a glance at the
    :ref:`TESPy components section <tespy_modules_components_label>`.

Finally, we can change the heat demand and run several offdesign calculations
to calculate the partload COP value.

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_18]
    :end-before: [sec_19]

After successfully modeling the heat pump in design and offdesign cases, you
can now start using your model for further calculations. For example, if you
have a time series of required heat flow of your consumer, you can loop over
the series and perform offdesign calculation adjusting the heat flow every
time. Of course, this is possible with every offdesign parameter.

Have fun working with TESPy!
