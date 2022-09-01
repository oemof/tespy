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
the plants topology in the figure. Also, you will find a fully working model in
the last chapter of this tutorial.

The main purpose of the heat pump is to deliver heat e.g. for the consumers of
a heating system. Thus, the heat pump's parameters will be set in a way, which
supports this target.
Generally, if systems are getting more complex, it is highly recommended to set
up your plant in incremental steps. This tutorial divides the plant in three
sections: The consumer part, the valve and the evaporator and the compressor as
last element. Each new section will be appended to the existing ones.


Set up a Network
^^^^^^^^^^^^^^^^
First, we have to create an instance of the
:py:class:`tespy.networks.network.Network` class. The network is the main
container of the model and will be required in all following sections. First,
it is necessary to specify a list of the fluids used in the plant. In this
example we will work with water (H\ :sub:`2`\O) and ammonia (NH\ :sub:`3`\).
Water is used for the cold side of the heat exchanger, for the consumer and
for the hot side of the environmental temperature. Ammonia is used as coolant
within the heat pump circuit. If you don’t specify the unit system, the
variables are set to SI-Units.

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

    We label the sink for the coolant :code:`"valve"`, as for our next
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
For the condenser we set pressure ratios on hot and cold side and additionally
we set a value for the upper terminal temperature difference as design
parameter and the heat transfer coefficient as offdesign parameter. The
consumer will have pressure losses, too. Further we set the isentropic
efficiency for the pump, the offdesign efficiency is calculated with a
characteristic function. Thus, we set the efficiency as design parameter and
the characteristic function as offdesign parameter. In offdesign calculation
the consumer's pressure ratio will be a function of the mass flow, thus as
offdesign parameter we select zeta. The most important parameter is the
consumer's heat demand since it decides the overall mass flow in the systems.
We marked this setting as :code:`"key parameter"`.

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_4]
    :end-before: [sec_5]

In order to calculate this network further parametrization is necessary, as
e.g. the fluids are not determined yet: At the hot inlet of the condenser we
define the temperature and the fluid vector. In order to fully determine the
fluid's state at this point, an information on the pressure is required. This
is achieved by setting the terminal temperature difference (see above). The
same needs to be done for the consumer cycle. We suggest to set the parameters
at the pump's inlet. On top, we assume that the consumer requires a constant
inlet temperature. The :code:`CycleCloser` automatically makes sure, that the
fluid's state at the consumer's outlet is the same as at the pump's inlet.

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_5]
    :end-before: [sec_6]

.. note::

    In TESPy there are two different types of calculations: design point and
    offdesign calculation. All parameters specified in the design attribute of
    a component or connection, will be unset in a offdesign calculation, all
    parameters specified in the offdesign attribute of a component or
    connection will be set for the offdesign calculation. The value for these
    parameters is the value derived from the design-calculation.

    Generally, the design calculation is used for designing your system in the
    way you want it to look like. This means, that you might want to specify a
    design point isentropic efficiency, pressure loss or terminal temperature
    difference. After you have designed your system, you are able to make
    offdesign calculations with TESPy. The offdesign calculation is used to
    predict the system's behavior at different points of operation. For this
    case, this might be different ambient temperature, different feed flow
    temperature, or partial load.

Solve
+++++
After creating the system, we want to solve our network. First, we calculate
the design case and directly after we can perform the offdesign calculation.
Before finishing the full network, we do this without changing the value of
any parameters to test, if the calculation actually works. All results have
to be identical to the design case results, if no value was modified.

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
pressure ratios on hot side as well as the lower terminal temperature
difference. This specification is similar to pinch point layout for waste heat
steam generators. The pressure ratio of the cold side *MUST NOT* be specified
in this setup as the drum imposes pressure equality for all inlets and outlets.

We use the hot side pressure ratio and the lower terminal temperature
difference as design parameters and choose zeta as well as a characteristic
function referring to the area independent heat transfer coefficient as its
offdesign parameters.

On top of that, the characteristic function of the evaporator should follow
the default characteristic line of 'EVAPORATING FLUID' on the cold side and
the default line 'DEFAULT' on the hot side. These lines are defined in the
:ref:`tespy_data_label`.

The superheater will also use the pressure ratios on hot and cold side.
Further we set a value for the upper terminal temperature difference. For
offdesign and design parameter specification the same logic as for the
evaporator and the already existing part of the network is applied. The system
designer has to answer the question: Which parameters are design point
parameters and how does the component perform at a different operation point.

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_9]
    :end-before: [sec_10]

.. seealso::

    If you want to learn more about handling characteristic functions you
    should have a glance at the
    :ref:`TESPy components section <tespy_modules_components_label>`.

Next step is the connection parametrization: The pressure in the drum and the
enthalpy of the wet steam reentering the drum need to be determined. For the
enthalpy we can specify a reference of the circulating mass flow to the main
cycle mass flow. The pressure is achieved through the given lower terminal
temperature difference of the evaporator and its hot side outlet temperature.
As we have specified a terminal temperature difference at the evaporator's
cold side inlet (:code:`ttd_l`), it might be necessary to state a starting
value for the pressure or the state of the fluid (gaseous), as we are near to
the two-phase region. On the hot side inlet of the superheater we define the
temperature, pressure and the fluid. Since the pressure between superheater
and first compressor will be a result of the pressure losses in the
superheater and we set the terminal temperature difference there, bad starting
values will lead to a linear dependency, as a temperature and a pressure are
set while the fluid's state could be within the two phase region. Thus, we
choose to specify :code:`state="g"`, so the solver will keep the fluid in
gaseous state at all times. At last we have to fully determine the state of
the incoming fluid at the superheater's hot side.

.. attention::

    Do only use the :code:`state` keyword if you know the fluid's state prior
    to the simulation. If you specify the fluid to be gaseous but the correct
    result of the simulation would be within the two-phase region, your
    calculation most likely will not converge.

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
network. This requires to change the connections 0, 6 and 13. The connection 6
has to be changed to include the compressor. After the last compressor stage,
connection 0 has to redefined, since we need to include the CycleCloser of the
working fluid's cycle. The connection 13 has to be connected to the heat
exchanger for intermittent cooling and finally with a pump to make the water
from the ambient heat source flow through the heat exchangers.

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
the superheater downstream. We will also replace the source for the coolant
of :code:`c0` at the condenser with another CycleCloser (:code:`cool_closer`),
to make sure the fluid properties after the second compressor are identical to
the fluid properties at the condenser's inlet.

.. tip::

    The intermittent cooling extracts heat from the cycle. As this heat is
    however used to increase the evaporation pressure of the working fluid due
    to the higher temperature level of the heat source, the reduction is very
    limited. We use a two stage compressor, because in a single stage
    compression, the outlet temperature of the coolant might violate technical
    boundary conditions of the real-world component.

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_12]
    :end-before: [sec_13]

Connections
+++++++++++
Consequently to the addition of the cycle closer we have to adjust the
connection definition touching the new cycle closer.

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_13]
    :end-before: [sec_14]

Parametrization
+++++++++++++++
For the two compressor we defined an isentropic efficiency and for the
offdesign calculation a generic characteristic line for the isentropic
efficiency will be applied. The first compressor has a fixed pressure ratio,
the seconds compressor pressure ratio will result from the required pressure
at the condenser. The heat exchanger comes with pressure ratios on both sides.
The parametrization of all other components remains identical.

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_14]
    :end-before: [sec_15]

.. code-block:: python

    cp1.set_attr(eta_s=0.8, design=['eta_s'], offdesign=['eta_s_char'])
    cp2.set_attr(eta_s=0.8, pr=5, design=['eta_s'], offdesign=['eta_s_char'])
    he.set_attr(pr1=0.99, pr2=0.98, design=['pr1', 'pr2'],
                offdesign=['zeta1', 'zeta2', 'kA_char'])

Regarding the connections, on the hot side after the intercooler we set the
temperature. For the cold side of the heat exchanger we set the temperature,
the pressure and the fluid on the inlet flow, at the outlet we specify the
temperature as a design parameter. In offdesign calculation, this will be a
result from the given heat transfer coefficient (see parametrisation of
intercooler, kA_char is an offdesign parameter). Last, make sure the fluid
properties after the compressor outlet are identical to those at the condenser
inlet using the references.

The last step leads to a necessary redefinition of the parametrization of the
existing model: As the enthalpy at the outlet of the second compressor is a
result of the given pressure ratio and the isentropic efficiency, it is not
allowed to set the temperature at the condenser's hot inlet anymore.

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_15]
    :end-before: [sec_16]

.. code-block:: python

    # condenser system

    c_in_cd.set_attr(fluid={'water': 0, 'NH3': 1})

    # compressor-system

    he_cp2.set_attr(T=40, p0=10)
    ic_in_he.set_attr(p=5, T=20, fluid={'water': 1, 'NH3': 0})
    he_ic_out.set_attr(T=30, design=['T'])

Solve
+++++

Here again, using the saved results from previous calculations is always
favorable, but with manually adjusted starting values and the :code:`state`
specifier, the calculation should still converge. If you want to use the
previous part to initialise start the solver with

.. literalinclude:: /../tutorial/advanced/stepwise.py
    :language: python
    :start-after: [sec_16]
    :end-before: [sec_17]

.. code-block:: python

    nw.solve('design', init_path='condenser')


Explore Partload Performance
++++++++++++++++++++++++++++

Further tasks
^^^^^^^^^^^^^
After successfully modeling the heat pump in design and offdesign cases, you
can now start using your model for further calculations. For example, if you
have a time series of required heat flow of your consumer, you can loop over
the series and perform offdesign calculation adjusting the heat flow every
time. Of course, this is possible with every offdesign parameter.

Have fun working with TESPy!
