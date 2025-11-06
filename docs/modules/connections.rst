.. _tespy_modules_connections_label:

Connections
===========

This section provides an overview of the parametrisation of connections, how to
use referencing of variables in different locations of your problem and an
overview on the :ref:`PowerConnection <tespy_powerconnections_label>`,
which is a connection to represent non-material energy flows.

Parametrisation
---------------

As mentioned in the introduction, for each connection you can specify the
following parameters:

* mass flow :code:`m`,
* volumetric flow :code:`v`,
* pressure :code:`p`,
* enthalpy :code:`h`,
* temperature :code:`T`,
* temperature difference to bubble line :code:`td_bubble`
* vapor mass fraction for pure fluids :code:`x`,
* a fluid vector :code:`fluid` and
* a balance closer for the fluid vector :code:`fluid_balance`.

Setting and unsetting values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A :code:`Connection` always connects two components, we will create a simple
problem using a :code:`SimpleHeatExchanger` to showcase specification options.

.. code-block:: python

    >>> from tespy.networks import Network
    >>> from tespy.connections import Connection
    >>> from tespy.components import Sink, Source, SimpleHeatExchanger
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(
    ...     temperature="°C", power="kW", pressure="bar", enthalpy="kJ/kg"
    ... )
    >>> source = Source('source')
    >>> heatexchanger = SimpleHeatExchanger('heat exchanger')
    >>> sink = Sink('sink')
    >>> c1 = Connection(source, "out1", heatexchanger, "in1", label="c1")
    >>> c2 = Connection(heatexchanger, "out1", sink, "in1", label="c2")
    >>> nw.add_conns(c1, c2)

It is possible to set simple values and then solve, e.g.:

- mass flow, pressure and temperature (and fluid) at inlet
- enthalpy at outlet
- (pressure drop in heat exchanger)

.. code-block:: python

    >>> c1.set_attr(fluid={"R290": 1}, m=5, p=3, T=50)
    >>> c2.set_attr(h=500)
    >>> heatexchanger.set_attr(dp=0)
    >>> nw.solve("design")

Both connections will have all results available, these can be accessed

- as SI value
- as value in the network's default unit of the respective quantity
- as value with the corresponding quantity (:code:`pint.Quantity`)

The results include mass flow, pressure, enthalpy, temperature, temperature,
vapor quality, specific volume :code:`vol` and volumetric flow.

.. code-block:: python

    >>> round(c1.h.val, 1)  # value in kJ/kg but without unit attached
    668.9
    >>> round(c1.vol.val_with_unit, 3)  # will be in m3/kg
    <Quantity(0.195, 'm3 / kilogram')>
    >>> round(c2.T.val_SI, 1)  # SI value
    259.0
    >>> round(c2.T.val_with_unit, 1)  # will be in °C
    <Quantity(-14.2, 'degree_Celsius')>

You can also provide quantities to a specific parameter to individually specify
a unit to a parameter, e.g. inlet mass flow. Note, that units are retained when
set with individual quantity.

.. code-block:: python

    >>> Q = nw.units.ureg.Quantity
    >>> c1.set_attr(m=Q(2, "t/h"))
    >>> nw.solve("design")
    >>> c1.m.val_with_unit
    <Quantity(2, 'metric_ton / hour')>
    >>> round(c2.m.val_with_unit, 2)
    <Quantity(0.56, 'kilogram / second')>

For pure fluids or CoolProp/REFPROP mixtures we can also specify two-phase
properties:

- vapor mass fraction/quality :code:`x`
- dew line temperature difference for superheating :code:`td_dew`
- bubble line temperature difference for subcooling :code:`td_bubble`

We can replace the inlet temperature specification e.g. with superheating. The
unit of :code:`temperature_difference` is different from the unit for
:code:`temperature`. Unsetting a value is simple: Just set it to :code:`None`.

.. code-block:: python

    >>> c1.set_attr(T=None)  # unset the value
    >>> nw.units.default["temperature"]
    '°C'
    >>> nw.units.default["temperature_difference"]
    'delta_degC'
    >>> c1.set_attr(td_dew=20)
    >>> nw.solve("design")
    >>> round(c1.T.val, 2)
    5.82

Setting starting values
^^^^^^^^^^^^^^^^^^^^^^^

Setting starting values for the variables can be helpful in some situations.
You can do this for the following properties:

- mass flow
- pressure
- enthalpy

.. code-block:: python

    >>> c1.set_attr(m0=4)
    >>> c2.set_attr(h0=300, p0=4)

Linear relationships
^^^^^^^^^^^^^^^^^^^^

It is also possible to set up linear relationships between different
specifications in your system in the form of:

.. math::

    x_0 = a * x_1 + b

It is possible to specify these for

- mass flow, pressure and enthalpy as well as
- temperature and volumetric flow.

For example, instead of h we can specify a reference to the temperature at
c1. The factor is always based on SI value, the delta is in the default unit of
the respective property. The starting value for h is required in this context
because in the previous calculation the fluid was in two-phase state, meaning
the partial derivative of the temperature of c2 with respect to enthalpy would
be zero otherwise.

.. code-block:: python

    >>> from tespy.connections import Ref
    >>> factor = 1
    >>> delta = 25
    >>> c2.set_attr(h=None, T=Ref(c1, factor, delta), h0=1000)
    >>> nw.solve("design")
    >>> round(c1.T.val_SI * factor + delta - 273.15, 2)
    30.82
    >>> round(c2.T.val, 2)
    30.82

Instead we could reference volumetric flow at outlet to find the temperature at
outlet.

.. code-block:: python

    >>> c2.set_attr(T=None)
    >>> factor = 1.2
    >>> delta = 0
    >>> c2.set_attr(v=Ref(c1, factor, delta))
    >>> nw.solve("design")
    >>> round(c1.v.val_SI * factor + delta, 2)
    0.11
    >>> round(c2.v.val_SI, 2)
    0.11
    >>> round(c2.T.val, 2)
    52.91

For more complex (and arbitrary) relationships between variables of the system
use the :code:`UserDefinedEquation` class. Some examples can be found in
:ref:`this section <tespy_ude_label>`.

Fluid specification
^^^^^^^^^^^^^^^^^^^

This sections shows some details on the specification of fluids.

.. code-block:: python

    # set both elements of the fluid vector
    >>> c1.set_attr(fluid={'water': 1})

    # same thing, but using data container
    >>> c1.fluid.set_attr(_val={'water': 1}, _is_set={'water'})
    >>> c1.fluid.is_set
    {'water'}

    # set starting values
    >>> c1.set_attr(fluid0={'water': 1})

    # unset full fluid vector
    >>> c1.set_attr(fluid={'water': None})
    >>> c1.fluid.is_set
    set()

    # unset part of fluid vector
    >>> c1.set_attr(fluid={'N2': 0.7, "O2": 0.3})
    >>> c1.fluid.is_set.remove('N2')
    >>> c1.fluid.is_set
    {'O2'}

CoolProp and REFPROP
++++++++++++++++++++

It is possible to specify the fluid property back end of the fluids by adding
the name of the back end in front of the fluid's name. For incompressible binary
mixtures, you can append the water volume/mass fraction to the fluid's name, for
example:

.. code-block:: python

    >>> c1.set_attr(fluid={'water': 1})  # HEOS back end
    >>> c1.set_attr(fluid={'INCOMP::water': 1})  # incompressible fluid
    >>> c1.set_attr(fluid={'BICUBIC::air': 1})  # bicubic back end
    >>> c1.set_attr(fluid={'INCOMP::MPG[0.5]|mass': 1})  # binary incompressible mixture

You can also specify REFPROP based fluids, e.g. R513A, which is a mass based
mixture of R134a and R1234yf:

.. code-block:: python

    >>> c1.set_attr(fluid={'REFPROP::R134A[0.44]&R1234yf[0.56]|mass': 1})  # REFPROP back end

.. note::

    Without further specifications CoolProp will be used as fluid property
    database. If you do not specify a back end, the **default back end**
    :code:`HEOS` will be used. For an overview of the back ends available please
    refer to the :ref:`fluid property section <tespy_fluid_properties_label>`.

Other Backends
++++++++++++++

You can also change the engine, for example to the iapws library. It is even
possible, that you define your own custom engine, e.g. using polynomial
equations. Please check out the fluid properties' section in the docs on how to
do this.

.. code-block:: python

    >>> from tespy.tools.fluid_properties.wrappers import IAPWSWrapper
    >>> c1.set_attr(fluid={'H2O': 1}, fluid_engines={"H2O": IAPWSWrapper})

Please also check out the section on
:ref:`custom fluid properties <tespy_fluid_properties_label>` for more
information.

Access from the :code:`Network` object
--------------------------------------

You may want to access the network's connections other than using the variable
names, for example in an imported network or connections from a subsystem. It
is possible to access these using the connection's label. By default, the label
is generated by this logic:

:code:`source:source_id_target:target_id`, where

- :code:`source` and :code:`target` are the  labels of the components that are
  connected.
- :code:`source_id` and :code:`target_id` are e.g. :code:`out1` and
  :code:`in2` respectively.

.. code-block:: python

    >>> conn = nw.get_conn('c1')
    >>> conn.label
    'c1'
    >>> conn.set_attr(p=7)
    >>> conn.p.val
    7.0

.. note::

    The label can only be specified on creation of the connection. Changing the
    label after might break this access method.

.. _tespy_powerconnections_label:

PowerConnections
================

PowerConnections can be used to represent non-material energy flow, like power
or heat. You can make use of generators, motors and buses.

Different use-cases for the implementation of :code:`PowerConnection` with the
respective power components can be:

- apply motor or generator efficiencies
- connect multiple turbomachines on a single shaft
- collect all electricity production and own consumption to calculate net
  power

The handling of the :code:`PowerConnection` and the respective components is
identical to standard components. The following components are available:

- :py:class:`tespy.components.power.generator.Generator`: generate electricity from mechanical energy
- :py:class:`tespy.components.power.motor.Motor`: generate mechanical energy from electricity
- :py:class:`tespy.components.power.bus.PowerBus`: balance all inflows and outflows of power into a bus
- :py:class:`tespy.components.power.sink.PowerSink`: e.g. represent power fed into the electricity grid
- :py:class:`tespy.components.power.source.PowerSource`: e.g. represent power drawn from the electricity grid

For more details on the components please go to the respective section of the
:ref:`documentaton <tespy_modules_components_label>` and the respective API
documentation linked in the list above.

Parameters
----------

The :code:`PowerConnection` only holds a single parameter, namely the power
flow :code:`E` (:math:`\dot E`), which is measured in Watts. You can create a
:code:`PowerConnection` instance by connecting to a component that has a
respective inlet or outlet. For example, consider a turbine generating
electricity. First we can set up a system as we are used to do without any
:code:`PowerConnections`:

.. code-block:: python

    >>> from tespy.components import Source, Sink, Turbine
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(temperature="degC", pressure="bar")
    >>> so = Source("source")
    >>> turbine = Turbine("turbine")
    >>> si = Sink("sink")
    >>> c1 = Connection(so, "out1", turbine, "in1", label="c1")
    >>> c2 = Connection(turbine, "out1", si, "in1", label="c2")
    >>> nw.add_conns(c1, c2)

We can parametrize the model, e.g. consider the turbine part of a gas turbine,
which expands hot flue gases:

.. code-block:: python

    >>> c1.set_attr(fluid={"air": 1}, p=10, T=1000, m=1)
    >>> c2.set_attr(p=1)
    >>> turbine.set_attr(eta_s=0.9)
    >>> nw.solve("design")
    >>> round(turbine.P.val / 1e3)
    -577

We can add a connection between the turbine and the grid. This will add one
extra variable to our problem (the energy flow :code:`E`) but also one extra
equation, namely the turbine energy balance. Therefore, after adding the new
connection, there is nothing to change to make the model solve.

.. code-block:: python

    >>> from tespy.connections import PowerConnection
    >>> from tespy.components import PowerSink
    >>> grid = PowerSink("grid")
    >>> e1 = PowerConnection(turbine, "power", grid, "power", label="e1")
    >>> nw.add_conns(e1)
    >>> nw.solve("design")
    >>> round(e1.E.val / 1e3)
    577

.. note::

    Note that the value of the energy flow of a :code:`PowerConnection` will
    always be positive in the defined direction (from one component's outlet
    to another component's inlet).

To learn what power connections are available in each of the component classes
see the respective API documentation. Below you will find more examples
utilizing the :code:`PowerConnection`.

Examples
--------

Single shaft gas turbine
^^^^^^^^^^^^^^^^^^^^^^^^

To make a more elaborate example, we will implement an open gas turbine
system using air as working fluid and a heater. You can also model gas
turbines with combustion, for this example, the focus is on modeling the
single shaft gas turbine system.

First, we import the necessary components and set up the material flow
system connecting the compressor to the heater and to the turbine.

.. code-block:: python

    >>> from tespy.connections import Connection, PowerConnection
    >>> from tespy.components import (
    ...     Turbine, Source, Sink, Compressor, SimpleHeatExchanger, PowerBus,
    ...     Generator, PowerSink
    ... )
    >>> from tespy.networks import Network
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(temperature="degC", pressure="bar")
    >>> so = Source("source")
    >>> heater = SimpleHeatExchanger("heater")
    >>> compressor = Compressor("compressor")
    >>> turbine = Turbine("turbine")
    >>> si = Sink("sink")
    >>> c1 = Connection(so, "out1", compressor, "in1", label="c1")
    >>> c2 = Connection(compressor, "out1", heater, "in1", label="c2")
    >>> c3 = Connection(heater, "out1", turbine, "in1", label="c3")
    >>> c4 = Connection(turbine, "out1", si, "in1", label="c4")

Next, we can set up the energy flows. Since the turbine and the compressor
are rotating on the same shaft, the turbine powers the compressor and the
generator at the same time. For this, we can use a PowerBus to represent
the shaft, which gets powered by the turbine. The turbine's power connector
is connected to the 'power_in1' connector of the shaft. The shaft
connects to the compressor's connector 'power' and to the generator's
connector 'power_in'. The generator then is connected at its outlet
'power_out' to the grid representation at the connector 'power'.

.. code-block:: python

    >>> shaft = PowerBus("shaft", num_in=1, num_out=2)
    >>> generator = Generator("generator")
    >>> grid = PowerSink("grid")
    >>> e1 = PowerConnection(turbine, "power", shaft, "power_in1", label="e1")
    >>> e2 = PowerConnection(shaft, "power_out1", compressor, "power", label="e2")
    >>> e3 = PowerConnection(shaft, "power_out2", generator, "power_in", label="e3")
    >>> e4 = PowerConnection(generator, "power_out", grid, "power", label="e4")
    >>> nw.add_conns(c1, c2, c3, c4, e1, e2, e3, e4)

We can parametrize the system, in this example, we fix the ambient air
temperature, pressure and mass flow, the turbine inlet temperature and the
turbine outlet pressure (to be equal to the ambient pressure).

.. code-block:: python

    >>> c1.set_attr(fluid={"air": 1}, m=1, p=1, T=25)
    >>> c3.set_attr(T=1000)
    >>> c4.set_attr(p=1)

In the components the turbine's and the compressor's efficiency are set as
well as the compressor's pressure ratio and the heater's pressure drop.

.. code-block:: python

    >>> turbine.set_attr(eta_s=0.9)
    >>> compressor.set_attr(eta_s=0.9, pr=15)
    >>> heater.set_attr(dp=0)

With the four power connections we have four additional variables in our
system. The compressor, the turbine and the shaft all deliver one equation
(their energy balance equation), meaning, one parameter is missing to fully
set up our problem. This could be the generator efficiency. With that, we
can solve the system, and check what amount of electricity is generated
through the generator.

.. code-block:: python

    >>> generator.set_attr(eta=0.95)
    >>> nw.solve("design")
    >>> round(e4.E.val_SI / 1e3, 1)
    246.9

Alternatively, we could also fix the electricity output to a specific
target value and unset the air mass flow. This will calculate the required
air mass flow to generate the desired amount of electricity.

.. code-block:: python

    >>> e4.set_attr(E=3e5)
    >>> c1.set_attr(m=None)
    >>> nw.solve("design")
    >>> round(c1.m.val, 3)
    1.215

Single shaft feed water pump powered by a turbine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a pump that is powered by a turbine. The turbine's :code:`turbine_fwp`
power output must therefore be equal to the pump's :code:`fwp` power
consumption.

.. code-block:: python

    >>> from tespy.networks import Network
    >>> from tespy.components import Pump, Turbine, Source, Sink
    >>> from tespy.connections import Connection, PowerConnection

    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(temperature="degC", pressure="bar")
    >>> cond = Source("condensate")
    >>> fwp = Pump("feed water pump")
    >>> feedwater = Sink("feedwater")
    >>> c1 = Connection(cond, "out1", fwp, "in1")
    >>> c2 = Connection(fwp, "out1", feedwater, "in1")
    >>> ls = Source("live steam")
    >>> turbine_fwp = Turbine("turbine fwp")
    >>> ws = Sink("Waste steam")
    >>> c11 = Connection(ls, "out1", turbine_fwp, "in1")
    >>> c12 = Connection(turbine_fwp, "out1", ws, "in1")
    >>> e1 = PowerConnection(turbine_fwp, "power", fwp, "power")
    >>> nw.add_conns(c1, c2, c11, c12, e1)

We can set up the system in a way, that calculates the required mass flow
of steam through the turbine to power the feed water pump and find the
power flow by accessing the respective attribute of the power connection.

.. code-block:: python

    >>> c1.set_attr(fluid={"water": 1}, p=0.5, x=0, m=10)
    >>> c2.set_attr(p=50)
    >>> fwp.set_attr(eta_s=0.75)
    >>> c11.set_attr(fluid={"water": 1}, p=40, T=500)
    >>> c12.set_attr(p=0.55)
    >>> turbine_fwp.set_attr(eta_s=0.9)
    >>> nw.solve("design")
    >>> nw.assert_convergence()
    >>> round(e1.E.val_SI / 1e3)
    68

Logic to force same power of two compressors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this example we combine the PowerConnection with a UserDefinedEquation.
Two air compressors should run in series and at identical power. For this the
intermediate pressure is variable.

    >>> from tespy.components import Source, Sink, Compressor, PowerSource
    >>> from tespy.connections import Connection, PowerConnection
    >>> from tespy.networks import Network
    >>> from tespy.tools import UserDefinedEquation
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(temperature="degC", pressure="bar")
    >>> so = Source("air source")
    >>> compressor1 = Compressor("compressor 1")
    >>> compressor2 = Compressor("compressor 2")
    >>> si = Sink("compressed air")
    >>> grid1 = PowerSource("grid compressor 1")
    >>> grid2 = PowerSource("grid compressor 2")
    >>> c1 = Connection(so, "out1", compressor1, "in1", label="c1")
    >>> c2 = Connection(compressor1, "out1", compressor2, "in1", label="c2")
    >>> c3 = Connection(compressor2, "out1", si, "in1", label="c3")
    >>> e1 = PowerConnection(grid1, "power", compressor1, "power")
    >>> e2 = PowerConnection(grid2, "power", compressor2, "power")
    >>> nw.add_conns(c1, c2, c3, e1, e2)
    >>> c1.set_attr(fluid={"air": 1}, m=1, p=1, T=25)
    >>> c3.set_attr(p=5)
    >>> compressor1.set_attr(eta_s=0.85)
    >>> compressor2.set_attr(eta_s=0.85)
    >>> def same_power_ude(ude):
    ...     e1, e2 = ude.conns
    ...     return e1.E.val_SI - e2.E.val_SI
    >>> def same_power_dependents(ude):
    ...     e1, e2 = ude.conns
    ...     return [c.E for c in ude.conns]
    >>> ude = UserDefinedEquation(
    ...     "power equality ude",
    ...     func=same_power_ude,
    ...     dependents=same_power_dependents,
    ...     conns=[e1, e2]
    ... )
    >>> nw.add_ude(ude)
    >>> nw.solve("design")
    >>> nw.assert_convergence()
    >>> round(e1.E.val / 1e3) == round(e2.E.val / 1e3)
    True
    >>> round(e1.E.val / 1e3)
    105

Including part load model for motor efficiency
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This example how a partload efficiency curve can be applied to a motor. For
this, let's assume the motor powers a refrigeration compressor. We can set up
the model by connecting the compressor to the refrigerant flows as usual and
add the :code:`PowerConnection` to the electricity grid via a :code:`Motor`
instance.

.. code-block:: python

    >>> from tespy.components import Source, Sink, Compressor, PowerSource, Motor
    >>> from tespy.connections import Connection, PowerConnection
    >>> from tespy.networks import Network
    >>> from tespy.tools import CharLine

    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(temperature="degC", pressure="bar")
    >>> so = Source("evaporated refrigerant")
    >>> compressor = Compressor("compressor")
    >>> si = Sink("compressed refrigerant")
    >>> grid = PowerSource("grid")
    >>> motor = Motor("motor")
    >>> c1 = Connection(so, "out1", compressor, "in1", label="c1")
    >>> c2 = Connection(compressor, "out1", si, "in1", label="c2")
    >>> e1 = PowerConnection(grid, "power", motor, "power_in", label="e1")
    >>> e2 = PowerConnection(motor, "power_out", compressor, "power", label="e2")
    >>> nw.add_conns(c1, c2, e1, e2)

The design efficiency is 0.98, the compressor's design efficiency is 0.85. On
top we fix the inlet state and mass flow as well as the compressor's pressure
ratio. For the characteristics of the motor's efficiency we can pass data to a
:code:`CharLine` instance, which is set to be used for the :code:`eta_char`
method in the model of the motor.

.. code-block:: python

    >>> c1.set_attr(fluid={"R290": 1}, m=1, td_dew=10, T=10)
    >>> compressor.set_attr(pr=3, eta_s=0.85, design=["eta_s"], offdesign=["eta_s_char"])
    >>> motor.set_attr(eta_char=CharLine(x=[0.5, 0.75, 1, 1.25], y=[0.9, 0.975, 1, 0.975]))
    >>> motor.set_attr(eta=0.98, design=["eta"], offdesign=["eta_char"])
    >>> nw.solve("design")
    >>> nw.save("design.json")
    >>> nw.assert_convergence()

After performing the design simulation we can change the fluid mass flow and
observe the change in efficiency of the motor:

.. code-block:: python

    >>> c1.set_attr(m=0.8)
    >>> nw.solve("offdesign", design_path="design.json", init_path="design.json")
    >>> nw.assert_convergence()
    >>> round(motor.eta.val, 3)
    0.966

.. note::

    As mentioned in the component section: It is also possible to import your
    custom characteristics from the :code:`HOME/.tespy/data` folder. Read more
    about this :ref:`here <tespy_modules_characteristics_label>`.
