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

* mass flow* (m),
* volumetric flow (v),
* pressure* (p),
* enthalpy* (h),
* temperature* (T),
* vapor mass fraction for pure fluids (x),
* a fluid vector (fluid) and
* a balance closer for the fluid vector (fluid_balance).

It is possible to specify values, starting values and references. The data
containers for connections are `dc_prop` for fluid properties (mass flow,
pressure, enthalpy, temperature, etc.) and `dc_flu` for fluid composition. If
you want to specify information directly to these do it with caution, data types
are not enforced.

In order to create the connections we create the components to connect first.

.. code-block:: python

    >>> from tespy.connections import Connection, Ref
    >>> from tespy.components import Sink, Source

    # create components
    >>> source1 = Source('source 1')
    >>> source2 = Source('source 2')
    >>> sink1 = Sink('sink 1')
    >>> sink2 = Sink('sink 2')

    # create connections
    >>> myconn = Connection(source1, 'out1', sink1, 'in1')
    >>> myotherconn = Connection(source2, 'out1', sink2, 'in1')

    # set pressure and vapor mass fraction by value, temperature and enthalpy
    # analogously
    >>> myconn.set_attr(p=7, x=0.5)
    >>> myconn.p.val, myconn.x.val
    (7, 0.5)

    # set starting values for mass flow, pressure and enthalpy (has no effect
    # on temperature and vapor mass fraction!)
    >>> myconn.set_attr(m0=10, p0=15, h0=100)
    >>> myconn.m.val0, myconn.p.val0, myconn.h.val0
    (10, 15, 100)

    # do the same directly on the data containers
    >>> myconn.p.set_attr(val=7, is_set=True)
    >>> myconn.x.set_attr(val=0.5, is_set=True)

    >>> myconn.m.set_attr(val0=10)
    >>> myconn.p.set_attr(val0=15)
    >>> myconn.h.set_attr(val0=100)
    >>> myconn.m.val0, myconn.p.val0, myconn.h.val0
    (10, 15, 100)

    # specify a referenced value: pressure of myconn is 1.2 times pressure at
    # myotherconn minus 5 (unit is the network's corresponding unit)
    >>> myconn.set_attr(p=Ref(myotherconn, 1.2, -5))

    # specify value and reference at the same time
    >>> myconn.p_ref.set_attr(ref=Ref(myotherconn, 1.2, -5), is_set=True)
    >>> myconn.p.set_attr(val=7, is_set=True)

    # possibilities to unset values
    >>> myconn.set_attr(p=None)
    >>> myconn.p.is_set
    False

    >>> myconn.set_attr(p=10)
    >>> myconn.p.set_attr(is_set=False)
    >>> myconn.p.is_set
    False

    >>> myconn.p_ref.set_attr(is_set=False)  # for referenced values
    >>> myconn.p_ref.is_set
    False

Fluid specification
^^^^^^^^^^^^^^^^^^^

If you want to specify the fluid vector you can do it in the following way.

.. code-block:: python

    # set both elements of the fluid vector
    >>> myconn.set_attr(fluid={'water': 1})

    # same thing, but using data container
    >>> myconn.fluid.set_attr(_val={'water': 1}, _is_set={'water'})
    >>> myconn.fluid.is_set
    {'water'}

    # set starting values
    >>> myconn.set_attr(fluid0={'water': 1})

    # same thing, but using data container
    >>> myconn.fluid.set_attr(val0={'water': 1})

    # unset full fluid vector
    >>> myconn.set_attr(fluid={'water': None})
    >>> myconn.fluid.is_set
    set()

    # unset part of fluid vector
    >>> myconn.set_attr(fluid={'water': 1})
    >>> myconn.fluid.is_set.remove('water')
    >>> myconn.fluid.is_set
    set()

.. note::

    References can not be used for fluid composition at the moment!

It is possible to specify the fluid property back end of the fluids by adding
the name of the back end in front of the fluid's name. For incompressible binary
mixtures, you can append the water volume/mass fraction to the fluid's name, for
example:

.. code-block:: python

    >>> myconn.set_attr(fluid={'water': 1})  # HEOS back end
    >>> myconn.set_attr(fluid={'INCOMP::water': 1})  # incompressible fluid
    >>> myconn.set_attr(fluid={'BICUBIC::air': 1})  # bicubic back end
    >>> myconn.set_attr(fluid={'INCOMP::MPG[0.5]|mass': 1})  # binary incompressible mixture

.. note::

    Without further specifications CoolProp will be used as fluid property
    database. If you do not specify a back end, the **default back end**
    :code:`HEOS` will be used. For an overview of the back ends available please
    refer to the :ref:`fluid property section <tespy_fluid_properties_label>`.

You can also change the engine, for example to the iapws library. It is even
possible, that you define your own custom engine, e.g. using polynomial
equations. Please check out the fluid properties' section in the docs on how to
do this.

.. code-block:: python

    >>> from tespy.tools.fluid_properties.wrappers import IAPWSWrapper
    >>> myconn.set_attr(fluid={'H2O': 1}, fluid_engines={"H2O": IAPWSWrapper})

Access from the :code:`Network` object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

    >>> from tespy.networks import Network

    >>> mynetwork = Network()
    >>> myconn = Connection(source1, 'out1', sink1, 'in1', label='myconnlabel')
    >>> mynetwork.add_conns(myconn)
    >>> mynetwork.get_conn('myconnlabel').set_attr(p=1e5)
    >>> myconn.p.val
    100000.0

.. note::

    The label can only be specified on creation of the connection. Changing the
    label after might break this access method.

.. _tespy_powerconnections_label:

PowerConnections
================

PowerConnections can be used to represent non-material energy flow, like power
or heat. You can make use of generators, motors and buses.

Different use-cases for the implementation of powerconnections with the
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
- :py:class:`tespy.components.power.sink.PowerSink`: e.g. represet power fed into the electricity grid
- :py:class:`tespy.components.power.source.PowerSource`: e.g. represent power drawn from the electricity grid

For more details on the components please go to the respective section of the
:ref:`documentaton <tespy_modules_components_label>` and the respective API
documentation linked in the list above.

Parameters
^^^^^^^^^^

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
    >>> nw = Network(p_unit="bar", T_unit="C", iterinfo=False)
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
    >>> round(e1.e.val / 1e3)
    577

.. note::

    Note that the value of the energy flow of a :code:`PowerConnection` will
    always be positive in the defined direction (from one component's outlet
    to another component's inlet).

To learn what power connections are available in each of the component classes
see the respective API documentation. Below you will find more examples
utilizing the :code:`PowerConnection`.

Example: Single shaft gas turbine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
    >>> nw = Network(T_unit="C", p_unit="bar", iterinfo=False)
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
    >>> round(e4.e.val_SI / 1e3, 1)
    246.9

Alternatively, we could also fix the electricity output to a specific
target value and unset the air mass flow. This will calculate the required
air mass flow to generate the desired amount of electricity.

.. code-block:: python

    >>> e4.set_attr(e=3e5)
    >>> c1.set_attr(m=None)
    >>> nw.solve("design")
    >>> round(c1.m.val, 3)
    1.215

Example: Single shaft feed water pump powered by a turbine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a pump that is powered by a turbine. The turbine's :code:`turbine_fwp`
power output must therefore be equal to the pump's :code:`fwp` power
consumption.

.. code-block:: python

    >>> from tespy.networks import Network
    >>> from tespy.components import Pump, Turbine, Source, Sink
    >>> from tespy.connections import Connection, PowerConnection

    >>> nw = Network(p_unit="bar", T_unit="C", iterinfo=False)
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
    >>> round(e1.e.val_SI / 1e3)
    68

Example: Logic to force same power of two turbines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Example: Including part load model for motor efficiency
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As mentioned in the component section: It is also possible to import your
custom characteristics from the :code:`HOME/.tespy/data` folder. Read more
about this :ref:`here <tespy_modules_characteristics_label>`.
