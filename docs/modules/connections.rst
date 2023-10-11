.. _tespy_modules_connections_label:

Connections
===========

This section provides an overview of the parametrisation of connections, the
usage of references and busses (connections for energy flow).

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


If you want to specify the fluid vector you can do it in the following way.

.. code-block:: python

    # set both elements of the fluid vector
    >>> myconn.set_attr(fluid={'water': 1})

    # same thing, but using data container
    >>> myconn.fluid.set_attr(val={'water': 1}, is_set={'water'})
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
    >>> myconn.set_attr(fluid={'INCOMP::MPG[0.5]': 1})  # binary incompressible mixture

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

.. _tespy_busses_label:

Busses
------

Busses are energy flow connectors. You can sum the energy flow of different
components and create relations between components regarding mass independent
energy transport.

Different use-cases for busses could be:

- post-processing
- introduce motor or generator efficiencies
- create relations of different components

The handling of busses is very similar to connections and components. You need
to add components to your busses as a dictionary containing at least the
instance of your component. Additionally, you may provide a characteristic line,
linking the ratio of actual value to a referenced value (design case value) to
an efficiency factor the component value of the bus is multiplied with. For
instance, you can provide a characteristic line of an electrical generator or
motor for a variable conversion efficiency. The referenced value is retrieved
by the design point of your system. Offdesign calculations use the referenced
value from your system's design point for the characteristic line. In design
case, the ratio will always be 1.

After a simulation, it is possible to output the efficiency of a component on
a bus and to output the bus value of the component using

- :code:`mycomponent.calc_bus_efficiency(mybus)`
- :code:`mycomponent.calc_bus_value(mybus)`

These data are also available in the network's results dictionary and contain

- the bus value,
- the component value,
- the efficiency value and
- the design value of the bus.

.. code-block:: python

    bus_results = mynetwork.results['power output']

.. note::

    The available keywords for the dictionary are:

    - 'comp' for the component instance.
    - 'param' for the parameter (e.g. the combustion engine has various
      parameters)
    - 'char' for the characteristic line
    - 'base' the base for efficiency definition
    - 'P_ref' for the reference value of the component

    There are different specification possibilities:

    - If you specify the component only, the parameter will be default and the
      efficiency factor of the characteristic line will be 1 independent of
      the load.
    - If you specify a numeric value for char, the efficiency factor will be
      equal to that value independent of the load.
    - If you want to specify a characteristic line, provide
      a :py:class:`CharLine <tespy.tools.characteristics.CharLine>`
      object.
    - Specify :code:`'base': 'bus'` if you want to change from the default base
      to the bus as base. This means, that the definition of the efficiency
      factor will change according to your specification.

      .. math ::

          \eta = \begin{cases}
          \frac{\dot{E}_\mathrm{component}}{\dot{E}_\mathrm{bus}} &
          \text{'base': 'bus'}\\
          \frac{\dot{E}_\mathrm{bus}}{\dot{E}_\mathrm{component}} &
          \text{'base': 'component'}
          \end{cases}

      This applies to the calculation of the bus value analogously.

      .. math::

          \dot{E}_\mathrm{bus} = \begin{cases}
          \frac{\dot{E}_\mathrm{component}}{f\left(
          \frac{\dot{E}_\mathrm{bus}}{\dot{E}_\mathrm{bus,design}}\right)} &
          \text{'base': 'bus'}\\
          \dot{E}_\mathrm{component} \cdot f\left(
          \frac{\dot{E}_\mathrm{component}}
          {\dot{E}_\mathrm{component,design}}\right) &
          \text{'base': 'component'}
          \end{cases}

The examples below show the implementation of busses in your TESPy simulation.

Create a pump that is powered by a turbine. The turbine's :code:`turbine_fwp`
power output must therefore be equal to the pump's :code:`fwp` power
consumption.

.. code-block:: python

    >>> from tespy.networks import Network
    >>> from tespy.components import Pump, Turbine, CombustionEngine
    >>> from tespy.connections import Bus

    >>> my_network = Network()
    >>> fwp = Pump("feed water pump")
    >>> turbine_fwp = Turbine("turbine fwp")

    # the total power on this bus must be zero
    # this way we can make sure the power of the turbine has the same value as
    # the pump's power but with negative sign
    >>> fwp_bus = Bus("feed water pump bus", P=0)
    >>> fwp_bus.add_comps({"comp": turbine_fwp}, {"comp": fwp, "base": "bus"})
    >>> my_network.add_busses(fwp_bus)

Create two turbines :code:`turbine1` and :code:`turbine2` which have the same
power output.

.. code-block:: python

    # the total power on this bus must be zero, too
    # we make sure the two turbines yield the same power output by adding the char
    # parameter for the second turbine and using -1 as char
    >>> turbine_1 = Turbine("turbine 1")
    >>> turbine_2 = Turbine("turbine 2")

    >>> turbine_bus = Bus('turbines', P=0)
    >>> turbine_bus.add_comps({'comp': turbine_1}, {'comp': turbine_2, 'char': -1})
    >>> my_network.add_busses(turbine_bus)

Create a bus for post-processing purpose only. Include a characteristic line
for a generator and add two turbines :code:`turbine_1` and :code:`turbine_2`
to the bus.

.. code-block:: python

    >>> import numpy as np
    >>> from tespy.tools.characteristics import CharLine

    # bus for postprocessing, no power (or heat flow) specified but with variable
    # conversion efficiency
    >>> power_bus = Bus("power output")
    >>> x = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.1])
    >>> y = np.array([0.85, 0.93, 0.95, 0.96, 0.97, 0.96])

    # create a characteristic line for a generator
    >>> gen1 = CharLine(x=x, y=y)
    >>> gen2 = CharLine(x=x, y=y)
    >>> power_bus.add_comps(
    ...     {'comp': turbine_1, 'char': gen1},
    ...     {'comp': turbine_2, 'char': gen2}
    ... )
    >>> my_network.add_busses(power_bus)

Create a bus for the electrical power output of a combustion engine
:code:`comb_engine`. Use a generator for power conversion and specify the total
power output.

.. code-block:: python

    >>> comb_engine = CombustionEngine("engine")

    # bus for combustion engine power
    >>> el_power_bus = Bus('combustion engine power', P=-10e6)
    >>> el_power_bus.add_comps({'comp': comb_engine, 'param': 'P', 'char': gen1})

Create a bus for the electrical power input of a pump :code:`pu` with
:code:`'bus'` and with :code:`'component'` as base. In both cases, the value of
the component power will be identical. Due to the different efficiency
definitions the value of the bus power will differ in part load.

.. code-block:: python

    >>> import numpy as np
    >>> from tespy.components import Pump, Sink, Source
    >>> from tespy.connections import Bus, Connection
    >>> from tespy.networks import Network
    >>> from tespy.tools.characteristics import CharLine

    >>> nw = Network(iterinfo=False, p_unit='bar', T_unit='C')

    >>> si = Sink('sink')
    >>> so = Source('source')
    >>> pu = Pump('pump')

    >>> so_pu = Connection(so, 'out1', pu, 'in1')
    >>> pu_si = Connection(pu, 'out1', si, 'in1')

    >>> nw.add_conns(so_pu, pu_si)

    # bus for combustion engine power
    >>> x = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.1])
    >>> y = np.array([0.85, 0.93, 0.95, 0.96, 0.97, 0.96])

    # create a characteristic line for a generator
    >>> mot_bus_based = CharLine(x=x, y=y)
    >>> mot_comp_based = CharLine(x=x, y=1 / y)
    >>> bus1 = Bus('pump power bus based')
    >>> bus1.add_comps({'comp': pu, 'char': mot_bus_based, 'base': 'bus'})

    # the keyword 'base': 'component' is the default value, therefore it does
    # not need to be passed
    >>> bus2 = Bus('pump power component based')
    >>> bus2.add_comps({'comp': pu, 'char': mot_comp_based})

    >>> nw.add_busses(bus1, bus2)

    >>> so_pu.set_attr(fluid={'H2O': 1}, m=10, p=5, T=20)
    >>> pu_si.set_attr(p=10)

    >>> pu.set_attr(eta_s=0.75)

    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> print('Bus based efficiency:', round(pu.calc_bus_efficiency(bus1), 2))
    Bus based efficiency: 0.97

    >>> print('Component based efficiency:', round(1 / pu.calc_bus_efficiency(bus2), 2))
    Component based efficiency: 0.97

    >>> print('Bus based bus power:', round(pu.calc_bus_value(bus1)))
    Bus based bus power: 6883

    >>> print('Component based bus power:', round(pu.calc_bus_value(bus2)))
    Component based bus power: 6883

    >>> so_pu.set_attr(m=8)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> print('Bus based efficiency:', round(pu.calc_bus_efficiency(bus1), 2))
    Bus based efficiency: 0.96

    >>> print('Component based efficiency:', round(1 / pu.calc_bus_efficiency(bus2), 2))
    Component based efficiency: 0.96

    >>> print('Bus based bus power:', round(pu.calc_bus_value(bus1)))
    Bus based bus power: 5562

    >>> print('Component based bus power:', round(pu.calc_bus_value(bus2)))
    Component based bus power: 5564

    # get DataFrame with the bus results
    >>> bus_results = nw.results['pump power bus based']

.. note::

    The x-values of the characteristic line represent the relative load of the
    component: actual value of the bus divided by the reference/design point
    value. In design-calculations the x-value used in the function evaluation
    will always be at 1.

As mentioned in the component section: It is also possible to import your
custom characteristics from the :code:`HOME/.tespy/data` folder. Read more
about this :ref:`here <tespy_modules_characteristics_label>`.
