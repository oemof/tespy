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
 * vapour mass fraction for pure fluids (x),
 * a fluid vector (fluid) and
 * a balance closer for the fluid vector (fluid_balance).

It is possible to specify values, starting values, references and data
containers. The data containers for connections are dc_prop for fluid
properties (mass flow, pressure, enthalpy, temperature and vapour mass
fraction) and dc_flu for fluid composition. If you want to specify
data_containers, you need to import them from the :code:`tespy.tools` module.

In order to create the connections we create the components to connect first.

.. code-block:: python

    from tespy.tools import dc_prop
    from tespy.connections import connection, ref
    from tespy.components import sink, source

    # create components
    source1 = source('source 1')
    source2 = source('source 2')
    sink1 = sink('sink 1')
    sink2 = sink('sink 2')

    # creat connections
    myconn = connection(source1, 'out1', sink1, 'in1')
    myotherconn = connection(source2, 'out1', sink2, 'in1')

    # set pressure and vapour mass fraction by value, temperature and enthalpy
    # analogously
    myconn.set_attr(p=7, x=0.5)

    # set starting values for mass flow, pressure and enthalpy (has no effect
    # on temperature and vapour mass fraction!)
    myconn.set_attr(m0=10, p0=15, h0=100)

    # do the same with a data container
    myconn.set_attr(p=dc_prop(val=7, val_set=True),
                    x=dc_prop(val=0.5, val_set=True))
    myconn.set_attr(m=dc_prop(val0=10), p=dc_prop(val0=15),
                    h=dc_prop(val0=100))

    # specify a value in a different unit for a specific parameter
    myconn.set_attr(p=dc_prop(val=7, val_set=True, unit='MPa', unit_set=True)

    # specify a referenced value: pressure of myconn is 1.2 times pressure at
    # myotherconn minus 5 Pa (always SI unit here)
    myconn.set_attr(p=ref(myotherconn, 1.2, -5))

    # specify value and reference at the same time
    myconn.set_attr(p=dc_prop(val=7, val_set=True,
                    ref=ref(myotherconn, 1.2, -5), ref_set=True))

    # unset value and reference
    myconn.set_attr(p=np.nan)
    myconn.p.set_attr(val_set=False, ref_set=False)

If you want to specify the fluid vector you can do it in the following way:

.. code-block:: python

    from tespy.tools import dc_flu

    # set both elements of the fluid vector
    myconn.set_attr(fluid={'water': 1, 'air': 0})
    # same thing, but using data container
    myconn.set_attr(fluid=dc_flu(val={'water': 1, 'air': 0},
                    val_set:{'water': True, 'air': True}))

    # set starting values
    myconn.set_attr(fluid0={'water': 1, 'air': 0})
    # same thing, but using data container
    myconn.set_attr(fluid=dc_flu(val0={'water': 1, 'air': 0}))

    # unset values
    myconn.fluid.set_attr(val_set={'water': False, 'air': False})

References can not be used for fluid composition at the moment!


.. _tespy_busses_label:

Busses
------

Busses are energy flow connectors. You can sum the energy flow of different
components and create relations between components regarding mass free energy
transport.

Different use-cases for busses could be:

- Easy post-processing.
- Introduce motor or generator efficiencies.
- Create relations of different components.

The handling of busses is very similar to connections and components. You need
to add components to your busses as a dictionary containing at least the
instance of your component. Additionally you may provide a characteristic line,
linking the ratio of actual value to referenced value (design case value) to a
factor the actual value of the component is multiplied with on the bus. For
instance, you can provide a characteristic line of an electrical generator or
motor for a variable conversion efficiency. The referenced value is retrieved
by the design point of your system. Offdesign calculations use the referenced
value from your system's design point for the characteristic line. In design
case, the ratio will always be 1.

.. note::

    The available keywords for the dictionary are
    - 'c' for the component instance.
    - 'p' for the parameter (the combustion engine has various parameters,
      have a look at the
      :ref:`combustion engine example <combustion_engine_label>`).
    - 'P_ref' for the reference value of the component.
    - 'char' for the characteristic line.

    There are different specification possibilites:
    - If you specify the component only, the parameter will be default
      (not working with cogeneration unit) and the conversion factor of the
      characteristic line will be 1 for every load.
    - If you specify a numeric value for char, the conversion factor will be
      equal to that value for every load.
    - If you want to specify a characteristic line, provide a
      :py:class:` <tespy.components.characteristics.char_line>` object.

The examples below shows the implementation of busses in your TESPy simulation.

Create a pump that is powered by a turbine. The turbine's power output must
therefore be equal to the pump's power consumption.

.. code-block:: python

    from tespy.networks import network
    from tespy.components import pump, turbine, combustion_engine
    from tespy.connections import bus

    # the total power on this bus must be zero
    # this way we can make sure the power of the turbine has the same value as
    # the pump's power but with negative sign
    fwp_bus = bus('feed water pump bus', P=0)
    fwp_bus.add_comps({'c': turbine_fwp}, {'c': fwp})
    my_network.add_busses(fwp_bus)

Create two turbines which have the same power output.

.. code:: python

    # the total power on this bus must be zero, too
    # we make sure the two turbines yield the same power output by adding the char
    # parameter for the second turbine and using -1 as char
    turbine_bus = bus('turbines', P=0)
    turbine_bus.add_comps({'c': turbine_1}, {'c': turbine_2, 'char': -1})
    my_network.add_busses(turbine_bus)

Create a bus for post-processing purpose only. Include a characteristic line
of a generator.

.. code:: python

    # bus for postprocessing, no power (or heat flow) specified but with variable
    # conversion efficiency
    power_bus = bus('power output')
    x = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.1])
    y = np.array([0.85, 0.93, 0.95, 0.96, 0.97, 0.96])
    # createa characteristic line for a generator
    gen1 = char_line(x=x, y=y)
    gen2 = char_line(x=x, y=y)
    power.add_comps({'c': turbine_hp, 'char': gen1}, {'c': turbine_lp, 'char': gen2})
    my_network.add_busses(power_bus)

Create a bus for the electrical power output of a combustion engine. Use a
generator for power conversion an specify the total power output.

.. code:: python

    # bus for cogeneration unit power
    x = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.1])
    y = np.array([0.85, 0.93, 0.95, 0.96, 0.97, 0.96])
    # createa characteristic line for a generator
    gen = char_line(x=x, y=y)
    el_power_bus = bus('combustion engine power', P=10e6)
    el_power_bus.add_comps({'c': comb_engine, 'p': 'P', 'char': gen})


.. note::

    The x-values of the characteristic line represent the relative load of the
    component: actual value of the bus divided by the reference/design point
    value. In design-calculations the x-value used in the function evaluation
    will always be at 1.

As mentioned in the component section: If you want to learn more about the
usage of characteristic funcitons in TESPy have a look at
:ref:`this part <using_tespy_characteristics_label>`.
