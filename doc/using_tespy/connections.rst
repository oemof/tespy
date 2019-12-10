TESPy connections
=================

This section provides an overview of the parametrisation of connections, the usage of references and busses (connections for energy flow).

Parametrisation
---------------

As mentioned in the introduction, for each connection you can specify the following parameters:

 * mass flow* (m),
 * volumetric flow (v),
 * pressure* (p),
 * enthalpy* (h),
 * temperature* (T),
 * vapour mass fraction for pure fluids (x),
 * a fluid vector (fluid) and
 * a balance closer for the fluid vector (fluid_balance).

It is possible to specify values, starting values, references and data containers. The data containers for connections are dc_prop for fluid properties (mass flow, pressure, enthalpy, temperature and vapour mass fraction)
and dc_flu for fluid composition. If you want to specify data_containers, you need to import them from the :code:`tespy.tools` module.

.. code-block:: python

    from tespy.tools import dc_prop
    from tespy.connections import connection, ref

    # creat connections
    myconn = connection(example_comp1, 'out1', example_comp2, 'in1')
    myotherconn = connection(example_comp3, 'out', example_comp4, 'in1')
	
    # set pressure and vapour mass fraction by value, temperature and enthalpy analogously
    myconn.set_attr(p=7, x=0.5)

    # set starting values for mass flow, pressure and enthalpy (has no effect on temperature and vapour mass fraction!)
    myconn.set_attr(m0=10, p0=15, h0=100)

    # do the same with a data container
    myconn.set_attr(p=dc_prop(val=7, val_set=True), x=dc_prop(val=0.5, val_set=True))
    myconn.set_attr(m=dc_prop(val0=10), p=dc_prop(val0=15), h=dc_prop(val0=100))

    # specify a value in a different unit for a specific parameter
    myconn.set_attr(p=dc_prop(val=7, val_set=True, unit='MPa', unit_set=True)

    # specify a referenced value: pressure of myconn is 1.2 times pressure at myotherconn minus 5 Pa (always SI unit here)
    myconn.set_attr(p=ref(myotherconn, 1.2, -5))

    # specify value and reference at the same time
    myconn.set_attr(p=dc_prop(val=7, val_set=True, ref=ref(myotherconn, 1.2, -5), ref_set=True))

    # unset value and reference
    myconn.set_attr(p=np.nan)
    myconn.p.set_attr(val_set=False, ref_set=False)

If you want to specify the fluid vector you can do it in the following way:

.. code-block:: python

    from tespy.tools import dc_flu

    # set both elements of the fluid vector
    myconn.set_attr(fluid={'water': 1, 'air': 0})
    # same thing, but using data container
    myconn.set_attr(fluid=dc_flu(val={'water': 1, 'air': 0}, val_set:{'water': True, 'air': True}))

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

Busses can be used to add up the power of different turbomachinery or to add up heat flow of different heat exchangers within your network.
The handling is very similar to connections and components. You need to add components to your busses as a dictionary containing at least the instance of your component.
Additionally you may provide a characteristic line, linking the ratio of actual heat flow/power to referenced heat flow/power to a factor the actual heat flow/power of the component is multiplied with on the bus.
For instance, you can provide a characteristic line of an electrical generator or motor for a variable conversion efficiency. The referenced value (P_ref) is retrieved by the design point of your system.
Offdesign calculations use the referenced value from your system design point for the characteristic line. In design case, the heat flow/power ratio thus will be equal to 1.

.. note::
	The available keywords for the dictionary are

	- 'c' for the component instance,
	- 'p' for the parameter (the cogeneration unit has various parameters, have a look at the :ref:`cogeneration unit example <cogeneration_unit_label>`),
	- 'P_ref' for the reference heat flow/power value of the component and
	- 'char' for the characteristic line.

	There are different specification possibilites:

	- If you specify the component only, the parameter will be default (not working with cogeneration unit) and the conversion factor of the characteristic line will be 1 for every load.
	- If you specify a numeric value for char, the conversion factor will be that value for every load.
	- If you want to specify a characteristic line, you need to provide a :py:class:`TESPy characteristics <tespy.components.characteristics.characteristics>` object.

This can be used for easy post processing, e. g. to calculate thermal efficiency or you can build up relations between components in your network.
If you want to use the busses for postprocessing only, you must not specify the sum of the power or heat flow on your bus.
If you set a value for P (equal parameter for heat flow or power), an additional equation will be added to your network. This way the total heat flow/power of the bus will equal to the specified value.
This could be useful, e. g. for establishing relations between different components, for instance when using a steam turbine powered feed water pump.
In the code example the power of the turbine and the feed water pump is added up and set to zero, as the turbines and feed water pumps power have to be equal in absolute value but have different sign.
The sign can be manipulated, e. g. in order to design two turbines with equal power output.
Do not forget to add the busses to you network.

.. code-block:: python

    from tespy.networks import network
    from tespy.connections import bus
    from tespy.characteristics import characteristics

    ...

    fwp_bus = bus('feed water pump', P=0) # set a value for the total power on this bus.
    fwp_bus.add_comps({'c': turbine_fwp}, {'c': fwp})

    turbine_bus = bus('turbines', P=0) # set a value for the total power on this bus
    turbine_bus.add_comps({'c': turbine_hp}, {'c': turbine_hp, 'char': -1})
    # the values for the busses power can be altered by using .set_attr()

    power = con.bus('power output') # bus for postprocessing, no power (or heat flow) specified but with variable conversion efficiency
    x = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.1])
    y = np.array([0.85, 0.93, 0.95, 0.96, 0.97, 0.96])
    gen = characteristics(x=x, y=y) # characteristic line for a generator
    power.add_comps({'c': turbine_hp, 'char': gen}, {'c': turbine_lp, 'char': gen})

    chp = bus('chp power') # bus for cogeneration unit power
    chp.add_comps({'c': cog_unit, 'p': 'P', 'char': gen})

    my_network.add_busses(fwp_bus, turbine_bus, power)
	
.. note::

	The x-values of the characteristic line represent the relative load of the component: actual value of the bus divided by the reference/design point value.
	In design-calculations the x-value used in the function evaluation will always be at 1.
