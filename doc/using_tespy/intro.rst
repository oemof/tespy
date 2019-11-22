Introduction
============

Set up a plant
--------------

In order to simulate a plant you will have to create a tespy.network first. The network is the main container for the model.

You need to specify a list of the fluids you need for the calculation in your plant. For more information on the fluid properties jump to the :ref:`bottom of this page <tespy_fluid_properties_label>`.

.. code-block:: python

	from tespy import nwk
	# create a network object with air and water as fluids
	fluid_list = ['air', 'water']
	my_plant = nwk.network(fluids=fluid_list)

On top of that, it is possible to specify a unit system and value ranges for the networks variables. If you do not specify these, TESPy will use SI-units.
The specification of the **value range** is used to **improve convergence stability**, in case you are dealing with **fluid mixtures**, e. g. using a combustion chamber.

.. code-block:: python

	from tespy import nwk

	# set the unitsystem for temperatures to °C, for pressure to bar and enthalpy to kJ / kg
	my_plant.set_attr(T_unit='C', p_unit='bar', h_unit='kJ / kg')
	my_plant.set_attr(T_range=[10, 700], p_range=[0.05, 150], h_range=[15, 4000])

Now you can start to create the components of the network.

Set up components
-----------------

Available components can be found :ref:`here <using_tespy_components_label>`. If you set up a component you have to specify a (within one network) unique label.
Moreover, it is possible to specify parameters for the component, for example power P for a turbine or upper terminal temperature difference ttd_u of a heat exchanger.
The full list of parameters for a specific component (e. g. a valve) is stated in the classes documentation.

.. note::
	Parameters for components are generally optional. Only the components label and in case you want to use a stoichiometric combustion chamber, its fuel and air composition, are mandatory parameters to provide.
	If an optional parameter is not specified by the user, it will be a result of the plants simulation. In this way, the set of equations a component returns is determined by which parameters you specify.
	You can find all equations in the :ref:`components documentation <using_tespy_components_label>` as well. The example below shows how to create a component with specific parameters, set or reset and how to unset a parameter:

.. code-block:: python

	from tespy import cmp
	import numpy as np

	feed_water_pump = cmp.pump(label='hp pump', P=1e3) # create pump labeled 'hp pump'
	feed_water_pump.set_attr(P=2e3, eta_s=0.9) # set the power to 2000 W, set isentropic efficiency to 90 %
	feed_water_pump.set_attr(P=np.nan) # unset power

After setting up the components the next step is to connect the components in your network.

Establish connections
---------------------

Connections are used to link two components (outlet of component 1 to inlet of component 2, source to target).
If two components are connected to each other the fluid properties at the source will be equal to the properties at the target.
It is possible to set the properties on each connection in a similar way as parameters are set for components. You may specify:

 * mass flow* (m),
 * volumetric flow (v),
 * pressure* (p),
 * enthalpy* (h),
 * temperature* (T),
 * vapour mass fraction for pure fluids (x),
 * temperature difference to boiling point for pure fluids (Td_bp),
 * fluids state for pure fluids (state='l' for liquid or state='g' for gaseous),
 * a fluid vector (fluid) and
 * a balance closer for the fluid vector (fluid_balance).

All parameters but the fluid vector, state and balance have to be numeric values. The fluid vector has to be specified as dictonary, see the example below.
The parameter :code:`fluid_balance` can only be :code:`True` or :code:`False`, the parameter :code:`state` can only be :code:`'l'` (liquid) or :code:`'g'` (gaseous).
For the properties marked with * it is possible to use references instead of numeric values.
This can be used for example if you want to have the pressure in two parts of your network related in a specific way but you do not know the values prior to the plant simulation.

.. code-block:: python

	from tespy import con

	ws_cond = con.connection(waste_steam_source, 'out1', condenser, 'in1', x=0.97) # waste steam source to condenser hot side inlet and setting vapour mass fraction
	cond_cp = con.connection(condenser, 'out1', condensate_pump, 'in1', fluid={'water': 1, 'air': 0}, Td_bp=-3) # setting a fluid vector: {'fluid i': mass fraction i}, subcooling to 3 K (15/9 K if temperature unit is Fahrenheit)
	cp_fwt = con.connection(condensate_pump, 'out1', feed_water_tank, 'in1', state='l') # enthalpy values will be manipulated in calculation process in a way, that the fluids state is liquid all the time
	fwt_fwp = con.connection(feed_water_tank, 'out1', feed_water_pump, 'in1') # connection without parameter specification
	fwp_eco = con.connection(feed_water_pump, 'out1', economiser, 'in2', v=10) #  setting volumetric flow
	eco_drum = con.connection(economiser, 'out2', drum, 'in1', T=320, p=con.ref(fwp_eco, 0.98, 0)) # setting temperature and pressure via reference object (pressure at this point is 0.98 times of pressure at connection fwp_eco)
	eva_eco = con.connection(evaporator, 'out1', economiser, 'in1', T=350, m=100) # setting temperature and mass flow
	eco_fgs = con.connection(economiser, 'out1', flue_gas_sink, 'in1', fluid_balance=True, fluid={'air': 1}, p=1) # setting fluid vector partially as well as the fluid balance parameter and pressure

	# this line is crutial, you have to add all connections to your network!
	my_plant.add_conns(ws_cond, cond_cp, cp_fwt, fwt_fwp, fwp_eco, eco_drum, eva_eco, eco_fgs)

.. figure:: api/_images/intro_connections.svg
    :align: center

    Figure 2: Topology after defining the above connections.

If you want to set, reset or unset a connection parameter the same logic as for the components is applied.

.. code-block:: python

	ws_cond.set_attr(x=0.95, p=0.05) # reset vapour mass fraction, set pressure
	eco_drum.set_attr(p=np.nan) # unset pressure

Start your calculation
----------------------

After building your network, the components and the connections, add the following line at the end of your script and off you go:

.. code-block:: python

	my_plant.solve(mode='design')

Please be aware, that the execution of the lines of code above will not create a solvable TESPy network. For good first examples jump to the :ref:`TESPy examples <tespy_examples_label>`.

In order to get a good overview of the TESPy functionalities, the following sections will walk you through the different TESPy modules in detail.
