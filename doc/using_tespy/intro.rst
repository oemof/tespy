Introduction
============

Set up a plant
--------------

In order to simulate a plant you will have to create a tespy.network first. The network is the main container for the model.
The model used in this introduction is shown in figure 2. It consists of a central heating plant and a consumer, represented by a heat exchanger with a valve.

.. figure:: api/_images/dhs_tut_scheme.svg

    :align: center

    Figure 2: Topology of the simplest district heating system

.. _using_tespy_introduction_label:


You need to specify a list of the fluids you need for the calculation in your plant. For more information on the fluid properties jump to the :ref:`bottom of this page <tespy_fluid_properties_label>`.

.. code-block:: python

	from tespy import nwk
	# create a network object with water as fluid
	fluid_list = ['water']
	my_plant = nwk.network(fluids=fluid_list)

On top of that, it is possible to specify a unit system and value ranges for the networks variables. If you do not specify these, TESPy will use SI-units.
The specification of the **value range** is used to **improve convergence stability**, in case you are dealing with **fluid mixtures**, e. g. using a combustion chamber.

.. code-block:: python

	from tespy import nwk

	# set the unitsystem for temperatures to °C, for pressure to bar and enthalpy to kJ / kg
	my_plant.set_attr(T_unit='C', p_unit='bar', h_unit='kJ / kg')
	my_plant.set_attr(T_range=[0, 100], p_range=[0.05, 150])

Now you can start to create the components of the network.

Set up components
-----------------

Available components can be found :ref:`here <using_tespy_components_label>`. If you set up a component you have to specify a (within one network) unique label.
Moreover, it is possible to specify parameters for the component, for example power P for a turbine or upper terminal temperature difference ttd_u of a heat exchanger.
The full list of parameters for a specific component (e. g. a valve) is stated in the classes documentation.

.. note::
	Parameters for components are generally optional. Only the components label and in case you want to use a stoichiometric combustion chamber, its fuel and air composition, are mandatory parameters to provide.
	If an optional parameter is not specified by the user, it will be a result of the plants simulation. In this way, the set of equations a component returns is determined by which parameters you specify.
	You can find all equations in the :ref:`components documentation <using_tespy_components_label>` as well. The example below shows how to create a component with specific parameters, set or reset and how to unset a parameter.
The definition of the other parameters can be found here: :py:class:`Pipe <tespy.components.piping.pipe>` and :py:class:`Heat exchanger simple <tespy.components.heat_exchangers.heat_exchanger_simple>`. 


.. code-block:: python
	
    from tespy import cmp
	import numpy as np
    
	# sources & sinks (central heating plant)
    
    so = cmp.source('heat source output')
    si = cmp.sink('heat source input')
    
    
    # consumer
    
    cons = cmp.heat_exchanger_simple(label='consumer')
    cons.set_attr(Q=-10000, pr=1)  # Q in W
    val = cmp.valve('valve')
    val.set_attr(pr=1)  # pr - pressure ratio (input/output) in per unit
    
    # pipes
    
    pipe_feed = cmp.pipe('pipe_feed')
    pipe_back = cmp.pipe('pipe_back')
    
    pipe_feed.set_attr(ks=0.0005,  # roughness in meters
                      L=100,  # length in m
                      D=0.06,  # diameter in m
                      kA=10,  # kA value - area and length independent! in W/K
                      Tamb=10)  # ambient temperature of the pipe environment (ground temperature)
    pipe_back.set_attr(ks=0.0005,
                      L=100,
                      D=0.06,
                      kA=10,
                      Tamb=10)


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
In this case, we just set input and output temperature of the system, as well as the input pressure.

.. code-block:: python

	from tespy import con
    import numpy as np
    
    # connections in the dhs

	so_pif = con.connection(so, 'out1', pipe_feed, 'in1')
    so_pif.set_attr(T=90, p=15, fluid={'water': 1})

	pif_cons = con.connection(pipe_feed, 'out1', cons, 'in1')
	cons_val = con.connection(cons, 'out1', val, 'in1', T=60, p=5)

	val_pib = con.connection(val, 'out1', pipe_back, 'in1')
	pib_si = con.connection(pipe_back, 'out1', si, 'in1')

    # this line is crutial: you have to add all connections to your network!
	my_plant.add_conns(so_pif, pif_cons, cons_val, val_pib, pib_si)


If you want to set, reset or unset a connection parameter the same logic as for the components is applied. In this example, we want to unset the output pressure, because it would overdetermine our system.

.. code-block:: python

    cons_val.set_attr(p=np.nan)  # unset pressure

Start your calculation
----------------------

After building your network, the components and the connections, add the following line at the end of your script and off you go:

.. code-block:: python

	my_plant.solve(mode='design')
    my_plant.print_results()

For further examples, that go deeper into TESPy, jump to the :ref:`TESPy examples <tespy_examples_label>`.

In order to get a good overview of the TESPy functionalities, the following sections will walk you through the different TESPy modules in detail.
