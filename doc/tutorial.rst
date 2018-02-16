Tutorial – heat pump
====================

.. contents::
    :depth: 1
    :local:
    :backlinks: top
	
1 Task
------
This tutorial deals with the creation of a heat pump model. You can see the component plan in figure 1.

.. figure:: api/_images/tutorial_heat_pump.svg
    :align: center
	
    Figure 1: Topology of the heat pump.

You have two possibilities to create the model. Either you create the whole circuit in one term or you divide it in sections, if the system gets too complex. In this tutorial the model is built up in four sections.
In the style of “How can I use TESPy?” it is necessary to implement network, components and connections. 


2 Set up a plant
----------------

In order to simulate a plant you have to create a tespy.network. The network is the main container for the model and is needed for all sections.
You need to specify a list of the fluids you need for the calculation in your plant. For the heat pump you need water (H2O) and ammonia (NH3). Water is used for the cold side of the heat exchanger, for the consumer and for the hot side of the environmental temperature. In contrast the ammonia is used for the circuit of heat pump.
Further it is possible to choose a unit system for variables for example mass flow, temperature, pressure and enthalpy. If you don’t specify the unit system, the variables are set to SI-Units. 
Closing, it is necessary to set a range for the used variables, in this simulation these are temperature, pressure and enthalpy.

Try to set a plant in your model. If you have problems, you will look at the solution in the next passage.

.. code-block:: python

    from tespy import nwk
	# create a network object with water and ammonia as fluids
	# set the unitsystem for temperatures to °C, for pressure to bar, for massflow to kg / s and for enthalpy to kJ / kg

    my_plant = nwk.network(fluids=['water' , 'NH3'], T='C', p='bar', h='kJ / kg', m='kg / s', p_range=[0.1, 100], T_range=[1, 500], h_range=[10, 10000])

Because of the circuit it is unnecessary where you start to create the model. Let´s start with the compressor-system.


3 Compressor-system
-------------------

3.1 Set up the components
^^^^^^^^^^^^^^^^^^^^^^^^^

In figure 1 you can see that the compressor-system is built of two compressors and one heat-exchanger. Because of a starting and an ending point, you need to create a source and a sink. Additional you need also a sink and a source on the cold side of the heat exchanger. 
Parameters for components are generally optional. On the one hand experiences offer that compressor 1 has an isentropic efficiency to 80 % and a pressure ratio form inlet to outlet to 3. On the other hand, compressor 2 has also an isentropic efficiency to 80 %, but a pressure ratio form inlet to outlet to 4. At least you need to parametrize the heat-exchanger. On the hot side there is a pressure ratio of 0.99 and on the cold side of 0.98.
Try to set components in your model. If you have problems, you will look at the solution in the next passage.


.. code-block:: python

	from tespy import cmp
	# %% components

	# source & sink

	source	= cmp.source('source')
	sink		= cmp.sink('sink')
	source1	= cmp.source('source1')
	sink1		= cmp.sink('sink1')

	# compressor-system

	compressor1    = cmp.compressor('compressor1', eta_s=0.8, dp=3)
	compressor2    = cmp.compressor('compressor2', eta_s=0.8, dp=4)
	heat_exchanger = cmp.heat_exchanger('heat_exchanger’, dp1=0.99, dp2=0.98)


3.2 Establish connections
^^^^^^^^^^^^^^^^^^^^^^^^^

Connections are used to link two components (outlet of component 1 to inlet of component 2, source to target). If two components are connected to each other the fluid properties at the source will be equal to the properties at the target. It is possible to set the properties on each connection in a similar way as parameters are set for components.
It is necessary that every connection has a fluid and the right amount of variables. For example pressure, temperature or mass flow. In this plant the mass flow of the heat pump circuit is 0.2.
Try to set connections in your model. If you have problems, you will look at the solution in the next passage.


.. code-block:: python
	
	from tespy import con
	a = con.connection(source, 'out1', compressor1, 'in1' , T=10, m=0.2, fluid={'water':0,'NH3':1}) # source to compressor setting temperature, mass flow and ammonia as fluid vector
	b = con.connection(compressor1, 'out1', heat_exchanger, 'in1', p=15) # setting pressure 
	c = con.connection(heat_exchanger, 'out1', compressor2, 'in1') # connection without parameter specification
	d = con.connection(compressor2, 'out1', sink, 'in1', T=170) #  setting temperature
	e = con.connection(source1, 'out1', heat_exchanger, 'in2', p=1, T=20, m=5, fluid={'water':1,'NH3':0} # setting temperature, pressure, mass flow and water as fluid vector
	f = con.connection(heat_exchanger,'out2',sink1,'in1') # connection without parameter specification
	my_plant.add_conns(a, b, c, d, e, f)


3.3 Simulate your plant
^^^^^^^^^^^^^^^^^^^^^^^

Now you need to create the solver for your network. To simulate your plant, follow the steps of “How can I use TESPy?”. For the first simulation it is enough to create a design simulation.

.. code-block:: python
	
	network.solve('design')
	network.process_components('post')
	network.save('compressor-system')
	

4 Condenser/consumer
--------------------

4.1 Set up the components
^^^^^^^^^^^^^^^^^^^^^^^^^

In figure 1 you can see that the condenser is built of one condenser (heat-exchanger), one pump and one consumer (heat-exchanger-simple). 
Parameters for components are generally optional. Experiences offer that the condenser has on the hot and cold side a pressure ratio of 0.99. At least experiences show that the pump has an isentropic efficiency to 80 % and a pressure ratio form inlet to outlet to 1.01.

Try to set components in your model. If you have problems, you will look at the solution in the next passage.


.. code-block:: python

	from tespy import cmp
	# %% components

	# source & sink

	source 	= cmp.source('source')
	sink		= cmp.sink('sink')
	source1	= cmp.source('source1')
	sink1		= cmp.sink('sink1')

	# condenser/consumer

	heat_exchanger        = cmp.heat_exchanger('condenser', dp1=0.99, dp2=0.99)
	pump                  = cmp.pump('pump', eta_s=0.8, dp=1.01)	
	heat_exchanger_simple = cmp.heat_exchanger_simple('consumer')


4.2 Establish connections
^^^^^^^^^^^^^^^^^^^^^^^^^

Try to set connections in your model. If you have problems, you will look at the solution in the next passage.


.. code-block:: python
	
	from tespy import con
	g = con.connection(source, 'out1', heat_exchanger, 'in1' , p=60, T=170, m=0.2, fluid={'water':0,'NH3':1}) # source to condenser setting pressure, temperature, mass flow and ammonia as fluid vector
	h = con.connection(heat_exchanger, 'out1', sink, 'in1') # connection without parameter specification
	i = con.connection(source1, 'out1', pump, 'in1', T=60, m=1.2, fluid={'water':1,'NH3':0}) # setting temperature, mass flow and water as fluid vector
	j = con.connection(pump, 'out1', heat_exchanger, 'in2', p=10) #  setting pressure
	k = con.connection(heat_exchanger, 'out2', heat_exchanger_simple, 'in1', T=105) # setting temperature
	l = con.connection(heat_exchanger_simple,'out1',sink1,'in1', h=con.ref(i,1,0), T=con.ref(i,1,0)) # setting the same temperature and enthalpy as you find at connection i
	my_plant.add_conns(g, h, i, j, k, l)


4.3 Simulate your plant
^^^^^^^^^^^^^^^^^^^^^^^

See section 3.3.

5 Vessel
--------

5.1 Set up the components
^^^^^^^^^^^^^^^^^^^^^^^^^

In figure 1 you can see that you only need a vessel. Because of a starting and an ending point, you need to create a source and a sink. 

Parameters for components are generally optional. Experiences offer that the vessel has a pressure ratio of 0.085.

Try to set components in your model. If you have problems, you will look at the solution in the next passage.


.. code-block:: python

	from tespy import cmp
	# %% components

	# source & sink

	source 	= cmp.source('source')
	sink		= cmp.sink('sink')

	# vessel

	vessel = cmp.vessel('vessel', dp=0.086)


5.2 Establish connections
^^^^^^^^^^^^^^^^^^^^^^^^^

Try to set connections in your model. If you have problems, you will look at the solution in the next passage.


.. code-block:: python
	
	from tespy import con
	m = con.connection(source, 'out1', vessel, 'in1', m=0.2, fluid={'water':0,'NH3':1}) # source to vessel setting mass flow and ammonia as fluid vector
	n = con.connection(vessel, 'out1', sink, 'in1', T=5, p=5.157) # setting pressure and temperature
	my_plant.add_conns(m, n)


5.3 Simulate your plant
^^^^^^^^^^^^^^^^^^^^^^^

See section 3.3.

6 Evaporator-system
-------------------

6.1 Set up the components
^^^^^^^^^^^^^^^^^^^^^^^^^

In figure 1 you can see that the evaporator is the most complex section. Besides the evaporator (heat-exchanger) the system is built up of a drum, a pump and a superheater (heat-exchanger).
 
Parameters for components are generally optional. Experiences offer that the evaporator has on hot and cold side a pressure ratio of 0.99. The superheater has the same pressure ratio on hot and cold side. At least it is known that the heat pump has an isentropic efficiency to 80 %.

Try to set components in your model. If you have problems, you will look at the solution in the next passage.


.. code-block:: python

	from tespy import cmp
	# %% components

	# source & sink

	source 	= cmp.source('source')
	sink		= cmp.sink('sink')
	source1	= cmp.source('source1')
	sink1		= cmp.sink('sink1')

	# evaporator-system

	drum			= cmp.drum('drum')
	evaporator	= cmp.heat_exchanger ('evaporator', dp1=0.99, dp2=0.99)
	superheater	= cmp.heat_exchanger('superheater', dp1=0.99, dp2=0.99)
	pump			= cmp.pump('pump', eta_s=0.8)


6.2 Establish connections
^^^^^^^^^^^^^^^^^^^^^^^^^

Setting the right parametrization is very difficult for the evaporator-system. It would therefore recommend that you set a four-time higher mass flow from pump to evaporator than from source to pump. //Erklären warum// Know your components, for example in the drum the pressure is everywhere the same.
Try to set connections in your model. If you have problems, you will look at the solution in the next passage. 


.. code-block:: python
	
	from tespy import con
	o = con.connection(source, 'out1', drum, 'in1' , p=5.157, T=5,  m=0.2, fluid={'water':0,'NH3':1}) # source to drum setting pressure, temperature, mass flow and ammonia as fluid vector
	p = con.connection(drum, 'out1', pump, 'in1') # connection without parameter specification
	q = con.connection(pump, 'out1', evaporator, 'in2', m=con.ref(o,4,0)) # setting a four-time higher mass flow as you find at connection o
	r = con.connection(evaporator, 'out2', drum, 'in2') # connection without parameter specification
	s = con.connection(drum, 'out2', superheater, 'in2') # connection without parameter specification
	t = con.connection(superheater, 'out2', sink, 'in1') # connection without parameter specification
	u = con.connection(source1, 'out1', superheater, 'in1', p=1, T=12, m=20, fluid={'water':1,'NH3':0}) # setting pressure, temperature, mass flow and water as fluid vector
	v = con.connection(superheater, 'out1', evaporator, 'in1', T=11.967) # setting temperature
	w = con.connection(evaporator, 'out1' ,sink1, 'in1') # connection without parameter specification
	my_plant.add_conns(o, p, q, r, s, t, u, v, w)


6.3 Simulate your plant
^^^^^^^^^^^^^^^^^^^^^^^

See section 3.3.
	
7 Linking sections to heat pump
-------------------------------

Now if every separated system works, you need to link all section to one model. You are going to create your heat pump. But you need to consider following things:
-	Set all components in one list. Attention: It isn´t allowed to label components same. (e.g.: source, source1, source2, ...)
-	You need also a starting and ending point. (TESpy can´t simulate a fleeting circuit)
-	You only need one-time ammonia as fluid vector.
-	You have to set less parametrization in connection. Besser erklärenXX
It would therefore recommend that you start to link the sections, one by one. (e.g. compressor-system to condenser). As soon as your two-section model work, you can add the vessel and as soon as your three-section model work, you can add the evaporator-system. 

Try to set components and connections in your model. If you have prob-lems, you can look up the complete solution in the next passage.

from tespy import cmp, con, nwk

	#%% network

	my_plant = nwk.network(fluids=['water', 'NH3'], T='C', p='bar', h='kJ / kg', m='kg / s', p_range=[0.001, 100], T_range=[1, 500], h_range=[10, 10000])

	# %% components

	# source & sink

	source	= cmp.source('source')
	sink		= cmp.sink('sink')
	source1	= cmp.source('source1')
	sink1		= cmp.sink('sink1')
	source2	= cmp.source('source2')
	sink2		= cmp.sink('sink2')
	source3	= cmp.source('source3')
	sink3		= cmp.sink('sink3')

	# compressor-system

	compressor1		= cmp.compressor('compressor1', eta_s=0.8, dp=3)
	compressor2		= cmp.compressor('compressor2', eta_s=0.8, dp=4)
	heat_exchanger	= cmp.heat_exchanger('heat_exchanger', dp1=0.99, dp2=0.98)

	# condenser

	heat_exchanger1			= cmp.heat_exchanger('condenser', dp1=0.99, dp2=0.99)
	pump							= cmp.pump('pump', eta_s=0.8, dp=1.01)
	heat_exchanger_simple	= cmp.heat_exchanger_simple('consumer')

	# vessel

	vessel	= cmp.vessel('vessel')

	# evaporator

	drum			= cmp.drum('drum')
	evaporator	= cmp.heat_exchanger ('evaporator', dp1=0.99, dp2=0.99)
	superheater	= cmp.heat_exchanger('superheater', dp1=0.99, dp2=0.99)
	pump1			= cmp.pump('pump1', eta_s=0.8)


	#%% connections

	# compressor-system

	a = con.connection(source, 'out1', compressor1, 'in1', T=10, m=0.2, fluid={'water':0,'NH3':1})
	b = con.connection(compressor1, 'out1', heat_exchanger, 'in1', p=15)
	c = con.connection(heat_exchanger, 'out1', compressor2, 'in1')
	d = con.connection(source1, 'out1', heat_exchanger, 'in2', p=1, T=20, m=5, fluid={'water':1,'NH3':0})
	e = con.connection(heat_exchanger, 'out2', sink1, 'in1')
	f = con.connection(compressor2, 'out1', heat_exchanger1, 'in1', T=170)
	my_plant.add_conns(a, b, c, d, e, f)

	# condenser

	g = con.connection(heat_exchanger1, 'out1', vessel, 'in1')
	h = con.connection(source2, 'out1', pump, 'in1', T=60, flu-id={'water':1,'NH3':0})
	i = con.connection(pump, 'out1', heat_exchanger1, 'in2', p=10)
	j = con.connection(heat_exchanger1, 'out2', heat_exchanger_simple, 'in1', T=105)
	k = con.connection(heat_exchanger_simple, 'out1', sink2, 'in1', h=con.ref(so2_pu1,1,0), T=con.ref(g,1,0))
	my_plant.add_conns(g, h, i, j, k)

	# vessel

	l = con.connection(vessel, 'out1', drum, 'in1')
	my_plant.add_conns(l)

	# evaporator

	m = con.connection(drum, 'out1', pump1, 'in1')
	n = con.connection(pump1, 'out1', evaporator, 'in2',  m=con.ref(ve_dr,4,0))
	o = con.connection(evaporator, 'out2', drum, 'in2')
	p = con.connection(drum, 'out2', superheater, 'in2')
	q = con.connection(superheater, 'out2', sink, 'in1', T=con.ref(a,1,0)) # setting the same pressure and temperature as you find at the starting connection (simulating the  fleeting circuit)
	my_plant.add_conns(m, n, o, p, q)

	r = con.connection(source3, 'out1', superheater, 'in1', T=12, p=1, m=20, fluid={'water':1,'NH3':0})
	s = con.connection(superheater, 'out1', evaporator, 'in1')
	t = con.connection(evaporator, 'out1', sink3, 'in1', T=10)
	my_plant.add_conns(r, s, t)


	#%% Solver

	nw.solve('design')
	nw.process_components('post')
	nw.save('heat pump')


In figure 2 you can see the heat pump including pressure, temperature, enthalpy and mass flow.	