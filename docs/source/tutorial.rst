.. _ppsim_label:

~~~~~~~~~~~
oemof-ppsim
~~~~~~~~~~~

Ppsim is an oemof-package for the simulation of thermal engeneering plants, such as thermal power plants or heat pumps. The package can be used to calculate stationary plant operation for design and offdesign cases. It is possible to derive plant characteristics which can be used in energy network simulations. In the :ref:`ppsim_examples_label` section you can find a short introduction how to use ppsim.

.. contents::
    :depth: 2
    :local:
    :backlinks: top


How can I use ppsim?
--------------------

To use ppsim you have to install oemof-ppsim. The installation provides some examples for testing purposes. Once the examples work you are close to your first power plant simulation.

Set up a plant
^^^^^^^^^^^^^^

In order to simulate a plant you will have to create a ppsim.network. The network is the main container for the model.

You need to specify a list of the fluids you need for the calculation in your plant. The whole list of available fluids can be found in the docs. If you want to work with fluid mixtures, for now, it is only possible to use ideal mixtures of real gases. Mixtures of liquids and gases or two liquids are not available. For more information on the fluid properties see the CoolProp documentation on `Pure and Pseudo-Pure fluid properties <http://www.coolprop.org/fluid_properties/PurePseudoPure.html>`_. Further it is possible to choose a unit system for mass flow, temperature and pressure.

.. code-block:: python

    import oemof.ppsim as pp
    my_plant = pp.network(fluids=['air', 'water'], T='C', p='bar')

Now you can start to create the components of the network.


Set up the components
^^^^^^^^^^^^^^^^^^^^^

All components can be found in the package components.components. If you set up a component you have to specify a (within one network) unique label. Moreover, it is possible to specify parameters for the component, for example power P for a turbine or upper terminal temperature difference dT_G_u of a heat exchanger. You get a full list of parameters for a specific component (e. g. a vessel) by typing:

.. code-block:: python

    import oemof.ppsim.components.components as comp
    print(comp.vessel.variables())
	
All parameters besides the label and the fuel for combustion chambers are optional, thus if not set, they will be a result of the plant simulation. In this way, the set of equations a component returns is determined by which parameters you specify. The example below shows how to create a component with specific parameters, set or reset and how to unset a parameter:

.. _pump-parametrisation:
.. code-block:: python

	my_pump = comp.pump(label='pump', P=1e3) # create pump
	my_pump.set_attr(P=2e3, eta_s=0.9) # reset power, set isentropic efficiency
	my_pump.set_attr(P=math.nan) # unset power
	
After setting up the components the next step is to connect the components in our network.

Establish connections
^^^^^^^^^^^^^^^^^^^^^

Connections are used to link two components (outlet of component 1 to inlet of component 2, source to target). If two components are connected to each other the fluid properties at the source will be equal to the properties at the target. It is possible to set the properties on each connection in a similar way as parameters are set for components. It is possible to specify:

 * mass flow*,
 * pressure*,
 * enthalpy*,
 * temperature*,
 * vapour mass fraction for pure fluids and
 * a fluid vector.

All parameters but the fluid vector have to be numeric values. The fluid vector has to be specified as dictonary, see the example below. For the properties marked with * it is possible to use references instead of numeric values. This can be used for example if you want to have the same pressure in two parts of your network but you do not know the pressure prior to the plant simulation.

.. code-block:: python
	
	import oemof.ppsim.connections.connections as conn
	import oemof.ppsim.connections.references as ref
	a = conn(waste_steam_source, 'out1', condenser, 'in1', x=0.97) # waste steam source to condenser hot side inlet and setting vapour mass fraction
	b = conn(condenser, 'out1', feed_water_pump, 'in1', fluid={'water': 1, 'air': 0}) # setting a fluid vector: {'fluid i': mass fraction i}
	c = conn(feed_water_tank, 'out1', feed_water_pump, 'in1') # connection without parameter specification
	d = conn(feed_water_pump, 'out1', economiser, 'in2', p=150) #  setting pressure
	e = conn(economiser, 'out2', drum, 'in1', T=320, p=ref.ref(d, 0.98, 0)) # setting temperature and pressure via reference object
	f = conn(evaporator, 'out1', economiser, 'in1', T=350, m=100) # setting temperature and mass flow
	g = conn(economiser, 'out1', flue_gas_sink, 'in1', fluid={'water': 0, 'air': 1}, p=1.013) # setting fluid vector and pressure

If you want to set, reset or unset a connection parameter the same logic as for the components is applied.

.. code-block:: python

	a.set_attr(x=0.95, p=0.05) # reset vapour mass fraction, set pressure
	d.set_attr(p=math.nan) # unset pressure
	
Busses: power connections
^^^^^^^^^^^^^^^^^^^^^^^^^

Busses can be used to add up the power of different turbomachinery or to add up heat flux of different heat exchangers within your network. This can be used either for easy post processing, e. g. to calculate thermal efficiency or you can build up relations between components in your network. If you want to use the busses for postprocessing only, you do not specify the sum of the power or heat flux on your bus. For establishing relations between different components, for instance when using a steam turbine powered feed water pump, you have to set the total power on this bus. In the code example the power of the turbine and the feed water pump is added up and set to zero, as the turbines and feed water pumps power have to be equal in absolute value but have different sign. The sign can be manipulated, e. g. in order to design two turbines with equal power output.

.. code-block:: python
	
	import oemof.ppsim.connections.busses as bus
	p = bus('feed water pump', P=0)
	p.add_comp([turbine_fwp, 1], [fwp, 1])
	p = bus('turbines', P=0)
	p.add_comp([turbine_hp, 1], [turbine_lp, -1])
	
Two labels for busses have a predefined function in the postprocessing analysis: 'P' and 'Q_diss'. If you specify these labels for your busses, 'P' will be interpreted as the total power of your process and 'Q_diss' as total amount of dissipated heat flux (from the process, not internally). Given these key figures, thermal efficiency and COP will be calculated and an entropy analysis for your systems components will be performed.

Subsystems
^^^^^^^^^^

Subsystems are an easy way to add frequently used component groups such as a drum with evaporator or a preheater with desuperheater to your system. You can use the predefined subsystems or create a subsytem yourself from a network object. Every subsystem must have two interfaces, an inlet interface and an outlet interface. These interfaces have a variable number of connections, which can be connected with the rest of your network. The example below uses the predefined subsystem preheater with desuperheater. The subsystems interfaces are subsys.inlet and subsys.outlet, both with two connections. All connections (and components) of the subsystem have to be added to the network in order to start a simulation. This can easily be done by adding the whole subsystem object to your network.

.. code-block:: python

	source = source(label='source1')
	sink = sink(label='sink1')
	source2 = source(label='source2')
	sink2 = sink(label='sink2')

	subsys = ph_desup(label='sub1', dT_G=8, dp1_desup=1, dp2_desup=1, dp1_cond=1, dp2_cond=1)

	a = connection(source, 'out1', subsys.inlet, 'in1', m=5, p=4, h=29e5, fluid={'water': 1})
	b = connection(subsys.outlet, 'out1', sink, 'in1')
	c = connection(source2, 'out1',subsys.inlet,'in2', p=50, h=3e5, fluid={'water': 1})
	d = connection(subsys.outlet, 'out2', sink2, 'in1', p0=50)

	nw = network(fluids=['water'], T='C')
	nw.add_conn(a, b, c, d)
	nw.add_subsys(subsys)

Simulate your plant
^^^^^^^^^^^^^^^^^^^

Before learning how to start the simulation a short introduction on how the solution process works is provdided below.

Introduction
++++++++++++

A ppsim.network can be represented as a linear system of non-linear equations, consequently the solution is obtained with numerical methods. ppsim uses the n-dimensional newton algorithm to find the systems solution, which may only be found, if the network is parameterized correctly. The variables of the system are mass flow, pressure, enthalpy and the fluid components on each connection of the network. Thus, the number of fluids you specify in the fluid list for the network and the number of connections determine the number of variables in the system:

.. math:: num_{var} = num_{conn} \cdot (3 + num_{fluids}).

The newton algorithm requires the calculation of residual values for the equations and partial derivatives of all variables (jacobian matrix). In the next step the matrix has to be inverted and multiplied with the residual vector to calculate the increment for the systems variables. This process is repeated until every equations result in the system is correct, thus the residual values are smaller than a specified error tolerance.

jacobian matrix J

.. math::
	J(\vec{x})=\left(\begin{array}{cccc}
	\frac{\partial f_1}{\partial x_1} & \frac{\partial f_1}{\partial x_2} & \cdots & \frac{\partial f_1}{\partial x_n} \\ 
	\frac{\partial f_2}{\partial x_1} & \frac{\partial f_2}{\partial x_2} & \cdots & \frac{\partial f_2}{\partial x_n} \\ 
	\vdots & \vdots & \ddots & \vdots \\
	\frac{\partial f_n}{\partial x_1} & \frac{\partial f_n}{\partial x_2} & \cdots & \frac{\partial f_n}{\partial x_n}
	\end{array}\right)
	
calculate increment

.. math::
	\vec{x}_{i+1}=\vec{x}_i-J(\vec{x}_i)^{-1}\cdot f(\vec{x}_i)
	
stop when

.. math::
	||f(\vec{x}_i)|| \leq \epsilon

This means that you have to provide the exact amount of required parameters (neither less nor more) and the parametrisation must not lead to linear dependencies. Each parameter you set for a connection or each power respectively heat flux you set for a bus will add one equation. On top, each component provides a different amount of basic equations plus the equations provided by your component specification. For example, setting the power of the pump above results in an additional equation compared to a pump without specified power:

.. math::
	\forall i \in \mathrm{network.fluids} \, &0 = fluid_{i,in} - fluid_{i,out}\\
											 &0 = \dot{m}_{in} - \dot{m}_{out}\\
					 \mathrm{additional:} \, &0 = 1000 - \dot{m}_{in} (\cdot {h_{out} - h_{in}})


Initialise the calculation
++++++++++++++++++++++++++

The newton algorithm requires starting values for all variables of the system. A high quality of initial values (low deveiation from solution) improves convergence speed and stability, whereas bad starting values might lead to instabilty and diverging calculation can be the result. In order to provide good initial values you can choose between three different initialisation options:

* initialise with standard values,
* provide starting values on your connections (see connection d in the subsystem example, usage: :code:`m0, p0, h0`) and
* provide a .csv-file of a previously calculated network.

The last option usually yields the best performance and is highly receommended. In order to initialise your calculation from a *.csv-file, you need to provide the filename *. The file does not need to contain all connections of your network, thus you can build up your network bit by bit and initialise the untouched part of your network from the .csv-file.

Solve the network
+++++++++++++++++

Starting with the subsystem example, in order to start your calculation you need to add the following line to your code:

.. code-block:: python

	solve.loop(nw, init_file=None, design_file=None, mode='design')
	
This starts the initialisation of your network and proceeds to its calculation.

* :code:`nw` is the network object,
* :code:`init_file` is the .csv-file you want to use for initialisation,
* :code:`design_file` is the .csv-file which holds the information of your plants design point and
* :code:`mode` is the calculation mode (design-calculation or offdesign-calculation).

There are two modes available (:code:`'design'` and :code:`'offdesign'`). If you choose :code:`offdesign` as calculation mode a design file must be specified. The initialisation file is always optional, if you specify it to be :code:`None`, the initialisation from .csv-file will be skipped.

Postprocessing
++++++++++++++

The preprocessing has three functions you can apply to your calculation:

* plot the convergence history,
* print the results to prompt and
* save the results in a .csv-file.

The plotting function is designed to use for trouble shooting when your calculation does not converge. Therefore you can specify a maximum number of iterations for the newton algorithm before calculation will be canceled. As a result you get a plot of mass flow, pressure and enthalpy on all connections of your network. From there it might be possible to identify e. g. oscillating values, which might be stabilised with improved initialisation parameters.

Offdesign calculations
++++++++++++++++++++++
	

After designing your process you might want to gain information on offdesign behaviour. By stating :code:`'offdesing'` as calculation mode, you can switch the component behaviour to offdesign. For example, this means that pressure drop in a pipe will be the result of reynolds number and the pipes dimensions. The table below shows all offdesign parameters available.

=======================	======================	===================================================
 component             	 parameter            	 affects
=======================	======================	===================================================
 vessel                	 zeta                  	 pressure drop
-----------------------	----------------------	---------------------------------------------------
 pipe                  	 | zeta :sup:`1`       	 | pressure drop
                       	 | dimensions :sup:`1` 	 | pressure drop
-----------------------	----------------------	---------------------------------------------------
 simple heat exchanger 	 zeta                 	 pressure drop
-----------------------	----------------------	---------------------------------------------------
 heat exchanger        	 | zeta1              	 | pressure drop hot side
                       	 | zeta2              	 | pressure drop cold side
                       	 | kA                 	 | heat flux
-----------------------	----------------------	---------------------------------------------------
 pump                  	 characteristic       	 isentropic efficiency
-----------------------	----------------------	---------------------------------------------------
 turbine               	 | cone law           	 | pressure drop, volumetric flow
                       	 | characteristic     	 | isentropic efficiency
-----------------------	----------------------	---------------------------------------------------
 compressor            	 | characteristic     	 | mass flow, pressure rise, isentropic efficiency
                       	 | vigv angle :sup:`2` 	 | see above, one arbitrary parameter less
=======================	======================	===================================================

1: If you set both parameters the length or the diameter must be a free parameter.

2: When setting the vigv angle the characteristic map will be used for a specific vigv angle. The vigv angle is a result of the calculation, if you use the characteristic map only

ppsim examples
--------------



(:download:`source file <../examples/solph/variable_chp/variable_chp.py>`, :download:`data file <../examples/solph/variable_chp/variable_chp.csv>`)

