.. _using_tespy_label:

###########
Using TESPy
###########
	
TESPy provides a simulation pacakge for component based thermal engineering containing the most important
basic components of such plants. In the introduction you will learn the basics of modelling component based
plants in TESPy.

We give an overview on the available components, introduce you to creating you own components and component
groups and give a short introduction on how TESPys solver works and how to handle different calculations modes.
Information on handling of fluid properties can be found at the end of this page.

On top of a ` step by step tutorial <http://tespy.readthedocs.io/en/latest/tutorial.html>`_ on how to
set up a heat pump in TESPy, we provide two basic examples in the `examples section
<http://tespy.readthedocs.io/en/latest/examples.html>`_.

.. figure:: api/_images/tutorial_heat_pump.svg
    :align: center
	
    Figure 1: Topology of a heat pump.

Introduction
============

Set up a plant
--------------

In order to simulate a plant you will have to create a tespy.network. The network is the main container for the model.

You need to specify a list of the fluids you need for the calculation in your plant. For more information on the fluid properties jump to the `bottom of this page <http://tespy.readthedocs.io/en/dev/using_tespy.html#fluid-properties-in-tespy>`_.

.. code-block:: python

    from tespy import nwk
	# create a network object with air and water as fluids
	fluid_list = ['air', 'water']
    my_plant = nwk.network(fluids=fluid_list)

On top of that, it is possible to specify a unit system and value ranges for the networks variables:

.. code-block:: python

    from tespy import nwk
	
	# set the unitsystem for temperatures to °C and for pressure to bar
	my_plant.set_attr(T_unit='C', p_unit='bar', h_unit='kJ / kg')
	my_plant.set_attr(T_range=[10, 700], p_unit=[0.05, 150], h_unit=[15, 4000])

Now you can start to create the components of the network.


Set up components
-----------------

Available components can be found `here <http://tespy.readthedocs.io/en/dev/using_tespy.html#available-components>`_. If you set up a component you have to specify a (within one network) unique label. Moreover, it is possible to specify parameters for the component, for example power P for a turbine or upper terminal temperature difference ttd_u of a heat exchanger. The full list of parameters for a specific component (e. g. a vessel) is stated in the classes documentation.

Parameters for components are generally optional. Only the components label and in case you want to use a combustion chamber, the combustion chambers fuel, are mandatory parameters to provide. If an optional parameter is not specified by the user, it will be a result of the plants simulation. In this way, the set of equations a component returns is determined by which parameters you specify. You can find all equations in the `components documentation <http://tespy.readthedocs.io/en/dev/using_tespy.html#available-components>`_ as well. The example below shows how to create a component with specific parameters, set or reset and how to unset a parameter:

.. _pump-parametrisation:
.. code-block:: python

	from tespy import cmp
	my_pump = cmp.pump(label='hp pump', P=1e3) # create pump labeled 'hp pump'
	my_pump.set_attr(P=2e3, eta_s=0.9) # set the power to 2000 W, set isentropic efficiency to 90 %
	my_pump.set_attr(P=math.nan) # unset power
	
After setting up the components the next step is to connect the components in your network.

Establish connections
---------------------

Connections are used to link two components (outlet of component 1 to inlet of component 2, source to target). If two components are connected to each other the fluid properties at the source will be equal to the properties at the target. It is possible to set the properties on each connection in a similar way as parameters are set for components. You may specify:

 * mass flow*,
 * pressure*,
 * enthalpy*,
 * temperature*,
 * vapour mass fraction for pure fluids,
 * a fluid vector and
 * a balance closer for the fluid vector.

All parameters but the fluid vector have to be numeric values. The fluid vector has to be specified as dictonary, see the example below. The parameter :code:`fluid_balance` ca only be :code:`True` or :code:`False`. For the properties marked with * it is possible to use references instead of numeric values. This can be used for example if you want to have the pressure in two parts of your network related in a specific way but you do not know the values prior to the plant simulation.

.. code-block:: python
	
	from tespy import con
	ws_cond = con.connection(waste_steam_source, 'out1', condenser, 'in1', x=0.97) # waste steam source to condenser hot side inlet and setting vapour mass fraction
	cond_fwp = con.connection(condenser, 'out1', feed_water_pump, 'in1', fluid={'water': 1, 'air': 0}) # setting a fluid vector: {'fluid i': mass fraction i}
	fwt_fwp = con.connection(feed_water_tank, 'out1', feed_water_pump, 'in1') # connection without parameter specification
	fwp_eco = con.connection(feed_water_pump, 'out1', economiser, 'in2', p=150) #  setting pressure
	eco_drum = con.connection(economiser, 'out2', drum, 'in1', T=320, p=con.ref(d, 0.98, 0)) # setting temperature and pressure via reference object
	eva_eco = con.connection(evaporator, 'out1', economiser, 'in1', T=350, m=100) # setting temperature and mass flow
	eco_fgs = con.connection(economiser, 'out1', flue_gas_sink, 'in1', fluid_balance=True, fluid={'air': 1}, p=1) # setting fluid vector partially as well as the fluid balance parameter and pressure

If you want to set, reset or unset a connection parameter the same logic as for the components is applied.

.. code-block:: python

	ws_cond.set_attr(x=0.95, p=0.05) # reset vapour mass fraction, set pressure
	fwp_eco.set_attr(p=math.nan) # unset pressure
	
Busses: power connections
-------------------------

Busses can be used to add up the power of different turbomachinery or to add up heat flux of different heat exchangers within your network. The handling is very similar to connections and components. You need to add components to your busses as a list containing the component object and a factor, the power of said component will be multiplied with. Do not forget to add the busses to you network.

This can be used for easy post processing, e. g. to calculate thermal efficiency or you can build up relations between components in your network. If you want to use the busses for postprocessing only, you do not specify the sum of the power or heat flux on your bus. If you set a value for P (equal parameter for heat flux or power), an additional equation will be added to your network. This could be useful, e. g. for establishing relations between different components, for instance when using a steam turbine powered feed water pump. In the code example the power of the turbine and the feed water pump is added up and set to zero, as the turbines and feed water pumps power have to be equal in absolute value but have different sign. The sign can be manipulated, e. g. in order to design two turbines with equal power output.

.. code-block:: python
	
	from tespy import nwk, con
	
	...
	
	fwp_bus = con.bus('feed water pump', P=0) # set a value for the total power on this bus.
	fwp_bus.add_comps([turbine_fwp, 1], [fwp, 1])
	
	turbine_bus = con.bus('turbines', P=0) # set a value for the total power on this bus
	turbine_bus.add_comps([turbine_hp, 1], [turbine_lp, -1])
	# the values for the busses power can be altered by using .set_attr()
	
	power = con.bus('power output') # bus for postprocessing, no power (or heat flux) specified
	power.add_comps([turbine_hp, 1], [turbine_lp, 1])
	
	my_network.add_busses(fwp_bus, turbine_bus, power)
	
Two labels for busses have a predefined function in the postprocessing analysis: 'P_res' and 'Q_diss'. If you specify these labels for your busses, 'P_res' will be interpreted as the total power of your process and 'Q_diss' as total amount of dissipated heat flux (from the process, not internally). Given these key figures, thermal efficiency or COP will be calculated and an entropy analysis for your systems components will be performed.*

*Planned feature, not implemented yet!

Subsystems/Component groups
---------------------------

Subsystems are an easy way to add frequently used component groups such as a drum with evaporator or a preheater with desuperheater to your system. You can use the predefined subsystems or `create a subsytem yourself <http://tespy.readthedocs.io/en/dev/using_tespy.html#tespy-subsystems-component-groups>`_. Every subsystem must have two interfaces, an inlet interface and an outlet interface. These interfaces have a variable number of connections, which can be connected with the rest of your network. The example below uses the predefined subsystem preheater with desuperheater (:code:`ph_desup()`). The subsystems interfaces are subsystem.inlet and subsystem.outlet, both with two connections. All connections (and components) of the subsystem have to be added to the network in order to start a simulation. This can easily be done by adding the whole subsystem object to your network.

.. code-block:: python

	from tespy import subsys, cmp
	ext = cmp.source(label='extraction steam')
	cond = cmp.sink(label='condensate')
	fwc = cmp.source(label='feed water cold')
	fww = cmp.sink(label='feed water warm')

	# a preheater with desuperheater part
	preheater = subsys.ph_desup(label='sub1')

	# connections into the subsystem are attached to subsystem.inlet, connections out of the subsystem to subsystem.outlet
	ext_pre = connection(source, 'out1', subsystem.inlet, 'in1')
	pre_cond = connection(subsystem.outlet, 'out1', sink, 'in1')
	fwc_pre = connection(source2, 'out1',subsystem.inlet,'in2')
	pre_fwc = connection(subsystem.outlet, 'out2', sink2, 'in1')
	
	# parametrisation
	preheater.set_attr(ttd=8, pr1_desup=1, pr2_desup=1, pr1_cond=1, pr2_cond=1)
	
	ext_pre.set_attr(m=5, p=4, h=29e5, fluid={'water': 1})
	fwc_pre.set_attr(p=50, h=3e5, fluid={'water': 1})
	pre_fwc.set_attr(p0=50)

	# create the network and connections and subsystems
	my_plant.add_conns(ext_pre, pre_cond, fwc_pre, pre_fwc)
	my_plant.add_subsys(subsys)

Start your calculation
----------------------

At the bottom of your script add the following line and off you go! Additional/advanced information on the solving process and which options are available are found `here <http://tespy.readthedocs.io/en/dev/using_tespy.html#solving-a-tespy-network>`_.

.. code-block:: python

	my_plant.solve(mode='design')
	
How can TESPy contribute to your energy system calculations?
------------------------------------------------------------

In this part you learn how you can use TESPy for your energy system calculations: In energy system calculations, for instance in oemof-solph, plants are usually modelled as abstract components on a much lower level of detail. In order to represent a plant within an abstract component it is possible to supply characteristics establishing a connection between your energy system model and a specific plant model. Thus the characteristics are a representation of a specific plant layout in terms of topology and process parameters.

The following part will show how to generate characteristics for a CHP unit. There are various technologies and concepts, for this example we will generate characteristics for a simple CHP with a backpressure steam turbine and a regenerative reheating unit as shown in the figure below. We want the characteristics to provide a correlation between output power and output heat flux at different temperatures of flow into a district heating system.

.. figure:: api/_images/CHP.svg
    :align: center
	
    Topology of the power plant.

Important design information can be obtained from the table below, the locations are indicated in the figure. After designing the plant, the mass flow in the main steam cycle has been changed stepwise from a slight overload of 50 kg/s to lower part loads (30 kg/s) with a stepwidth of 5 kg/s. Further the required temperature for the heating system was changed from 80 °C to 120 °C in steps of 10 K.

=========== =============== ======= ========
 location    parameter       value   unit
=========== =============== ======= ========
 fs          | pressure      | 100   | bar
             | temperature   | 550   | °C
             | mass flow     | 47    | kg/s
----------- --------------- ------- --------
 extr        pressure        10      bar
----------- --------------- ------- --------
 condenser   ttd_u :sup:`2`  10       K
----------- --------------- ------- --------
 reheater    ttd_u :sup:`2`  7       K
----------- --------------- ------- --------
 from_hs     | pressure      | 10    | bar
             | temperature   | 60    | °C
----------- --------------- ------- --------
 to_hs       temperature     110     °C
=========== =============== ======= ========

2: ttd_u is the upper terminal temperature difference, defined as temperature difference between hot side inlet and cold side outlet.

As a result we get the PQ-diagram of this power plant containing the characteristics at different temperatures in the heating system. Within your oemof-solph energy system it is now possible to implement the characteristic lines as a function of the temperature level in the heating system.

.. figure:: api/_images/PQ_diagram.svg
    :align: center
	
    PQ-diagram for a CHP unit.
	
Download the :download:`source file <../examples/chp.py>` of this example.
	
Solving a TESPy Network
=======================

Before learning how solve your TESPy network a short introduction on how the solution process works is provdided below.

Algorithm
---------

A TESPy Network can be represented as a linear system of non-linear equations, consequently the solution is obtained with numerical methods. TESPy uses the n-dimensional newton algorithm to find the systems solution, which may only be found, if the network is parameterized correctly. The variables of the system are mass flow, pressure, enthalpy and the fluid components on each connection of the network. Thus, the number of fluids you specify in the fluid list for the network and the number of connections determine the number of variables in the system:

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

This means that you have to provide the exact amount of required parameters (neither less nor more) and the parametrisation must not lead to linear dependencies. Each parameter you set for a connection or each power respectively heat flux you set for a bus will add one equation. On top, each component provides a different amount of basic equations plus the equations provided by your component specification. For example, setting the power of a pump results in an additional equation compared to a pump without specified power:

.. math::
	\forall i \in \mathrm{network.fluids} \, &0 = fluid_{i,in} - fluid_{i,out}\\
											 &0 = \dot{m}_{in} - \dot{m}_{out}\\
					 \mathrm{additional:} \, &0 = 1000 - \dot{m}_{in} (\cdot {h_{out} - h_{in}})
					 
Solving					 
-------

After you added all of your connections, subsystems and busses to your network, you can start the calculation with the following command.

.. code-block:: python

	nw.solve(init_file=None, design_file=None, mode='design', dec='.', max_iter=50)
	
This starts the initialisation of your network and proceeds to its calculation.

* :code:`nw` is the network object,
* :code:`init_file` is the .csv-file you want to use for initialisation,
* :code:`design_file` is the .csv-file which holds the information of your plants design point,
* :code:`mode` is the calculation mode (design-calculation or offdesign-calculation) and
* :code:`max_iter` is the maximum amount of iterations performed by the solver.

There are two modes available (:code:`'design'` and :code:`'offdesign'`). If you choose :code:`offdesign` as calculation mode a design file must be specified. The initialisation file is always optional but very valuable, if you specify it to be :code:`None`, the initialisation from .csv-file will be skipped.

Initialisation
^^^^^^^^^^^^^^

The newton algorithm requires starting values for all variables of the system. A high quality of initial values (low deveiation from solution) improves convergence speed and stability, whereas bad starting values might lead to instabilty and diverging calculation can be the result. In order to provide good initial values you can choose between three different initialisation options:

* initialise with standard values,
* provide starting values on your connections (see connection d in the subsystem example, usage: :code:`m0, p0, h0`) and
* provide a .csv-file of a previously calculated network.

The last option usually yields the best performance and is highly receommended. In order to initialise your calculation from a .csv-file, you need to provide the filename (e. g. myfile_results.csv). The file does not need to contain all connections of your network, thus you can build up your network bit by bit and initialise the existing parts of your network from the .csv-file. Be aware that a change within the fluid vector does not allow this practice. Thus, if you plan to use additional fluids in parts of the network you have not touched until now, you will need to state all fluids from the beginning.

Postprocessing
^^^^^^^^^^^^^^

The postprocessing has three functions you can apply to your calculation:

* plot the convergence history (:code:`nw.plot_convergence()`),
* print the results to prompt (:code:`nw.print_results()`) and
* save the results in a .csv-file (:code:`nw.save(filename, dec='.')`).

The main purpose of the plotting function is trouble shooting when your calculation does not converge. Therefore you specify a maximum number of iterations for the solver (:code:`max_iter`). As a result you get a plot of mass flow, pressure and enthalpy on all connections of your network. From there it might be possible to identify e. g. oscillating values or values that stay beyond the specified bounds of the fluid properties.

You can print the components and its properties to the prompt or, if you choose to save your results in a .csv-file, open the file and look up the components results in the file 'filename_comp.csv'. The mass flows and fluid properties of all connections are stored in the file 'filename_conn.csv'. On top, you can specify the decimal separator with :code:`nw.save(filename, dec='.')`.

Offdesign calculation
^^^^^^^^^^^^^^^^^^^^^
	
After designing your process you might want to gain information on offdesign behaviour. By stating :code:`'offdesing'` as calculation mode, you can auto-switch the components and connections to offdesign mode. This means, that all parameters given in :code:`component.design` will be unset and instead all parameters provided in :code:`component.offdesign` will be set. The same action is performed for the connections.

The default design and offdesign parameters for components can be found in the components documentation. For connections, there are no default design and offdesign parameters. For example, in order to specify custom design and offdesign parameters for a turbine use

.. code-block:: python

	turbine.set_attr(design=['P', 'eta_s'], offdesign=['cone', 'char'])
	
and for connections it works in the same way.

.. code-block:: python

	connection.set_attr(design=['h'], offdesign=['T'])
	
The table below contains frequently used offdesign parameters of the components.

=======================	======================	===================================================
 component             	 parameter            	 affects
=======================	======================	===================================================
 vessel                	 zeta                  	 pressure drop
-----------------------	----------------------	---------------------------------------------------
 pipe                  	 | zeta                	 | pressure drop
                       	 | k_s, D, L           	 | pressure drop (via dimensions and roughness)
                       	 | kA, t_a             	 | heat flux (using constant ambient temperature)
-----------------------	----------------------	---------------------------------------------------
 simple heat exchanger 	 see pipe              	  
-----------------------	----------------------	---------------------------------------------------
 heat exchanger        	 | zeta1              	 | pressure drop hot side
                       	 | zeta2              	 | pressure drop cold side
                       	 | kA                 	 | heat flux
-----------------------	----------------------	---------------------------------------------------
 pump                  	 char                  	 isentropic efficiency
-----------------------	----------------------	---------------------------------------------------
 turbine               	 | cone               	 | pressure drop, volumetric flow
                       	 | char                	 | isentropic efficiency
-----------------------	----------------------	---------------------------------------------------
 compressor            	 | char                	 | mass flow, pressure rise, isentropic efficiency
                       	 | vigv :sup:`1`         | see above, one arbitrary parameter less
=======================	======================	===================================================

1: When setting the vigv angle the characteristic map will be used for a specific vigv angle. The vigv angle is a result of the calculation, if you use the characteristic map only.

TESPy components
================

Available components
--------------------

More information on the components can be gathered from the code documentation. We have linked the base class containing a figure and basic informations as well as the equations.

- `Source <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.source>`_ (no equations)
- `Sink <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.sink>`_ (no equations)
- `Merge <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.merge>`_ (`equations <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.merge.equations>`_)
- `Splitter <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.splitter>`_ (`equations <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.splitter.equations>`_)
- `Vessel <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.vessel>`_ (`equations <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.vessel.equations>`_)
- Turbomachines
	* `Pump <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.pump>`_ (`equations <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.turbomachine.equations>`_)
	* `Compressor <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.compressor>`_ (`equations <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.turbomachine.equations>`_)
	* `Turbine <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.turbine>`_ (`equations <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.turbomachine.equations>`_)
- `Drum <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.combustion_chamber>`_ (`equations <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.combustion_chamber.equations>`_)
- Heat exchangers
	* `Heat exchanger <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.heat_exchanger>`_ (`equations <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.heat_exchanger.equations>`_)
	* `Equations <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.condenser>`_ (`equations <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.condenser.equations>`_)
	* `Desuperheater <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.desuperheater>`_ (`equations <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.desuperheater.equations>`_)
	* `Heat exchanger simple <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.heat_exchanger_simple>`_ (`equations <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.heat_exchanger_simple.equations>`_)
	* `Pipe <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.pipe>`_ (`equations <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.pipe.equations>`_)
- `Drum <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.drum>`_ (`equations <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.drum.equations>`_)

Custom components
-----------------

If required, you can add custom components. These components should inherit from tespy.components.components class or its children. In order to do that, create a python file in your working directory and import the tespy.components.components module. The most important functions are

- :code:`attr(self)`,
- :code:`inlets(self)`,
- :code:`outlets(self)`,
- :code:`equations(self, nw)`,
- :code:`derivatives(self, nw)` and
- :code:`calc_parameters(self, nw)`,

where :code:`nw` is a tespy.networks.network object.

The starting lines of your file would look like this:

.. code:: python
	
	from tespy import cmp
	
	
	class my_custom_component(cmp.component):
	
	
Attributes
^^^^^^^^^^

:code:`attr(self)` must return a list with strings in it. These are the attributes you can specify when you want to parametrize your component.

Inlets and outlets
^^^^^^^^^^^^^^^^^^

:code:`inlets(self)` and :code:`outlets(self)` respectively must return a list of strings. The list may look like this:

.. code:: python

	def inlets(self):
		return ['in1', 'in2']

	def outlets(self):
		return ['out1', 'out2']

The number of inlets and outlets might even be generic, e. g. if you have added an attribute :code:`'num_in'` in :code:`attr(self)`:

.. code:: python

    def inlets(self):
        if self.num_in_set:
            return ['in' + str(i + 1) for i in range(self.num_in)]
        else:
            self.set_attr(num_in=2)
            return self.inlets()

Equations
^^^^^^^^^

The equations contain the information on the changes to the fluid properties within the component. Each equations must formulated in a way, that the correct result will be zero, e. g.:

.. math::

	0 = \dot{m}_{in} - \dot{m}_{out}
	
The equations method requires a tespy.networks.network object as parameter. You can aquire a list of the ingoing and outgoing equations by the following command:

.. code:: python

    def inlets(self):
        if self.num_in_set:
            return ['in' + str(i + 1) for i in range(self.num_in)]
        else:
            self.set_attr(num_in=2)
            return self.inlets()

The equations are added to a list one after another, which will be returned at the end.

Derivatives
^^^^^^^^^^^
	
You need to calculate the partial derivatives of the equations to all variables of the network. This means, that you have to calculate the partial derivatives to mass flow, pressure, enthalpy and all fluids in the fluid vector on each incomming or outgoing connection of the component.

Add all derivatives to a list (in the same order as the equations) and return the list as numpy array (:code:`np.asarray(list)`). The derivatives can be calculated analytically or numerically by using the inbuilt function :code:`ddx_func(self, inlets, outlets, func, dx, pos)`.

- :code:`inlets` and :code:`outlets` are a list of the connections at the inlets and the outlets,
- :code:`func` is the function you want to calculate the derivatives for,
- :code:`dx` is the variable you want to calculate the derivative to and
- :code:`pos` indicates the connection you want to calculate the derivative for, e. g. :code:`pos=1` means, that counting your inlets and outlets from low index to high index (first inlets, then outlets), the connection to be used is the second connection in that list.

For a good start just look into the source code of the inbuilt components. If you have further questions feel free to contact us.

TESPy subsystems/component groups
=================================

You can use subsystems in order to represent groups of different components. These are highly customizable and thus a very powerful tool, if you require to use specific component groups frequently. You will learn how to create your own subsytems. Create a .py file in your working-directory with the class-definition of your custom subsystem. This usually includes the following methods:

- :code:`attr`: list of subsystem attributes,
- :code:`create_comps`: define the number of interfaces and create the necessary components,
- :code:`set_comps`: parametrize the components with the defined attributes from :code:`attr`,
- :code:`create_conns`: create the subsystems topology and
- :code:`set_conns`: parametrize them.

The following section shows, how the different functions of a subsystem can be defined. The code is taken from the subsystem drum with evaporator and natural flow.

Your file will start with the following lines:

.. code-block:: python

	from tespy import con, cmp, subsys
	
	
	class dr_eva_natural (subsys.subsystem):

Add the attr method:

.. code-block:: python
	
	def attr(self):
		# define available attributes for subsystem
		# num_i and num_o are excluded, as they are predefined in this subsystem
		return ([n for n in subsys.subsystem.attr(self) if
				 n != 'num_i' and n != 'num_o'] +
				['dp1_eva', 'PP', 'circ_num'])

Create the components
---------------------

The inlet and the outlet of the subsystem must be an attribute of the subsystem in order to reference to these when you are creating a network and want to connect the subsystem to the rest of the network.

.. code-block:: python

	def create_comps(self):
		# create the components

		self.num_i = 2
		self.num_o = 2
		self.inlet = cmp.subsys_interface(label=self.label + '_inlet',
										  num_inter=self.num_i)
		self.outlet = cmp.subsys_interface(label=self.label + '_outlet',
										   num_inter=self.num_o)
		self.drum = cmp.drum(label=self.label + '_drum')
		self.evaporator = cmp.heat_exchanger(label=self.label + '_evaporator',
											 mode='man')

As specific attributes refer to specific components in the subsystem, it is necessery, that the evaporator is stored as attribute of the subsystem as well. Else it would not be possible to set values for the parametrization.


Parametrize the components
--------------------------

.. code-block:: python

	def set_comps(self):
		# set component parameters

		self.evaporator.set_attr(ttd_l=self.PP)
		self.evaporator.set_attr(pr1=self.pr1_eva)

Create the connections
----------------------

In this example the components are saved in a list which is an attribute of the subsystem. As only the fourth and the last connections must be referenced in :code:`set_conns` it would be sufficient to store those connection as attributes of the subsystem.

.. code-block:: python

	def create_conns(self):
		# create the connections

		self.conns = []

		self.conns += [con.connection(self.inlet, 'out1', self.evaporator, 'in1')]
		self.conns += [con.connection(self.evaporator, 'out1', self.outlet, 'in1')]
		self.conns += [con.connection(self.inlet, 'out2', self.drum, 'in1')]
		self.conns += [con.connection(self.drum, 'out1', self.evaporator, 'in2')]
		self.conns += [con.connection(self.evaporator, 'out2', self.drum, 'in2')]
		self.conns += [con.connection(self.drum, 'out2', self.outlet, 'in2')]

Parametrize the connections
---------------------------

The connection gets a ref object as attribute, thus it is necessary to look, if the subsystems attribute is set or not. For parametrization with specific values simply use :code:`self.conns[3].set_attr(m=self.mass_flow)`. :code:`self.mass_flow` must be a subsystem attribute in this example.

.. code-block:: python
	def set_conns(self):
		# set connection parameters

		if self.circ_num_set:
			self.conns[3].set_attr(m=con.ref(self.conns[-1], self.circ_num, 0))
		else:
			self.conns[3].set_attr(m=np.nan)

Add more felxibility
--------------------

If you want to add even more flexibility, you might need to manipulate the :code:`__init__()` method. For example, if you want a variable number of inlets and outlets because you have a variable number of components groups within your subsystem, you may introduce an attribute which is set on initialisation and lets you create and parametrize components and connections generically. This might be very interesting for district heating systems, turbines with several sections of equal topology, etc..

Fluid properties in TESPy
=========================

The basic fluid properties are handled by `CoolProp <http://www.coolprop.org/>`_. All available fluids can be found on their homepage. 

Pure and pseudo-pure fluids
---------------------------

If you use pure fluids, TESPy directly uses CoolProp functions to gather all fluid properties. CoolProp covers the most important fluids such as water, air as a pseudo-pure fluid as well as its components, several fuels and refrigerants etc.. Look for the aliases in the `list of fluids http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids>`_. All fluids provided in this list cover liquid and gaseous state and the two-phase region.

Incompressible fluids
---------------------

If you are looking for heat transer fluids, the `list of incompressible fluids <http://www.coolprop.org/fluid_properties/Incompressibles.html>`_ might be interesting for you. In contrast to the pure fluids, the properties cover liquid state only.

Fluid mixtures
--------------

CoolProp provides fluid properties for two component mixtures. BUT: These are NOT integrated in TESPy! Nevertheless, you can use fluid mixtures for gases:

Ideal mixtures of gaseous fluids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TESPy can handle mixtures of gaseous fluids, by using the single fluid properties from CoolProp together with corresponding equations for mixtures. The equations can be found in the `tespy.helpers module <http://tespy.readthedocs.io/en/dev/api/tespy.html#module-tespy.helpers>`_ and are applied automatically to the fluid vector.

Other mixtures
^^^^^^^^^^^^^^

It is NOT POSSIBLE to use mixtures of liquid and other liquid or gaseous fluids AT THE MOMENT! If you try to use a mixture of two liquid or gaseous fluids and liquid fluids, e. g. water and methanol or liquid water and air, the equations will still be applied, but obviously return bad values. If you have ideas for the implementation of new kinds of mixtures we appreciate you contacting us.
