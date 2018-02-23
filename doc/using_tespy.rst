.. _using_tespy_label:

###########
Using TESPy
###########

.. figure:: api/_images/tutorial_heat_pump.svg
    :align: center
	
    Figure 1: Topology of a heat pump.
	
TESPy provides a simulation pacakge for component based thermal engineering containing the most important
basic components of such plants. The following parts show an overview on the available components,
component groups and on the handling of fluid properties. Furthermore you can learn more about the solution
process within TESPy and how to handle different calculations modes.

We provide a `tutorial <http://tespy.readthedocs.io/en/latest/tutorial.html>`_ on how to use TESPy.

TESPy components
================

Available components
--------------------

More information on the components can be gathered from the code documentation. We have linked the base class containing a figure and basic informations as well as the equations.

- `Source <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.source>`_ (`equations <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.source.equations>`_)
- `Sink <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.sink>`_ (`equations <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.sink.equations>`_)
- `Merge <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.merge>`_ (`equations <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.merge.equations>`_)
- `Splitter <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.splitter>`_ (`equations <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.splitter.equations>`_)
- `Vessel <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.vessel>`_ (`equations <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.vessel.equations>`_)
- Turbomachines
	* `Pump <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.pump>`_ (`equations <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.pump.equations>`_)
	* `Compressor <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.compressor>`_ (`equations <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.compressor.equations>`_)
	* `Turbine <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.turbine>`_ (`equations <http://tespy.readthedocs.io/en/dev/api/tespy.components.html#tespy.components.components.turbine.equations>`_)
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

* plot the convergence history (:code:`nw.plot_convergence(mode)`),
* print the results to prompt (:code:`nw.process_components(mode)`) and
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

TESPy can handle mixtures of gaseous fluids, by using the single fluid properties from CoolProp together with corresponding equations for mixtures. The equations can be found in the tespy.helpers module and are applied automatically to the fluid vector.

Other mixtures
^^^^^^^^^^^^^^

It is NOT POSSIBLE to use mixtures of liquid and other liquid or gaseous fluids AT THE MOMENT! If you try to use a mixture of two liquid or gaseous fluids and liquid fluids, e. g. water and methanol or liquid water and air, the equations will still be applied, but obviously return bad values. If you have ideas for the implementation of new kinds of mixtures we appreciate you contacting us.
