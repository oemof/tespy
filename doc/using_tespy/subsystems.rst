.. _tespy_subsystems_label:

TESPy subsystems/component groups
=================================

Usage
-----

Subsystems are an easy way to add frequently used component groups such as a drum with evaporator or a preheater with desuperheater to your system.
You can use the predefined subsystems or :ref:`create a subsytem yourself <tespy_subsystems_label>`. Every subsystem must have two interfaces, an inlet interface and an outlet interface.
These interfaces have a variable number of connections, which can be connected with the rest of your network. The example below uses the predefined subsystem preheater with desuperheater (:code:`ph_desup_cond()`).
The subsystems interfaces are subsystem.inlet and subsystem.outlet, both with two connections. All connections (and components) of the subsystem have to be added to the network in order to start a simulation.
This can easily be done by adding the whole subsystem object to your network.

.. code-block:: python

	from tespy import subsys, cmp
	ext = cmp.source(label='extraction steam')
	cond = cmp.sink(label='condensate')
	fwc = cmp.source(label='feed water cold')
	fww = cmp.sink(label='feed water warm')

	# a preheater with desuperheater part
	preheater = subsys.ph_desup(label='sub1')

	# connections into the subsystem are attached to subsystem.inlet, connections out of the subsystem to subsystem.outlet
	ext_pre = connection(ext, 'out1', subsystem.inlet, 'in1')
	pre_cond = connection(subsystem.outlet, 'out1', cond, 'in1')
	fwc_pre = connection(fwc, 'out1',subsystem.inlet,'in2')
	pre_fwc = connection(subsystem.outlet, 'out2', fww, 'in1')

	# parametrisation
	preheater.set_attr(ttd=8, pr1_desup=1, pr2_desup=1, pr1_cond=1, pr2_cond=1)

	ext_pre.set_attr(m=5, p=4, h=29e5, fluid={'water': 1})
	fwc_pre.set_attr(p=50, h=3e5, fluid={'water': 1})
	pre_fwc.set_attr(p0=50)

	# create the network and connections and subsystems
	my_plant.add_conns(ext_pre, pre_cond, fwc_pre, pre_fwc)
	my_plant.add_subsys(subsys)


.. figure:: api/_images/intro_subsys.svg
    :align: center

    Figure 3: Topology of the subsystem.

Custom subsystems
-----------------

You can use subsystems in order to represent groups of different components. These are highly customizable and thus a very powerful tool, if you require to use specific component groups frequently.
You will learn how to create your own subsystems. Create a .py file in your working-directory with the class-definition of your custom subsystem. This usually includes the following methods:

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
^^^^^^^^^^^^^^^^^^^^^

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
		self.evaporator = cmp.heat_exchanger(label=self.label + '_evaporator')

As specific attributes refer to specific components in the subsystem, it is necessery, that the evaporator is stored as attribute of the subsystem as well. Else it would not be possible to set values for the parametrization.

Parametrize the components
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

	def set_comps(self):
		# set component parameters

		self.evaporator.set_attr(ttd_l=self.PP)
		self.evaporator.set_attr(pr1=self.pr1_eva)

Create the connections
^^^^^^^^^^^^^^^^^^^^^^

Create a list called :code:`self.conns` and add the connections to that list.

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The connection gets a ref object as attribute, thus it is necessary to look, if the subsystems attribute is set or not.
For parametrization with specific values simply use :code:`self.conns[3].set_attr(m=self.mass_flow)`. :code:`self.mass_flow` must be a subsystem attribute in this example.

.. code-block:: python

	def set_conns(self):
		# set connection parameters

		if self.circ_num_set:
			self.conns[3].set_attr(m=con.ref(self.conns[-1], self.circ_num, 0))
		else:
			self.conns[3].set_attr(m=np.nan)

Add more felxibility
^^^^^^^^^^^^^^^^^^^^

If you want to add even more flexibility, you might need to manipulate the :code:`__init__()` method.
For example, if you want a variable number of inlets and outlets because you have a variable number of components groups within your subsystem,
you may introduce an attribute which is set on initialisation and lets you create and parametrize components and connections generically.
This might be very interesting for district heating systems, turbines with several sections of equal topology, etc..
For a good start, you can have a look into the sub_consumer.py at the `tespy_examples repository <https://github.com/oemof/oemof-examples/tree/master/oemof_examples/tespy/district_heating>`_.
