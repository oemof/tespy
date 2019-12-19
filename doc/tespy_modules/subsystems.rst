.. _tespy_subsystems_label:

Subsystems and component groups
===============================

Usage
-----

REWORK USAGE PART 

Subsystems are an easy way to add frequently used component groups such as a
drum with evaporator or a preheater with desuperheater to your system. You can
:ref:`create a subsytem yourself <tespy_subsystems_label>`. Every subsystem
should have two interfaces, an inlet interface and an outlet interface. These
interfaces have a variable number of connections, which can be connected with
the rest of your network.

Custom subsystems
-----------------

You can use subsystems in order to represent groups of different components.
These are highly customizable and thus a very powerful tool, if you require to
use specific component groups frequently. You will learn how to create your own
subsystems. Create a .py file in your working-directory with the
class-definition of your custom subsystem. This includes the following
methods:

- :code:`__init__`: Initialisation of subsystem object. Do not override this
  method, if you do not need additional input parameters regarding the 
  subsystem's topology. However, if you need additional parameters, e. g. the
  number of components in a subsystem should be determined on creation, take 
  the standard :code:`__init__` method and add your code between the label 
  declaration and the components and connection creation.
- :code:`create_comps`: Create the components of your subsystem and save them
  in the :code:`subsystem.comps` attribute (dictionary).
- :code:`create_conns`: Create the connections of your subsystem and save them
  in the :code:`subsystem.conns` attribute (dictionary).

Example: Waste heat steam generator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following section shows, how the different functions of a subsystem can be
defined. We create a waste heat steam generator containing a superheater, an
evaporator and a drum as well as an economizer.

Start the file with the following lines (or add your class in at the end of the
:py:mod:`subsystems module<tespy.components.subsystems>`.

.. code-block:: python

    from tespy.components import subsystem, heat_exchanger, drum
    from tespy.connections import connection

    class waste_heat_steam_generator(subsystem):
		"""Class documentation"""
		
		def create_comps(self):
			"""Method documentation"""
			self.comps['eco'] = heat_exchanger('economizer')
			self.comps['eva'] = heat_exchanger('evaporator')
			self.comps['sup'] = heat_exchanger('superheater')
			self.comps['drum'] = drum('drum')
		
		def create_conns(self):
			"""Method documentation"""
			self.conns['eco_dr'] = connection(self.comps['eco'], 'out2',
			                                  self.comps['drum'], 'in1')
			self.conns['dr_eva'] = connection(self.comps['drum'], 'out1',
			                                  self.comps['eva'], 'in2')
			self.conns['eva_dr'] = connection(self.comps['eva'], 'out2',
			                                  self.comps['drum'], 'in2')
			self.conns['dr_sup'] = connection(self.comps['drum'], 'out2',
			                                  self.comps['sup'], 'in2')
			self.conns['sup_eva'] = connection(self.comps['sup'], 'out1',
			                                   self.comps['eva'], 'in1')
			self.conns['eva_eco'] = connection(self.comps['eva'], 'out1',
			                                   self.comps['eco'], 'in1')


Add more felxibility
^^^^^^^^^^^^^^^^^^^^

If you want to add even more flexibility, you might need to manipulate the
:code:`__init__()` method. For example, if you want a variable number of inlets
and outlets because you have a variable number of components groups within your
subsystem, you may introduce an attribute which is set on initialisation and
lets you create and parametrize components and connections generically. This
might be very interesting for district heating systems, turbines with several
sections of equal topology, etc.. For a good start, you can have a look into the
sub_consumer.py at the `tespy_examples repository <https://github.com/oemof/oemof-examples/tree/master/oemof_examples/tespy/district_heating>`_.
