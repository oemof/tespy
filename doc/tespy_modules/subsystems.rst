.. _tespy_subsystems_label:

Subsystems and component groups
===============================

Subsystems are an easy way to add frequently used component groups such as a
drum with evaporator or a preheater with desuperheater to your system. In this
section you will learn how to create a subsystem and implement it in your work.
The subsystems are highly customizable and thus a very powerful tool, if you
require to use specific component groups frequently. We provide an example, of
how to create a simple subsystem and use it in a simulation.

Custom subsystems
-----------------

Create a :code:`.py` file in your working-directory. This file contains the
class definition of your subsystem and at minimum two methods:

- :code:`create_comps`: Method to create the components of your subsystem and
  save them in the :code:`subsystem.comps` attribute (dictionary).
- :code:`create_conns`: Method to create the connections of your subsystem and
  save them in the :code:`subsystem.conns` attribute (dictionary).

All other functionalities are inherited by the parent class of the
:py:class:`subsystem <tespy.components.subsystems.subsystem>` object.

Example
-------

Create the subsystem
^^^^^^^^^^^^^^^^^^^^

We create a subsystem for the usage of a waste heat steam generator. The
subsystem is built up of a superheater, an evaporator, a drum and an economizer
as seen in the figure below.

.. figure:: api/_images/subsystem_waste_heat_generator.svg
    :align: center

    Figure: Topology of the waste heat steam generator.

Create a file, e. g. :code:`mysubsystems.py` and add the following lines:
- Imports of the necessary classes from tespy.
- Class definition of the subsystem (inheriting from subsystem class).
- Methods for component and connection creation. Both, components and
  connections, are stored in a dictionary for easy access by their respective
  label.

.. code-block:: python

    from tespy.components import subsystem, heat_exchanger, drum
    from tespy.connections import connection

    class waste_heat_steam_generator(subsystem):
        """Class documentation"""

        def create_comps(self):
            """Create the subsystem's components."""
            self.comps['eco'] = heat_exchanger('economizer')
            self.comps['eva'] = heat_exchanger('evaporator')
            self.comps['sup'] = heat_exchanger('superheater')
            self.comps['drum'] = drum('drum')

        def create_conns(self):
            """Define the subsystem's connections."""
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

Import your subsystem
^^^^^^^^^^^^^^^^^^^^^

In a different script, we create a network and import the subsystem we just
created along with the different tespy classes required. The location of the
:code:`mysubsystems.py` file must be known by your python installation or lie
within the same folder as your script.

.. code-block:: python

    from tespy.components import source, sink
    from tespy.connections import connection
    from tespy.networks import network
    
    from mysubsystems import waste_heat_steam_generator as whsg


Add more felxibility
--------------------

- :code:`__init__`: Initialisation of subsystem object. Do not override this
  method, if you do not need additional input parameters regarding the
  subsystem's topology. However, if you need additional parameters, e. g. the
  number of components in a subsystem should be determined on creation, take
  the standard :code:`__init__` method and add your code between the label
  declaration and the components and connection creation.

If you want to add even more flexibility, you might need to manipulate the
:code:`__init__()` method. For example, if you want a variable number of inlets
and outlets because you have a variable number of components groups within your
subsystem, you may introduce an attribute which is set on initialisation and
lets you create and parametrize components and connections generically. This
might be very interesting for district heating systems, turbines with several
sections of equal topology, etc.. For a good start, you can have a look into the
sub_consumer.py at the `tespy_examples repository <https://github.com/oemof/oemof-examples/tree/master/oemof_examples/tespy/district_heating>`_.
