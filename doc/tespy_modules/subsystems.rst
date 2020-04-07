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

    from tespy.networks import network
    from tespy.components import source, sink
    from tespy.connections import connection
    import numpy as np

    from mysubsystems import waste_heat_steam_generator as whsg

    # %% network definition

    fluid_list = ['air', 'water']
    nw = network(fluid_list, p_unit='bar', T_unit='C')

    # %% component definition

    feed_water = source('feed water inlet')
    steam = sink('live steam outlet')

    waste_heat = source('waste heat inlet')
    chimney = sink('waste heat chimney')

    sg = whsg('waste heat steam generator')

    # %% connection definition

    fw_sg = connection(feed_water, 'out1', sg.comps['eco'], 'in2')
    sg_ls = connection(sg.comps['sup'], 'out2', steam, 'in1')
    fg_sg = connection(waste_heat, 'out1', sg.comps['sup'], 'in1')
    sg_ch = connection(sg.comps['eco'], 'out1', chimney, 'in1')

    nw.add_conns(fw_sg, sg_ls, fg_sg, sg_ch)
    nw.add_subsys(sg)

    # %% connection parameters

    fw_sg.set_attr(fluid={'air': 0, 'water': 1}, T=25)
    fg_sg.set_attr(fluid={'air': 1, 'water': 0}, T=650, m=100)

    sg_ls.set_attr(p=130)
    sg_ch.set_attr(p=1)

    sg.conns['eva_dr'].set_attr(x=0.6)

    # %% component parameters

    sg.comps['eco'].set_attr(pr1=0.999,  pr2=0.97, ttd_u=25,
                             design=['pr1', 'pr2', 'ttd_u'],
                             offdesign=['zeta1', 'zeta2', 'kA'])

    sg.comps['eva'].set_attr(pr1=0.999, ttd_l=20, design=['pr1', 'ttd_l'],
                             offdesign=['zeta1', 'kA'])

    sg.comps['sup'].set_attr(pr1=0.999,  pr2=0.99, ttd_u=50,
                             design=['pr1', 'pr2', 'ttd_u'],
                             offdesign=['zeta1', 'zeta2', 'kA'])

    # %% solve

    # solve design case
    nw.solve('design')
    nw.print_results()
    nw.save('tmp')

    # offdesign test
    nw.solve('offdesign', design_path='tmp')


Add more flexibility
--------------------

If you want to add even more flexibility, you might need to manipulate the
:code:`__init__` method of your custom subsystem class. Usually, you do not
need to override this method. However, if you need additional parameters, e. g.
in order to alter the subsystem's topology or specify additional information,
take a look at the standard
:py:meth:`__init__ <tespy.components.subsystems.subsystem>` method and add your
code between the label declaration and the components and connection creation.

For example, if you want a variable number of inlets and outlets because you
have a variable number of components groups within your subsystem, you may
introduce an attribute which is set on initialisation and lets you create and
parameterize components and connections generically. This might be very
interesting for district heating systems, turbines with several sections of
equal topology, etc.. For a good start, you can have a look at the
:code:`sub_consumer.py` of the district heating network in the
`oemof_examples <https://github.com/oemof/oemof-examples/tree/master/oemof_examples/tespy/district_heating>`_
repository.
