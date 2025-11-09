.. _tespy_modules_subsystems_label:

Subsystems: Reusable Component Groups
=====================================

Subsystems are an easy way to add frequently used component groups such as a
drum with evaporator or a preheater with desuperheater to your system. In this
section you will learn how to create a subsystem and implement it in your work.
The subsystems are highly customizable and thus a very powerful tool, if you
require using specific component groups frequently. We provide an example, of
how to create a simple subsystem and use it in a simulation.

Custom subsystems
-----------------

Create a :code:`.py` file in your working-directory. This file contains the
class definition of your subsystem and at minimum one method:

- :code:`create_network`: Method to create the network of your subsystem.

On top of that you need to add methods to define the available interfaces of
your subsystem to the remaining network through specifying the number of inlets
and outlets in the :code:`__init__` method of your class as seen in the code
example below.

All other functionalities are inherited by the parent class of the
:py:class:`subsystem <tespy.components.subsystem.Subsystem>` object.

Example
-------

Create the subsystem
^^^^^^^^^^^^^^^^^^^^

We create a subsystem for the usage of a waste heat steam generator. The
subsystem is built up of a superheater, an evaporator, a drum and an economizer
as seen in the figure below.

.. figure:: /_static/images/modules/subsystem_waste_heat_generator.svg
    :align: center
    :alt: Topology of the waste heat steam generator
    :figclass: only-light

    Figure: Topology of the waste heat steam generator

.. figure:: /_static/images/modules/subsystem_waste_heat_generator_darkmode.svg
    :align: center
    :alt: Topology of the waste heat steam generator
    :figclass: only-dark

    Figure: Topology of the waste heat steam generator

Create a file, e.g. :code:`mysubsystems.py` and add the following lines:

- Imports of the necessary classes from tespy.
- Class definition of the subsystem (inheriting from subsystem class).
- Methods for component and connection creation. Both, components and
  connections, are stored in a dictionary for easy access by their respective
  label.

.. code-block:: python

    >>> from tespy.components import Subsystem, HeatExchanger, Drum
    >>> from tespy.connections import Connection

    >>> class WasteHeatSteamGenerator(Subsystem):
    ...     """Class documentation"""
    ...     def __init__(self, label):
    ...         self.num_in = 2
    ...         self.num_out = 2
    ...         super().__init__(label)
    ...
    ...     def create_network(self):
    ...         """Define the subsystem's connections."""
    ...         eco = HeatExchanger('economizer')
    ...         eva = HeatExchanger('evaporator')
    ...         sup = HeatExchanger('superheater')
    ...         drum = Drum('drum')
    ...
    ...         inlet_eco = Connection(self.inlet, 'out2', eco, 'in2', label='1')
    ...         eco_dr = Connection(eco, 'out2', drum, 'in1', label='2')
    ...         dr_eva = Connection(drum, 'out1', eva, 'in2', label='3')
    ...         eva_dr = Connection(eva, 'out2', drum, 'in2', label='4')
    ...         dr_sup = Connection(drum, 'out2', sup, 'in2', label='5')
    ...         sup_outlet = Connection(sup, 'out2', self.outlet, 'in2', label='6')
    ...
    ...         self.add_conns(inlet_eco, eco_dr, dr_eva, eva_dr, dr_sup, sup_outlet)
    ...
    ...         inlet_sup = Connection(self.inlet, 'out1', sup, 'in1', label='11')
    ...         sup_eva = Connection(sup, 'out1', eva, 'in1', label='12')
    ...         eva_eco = Connection(eva, 'out1', eco, 'in1', label='13')
    ...         eco_outlet = Connection(eco, 'out1', self.outlet, 'in1', label='14')
    ...
    ...         self.add_conns(inlet_sup, sup_eva, eva_eco, eco_outlet)

Make use of your subsystem
^^^^^^^^^^^^^^^^^^^^^^^^^^

We create a network and use the subsystem we just created along with the
different tespy classes required.

.. code-block:: python

    >>> from tespy.networks import Network
    >>> from tespy.components import Source, Sink
    >>> from tespy.connections import Connection
    >>> import numpy as np

    >>> # %% network definition
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(temperature="degC", pressure="bar")

    >>> # %% component definition
    >>> feed_water = Source('feed water inlet')
    >>> steam = Sink('live steam outlet')
    >>> waste_heat = Source('waste heat inlet')
    >>> chimney = Sink('waste heat chimney')

    >>> sg = WasteHeatSteamGenerator('waste heat steam generator')

    >>> # %% connection definition
    >>> fw_sg = Connection(feed_water, 'out1', sg, 'in2')
    >>> sg_ls = Connection(sg, 'out2', steam, 'in1')
    >>> fg_sg = Connection(waste_heat, 'out1', sg, 'in1')
    >>> sg_ch = Connection(sg, 'out1', chimney, 'in1')

    >>> nw.add_conns(fw_sg, sg_ls, fg_sg, sg_ch)
    >>> nw.add_subsystems(sg)

    >>> # %% connection parameters
    >>> fw_sg.set_attr(fluid={'water': 1}, T=25, m0=15)
    >>> fg_sg.set_attr(fluid={'air': 1}, T=650, m=100)
    >>> sg_ls.set_attr(p=130, T=600, design=['T'])
    >>> sg_ch.set_attr(p=1)

    >>> sg.get_conn('4').set_attr(x=0.6)

    >>> # %% component parameters
    >>> sg.get_comp('economizer').set_attr(
    ...     pr1=0.999,  pr2=0.97, design=['pr1', 'pr2'],
    ...     offdesign=['zeta1', 'zeta2', 'kA_char']
    ... )

    >>> sg.get_comp('evaporator').set_attr(
    ...     pr1=0.999, ttd_l=20, design=['pr1', 'ttd_l'],
    ...     offdesign=['zeta1', 'kA_char']
    ... )

    >>> sg.get_comp('superheater').set_attr(
    ...     pr1=0.999,  pr2=0.99, design=['pr1', 'pr2'],
    ...     offdesign=['zeta1', 'zeta2', 'kA_char']
    ... )

    >>> sg.get_conn('2').set_attr(td_bubble=5, design=['td_bubble'])

    >>> # %% solve
    >>> # solve design case
    >>> nw.solve('design')
    >>> nw.assert_convergence()
    >>> nw.save('tmp.json')

    >>> # offdesign test
    >>> nw.solve('offdesign', design_path='tmp.json')
    >>> nw.assert_convergence()
