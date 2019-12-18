.. _tespy_modules_label:

~~~~~~~~~~~~~~~~~
The TESPy modules
~~~~~~~~~~~~~~~~~

TESPy provides a simulation package for component based thermal engineering
containing the most important components of such plants. As reference
literature, the following books can be recommended:

 * Leopold Böswirth, Sabine Bschorer, Thomas Buck: Technische Strömungslehre,
   Wiesbaden, 2012. (german)
 * Hans Dieter Baehr, Stephan Kabelac: Thermodynamik, Berlin, 2016. (german)

In the introduction you will learn the basics of modelling component based
plants in TESPy. We then give an overview on the main TESPy modules, which will
be demonstrated on a very simple example:

 * tespy.networks,
 * tespy.components,
 * tespy.connections,
 * tespy.network_reader and
 * tespy.tools.

At the end of this page we give a brief overview how TESPy handles fluid
properties.

.. contents:: `Contents`
    :depth: 1
    :local:
    :backlinks: top

.. include:: using_tespy/networks.rst
.. _using_tespy_components_label:
.. include:: using_tespy/components.rst
.. include:: using_tespy/subsystems.rst
.. _using_tespy_connections_label:
.. include:: using_tespy/connections.rst
.. _using_tespy_characteristics_label:
.. include:: using_tespy/characteristics.rst
.. include:: using_tespy/fluid_properties.rst
.. include:: using_tespy/other.rst
