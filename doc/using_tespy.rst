.. _using_tespy_label:

~~~~~~~~~~~
Using TESPy
~~~~~~~~~~~

TESPy provides a simulation package for component based thermal engineering containing the most important
components of such plants. In the introduction you will learn the basics of modelling component based
plants in TESPy. We then give an overview on the main TESPy modules:

 * tespy.networks,
 * tespy.components,
 * tespy.connections,
 * tespy.network_reader and
 * tespy.tools.

At the end of this page we give a brief overview how TESPy handles fluid properties.

.. figure:: api/_images/tutorial_heat_pump.svg
    :align: center

    Figure 1: Topology of a heat pump.

.. _using_tespy_introduction_label:

We highly recommend to check our :ref:`step by step tutorial <heat_pump_tutorial_label>` on how to
set up a heat pump in TESPy. You will learn, how to set up and design a plant as well as calculate offdesign/partload performance.
Additionally we provide basic examples in the :ref:`examples section <tespy_examples_label>`.

.. contents:: `Contents`
    :depth: 1
    :local:
    :backlinks: top

.. include:: using_tespy/intro.rst
.. include:: using_tespy/networks.rst
.. include:: using_tespy/components.rst
.. include:: using_tespy/subsystems.rst
.. include:: using_tespy/connections.rst
.. include:: using_tespy/fluid_properties.rst
.. include:: using_tespy/other.rst