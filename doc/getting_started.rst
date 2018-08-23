~~~~~~~~~~~~~~~
Getting started
~~~~~~~~~~~~~~~

.. _tespy_examples_label:

Examples
========

In the example section we provide a variety of TESPy applications:

* a very basic model of the clausius rankine process,
* the calculation of backpressure lines of a chp at different loads and feed flow temperature levels,
* modeling approach for a district heating system with various consumers and system infrastructure and
* the COP of a heat pump dependent on load and ambient temperature.

You can find all examples in the `tespy examples github repository <https://github.com/fwitte/tespy_examples/>`_.

.. contents:: `Examples`
    :depth: 1
    :local:
    :backlinks: top

.. include:: getting_started/basic.rst
.. include:: getting_started/chp.rst
.. include:: getting_started/ccbp.rst
.. include:: getting_started/district_heating.rst
.. include:: getting_started/heat_pump.rst

.. _tespy_tutorial_label:

Tutorials
=========

We provide two different tutorials for you to better understand how to work with TESPy.
You will learn how to create basic models and get the idea of designing a plant and simulating the offdesign behaviour in the heat pump tutorial
(this is the bases for the COP calculation of the heat pump in the examples).
On top of that, we created a tutorial for the usage of the combustion chamber: It is an important component for thermal power plants while beeing a source for many errors in the calculation.

.. contents:: `Tutorials`
    :depth: 1
    :local:
    :backlinks: top

.. include:: getting_started/tutorial_heat_pump.rst
.. include:: getting_started/tutorial_combustion_chamber.rst
