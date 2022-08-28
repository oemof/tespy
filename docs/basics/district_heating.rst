.. _tespy_basics_district_heating_label:

District Heating Network
========================

.. figure:: /_static/images/basics/district_heating.svg
    :align: center
    :alt: Topology of the district heating network
    :figclass: only-light

    Figure: Topology of the district heating network

.. figure:: /_static/images/basics/district_heating_darkmode.svg
    :align: center
    :alt: Topology of the district heating network
    :figclass: only-dark

    Figure: Topology of the district heating network

In this section we provide you with a very simple example as firsts steps in
using TESPy. The model used in this introduction is shown the figure. It
consists of a central heating plant and a consumer, represented by a heat
exchanger with a control valve.

Setting up the System
^^^^^^^^^^^^^^^^^^^^^


Fixed Pipe Dimensions
^^^^^^^^^^^^^^^^^^^^^
- Set pipe parameters, calculate heat loss and pressure loss

Calculate Pipe Dimensions
^^^^^^^^^^^^^^^^^^^^^^^^^
- Set pipe length, pressure loss, calculate diameter
- Set max temperature change per length, calculate kA -> data for insulation
