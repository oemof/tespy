.. _basic_example_label:

Clausius rankine cycle
----------------------

This example provides a model for a basic clausius rankine cycle.
The process flow diagram is shown in the image below, the source code can be found at the `tespy_examples repository <https://github.com/fwitte/tespy_examples/blob/master/basic/clausius_rankine.py>`_.

.. figure:: api/_images/basic.svg
    :align: center
	
    Figure: Topology of the basic clausius rankine cycle.

The basic clausius rankine cycle is built up of a steam turbine, a condenser, the feed water pump and the steam generator. The ideal process' isentropic efficiencies of the steam turbine and the pump are at a value of 100 %, and pressure losses in the condenser and the steam generator are non-existent, which would result in the thermal efficiency beeing equal to the carnot efficiency. For this example realistic figures have been chosen.
After the plant design an offdesign calculation with 90 % rated power is performed. The inline-comments give you hints which and why offdesign parameters have been choosen. Additionally we added a calculation of thermal efficiency with the standard method as well as with the entropy method:

.. math::

    \eta_{th} = \frac{|\sum P|}{\sum \dot{Q}_{in}}

    \eta_{th} = \eta_c - \sum_i \Delta \eta_i

    T_{m,in} = \frac{\sum \dot{Q}_{in}}{\sum \dot{S}_{in}}

    T_{m,out} = \frac{\sum \dot{Q}_{out}}{\sum \dot{S}_{out}}

    \eta_c = 1 - \frac{T_{m,out}}{T_{m,in}}

    \Delta \eta_i = \frac{T_{m,out} \cdot \dot{S}_{irr,i}}
    {\sum \dot{Q}_{in}}

    \text{indices: in = input, out = output}
    
.. figure:: api/_images/basic_efficiency.svg
    :align: center
	
    Figure: Efficiency of the basic clausius rankine cycle in design case.    
    
For a good start you can try to modify or exchange parameters. E. g. adjust the value of the upper terminal temperature difference at the condenser, or replace this parameter with a pressure at the turbine's outlet. In oder to get more familiar with how TESPy works you could try to insert more components, maybe add an extraction of the steam turbine for preheating the feed water. It is strongly recommended to add new components bit per bit, troubleshooting is much easier this way.
