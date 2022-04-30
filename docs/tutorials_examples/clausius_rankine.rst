Clausius rankine cycle
----------------------

This example provides a model for a basic clausius rankine cycle.
The process flow diagram is shown in the image below, the source code can be
found at the TESPy `examples repository
<https://github.com/oemof/oemof-examples/tree/master/oemof_examples/tespy/clausius_rankine>`__.

.. figure:: api/_images/basic.svg
    :align: center

    Figure: Topology of the basic clausius rankine cycle.

The basic clausius rankine cycle is built up of a steam turbine, a condenser,
the feed water pump and the steam generator. The ideal process' isentropic
efficiencies of the steam turbine and the pump are at a value of 100 %, and
pressure losses in the condenser and the steam generator are non-existent,
which would result in the thermal efficiency being equal to the Carnot
efficiency. For this example realistic figures have been chosen.
After the plant design an offdesign calculation with 90 % rated power is
performed. The inline-comments give you hints which and why offdesign
parameters have been chosen.

For a good start you can try to modify or exchange parameters. E.g. adjust the
value of the upper terminal temperature difference at the condenser, or replace
this parameter with a pressure at the turbine's outlet. In order to get more
familiar with how TESPy works you could try to insert more components, maybe
add an extraction of the steam turbine for preheating the feed water. It is
strongly recommended to add new components bit per bit, troubleshooting is much
easier this way.
