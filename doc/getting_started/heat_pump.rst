.. _heat_pump_cop_label:

~~~~~~~~~~~~~~~~~~
COP of a heat pump
~~~~~~~~~~~~~~~~~~

This example is based on the :ref:`heat pump tutorial <heat_pump_tutorial_label>` and shows how to calculate the COP of a heat pump at different ambient temperatures and different loads of the plant.
The idea is very similar to the :ref:`CHP example <chp_example_label>`, thus you should have a look at the tutorial and the CHP example first.
You will find the source code `in this repository <https://github.com/fwitte/tespy_examples/blob/master/heat_pump>`_.

.. figure:: https://github.com/fwitte/tespy_examples/tree/master/heat_pump/flow_diagram.svg
    :align: center
    Figure: Topology of the heat pump unit.

After the plant has been designed, the consumer's heat demand and the ambient temperature are modified within defined ranges.
Generally, if you are performing offdesign calculation, keep in mind, that a good initial guess/solution is the key to good convergence progress. This is why, we initialise the calculation at a higher ambient temperature with the results
from the calculation of the same load and the nearest ambient temperature possible (in this case always the calculation one step before). This helps the algorithm to stabilize and find a solution.
If you skip out on a large range of temperature or power, you might run into convergence issues. The figure below shows the COP of the heat pump for the different temperature levels and at different loads.
    
.. figure:: https://github.com/fwitte/tespy_examples/tree/master/heat_pump/COP.svg
    :align: center
    Figure: COP of the heat pump.

