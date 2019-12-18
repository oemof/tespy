COP of a heat pump
------------------

This example is based on the :ref:`heat pump tutorial <heat_pump_tutorial_label>` and shows how to calculate the COP of a heat pump at different ambient temperatures and different loads of the plant.
The idea is very similar to the :ref:`CHP example <chp_example_label>`, thus you should have a look at the tutorial and the CHP example first.
Be aware, that there are slight changes regarding the topology of the system from the tutorial to this example.
You will find the source code `in this repository <https://github.com/oemof/oemof-examples/tree/master/oemof_examples/tespy/heat_pump>`_.

.. figure:: api/_images/heat_pump_example.svg
    :align: center
	
    Figure: Topology of the heat pump unit.

After the plant has been designed, the consumer's heat demand and the ambient temperature are modified within defined ranges.
Generally, if you are performing offdesign calculation, keep in mind, that a good initial guess/solution is the key to good convergence progress. This is why, we initialise the calculation at a higher ambient temperature with the results
from the calculation of the same load and the nearest ambient temperature possible (in this case always the calculation one step before). This helps the algorithm to stabilize and find a solution.
If you skip out on a large range of temperature or power, you might run into convergence issues. The figures below show the COP of the heat pump for two different heat sources at different temperature levels and at different loads:
In the first figure water is used as heat source, in the second one air. Obviously, the heat pump using air performs much worse. This is mainly to the high power consumption of fan, as for the same amount of heat to be transferred, a much higher volume has to be moved.
    
.. figure:: api/_images/heat_pump_COP_water.svg
    :align: center
	
    Figure: COP of the heat pump using water as heat source. 
	
.. figure:: api/_images/heat_pump_COP_air.svg
    :align: center
	
    Figure: COP of the heat pump using air as heat source.
