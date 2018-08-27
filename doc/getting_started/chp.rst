.. _chp_example_label:


CHP with backpressure turbine
-----------------------------

We have set up a simple combined heat and power unit for this example. A backpressure steam turbine is operated with steam extraction for preheating purposes.
You will find the source code `here <https://github.com/fwitte/tespy_examples/blob/master/chp/chp.py>`_.

.. figure:: api/_images/chp.svg
    :align: center
		
    Figure: Topology of the chp unit.

At first, a plant design is chosen: The unit provides a total power of 5 MW and heating at a temperature of 110 °C for the feed flow.
After that, the temperature at feed flow and live steam mass flow are altered (70 °C to 120 °C and 60 % to 105 % in respect to design mass flow) to cover the unit's range of operation.
Thus, the calculation mode is switched to offdesign and the temperature and mass flow are altered in two embedded loops.
The latest calculated case of the same mass flow is selected for initialisation in order to archieve better and faster convergence.
The results are saved to .csv-files and the following plot of backpressure lines will be created.

.. figure:: api/_images/chp_PQ.svg
    :align: center	
	
    Figure: Backpressure lines of a CHP unit.
