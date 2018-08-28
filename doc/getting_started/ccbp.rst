.. _chp_example_label:

Combined cycle with backpressure turbine
----------------------------------------

Another example for chp units, this one is a combined cycle power plant using a backpressure steam turbine. Additionally to the heat extraction at the steam turbine condenser,
the distrcit heating water extracts energy from the waste heat of the waste heat steam generator. You will find the source code `here <https://github.com/fwitte/tespy_examples/blob/master/ccbp/cc_bp.py>`_.

.. figure:: api/_images/cc_bp.svg
    :align: center
		
    Figure: Topology of the chp unit.

The example file handles the plant's design as well as the offdesign performance.
