.. _solar_collector_example_label:

Solar collector
---------------

This example shows how you can use a solarthermal collector in TESPy.
The process flow diagram is shown in the image below, the source code can be found at the `tespy_examples repository <https://github.com/fwitte/tespy_examples/blob/master/solar_collector/solar_collector.py>`_.

.. figure:: api/_images/solar_collector.svg
    :align: center
	
    Figure: Topology of the solar collector.

The solarthermal collector is used to transfer heat from the solar radiation to the collector fluid.
The TESPy component :py:class:`solar_collector <tespy.components.components.solar_collector>` inherits from the :py:class:`simple_heat_exchanger <tespy.components.components.simple_heat_exchanger>` component.
An energy balance is applied according to the :py:class:`solar collector energy func <tespy.components.components.solar_collector.energy_func>` method, which takes the collector's

- surface area :code:`A`,
- loss key figures :code:`lkf_lin` (linear) and :code:`lkf_quad` (quadratic),
- ambient temperature :code:`t_a` as well as
- area independent absorped energy :code:`E` (radiation on inclined surface minus optical losses)

into account.

In the script different ways of parametrisation are shown. In the last part a collector is designed and the offdesign performance at different rates of absorption and ambient temperatures is calculated subsequently.
Assuming a constant mass flow through the collector, the outlet temperature and the pressure losses of the collector are calculated.

E. g., if you want to calculate the performance of the collector within a specifc period of time, you could have the absorped energy and the ambient temperature as input time series and iterate over said series.
As the absorped energy of the collector is a function of the global radiation on the inclined surface, datetime and location only (optiacal losses are not temperature dependet), you could calculate the absorption in a preprocessing script.
If you are to create such a script, we would appreciate you sharing and adding it to TESPy!
