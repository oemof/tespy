Solar collector
---------------

This example shows how you can use a solarthermal collector in TESPy.
The process flow diagram is shown in the image below, the source code can be
found at the `tespy_examples repository
<https://github.com/oemof/oemof-examples/tree/master/oemof_examples/tespy/solar_collector>`_.

.. figure:: api/_images/solar_collector.svg
    :align: center

    Figure: Topology of the solar collector.

The solarthermal collector is used to transfer heat from the solar radiation to
the collector fluid. The TESPy component
:py:class:`tespy.components.heat_exchangers.solar_collector` inherits from
the :py:class:`tespy.components.heat_exchangers.heat_exchanger_simple`
component. An energy balance is applied according to the
:py:meth:`tespy.components.heat_exchangers.solar_collector.energy_func`
method, which takes the collector's

- surface area :code:`A`,
- loss key figures :code:`lkf_lin` (linear) and :code:`lkf_quad` (quadratic),
- ambient temperature :code:`Tamb`,
- optical efficiency :code:`eta_opt` as well as
- incoming radiation :code:`E` (W/m^2)

into account.

In the script different ways of parametrisation are shown. In the last part a
collector is designed and the offdesign performance at different rates of
absorption and ambient temperatures is calculated subsequently. Assuming a
constant mass flow through the collector, the outlet temperature and the
pressure losses of the collector are calculated.

For example, if you want to calculate the performance of the collector within
a specific period of time, you could have the absorbed energy and the ambient
temperature as input time series and iterate over said series. As the absorbed
energy of the collector is a function of the global radiation on the inclined
surface, datetime and location only (optical losses are not temperature
dependent), you could calculate the absorption in a preprocessing script. If
you are to create such a script, we would appreciate you sharing and adding it
to TESPy!
