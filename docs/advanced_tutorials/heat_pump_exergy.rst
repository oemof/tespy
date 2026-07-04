.. _tutorial_heat_pump_exergy_label:

Exergy Analysis of a Ground-Coupled Heat Pump
---------------------------------------------

.. note::

    The exergy analysis in this tutorial uses
    `exerpy <https://github.com/oemof/exerpy>`__, the dedicated external
    library for exergy analysis. Exerpy is fully compatible with TESPy models -
    build your network as usual and pass the system boundary crossing streams to
    exerpy for automatic analysis. On top of physical exergy, exerpy offers
    advanced features like exergoeconomic methods. More examples are available
    in the `exerpy documentation <https://exerpy.readthedocs.io/en/latest/examples.html>`__.

Task
^^^^

This tutorial shows how to set up and carry out an exergy analysis for a
ground-coupled heat pump (GCHP). In addition, various post-processing options
are presented. To investigate the impact of refrigerant choice on COP and
exergetic efficiency, the same :code:`HeatPumpModel` class is evaluated with
different refrigerants (NH3 and R290). Finally, the influence of varying
different parameters on COP and exergy efficiency is investigated and plotted.

.. note::

    Please note, currently this tutorial is intended to show the user, how to
    carry out an exergy analysis for a simple system and how to use this
    toolbox in several investigations of a specific system. While there is a
    very short description of the setup, methodology and results, an in-depth
    discussion of the method and the results is not yet provided. If you would
    like to add this to the documentation you are welcome to contact us via our
    GitHub.

Since there is an existing tutorial for
:ref:`creating a heat pump <tutorial_heat_pump_label>`, this tutorial
starts with the explanations for setting up the exergy analysis. Note however,
that the heat pump model differs slightly in structure from the model in the
previous tutorial. All related Python scripts of the fully working GCHP-model
are listed in the following:

- GCHP model class:
  :download:`heat_pump_model.py </../tutorial/heat_pump_exergy/heat_pump_model.py>`
- Single-fluid evaluation (NH3):
  :download:`NH3_example.py </../tutorial/heat_pump_exergy/NH3_example.py>`
- Full calculations for all refrigerants with parametric studies:
  :download:`all_calculations.py </../tutorial/heat_pump_exergy/all_calculations.py>`
- Plots of the results of the parameter variations:
  :download:`plots.py </../tutorial/heat_pump_exergy/plots.py>`

The figure below shows the topology of the GCHP. In this model, a
ground-coupled heat pump is modeled, which is for instance connected to a
single-family house with underfloor heating. The heating system represents the
heat demand of the house. The geothermal heat collector is represented by a
ground heat feed flow (Source) and return flow (Sink). The heat pump circuit
consists of the basic components: condenser, expansion valve, evaporator and
compressor.

.. figure:: /_static/images/tutorials/heat_pump_exergy/flowsheet.svg
    :align: center
    :alt: Topology of the Ground-Couped Heat Pump (GCHP)
    :figclass: only-light

    Figure: Topology of the Ground-Couped Heat Pump (GCHP).

.. figure:: /_static/images/tutorials/heat_pump_exergy/flowsheet_darkmode.svg
    :align: center
    :alt: Topology of the Ground-Couped Heat Pump (GCHP)
    :figclass: only-dark

    Figure: Topology of the Ground-Couped Heat Pump (GCHP).

The input data of the model are based on different literature. In general, the
model of the GCHP is based on a data sheet of a real heat pump
(`Viessmann Vitocal 300-G <https://www.viessmann.de/de/wohngebaeude/waermepumpe/sole-wasser-waermepumpen/vitocal-300-g.html>`_ ).
However, the data are used as approximate values to create a model that works
with both NH3 and R290, although the mentioned heat pump is designed to use
R410A. The range of the underfloor heating system temperature and the range of
the geothermal temperature are assumptions based on measured data from the
research project
`WPsmart <https://wp-monitoring.ise.fraunhofer.de/wp-smart-im-bestand/german/index/index.html>`_
and :cite:`Chen2015`. The average outdoor temperature is taken from
:cite:`Chen2015`.

TESPy model
^^^^^^^^^^^

In principle, the GCHP-model corresponds to the flowsheet shown above.
The heating system and the geothermal heat collector can be modeled as sources
and sinks, which represent the feed and the return flow in both cases.
The condenser is modeled as :code:`MovingBoundaryHeatExchanger` instance, while
the evaporator is modeled using a :code:`HeatExchanger` instance. In total, the
TESPy model consists of 11 components.

The model is encapsulated in a :code:`HeatPumpModel` class that inherits from
:code:`ModelTemplate`. This design allows the same network to be instantiated
for different refrigerants by simply passing the fluid name as a constructor
argument, and provides a high-level interface for parametric studies through
the :code:`sensitivity_analysis` method. The constructor stores the design
parameters and delegates network creation to the :code:`_create_network` method
via the :code:`ModelTemplate` base class:

.. literalinclude:: /../tutorial/heat_pump_exergy/heat_pump_model.py
   :language: python
   :pyobject: HeatPumpModel.__init__

The :code:`_parameter_lookup` method defines the named parameters that can be
read and written through the :code:`set_parameters` and
:code:`sensitivity_analysis` interfaces. For parameters that require custom
getter or setter logic - such as :code:`T_geo`, which maps to the feed
temperature with a fixed offset, or :code:`COP` and :code:`epsilon` which are
derived quantities - a dictionary with :code:`"get"` and/or :code:`"set"` keys
is used. Simple connection or component attributes can be given as a path list
directly:

.. literalinclude:: /../tutorial/heat_pump_exergy/heat_pump_model.py
   :language: python
   :pyobject: HeatPumpModel._parameter_lookup

In real systems, the circulating brine in the geothermal collector usually
consists of a mixture of water and antifreeze. Pure water is used as the
circulating fluid in this example. In fact, some geothermal collectors are
filled with water, provided that the ground temperature is high enough
throughout the year, such as in :cite:`Chen2015`.

The following parameter specifications were made for the design case
calculation:

- isentropic efficiency values
- electrical conversion efficiencies of compressor and pumps
- terminal temperature difference values at condenser and evaporator
- pressure losses in condenser and evaporator
- hot and cold side heat transfer coefficients of evaporator
- temperature difference to boiling point of refrigerant at compressor inlet
- temperatures and pressure of heating system feed and return flow
- temperatures and pressure of geothermal heat collector feed and return flow
- condenser heat output

The network is built, parametrized and solved within the :code:`_create_network`
method. The motor and power distribution objects are also set up here, so that
power connections cross the system boundary in a way that is compatible with
exerpy. The geothermal return temperature :code:`c13` tracks the feed
temperature :code:`c11` via a :code:`Ref` object with a fixed offset of 3 °C.
Similarly, the heating system return temperature :code:`c21` tracks the feed
temperature :code:`c23` with a fixed offset of 5 °C.

.. literalinclude:: /../tutorial/heat_pump_exergy/heat_pump_model.py
   :language: python
   :pyobject: HeatPumpModel._create_network

The units used are temperature in °C, pressure in bar and enthalpy in kJ/kg.
The ambient state (:code:`Tamb`, :code:`pamb`) is stored on the instance and
passed to exerpy when the exergy analysis is run.

h-log(p)-diagram
^^^^^^^^^^^^^^^^

At first, we will have a short look at the h-log(p)-diagram of the process,
exemplary for NH3 as working fluid. Such diagrams are useful to better
understand a process, therefore we will quickly present how to generate it
using TESPy. The :code:`plot_logph_diagram_matplotlib` method wraps the
`fluprodia <https://fluprodia.readthedocs.io/en/latest/>`_ library to generate
the diagram directly from the network state. The cycle-closing connection
:code:`'c1'` is passed as the starting point so that the method can trace the
full refrigerant cycle:

.. literalinclude:: /../tutorial/heat_pump_exergy/NH3_example.py
   :language: python
   :start-after: [logph]
   :end-before: [exergy]

.. note::

    For more information on fluprodia integration also see
    :ref:`here <fluprodia_label>`.

.. figure:: /_static/images/tutorials/heat_pump_exergy/NH3_logph.svg
    :align: center
    :alt: Fluid Property Diagram h-log(p) of the GCHP
    :figclass: only-light

    Figure: h-log(p) diagram of the NH3 GCHP.

.. figure:: /_static/images/tutorials/heat_pump_exergy/NH3_logph_darkmode.svg
    :align: center
    :alt: Fluid Property Diagram h-log(p) of the GCHP
    :figclass: only-dark

    Figure: h-log(p) diagram of the NH3 GCHP.

The resulting fluid property diagram is shown in the figure above. It can
easily be seen, that the evaporator slightly overheats the working fluid, while
it leaves the condenser in saturated liquid state. The working fluid temperature
after leaving the compressor is quite high with far more than 100 °C given the
heat sink only requires a temperature of only 40 °C. In comparison, R290
leaves the compressor at a lower temperature.

More examples of creating fluid property diagrams can be found in the fluprodia
documentation referenced above.

Exergy analysis
^^^^^^^^^^^^^^^
Following, the main tasks of this tutorial are presented. First, the exergy
analysis is set up for the respective network and carried out for the base
case. Subsequently, the influence of different parameters such as temperature
of the heat source and sink as well as ambient temperature and part load
operation of the heat pump regarding exergetic efficiency are investigated.

Analysis setup
++++++++++++++

After the network has been solved, the exergy analysis is carried out via the
:code:`run_exergy_analysis` method, which internally creates an
:py:class:`exerpy.ExergyAnalysis` instance from the TESPy network. All exergy
streams crossing the system boundary must be classified as:

- fuel exergy :code:`E_F` - resources supplied to the system
- product exergy :code:`E_P` - desired output of the system
- exergy loss streams :code:`E_L` - exergy discarded to the environment

In exerpy, each of these is a dictionary with :code:`"inputs"` and
:code:`"outputs"` keys containing the labels of the boundary-crossing
connections.

For the GCHP the electrical power is supplied via a :code:`PowerConnection`
labelled :code:`'e1'` (grid side). The geothermal heat boundary is represented
by the material connections :code:`'c11'` (inlet from ground) and
:code:`'c13'` (outlet to ground). The heating system boundary is represented
by :code:`'c23'` (feed flow to house) and :code:`'c21'` (return flow from
house). In the example of the GCHP, only :code:`E_F` and :code:`E_P` are
defined. Ambient temperature and pressure are passed in the network units
(°C and bar); the method converts them to SI units (K, Pa) before calling
exerpy:

.. literalinclude:: /../tutorial/heat_pump_exergy/NH3_example.py
   :language: python
   :start-after: [exergy]

Results
+++++++

The results can be printed and retrieved as DataFrames using the
:py:meth:`exerpy.ExergyAnalysis.exergy_results` method:

.. code-block:: python

   df_comp, df_material, df_power = ean.exergy_results()

The overall system results (total :code:`E_F`, :code:`E_P`, :code:`E_D` and
:code:`epsilon`) are available directly as attributes:

.. code-block:: python

   print(f"E_F = {ean.E_F:.1f} W")
   print(f"E_P = {ean.E_P:.1f} W")
   print(f"epsilon = {ean.epsilon:.3f}")

An exergy destruction waterfall diagram can be generated with:

.. code-block:: python

   ean.plot_exergy_waterfall(title='NH3 Heat Pump Exergy Analysis')

Parametric analysis
^^^^^^^^^^^^^^^^^^^
Below, different parametric analyses will be presented considering the
following issues:

- plot exergy destruction
- varying ambient and geothermal temperature
- varying geothermal and heating system temperature
- varying heating load and geothermal temperature

In order to be able to compare the results of the two refrigerants NH3 and
R290, all calculations are collected in a single script
:download:`all_calculations.py </../tutorial/heat_pump_exergy/all_calculations.py>`
that loops over both fluids. The :code:`HeatPumpModel` class makes it
straightforward to switch refrigerants - only the fluid name changes. The plots
in this tutorial are created with `Matplotlib <https://matplotlib.org/>`_ in a
separate script :download:`plots.py </../tutorial/heat_pump_exergy/plots.py>`.
For installation instructions or further documentation please see the Matplotlib
documentation.

For the post-processing, the following additional packages
are required:

.. code-block:: python

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

The overall structure of the calculations script is a loop over both working
fluids. For each fluid the model is created, the design case is solved and
saved, the h-log(p) diagram is generated and the exergy analysis is run at the
design point before the parametric studies begin:

.. literalinclude:: /../tutorial/heat_pump_exergy/all_calculations.py
   :language: python
   :end-before: [ed_export]

.. note::

    All code excerpts shown in the following subsections are **continuations
    of the same** :code:`for fluid in ["NH3", "R290"]:` loop introduced
    above. They are not standalone scripts.

Plot exergy destruction
+++++++++++++++++++++++
In order to visualize how much exergy of the fuel exergy :code:`E_F` the
individual components of the GCHP destroy, the exergy destruction :code:`E_D`
can be displayed in a bar chart as shown at the end of this section. The
waterfall diagram is created directly from the :code:`ExergyAnalysis` instance
by calling :py:meth:`exerpy.ExergyAnalysis.plot_exergy_waterfall`.

.. figure:: /_static/images/tutorials/heat_pump_exergy/NH3_waterfall.svg
    :align: center
    :alt: Waterfall plot for ammonia heat pump
    :figclass: only-light

    Figure: Waterfall diagram for the ammonia heat pump

.. figure:: /_static/images/tutorials/heat_pump_exergy/NH3_waterfall_darkmode.svg
    :align: center
    :alt: Waterfall plot for ammonia heat pump
    :figclass: only-dark

    Figure: Waterfall diagram for the ammonia heat pump

In addition, the component-level results are exported to a :code:`.csv` file
so that they can be combined across refrigerants in the separate plot script.
The code below extracts :code:`E_D` and the running fuel-exergy remainder for
each component that destroys more than 1 W, and saves a compact DataFrame:

.. literalinclude:: /../tutorial/heat_pump_exergy/all_calculations.py
   :language: python
   :start-after: [ed_export]
   :end-before: [parametric]

.. note::

    In order to be able to use the data from the data frames in a separate
    script for plot creation, all data frames must be saved as a file with
    their own individual name.

In the separate plot script
(:download:`plots.py </../tutorial/heat_pump_exergy/plots.py>`) the
:code:`.csv` files can now be re-imported to create plots with Matplotlib. The
Python code for creating the bar chart is included in the previously
referenced plot script and can be found there. For more information on
creating plots with Matplotlib, please check the
`Matplotlib documentation <https://matplotlib.org/>`_. The resulting bar chart
is shown below.

.. figure:: /_static/images/tutorials/heat_pump_exergy/diagram_E_D.svg
    :align: center
    :alt: Comparison of exergy destruction and exergy efficiency
    :figclass: only-light

    Figure: Comparison of exergy destruction and exergy efficiency of both
    working fluids in design case.

.. figure:: /_static/images/tutorials/heat_pump_exergy/diagram_E_D_darkmode.svg
    :align: center
    :alt: Comparison of exergy destruction and exergy efficiency
    :figclass: only-dark

    Figure: Comparison of exergy destruction and exergy efficiency of both
    working fluids in design case.

The bar chart shows how much exergy the individual components of the GCHP
destroy in absolute terms and as a percentage of the fuel exergy :code:`E_F`.
After deducting the destroyed exergy :code:`E_D`, the product exergy
:code:`E_P` remains. Overall, it is noticeable that the GCHP with NH3 requires
less fuel exergy than the GCHP with R290, with the same amount of product
exergy. Furthermore, with NH3 the condenser has the highest exergy destruction,
whereas with R290 the valve destroys the largest amount of exergy.

Varying ambient and geothermal temperature
++++++++++++++++++++++++++++++++++++++++++
In order to consider the influence of a change in ambient temperature or
geothermal temperature on the exergetic efficiency, parametric studies are
performed with different values of these parameters.

For the variation of the ambient temperature :code:`Tamb`, only the exergy
analysis is re-executed without re-solving the network - the thermodynamic
state is unchanged and only the reference temperature shifts. The ambient
temperature is varied between 4°C and 20°C.

The mean geothermal temperature :code:`Tgeo` is varied between 14°C and 8°C
via offdesign calculations using the
:py:meth:`~tespy.models.template.ModelTemplate.sensitivity_analysis` method.
The method accepts a :code:`param_dict` with the parameter name and a list of
values, a list of result quantities to collect and the solve mode. When a
:code:`postproc_func` is provided, it is called after each successful solve -
here it runs the exergy analysis so that :code:`epsilon` is up to date before
the results are recorded.

.. literalinclude:: /../tutorial/heat_pump_exergy/all_calculations.py
   :language: python
   :start-after: [parametric]
   :end-before: [tgeo_ths]

The results of the calculation can be plotted as shown in the following
figure. The related Python code to create this plot can be found in the plot
script (:download:`plots.py </../tutorial/heat_pump_exergy/plots.py>`). For
further documentation please see the `Matplotlib <https://matplotlib.org/>`__
documentation.

.. figure:: /_static/images/tutorials/heat_pump_exergy/diagram_eps_Tamb_Tgeo.svg
    :align: center
    :alt: Varying Tamb and Tgeo of the GCHP
    :figclass: only-light

    Figure: Varying ambient and geothermal temperature.

.. figure:: /_static/images/tutorials/heat_pump_exergy/diagram_eps_Tamb_Tgeo_darkmode.svg
    :align: center
    :alt: Varying Tamb and Tgeo of the GCHP
    :figclass: only-dark

    Figure: Varying ambient and geothermal temperature.

It can be recognized that the specified ambient temperature :code:`Tamb` used
in the exergy analysis has a considerable influence on the exergetic efficiency
epsilon. The closer the ambient temperature is to the temperature of the
heating system, the lower the exergetic efficiency. This can be argued from
the fact that while :code:`E_F` and :code:`E_P` both decrease with increasing
:code:`Tamb`, :code:`E_P` decreases proportionally more than :code:`E_F`. In
comparison, it can be seen on the right that with increasing :code:`Tgeo`, and
thus decreasing temperature difference between geothermal heat collector and
heating system, epsilon increases. This can be explained by the resulting
decrease in :code:`E_F` with :code:`E_P` remaining constant.

Varying geothermal and heating system temperature
+++++++++++++++++++++++++++++++++++++++++++++++++
Another relation that can be investigated is the influence of a change in the
geothermal and the heating system temperatures on the exergetic efficiency and
the COP of the GCHP. In this calculation :code:`Tgeo` is varied between 14°C
and 10°C. The heating system temperature :code:`Ths` is varied between 45°C
and 35°C. All temperature values are mean values of the feed and return flow
temperatures.

The full Cartesian product of :code:`Tgeo_range` and :code:`Ths_range` is
assembled with :func:`itertools.product` and passed to
:py:meth:`~tespy.models.template.ModelTemplate.sensitivity_analysis` as
parallel lists. Both :code:`T_geo` and :code:`T_hs` are varied simultaneously
within a single offdesign loop, and the results DataFrame is then pivoted to
obtain COP and exergetic efficiency as functions of both temperatures:

.. literalinclude:: /../tutorial/heat_pump_exergy/all_calculations.py
   :language: python
   :start-after: [tgeo_ths]
   :end-before: [tgeo_q]

The results of this calculation are shown in the following figure. The
corresponding Python code can likewise be found in the plot script
(:download:`plots.py </../tutorial/heat_pump_exergy/plots.py>`).

.. figure:: /_static/images/tutorials/heat_pump_exergy/diagram_cop_eps_Tgeo_Ths.svg
    :align: center
    :alt: Varying Tgeo and Ths of the GCHP
    :figclass: only-light

    Figure: Varying geothermal and heating system temperature.

.. figure:: /_static/images/tutorials/heat_pump_exergy/diagram_cop_eps_Tgeo_Ths_darkmode.svg
    :align: center
    :alt: Varying Tgeo and Ths of the GCHP
    :figclass: only-dark

    Figure: Varying geothermal and heating system temperature.

It can be seen that the GCHP with NH3 has a better exergetic efficiency than
with R290. As in the prior investigation, an increasing geothermal heat
collector temperature also has a favorable effect on epsilon. The opposite
behavior of epsilon and COP for both refrigerants is remarkable. The COP drops
while the exergetic efficiency rises. This can be explained by the fact that at
constant heating load :code:`Q`, the required electrical power input increases
as the heating system temperature rises. However regarding exergetic
efficiency, :code:`E_F` and :code:`E_P` both increase with increasing heating
system temperature. The ratio between these two parameters is such that
the exergetic efficiency improves as the heating system temperature rises.

Varying geothermal temperature and heating load
+++++++++++++++++++++++++++++++++++++++++++++++
Finally, the influence of the simultaneous variation of the geothermal
temperature :code:`Tgeo` and the heating load :code:`Q` on the exergetic
efficiency and the COP of the GCHP is examined. The investigation is carried
out in the same way as the variation of :code:`Tgeo` and :code:`Ths` described
above. In contrast to the previous investigation, :code:`Q` is varied here
instead of :code:`Ths`. The range of :code:`Q` varies between 4.3 and 2.8 kW.
The rated load was previously set at 4 kW in the design calculation.

.. literalinclude:: /../tutorial/heat_pump_exergy/all_calculations.py
   :language: python
   :start-after: [tgeo_q]

.. figure:: /_static/images/tutorials/heat_pump_exergy/diagram_cop_eps_Tgeo_Q.svg
    :align: center
    :alt: Varying Tgeo and Q of the GCHP
    :figclass: only-light

    Figure: Varying geothermal temperature and heat load.

.. figure:: /_static/images/tutorials/heat_pump_exergy/diagram_cop_eps_Tgeo_Q_darkmode.svg
    :align: center
    :alt: Varying Tgeo and Q of the GCHP
    :figclass: only-dark

    Figure: Varying geothermal temperature and heat load.

The results are shown in the figure above. As before, the Python code for
creating the plot can be found in the plot script
(:download:`plots.py </../tutorial/heat_pump_exergy/plots.py>`).
The partial load behavior of the GCHP, which results from the characteristic
lines of the efficiencies of the individual components, can be recognized
in the curves shown.

Conclusion
^^^^^^^^^^
This tutorial provides an exemplary insight into post-processing with the
TESPy exergy analysis tool. Of course, other parameters can also be examined
and varied. Feel free to try out different parameter variations. But make sure
that the data ranges are not only adjusted in the Python script of the model,
but also in the Python script of the plots, if a plot is created with the
stand-alone plot script.

More examples of exergy analysis can be found in the
:ref:`TESPy analysis section <advanced_exergy_label>` and in the
`exerpy documentation <https://exerpy.readthedocs.io/en/latest/examples.html>`__.
If you are interested in contributing or have questions and remarks on this
tutorial, you are welcome to file an issue at our GitHub page.
