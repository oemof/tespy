.. _tespy_tutorial_heat_pump_exergy_label:

Exergy Analysis of a Ground-Coupled Heat Pump
---------------------------------------------

Task
^^^^

This tutorial shows how to set up and carry out an exergy analysis for a
ground-coupled heat pump (GCHP). In addition, various post-processing options
are presented. To investigate the impact of refrigerant choice on COP and
exergetic efficiency, two Python scripts of the same network with different
refrigerants (NH3 and R410A) are created. Finally, the influence of varying
different parameters on COP and exergetic efficiency is investigated and
plotted.

.. note::

    Please note, currently this tutorial is intended to show the user, how to
    carry out an exergy analysis for a simple system and how to use this
    toolbox in several investigations of a specific system. While there is a
    very short description of the setup, methodology and results, an in-depth
    discussion of the method and the results is not yet provided. If you would
    like to add this to the documentation you are welcome to contact us via our
    GitHub.

Since there is an existing tutorial for
:ref:`creating a heat pump <heat_pump_tutorial_label>`, this tutorial starts
with the explanations for setting up the exergy analysis. Note, however, that
the heat pump model differs slightly in structure from the model in the
previous tutorial. All related Python scripts of the fully working GCHP-model
are listed in the following:

- GCHP with NH3 (the model only):
  :download:`NH3.py </../tutorial/heat_pump_exergy/NH3.py>`
- GCHP with R410A (the model only):
  :download:`R410A.py </../tutorialheat_pump_exergy/R410A.py>`
- GCHP with NH3 (model and post-processing):
  :download:`NH3_calculations.py </../tutorial/heat_pump_exergy/NH3_calculations.py>`
- GCHP with R410A (model and post-processing):
  :download:`R410A_calculations.py </../tutorial/heat_pump_exergy/R410A_calculations.py>`
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
with both NH3 and R410A, although the mentioned heat pump is designed to use
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
The condenser is modeled as :code:`Condenser` instance, while the evaporator
is modeled using :code:`HeatExchanger` instance. In total, the TESPy model
consists of 11 components.

In real systems, the circulating brine in the geothermal collector usually
consists of a mixture of water and antifreeze. Since the calculation with
mixtures of incompressible fluids is not yet fully implemented in TESPy, pure
water is used as the circulating fluid in this network. In fact, some
geothermal collectors are filled with water, provided that the ground
temperature is high enough throughout the year, such as in :cite:`Chen2015`.

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

The model using NH3 as refrigerant and the model using R410A as refrigerant
differ in the fluid definition, the naming of the stored files and the
specification of the starting values only. The definition of the starting
values is necessary to obtain a numerical solution for the first calculation.
In this tutorial, the given code examples are shown exemplary for the model
with NH3 as refrigerant only.

The units used and the ambient state are defined as follows:

.. code-block:: python

    nw = Network(
        fluids=['water', 'NH3'], T_unit='C', p_unit='bar',
        h_unit='kJ / kg', m_unit='kg / s'
    )

    pamb = 1.013
    Tamb = 2.8

For the model using R410A as refrigerant, the fluid definition is accordingly
:code:`'R410A'` instead of :code:`'NH3'`.

The temperature of the heating system feed flow is set to 40°C in design
calculation. The difference between feed and return flow temperature is kept
constant at 5°C. Therefore the return flow is set to 35°C.

The geothermal heat collector temperature is defined as follows:

.. code-block:: python

    Tgeo = 9.5

:code:`Tgeo` is the mean geothermal temperature. The difference between
feed and return flow temperature is kept constant at 3°C. Therefore, the feed
flow temperature in the design calculation is set to :code:`Tgeo + 1.5°C` and
the return flow temperature is set to :code:`Tgeo - 1.5°C`.

The complete Python code of the TESPy models is available in the scripts
:download:`NH3.py </../tutorial/NH3.py>` with NH3 as refrigerant and
:download:`R410A.py </../tutorial/R410A.py>` with R410A as refrigerant. All
other specified values of the component and connection parameters can be found
in these Python scripts.

In the scripts
:download:`NH3_calculations.py </../tutorial/NH3_calculations.py>` and
:download:`R410A_calculations.py </../tutorial/R410A_calculations.py>`,
the Python code of the TESPy models of the GCHP is extended to handle the
different tasks mentioned in the introduction. In these two scripts you can
find the corresponding Python code for all calculations that will be presented
in the next sections of the tutorial. As previously mentioned, the given code
examples in the following are only shown exemplary for the GCHP with NH3 as
refrigerant. If the scripts differ beyond the mentioned points, it will be
pointed out at the respective place of the tutorial.

h-log(p)-diagram
^^^^^^^^^^^^^^^^

At first, we will have a short look at the h-log(p)-diagram of the process,
exemplary for NH3 as working fluid. Such diagrams are useful to better
understand a process, therefore we will quickly present how to generate it
using TESPy with fluprodia. For more information and installation
instructions for fluprodia please have a look at the
`online documentation <https://fluprodia.readthedocs.io/en/latest/>`_.

The data for the diagram are first saved in a dictionary :code:`result_dict`
using the :code:`get_plotting_data` method of each component that is to be
visualized.

.. code-block:: python

    from fluprodia import FluidPropertyDiagram

    result_dict = {}
    result_dict.update({ev.label : ev.get_plotting_data()[2]})
    result_dict.update({cp.label : cp.get_plotting_data()[1]})
    result_dict.update({cd.label : cd.get_plotting_data()[1]})
    result_dict.update({va.label : va.get_plotting_data()[1]})

.. note::

    The first level key of the nested dictionary returned from the
    :code:`get_plotting_data` method contains the connection id of the state
    change. Make sure you specify the correct id for the components to be
    displayed. A table of the state change and the respective id can be found
    :ref:`here <FluProDia_label>`.

Next, a :code:`FluidPropertyDiagram` instance is created and the units of the
diagram are specified.

.. code-block:: python

    diagram = FluidPropertyDiagram('NH3')
    diagram.set_unit_system(T='°C', p='bar', h='kJ/kg')

Afterwards, the dictionary can be passed to the :code:`calc_individual_isoline`
method of the :code:`FluidPropertyDiagram` object. In addition, the axis
limits are set. The :code:`calc_isolines` method calculates all isolines of the
diagram and the :code:`draw_isolines` method draws the isolines of the
specified type. Finally, the results can be plotted and the diagram can be
saved with the code shown below.

.. code-block:: python

    for key, data in result_dict.items():
            result_dict[key]['datapoints'] = diagram.calc_individual_isoline(**data)

    diagram.set_limits(x_min=0, x_max=2100, y_min=1e0, y_max=2e2)
    diagram.calc_isolines()
    diagram.draw_isolines('logph')

    for key in result_dict.keys():
        datapoints = result_dict[key]['datapoints']
        diagram.ax.plot(datapoints['h'],datapoints['p'], color='#ff0000')
        diagram.ax.scatter(datapoints['h'][0],datapoints['p'][0], color='#ff0000')

    diagram.save('NH3_logph.svg')

.. figure:: /_static/images/tutorials/heat_pump_exergy/NH3_logph.svg
    :align: center
    :alt: Fluid Property Diagram h-log(p) of the GCHP

    Figure: h-log(p) diagram of the NH3 GCHP.

The resulting fluid property diagram is shown in the figure above. It can
easily be seen, that the evaporator slightly overheats the working fluid, while
the it leaves the condenser in saturated liquid state. The working fluid
temperature after leaving the compressor is quite high with far more than
100 °C given the heat sink only requires a temperature of only 40 °C. In
comparison, the R410A leaves the compressor at about 75 °C.

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

After the network has been built, the exergy analysis can be set up. For this
purpose, all exergy flows entering and leaving the network must be defined.
The exergy flows are defined as a list of busses as follows:

- fuel exergy :code:`E_F`
- product exergy :code:`E_P`
- exergy loss streams :code:`E_L`
- internal exergy streams not bound to connections :code:`internal_busses`

First, the busses for the exergy analysis must be defined. The first bus is
for the electrical energy supply of the compressor and the pumps. The motor
efficiency is calculated by a characteristic line. This power input bus
represents fuel exergy.

The product exergy is the heat supply of the condenser to the heating system,
which is represented by the heating system bus. The bus consists of the
streams :code:`hs_ret` and :code:`hs_feed`. Note that the :code:`base`
keyword of the stream entering the network :code:`hs_ret` must be set to
:code:`bus`.

Lastly, the geothermal heat bus represents the heat that is transferred from
the geothermal heat collector to the evaporator. The bus consists of the
streams :code:`gh_in` and :code:`gh_out`. Here, the :code:`base` of the stream
:code:`gh_in` is set to :code:`bus`, because this stream represents an energy
input from outside of the network. In this example, the geothermal heat bus is
defined as fuel exergy, because the ambient temperature :code:`Tamb` is set at
a lower temperature than the temperature of the geothermal heat collector.

.. code-block:: python

    x = np.array([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2])
    y = np.array([0, 0.86, 0.9, 0.93, 0.95, 0.96, 0.95])

    char = CharLine(x=x, y=y)
    power = Bus('power input')
    power.add_comps({'comp': cp, 'char': char, 'base': 'bus'},
                    {'comp': ghp, 'char': char, 'base': 'bus'},
                    {'comp': hsp, 'char': char, 'base': 'bus'})

    heat_cons = Bus('heating system')
    heat_cons.add_comps({'comp': hs_ret, 'base': 'bus'}, {'comp': hs_feed})

    heat_geo = Bus('geothermal heat')
    heat_geo.add_comps({'comp': gh_in, 'base': 'bus'},
                       {'comp': gh_out})

    nw.add_busses(power, heat_cons, heat_geo)

In order to carry out the exergy analysis an :code:`ExergyAnalysis` instance
passing the network to analyse as well as the respective busses is created.
The product exergy is represented by the bus :code:`power`. The busses
:code:`heat_cons` and :code:`heat_geo` are passed as fuel exergy.
In the example of the GCHP, only :code:`E_F` and :code:`E_P` are defined.
Other examples of exergy analysis setup can be found in the
:ref:`TESPy analysis <tespy_advanced_exergy_label>` page and in the API
documentation of class :py:class:`tespy.tools.analyses.ExergyAnalysis`.

.. code-block:: python

   ean = ExergyAnalysis(network=nw,
                        E_F=[power, heat_geo],
                        E_P=[heat_cons])

   ean.analyse(pamb, Tamb)

The :py:meth:`tespy.tools.analyses.ExergyAnalysis.analyse` method will run the
exergy analysis automatically. This method expects information about the
ambient pressure and ambient temperature. Additionally, an automatic check of
consistency is performed by the analysis as further described in
:ref:`TESPy analysis <tespy_advanced_exergy_label>`.

Results
+++++++

The results can be printed by using the
:py:meth:`tespy.tools.analyses.ExergyAnalysis.print_results` method.

.. code-block:: python

   ean.print_results()

Further descriptions of which tables are printed and how to select what is
printed can be found in the :ref:`TESPy analysis section <tespy_advanced_exergy_label>`.
There you can also find more detailed descriptions of how to access the
underlying data for the tabular printouts, which are stored in
`pandas DataFrames <https://pandas.pydata.org/pandas-docs/stable/user_guide/dsintro.html>`_.

With the `plotly <https://plotly.com/>`_ library installed, the results can
also be displayed in a `sankey diagram <https://plotly.com/python/sankey-diagram/>`_.
The :py:meth:`tespy.tools.analyses.ExergyAnalysis.generate_plotly_sankey_input`
method returns a dictionary containing links and nodes for the sankey diagram.

.. code-block:: python

   links, nodes = ean.generate_plotly_sankey_input()
    fig = go.Figure(go.Sankey(
        arrangement="snap",
        node={
            "label": nodes,
            'pad': 11,
            'color': 'orange'},
        link=links))
    plot(fig, filename='NH3_sankey')


.. figure:: /_static/images/tutorials/heat_pump_exergy/NH3_sankey.svg
    :align: center
    :alt: Sankey diagram of the Ground-Coupled Heat Pump (GCHP)

    Figure: Sankey diagram of the GCHP (open in
    new tab to enlarge).

In the figure above you can see the sankey diagram which is created by running
the script of the GCHP with NH3 as refrigerant. Information about, for example,
the colors used or the node order can be found in the
:ref:`TESPy analysis section <tespy_advanced_exergy_label>`.

Post-Processing
^^^^^^^^^^^^^^^
Below, different possibilities of post-processing and visualization of the
exergy analysis results will be presented. The following issues will be
considered:

- plot exergy destruction
- varying ambient and geothermal temperature
- varying geothermal and heating system temperature
- varying heating load and geothermal temperature

In order to be able to compare the results of the two refrigerants NH3 and
R410A, plots of the results of the mentioned issues are created in a separate
plot script :download:`plots.py </../tutorial/plots.py>`. The plots in this
tutorial are created with `Matplotlib <https://matplotlib.org/>`_. For
installation instructions or further documentation please see the Matplotlib
documentation.

For the post-processing, the following additional packages
are required:

.. code-block:: python

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

Plot exergy destruction
+++++++++++++++++++++++
In order to visualize how much exergy of the fuel exergy :code:`E_F` the
individual components of the GCHP destroy, the exergy destruction :code:`E_D`
can be displayed in a bar chart as shown at the end of this section.

To create this diagram, the required data for the diagram must first be
handled. As shown below, the three lists :code:`comps`, :code:`E_D` and
:code:`E_P` are created and first filled with the values for the top bar. A
loop is then used to add all component labels to the list :code:`comps` that
destroy a noticeable amount of exergy (> 1W).  The list :code:`E_D` contains
the corresponding values of the destroyed exergy. List :code:`E_P`, in turn,
contains the value of the exergy that remains after subtracting the destroyed
exergy from the fuel exergy.

.. code-block:: python

    comps = ['E_F']
    E_F = ean.network_data.E_F
    E_D = [0]
    E_P = [E_F]
    for comp in ean.component_data.index:
        # only plot components with exergy destruction > 1 W
        if ean.component_data.E_D[comp] > 1 :
            comps.append(comp)
            E_D.append(ean.component_data.E_D[comp])
            E_F = E_F-ean.component_data.E_D[comp]
            E_P.append(E_F)
    comps.append("E_P")
    E_D.append(0)
    E_P.append(E_F)

With regard to the bar chart to be created, the filled lists are then saved in
a panda DataFrame and exported to a :code:`.csv` file. Exporting the data is
necessary in order to be able to use the results of the two scripts of the
different refrigerants NH3 and R410A in a separate script.

.. code-block:: python

    df_comps = pd.DataFrame(columns= comps)
    df_comps.loc["E_D"] = E_D
    df_comps.loc["E_P"] = E_P
    df_comps.to_csv('NH3_E_D.csv')

.. note::

    In order to be able to use the data from the data frames in a separate
    script for plot creation, all data frames must be saved as a file with
    their own individual name.

In the separate plot script (:download:`plots.py </../tutorial/plots.py>`) the :code:`.csv` files can
now be re-imported to create plots with Matplotlib. The Python code for
creating the bar chart is included in the previously referenced plot script
and can be found there. For more information on creating plots with
Matplotlib, please check the
`Matplotlib documentation <https://matplotlib.org/>`_. The resulting bar chart
is shown below.

.. figure:: /_static/images/tutorials/heat_pump_exergy/diagram_E_D.svg
    :align: center
    :alt: Comparison of exergy destruction and exergy efficiency

    Figure: Comparison of exergy destruction and exergy efficiency of both
    working fluids in design case.

The bar chart shows how much exergy the individual components of the GCHP
destroy in absolute terms and as a percentage of the fuel exergy :code:`E_F`.
After deducting the destroyed exergy :code:`E_D`, the product exergy
:code:`E_P` remains. Overall, it is noticeable that the GCHP with NH3 requires
less fuel exergy than the GCHP with R410A, with the same amount of product
exergy. Furthermore, with NH3 the condenser has the highest exergy destruction,
whereas with R410A the valve destroys the largest amount of exergy.

Varying ambient and geothermal temperature
++++++++++++++++++++++++++++++++++++++++++
In order to consider the influence of a change in ambient temperature or
geothermal temperature on the exergetic efficiency, offdesign calculations are
performed with different values of these parameters. The first step is to
create dataframes as shown below. The ambient temperature :code:`Tamb`
is varied between 1°C and 20°C. The mean geothermal temperature :code:`Tgeo`
is varied between 11.5°C and 6.5°C. Note that the geothermal temperature
:code:`Tgeo` is given as a mean value of the feed an return flow temperatures,
as described in the beginning of this tutorial.

.. code-block:: python

    Tamb_design = Tamb
    Tgeo_design = Tgeo
    i = 0

    # create data ranges and frames
    Tamb_range = np.array([1,4,8,12,16,20])
    Tgeo_range = np.array([11.5, 10.5, 9.5, 8.5, 7.5, 6.5])
    df_eps_Tamb = pd.DataFrame(columns= Tamb_range)
    df_eps_Tgeo = pd.DataFrame(columns= Tgeo_range)

Next, the exergetic efficiency epsilon can be calculated for the different
values of :code:`Tamb` in :code:`Tamb_range` by calling the
:py:meth:`tespy.tools.analyses.ExergyAnalysis.analyse` method in a loop. The
results are saved in the created dataframe and exported to a .csv file.

.. code-block:: python

    # calculate epsilon depending on Tamb
    eps_Tamb = []
    print("Varying ambient temperature:\n")
    for Tamb in Tamb_range:
        i += 1
        ean.analyse(pamb, Tamb)
        eps_Tamb.append(ean.network_data.epsilon)
        print("Case %d: Tamb = %.1f °C"%(i,Tamb))

    # save to data frame
    df_eps_Tamb.loc[Tgeo_design] = eps_Tamb
    df_eps_Tamb.to_csv('NH3_eps_Tamb.csv')

.. note::

    If only the ambient state (temperature or pressure) changes, there is no
    need to create a new :code:`ExergyAnalysis` instance. Instead, you can
    simply call the :py:meth:`tespy.tools.analyses.ExergyAnalysis.analyse`
    method with the new ambient state. A new instance only needs to be created
    when there are changes in the topology of the network.

The following calculation of the network with different geothermal mean
temperatures is carried out as an offdesign calculation. Again, no new
:code:`ExergyAnalysis` instance needs to be created. The ambient temperature
:code:`Tamb` is reset to the design value.

.. code-block:: python

    # calculate epsilon depending on Tgeo
    eps_Tgeo = []
    print("\nVarying mean geothermal temperature:\n")
    for Tgeo in Tgeo_range:
        i += 1
        # set feed and return flow temperatures around mean value Tgeo
        gh_in_ghp.set_attr(T=Tgeo+1.5)
        ev_gh_out.set_attr(T=Tgeo-1.5)
        nw.solve('offdesign', init_path=path, design_path=path)
        ean.analyse(pamb, Tamb_design)
        eps_Tgeo.append(ean.network_data.epsilon)
        print("Case %d: Tgeo = %.1f °C"%(i,Tgeo))

    # save to data frame
    df_eps_Tgeo.loc[Tamb_design] = eps_Tgeo
    df_eps_Tgeo.to_csv('NH3_eps_Tgeo.csv')

The results of the calculation can be plotted as shown in the following
figure. The related Python code to create this plot can be found in the plot
script (:download:`plots.py </../tutorial/plots.py>`). For further documentation
please see the `Matplotlib documentation <https://matplotlib.org/>`_.

.. figure:: /_static/images/tutorials/heat_pump_exergy/diagram_eps_Tamb_Tgeo.svg
    :align: center
    :alt: Varying Tamb and Tgeo of the GCHP

    Figure: Varying ambient and geothermal temperature.

It can be recognized that the specified ambient temperature :code:`Tamb` used
in the :code:`analyse` method of the :code:`ExergyAnalysis` instance has a
considerable influence on the exergetic efficiency epsilon. The closer the
ambient temperature is to the temperature of the heating system, the lower the
exergetic efficiency. This can be argued from the fact that while :code:`E_F`
and :code:`E_P` both decrease with increasing :code:`Tamb`, :code:`E_P`
decreases proportionally more than :code:`E_F`. In comparison, it can be seen
on the right that with increasing :code:`Tgeo`, and thus decreasing
temperature difference between geothermal heat collector and heating system,
epsilon increases. This can be explained by the resulting decrease in
:code:`E_F` with :code:`E_P` remaining constant.

Varying geothermal and heating system temperature
+++++++++++++++++++++++++++++++++++++++++++++++++
Another relation that can be investigated is the influence of a change in the
geothermal and the heating system temperatures on the exergetic efficiency and
the COP of the GCHP. Again, the first step is to create data frames. In this
calculation :code:`Tgeo` is varied between 10.5°C and 6.5°C. The heating
system temperature :code:`Ths` is varied between 42.5°C and 32.5°C. As before,
all temperature values are mean values of the feed and return flow
temperatures.

.. code-block:: python

    # create data ranges and frames
    Tgeo_range = [10.5, 8.5, 6.5]
    Ths_range = [42.5, 37.5, 32.5]
    df_eps_Tgeo_Ths = pd.DataFrame(columns= Ths_range)
    df_cop_Tgeo_Ths = pd.DataFrame(columns= Ths_range)

The values of :code:`Tgeo` and :code:`Ths` are varied simultaneously within
the specified range and again the exergetic efficiency is calculated. In
addition, the COP is calculated for each parameter combination. The data is
stored in two dataframes with the range of :code:`Tgeo` as rows and the range
of :code:`Ths` as columns.

.. code-block:: python

    # calculate epsilon and COP
    print("\nVarying mean geothermal temperature and "+
          "heating system temperature:\n")
    for Tgeo in Tgeo_range:
        # set feed and return flow temperatures around mean value Tgeo
        gh_in_ghp.set_attr(T=Tgeo+1.5)
        ev_gh_out.set_attr(T=Tgeo-1.5)
        epsilon = []
        cop = []
        for Ths in Ths_range:
            i += 1
            cd_hs_feed.set_attr(T=Ths+2.5)
            hs_ret_hsp.set_attr(T=Ths-2.5)
            if Ths == Ths_range[0]:
                nw.solve('offdesign', init_path=path, design_path=path)
            else:
                nw.solve('offdesign', design_path=path)
            ean.analyse(pamb, Tamb_design)
            epsilon.append(ean.network_data.epsilon)
            cop += [abs(cd.Q.val) / (cp.P.val + ghp.P.val + hsp.P.val)]
            print("Case %d: Tgeo = %.1f °C, Ths = %.1f °C"%(i,Tgeo,Ths))

        # save to data frame
        df_eps_Tgeo_Ths.loc[Tgeo] = epsilon
        df_cop_Tgeo_Ths.loc[Tgeo] = cop

    df_eps_Tgeo_Ths.to_csv('NH3_eps_Tgeo_Ths.csv')
    df_cop_Tgeo_Ths.to_csv('NH3_cop_Tgeo_Ths.csv')


The results of this calculation are shown in the following figure. The
corresponding Python code can likewise be found in the plot script
(:download:`plots.py </../tutorial/plots.py>`).

.. figure:: /_static/images/tutorials/heat_pump_exergy/diagram_cop_eps_Tgeo_Ths.svg
    :align: center
    :alt: Varying Tgeo and Ths of the GCHP

    Figure: Varying geothermal and heating system temperature.

It can be seen that the GCHP with NH3 has a better exergetic efficiency than
with R410A. As in the prior investigation, an increasing geothermal heat
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
The rated load was previously set at 4 kW in the design calculation. Due to the
similarity to the previous parameter variation, the corresponding Python code
is not presented, but can be found in the scripts linked at the beginning
instead.

.. figure:: /_static/images/tutorials/heat_pump_exergy/diagram_cop_eps_Tgeo_Q.svg
    :align: center
    :alt: Varying Tgeo and Q of the GCHP

    Figure: Varying geothermal temperature and heat load.

The results are shown in the figure above. As before, the Python code for
creating the plot can be found in the plot script
(:download:`plots.py </../tutorial/plots.py>`).
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
:ref:`TESPy analysis section <tespy_advanced_exergy_label>` and in the API
documentation of the :py:class:`tespy.tools.analyses.ExergyAnalysis` class. If
you are interested in contributing or have questions and remarks on this
tutorial, you are welcome to file an issue at our GitHub page.
