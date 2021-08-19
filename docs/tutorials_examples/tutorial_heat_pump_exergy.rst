Exergy Analysis of a Ground-Coupled Heat Pump
---------------------------------------------

.. contents::
    :depth: 1
    :local:
    :backlinks: top
     
          
Task
^^^^

This tutorial shows how to set up and carry out an exergy analysis for a 
ground-coupled heat pump (GCHP). In addition, various post-processing options 
are presented. To investigate the impact of refrigerant choice on COP and 
exergetic efficiency, two Python scripts of the same network with different 
refrigerants (NH3 and R410A) are created. Finally, the influence of varying 
different parameters on COP and exergetic efficiency is investigated and 
plotted. 

Since there is an existing tutorial for creating a heat pump (Referenz Heat Pump Tutorial), 
this tutorial starts with the explanations for setting up the exergy analysis. 
Note, however, that the heat pump model differs slightly in structure from the
model in the previous tutorial. All related Python scripts of the fully 
working GCHP-model are listed in the following:

- GCHP with NH3 (the model only): 
    :download:`NH3 <../tutorial/NH3.py>`
- GCHP with R410A (the model only): 
    :download:`R410A <../tutorial/R410A.py>`
- GCHP with NH3 (model and post-processing): 
    :download:`NH3_calculations <../tutorial/NH3_calculations.py>`
- GCHP with R410A (model and post-processing): 
    :download:`R410A_calculations <../tutorial/R410A_calculations.py>`
- Plots of the results of the parameter variations: 
    :download:`plots <../tutorial/plots.py>`


.. figure:: api/_images/heat_pump_exergy_flowsheet.svg
    :align: center
    :alt: Topology of the Ground-Coupled Heat Pump (GCHP)
    
    
The figure above shows the topology of the GCHP. In this model, a 
ground-coupled heat pump is modeled, which is for instance connected to a 
single-family house with underfloor heating. The heating system represents the 
consumed heat of the house. The geothermal heat collector is represented by a 
ground heat feed flow (source) and return flow (sink). The heat pump circuit 
consists of the basic components: condenser, expansion valve, evaporator and 
compressor.  

The input data of the model are based on different literature. In general, 
the model of the GCHP is based on a data sheet of a real heat pump 
(`Viessmann Vitocal 300-G <https://www.viessmann.de/de/wohngebaeude/waermepumpe/sole-wasser-waermepumpen/vitocal-300-g.html>`_ ). 
However, the data are used as approximate values to create a model 
that works with both NH3 and R410A, although the mentioned heat pump is 
designed to use R410A. 
The range of the underfloor heating system temperature and the range of the 
geothermal temperature are assumptions based on measured data from the 
research project 
`WPsmart <https://wp-monitoring.ise.fraunhofer.de/wp-smart-im-bestand/german/index/index.html>`_ 
and :cite:`Chen2015`. The average outdoor temperature is 
taken from :cite:`Chen2015`.



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

The model using NH3 as refrigerant (Referenz NH3 Modell) and the model using 
R410A as refrigerant (Referenz R410A Modell) differ in the fluid definition 
and the specification of the starting values only. The definition of the 
starting values is necessary to obtain a numerical solution for the first 
calculation without an inital path. In this tutorial, the given code examples 
are only given for the example with NH3 as refrigerant.

The following parameter specifications were made for the design case 
calculation:

- isentropic efficiency values 
- electrical conversion efficiencies of compressor and pumps
- terminal temperature difference values at condenser and evaporator
- pressure losses in condenser and evaporator 
- hot and cold side heat transfer coefficents of evaporator
- temperature difference to boiling point of refrigerant at compressor inlet
- temperatures and pressure of heating sytem feed and return flow
- temperatures and pressure of geothermal heat collector feed and return flow 
- condenser heat output
 
The units used and the ambient state are defined as follows:

.. code-block:: python

    nw = Network(fluids=['water', 'NH3'], T_unit='C', p_unit='bar',
             h_unit='kJ / kg', m_unit='kg / s')

    pamb = 1.013
    Tamb = 2.8

The temperature of the heating system feed flow is set to 40°C in design 
calculation. The difference between feed and return flow temperature is kept 
constant at 5°C. It follows, that the return flow is set to 35°C. 

The geothermal heat collector temperature is defined as follows:

.. code-block:: python

    Tgeo = 9.5
    
:code:`Tgeo` is the mean geothermal temperature. The difference between 
feed and return flow temperature is kept constant at 3°C. Therefore, the feed 
flow temperature in the design calculation is set to :code:`Tgeo + 1.5°C` and 
the return flow temperature is set to :code:`Tgeo - 1.5°C`. 

All other specified values of the component and connection parameters can be 
found in the Python scripts referenced above.


Exergy analysis
^^^^^^^^^^^^^^^

Analysis setup
++++++++++++++

After the network has been built, the exergy analysis can be set up. For this, 
all exergy flows entering and leaving the network must be defined. The exergy 
flows are defined as a list of busses as follows: 
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
streams :code:`hs_ret` and :code:`hs_feed`. Note, that the :code:`base` 
keyword of the stream entering the network :code:`hs_ret` must be set to 
:code:`bus`. 

Lastly, the geothermal heat bus represents the heat that is transferred from 
the geothermal heat collector to the evaporator. The bus consists of the 
streams :code:`gh_in` and :code:`gh_out`. Here, the :code:`base` of the stream 
:code:`gh_in` is set to :code:`bus`, because this stream represents an energy 
input from outside of the network. 
In this example, the geothermal heat bus is defined as fuel exergy, because 
the ambient temperature :code:`Tamb` is set at a lower temperature 
than the temperature of the geothermal heat collector. 

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
Other examples of exergy analysis setup can be found (Referenz Tutroial Exergy) 
and in the API documentation of class 
:py:class:`tespy.tools.analyses.ExergyAnalysis`.  

.. code-block:: python

   ean = ExergyAnalysis(network=nw, 
                         E_F=[power, heat_geo], 
                         E_P=[heat_cons])
                                      
   ean.analyse(pamb, Tamb)

The :py:meth:`tespy.tools.analyses.ExergyAnalysis.analyse` method will run the 
exergy analysis automatically. This method expects information about the 
ambient pressure and ambient temperature. Additionally, an automatic check of
consistency is performed by the analysis as further described in 
(Referenz Tutorial Exergy). 


Results
+++++++

The results can be printed by using the 
:py:meth:`tespy.tools.analyses.ExergyAnalysis.print_results` method.

.. code-block:: python

   ean.print_results()

Further descriptions of which tables are printed and how to select what is 
printed can be found in the tutorial (Referenz Tutorial Exergie).
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
    

.. figure:: api/_images/NH3_sankey.html
    :align: center
    :alt: Sankey diagram of the Ground-Coupled Heat Pump (GCHP)

In the figure above you can see the sankey diagram which is created by running 
the script of the GCHP with NH3 as refrigerant (Skript Referenz einfügen). 
Information about, for example, the colors used or the node order can be found 
in the tutorial (Referenz Tutorial Exergie). 

The full Python code up to this step of the tutorial is available in the 
scripts (Referenz Skript NH3 Modell) with NH3 as refrigerant and
(Referenz Skript R410A Modell) with R410A as refrigerant. 
 

Post-Processing
^^^^^^^^^^^^^^^ 

Below, different possibilities of post-processing and visualization of
the exergy analysis results will be presented. The following issues will be 
considered: 
- create an h-log(p) diagram
- plot exergy destruction
- varying ambient and geothermal temperature
- varying geothermal and heating system temperature
- varying heating load and geothermal temperature

In the scripts (Referenz Skript NH3 Berechnungen) and (Referenz Skript R410A Berechnungen), 
the Phython code of the first steps of this tutorial is extended to handle the 
listed post-processing issues. Since these scripts differ almost only in the 
definition of the fluid, the specification of the starting values and the 
naming of the stored files, the lines of code from the scripts listed below 
are as before only shown using NH3 as an example. If the scripts differ beyond 
the mentioned points, it will be pointed out at the respective place of the 
tutorial. 

In addition, script (Referenz Skript Plots) includes the python code for 
creating the plots of the last three issues. The plots in this tutorial are 
created with `Matplotlib <https://matplotlib.org/>`_. For installation 
instructions or further documentation please see the Matplotlib documentation.  

For the post-processing issues considered, the following additional packages 
are required:

.. code-block:: python

    from fluprodia import FluidPropertyDiagram
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt


h-log(p)-diagram
++++++++++++++++
 


Plot exergy destruction
+++++++++++++++++++++++

By adding the code below, the exergy destruction of the components is 
displayed in a block diagram, as shown in the two figures below. Only exergy 
destruction of components higher than 1 W will be displayed.

.. code-block:: python

    comps = []
    x= []
    for comp in ean.component_data.index:
        # only plot components with exergy destruction > 1 W
        if ean.component_data.E_D[comp] > 1 : 
            comps.append(comp)
            x.append(ean.component_data.E_D[comp])
    y = (comps)
    y_pos = np.arange(len(comps))
         
    fig, ax = plt.subplots()
    hbars = ax.barh(y_pos, x, align='center')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(y)
    ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_xlabel('E_D')
    ax.set_title('Component Exergy Destruction "NH3"')
    ax.set_xlim(right=200)  # adjust xlim to fit labels
    
    plt.show()
    fig.savefig('NH3_E_D.svg', bbox_inches='tight')

.. figure:: api/_images/NH3_E_D.svg
    :align: center
    :alt: Exergy Destruction of the GCHP - NH3
    
.. figure:: api/_images/R410A_E_D.svg
    :align: center
    :alt: Exergy Destruction of the GCHP - R410A


Varying ambient and geothermal temperature
++++++++++++++++++++++++++++++++++++++++++

In order to consider the influence of a change in ambient temperature or 
geothermal temperature on the exergetic efficiency, offdesign calculations are 
performed with different values of these parameters. The first step is to 
create data frames as shown below. The ambient temperature :code:`Tamb` 
is varied between 1°C and 20°C. The mean geothermal temperature :conde:`Tgeo`
is varied between 11.5°C and 6.5°C. 
Note that the geothermal temperature :code:`Tgeo` is given as a mean value of 
the feed an return flow temperatures, as described in the beginning of this 
tutorial. 

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
results are saved in the created data frame. 

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

.. note::

    In order to be able to use the data from the data frames in a separate 
    script for plot creation, all data frames must be saved as a file with 
    their own individual name. 

The results of the calculation can be plotted as shown in the following 
figure. The related python code to create this plot can be found in the script 
(Referenz Plot Skript). For further documentation please see the 
`Matplotlib documentation <https://matplotlib.org/>`_. 

.. figure:: api/_images/diagram_eps_Tamb_Tgeo.svg
    :align: center
    :alt: Varying Tamb and Tgeo of the GCHP


Varying geothermal and heating system temperature
+++++++++++++++++++++++++++++++++++++++++++++++++

Another issue that can be investigated is the influence of a change in the 
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
stored in two data frames with the range of :code:`Tgeo` as rows and the range 
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
(Referenz Plot Skript).     
    
.. figure:: api/_images/diagram_cop_eps_Tgeo_Ths.svg
    :align: center
    :alt: Varying Tgeo and Ths of the GCHP


Varying geothermal temperature and heating load
+++++++++++++++++++++++++++++++++++++++++++++++

Finally, the influence of the simultaneous variation of the geothermal 
temperature :code:`Tgeo` and the heating load :code:`Q` on the exergetic 
efficiency and the COP of the GCHP is examined. The investigation is carried 
out in the same way as the variation of :code:`Tgeo` and :code:`Ths` described 
above. In contrast to the previous investigation, :code:`Q` is varied here 
instead of :code:`Ths`. The range of :code:`Q` varies between 4.3 and 2.8 kW.
The rated load was previously set at 4 W in the design calculation. Due to the 
similarity to the previous parameter variation, the corresponding Python code 
is not presented, but can be found in the scripts (Referenz NH3 calculations) 
and (Referenz R410A calculations) instead. 

.. figure:: api/_images/diagram_cop_eps_Tgeo_Q.svg
    :align: center
    :alt: Varying Tgeo and Q of the GCHP

The results are shown in the figure above. As before, the Python code for 
creating the plot can be found in the script (Referenz plot Skript).


Conclusion
^^^^^^^^^^

This tutorial provides an exemplary insight into post-processing with the
TESPy exergy analysis tool. Of course, other parameters can also be examined 
and varied. Feel free to try out different parameter variations. But make sure 
that the data ranges are not only adjusted in the Python script of the model, 
but also in the Python script of the plots, if a plot is created with the 
stand-alone plot script. 

More examples of exergy analysis can be found in the documentation of the 
exergy analysis (Referenz Tutorial Exergie) and in the API documentation of 
the :py:class:`tespy.tools.analyses.ExergyAnalysis` class. If you are 
interested in contributing, you are welcome to file an issue at our GitHub 
page. 




