.. _plotting_label:

Plotting
========

The :py:mod:`tespy.plotting` module draws property diagrams, Q-T diagrams,
characteristics and results of your models. Every kind of plot is available
with matplotlib and with Plotly:

- :code:`plot_*` functions draw with matplotlib and return the axes,
- :code:`plot_*_plotly` functions draw with Plotly and return the figure, in
  which hovering a point shows the values behind it,
- :code:`*_data` functions return the data of the plot without drawing.

Property diagrams are drawn by `CoolPlot <https://github.com/CoolProp/CoolPlot>`__.
Install the plotting dependencies with

.. code-block:: bash

    pip install tespy[plot]         # matplotlib and CoolPlot
    pip install tespy[interactive]  # additionally Plotly

All examples on this page use the following heat pump.

.. code-block:: python

    >>> from tespy.networks import Network
    >>> from tespy.connections import Connection
    >>> from tespy.components import (
    ...     CycleCloser, MovingBoundaryHeatExchanger, Compressor, Valve,
    ...     SimpleHeatExchanger, Source, Sink
    ... )

    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(
    ...     temperature="°C", pressure="bar", pressure_difference="bar",
    ...     heat="kW", power="kW", enthalpy="kJ/kg", entropy="kJ/kg/K"
    ... )

    >>> cp = Compressor("compressor")
    >>> cc = CycleCloser("cycle_closer")
    >>> cd = MovingBoundaryHeatExchanger("condenser")
    >>> va = Valve("expansion valve")
    >>> ev = SimpleHeatExchanger("evaporator")
    >>> so = Source("water source")
    >>> si = Sink("water sink")

    >>> c1 = Connection(cc, "out1", cd, "in1", label="c1")
    >>> c2 = Connection(cd, "out1", va, "in1", label="c2")
    >>> c3 = Connection(va, "out1", ev, "in1", label="c3")
    >>> c4 = Connection(ev, "out1", cp, "in1", label="c4")
    >>> c5 = Connection(cp, "out1", cc, "in1", label="c5")
    >>> a1 = Connection(so, "out1", cd, "in2", label="a1")
    >>> a2 = Connection(cd, "out2", si, "in1", label="a2")
    >>> nw.add_conns(c1, c2, c3, c4, c5, a1, a2)

    >>> cd.set_attr(dp1=0, dp2=0, Q=-1e3, td_pinch=5)
    >>> ev.set_attr(dp=0)
    >>> cp.set_attr(eta_s=0.8)
    >>> c1.set_attr(fluid={"R290": 1})
    >>> c2.set_attr(td_bubble=5)
    >>> c4.set_attr(td_dew=5, T=15)
    >>> a1.set_attr(fluid={"water": 1}, p=1, T=50)
    >>> a2.set_attr(T=65)
    >>> nw.solve("design")
    >>> nw.assert_convergence()

Property diagrams
-----------------

Pass the network and the label of any connection of the part of the network
you want to see. All connections with the same fluid that are linked by
components are drawn. The diagram type names the vertical axis first:
:code:`"ph"`, :code:`"Ts"`, :code:`"hs"`, :code:`"ps"`, :code:`"prho"`,
:code:`"Trho"` and :code:`"pT"`. The axes use the units of the network.

.. code-block:: python

    >>> from tespy.plotting import plot_property_diagram
    >>> ax = plot_property_diagram(nw, "c1")
    >>> ax = plot_property_diagram(nw, "c1", "Ts", title="R290 heat pump")

.. figure:: /_static/images/modules/plotting/logph_light.png
    :align: center
    :alt: log(p)-h diagram of a R290 heat pump
    :figclass: only-light

    Figure: log(p)-h diagram of a R290 heat pump

.. figure:: /_static/images/modules/plotting/logph_dark.png
    :align: center
    :alt: log(p)-h diagram of a R290 heat pump
    :figclass: only-dark

    Figure: log(p)-h diagram of a R290 heat pump

The state points are the results of TESPy. CoolPlot draws the saturation
dome, the isolines and the paths between the states. To change the diagram
before it is drawn, e.g. the isolines or the axis limits, pass a function
that receives the :code:`coolplot.PropertyDiagram`:

.. code-block:: python

    >>> def more_isotherms(diagram):
    ...     _ = diagram.set_lines("T", values=list(range(-20, 121, 10)), labels=True)
    >>> ax = plot_property_diagram(nw, "c1", diagram_hook=more_isotherms)

Property diagrams are available for fluids of the CoolProp back end, pure
fluids, predefined blends and CoolProp mixtures with a mass or molar based
composition. Other fluid property back ends and mixtures of several fluids
raise an error.

Q-T diagrams
------------

The Q-T diagram of a heat exchanger shows the temperatures of both sides over
the heat transferred, with the section and phase boundaries as markers and
the pinch. A list of heat exchangers is drawn next to each other, e.g. the
economizer, evaporator and superheater of a steam generator.

.. code-block:: python

    >>> from tespy.plotting import plot_QT_diagram
    >>> ax = plot_QT_diagram(cd)

.. figure:: /_static/images/modules/plotting/QT_light.png
    :align: center
    :alt: Q-T diagram of a condenser
    :figclass: only-light

    Figure: Q-T diagram of a condenser

.. figure:: /_static/images/modules/plotting/QT_dark.png
    :align: center
    :alt: Q-T diagram of a condenser
    :figclass: only-dark

    Figure: Q-T diagram of a condenser

Characteristics
---------------

Characteristic lines, maps and polynomials are drawn with their valid range.
Operating points outside the range of the data are marked, an extrapolated
part of a line is dashed. For a component of a solved network, pass the
component and the name of the characteristic to mark its operating point,
e.g. :code:`plot_characteristic(compressor, "eta_s_char")`. Polynomials have
no range of their own, give it with :code:`ranges`.

.. code-block:: python

    >>> from tespy.plotting import plot_characteristic
    >>> from tespy.tools.characteristics import CharMap, EN12900Polynomial
    >>> compressor_map = CharMap(
    ...     x=[0.9, 1.0, 1.1],
    ...     y=[[0.5, 0.8, 1.0], [0.6, 0.9, 1.1], [0.7, 1.0, 1.2]],
    ...     z=[[0.8, 0.84, 0.85], [0.82, 0.87, 0.88], [0.84, 0.89, 0.9]],
    ... )
    >>> ax = plot_characteristic(compressor_map, operating_points=[(0.95, 0.85)])
    >>> capacity = EN12900Polynomial(
    ...     [4000, 120, -40, 1.5, -0.6, 0.2, 0, 0, 0, 0], value_unit="W"
    ... )
    >>> ax = plot_characteristic(capacity, ranges=[(-25, 10), (25, 60)])

.. figure:: /_static/images/modules/plotting/char_map_light.png
    :align: center
    :alt: Characteristic map with an operating point
    :figclass: only-light

    Figure: Characteristic map with an operating point

.. figure:: /_static/images/modules/plotting/char_map_dark.png
    :align: center
    :alt: Characteristic map with an operating point
    :figclass: only-dark

    Figure: Characteristic map with an operating point

Results
-------

Further plots show the results of a parameter sweep and results on a grid
of two parameters. Both take a :code:`pandas.DataFrame` or a dictionary of
lists, one row per solve.

.. code-block:: python

    >>> from tespy.plotting import plot_operating_map, plot_sweep
    >>> results = {
    ...     "T_cond": [35, 45, 55, 35, 45, 55],
    ...     "T_evap": [0, 0, 0, 10, 10, 10],
    ...     "COP": [5.0, 4.0, 3.3, 6.2, 4.8, 3.9],
    ... }
    >>> ax = plot_sweep(results, "T_cond", "COP", group="T_evap")
    >>> ax = plot_operating_map(results, "T_cond", "T_evap", "COP")

.. figure:: /_static/images/modules/plotting/sweep_light.png
    :align: center
    :alt: COP of two working fluids over the condensation temperature
    :figclass: only-light

    Figure: Parameter sweep

.. figure:: /_static/images/modules/plotting/sweep_dark.png
    :align: center
    :alt: COP of two working fluids over the condensation temperature
    :figclass: only-dark

    Figure: Parameter sweep

Structure of the system of equations
------------------------------------

The incidence matrix shows which equation depends on which variable, for
the equations and variables that remain after the presolving of the network.
With :code:`block_order=True` equations and variables are sorted by the
blocks the solver works through: the matrix is then block lower triangular,
every block only depends on its own variables and on those of blocks solved
before it. Blocks of several equations that are solved together are
highlighted and outlined, dependencies on earlier blocks are lighter. It is the plot
counterpart of :code:`nw.print_incidence_matrix`.

.. code-block:: python

    >>> from tespy.plotting import plot_incidence_matrix
    >>> ax = plot_incidence_matrix(nw)
    >>> ax = plot_incidence_matrix(nw, block_order=True)

.. figure:: /_static/images/modules/plotting/block_form_light.png
    :align: center
    :alt: Block lower triangular form of the incidence matrix
    :figclass: only-light

    Figure: Block lower triangular form of the incidence matrix

.. figure:: /_static/images/modules/plotting/block_form_dark.png
    :align: center
    :alt: Block lower triangular form of the incidence matrix
    :figclass: only-dark

    Figure: Block lower triangular form of the incidence matrix

Interactive plots in the browser
--------------------------------

The :code:`_plotly` functions take the same arguments and return a Plotly
figure. Hovering a point shows its values with units, e.g. temperature,
pressure, enthalpy, entropy and mass flow of a state, or the temperature
difference along a heat exchanger.

.. code-block:: python

    >>> from tespy.plotting import plot_property_diagram_plotly, plot_QT_diagram_plotly
    >>> figure = plot_property_diagram_plotly(nw, "c1")
    >>> figure = plot_QT_diagram_plotly(cd)

Save the figure as a web page or show it in a notebook:

.. code-block:: python

    figure.write_html("heat_pump.html", include_plotlyjs="cdn")
    figure.show()

Every point knows the network object it belongs to. In the browser, the last
entry of the :code:`customdata` of a point is a reference such as
:code:`{"connection": "c1"}` or :code:`{"component": "condenser", "side":
"hot"}`, and :code:`figure.layout.meta["tespy"]` holds the data of the whole
plot.

Plot data
---------

Every plot is created from plain data, which you can retrieve without
drawing and export to JSON, e.g. to process it in another program. The
format is provisional and carries a version number.

.. code-block:: python

    >>> from tespy.plotting import property_diagram_data
    >>> description = property_diagram_data(nw, "c1")
    >>> description.x_axis.title
    'specific enthalpy in kJ/kg'
    >>> states = description.series[0]
    >>> states.references[0]
    {'connection': 'c1'}
    >>> text = description.to_json()

An exported description can be drawn again with
:py:func:`tespy.plotting.draw_matplotlib` or
:py:func:`tespy.plotting.draw_plotly`.

Styles
------

TESPy ships three styles:

- :code:`"tespy"` - the default for light backgrounds,
- :code:`"tespy-dark"` - the same colors stepped for dark backgrounds,
- :code:`"tespy-paper"` - black, gray and one blue, where dashes, markers
  and hatches tell series apart, e.g. for printing.

.. figure:: /_static/images/modules/plotting/logph_paper.png
    :align: center
    :alt: log(p)-h diagram in the paper style

    Figure: log(p)-h diagram in the :code:`"tespy-paper"` style

A style holds colors by their role: surface, text, hot and cold side, process
lines, categories, colormaps and isolines. Choose a style per plot or for the
session:

.. code-block:: python

    >>> from tespy.plotting import set_defaults
    >>> ax = plot_QT_diagram(cd, style="tespy-dark")
    >>> set_defaults(style="tespy-paper")
    >>> set_defaults(style="tespy")

To use your own colors, derive a style from a shipped one and register it:

.. code-block:: python

    >>> from tespy.plotting import get_style, register_style
    >>> company = get_style("tespy").derive(
    ...     "company",
    ...     categorical=("#00395B", "#EC6707", "#74ADC0", "#B54036"),
    ...     cycle="#00395B", hot="#B54036", cold="#00395B",
    ... )
    >>> register_style(company, overwrite=True)
    >>> ax = plot_QT_diagram(cd, style="company")

A style can also be read from a TOML file. Only the roles that differ from
the style named in :code:`based_on` are needed:

.. code-block:: toml

    name = "company"
    based_on = "tespy"
    categorical = ["#00395B", "#EC6707", "#74ADC0", "#B54036"]
    cycle = "#00395B"

.. code-block:: python

    from tespy.plotting import Style, register_style
    register_style(Style.from_toml("company.toml"))

Logo and information block
--------------------------

Every plot shows the TESPy logo and an information block with the TESPy
version. Property diagrams add the reference of the equation of state and
the versions of CoolProp and CoolPlot. In matplotlib plots the legend is
part of the same block, which is placed in the corner with the least data.
Add your own lines with :code:`info`, fix the corner with
:code:`logo_location` or switch logo and information off with
:code:`logo=False`, per plot or for the session:

.. code-block:: python

    >>> ax = plot_QT_diagram(cd, info=["R290 heat pump"], logo_location="upper left")
    >>> set_defaults(logo=False)
    >>> set_defaults(logo=True)

Model templates
---------------

Models based on :py:class:`tespy.models.template.ModelTemplate` have the
methods :code:`plot_property_diagram`, :code:`plot_QT_diagram` and their
:code:`_plotly` variants, which take connection and component labels.

Migration from fluprodia
------------------------

The plotting methods of the model template based on fluprodia are
deprecated and will be removed in version 0.14.

.. list-table::
    :header-rows: 1

    * - Deprecated
      - New
    * - :code:`plot_logph_diagram_matplotlib(label)`
      - :code:`plot_property_diagram(label, "ph")`
    * - :code:`plot_Ts_diagram_matplotlib(label)`
      - :code:`plot_property_diagram(label, "Ts")`
    * - :code:`plot_QT_diagram_matplotlib(label)`
      - :code:`plot_QT_diagram(label)`
    * - :code:`register_diagram(fluid, diagram)`
      - :code:`diagram_hook` argument of :code:`plot_property_diagram`
    * - :code:`CharLine.plot(path, ...)`, :code:`CharMap.plot(path, ...)`
      - :py:func:`tespy.plotting.plot_characteristic`

.. code-block:: python

    >>> from matplotlib import pyplot as plt
    >>> plt.close("all")
