.. _modules_networks_label:

Networks
========
The network class handles preprocessing, solving and post-processing.
We will walk you through all the important steps.

Setup
-----
Network container
^^^^^^^^^^^^^^^^^
The TESPy network contains all data of your plant, which in terms of the
calculation is represented by a nonlinear system of equations. The system
variables of your TESPy network are a subset of the following:

* mass flow,
* pressure,
* enthalpy and
* the mass fractions of the fluids of every :code:`Connection` instance as well
  as
* energy flow (power or heat) for :code:`PowerConnection` instances and
* component variables (e.g. compressor rpm) if declared a variable.

The solver will simplify the variable space in a presolving step and then solve
for the remaining variables. In the following, you will find information on the
:code:`Network` setup, solving and debugging.

.. _units_label:

Unit specifications
+++++++++++++++++++
There are a couple of ways to impose boundary conditions with the respective
units to your models. The unit can be specified through two different ways:

- as default units per quantity (e.g. temperature, temperature difference,
  power, ...) for all properties of the components and connections, or
- as individual units of a specific parameter.

.. tip::

    TESPy implements :code:`pint` in the back-end to handle the units. If you
    want to learn more on how to work with :code:`pint` units, check out the
    respective `documentation <https://pint.readthedocs.io>`__.

The example below shows, how you can setup units for your :code:`Network`.

.. code-block:: python

    >>> from tespy.networks import Network
    >>> nw = Network()

Default units are SI units. We can check, what unit is the default unit like
this:

.. code-block:: python

    >>> nw.units.get_default("temperature")
    'kelvin'
    >>> nw.units.get_default("power")
    'W'
    >>> nw.units.get_default("specific_volume")
    'm3/kg'

We can change the default units as follows:

.. code-block:: python

    >>> nw.units.set_defaults(temperature="degC")  # celsius
    >>> nw.units.get_default("temperature")
    'degC'
    >>> nw.units.set_defaults(power="hp")  # horse power
    >>> nw.units.get_default("power")
    'hp'
    >>> nw.units.set_defaults(efficiency="%")  # percent
    >>> nw.units.set_defaults(pressure="bar", pressure_difference="bar")

The unit specification then applies to all parameters with the same quantity.
For example, let's set up a model of a compressor.

.. code-block:: python

    >>> from tespy.components import Compressor, Source, Sink
    >>> from tespy.connections import Connection
    >>> source = Source("gas inflow")
    >>> compressor = Compressor("compressor")
    >>> sink = Sink("gas discharge")
    >>> c1 = Connection(source, "out1", compressor, "in1", label="c1")
    >>> c2 = Connection(compressor, "out1", sink, "in1", label="c2")
    >>> nw.add_conns(c1, c2)

Now we can parametrize our problem. It will utilize the unit specifications we
created above.

.. code-block:: python

    >>> c1.set_attr(fluid={"air": 1}, m=1, p=1, T=25)  # p in bar, T in celsius
    >>> c2.set_attr(p=3)
    >>> compressor.set_attr(eta_s=80)  # efficiency in %

We can disable the iteration printouts by setting :code:`iterinfo=False`

    >>> nw.iterinfo = False
    >>> nw.solve("design")

Now we can check results, e.g. the power of the compressor, which is expected
to be in horse power:

.. code-block:: python

    >>> round(compressor.P.val, 0)
    185.0

We can also retrieve the value with the respective unit, and then use pint to
transform it into what ever unit we need:

.. code-block:: python

    >>> round(compressor.P.val_with_unit, 0)
    <Quantity(185.0, 'horsepower')>
    >>> round(compressor.P.val_with_unit.to("kW"), 0)
    <Quantity(138.0, 'kilowatt')>

Alternatively, we can specify an individual unit using the :code:`Quantity`
class of pint. For that you have to utilize the :code:`UnitRegistry` of
your :code:`Network.units`: :code:`ureg`.

.. code-block:: python

    >>> ureg = nw.units.ureg
    >>> c1.set_attr(m=ureg.Quantity(1, "t/h"))
    >>> c1.m.val_with_unit
    <Quantity(1, 'metric_ton / hour')>

.. caution::

    If you now modify this number, it will remove the individual unit! On top,
    by specifying a bare number, you will always remove the unit information.
    This information is only reconnected once, the :code:`Network` is
    initialized in context of a simulation!

    .. code-block:: python

        >>> c1.set_attr(m=5)
        >>> c1.m.val_with_unit
        5

    .. code-block:: python

        >>> nw.solve("design")
        >>> c1.m.val_with_unit
        <Quantity(5, 'kilogram / second')>

To understand, what quantity is associated with a specific parameter, you can
do the following:

.. code-block:: python

    >>> compressor.dp.quantity
    'pressure_difference'
    >>> c1.td_dew.quantity
    'temperature_difference'

It is also possible to use your own :code:`UnitRegistry`:

.. code-block:: python

    >>> from pint import UnitRegistry
    >>> ureg = UnitRegistry()
    >>> nw.units.set_ureg(ureg)

Changing the ureg will only have effect on future specifications. Existing
quantities are not changed.

.. attention::

    The quantities :code:`pressure` and :code:`pressure_difference` as well as
    :code:`temperature` and :code:`temperature_difference` need to be set
    individually!

.. _printout_logging_label:

Printouts and logging
+++++++++++++++++++++
TESPy comes with an inbuilt logger. If you want to keep track of
debugging-messages, general information, warnings or errors you should enable
the logger. At the beginning of your python script e.g. add the following
lines:

.. code-block:: python

    >>> from tespy.tools import logger
    >>> import logging
    >>> ();logger.define_logging(
    ...     logpath="myloggings", log_the_path=True, log_the_version=True,
    ...     screen_level=logging.ERROR, file_level=logging.DEBUG
    ... );()  # +doctest: ELIPSIS
    (...)

The log-file will be saved to :code:`~/.tespy/log_files/` by default. All
available options are documented in the
:py:func:`API <tespy.tools.logger.define_logging>`.

Prior to solving the network there are options regarding the **console**
**printouts for the calculation progress**. Specify, if you want to enable or
disable convergence progress printouts:

.. code-block:: python

    >>> nw.iterinfo
    False

    # enable iteration information printout
    >>> nw.iterinfo = True
    >>> nw.iterinfo
    True

    # disable iteration information printout
    >>> nw.iterinfo = False

Adding connections
++++++++++++++++++
As seen in the introduction, you will have to create your networks from the
components and the connections between them. You can add connections directly
or via subsystems using the corresponding methods:

.. code-block:: python

    >>> nw.add_conns()
    >>> nw.add_subsystems()

.. note::

    You do not need to add the components to the network, as they are inherited
    via the added connections. After having set up your network and added all
    required elements, you can start the calculation.

There are two types of connections, you can learn about them more in
:ref:`these sections <modules_connections_label>`.

Start calculation
^^^^^^^^^^^^^^^^^
You can start the solution process with the following line:

.. code-block:: python

    nw.solve(mode='design')

This starts the initialisation of your network and proceeds to its calculation.
The specification of the **calculation mode is mandatory**, This is the list of
available keywords:

- :code:`mode` is the calculation mode (:code:`'design'`-calculation or
  :code:`'offdesign'`-calculation).
- :code:`init_path` is the path to a saved network state (json file path or
  a :code:`dict` returned by :code:`nw.save(as_dict=True)`) to use for
  initialisation.
- :code:`design_path` is the path to a saved network state (json file path or
  a :code:`dict` returned by :code:`nw.save(as_dict=True)`) which holds the
  information of your plant's design point.
- :code:`max_iter` is the maximum amount of iterations performed by the
  solver.
- :code:`min_iter` is the minimum amount of iterations before a solution can
  be accepted (given the convergence criterion is satisfied).
- :code:`init_only` stop after initialisation (True/False).
- :code:`init_previous` use starting values from previous simulation
  (True/False).
- :code:`use_cuda` use cuda instead of numpy for matrix inversion, speeds up
  simulation in some cases by outsourcing calculation to graphics card. For
  more information please visit the
  `cupy documentation <https://docs.cupy.dev/en/stable/index.html>`_.

There are two calculation modes available (:code:`'design'` and
:code:`'offdesign'`), which are explained in the subsections below. If you
choose :code:`offdesign` as calculation mode the specification of a
:code:`design_path` is mandatory.

The usage of an initialisation path is always optional but highly recommended,
as the convergence of the solution process will be improved, if you provide
good starting values. If you do not specify an :code:`init_path`, the
initialisation from saved results will be skipped.
:code:`init_only=True` usually is used for debugging. Or, you could use this
feature to export a not solved network, if you want to do the parametrisation
in .csv-files rather than your python script.

The :code:`init_previous` parameter can be used in design and offdesign
calculations and works very similar to specifying an :code:`init_path`.
In contrast, starting values are taken from the previous calculation. Specifying
the :code:`init_path` overwrites :code:`init_previous`.

Design mode
+++++++++++
The design mode is used to design your system and is always the first
calculation of your plant. **The offdesign calculation is always based on a**
**design calculation!** Obviously as you are designing the plant the way you
want, you are flexible to choose the parameters to specify. However, you can
not specify parameters that are based on a design case, as for example the
isentropic efficiency characteristic function of a turbine or a pump.
Specifying a value for the efficiency is of course possible.

Offdesign mode
++++++++++++++
The offdesign mode is used to **calculate the performance of your plant, if**
**parameters deviate from the plant's design point**. This can be partload
operation, operation at different temperature or pressure levels etc.. Thus,
before starting an offdesign calculation you have to design your plant first.
By stating :code:`'offdesign'` as calculation mode, **components and**
**connections will switch to the offdesign mode.** This means that all
parameters provided as design parameters will be unset and all parameters
provided as offdesign parameters will be set instead. You can specify a
connection's or component's (off-)design parameters using the
:code:`set_attr` method.

For example, for a condenser you would usually design it to a maximum terminal
temperature difference, in offdesign the heat transfer coefficient is selected.
The heat transfer coefficient is calculated in the preprocessing of the
offdesign case based on the results from the design-case. Of course, this
applies to all other parameters in the same way. Also, the pressure drop is a
result of the geometry for the offdesign case, thus we swap the pressure ratios
with geometry independent zeta :math:`\frac{\zeta}{D^4}` values.

.. code-block:: python

    mycomponent.set_attr(
        design=['ttd_u', 'pr1', 'pr2'], offdesign=['UA', 'zeta1_d4', 'zeta2_d4']
    )

.. note::

    Some parameters come with characteristic functions based on the design case
    properties. This means, that e.g. the isentropic efficiency of a turbine
    is calculated as function of the actual mass flow to design mass flow
    ratio. You can provide your own (measured) data or use the already existing
    data from TESPy. All standard characteristic functions are available at
    :ref:`data_label`.

For connections it works in the same way, e.g. write

.. code-block:: python

    myconnection.set_attr(design=['h'], offdesign=['T'])

if you want to replace the enthalpy with the temperature for your offdesign.
The temperature is a result of the design calculation and that value is then
used for the offdesign calculation in this example.

To solve your offdesign calculation, use:

.. code-block:: python

    nw.solve(mode='offdesign', design_path='path/to/designpoint.json')

Component-level design references
*********************************

In some situations a single component - or a small group of components - has
been redesigned or replaced, so its individual design point differs from the
one captured in the network-wide :code:`design_path`. TESPy supports two
complementary mechanisms to handle this.

**Individual design path in an offdesign simulation**

You can assign a separate :code:`design_path` directly to any component.
During the offdesign preprocessing TESPy will then load that component's
design values (e.g. the heat-transfer coefficient :code:`UA` of a heat
exchanger or the isentropic efficiency of a compressor) from the component's
own path, while every other component keeps using the global
:code:`design_path`. The design values of all connections adjacent to that
component are also read from the component's individual path and stored
internally; they serve as the reference for characteristic functions that
scale with mass-flow or similar.

.. code-block:: python

    # Heat exchanger was replaced; new design case is in a separate file.
    heatex.set_attr(design_path='path/to/new_heatex_design.json')
    nw.solve(mode='offdesign', design_path='path/to/network_design.json')

Setting :code:`design_path=None` on a component reverts it to the global
reference on the next solve.

.. note::

    It is possible to use an individual design from an isolated sub-network
    that contains only the component of interest. Connection labels or even
    component labels in that file do not need to match the labels in the main
    network: if no label match is found TESPy tries to identify the correct
    entry via the port topology (inlet/outlet IDs) stored in the file.
    This topology information is only available in design files exported with
    TESPy 0.9.15 or higher. Older exports require matching labels.

**Local offdesign component in a design simulation**

It is also possible to keep most of the network in design mode, while one
component that has already been built should be treated in offdesign mode.
Setting :code:`local_offdesign=True` together with a :code:`design_path` on
that component achieves this. The component's offdesign equations become active
for that solve, and its design values and adjacent-connection references are
loaded from its own :code:`design_path`.

.. code-block:: python

    # Existing heat exchanger is frozen in offdesign; rest of network is sized.
    heatex.set_attr(
        local_offdesign=True,
        design_path='path/to/existing_heatex_design.json'
    )
    nw.solve(mode='design')

Solving and debugging
---------------------
Everything about the solution process - how the variable space is reduced,
how starting values are generated, how the equation system is decomposed
into blocks and solved, how to interact with a running solution process and
how to debug non-converging models - is described in the
:ref:`solver section <modules_solver_label>`.
