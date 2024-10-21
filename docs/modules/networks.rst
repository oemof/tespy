.. _tespy_modules_networks_label:

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
* the mass fractions of the fluids

of every connection.

The solver will simplify the variable space in a presolving step and then solve
for the remaining variables. If your **system includes fluid mixtures**, you
might want to **make use of the value ranges** for the system variables. This
improves the stability of the algorithm. Try to fit the boundaries as tight as
possible, for instance, if you know that the maximum pressure in the system will
be at 10 bar, use it as upper boundary.

.. note::

    Value ranges for pure fluids are not required as these are dealt with
    automatically.

.. code-block:: python

    >>> from tespy.networks import Network

    >>> my_plant = Network()
    >>> my_plant.p_unit
    'Pa'
    >>> my_plant.set_attr(p_unit='bar', h_unit='kJ / kg')
    >>> my_plant.set_attr(p_range=[0.05, 10], h_range=[15, 2000])
    >>> my_plant.p_unit
    'bar'
    >>> my_plant.p_range_SI
    [5000.0, 1000000.0]

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
    ...     screen_level=logging.INFO, file_level=logging.DEBUG
    ... );()  # +doctest: ELIPSIS
    (...)

The log-file will be saved to :code:`~/.tespy/log_files/` by default. All
available options are documented in the
:py:func:`API <tespy.tools.logger.define_logging>`.

Prior to solving the network there are options regarding the **console**
**printouts for the calculation progress**. Specify, if you want to enable or
disable convergence progress printouts:

.. code-block:: python

    >>> my_plant.iterinfo
    True

    # disable iteration information printout
    >>> my_plant.set_attr(iterinfo=False)
    >>> my_plant.iterinfo
    False

    # enable iteration information printout
    >>> my_plant.set_attr(iterinfo=True)

Adding connections
++++++++++++++++++
As seen in the introduction, you will have to create your networks from the
components and the connections between them. You can add connections directly
or via subsystems using the corresponding methods:

.. code-block:: python

    >>> my_plant.add_conns()
    >>> my_plant.add_subsys()

.. note::

    You do not need to add the components to the network, as they are inherited
    via the added connections. After having set up your network and added all
    required elements, you can start the calculation.

Busses: Energy Connectors
+++++++++++++++++++++++++
Another type of connection is the bus: Busses are connections for massless
transfer of energy e.g. in turbomachinery or heat exchangers. They can be used
to model motors or generators, too. Add them to your network with the following
method:

.. code-block:: python

    >>> my_plant.add_busses()

You will learn more about busses and how they work in
:ref:`this part <tespy_busses_label>`.

Start calculation
^^^^^^^^^^^^^^^^^
You can start the solution process with the following line:

.. code-block:: python

    my_plant.solve(mode='design')

This starts the initialisation of your network and proceeds to its calculation.
The specification of the **calculation mode is mandatory**, This is the list of
available keywords:

- :code:`mode` is the calculation mode (:code:`'design'`-calculation or
  :code:`'offdesign'`-calculation).
- :code:`init_path` is the path to the network folder you want to use for
  initialisation.
- :code:`design_path` is the path to the network folder which holds the
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
with zeta values.

.. code-block:: python

    mycomponent.set_attr(
        design=['ttd_u', 'pr1', 'pr2'], offdesign=['kA', 'zeta1', 'zeta2']
    )

.. note::

    Some parameters come with characteristic functions based on the design case
    properties. This means, that e.g. the isentropic efficiency of a turbine
    is calculated as function of the actual mass flow to design mass flow
    ratio. You can provide your own (measured) data or use the already existing
    data from TESPy. All standard characteristic functions are available at
    :ref:`tespy_data_label`.

For connections it works in the same way, e.g. write

.. code-block:: python

    myconnection.set_attr(design=['h'], offdesign=['T'])

if you want to replace the enthalpy with the temperature for your offdesign.
The temperature is a result of the design calculation and that value is then
used for the offdesign calculation in this example.

To solve your offdesign calculation, use:

.. code-block:: python

    my_plant.solve(mode='offdesign', design_path='path/to/network_designpoint')

Solving
-------
A TESPy network can be represented as a linear system of nonlinear equations,
consequently the solution is obtained with numerical methods. TESPy uses the
n-dimensional Newton-Raphson method to find the system's solution, which may
only be found, if the network is parameterized correctly. **The number of**
**variables n changes depending on your system's topology and your**
**specifications**. On top of that, the presolver reduces the number of
variables based on your model structure and your specifications.

**General preprocessing**

* check network consistency and initialise components (if network topology is
  changed to a prior calculation only).
* create a topology representation of the components and the connections.
* simplify the variable space based on the plant's topology.
* perform design/offdesign switch (for offdesign calculations only).
* preprocessing of offdesign case using the information from the
  :code:`design_path` argument.

The network check is used to find errors in the network topology, the
calculation can not start without a successful check. The design/offdesign
switch is described in the network setup section. For offdesign calculation the
:code:`design_path` argument is required. The design point information is
extracted from that path in preprocessing. For this, you will need to save
your network's design point information using:

.. code-block:: python

    my_plant.save('path/for/savestate')

**Simplifying the variable space**

To reduce the size of the system of equations a reduction of the variable space
is performed in the initialisation of a calculation. For every of the primary
variables (mass flow, pressure, enthalpy and fluid mass fractions), if a value
is directly specified by the user, the respective variable is removed from the
variable space, because it does not need to be solved.

Furthermore, there are three steps to simplify the variable space, i.e.
regarding mass flow, regarding the fluid composition and regarding pressure and
enthalpy.

First, based on the topology of the network, different branches are created.
These are:

- branches, in which the mass flow is equal in every of its connections and
- branches, in which the fluid composition is equal in every of its connections.

For every mass flow branch, the variable space is reduced to a single mass flow.
For example, in a simple Clausius Rankine cycle there will only be a single
mass flow in the variable space. Analogously in every fluid composition branch,
the variable space is reduced to a single vector containing the variable fluids
of that branch. For example, if a mass flow is split in two streams using a
splitter, the fluid composition remains constant downstream of the splitter.
Therefore, all connections downstream of the splitter share the same fluid
composition as upstream of the splitter.

The next step is a reduction of the fluid vector specifications: Consider a case
with a couple of potential fluids on a fluid branch, e.g. oxygen, nitrogen,
argon, carbon dioxide and water at the outlet of a combustion chamber. All fluid
mass fractions specified by the user will be fixed and removed from the variable
space. If then, only a single fluid remains with "unknown" mass fraction, we can
assign a mass fraction to that fluid, which is equal to 1 minus the sum of all
other fluids' mass fractions.

Finally, presolving is applied to pressure and enthalpy, whenever the fluid
composition is fixed. If either pressure or enthalpy is specified by the user
and on top of that temperature, vapor quality or temperature difference to
saturation temperature, the respective variable (enthalpy or pressure) can
directly be calculated. Similarly, if temperature and temperature difference to
saturation temperature or vapor quality are specified, both pressure and
enthalpy can be deducted.

**Finding starting values**

The algorithm requires starting values for all variables of the system, thus an
initialisation of the system is run prior to calculating the solution. **High**
**quality initial values are crucial for convergence speed and stability**, bad
starting values might lead to instability and diverging calculation can be the
result. The following steps are performed in finding starting values:

* fluid composition guessing.
* fluid property initialisation.
* initialisation from previous simulation run (:code:`init_previous`).
* initialisation from .csv (setting starting values from :code:`init_path`
  argument).

Starting value generation for your calculations starts with the fluid
composition guessing in case the fluid composition is not fixed. The available
fluids will be assigned the same mass fraction :math:`x`, if no starting value
is supplied. The mass fractions are distributed to 1 minus the sum of all user
specified mass fractions: :math:`x=\frac{1-\sum\text{x_spec}}{n}`. If you are
using combustion chambers these will be replaced by a generic flue gas
composition will be calculated prior to the propagation.

Next the fluid property initialisation uses user specified starting values or
the results from the previous simulation to set starting values for mass flow,
pressure and enthalpy. Otherwise, generic starting values are generated on basis
of which components a connection is linked to. If you **do not want** to use the
results of a previous calculation, you need to specify
:code:`init_previous=False` on the :code:`Network.solve` method call.

Last step in starting value generation is the initialisation from a saved
network state. In order to initialise your calculation with this method, you
need to provide the path to the saved network in the :code:`init_path` argument
of the `solve` method. TESPy searches through the connections.csv file. If a
connection with the respective label is found, the starting values for the
system variables are taken over from that file.

.. note::

    The files do not need to contain all connections of your network. You can
    build your network step by step and initialise the existing parts of your
    network from the :code:`init_path`. Be aware that a change within the fluid
    vector does not allow this practice! If you plan to use additional fluids
    in parts of the network you have not touched until now, you will need to
    state all fluids from the beginning.

Algorithm
^^^^^^^^^
In this section we will give you an introduction to the solving algorithm
implemented.

Newton-Raphson method
+++++++++++++++++++++
The Newton-Raphson method requires the calculation of residual values for the
equations and of the partial derivatives to all system variables (Jacobian
matrix). In the next step the matrix is inverted and multiplied with the
residual vector to calculate the increment for the system variables. This
process is repeated until every equation's result in the system is "correct",
thus the residual values are smaller than a specified error tolerance. All
equations are of the same structure:

.. math::

    0 = \text{expression}

calculate the residuals

.. math::

    f(\vec{x}_i)

Jacobian matrix J

.. math::

    J(\vec{x})=\left(\begin{array}{cccc}
    \frac{\partial f_1}{\partial x_1} & \frac{\partial f_1}{\partial x_2} &
    \cdots & \frac{\partial f_1}{\partial x_n} \\
    \frac{\partial f_2}{\partial x_1} & \frac{\partial f_2}{\partial x_2} &
    \cdots & \frac{\partial f_2}{\partial x_n} \\
    \vdots & \vdots & \ddots & \vdots \\
    \frac{\partial f_n}{\partial x_1} & \frac{\partial f_n}{\partial x_2} &
    \cdots & \frac{\partial f_n}{\partial x_n}
    \end{array}\right)

derive the increment

.. math::

    \vec{x}_{i+1}=\vec{x}_i-J(\vec{x}_i)^{-1}\cdot f(\vec{x}_i)

while

.. math::

    ||f(\vec{x}_i)|| > \epsilon

.. note::

    You have to provide the exact amount of required parameters (neither less
    nor more) and the parametrisation must not lead to linear dependencies.
    Each parameter you set for a connection and each energy flow you specify
    for a bus will add one equation to your system. On top, each component
    provides a different amount of basic equations plus the equations provided
    by your component specification.

For example, consider a pump: Total mass flow as well as the fluid mass
fractions of the mixture entering the pump will be identical at the outlet. The
pump delivers two mandatory equations. If you additionally specify, e.g. the
power :math:`P` to be 1000 W, the set of equations will look like this:

.. math::

    \forall i \in \mathrm{network.fluids} \, &0 = fluid_{i,in} -fluid_{i,out}\\
    &0 = \dot{m}_{in} - \dot{m}_{out}\\
    \mathrm{additional:} \, &0 = 1000 - \dot{m}_{in} (\cdot {h_{out} - h_{in}})

.. _tespy_modules_convergence_check_label:

Convergence stability
+++++++++++++++++++++
One of the main downsides of the Newton-Raphson method is that the initial
step width is very large and that it does not know physical boundaries, for
example mass fractions smaller than 0 and larger than 1 or negative pressure.
Also, the large step width can adjust enthalpy or pressure to quantities that
are not covered by the fluid property databases. This would cause an inability
e.g. to calculate a temperature from pressure and enthalpy in the next
iteration of the algorithm. In order to improve convergence stability, we have
added a convergence check.

**The convergence check manipulates the system variables after the increment**
**has been added**. This manipulation has four steps, the first two are always
applied:

* Cut off fluid mass fractions smaller than 0 and larger than 1. This way a
  mass fraction of a single fluid component never exceeds these boundaries.
* Check, whether the fluid properties of pure fluids are within the available
  ranges of CoolProp and readjust the values if not.

The next two steps are applied, if the user did not specify an
:code:`init_path` and the iteration count is lower than 3, thus in the first
three iteration steps of the algorithm only. In other cases this convergence
check is skipped.

* Fox mixtures: check, if the fluid properties (pressure, enthalpy and mass
  flow) are within the user specified boundaries
  (:code:`p_range, h_range, m_range`) and if not, cut off higher/lower values.
* Check the fluid properties of the connections based on the components they
  are connecting. For example, check if the pressure at the outlet of a turbine
  is lower than the pressure at the inlet or if the flue gas composition at a
  combustion chamber's outlet is within the range of a "typical" flue gas
  composition. If there are any violations, the corresponding variables are
  manipulated. If you want to look up, what exactly the convergence check for a
  specific component does, look out for the :code:`convergence_check` methods
  in the :py:mod:`tespy.components module <tespy.components>`.

In a lot of different tests the algorithm has found a near enough solution
after the third iteration, further checks are usually not required.

.. tip::

    To check if the solver successfully found a solution for your model you can
    check the `.converged` attribute of the Network class after calling the
    `solve` method. It will be `True` in case no linear dependency was and the
    residual value of all equations is below the minimum threshold.

Calculation speed improvement
+++++++++++++++++++++++++++++
For improvement of calculation speed, the calculation of specific derivatives
is skipped, if the change of the corresponding variable was below a
threshold of :code:`1e-12` in the iteration before.

As a user you can take two more measures to improve calculation speed: Specify
primary variables whenever possible/reasonable. This will not only reduce the
variable space but also remove the necessity to calculate partial derivatives
towards them.

Troubleshooting
+++++++++++++++
In this section we show you how you can troubleshoot your calculation and list
up common mistakes. If you want to debug your code, make sure to enable the
logger and have a look at the log-file at :code:`~/.tespy/` (or at your
specified location).

First, make sure your network topology is set up correctly, TESPy will prompt
an Error, if not. TESPy will prompt an error, too, if you did not provide
enough or if you provide too many parameters for your calculation, but you will
not be given an information which specific parameters are under- or
over-determined.

.. note::

    Always keep in mind, that the system has to find a value for mass flow,
    pressure, enthalpy and the fluid mass fractions. Try to build up your
    network step by step and have in mind, what parameters will be determined
    by adding a component without any parametrisation. This way, you can easily
    determine, which parameters are still to be specified.

If you are modeling a cycle, e.g. the Clausius Rankine cylce, you need to make
a cut in the cycle using the :code:`CycleCloser` or a :code:`Sink` and a
:code:`Source` not to over-determine the system. Have a look in the
:ref:`tutorial section <tespy_basics_label>` to understand why this is
important and how it can be implemented.

If you have provided the correct number of parameters in your system and the
calculations stops after or even before the first iteration, there might be a
couple of reasons for that:

- Sometimes, the fluid property database does not find a specific fluid
  property in the initialisation process, have you specified the values in the
  correct unit?
- A linear dependency in the Jacobian matrix due to bad parameter settings
  stops the calculation (over-determining one variable, while missing out on
  another).
- A linear dependency in the Jacobian matrix due to bad starting values stops
  the calculation.

The first reason can be eliminated by carefully choosing the parametrisation.
**A linear dependency due to bad starting values is often more difficult to**
**resolve, and it may require some experience.** In many cases, the linear
dependency is caused by equations, that require the **calculation of a**
**temperature**, e.g. specifying a temperature at some point of the network,
terminal temperature differences at heat exchangers, etc.. In this case,
**the starting enthalpy and pressure should be adjusted in a way, that the**
**fluid state is not within the two-phase region:** The specification of
temperature and pressure in a two-phase region does not yield a distinct value
for the enthalpy. Even if this specific case appears after some iterations,
better starting values often do the trick.

Another frequent error is that fluid properties move out of the bounds given by
the fluid property database. The calculation will stop immediately.
**Adjusting pressure and enthalpy ranges for the convergence check** might help
in this case.

.. note::

    If you experience slow convergence or instability within the convergence
    process, it is sometimes helpful to have a look at the iteration
    information. This is printed by default and provides information on the
    residuals of your systems' equations and on the increments of the systems'
    variables. Maybe it is only one variable causing the instability, its
    increment is much larger than the increment of the other variables?

If you run multiple simulations and a simulation crashed due to an internal
error (e.g. fluid property related), and you still want the next simulations to
perform correctly, you have to call the
:py:meth:`tespy.networks.network.Network.reset_topology_reduction_specifications`
method after the failed simulation and before you run the next simulation.
Usually, this happens automatically as part of the post-processing, but in case
the simulation crashed before that, this step cannot be executed. Then,
restarting the simulation is not possible.

Did you experience other errors frequently and have a workaround/tips for
resolving them? You are very welcome to contact us and share your experience
for other users!

Post-processing
---------------
A post-processing is performed automatically after the calculation finished. You
have further options:

- Automatically create a documentation of your model.
- Print the results to prompt (:code:`print_results()`).
- Save the results in structure of .csv-files (:code:`save()`).
- Generate fluid property diagrams with an external tool.

Automatic model documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Using the automatic TESPy model documentation you can create an overview of
all input parameters, specifications and equations as well as characteristics
applied in LaTeX format. This enables high

- **transparency**,
- **readability** and
- **reproducibility**.

In order to use the model documentation, you need to import the corresponding
method and pass your network information. At the moment, you can the following
optional arguments to the method:

- :code:`path`: Basepath, where the LaTeX data and figures are exported to.
- :code:`filename`: Filename of the report.
- :code:`fmt`: A formatting dictionary, for a sample see below.

.. code-block:: python

    from tespy.tools import document_model

    fmt = {
        'latex_body': True,  # adds LaTeX body to compile report out of the box
        'include_results': True,  # include parameter specification and results
        'HeatExchanger': {  # for components of class HeatExchanger
            'params': ['Q', 'ttd_l', 'ttd_u', 'pr1', 'pr2']},  # change columns displayed
        'Condenser': {  # for components of class HeatExchanger
            'params': ['Q', 'ttd_l', 'ttd_u', 'pr1', 'pr2']
            'float_fmt': '{:,.2f}'},  # change float format of data
        'Connection': {  # for Connection instances
            'p': {'float_fmt': '{:,.4f}'},  # change float format of pressure
            's': {'float_fmt': '{:,.4f}'},
            'h': {'float_fmt': '{:,.2f}'},
            'params': ['m', 'p', 'h', 's']  # list results of mass flow, ...
            'fluid': {'include_results': False}  # exclude results of fluid composition
        },
        'include_results': True,  # include results
        'draft': False  # disable draft mode
    }
    document_model(mynetwork, fmt=fmt)

.. note::

    Specified values are displayed in any case. The selection of which
    parameters to show and which to exclude only applies to results.

After having exported the LaTeX code, you can simply use :code:`\input{}`
in your main LaTeX document to include the documentation of your model. In
order to compile correctly you need to load the following LaTeX packages:

* graphicx
* float
* hyperref
* booktabs
* amsmath
* units
* cleveref
* longtable

For generating different file formats, like markdown, html or
restructuredtext, you could try the `pandoc <https://pandoc.org/>`_ library.
For examples, of how the reports look you can have a look at the
`examples <https://github.com/oemof/oemof-examples/tree/master/oemof_examples/tespy>`_
repository, or just try it yourself :).

This feature is introduced in version 0.4.0 and still subject to changes. If
you have any suggestions, ideas or feedback, you are very welcome to submit an
issue on our GitHub or even open a pull request.

Results printing
^^^^^^^^^^^^^^^^
To print the results in your console use the :code:`print_results()` method.
It will print tables containing the component, connection and bus properties.
Some results will be colored, the colored results indicate

* if a parameter was specified as value before calculation.
* if a parameter is out of its predefined value bounds (e.g. efficiency > 1).
* if a component parameter was set to :code:`'var'` in your calculation.

The color for each of those categories is different and might depend on the
console settings of your machine. If you do not want the results to be colored
you can instead call the method the following way:

.. code-block:: python

    my_plant.print_results(colored=False)

If you want to limit your printouts to a specific subset of components,
connections and busses, you can specify the :code:`printout` parameter to block
individual result printout.

.. code-block:: python

    mycomp.set_attr(printout=False)
    myconn.set_attr(printout=False)
    mybus.set_attr(printout=False)

If you want to prevent all printouts of a subsystem, add something like this:

.. code-block:: python

    # connections
    for c in mysubsystem.conns.values():
        c.set_attr(printout=False)

    # components
    for c in mysubsystem.comps.values():
        c.set_attr(printout=False)

Save your results
^^^^^^^^^^^^^^^^^
If you choose to save your results the specified folder will be created
containing information about the network, all connections, busses, components
and characteristics.

In order to perform calculations based on your results, you can access all
components' and connections' parameters:

The easiest way to access the results of one specific component looks like this

.. code:: python

    eff = mycomp.eta_s.val  # isentropic efficiency of mycomp
    P = mycomp.P.val

and similar for connection parameters:

.. code:: python

    mass_flow = myconn.m.val  # value in specified network unit
    mass_flow_SI = myconn.m.val_SI  # value in SI unit
    mass_fraction_oxy = myconn.fluid.val['O2']  # mass fraction of oxygen
    specific_volume = myconn.vol.val  # value in specified network unit
    specific_entropy = myconn.s.val  # value in specified network unit
    volumetric_flow = myconn.v.val  # value in specified network unit
    specific_exergy = myconn.ex_physical  # SI value only

On top of that, you can access pandas DataFrames containing grouped results
for the components, connections and busses. The instance of class Network
provides a results dictionary.

.. code:: python

    # key for connections is 'Connection'
    results_for_conns = my_plant.results['Connection']
    # keys for components are the respective class name, e.g.
    results_for_turbines = my_plant.results['Turbine']
    results_for_heat_exchangers = my_plant.results['HeatExchanger']
    # keys for busses are the labels, e.g. a Bus labeled 'power input'
    results_for_mybus = my_plant.results['power input']

The index of the DataFrames is the connection's or component's label.

.. code:: python

    results_for_specific_conn = my_plant.results['Connection'].loc['myconn']
    results_for_specific_turbine = my_plant.results['Turbine'].loc['turbine 1']
    results_for_component_on_bus = my_plant.results['power input'].loc['turbine 1']

The full list of connection and component parameters can be obtained from the
respective API documentation.

Network reader
==============
The network reader is a useful tool to import networks from a data structure
using .csv-files. In order to re-import an exported TESPy network, you must
save the network first.

.. code:: python

    my_plant.export('mynetwork')

This generates a folder structure containing all relevant files defining your
network (general network information, components, connections, busses,
characteristics) holding the parametrization of that network. You can re-import
the network using following code with the path to the saved documents. The
generated network object contains the same information as a TESPy network
created by a python script. Thus, it is possible to set your parameters in the
.csv-files, too. The imported network is handled identically as a manually
created network.

.. code:: python

    from tespy.networks import load_network
    imported_plant = load_network('path/to/mynetwork')
    imported_plant.solve('design')

.. note::

    Imported busses, components and connections are accessible by their label,
    e.g. :code:`imported_plant.busses['total heat output']`,
    :code:`imported_plant.get_comp('condenser')` and
    :code:`imported_plant.get_conn('myconnectionlabel')` respectively. If
    you did not provide labels for your connections, by default, the
    connection's label will be according to this principle:
    :code:`'source-label:source-id_target-label:target-id'`, where source and
    target are the labels of the connected components.
