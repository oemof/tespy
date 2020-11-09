.. _using_tespy_networks_label:

Networks
========
The network class handles preprocessing, solving and postprocessing.
We will walk you through all the important steps.

Setup
-----
Network container
^^^^^^^^^^^^^^^^^
The TESPy network contains all data of your plant, which in terms of the
calculation is represented by a nonlinear system of equations. The system
variables of your TESPy network are:

 * mass flow,
 * pressure,
 * enthalpy and
 * the mass fractions of the network's fluids.

The solver will solve for these variables. As stated in the introduction the
list of fluids is passed to your network on creation. If your **system**
**includes fluid mixtures**, you should **always make use of the value ranges**
for the system variables. This improves the stability of the algorithm. Try to
fit the boundaries as tight as possible, for instance, if you know that the
maximum pressure in the system will be at 10 bar, use it as upper boundary.

.. note::

    Value ranges for pure fluids are not required as these are dealt with
    automatically.

.. code-block:: python

    from tespy.networks import network

    fluid_list = ['CO2', 'H2O', 'N2', 'O2', 'Ar']
    my_plant = network(fluids=fluid_list)
    my_plant.set_attr(p_unit='bar', h_unit='kJ / kg')
    my_plant.set_attr(p_range=[0.05, 10], h_range=[15, 2000])

.. note::

    It is possible to specify the fluid property back-end of the network fluids
    by adding the name of the back-end in front of the fluid's name within the
    list of fluids of your network. For example:

    .. code-block:: python

        from tespy.networks import network

        fluid_list = ['CO2', 'BICUBIC::H2O', 'INCOMP::DowQ']
        network(fluids=fluid_list)

    If you do not specify a back-end, the **default back-end** :code:`HEOS`
    will be used (as for :code:`CO2`). In this example, :code:`H2O` will be
    using the :code:`BICUBIC` back-end and :code:`DowQ` the back-end for
    incompressible fluids :code:`INCOMP`. For an overview of the back-ends
    available please refer to the
    :ref:`fluid property section <tespy_fluid_properties_label>`.

    **It is not possible to use the same fluid in different back-ends!**

.. _printout_logging_label:

Printouts and logging
+++++++++++++++++++++
TESPy comes with an inbuilt logger. If you want to keep track of
debugging-messages, general information, warnings or errors you should enable
the logger. At the beginning of your python script e.g. add the following
lines:

.. code-block:: python

    from tespy.tools import logger
    import logging
    logger.define_logging(
        log_path=True, log_version=True,
        screen_level=logging.INFO, file_level=logging.DEBUG
    )

The log-file will be saved to :code:`~/.tespy/log_files/` by default. All
available options are documented in the
:py:func:`API <tespy.tools.logger.define_logging>`.

Prior to solving the network there are options regarding the **console**
**printouts for the calculation progress**. Specify, if you want to enable or
disable convergence progress printouts:

.. code-block:: python

    # disable iteration information printout
    myplant.set_attr(iterinfo=False)
    # enable iteration information printout
    myplant.set_attr(iterinfo=True)

Adding connections
++++++++++++++++++
As seen in the introduction, you will have to create your networks from the
components and the connections between them. You can add connections directly
or via subsystems and networks holding them by using the appropriate methods:

.. code-block:: python

    myplant.add_conns()
    myplant.add_subsys()
    myplant.add_nwks()

.. note::

    You do not need to add the components to the network, as they are inherited
    via the added connections. After having set up your network and added all
    required elements, you can start the calculation.

Busses: power connections
+++++++++++++++++++++++++
Another type of connection is the bus: Busses are power connections for e.g.
turbomachines or heat exchangers. They can be used to model motors or
generators, too. Add them to your network with the following method:

.. code-block:: python

    myplant.add_busses()

You will learn more about busses and how they work in
:ref:`this part <tespy_busses_label>`.

Start calculation
^^^^^^^^^^^^^^^^^
You can start the solution process with the following line:

.. code-block:: python

    myplant.solve(mode='design')

This starts the initialisation of your network and proceeds to its calculation.
The specification of the **calculation mode is mandatory**, This is the list of
available keywords:

 * :code:`mode` is the calculation mode (:code:`"design"`-calculation or
   :code:`"offdesign"`-calculation).
 * :code:`init_path` is the path to the network folder you want to use for
   initialisation.
 * :code:`design_path` is the path to the network folder which holds the
   information of your plant's design point.
 * :code:`max_iter` is the maximum amount of iterations performed by the
   solver.
 * :code:`min_iter` is the minimum amount of iterations before a solution can
   be accepted (given the convergence criterion is satisfied).
 * :code:`init_only` stop after initialisation (True/False).
 * :code:`init_previous` use starting values from previous simulation
   (True/False).

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
the :code:`ìnit_path` overwrites :code:`init_previous`.

Design mode
+++++++++++
The design mode is used to design your system and is always the first
calculation of your plant. **The offdesign calculation is always based on a**
**design calculation!**. Obviously as you are designing the plant the way you
want, you are flexible to choose the parameters to specify. However, you can
not specify parameters that are based on a design case, as for example the
isentropic efficiency characteristic function of a turbine or a pump.
Specifying a value for the efficiency is of course possible.

Offdesign mode
++++++++++++++
The offdesign mode is used to **calculate the performance of your plant, if**
**parameters deviate from the plant's design point**. This can be partload
operation, operation at different temperature or pressure levels etc.. Thus,
before starting an offdesing calculation you have to design your plant first.
By stating :code:`'offdesign'` as calculation mode, **components and**
**connections will switch to the offdesign mode.** This means that all
parameters provided as design parameters will be unset and all parameters
provided as offdesign parameters will be set instead. You can specify a
connection's or component's (off-)design parameters using the set_attr method.

You can specify design and offdesign parameters for components and connections.
For example, for a condenser you would usually design it to a maximum terminal
temperature difference, in offdesign the heat transfer coefficient is selected.
The heat transfer coefficient is calculated in the preprocessing of the
offdesign case based on the results from the design-case. Of course, this
applies to all other parameters in the same way. Also, the pressure drop is a
result of the geometry for the offdesign case, thus we swap the pressure ratios
with zeta values.

.. code-block:: python

    mycomponent.set_attr(design=['ttd_u', 'pr1', 'pr2'],
                         offdesign=['kA', 'zeta1', 'zeta2'])

.. note::

    Some parameters come with characteristic functions based on the design case
    properties. This means, that e.g. the isentropic efficiency of a turbine
    is calculated as function of the actual mass flow to design mass flow
    ratio. You can provide your own (measured) data or use the already existing
    data from TESPy. All standard characteristic functions are available at
    :py:mod:`tespy.data`.

For connections it works in the same way, e.g. write

.. code-block:: python

    myconnection.set_attr(design=['h'], offdesign=['T'])

if you want to replace the enthalpy with the temperature for your offdesign.
The temperature is a result of the design calculation and that value is then
used for the offdesign calculation in this example.

To solve your offdesign calculation, use:

.. code-block:: python

    myplant.solve(mode='offdesign', design_path='path/to/network_designpoint')

Solving
-------
A TESPy network can be represented as a linear system of nonlinear equations,
consequently the solution is obtained with numerical methods. TESPy uses the
n-dimensional Newton–Raphson method to find the systems solution, which may
only be found, if the network is parameterized correctly. **The number of
variables n** is :math:`n = num_{conn} \cdot (3 + num_{fluids})`.

The algorithm requires starting values for all variables of the system, thus an
initialisation of the system is run prior to calculating the solution. **High**
**quality initial values are crutial for convergence speed and stability**, bad
starting values might lead to instability and diverging calculation can be the
result. There are different levels for the initialisation.

Initialisation
^^^^^^^^^^^^^^
The initialisation is performed in the following steps.

**General preprocessing:**

 * check network consistency and initialise components (if network topology is
   changed to a prior calculation only).
 * perform design/offdesign switch (for offdesign calculations only).
 * preprocessing of offdesign case using the information from the
   :code:`design_path` argument.

**Finding starting values:**

 * fluid propagation.
 * fluid property initialisation.
 * initialisation from previous simulation run (:code:`ìnit_previous`).
 * initialisation from .csv (setting starting values from :code:`init_path`
   argument).

The network check is used to find errors in the network topology, the
calculation can not start without a successful check. For components, a
preprocessing of some parameters is necessary. It is performed by the
:code:`comp_init` method of the components. You will find the methods in the
:py:mod:`components module <tespy.components>`. The design/offdesign switch is
described in the network setup section. For offdesign calculation the
:code:`design_path` argument is required. The design point information is
extracted from that path in preprocessing. For this, you will need to export
your network's design point information using:

.. code-block:: python

    myplant.save('path/for/export')

Starting value generation for your calculations starts with the fluid
propagation. **The fluid propagation is a very important step in the**
**initialisation.** Often, you will specify the fluid at one point of the
network only, all other connections are missing an initial information on the
fluid, if you are not using an :code:`init_path`. The fluid propagation will
push/pull the specified fluid through the network. If you are using combustion
chambers these will be starting points and a generic flue gas composition will
be calculated prior to the propagation. You do not necessarily need to state a
starting value for the fluid at every point of the network.

.. note::
    If the fluid propagation fails, you often experience an error, where the
    fluid property database can not find a value, because the fluid is 'nan'.
    Providing starting values manually can fix this problem.

If available, the fluid property initialisation uses the user specified starting
values or the results from the previous simulation. Otherwise generic starting
values are generated on basis of which components a connection is linked to.
If you do not want to use the results of a previous calculation, you need to
specify :code:`init_previous=False` on the :code:`network.solve` method call.

Last step in starting value generation is the initialisation from a saved
network structure. In order to initialise your calculation from the
:code:`init_path`, you need to provide the path to the saved/exported network.
If you specify an :code:`init_path` TESPy searches through the connections file
for the network topology and if the corresponding connection is found, the
starting values for the system variables are extracted from the connections
file.

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

Newton–Raphson method
+++++++++++++++++++++
The Newton–Raphson method requires the calculation of residual values for the
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

jacobian matrix J

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

.. _using_tespy_convergence_check_label:

Convergence stability
+++++++++++++++++++++
One of the main downsides of the Newton–Raphson method is that the initial
stepwidth is very large and that it does not know physical boundaries, for
example mass fractions smaller than 0 and larger than 1 or negative pressure.
Also, the large stepwidth can adjust enthalpy or pressure to quantities that
are not covered by the fluid property databases. This would cause an inability
e.g. to calculate a temperature from pressure and enthalpy in the next
iteration of the algorithm. In order to improve convergence stability, we have
added a convergence check.

**The convergence check manipulates the system variables after the increment**
**has been added**. This manipulation has four steps, the first two are always
applied:

 * Cut off fluid mass fractions smaller than 0 and larger than 1. This way a
   mass fraction of a single fluid components never exceeds these boundaries.
 * Check, whether the fluid properties of pure fluids are within the available
   ranges of CoolProp and readjust the values if not.

The next two steps are applied, if the user did not specify an
:code:`init_path` and the iteration count is lower than 3, thus in the first
three iteration steps of the algorithm only. In other cases this convergence
check is skipped.

 * Fox mixtures: check, if the fluid properties (pressure, enthalpy and
   mass flow) are within the user specified boundaries
   (:code:`p_range, h_range, m_range`) and if not, cut off higher/lower
   values.
 * Check the fluid properties of the connections based on the components they
   are connecting. For example, check if the pressure at the outlet of a
   turbine is lower than the pressure at the inlet or if the flue gas
   composition at a combustion chamber's outlet is within the range of a
   "typical" flue gas composition. If there are any violations, the
   corresponding variables are manipulated. If you want to look up, what
   exactly the convergence check for a specific component does, look out for
   the :code:`convergence_check` methods in the
   :py:mod:`tespy.components module <tespy.components>`.

In a lot of different tests the algorithm has found a near enough solution
after the third iteration, further checks are usually not required.

.. note::

    It is possible to improve the convergence stability manually when using
    pure fluids. If you know the fluid's state is liquid or gaseous prior to
    the calculation, you may provide the according value for the keyword e. g.
    :code:`myconn.set_attr(state='l')`. The convergence check manipulates the
    enthalpy values so that the fluid is always in the desired state at that
    point. For an example see the release information of
    :ref:`version 0.1.1 <whats_new_011_example_label>`.

Calculation speed improvement
+++++++++++++++++++++++++++++
For improvement of calculation speed, the calculation of specific equations and
derivatives is skipped if possible. There are two criteria for equations and
one criterion for derivatives that are checked for calculation intensive
operations, e. g. whenever fluid property library calls are necessary:

For component equations the recalculation of the residual value is skipped,

- if the absolute of the residual value of that equations is lower than the
  threshold of :code:`1e-12` in the iteration before and simultaneously
- the iteration count is not a multiple of 4.

Connection equations are evaluated at least in every second iteration. If a
temperature value has been specified, the equation will be evaluated in every
iteration.

The calculation of derivatives is skipped, if the change of the corresponding
variable was below a threshold of :code:`1e-12` in the iteration before.
Again, this does not apply to temperature value specification, as especially
when using fluid mixtures, the convergence stability is very sensitive to
these equations and derivatives.

.. note::

    In order to make sure, that every equations is evaluated at least twice,
    the minimum amount of iterations before convergence can be accepted is at
    4.

Troubleshooting
+++++++++++++++
In this section we show you how you can troubleshoot your calculation and list
up common mistakes. If you want to debug your code, make sure to enable the
logger and have a look at the log-file at :code:`~/.tespy/` (or at your
specified location).

First of all, make sure your network topology is set up correctly, TESPy will
prompt an Error, if not. TESPy will prompt an error, too, if you did not
provide enough or if you provide too many parameters for your calculation, but
you will not be given an information which specific parameters are under- or
overdetermined.

.. note::
    Always keep in mind, that the system has to find a value for mass flow,
    pressure, enthalpy and the fluid mass fractions. Try to build up your
    network step by step and have in mind, what parameters will be determined
    by adding an additional component without any parametrisation. This way,
    you can easily determine, which parameters are still to be specified.

When using multiple fluids in your network, e. g.
:code:`fluids=['water', 'air', 'methane']` and at some point you want to have
water only, you still need to specify the mass fractions for both air and
methane (although beeing zero) at that point
:code:`fluid={'water': 1, 'air': 0, 'methane': 0}`. Also, setting
:code:`fluid={water: 1}, fluid_balance=True` will still not be sufficient, as
the fluid_balance parameter adds only one equation to your system.

If you are modeling a cycle, e. g. the clausius rankine cylce, you need to make
a cut in the cycle using the cycle_closer or a sink and a source not to
overdetermine the system. Have a look in the
:ref:`heat pump tutorial <heat_pump_tutorial_label>` to understand why this is
important and how it can be implemented.

If you have provided the correct number of parameters in your system and the
calculations stops after or even before the first iteration, there are four
frequent reasons for that:

 * Sometimes, the fluid property database does not find a specific fluid
   property in the initialisation process, have you specified the values in the
   correct unit?
 * Also, fluid property calculation might fail, if the fluid propagation
   failed. Provide starting values for the fluid composition, especially, if
   you are using drums, merges and splitters.
 * A linear dependency in the jacobian matrix due to bad parameter settings
   stops the calculation (overdetermining one variable, while missing out on
   another).
 * A linear dependency in the Jacobian matrix due to bad starting values stops
   the calculation.

The first reason can be eliminated by carefully choosing the parametrization.
**A linear dependency due to bad starting values is often more difficult to**
**resolve and it may require some experience.** In many cases, the linear
dependency is caused by equations, that require the **calculation of a**
**temperature**, e. g. specifying a temperature at some point of the network,
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

Did you experience other errors frequently and have a workaround/tips for
resolving them? You are very welcome to contact us and share your experience
for other users!

Postprocessing
--------------
A postprocessing is performed automatically after the calculation finished. You
have two further options:

 * print the results to prompt (:code:`myplant.print_results()`) and
 * save the results in a .csv-file (:code:`myplant.save('savename')`).

To print the results in your console use the :code:`print_results()` method.
It will print tables containing the component, connection and bus properties.
Some of the results will be colored, the colored results indicate

 * if a parameter was specified as value before calculation.
 * if a parameter is out of its predefined value bounds (e.g. efficiency > 1).
 * if a component parameter was set to :code:`'var'` in your calculation.

The color for each of those categories is different and might depend on the
console settings of your machine. You can change the colors using
`ANSI escape codes <https://en.wikipedia.org/wiki/ANSI_escape_code#Colors>`_,
for example changing the color for user specified values to blue:

.. code-block:: python

      from tespy.tools.global_vars import coloring
      coloring['set'] = '\034[0m'

If you do not want the results to be colored you can instead call the method
the following way:

.. code-block:: python

    myplant.print_results(colored=False)

If you want to limit your prinouts to a specific subset of components,
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

If you choose to save your results the specified folder will be created
containing information about the network, all connections, busses, components
and characteristics.

In order to perform calculations based on your results, you can access all
components' and connections' parameters:

For the components this is the way to go

.. code:: python

    eff = mycomp.eta_s.val  # isentropic efficiency of mycomp
    P = mycomp.P.val

Use this code for connection parameters:

.. code:: python

    mass_flow = myconn.m.val  # value in specified network unit
    mass_flow_SI = myconn.m.val_SI  # value in SI unit
    mass_fraction_oxy = myconn.fluid.val['O2']  # mass fraction of oxygen
    specific_volume = myconn.vol.val  # value in specified network unit
    specific_entropy = myconn.s.val  # value in specified network unit

.. _FluProDia_label:

Creating fluid property diagrams
--------------------------------

.. image:: //api/_images/logph_diagram_states.svg
    :align: center

Figure: logph diagram of NH3 with a simple heat pump cycle.

.. image:: //_images/Ts_diagram_states.svg

Figure: Ts diagram of NH3 with a simple heat pump cycle.

CoolProp has an inbuilt feature for creating fluid property diagrams.
Unfortunately, the handling is not very easy at the moment. We recommend using
fluprodia (Fluid Property Diagram) instead. You can create and customize
different types of diagrams for all pure and pseudo-pure fluids available in
CoolProp. In order to plot your process data into a diagram, you can use the
:code:`get_plotting_data` method of each component. The method returns a
dictionary, that can be passed as **kwargs to the
:code:`calc_individual_isoline` method of a fluprodia
:code:`FluidPropertyDiagram` object. The fluprodia documentation provides
examples of how to plot a process into different diagrams, too. For more
information on fluprodia have a look at the
`online documentation <https://fluprodia.readthedocs.io/en/latest/>`_. You can
install the package with pip.

.. code-block:: bash

    pip install fluprodia

.. note::

    The plotting data a returned from the :code:`get_plotting_data` as a
    nested dictionary. The first level key contains the connection id of the
    state change (change state from incoming connection to outgoing
    connection). The table below shows the state change and the respective id.

    .. list-table:: State change and respective ids of dictionary
       :widths: 60 10 10 10
       :header-rows: 1

       * - component
         - state from
         - state to
         - id
       * - components with one mass flow
         - :code:`in1`
         - :code:`out1`
         - :code:`1`
       * - class heat_exchanger and subclasses
         - :code:`in1`
         - :code:`out1`
         - :code:`1`
       * -
         - :code:`in2`
         - :code:`out2`
         - :code:`2`
       * - class orc_evaporator
         - :code:`in1`
         - :code:`out1`
         - :code:`1`
       * -
         - :code:`in2`
         - :code:`out2`
         - :code:`2`
       * -
         - :code:`in3`
         - :code:`out3`
         - :code:`3`
       * - class merge
         - :code:`in1`
         - :code:`out1`
         - :code:`1`
       * -
         - :code:`in2`
         - :code:`out1`
         - :code:`2`
       * -
         - ...
         - ...
         - ...
       * - class node
         - :code:`in1` if :math:`\dot{m}>0`
         - :code:`out1` if :math:`\dot{m}>0`
         - :code:`1`
       * -
         - ...
         - ...
         - ...

    - All other components do not return any information as either there is no
      change in state or the state change is accompanied by a change in fluid
      composition.
    - For class node the state change is from connections with incoming mass
      flow connections with outgoing mass flow.

Network reader
==============
The network reader is a useful tool to import networks from a datastructure
using .csv-files. In order to reimport an exported TESPy network, you must save
the network first.

.. code:: python

    myplant.save('mynetwork')

This generates a folder structure containing all relevant files defining your
network (general network information, components, connections, busses,
characteristics) holding the parametrization of that network. You can reimport
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
    e. g. :code:`imported_plant.busses['total heat output']`,
    :code:`imported_plant.components['condenser']` and
    :code:`imported_plant.connections['myconnectionlabel']` respectively. If
    you did not provide labels for your connections, by default, the
    connection's label will be according to this principle:
    :code:`'source-label:source-id_target-label:target-id'`, where source and
    target are the labels of the connected components.
