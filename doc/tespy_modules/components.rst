Components
==========

In this section we will introduce you into the details of component
parametrisation and component characteristics. At the end of the section we
show you, how to create custom components.

List of components
------------------

More information on the components can be gathered from the code documentation.
We have linked the base class containing a figure and basic informations as
well as the equations.

- Basics
    * :py:class:`Cycle closer<tespy.components.basics.cycle_closer>`
    * :py:class:`Sink <tespy.components.basics.sink>`
    * :py:class:`Source <tespy.components.basics.source>`
    * :py:class:`Subsystem interface <tespy.components.basics.subsystem_interface>`
- Combustion
    * :py:class:`Combustion chamber <tespy.components.combustion.combustion_chamber>`
    * :py:class:`Combustion chamber stoichiometric <tespy.components.combustion.combustion_chamber_stoich>`
    * :py:class:`Combustion engine <tespy.components.combustion.combustion_engine>`
- Heat exchangers
    * :py:class:`Heat exchanger simple <tespy.components.heat_exchangers.heat_exchanger_simple>`
    * :py:class:`Heat exchanger <tespy.components.heat_exchangers.heat_exchanger>`
    * :py:class:`Condenser <tespy.components.heat_exchangers.condenser>`
    * :py:class:`Desuperheater <tespy.components.heat_exchangers.desuperheater>`
- Nodes
    * :py:class:`Node <tespy.components.nodes.node>`
    * :py:class:`Drum <tespy.components.nodes.drum>`
    * :py:class:`Merge <tespy.components.nodes.merge>`
    * :py:class:`Separator <tespy.components.nodes.separator>`
    * :py:class:`Splitter <tespy.components.nodes.splitter>`
- Piping
    * :py:class:`Pipe <tespy.components.piping.pipe>`
    * :py:class:`Valve <tespy.components.piping.valve>`
- Reactors
    * :py:class:`Water electrolyzer <tespy.components.reactors.water_electrolyzer>`
- Turbomachinery
    * :py:class:`Compressor <tespy.components.turbomachinery.compressor>`
    * :py:class:`Pump <tespy.components.turbomachinery.pump>`
    * :py:class:`Turbine <tespy.components.turbomachinery.turbine>`

.. _using_tespy_components_parametrisation_label:

Component parametrisation
-------------------------

Component parameters can be set and accessed in various ways. All parameters of
components are objects of a :code:`data_container` class. The data container
for component parameters it is called :code:`dc_cp`, :code:`dc_cc` for
component characteristics and :code:`dc_cm` for characteristic maps. The main
purpose of having a data container for the parameters (instead of pure
numbers), is added flexibility for the user. There are different ways for you
to specify a component parameter.

Component parameters
^^^^^^^^^^^^^^^^^^^^

The example shows different ways to specify the heat transfer coefficient of an
evaporator and how to unset the parameter again.

.. code-block:: python

    from tespy.components import heat_exchanger
    from tespy.tools import dc_cp
    import numpy as np

    he = heat_exchanger('evaporator')

    # specify the value
    he.set_attr(kA=1e5)
    # create a data container
    he.set_attr(kA=dc_cp(val=1e5, is_set=True))
    # set data container parameters
    he.kA.set_attr(val=1e5, is_set=True)

    # unset value
    he.set_attr(kA=np.nan)
    he.kA.set_attr(is_set=False)


Grouped parameters
^^^^^^^^^^^^^^^^^^

Grouped parameters are used whenever a component property depends on multiple
parameters. For instance, the pressure loss calculation via Darcy-Weissbach
requires information about the length, diameter and roughness of the pipe.
The solver will prompt a warning, if you do not specify all parameters required
by a parameter group. If parameters of the group are missing, the equation will
not be implemented by the solver.

.. code-block:: python

    from tespy.components import pipe
    import numpy as np

    my_pipe = pipe('pipe')

    # specify grouped parameters
    my_pipe.set_attr(D=0.1, L=20, ks=0.00005)

    # the solver will not use the Darcy-Weissbach-equation in this case
    my_pipe.set_attr(D=0.1, ks=0.00005)

There are three components using parameter groups:

- heat_exchanger_simple and pipe
    * :code:`hydro_group` (:code:`D`, :code:`L`, :code:`ks`)
    * :code:`kA_group` (:code:`kA`, :code:`Tamb`)
- solar_collector
    * :code:`hydro_group` (:code:`D`, :code:`L`, :code:`ks`)
    * :code:`energy_group` (:code:`E`, :code:`eta_opt`, :code:`lkf_lin`,
      :code:`lkf_quad`, :code:`A`, :code:`Tamb`)

Custom variables
^^^^^^^^^^^^^^^^

It is possible to use component parameters as variables of your system of
equations. For example, give a pressure ratio :code:`pr`, length :code:`L` and
roughness :code:`ks` of a pipe you want to calculate the pipe's diameter
:code:`D` required to achieve the specified pressure ratio. In this case you
need to specify the diameter the following way.

.. code-block:: python

    from tespy.components import pipe
    import numpy as np

    # custom variables
    my_pipe = pipe('my pipe')

    # make diameter variable of system
    my_pipe.set_attr(pr=0.98, L=100, ks=0.00002, D='var')

    # a second way of specifying this is similar to the
    # way used in the component parameters section
    # the benefit is, that val will be the starting value
    my_pipe.set_attr(pr=0.98, L=100, ks=0.00002)
    my_pipe.set_attr(D=dc_cp(val=0.2, is_set=True, is_var=True))

It is also possible to set value boundaries for you custom variable. You can do
this, if you expect the result to be within a specific range. But beware: This
might result in a non converging simulation, if the actual value is out of your
specified range.

.. code-block:: python

    # data container specification with identical result,
    # benefit: specification of bounds will increase stability
    my_pipe.set_attr(D=dc_cp(val=0.2, is_set=True, is_var=True,
                             min_val=0.1, max_val=0.3))

Component characteristics
^^^^^^^^^^^^^^^^^^^^^^^^^

Several components integrate parameters using a characteristic function. These
parameters come with default characteristics. As the user might not notice
this, TESPy prompts warnings in the preprocessing, if the default
characteristics are used. The default characteristics available can be found in
the :py:mod:`data <tespy.data>` module. Of course, it is possible to specify
your own characteristic functions.

.. note::

    **There are two different characteristics specifications**

    The characteristic function can be an auxiliary parameter of a different
    component property. This is the case for :code:`kA_char1`
    and :code:`kA_char2` of heat exchangers, :code:`kA_char` of simple
    heat exchangers and pipes as well as the characteristics of a combustion
    engine: :code:`tiP_char`, :code:`Q1_char`, :code:`Q2_char`
    and :code:`Qloss_char`.

    The characteristic function is an individual parameter of the component.
    This is the case for all other components!

    **What does this mean?**

    For the auxiliary functionality the main parameter,
    e. g. :code:`kA` of a heat exchanger must be set :code:`.kA.is_set=True`.

    For the other functionality the characteristics parameter must be
    set e. g. :code:`.eta_s_char.is_set=True`.

For example, :code:`kA` specification for heat exchangers:

.. code-block:: python

    from tespy.components import heat_exchanger
    from tespy.tools import dc_cc
    from tespy.tools.characteristics import load_default_char as ldc
    from tespy.tools.characteristics import char_line
    import numpy as np

    he = heat_exchanger('evaporator', kA=1e5)

    # use a characteristic line from the defaults: specify the component, the
    # parameter and the name of the characteristc function. Also, specify, what
    # type of characteristic function you want to use.
    kA_char1 = ldc('heat exchanger', 'kA_char1', 'EVAPORATING FLUID', char_line)
    kA_char2 = ldc('heat exchanger', 'kA_char2', 'EVAPORATING FLUID', char_line)
    he.set_attr(kA_char1=kA_char1, kA_char2=kA_char2)

    # specification of a data container yields same result. It is aditionally
    # possible to specify the characteristics parameter, mass flow in this case
    # the specification parameters available are stated in the components
    # class documentation
    he.set_attr(kA_char1=dc_cc(param='m', func=kA_char1),
                kA_char2=dc_cc(param='m', func=kA_char2))

    # use custom specification parameters
    x = np.array([0, 0.5, 1, 2])
    y = np.array([0, 0.8, 1, 1.2])
    kA_char1 = char_line(x, y)
    he.set_attr(kA_char1=kA_char1)

For example, :code:`eta_s_char` specification for a pump.

.. code-block:: python

    from tespy.components import pump
    from tespy.tools import dc_cc
    from tespy.tools.characteristics import load_default_char as ldc
    from tespy.tools.characteristics import char_line
    import numpy as np

    pu = pump('pump')

    # use a characteristic line from the defaults
    # CAUTION: this example does only specify the function to follow
    # the given default line. The parameter will not be used in a
    # simulation!
    eta_s_char = ldc('pump', 'eta_s_char', 'DEFAULT', char_line)
    pu.set_attr(eta_s_char=eta_s_char)

    # If we want to use the parameter in the simulation:
    eta_s_char = dc_cc(func=ldc('pump', 'eta_s_char', 'DEFAULT', char_line),
                       is_set=True)
    pu.set_attr(eta_s_char=eta_s_char)

Instead of writing your custom characteristic line information directly into
your Python script, TESPy provides a second method of implementation: It is
possible to store your data in the :code:`HOME/.tespy/data` folder and import
from there. For additional information on formatting and usage, look into
:ref:`this part <using_tespy_characteristics_label>`.

.. code-block:: python

    from tespy.tools.characteristics import load_custom_char as lcc

    eta_s_char = dc_cc(func=lcc('my_custom_char', char_line), is_set=True)
    pu.set_attr(eta_s_char=eta_s_char)


Characteristics are available for the following components and parameters:

- combustion engine
    * :py:meth:`tiP_char <tespy.components.combustion.combustion_engine.tiP_char_func>`: thermal input vs. power ratio.
    * :py:meth:`Q1_char <tespy.components.combustion.combustion_engine.Q1_char_func>`: heat output 1 vs. power ratio.
    * :py:meth:`Q2_char <tespy.components.combustion.combustion_engine.Q2_char_func>`: heat output 2 vs. power ratio.
    * :py:meth:`Qloss_char <tespy.components.combustion.combustion_engine.Qloss_char_func>`: heat loss vs. power ratio.
- compressor
    * :py:meth:`char_map <tespy.components.turbomachinery.compressor.char_map_func>`: component map for isentropic efficiency and pressure rise.
    * :py:meth:`eta_s_char <tespy.components.turbomachinery.compressor.eta_s_char_func>`: isentropic efficiency vs. pressure ratio.
- heat exchangers:
    * :py:meth:`kA1_char, kA2_char <tespy.components.heat_exchangers.heat_exchanger.kA_func>`: heat transfer coefficient vs. mass flow.
- pump
    * :py:meth:`eta_s_char <tespy.components.turbomachinery.pump.eta_s_char_func>`: isentropic efficiency vs. volumetric flow rate.
    * :py:meth:`flow_char <tespy.components.turbomachinery.pump.flow_char_func>`: pressure rise vs. volumetric flow.
- simple heat exchangers
    * :py:meth:`kA_char <tespy.components.heat_exchangers.heat_exchanger_simple.kA_func>`: heat transfer coefficient vs. mass flow.
- turbine
    * :py:meth:`eta_s_char <tespy.components.turbomachinery.turbine.eta_s_char_func>`: isentropic efficiency vs. isentropic enthalpy difference/pressure ratio/volumetric flow/mass flow.
- valve
    * :py:meth:`dp_char <tespy.components.piping.valve.dp_char_func>`: pressure drop vs. flow rate.
- water electrolyzer
    * :py:meth:`eta_char <tespy.components.reactors.water_electrolyzer.eta_char_func>`: efficiency vs. load ratio.

For more information to how the characteristic functions work
:ref:`click here <using_tespy_characteristics_label>`.

Custom components
-----------------

You can add own components. The class should inherit from the
:py:class:`component <tespy.components.components.component>` class or its
children. In order to do that, create a python file in your working directory
and import the base class for your custom component. Now create a class for
your component and at least add the following methods.

- :code:`component(self)`,
- :code:`attr(self)`,
- :code:`inlets(self)`,
- :code:`outlets(self)`,
- :code:`comp_init(self, nw)`,
- :code:`equations(self)`,
- :code:`derivatives(self)` as well as
- :code:`calc_parameters(self)`.

The starting lines of your file should look like this:

.. code:: python

    from tespy.components import component
    from tespy.tools import dc_cc, dc_cp

    class my_custom_component(component):
        """
        This is a custom component.

        You can add your documentation here. From this part, it should be clear
        for the user, which parameters are available, which mandatory equations
        are applied and which optional equations can be applied using the
        component parameters.
        """

        def component(self):
            return 'name of your component'

Attributes
^^^^^^^^^^

The attr method must return a dictionary with the attributes you want to use
for your component. The keys represent the attributes and the respective values
the type of data container used for this attribute. Using the data container
attributes it is possible to add defaults. Defaults for characteristic lines or
characteristic maps are loaded automatically by the component initialisation
method of class :code:`component`. For more information on the default
characteristics consider this
:ref:`chapter <using_tespy_characteristics_label>`.

.. code:: python

    def attr(self):
        return {'par1': dc_cp(min_val=0, max_val=1),
                'par2': dc_cc(param='m')}


Inlets and outlets
^^^^^^^^^^^^^^^^^^

:code:`inlets(self)` and :code:`outlets(self)` respectively must return a list
of strings. The list may look like this:

.. code:: python

    def inlets(self):
        return ['in1', 'in2']

    def outlets(self):
        return ['out1', 'out2']

The number of inlets and outlets might even be generic, e. g. if you have added
an attribute :code:`'num_in'` your code could look like this:

.. code:: python

    def inlets(self):
        if self.num_in.is_set:
            return ['in' + str(i + 1) for i in range(self.num_in.val)]
        else:
            # default number is 2
            return ['in1', 'in2']

Component initialisation
^^^^^^^^^^^^^^^^^^^^^^^^
In the component initialisation you need to determine the total number of
equations and set up the residual value vector as well as the matrix of patial
derivates. The method
:py:meth:`tespy.components.components.component.comp_init` already handles
counting the custom variables and setting up default characteristic lines for
you. The :code:`comp_init()` method of your new component should use call that
method. In order to determine the total number of equations, determine
the number of mandatory equations and the number of optional equations applied.

Then set up the residual value vector and the matrix of partial derivatives.
If the component delivers derivates that are constant, you can paste those
values into the matrix already. The code example shows the implementation of
the :py:meth:`tespy.components.turbomachinery.turbine.comp_init` method.

.. code:: python

    def comp_init(self, nw):

        component.comp_init(self, nw)

        # number of mandatroy equations for
        # fluid balance: num_fl
        # mass flow: 1
        self.num_eq = self.num_nw_fluids + 1

        # number of optional equations, count which parameters are applied
        for var in [self.P, self.pr, self.eta_s, self.eta_s_char, self.cone]:
        if var.is_set is True:
            self.num_eq += 1

        self.mat_deriv = np.zeros((
        self.num_eq,
        self.num_i + self.num_o + self.num_vars,
        self.num_nw_vars))

        self.vec_res = np.zeros(self.num_eq)
        pos = self.num_nw_fluids
        self.mat_deriv[0:pos] = self.fluid_deriv()
        self.mat_deriv[pos:pos + 1] = self.mass_flow_deriv()

Equations
^^^^^^^^^

The equations contain the information on the changes to the fluid properties
within the component. Each equation must be defined in a way, that the correct
result is zero, e. g. if the mass flow at the inlet :math:`\dot{m}_{in}` should
be equal to the mass flow at the outlet :math:`\dot{m}_{out}` and the pressure
at the outlet :math:`p_{out}` is smaller than the pressure at the inlet
:math:`p_{in}` by a specified pressure difference :math:`\Delta p`.

.. math::

    0 = \dot{m}_{in} - \dot{m}_{out}\\
    0 = p_{in} - p_{out} - \Delta p

The connections connected to your component are available as a list in
:code:`self.inl` and :code:`self.outl` respectively. Optional equations should
only be applied, if the paramter has been specified by the user.

.. code:: python

    def equations(self):

        k = 0
        self.vec_res[k] = self.inl[0].m.val_SI - self.outl[0].m.val_SI
        k += 1

        if self.dp.is_set:
            self.vec_res[k] = (
                self.inl[0].p.val_SI - self.outl[0].p.val_SI - self.dp.val
            k += 1

Derivatives
^^^^^^^^^^^

You need to calculate the partial derivatives of the equations to all variables
of the network. This means, that you have to calculate the partial derivatives
to mass flow, pressure, enthalpy and all fluids in the fluid vector on each
incomming or outgoing connection of the component.

Add all derivatives to the matrix (*in the same order as the equations!*).
The derivatives can be calculated analytically or numerically by using the
inbuilt function :code:`numeric_deriv(self, func, dx, pos, **kwargs)`.

- :code:`func` is the function you want to calculate the derivatives for,
- :code:`dx` is the variable you want to calculate the derivative to and
- :code:`pos` indicates the connection you want to calculate the derivative
  for, e. g. :code:`pos=1` means, that counting your inlets and outlets from
  low index to high index (first inlets, then outlets), the connection to be
  used is the second connection in that list.
- :code:`kwargs` are additional keyword arguments required for the function.

For a good start just look into the source code of the inbuilt components. If
you have further questions do not hesitate to contact us. The example above
would look like this:

.. code:: python

    def derivatives(self):

        k = 0
        self.mat_deriv[k, 0, 0] = 1
        self.mat_deriv[k, 1, 0] = -1
        k += 1

        if self.dp.is_set:
            self.mat_deriv[k, 0, 1] = 1
            self.mat_deriv[k, 1, 1] = -1
            k += 1
