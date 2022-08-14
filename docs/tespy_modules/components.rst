Components
==========

In this section we will introduce you to the details of component
parametrisation and component characteristics. At the end of the section we
show you how to create custom components.

List of components
------------------
More information on the components can be gathered from the code documentation.
We have linked the base class containing a figure and basic information as
well as the equations.

- Basics

    * :py:class:`Cycle closer<tespy.components.basics.cycle_closer.CycleCloser>`
    * :py:class:`Sink <tespy.components.basics.sink.Sink>`
    * :py:class:`Source <tespy.components.basics.source.Source>`
    * :py:class:`Subsystem interface <tespy.components.basics.subsystem_interface.SubsystemInterface>`

- Combustion

    * :py:class:`Combustion chamber <tespy.components.combustion.base.CombustionChamber>`
    * :py:class:`Diabatic combustion chamber <tespy.components.combustion.diabatic.DiabaticCombustionChamber>`
      (Advanced version of combustion chamber, featuring heat losses and
      pressure drop)
    * :py:class:`Combustion engine <tespy.components.combustion.engine.CombustionEngine>`

- Heat exchangers

    * :py:class:`Simplified heat exchanger <tespy.components.heat_exchangers.simple.HeatExchangerSimple>`
    * :py:class:`Solar collector <tespy.components.heat_exchangers.solar_collector.SolarCollector>`
    * :py:class:`Parabolic trough <tespy.components.heat_exchangers.parabolic_trough.ParabolicTrough>`
    * :py:class:`Heat exchanger <tespy.components.heat_exchangers.base.HeatExchanger>`
    * :py:class:`Condenser <tespy.components.heat_exchangers.condenser.Condenser>`
    * :py:class:`Desuperheater <tespy.components.heat_exchangers.desuperheater.Desuperheater>`

- Nodes

    * :py:class:`Droplet separator <tespy.components.nodes.droplet_separator.DropletSeparator>`
    * :py:class:`Drum <tespy.components.nodes.drum.Drum>`
    * :py:class:`Merge <tespy.components.nodes.merge.Merge>`
    * :py:class:`Separator <tespy.components.nodes.separator.Separator>`
    * :py:class:`Splitter <tespy.components.nodes.splitter.Splitter>`

- Piping

    * :py:class:`Pipe <tespy.components.piping.pipe.Pipe>`
    * :py:class:`Valve <tespy.components.piping.valve.Valve>`

- Reactors

    * :py:class:`Fuel cell <tespy.components.reactors.fuel_cell.FuelCell>`
    * :py:class:`Water electrolyzer <tespy.components.reactors.water_electrolyzer.WaterElectrolyzer>`

- Turbomachinery

    * :py:class:`Compressor <tespy.components.turbomachinery.compressor.Compressor>`
    * :py:class:`Pump <tespy.components.turbomachinery.pump.Pump>`
    * :py:class:`Turbine <tespy.components.turbomachinery.turbine.Turbine>`

List of custom components
-------------------------
Here we list the components integrated in the customs module.

- :py:class:`Evaporator for two-phase geothermal organic rankine cycle <tespy.components.customs.orc_evaporator.ORCEvaporator>`

.. _using_tespy_components_parametrisation_label:

Component parametrisation
-------------------------

All parameters of components are objects of a :code:`DataContainer` class. The
data container for component parameters is called
:code:`ComponentProperties`, :code:`ComponentCharacteristics` for component
characteristics, and :code:`ComponentCharacteristicMaps` for characteristic
maps. The main purpose of having a data container for the parameters (instead
of pure numbers), is added flexibility for the user. There are different ways
for you to specify and access component parameters.

Component parameters
^^^^^^^^^^^^^^^^^^^^

The example shows different ways to specify the heat transfer coefficient of an
evaporator and how to unset the parameter again.

.. code-block:: python

    from tespy.components import HeatExchanger
    from tespy.tools import ComponentProperties as dc_cp
    import numpy as np

    he = HeatExchanger('evaporator')

    # specify the value
    he.set_attr(kA=1e5)
    # specify via dictionary
    he.set_attr(kA={'val': 1e5, 'is_set': True})
    # set data container parameters
    he.kA.set_attr(val=1e5, is_set=True)

    # possibilities to unset a value
    he.set_attr(kA=np.nan)
    he.set_attr(kA=None)
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

    from tespy.components import Pipe
    import numpy as np

    my_pipe = Pipe('pipe')

    # specify grouped parameters
    my_pipe.set_attr(D=0.1, L=20, ks=0.00005)

    # the solver will not apply an equation, since the information of the
    # pipe's length is missing.
    my_pipe.set_attr(D=0.1, ks=0.00005)

There are several components using parameter groups:

- heat_exchanger_simple and pipe
    * :code:`hydro_group` (:code:`D`, :code:`L`, :code:`ks`)
    * :code:`kA_group` (:code:`kA`, :code:`Tamb`)
    * :code:`kA_char_group` (:code:`kA_char`, :code:`Tamb`)
- solar_collector
    * :code:`hydro_group` (:code:`D`, :code:`L`, :code:`ks`)
    * :code:`energy_group` (:code:`E`, :code:`eta_opt`, :code:`lkf_lin`,
      :code:`lkf_quad`, :code:`A`, :code:`Tamb`)
- parabolic_trough
    * :code:`hydro_group` (:code:`D`, :code:`L`, :code:`ks`)
    * :code:`energy_group` (:code:`E`, :code:`eta_opt`, :code:`aoi`,
      :code:`doc`, :code:`c_1`, :code:`c_2`, :code:`iam_1`, :code:`iam_2`,
      :code:`A`, :code:`Tamb`)
- compressor
    * :code:`char_map_eta_s_group` (:code:`char_map_eta_s`, :code:`igva`)
    * :code:`char_map_pr_group` (:code:`char_map_pr`, :code:`igva`)

Custom variables
^^^^^^^^^^^^^^^^
It is possible to use component parameters as variables of your system of
equations. In the component parameter list, if a parameter can be a string, it
is possible to specify this parameter as custom variable. For example, given
the pressure ratio :code:`pr`, length :code:`L` and roughness :code:`ks` of a
pipe you may want to calculate the pipe's diameter :code:`D` required to
achieve the specified pressure ratio. In this case you need to specify the
diameter the following way.

.. code-block:: python

    from tespy.components import Pipe
    from tespy.tools import ComponentProperties as dc_cp
    import numpy as np

    # custom variables
    my_pipe = Pipe('my pipe')

    # make diameter variable of system
    my_pipe.set_attr(pr=0.98, L=100, ks=0.00002, D='var')

    # a second way of specifying this is similar to the
    # way used in the component parameters section
    # val will be used as starting value
    my_pipe.set_attr(pr=0.98, L=100, ks=0.00002)
    my_pipe.set_attr(D={'val': 0.2, 'is_set': True, 'is_var': True})

It is also possible to set value boundaries for you custom variable. You can do
this, if you expect the result to be within a specific range. But beware: This
might result in a non converging simulation, if the actual value is out of your
specified range.

.. code-block:: python

    # data container specification with identical result,
    # benefit: specification of bounds will increase stability
    my_pipe.set_attr(D={
        'val': 0.2, 'is_set': True, 'is_var': True,
        'min_val': 0.1, 'max_val': 0.3})

.. _component_characteristic_specification_label:

Component characteristics
^^^^^^^^^^^^^^^^^^^^^^^^^

Several components integrate parameters using a characteristic function. These
parameters come with default characteristics. The default characteristics
available can be found in the :py:mod:`<tespy.data>` module. Of course, it is
possible to specify your own characteristic functions.

.. note::

    **There are two different characteristics specifications**

    The characteristic function can be an auxiliary parameter of a different
    component property. This is the case for :code:`kA_char1`
    and :code:`kA_char2` of heat exchangers as well as the characteristics of a
    combustion engine: :code:`tiP_char`, :code:`Q1_char`, :code:`Q2_char`
    and :code:`Qloss_char`.

    For all other components, the characteristic function is an individual
    parameter of the component.

    **What does this mean?**

    For the auxiliary functionality the main parameter, e.g. :code:`kA_char`
    of a heat exchanger must be set :code:`.kA_char.is_set=True`.

    For the other functionality the characteristics parameter must be
    set e.g. :code:`.eta_s_char.is_set=True`.

For example, :code:`kA_char` specification for heat exchangers:

.. code-block:: python

    from tespy.components import HeatExchanger
    from tespy.tools.characteristics import load_default_char as ldc
    from tespy.tools.characteristics import CharLine
    import numpy as np

    he = HeatExchanger('evaporator')

    # the characteristic function requires the design value of the property,
    # therefore the design value of kA must be set and additionally we set
    # the kA_char method. This is performed automatically, if you specify the
    # kA_char as offdesign parameter (usual case).
    he.set_attr(kA={'design': 1e5}, kA_char={'is_set': True})

    # use a characteristic line from the defaults: specify the component, the
    # parameter and the name of the characteristic function. Also, specify,
    # what type of characteristic function you want to use.
    kA_char1 = ldc('heat exchanger', 'kA_char1', 'DEFAULT', CharLine)
    kA_char2 = ldc('heat exchanger', 'kA_char2', 'EVAPORATING FLUID', CharLine)
    he.set_attr(kA_char2=kA_char2)

    # specification of a data container yields the same result. It is
    # additionally possible to specify the characteristics parameter, e.g. mass
    # flow for kA_char1 and volumetric flow for kA_char2
    he.set_attr(
        kA_char1={'char_func': kA_char1, 'param': 'm'},
        kA_char2={'char_func': kA_char2, 'param': 'v'})

    # or use custom values for the characteristic line e.g. kA vs volumetric
    # flow
    x = np.array([0, 0.5, 1, 2])
    y = np.array([0, 0.8, 1, 1.2])
    kA_char1 = CharLine(x, y)
    he.set_attr(kA_char1={'char_func': kA_char1, 'param': 'v'})

Full working example for :code:`eta_s_char` specification of a turbine.

.. code-block:: python

    from tespy.components import Sink, Source, Turbine
    from tespy.connections import Connection
    from tespy.networks import Network
    from tespy.tools.characteristics import CharLine
    import numpy as np

    fluid_list = ['water']
    nw = Network(fluids=fluid_list, p_unit='bar', T_unit='C', h_unit='kJ / kg')
    si = Sink('sink')
    so = Source('source')
    t = Turbine('turbine')
    inc = Connection(so, 'out1', t, 'in1')
    outg = Connection(t, 'out1', si, 'in1')
    nw.add_conns(inc, outg)

    # design value specification, cone law and eta_s characteristic as
    # offdesign parameters
    eta_s_design = 0.855
    t.set_attr(eta_s=eta_s_design, design=['eta_s'], offdesign=['eta_s_char','cone'])

    # Characteristics x as m/m_design and y as eta_s(m)/eta_s_design
    # make sure to cross the 1/1 point (design point) to yield the same
    # output in the design state of the system
    line = CharLine(
        x=[0.1, 0.3, 0.5, 0.7, 0.9, 1, 1.1],
        y=np.array([0.6, 0.65, 0.75, 0.82, 0.85, 0.855, 0.79]) / eta_s_design)

    # default parameter for x is m / m_design
    t.set_attr(eta_s_char={'char_func': line})
    inc.set_attr(fluid={'water': 1}, m=10, T=550, p=110, design=['p'])
    outg.set_attr(p=0.5)
    nw.solve('design')
    nw.save('tmp')
    # change mass flow value, e.g. 3 kg/s and run offdesign calculation
    inc.set_attr(m=3)
    nw.solve('offdesign', design_path='tmp')
    # isentropic efficiency should be at 0.65
    nw.print_results()

    # alternatively, we can specify the volumetric flow v / v_design for
    # the x lookup
    t.set_attr(eta_s_char={'param': 'v'})
    nw.solve('offdesign', design_path='tmp')
    nw.print_results()

Instead of writing your custom characteristic line information directly into
your Python script, TESPy provides a second method of implementation: It is
possible to store your data in the :code:`HOME/.tespy/data` folder and import
from there. For additional information on formatting and usage, look into
:ref:`this part <using_tespy_characteristics_label>`.

.. code-block:: python

    from tespy.tools.characteristics import load_custom_char as lcc

    eta_s_char = dc_cc(func=lcc('my_custom_char', CharLine), is_set=True)
    t.set_attr(eta_s_char=eta_s_char)

It is possible to allow value extrapolation at the lower and upper limit of the
value range at the creation of characteristic lines. Set the extrapolation
parameter to :code:`True`.

.. code-block:: python

    # use custom specification parameters
    x = np.array([0, 0.5, 1, 2])
    y = np.array([0, 0.8, 1, 1.2])
    kA_char1 = CharLine(x, y, extrapolate=True)
    he.set_attr(kA_char1=kA_char1)

    # set extrapolation to True for existing lines, e.g.
    he.kA_char1.func.extrapolate = True
    pu.eta_s_char.func.extrapolate = True

Characteristics are available for the following components and parameters:

- combustion engine

    * :py:meth:`tiP_char <tespy.components.combustion.engine.CombustionEngine.tiP_char_func>`: thermal input vs. power ratio.
    * :py:meth:`Q1_char <tespy.components.combustion.engine.CombustionEngine.Q1_char_func>`: heat output 1 vs. power ratio.
    * :py:meth:`Q2_char <tespy.components.combustion.engine.CombustionEngine.Q2_char_func>`: heat output 2 vs. power ratio.
    * :py:meth:`Qloss_char <tespy.components.combustion.engine.CombustionEngine.Qloss_char_func>`: heat loss vs. power ratio.

- compressor

    * :py:meth:`char_map <tespy.components.turbomachinery.compressor.Compressor.char_map_func>`: pressure ratio vs. non-dimensional mass flow.
    * :py:meth:`char_map <tespy.components.turbomachinery.compressor.Compressor.char_map_func>`: isentropic efficiency vs. non-dimensional mass flow.
    * :py:meth:`eta_s_char <tespy.components.turbomachinery.compressor.Compressor.eta_s_char_func>`: isentropic efficiency.

- heat exchangers:

    * :py:meth:`kA1_char, kA2_char <tespy.components.heat_exchangers.base.HeatExchanger.kA_char_func>`: heat transfer coefficient.

- pump

    * :py:meth:`eta_s_char <tespy.components.turbomachinery.pump.Pump.eta_s_char_func>`: isentropic efficiency.
    * :py:meth:`flow_char <tespy.components.turbomachinery.pump.Pump.flow_char_func>`: absolute pressure change.

- simple heat exchangers

    * :py:meth:`kA_char <tespy.components.heat_exchangers.simple.HeatExchangerSimple.kA_char_group_func>`: heat transfer coefficient.

- turbine

    * :py:meth:`eta_s_char <tespy.components.turbomachinery.turbine.Turbine.eta_s_char_func>`: isentropic efficiency.

- valve

    * :py:meth:`dp_char <tespy.components.piping.valve.Valve.dp_char_func>`: absolute pressure change.

- water electrolyzer

    * :py:meth:`eta_char <tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.eta_char_func>`: efficiency vs. load ratio.

For more information on how the characteristic functions work
:ref:`click here <using_tespy_characteristics_label>`.

Extend components with new equations
------------------------------------

You can easily add custom equations to the existing components. In order to do
this, you need to implement four changes to the desired component class:

- modify the :code:`get_variables(self)` method.
- add a method, that returns the result of your equation.
- add a method, that places the partial derivatives in the jacobian matrix of
  your component.
- add a method, that returns the LaTeX code of your equation for the automatic
  documentation feature.

In the :code:`get_variables(self)` method, add an entry for your new equation.
If the equation uses a single parameter, use the :code:`ComponentProperties`
type DataContainer (or the :code:`ComponentCharacteristics` type in case you
only apply a characteristic curve). If your equations requires multiple
parameters, add these parameters as :code:`ComponentProperties` or
:code:`ComponentCharacteristics` respectively and add a
:code:`GroupedComponentProperties` type DataContainer holding the information,
e.g. like the :code:`hydro_group` parameter of the
:py:class:`tespy.components.heat_exchangers.simple.HeatExchangerSimple`
class shown below.

.. code:: python

    # [...]
    'D': dc_cp(min_val=1e-2, max_val=2, d=1e-4),
    'L': dc_cp(min_val=1e-1, d=1e-3),
    'ks': dc_cp(val=1e-4, min_val=1e-7, max_val=1e-3, d=1e-8),
    'hydro_group': dc_gcp(
        elements=['L', 'ks', 'D'], num_eq=1,
        latex=self.hydro_group_func_doc,
        func=self.hydro_group_func, deriv=self.hydro_group_deriv),
    # [...]

:code:`latex`, :code:`func` and :code:`deriv` are pointing to the method that
should be applied for the corresponding purpose. For more information on
defining the equations, derivatives and the LaTeX equation you will find the
information in the next section on custom components.


Custom components
-----------------

You can add own components. The class should inherit from the
:py:class:`component <tespy.components.component.Component>` class or its
children. In order to do that, you can use the customs module or create a
python file in your working directory and import the base class for your
custom component. Now create a class for your component and at least add the
following methods.

- :code:`component(self)`,
- :code:`get_variables(self)`,
- :code:`get_mandatory_constraints(self)`,
- :code:`inlets(self)`,
- :code:`outlets(self)` and
- :code:`calc_parameters(self)`.

The starting lines of your file should look like this:

.. code:: python

    from tespy.components.component import Component
    from tespy.tools import ComponentCharacteristics as dc_cc
    from tespy.tools import ComponentProperties as dc_cp

    class my_custom_component(Component):
        """
        This is a custom component.

        You can add your documentation here. From this part, it should be clear
        for the user, which parameters are available, which mandatory equations
        are applied and which optional equations can be applied using the
        component parameters.
        """

        def component(self):
            return 'name of your component'

Mandatory Constraints
^^^^^^^^^^^^^^^^^^^^^

The :code:`get_mandatory_constraints()` method must return a dictionary
containing the information for the mandatory constraints of your component.
The corresponding equations are applied independently of the user
specification. Every key of the mandatory constraints represents one set of
equations. It holds another dictionary with information on

- the equations,
- the number of equations for this constraint,
- the derivatives,
- whether the derivatives are constant values or not (:code:`True/False`) and
- the LaTeX code for the model documentation.

For example, the mandatory equations of a valve look are the following:

.. math::

    0=\dot{m}_{\mathrm{in,1}}-\dot{m}_{\mathrm{out,1}}

    0=x_{fl\mathrm{,in,1}}-x_{fl\mathrm{,out,1}}\;\forall fl
    \in\text{network fluids}

    0=h_{\mathrm{in,1}}-h_{\mathrm{out,1}}

The corresponding method looks like this. The equations, derivatives and
LaTeX string generation are individual methods you need to define
(see next sections).

.. code-block:: python

    def get_mandatory_constraints(self):
        return {
            'mass_flow_constraints': {
                'func': self.mass_flow_func, 'deriv': self.mass_flow_deriv,
                'constant_deriv': True, 'latex': self.mass_flow_func_doc,
                'num_eq': 1},
            'fluid_constraints': {
                'func': self.fluid_func, 'deriv': self.fluid_deriv,
                'constant_deriv': True, 'latex': self.fluid_func_doc,
                'num_eq': self.num_nw_fluids},
            'enthalpy_equality_constraints': {
                'func': self.enthalpy_equality_func,
                'deriv': self.enthalpy_equality_deriv,
                'constant_deriv': True,
                'latex': self.enthalpy_equality_func_doc,
                'num_eq': 1}
        }

- :code:`func`: Method to be applied (returns residual value of equation).
- :code:`deriv`: Partial derivatives of equation to primary variables.
- :code:`latex`: Method returning the LaTeX string of the equation.

Attributes
^^^^^^^^^^

The :code:`get_variables()` method must return a dictionary with the attributes
you want to use for your component. The keys represent the attributes and the
respective values the type of data container used for this attribute. By using
the data container attributes, it is possible to add defaults. Defaults for
characteristic lines or characteristic maps are loaded automatically by the
component initialisation method of class
:py:class:`tespy.components.component.Component`. For more information on the
default characteristics consider this
:ref:`chapter <using_tespy_characteristics_label>`.

The structure is very similar to the mandatory constraints, using
DataContainers instead of dictionaries, e.g. for the Valve:

.. code:: python

    def get_variables(self):
        return {
            'pr': dc_cp(
                min_val=1e-4, max_val=1, num_eq=1,
                deriv=self.pr_deriv, func=self.pr_func,
                func_params={'pr': 'pr'}, latex=self.pr_func_doc),
            'zeta': dc_cp(
                min_val=0, max_val=1e15, num_eq=1,
                deriv=self.zeta_deriv, func=self.zeta_func,
                func_params={'zeta': 'zeta'}, latex=self.zeta_func_doc),
            'dp_char': dc_cc(
                param='m', num_eq=1,
                deriv=self.dp_char_deriv, func=self.dp_char_func,
                char_params={'type': 'abs'}, latex=self.dp_char_func_doc)
        }


Inlets and outlets
^^^^^^^^^^^^^^^^^^

:code:`inlets(self)` and :code:`outlets(self)` respectively must return a list
of strings. The list may look like this:

.. code:: python

    def inlets(self):
        return ['in1', 'in2']

    def outlets(self):
        return ['out1', 'out2']

The number of inlets and outlets might even be generic, e.g. if you have added
an attribute :code:`'num_in'` your code could look like this:

.. code:: python

    def inlets(self):
        if self.num_in.is_set:
            return ['in' + str(i + 1) for i in range(self.num_in.val)]
        else:
            # default number is 2
            return ['in1', 'in2']

Defining equations and derivatives
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Every equation required by the mandatory constraints and in the variables of
the component must be individual methods returning the residual value of the
equation applied. This logic accounts for the derivatives and the LaTeX
equation, too. The Valve's dp_char parameter methods are the following.

.. code:: python

    def dp_char_func(self):
        r"""
        Equation for characteristic line of difference pressure to mass flow.

        Returns
        -------
        residual : ndarray
            Residual value of equation.

            .. math::

                0=p_\mathrm{in}-p_\mathrm{out}-f\left( expr \right)
        """
        p = self.dp_char.param
        expr = self.get_char_expr(p, **self.dp_char.char_params)
        if not expr:
            msg = ('Please choose a valid parameter, you want to link the '
                   'pressure drop to at component ' + self.label + '.')
            logging.error(msg)
            raise ValueError(msg)

        return (
            self.inl[0].p.val_SI - self.outl[0].p.val_SI -
            self.dp_char.char_func.evaluate(expr))

    def dp_char_func_doc(self, label):
        r"""
        Equation for characteristic line of difference pressure to mass flow.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        p = self.dp_char.param
        expr = self.get_char_expr_doc(p, **self.dp_char.char_params)
        if not expr:
            msg = ('Please choose a valid parameter, you want to link the '
                   'pressure drop to at component ' + self.label + '.')
            logging.error(msg)
            raise ValueError(msg)

        latex = (
            r'0=p_\mathrm{in}-p_\mathrm{out}-f\left(' + expr +
            r'\right)')
        return generate_latex_eq(self, latex, label)

    def dp_char_deriv(self, increment_filter, k):
        r"""
        Calculate partial derivatives of difference pressure characteristic.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        if not increment_filter[0, 0]:
            self.jacobian[k, 0, 0] = self.numeric_deriv(
                self.dp_char_func, 'm', 0)
        if self.dp_char.param == 'v':
            self.jacobian[k, 0, 1] = self.numeric_deriv(
                self.dp_char_func, 'p', 0)
            self.jacobian[k, 0, 2] = self.numeric_deriv(
                self.dp_char_func, 'h', 0)
        else:
            self.jacobian[k, 0, 1] = 1

        self.jacobian[k, 1, 1] = -1

For the derivative, the partial derivatives to all variables of the network
are required. This means, that you have to calculate the partial derivatives
to mass flow, pressure, enthalpy and all fluids in the fluid vector on each
connection affecting the equation defined before.

The derivatives can be calculated analytically (preferred if possible) or
numerically by using the inbuilt method
:py:meth:`tespy.components.component.Component.numeric_deriv`, where

- :code:`func` is the function you want to calculate the derivatives for,
- :code:`dx` is the variable you want to calculate the derivative to and
- :code:`pos` indicates the connection you want to calculate the derivative
  for, e.g. :code:`pos=1` means, that counting your inlets and outlets from
  low index to high index (first inlets, then outlets), the connection to be
  used is the second connection in that list.
- :code:`kwargs` are additional keyword arguments required for the function.

LaTeX documentation
^^^^^^^^^^^^^^^^^^^
Finally, add a method that returns the equation as LaTeX string for the
automatic model documentation feature. Simple write the equation and return
it with the :py:meth:`tespy.tools.document_models.generate_latex_eq` method,
which automatically generates a LaTeX equation environment and labels the
equation, so you can reference it later. Therefore, the latex generation
methods needs the label as parameter.

Need assistance?
^^^^^^^^^^^^^^^^
You are very welcome to submit an issue on our GitHub!
