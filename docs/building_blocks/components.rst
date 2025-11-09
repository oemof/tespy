.. _tespy_modules_components_label:

Components
==========

In this section we will introduce you to the details of component
parametrisation and component characteristics. At the end of the section we
show you how to create custom components.

List of components
------------------

Below we have collected the different components available in tespy per module.
In the tabs you can view the available parameters per component and relevant
links to the underlying equations, the class documentation and example.

.. include:: _components_overview.rst

.. _tespy_modules_components_parametrisation_label:

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

    >>> from tespy.components import HeatExchanger
    >>> from tespy.tools import ComponentProperties as dc_cp
    >>> import numpy as np

    >>> he = HeatExchanger('evaporator')

    >>> # specify the value
    >>> he.set_attr(kA=1e5)
    >>> # specify via dictionary
    >>> he.set_attr(kA={'_val': 1e5, 'is_set': True})
    >>> # set data container parameters
    >>> he.kA.set_attr(_val=1e5, is_set=True)
    >>> he.kA.is_set
    True

    >>> # possibilities to unset a value
    >>> he.set_attr(kA=np.nan)
    >>> he.set_attr(kA=None)
    >>> he.kA.set_attr(is_set=False)
    >>> he.kA.is_set
    False

Grouped parameters
^^^^^^^^^^^^^^^^^^

Grouped parameters are used whenever a component property depends on multiple
parameters. For instance, the pressure loss calculation via Darcy-Weissbach
requires information about the length, diameter and roughness of the pipe.
The solver will prompt a warning, if you do not specify all parameters required
by a parameter group. If parameters of the group are missing, the equation will
not be implemented by the solver.

.. code-block:: python

    >>> from tespy.components import Pipe, Source, Sink
    >>> from tespy.networks import Network
    >>> from tespy.connections import Connection

    >>> nw = Network()
    >>> nw.units.set_defaults(temperature="degC", pressure="bar")

    >>> so = Source('source')
    >>> si = Sink('sink')
    >>> my_pipe = Pipe('pipe')

    >>> c1 = Connection(so, 'out1', my_pipe, 'in1')
    >>> c2 = Connection(my_pipe, 'out1', si, 'in1')
    >>> nw.add_conns(c1, c2)
    >>> c1.set_attr(fluid={"CH4": 1}, m=1, p=10, T=25)
    >>> c2.set_attr(p0=10, T=25)

    >>> # specify grouped parameters
    >>> my_pipe.set_attr(D=0.1, L=20, ks=0.00005)
    >>> nw.solve('design', init_only=True)
    >>> my_pipe.darcy_group.is_set
    True

    >>> # the solver will not apply an equation, since the information of the
    >>> # pipe's length is now missing (by removing it as follows).
    >>> c2.set_attr(p=10)
    >>> my_pipe.set_attr(L=None)
    >>> nw.solve('design', init_only=True)
    >>> my_pipe.darcy_group.is_set
    False

There are several components using parameter groups:

- heat_exchanger_simple and pipe

  * :code:`darcy_group` (:code:`D`, :code:`L`, :code:`ks`)
  * :code:`hw_group` (:code:`D`, :code:`L`, :code:`ks_HW`)
  * :code:`kA_group` (:code:`kA`, :code:`Tamb`)
  * :code:`kA_char_group` (:code:`kA_char`, :code:`Tamb`)

- solar_collector

  * :code:`darcy_group` (:code:`D`, :code:`L`, :code:`ks`)
  * :code:`hw_group` (:code:`D`, :code:`L`, :code:`ks_HW`)
  * :code:`energy_group` (:code:`E`, :code:`eta_opt`, :code:`lkf_lin`,
    :code:`lkf_quad`, :code:`A`, :code:`Tamb`)

- parabolic_trough

  * :code:`darcy_group` (:code:`D`, :code:`L`, :code:`ks`)
  * :code:`hw_group` (:code:`D`, :code:`L`, :code:`ks_HW`)
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

    >>> # make diameter variable of system
    >>> my_pipe.set_attr(pr=0.98, L=100, ks=0.00002, D='var')
    >>> c2.set_attr(p=None)
    >>> nw.solve("design", init_only=True)
    >>> my_pipe.darcy_group.is_set
    True

    >>> # a second way of specifying this is similar to the
    >>> # way used in the component parameters section
    >>> # val will be used as starting value
    >>> my_pipe.darcy_group.is_set = False
    >>> my_pipe.set_attr(pr=0.98, L=100, ks=0.00002)
    >>> my_pipe.set_attr(D={'_val': 0.2, 'is_set': True, '_is_var': True})
    >>> nw.solve("design", init_only=True)
    >>> my_pipe.darcy_group.is_set
    True

It is also possible to set value boundaries for you custom variable. You can do
this, if you expect the result to be within a specific range. But beware: This
might result in a non converging simulation, if the actual value is out of your
specified range.

.. code-block:: python

    >>> # data container specification with identical result,
    >>> # benefit: specification of bounds will increase stability
    >>> my_pipe.set_attr(D={
    ...     '_val': 0.2, 'is_set': True, '_is_var': True,
    ...     'min_val': 0.1, 'max_val': 0.3}
    ... )
    >>> round(my_pipe.D.max_val, 1)
    0.3

.. _component_characteristic_specification_label:

Component characteristics
^^^^^^^^^^^^^^^^^^^^^^^^^

Several components integrate parameters using a characteristic function. These
parameters come with default characteristics. The default characteristics
available can be found in the :ref:`tespy_data_label`. Of course, it is
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

    >>> from tespy.components import HeatExchanger
    >>> from tespy.tools.characteristics import load_default_char as ldc
    >>> from tespy.tools.characteristics import CharLine

    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(temperature="degC", pressure="bar")

    >>> he = HeatExchanger('evaporator')
    >>> cond = Source('condensate')
    >>> steam = Sink('steam')
    >>> gas_hot = Source('air inlet')
    >>> gas_cold = Sink('air outlet')

    >>> c1 = Connection(cond, "out1", he, "in2")
    >>> c2 = Connection(he, "out2", steam, "in1")
    >>> c3 = Connection(gas_hot, "out1", he, "in1")
    >>> c4 = Connection(he, "out1", gas_cold, "in1")

    >>> nw.add_conns(c1, c2, c3, c4)

    >>> c1.set_attr(fluid={'water': 1}, m=10, p=10, x=0)
    >>> c2.set_attr(p=10, x=1)
    >>> c3.set_attr(fluid={'air': 1}, T=250, p=1)
    >>> c4.set_attr(T=200, p=1)

    >>> nw.solve("design")
    >>> nw.save("design_case.json")
    >>> round(he.kA.val)
    503013

    >>> # the characteristic function is made for offdesign calculation.
    >>> he.set_attr(kA_char={'is_set': True})
    >>> c4.set_attr(T=None)
    >>> nw.solve("offdesign", design_path="design_case.json")
    >>> # since we did not change any property, the offdesign case yields the
    >>> # same value as the design kA value
    >>> round(he.kA.val)
    503013

    >>> c1.set_attr(m=9)
    >>> # use a characteristic line from the defaults: specify the component, the
    >>> # parameter and the name of the characteristic function. Also, specify,
    >>> # what type of characteristic function you want to use.
    >>> kA_char1 = ldc('HeatExchanger', 'kA_char1', 'DEFAULT', CharLine)
    >>> kA_char2 = ldc('HeatExchanger', 'kA_char2', 'EVAPORATING FLUID', CharLine)
    >>> he.set_attr(kA_char2=kA_char2)
    >>> nw.solve("offdesign", design_path="design_case.json")
    >>> round(he.kA.val)
    481745

    >>> # specification of a data container yields the same result. It is
    >>> # additionally possible to specify the characteristics parameter, e.g.
    >>> # mass flow for kA_char1 (identical to default case) and volumetric
    >>> # flow for kA_char2
    >>> he.set_attr(
    ...     kA_char1={'char_func': kA_char1, 'param': 'm'},
    ...     kA_char2={'char_func': kA_char2, 'param': 'v'}
    ... )
    >>> nw.solve("offdesign", design_path="design_case.json")
    >>> round(he.kA.val)
    481745

    >>> # or use custom values for the characteristic line e.g. kA vs volumetric
    >>> # flow
    >>> x = np.array([0, 0.5, 1, 2])
    >>> y = np.array([0, 0.8, 1, 1.2])
    >>> kA_char2 = CharLine(x, y)
    >>> he.set_attr(kA_char2={'char_func': kA_char2, 'param': 'v'})
    >>> nw.solve("offdesign", design_path="design_case.json")
    >>> round(he.kA.val)
    475107

Full working example for :code:`eta_s_char` specification of a turbine.

.. code-block:: python

    >>> from tespy.components import Sink, Source, Turbine
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> from tespy.tools.characteristics import CharLine
    >>> import numpy as np

    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(
    ...     temperature="degC", pressure="bar", enthalpy="kJ/kg"
    ... )
    >>> si = Sink('sink')
    >>> so = Source('source')
    >>> t = Turbine('turbine')
    >>> inc = Connection(so, 'out1', t, 'in1')
    >>> outg = Connection(t, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)

    >>> # design value specification, cone law and eta_s characteristic as
    >>> # offdesign parameters
    >>> eta_s_design = 0.855
    >>> t.set_attr(eta_s=eta_s_design, design=['eta_s'], offdesign=['eta_s_char','cone'])

    >>> # Characteristics x as m/m_design and y as eta_s(m)/eta_s_design
    >>> # make sure to cross the 1/1 point (design point) to yield the same
    >>> # output in the design state of the system
    >>> line = CharLine(
    ...     x=[0.1, 0.3, 0.5, 0.7, 0.9, 1, 1.1],
    ...     y=np.array([0.6, 0.65, 0.75, 0.82, 0.85, 0.855, 0.79]) / eta_s_design
    ... )

    >>> # default parameter for x is m / m_design
    >>> t.set_attr(eta_s_char={'char_func': line})
    >>> inc.set_attr(fluid={'water': 1}, m=10, T=550, p=110, design=['p'])
    >>> outg.set_attr(p=0.5)
    >>> nw.solve('design')
    >>> nw.save('tmp.json')
    >>> # change mass flow value, e.g. 3 kg/s and run offdesign calculation
    >>> inc.set_attr(m=3)
    >>> nw.solve('offdesign', design_path='tmp.json')
    >>> # isentropic efficiency should be at 0.65
    >>> round(t.eta_s.val, 2)
    0.65

    >>> # alternatively, we can specify the volumetric flow v / v_design for
    >>> # the x lookup
    >>> t.set_attr(eta_s_char={'param': 'v'})
    >>> nw.solve('offdesign', design_path='tmp.json')
    >>> round(t.eta_s.val, 2)
    0.84

Instead of writing your custom characteristic line information directly into
your Python script, TESPy provides a second method of implementation: It is
possible to store your data in the :code:`HOME/.tespy/data` folder and import
from there. For additional information on formatting and usage, look into
:ref:`this part <tespy_modules_characteristics_label>`.

.. code-block:: python

    from tespy.tools.characteristics import load_custom_char as lcc

    eta_s_char = dc_cc(func=lcc('my_custom_char', CharLine), is_set=True)
    t.set_attr(eta_s_char=eta_s_char)

It is possible to allow value extrapolation at the lower and upper limit of the
value range at the creation of characteristic lines. Set the extrapolation
parameter to :code:`True`.

.. code-block:: python

    # use custom specification parameters
    >>> x = np.array([0, 0.5, 1, 2])
    >>> y = np.array([0, 0.8, 1, 1.2])
    >>> kA_char1 = CharLine(x, y, extrapolate=True)
    >>> kA_char1.extrapolate
    True

    >>> # set extrapolation to True for existing lines, e.g.
    >>> he.kA_char1.char_func.extrapolate = True
    >>> he.kA_char1.char_func.extrapolate
    True

Characteristics are available for the following components and parameters:

- CombustionEngine

  * :py:meth:`tiP_char <tespy.components.combustion.engine.CombustionEngine.tiP_char_func>`: thermal input vs. power ratio.
  * :py:meth:`Q1_char <tespy.components.combustion.engine.CombustionEngine.Q1_char_func>`: heat output 1 vs. power ratio.
  * :py:meth:`Q2_char <tespy.components.combustion.engine.CombustionEngine.Q2_char_func>`: heat output 2 vs. power ratio.
  * :py:meth:`Qloss_char <tespy.components.combustion.engine.CombustionEngine.Qloss_char_func>`: heat loss vs. power ratio.

- Compressor

  * :py:meth:`char_map <tespy.components.turbomachinery.compressor.Compressor.char_map_pr_func>`: pressure ratio vs. non-dimensional mass flow.
  * :py:meth:`char_map <tespy.components.turbomachinery.compressor.Compressor.char_map_eta_s_func>`: isentropic efficiency vs. non-dimensional mass flow.
  * :py:meth:`eta_s_char <tespy.components.turbomachinery.compressor.Compressor.eta_s_char_func>`: isentropic efficiency.

- PolynomialCompressor

  * :py:meth:`<tespy.components.displacementmachinery.polynomial_compressor.PolynomialCompressor.eta_s_poly_group_func>`: isentropic efficiency based on EN12900 polynomial
  * :py:meth:`<tespy.components.displacementmachinery.polynomial_compressor.PolynomialCompressor.eta_vol_poly_group_func>`: volumetric efficiency based on EN12900 polynomial

- HeatExchanger:

  * :py:meth:`kA1_char, kA2_char <tespy.components.heat_exchangers.base.HeatExchanger.kA_char_func>`: heat transfer coefficient.

- Pump

  * :py:meth:`eta_s_char <tespy.components.turbomachinery.pump.Pump.eta_s_char_func>`: isentropic efficiency.
  * :py:meth:`flow_char <tespy.components.turbomachinery.pump.Pump.flow_char_func>`: absolute pressure change.

- SimpleHeatExchanger

  * :py:meth:`kA_char <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.kA_char_group_func>`: heat transfer coefficient.

- Turbine

  * :py:meth:`eta_s_char <tespy.components.turbomachinery.turbine.Turbine.eta_s_char_func>`: isentropic efficiency.

- Valve

  * :py:meth:`dp_char <tespy.components.piping.valve.Valve.dp_char_func>`: absolute pressure change.

- WaterElectrolyzer

  * :py:meth:`eta_char <tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.eta_char_func>`: efficiency vs. load ratio.

For more information on how the characteristic functions work
:ref:`click here <tespy_modules_characteristics_label>`.

Extend components with new equations
------------------------------------

You can easily add custom equations to the existing components. In order to do
this, you need to implement four changes to the desired component class:

- modify the :code:`get_parameters(self)` method.
- add a method, that returns the result of your equation.
- add a method, that returns the variables your equation depends on.

In the :code:`get_parameters(self)` method, add an entry for your new equation.
If the equation uses a single parameter, use the :code:`ComponentProperties`
type DataContainer (or the :code:`ComponentCharacteristics` type in case you
only apply a characteristic curve). If your equations requires multiple
parameters, add these parameters as :code:`ComponentProperties` or
:code:`ComponentCharacteristics` respectively and add a
:code:`GroupedComponentProperties` type DataContainer holding the information,
e.g. like the :code:`darcy_group` parameter of the
:py:class:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger`
class shown below.

.. code:: python

    # [...]
    'D': dc_cp(min_val=1e-2, max_val=2, d=1e-4, quantity="length"),
    'L': dc_cp(min_val=1e-1, d=1e-3, quantity="length"),
    'ks': dc_cp(val=1e-4, min_val=1e-7, max_val=1e-3, d=1e-8, quantity="length"),
    'darcy_group': dc_gcp(
        elements=['L', 'ks', 'D'], num_eq_sets=1,
        func=self.darcy_func,
        dependents=self.darcy_dependents
    ),
    # [...]

.. tip::

    With the :code:`quantity` keyword, tespy will automatically understand what
    unit conversion to apply to the respective parameter. E.g. in case you
    want to specify the roughness  :code:`ks` in millimeter, you can either
    set the default unit for length of your :code:`Network` to millimeter, or
    you can pass the :code:`ks` value as :code:`pint.Quantity` to your
    component using millimeter as unit. Then the conversion to the SI unit is
    taken care of automatically in the preprocessing and the respective
    equation will make use of the SI value.

:code:`func` and :code:`dependents` are pointing to the method that should be
applied for the corresponding purpose. For more information on defining the
equations and dependents see the next section on custom components. When
defining the dependents in a standalone way, the partial derivatives are
calculated automatically. If you want to insert the partial derivatives
manually, you can define another function and pass with the :code:`deriv`
keyword.
