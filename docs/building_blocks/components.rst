.. _modules_components_label:

Components
==========

This section provides an overview on all available component classes, how to
specify simple inputs to components, implement custom values in context of
characteristic lines or maps, and give hints as to how to implement custom
equations to existing components or create custom components yourself.

List of components
------------------

Below we have collected the different components available in tespy per module.
In the tabs you can view the available parameters per component and relevant
links to the underlying equations, the class documentation and example.

.. include:: _components_overview.rst

.. _modules_components_parametrisation_label:

Component parametrisation
-------------------------

All parameters of components can be set and unset through the :code:`set_attr`
method. As seen in the tables above, there are different types of parameters:

- Simple component parameters: If an equation/method is associated with the
  parameter, setting a value for the parameter will activate that equation.
- Grouped component parameters: Multiple simple parameters can be part of a
  parameter group. To activate the equation associated with the parameter
  group, values for all elements of the group have to be set by the user.
- Component parameters that serve as variable in the model: Specifically in the
  context of grouped parameters it is possible to set one of the parameters as
  a variable. If values for all other parameters of the group are specified as
  well, then the group equation is activated and the solver will include the
  one parameter marked as variable in the variables.
- Component characteristics and characteristic groups: Lookup-tables are
  implemented here to interpolate over a variable or an expression and modify
  the outcome of an equation.

Component parameters
^^^^^^^^^^^^^^^^^^^^

The example shows how to specify a component parameter and how to unset again
at the example of the heat transfer coefficient of an evaporator.

.. code-block:: python

    >>> from tespy.components import HeatExchanger
    >>> from tespy.tools import ComponentProperties as dc_cp
    >>> import numpy as np

    >>> he = HeatExchanger('evaporator')

    >>> # specify the value
    >>> he.set_attr(kA=1e5)
    >>> he.kA.val
    100000.0
    >>> he.kA.is_set
    True

    >>> # possibilities to unset a value
    >>> he.set_attr(kA=None)
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

.. _component_variables_label:

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

It is also possible to set value boundaries for you custom variable. You can do
this, if you expect the result to be within a specific range. But beware: This
might result in a non converging simulation, if the actual value is out of your
specified range.

.. code-block:: python

    >>> my_pipe.D.max_val = 0.3
    >>> round(my_pipe.D.max_val, 1)
    0.3

.. _component_characteristic_specification_label:

Component characteristics
^^^^^^^^^^^^^^^^^^^^^^^^^

Several components integrate parameters using a characteristic function. These
parameters come with default characteristics. The default characteristics
available can be found in the :ref:`this section <data_label>`. Of course, it
is possible to specify your own data for these characteristic functions.

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
    >>> he.set_attr(offdesign=["kA_char"])
    >>> c4.set_attr(design=["T"])
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
    >>> t.set_attr(
    ...     eta_s=eta_s_design,
    ...     design=['eta_s'], offdesign=['eta_s_char','cone']
    ... )

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

Finally, it is possible to allow value extrapolation at the lower and upper
limit of the value range at the creation of characteristic lines. Set the
extrapolation parameter to :code:`True`.

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

For more information on how the characteristic functions work
:ref:`click here <modules_characteristics_label>`.

Component customization
-----------------------

Please check the :ref:`advanced features <custom_components_label>`
on how to

- extend existing components with additional parameters and equations or
- create your own components from scratch
