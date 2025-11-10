.. _components_custom_components_label:

Custom components
-----------------

You can add own components. The class should inherit from the
:py:class:`component <tespy.components.component.Component>` class or its
children. In order to do that, you can use the customs module or create a
python file in your working directory and import the base class for your
custom component. Now create a class for your component and at least add the
following methods.

- :code:`component(self)`,
- :code:`get_parameters(self)`,
- :code:`get_mandatory_constraints(self)`,
- :code:`inlets(self)`,
- :code:`outlets(self)` and
- :code:`calc_parameters(self)`.

Optionally, you can add

- :code:`powerinlets(self)` and
- :code:`poweroutlets(self)`

in case your component should have methods to connect the material flows with
non-material flows associated with a :code:`PowerConnection`.

.. note::

  For more information on the :code:`PowerConnection` please check the
  respective :ref:`section in the docs <powerconnections_label>`.

The starting lines of your file should look like this:

.. code-block:: python

    from tespy.components.component import Component
    from tespy.tools import ComponentCharacteristics as dc_cc
    from tespy.tools import ComponentMandatoryConstraints as dc_cmc
    from tespy.tools import ComponentProperties as dc_cp

    class MyCustomComponent(Component):
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
- the number of equations for this constraint and
- the variables each equation depends on.

Furthermore more optional specifications can be made

- the partial derivatives,
- whether the derivatives are constant values or not (:code:`True/False`) and
- the structure_matrix keyword.

For example, the mandatory equations of the class :code:`Valve` look are the
following:

.. math::

    0=h_{\mathrm{in,1}}-h_{\mathrm{out,1}}

The corresponding method looks like this:

.. literalinclude:: /../src/tespy/components/piping/valve.py
    :pyobject: Valve.get_mandatory_constraints

The method inherits from the :code:`Component` base class and then adds the
enthalpy equality constraint on top of the mass flow equality and fluid
equality constraints.

.. note::

    In this simple case only the :code:`structure_matrix` has to be provided.
    It creates a mapping between linearly dependent pairs of variables and is
    utilized to simplify the problem during presolving. It is generally
    optional.

For equations, that depend on more than two variables, or that do not have
direct linear relationsships additional parameters have to be supplied, e.g.
see the respective method of the class :code:`HeatExchanger`.

.. literalinclude:: /../src/tespy/components/heat_exchangers/base.py
    :pyobject: HeatExchanger.get_mandatory_constraints

Here we have the following keywords:

- :code:`func`: Method to be applied (returns residual value of equation)
- :code:`dependents`: Method to return the variables :code:`func` depends on
- :code:`num_eq_sets`: Number of equations

.. note::

    In some cases the number of equations can depend on the length of the fluid
    vectors associated with the component. :code:`num_eq_sets` specifically
    points to the number of equation per all fluids in the fluid vector. Since
    this number is not necessarily known prior to solving the problem, there is
    a possibility to update the number of equations after presolving to
    determine the correct number. This update is only relevant for classes like
    :code:`Merge` and :code:`CombustionChamber` etc.. Feel free to reach out in
    the discussion forum, if you have any questions about it.

With the above mentioned specifications, tespy will apply the method to
calculate the residual value of your equation and automatically calculates its
partial derivatives towards all variables specified in the :code:`dependents`
list.

Finally, sometimes it is reasonable to not let tespy automatically calculate
all partial derivatives, because the calculation can be computationally
expensive. Instead you can additionally provide the following keyword:

- :code:`deriv`: A method that calculate the partial derivatives.

You will find more information and examples on this in the next sections.

You can also define mandatory constraints that are conditional, e.g. in context
of :code:`PowerConnections`. For example, the connection between the material
flow variables of the inlet and the outlet of the turbine to the non-material
energy output variable of the turbine should only be made, in case the turbine
is actually connected with a :code:`PowerConnection`:

.. literalinclude:: /../src/tespy/components/turbomachinery/turbine.py
    :pyobject: Turbine.get_mandatory_constraints

Attributes
^^^^^^^^^^

This part is very similar to the previous one. The :code:`get_parameters()`
method must return a dictionary with the attributes you want to use for your
component. The keys represent the attributes and the respective values the type
of data container used for this attribute. By using the data container
attributes, it is possible to add defaults. Defaults for characteristic lines
or characteristic maps are loaded automatically by the component initialisation
method of class
:py:class:`tespy.components.component.Component`. For more information on the
default characteristics consider this
:ref:`chapter <modules_characteristics_label>`.

The structure is very similar to the mandatory constraints, e.g. for the
class :code:`Valve`:

.. literalinclude:: /../src/tespy/components/piping/valve.py
    :pyobject: Valve.get_parameters

Inlets and outlets
^^^^^^^^^^^^^^^^^^

:code:`inlets(self)` and :code:`outlets(self)` respectively must return a list
of strings. The list may look like this (of class :code:`HeatExchanger`)

.. literalinclude:: /../src/tespy/components/heat_exchangers/base.py
    :pyobject: HeatExchanger.inlets

.. literalinclude:: /../src/tespy/components/heat_exchangers/base.py
    :pyobject: HeatExchanger.outlets

The number of inlets and outlets might even be variable, e.g. if you have added
an attribute :code:`'num_in'` your code could look like this (as in class
:code:`Merge`):

.. literalinclude:: /../src/tespy/components/nodes/merge.py
    :pyobject: Merge.inlets

Inlets and outlets for PowerConnections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If your component should incorporate :code:`PowerConnections` you can define
connctor ids in a similar way, for example power inlet for compressors or
power outlet for turbines. Here the methods are :code:`powerinlets` and
:code:`poweroutlets`.

.. literalinclude:: /../src/tespy/components/turbomachinery/compressor.py
    :pyobject: Compressor.powerinlets

.. literalinclude:: /../src/tespy/components/turbomachinery/turbine.py
    :pyobject: Turbine.poweroutlets

In a similar way, you can add flexibility with a dynamic number of inlets and
outlets:

.. literalinclude:: /../src/tespy/components/power/bus.py
    :pyobject: PowerBus.powerinlets

.. literalinclude:: /../src/tespy/components/power/bus.py
    :pyobject: PowerBus.poweroutlets

Define the required methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^
In the above section the concept of the component mandatory constraints and
their attributes was introduced. Now we need to fill the respective parts
with some life, i.e. how to define

- the :code:`structure_matrix` (optional),
- the :code:`func`,
- the :code:`dependents` and
- the :code:`deriv` (optional) methods.

Define the structure matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^
As mentioned, with the structure matrix you can make a mapping, in case two
variables are linked to each other with a linear relationship. The presolving
of a model will utilize this information to reduce the number of variables.
For example, for a specified pressure ratio :code:`pr` of a component, where
the inlet and the outlet pressure are linked through this equation:

.. math::

    p_\text{inlet} \cdot \text{pr} - p_\text{outlet} = 0

We can create a method and reference to it from the component mandatory
constraints or attribute dictionaries. In this method you have to

- place the partial derivatives towards both variables in the component's
  :code:`_structure_matrix` attribute.
- place any offset in the component's :code:`_rhs` attribute.

For the example above, the derivative to the inlet pressure is :code:`pr`, and
to the outlet pressure :code:`-1`. The offset/right hand side value of the
equation is 0.

.. literalinclude:: /../src/tespy/components/component.py
    :pyobject: Component.pr_structure_matrix

A different equation to simplify with this method could be the delta pressure
:code:`dp`. In this case, the :code:`_rhs` is not zero, it is the value of
:code:`dp`.

.. math::

    p_\text{inlet} - p_\text{outlet} = \text{dp}

.. literalinclude:: /../src/tespy/components/component.py
    :pyobject: Component.dp_structure_matrix

Define the equations
^^^^^^^^^^^^^^^^^^^^
The definition of an equation is quite straight forward: It must return its
residual value. For example, the equation of the :code:`dp_char` parameter
associated with the class :code:`Valve` is the following:

.. math::

    p_\text{inlet} - p_\text{outlet} - f_\text{dp}\left(x\right) = 0

.. literalinclude:: /../src/tespy/components/piping/valve.py
    :pyobject: Valve.dp_char_func

.. caution::

    Your equations should only use and access the SI values :code:`val_SI`
    associated with connection or component parameters in the back-end .

Define the dependents
^^^^^^^^^^^^^^^^^^^^^
Next, you have to define the list of variables the equation depends on, i.e.
towards which variables the partial derivatives should be calculated. In this
example, it is the inlet and the outlet pressure, as well as the mass flow and
in case the volumetric flow should be used to assess the characteristic
function, the inlet enthalpy.

.. literalinclude:: /../src/tespy/components/piping/valve.py
    :pyobject: Valve.dp_char_dependents

The solver will automatically determine, which of the variables returned by
this method are actual variables (have not been presolved) and the calculate
the derivative to the specified equation numerically using a central finite
difference. In the case of this method, this will be an extra 6 or 8 function
evaluations to determine the partial derivatives, if all of the indicated
variables are actually system variables (have not been presolved).

The only thing you have to do is, to make the method return a list of variables
the equation depends on.

It can be more complex than that when dealing with equations, which have
partial derivatives towards components of a fluid mixture. For example, the
energy balance of the :code:`CommbustionChamber` depends on the fuel's mass
fraction in the fluid mixtures of its inlets. To account for this in the
dependents specification, your method has to return a dictionary instead,
which uses the keys

- :code:`scalars` for all "standard" variables
- :code:`vectors` for all fluid mixture component variables

In the example below, the variable mixture components of the inlets are the
union of the set of fuels available in the :code:`CombustionChamber` and the
fluid components that are actually variable in the mixture. For this, a
subdictionary is created, which is a mapping of the fluid mixture container
:code:`c.fluid` to a set of fluid names
:code:`self.fuel_list & c.fluid.is_var`.

.. literalinclude:: /../src/tespy/components/combustion/base.py
    :pyobject: CombustionChamber.energy_balance_dependents

Define the derivatives
^^^^^^^^^^^^^^^^^^^^^^
The downside of the simple to use approach of defining the equation together
with its dependents is, that it can be computationally expensive to calculate
the partial derivatives. In this case, it may be reasonable to implement a
method specifically for the calculation of the partial derivatives.

For example, consider the isentropic efficiency equation of a :code:`Turbine`:

.. literalinclude:: /../src/tespy/components/turbomachinery/turbine.py
    :pyobject: Turbine.eta_s_func

The partial derivatives to the inlet and outlet pressure as well as the inlet
enthalpy can only be determined numerically. However, the partial derivative to
the outlet enthalpy can be obtained analytically, it is :code:`1`. To save the
extra evaluation of the equation in case the outlet enthalpy is a variable, we
can define the following method:

.. literalinclude:: /../src/tespy/components/turbomachinery/turbine.py
    :pyobject: Turbine.eta_s_deriv

To place the partial derivative you can use the :code:`_partial_derivative`
method and pass

- the variable
- the equation number (passed to your method through the argument k)
- the value of the partial derivative (a number or a callable)

  - in case you pass a number, it will put the value directly into the
    Jacobian
  - in case you pass a callable, the derivative will be determined numerically
    for the specified callable and the result will then be passed to the
    Jacobian

- the :code:`increment_filter`, which is a lookup for variables, that do not
  change anymore from one iteration to the next. In this case, the calculation
  of the derivative will be skipped.

.. attention::

    We cannot simply put down the derivatives for all variables in the Jacobian
    because we do not necessarily know (prior to solving) which variables will
    be mapped to a single variable because they are linearly dependent. Thus,
    we have to use the set of dependents, that is passed to our derivative
    method. Otherwise, the  calculation of the derivative, e.g. for outlet
    pressure may override the value for inlet pressure, even though both are
    pointing to the same variable. In case of numerical derivative calculation
    this is not an issue except for the extra computational effort. But if you
    have determined the derivatives analytically, then their value might change
    if two variables are mapped to a single one.

Need assistance?
^^^^^^^^^^^^^^^^
You are very welcome to open a discussion or submit an issue on the GitHub
repository!
