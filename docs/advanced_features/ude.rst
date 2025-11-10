.. _ude_label:

User defined equations
======================
User defined functions provide a powerful tool to the user as they enable
the definition of generic and individual equations that can be applied to your
TESPy model. In order to implement this functionality in your model you will
use the :py:class:`tespy.tools.helpers.UserDefinedEquation`. The API
documentation provides you with an interesting example application, too.

.. attention::

    The API of the :code:`UserDefinedEquation` has changed with version 0.9.
    Please consider the information in the changelog and on this page.

Getting started
---------------

For an easy start, let's consider two different streams. The mass flow of both
streams should be coupled within the model. There is already a possibility
covering simple relations, i.e. applying value referencing with the
:py:class:`tespy.connections.connection.Ref` class. This class allows
formulating simple linear relations:

.. math::

    0 = \dot{m}_1 - \left(\dot{m}_2 \cdot a + b\right)

Instead of this simple application, other relations could be useful. For
example, the mass flow of our first stream should be quadratic to the mass
flow of the second stream.

.. math::

    0 = \dot{m}_1 - \dot{m}_2^2

In order to apply this relation, we need to import the
:py:class:`tespy.tools.helpers.UserDefinedEquation` class into our model and
create an instance with the respective data. First, we set up the TESPy model.

.. code-block:: python

    >>> from tespy.networks import Network
    >>> from tespy.components import Sink, Source
    >>> from tespy.connections import Connection
    >>> from tespy.tools import UserDefinedEquation

    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(
    ...     pressure='bar', temperature='degC', enthalpy='kJ/kg'
    ... )

    >>> so1 = Source('source 1')
    >>> so2 = Source('source 2')
    >>> si1 = Sink('sink 1')
    >>> si2 = Sink('sink 2')

    >>> c1 = Connection(so1, 'out1', si1, 'in1')
    >>> c2 = Connection(so2, 'out1', si2, 'in1')

    >>> nw.add_conns(c1, c2)

    >>> c1.set_attr(fluid={'water': 1}, p=1, T=50)
    >>> c2.set_attr(fluid={'water': 1}, p=5, T=250, v=4)

In the model both streams are well-defined regarding pressure, enthalpy and
fluid composition. The second stream's mass flow is defined through
specification of the volumetric flow, we are missing the mass flow of the
connection :code:`c1`. As described, its value should be quadratic to the
(still unknown) mass flow of :code:`c2`. First, we now need to define the
equation in a function which returns the residual value of the equation.

.. code-block:: python

    >>> def my_ude(ude):
    ...     return ude.conns[0].m.val_SI - ude.conns[1].m.val_SI ** 2


.. note::

    The function must only take one parameter, i.e. the UserDefinedEquation
    class instance, and must be named :code:`ude`! It serves to access some
    important parameters of the equation:

    - connections or components required in the equation
    - automatic numerical derivatives
    - other (external) parameters (e.g. the CharLine in the API docs example of
      :py:class:`tespy.tools.helpers.UserDefinedEquation`)

.. attention::

    It is only possible to use the SI-values of the connection variables as
    these values are updated in every iteration. The values in the network's
    specified unit system are only updated after a simulation.

The second step is to define a function which returns on which variables the
equation depends. This is used to automatically determine the derivatives of
the equation to the system's variables.

.. code-block:: python

    >>> def my_ude_dependents(ude):
    ...     c1, c2 = ude.conns
    ...     return [c1.m, c2.m]

In theory, this is already sufficient information to use the equation in your
model. However, it is possible to additionally provide a function specifying
the derivatives. This is useful if the derivatives can be calculated
analytically. In order to do this, we create a function that updates the values
inside the Jacobian of the :code:`UserDefinedEquation`. We can use the
highlevel method :code:`partial_derivative` for this. In this case the partial
derivatives are easy to find:

- The derivative to mass flow of connection :code:`c1` is equal to :math:`1`
- The derivative to mass flow of connection :code:`c2` is equal to
  :math:`-2 \cdot \dot{m}_2`.

.. code-block:: python

    >>> def my_ude_deriv(increment_filter, k, dependents=None, ude=None):
    ...     c1 = ude.conns[0]
    ...     c2 = ude.conns[1]
    ...     ude.partial_derivative(c1.m, 1)
    ...     ude.partial_derivative(c2.m, -2 * ude.conns[1].m.val_SI)

.. attention::

    The function arguments have to look exactly as provided in the example!

Now we can create our instance of the :code:`UserDefinedEquation` and add it to
the network. The class requires four mandatory arguments to be passed:

- :code:`label` of type String.
- :code:`func` which is the function holding the equation to be applied.
- :code:`dependents` which is the function returning the dependent variables.
- :code:`deriv` (optional) which is the function holding the calculation of the
  Jacobian.
- :code:`conns` (optional) which is a list of the connections required by the
  equation. The order of the connections specified in the list is equal to the
  accessing order in the equation and derivative calculation.
- :code:`comps` (optional) which is a list of the components required by the
  equation. The order of the components specified in the list is equal to the
  accessing order in the equation and derivative calculation.
- :code:`params` (optional) which is a dictionary holding additional data
  required in the equation, dependents specification or derivative calculation.

.. code-block:: python

    >>> ude = UserDefinedEquation(
    ... 'my ude', my_ude, my_ude_dependents,
    ... deriv=my_ude_deriv, conns=[c1, c2]
    ... )
    >>> nw.add_ude(ude)
    >>> nw.solve('design')
    >>> round(c2.m.val_SI ** 2, 2) == round(c1.m.val_SI, 2)
    True
    >>> nw.del_ude(ude)

More examples
-------------

After warm-up let's create some more complex examples, e.g. the square root of
the temperature of the second stream should be equal to the logarithmic value of
the pressure squared divided by the mass flow of the first stream.

.. math::

    0 = \sqrt{T_2} - \ln\left(\frac{p_1^2}{\dot{m}_1}\right)

In order to access the temperature within the iteration process, we need to
calculate it with the respective method. We can import it from the
:py:mod:`tespy.tools.fluid_properties` module. Additionally, import numpy for
the logarithmic value.

.. code-block:: python

    >>> import numpy as np

    >>> def my_ude(ude):
    ...     return (
    ...         ude.conns[1].calc_T() ** 0.5
    ...         - np.log(abs(ude.conns[0].p.val_SI ** 2 / ude.conns[0].m.val_SI))
    ...     )

.. note::

    We use the absolute value inside the logarithm expression to avoid
    ValueErrors within the solution process as the mass flow is not restricted
    to positive values.

.. code-block:: python

    >>> def my_ude_dependents(ude):
    ...     c1 = ude.conns[0]
    ...     c2 = ude.conns[1]
    ...     return [c1.m, c1.p, c2.p, c2.h]

The derivatives can be determined analytically for the pressure and mass flow
of the first stream easily. For the temperature value, you can use the
predefined fluid property functions :code:`dT_mix_dph` and :code:`dT_mix_pdh`
respectively to calculate the partial derivatives.

.. code-block:: python

    >>> from tespy.tools.fluid_properties import dT_mix_dph
    >>> from tespy.tools.fluid_properties import dT_mix_pdh

    >>> def my_ude_deriv(increment_filter, k, dependents=None, ude=None):
    ...     c1 = ude.conns[0]
    ...     c2 = ude.conns[1]
    ...     ude.partial_derivative(c1.m, 1 / ude.conns[0].m.val_SI)
    ...     ude.partial_derivative(c1.p, - 2 / ude.conns[0].p.val_SI)
    ...     T = c2.calc_T()
    ...     # this API also works, it is not as convenient, but saves
    ...     # computational effort because the derivatives are only calculated
    ...     # on demand
    ...     if c2.p.is_var:
    ...         ude.partial_derivative(
    ...             c2.p,
    ...             dT_mix_dph(c2.p.val_SI, c2.h.val_SI, c2.fluid_data, c2.mixing_rule)
    ...             * 0.5 / (T ** 0.5)
    ...         )
    ...     if c2.h.is_var:
    ...         ude.partial_derivative(
    ...             c2.h,
    ...             dT_mix_pdh(c2.p.val_SI, c2.h.val_SI, c2.fluid_data, c2.mixing_rule)
    ...             * 0.5 / (T ** 0.5)
    ...         )

    >>> ude = UserDefinedEquation(
    ...     'ude numerical', my_ude, my_ude_dependents,
    ...     deriv=my_ude_deriv, conns=[c1, c2]
    ... )
    >>> nw.add_ude(ude)
    >>> nw.set_attr(m_range=[.1, 100])  # stabilize algorithm
    >>> nw.solve('design')
    >>> round(c1.m.val, 2)
    1.17

    >>> c1.set_attr(p=None, m=1)
    >>> nw.solve('design')
    >>> round(c1.p.val, 3)
    0.926

    >>> c1.set_attr(p=1)
    >>> c2.set_attr(T=None)
    >>> nw.solve('design')
    >>> round(c2.T.val, 1)
    257.0

But, what if the analytical derivative is not available? Then we can just
not specify the :code:`deriv` keyword to the :code:`UserDefinedEquation`:

.. code-block:: python

    >>> nw.del_ude(ude)
    >>> ude = UserDefinedEquation(
    ...     'ude numerical', my_ude, my_ude_dependents, conns=[c1, c2]
    ... )
    >>> nw.add_ude(ude)
    >>> c1.set_attr(p=None)
    >>> c2.set_attr(T=250)
    >>> nw.solve('design')
    >>> round(c1.p.val, 3)
    0.926
    >>> c1.set_attr(p=1)
    >>> c2.set_attr(T=None)
    >>> nw.solve('design')
    >>> round(c2.T.val, 1)
    257.0

Obviously, the downside is a slower performance of the solver, as for every
dependent the function will be evaluated fully twice (central finite
difference).

Last, we want to consider an example using additional parameters in the
UserDefinedEquation, where :math:`a` might be a factor between 0 and 1 and
:math:`b` is the steam mass fraction (also, between 0 and 1). The difference of
the enthalpy between the two streams multiplied with factor a should be equal
to the difference of the enthalpy of stream two and the enthalpy of saturated
gas at the pressure of stream 1. The definition of the UserDefinedEquation
instance must therefore be changed as below.

.. math::

    0 = a \cdot \left(h_2 - h_1 \right) -
    \left(h_2 - h\left(p_1, x=b \right)\right)

.. code-block:: python

    >>> from tespy.tools.fluid_properties import h_mix_pQ
    >>> from tespy.tools.fluid_properties import dh_mix_dpQ

    >>> def my_ude(ude):
    ...     a = ude.params['a']
    ...     b = ude.params['b']
    ...     c1 = ude.conns[0]
    ...     c2 = ude.conns[1]
    ...     return (
    ...         a * (c2.h.val_SI - c1.h.val_SI) -
    ...         (c2.h.val_SI - h_mix_pQ(c1.p.val_SI, b, c1.fluid_data))
    ...     )

    >>> def my_ude_dependents(ude):
    ...     c1 = ude.conns[0]
    ...     c2 = ude.conns[1]
    ...     return [c1.p, c1.h, c2.h]

    >>> def my_ude_deriv(ude):
    ...     a = ude.params['a']
    ...     b = ude.params['b']
    ...     c1 = ude.conns[0]
    ...     c2 = ude.conns[1]
    ...     ude.partial_derivative(c1.p, dh_mix_dpQ(c1.p.val_SI, b, c1.fluid_data))
    ...     ude.partial_derivative(c1.h, -a)
    ...     ude.partial_derivative(c2.p, a - 1)

    >>> ude = UserDefinedEquation(
    ...     'my ude', my_ude, my_ude_dependents,
    ...     deriv=my_ude_deriv, conns=[c1, c2], params={'a': 0.5, 'b': 1}
    ... )

One more example (using a CharLine for data point interpolation) can be found
in the API documentation of class
:py:class:`tespy.tools.helpers.UserDefinedEquation`.
