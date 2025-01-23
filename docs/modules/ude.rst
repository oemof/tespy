.. _tespy_ude_label:

User defined equations
======================
User defined functions provide a powerful tool to the user as they enable
the definition of generic and individual equations that can be applied to your
TESPy model. In order to implement this functionality in your model you will
use the :py:class:`tespy.tools.helpers.UserDefinedEquation`. The API
documentation provides you with an interesting example application, too.

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
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')

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
    class instance. The **name of the parameter is arbitrary**. We will use
    :code:`ude` in this example. It serves to access some important parameters
    of the equation:

    - connections required in the equation
    - Jacobian matrix to place the partial derivatives
    - automatic numerical derivatives
    - other (external) parameters (e.g. the CharLine in the API docs example of
      :py:class:`tespy.tools.helpers.UserDefinedEquation`)

.. note::

    It is only possible to use the SI-values of the connection variables as
    these values are updated in every iteration. The values in the network's
    specified unit system are only updated after a simulation.

The second step is to define the derivatives with respect to all primary
variables of the network, i.e. mass flow, pressure, enthalpy and fluid
composition of every connection. The derivatives have to be passed to the
Jacobian. In order to do this, we create a function that updates the values
inside the Jacobian of the :code:`UserDefinedEquation` and returns it:

- :code:`ude.jacobian` is a dictionary containing numpy arrays for every
  connection required by the :code:`UserDefinedEquation`.
- derivatives are referred to with the :code:`J_col` attribute of the
  variables of a :code:`Connection` object, i.e.
  :code:`c.m.J_col` for mass flow, :code:`c.p.J_col` for pressure,
  :code:`c.h.J_col` for enthalpy, and :code:`c.fluid.J_col[fluid_name]` for the
  derivative of the fluid composition towards a specific fluid
  :code:`fluid_name`.

If we calculate the derivatives of our equation, it is easy to find, that only
derivatives to mass flow are not zero.

- The derivative to mass flow of connection :code:`c1` is equal to :math:`1`
- The derivative to mass flow of connection :code:`c2` is equal to
  :math:`-2 \cdot \dot{m}_2`.

.. code-block:: python

    >>> def my_ude_deriv(ude):
    ...     c1 = ude.conns[0]
    ...     c2 = ude.conns[1]
    ...     if c1.m.is_var:
    ...         ude.jacobian[c1.m.J_col] = 1
    ...     if c2.m.is_var:
    ...         ude.jacobian[c2.m.J_col] = -2 * ude.conns[1].m.val_SI

Now we can create our instance of the :code:`UserDefinedEquation` and add it to
the network. The class requires four mandatory arguments to be passed:

- :code:`label` of type String.
- :code:`func` which is the function holding the equation to be applied.
- :code:`deriv` which is the function holding the calculation of the Jacobian.
- :code:`conns` which is a list of the connections required by the equation.
  The order of the connections specified in the list is equal to the accessing
  order in the equation and derivative calculation.
- :code:`params` (optional keyword argument) which is a dictionary holding
  additional data required in the equation or derivative calculation.

.. code-block:: python

    >>> ude = UserDefinedEquation('my ude', my_ude, my_ude_deriv, [c1, c2])
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

The derivatives can be determined analytically for the pressure and mass flow
of the first stream easily. For the temperature value, you can use the
predefined fluid property functions :code:`dT_mix_dph` and :code:`dT_mix_pdh`
respectively to calculate the partial derivatives.

.. code-block:: python

    >>> from tespy.tools.fluid_properties import dT_mix_dph
    >>> from tespy.tools.fluid_properties import dT_mix_pdh

    >>> def my_ude_deriv(ude):
    ...     c1 = ude.conns[0]
    ...     c2 = ude.conns[1]
    ...     if c1.m.is_var:
    ...         ude.jacobian[c1.m.J_col] = 1 / ude.conns[0].m.val_SI
    ...     if c1.p.is_var:
    ...         ude.jacobian[c1.p.J_col] = - 2 / ude.conns[0].p.val_SI
    ...     T = c2.calc_T()
    ...     if c2.p.is_var:
    ...         ude.jacobian[c2.p.J_col] = (
    ...             dT_mix_dph(c2.p.val_SI, c2.h.val_SI, c2.fluid_data, c2.mixing_rule)
    ...             * 0.5 / (T ** 0.5)
    ...         )
    ...     if c2.h.is_var:
    ...         ude.jacobian[c2.h.J_col] = (
    ...             dT_mix_pdh(c2.p.val_SI, c2.h.val_SI, c2.fluid_data, c2.mixing_rule)
    ...             * 0.5 / (T ** 0.5)
    ...         )

But, what if the analytical derivative is not available? You can make use of
generic numerical derivatives using the inbuilt method :code:`numeric_deriv`.
The methods expects the variable :code:`'m'`, :code:`'p'`, :code:`'h'` or
:code:`'fluid'` (fluid composition) to derive the function to as well as the
respective connection index from the list of connections. The "lazy" solution
for the above derivatives would therefore look like this:

.. code-block:: python

    >>> def my_ude_deriv(ude):
    ...     c1 = ude.conns[0]
    ...     c2 = ude.conns[1]
    ...     if c1.m.is_var:
    ...         ude.jacobian[c1.m.J_col] = ude.numeric_deriv('m', c1)
    ...     if c1.p.is_var:
    ...         ude.jacobian[c1.p.J_col] = ude.numeric_deriv('p', c1)
    ...     if c2.p.is_var:
    ...         ude.jacobian[c2.p.J_col] = ude.numeric_deriv('p', c2)
    ...     if c2.h.is_var:
    ...         ude.jacobian[c2.h.J_col] = ude.numeric_deriv('h', c2)

    >>> ude = UserDefinedEquation('ude numerical', my_ude, my_ude_deriv, [c1, c2])
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

Obviously, the downside is a slower performance of the solver, as for every
:code:`numeric_deriv` call the function will be evaluated fully twice
(central finite difference).

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

    >>> def my_ude_deriv(ude):
    ...     a = ude.params['a']
    ...     b = ude.params['b']
    ...     c1 = ude.conns[0]
    ...     c2 = ude.conns[1]
    ...     if c1.p.is_var:
    ...         ude.jacobian[c1.p.J_col] = dh_mix_dpQ(c1.p.val_SI, b, c1.fluid_data)
    ...     if c1.h.is_var:
    ...         ude.jacobian[c1.h.J_col] = -a
    ...     if c2.p.is_var:
    ...         ude.jacobian[c2.p.J_col] = a - 1

    >>> ude = UserDefinedEquation(
    ...     'my ude', my_ude, my_ude_deriv, [c1, c2], params={'a': 0.5, 'b': 1}
    ... )


One more example (using a CharLine for data point interpolation) can be found in
the API documentation of class
:py:class:`tespy.tools.helpers.UserDefinedEquation`.

Document your equations
-----------------------

For the automatic documentation of your models just pass the :code:`latex`
keyword on creation of the UserDefinedEquation instance. It should contain the
latex equation string. For example, the last equation from above:

.. code-block:: python

    latex = (
       r'0 = a \cdot \left(h_2 - h_1 \right) - '
       r'\left(h_2 - h\left(p_1, x=b \right)\right)'
    )

    ude = UserDefinedEquation(
       'my ude', my_ude, my_ude_deriv, [c1, c2], params={'a': 0.5, 'b': 1},
       latex={'equation': latex}
    )

The documentation will also create figures of :code:`CharLine` and
:code:`CharMap` objects provided. To add these, adjust the code like this.
Provide the :code:`CharLine` and :code:`CharMap` objects within a list.

.. code-block:: python

    ude = UserDefinedEquation(
       'my ude', my_ude, my_ude_deriv, [c1, c2], params={'a': 0.5, 'b': 1},
       latex={
           'equation': latex,
           'lines': [charline1, charline2],
           'maps': [map1]
       }
    )
