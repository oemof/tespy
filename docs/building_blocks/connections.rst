.. _modules_connections_label:

Connections
===========

This section provides an overview of the available :code:`Connection` classes
in the tabs below. Beyond that, it gives an introduction on how to
parametrize instances. Connections hold the variables that are solved for in
the system of equations of all your models.

.. include:: _connections_overview.rst

Connection Overview
-------------------

The tables above indicate, which specification parameters are available for
each class. To set/unset values the same logic applies as is used in
components.

Setting and unsetting values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We will create a simple problem using a :code:`SimpleHeatExchanger` to showcase
specification options. Further down there is a second simple problem
showcasing the :code:`PowerConnection`.

A :code:`Connection` always connects two components:

.. code-block:: python

    >>> from tespy.networks import Network
    >>> from tespy.connections import Connection
    >>> from tespy.components import Sink, Source, SimpleHeatExchanger
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(
    ...     temperature="°C", power="kW", pressure="bar", enthalpy="kJ/kg"
    ... )
    >>> source = Source('source')
    >>> heatexchanger = SimpleHeatExchanger('heat exchanger')
    >>> sink = Sink('sink')
    >>> c1 = Connection(source, "out1", heatexchanger, "in1", label="c1")
    >>> c2 = Connection(heatexchanger, "out1", sink, "in1", label="c2")
    >>> nw.add_conns(c1, c2)

It is possible to set simple values and then solve, e.g.:

- mass flow, pressure and temperature (and fluid) at inlet
- enthalpy at outlet
- (pressure drop in heat exchanger)

.. code-block:: python

    >>> c1.set_attr(fluid={"R290": 1}, m=5, p=3, T=50)
    >>> c2.set_attr(h=500)
    >>> heatexchanger.set_attr(dp=0)
    >>> nw.solve("design")

Both connections will have all results available, these can be accessed

- as SI value
- as value in the network's default unit of the respective quantity
- as value with the corresponding quantity (:code:`pint.Quantity`)

The results include mass flow, pressure, enthalpy, temperature, temperature,
vapor quality, specific volume :code:`vol` and volumetric flow.

.. code-block:: python

    >>> round(c1.h.val, 1)  # value in kJ/kg but without unit attached
    668.9
    >>> round(c1.vol.val_with_unit, 3)  # will be in m3/kg
    <Quantity(0.195, 'm3 / kilogram')>
    >>> round(c2.T.val_SI, 1)  # SI value
    259.0
    >>> round(c2.T.val_with_unit, 1)  # will be in °C
    <Quantity(-14.2, 'degree_Celsius')>

You can also provide quantities to a specific parameter to individually specify
a unit to a parameter, e.g. inlet mass flow. Note, that units are retained when
set with individual quantity.

.. code-block:: python

    >>> Q = nw.units.ureg.Quantity
    >>> c1.set_attr(m=Q(2, "t/h"))
    >>> nw.solve("design")
    >>> c1.m.val_with_unit
    <Quantity(2, 'metric_ton / hour')>
    >>> round(c2.m.val_with_unit, 2)
    <Quantity(0.56, 'kilogram / second')>

For pure fluids or CoolProp/REFPROP mixtures we can also specify two-phase
properties:

- vapor mass fraction/quality :code:`x`
- dew line temperature difference for superheating :code:`td_dew`
- bubble line temperature difference for subcooling :code:`td_bubble`
- dew line temperature :code:`T_dew` and bubble line temperature
  :code:`T_bubble` to impose the corresponding pressure to the model

We can replace the inlet temperature specification e.g. with superheating. The
unit of :code:`temperature_difference` is different from the unit for
:code:`temperature`. Unsetting a value is simple: Just set it to :code:`None`.

.. code-block:: python

    >>> c1.set_attr(T=None)  # unset the value
    >>> nw.units.default["temperature"]
    '°C'
    >>> nw.units.default["temperature_difference"]
    'delta_degC'
    >>> c1.set_attr(td_dew=20)
    >>> nw.solve("design")
    >>> round(c1.T.val, 2)
    5.82

Setting starting values
^^^^^^^^^^^^^^^^^^^^^^^

Setting starting values for the variables can be helpful in some situations.
You can do this for the following properties:

- mass flow :code:`m0`
- pressure :code:`p0`
- enthalpy :code:`h0`

These specifications are optional!

.. code-block:: python

    >>> c1.set_attr(m0=4)
    >>> c2.set_attr(h0=300, p0=4)

Referencing specifications
^^^^^^^^^^^^^^^^^^^^^^^^^^

It is also possible to set up linear relationships between parameters between
different instances of :code:`Connection` in your system in the following form:

.. math::

    x_0 = a * x_1 + b

It is possible to specify these for

- mass flow, pressure and enthalpy as well as
- temperature and volumetric flow.

For example, instead of h we can specify a reference to the temperature at
c1. The factor is always based on SI value, the delta is in the default unit of
the respective property. The starting value for h is required in this context
because in the previous calculation the fluid was in two-phase state, meaning
the partial derivative of the temperature of c2 with respect to enthalpy would
be zero otherwise.

.. code-block:: python

    >>> from tespy.connections import Ref
    >>> factor = 1
    >>> delta = 25
    >>> c2.set_attr(h=None, T=Ref(c1, factor, delta), h0=1000)
    >>> nw.solve("design")
    >>> round(c1.T.val_SI * factor + delta - 273.15, 2)
    30.82
    >>> round(c2.T.val, 2)
    30.82

Instead we could also reference volumetric flow at outlet to find the
temperature at outlet.

.. code-block:: python

    >>> c2.set_attr(T=None)
    >>> factor = 1.2
    >>> delta = 0
    >>> c2.set_attr(v=Ref(c1, factor, delta))
    >>> nw.solve("design")
    >>> round(c1.v.val_SI * factor + delta, 2)
    0.11
    >>> round(c2.v.val_SI, 2)
    0.11
    >>> round(c2.T.val, 2)
    52.91

For more complex (and arbitrary) relationships between variables of the system
use the :code:`UserDefinedEquation` class. Some examples can be found in
:ref:`this section <ude_label>`.

Fluid specification
^^^^^^^^^^^^^^^^^^^

This sections shows some details on the specification of fluids.

.. code-block:: python

    # set both elements of the fluid vector
    >>> c1.set_attr(fluid={'water': 1})

    # set starting values (might be necessary sometimes)
    >>> c1.set_attr(fluid0={'water': 1})

    # overwrite the existing fluid vector
    >>> c1.fluid.is_set
    {'water'}
    >>> c1.set_attr(fluid={'N2': 0.7, "O2": 0.3})
    >>> c1.fluid.is_set
    {'N2', 'O2'}

    # remove a single specification while keeping the other component inside
    >>> c1.fluid.is_set.remove('N2')
    >>> c1.fluid.is_set
    {'O2'}
    >>> c1.fluid.val
    {'N2': 0.7, 'O2': 0.3}

CoolProp and REFPROP
++++++++++++++++++++

It is possible to specify the fluid property back end of the fluids by adding
the name of the back end in front of the fluid's name. For incompressible binary
mixtures, you can append the water volume/mass fraction to the fluid's name, for
example:

.. code-block:: python

    >>> c1.set_attr(fluid={'water': 1})  # HEOS back end
    >>> c1.set_attr(fluid={'INCOMP::water': 1})  # incompressible fluid
    >>> c1.set_attr(fluid={'BICUBIC::air': 1})  # bicubic back end
    >>> c1.set_attr(fluid={'INCOMP::MPG[0.5]|mass': 1})  # binary incompressible mixture

You can also specify REFPROP based fluids, e.g. R513A, which is a mass based
mixture of R134a and R1234yf:

.. code-block:: python

    >>> c1.set_attr(fluid={'REFPROP::R134A[0.44]&R1234yf[0.56]|mass': 1})  # REFPROP back end

.. note::

    Without further specifications CoolProp will be used as fluid property
    database. If you do not specify a back end, the **default back end**
    :code:`HEOS` will be used. For an overview of the back ends available please
    refer to the :ref:`fluid property section <fluid_properties_label>`.

Other Backends
++++++++++++++

You can also change the engine, for example to the iapws library. It is even
possible, that you define your own custom engine, e.g. using polynomial
equations. Please check out the fluid properties' section in the docs on how to
do this.

.. code-block:: python

    >>> from tespy.tools.fluid_properties.wrappers import IAPWSWrapper
    >>> c1.set_attr(fluid={'H2O': 1}, fluid_engines={"H2O": IAPWSWrapper})

Please also check out the section on
:ref:`custom fluid properties <fluid_properties_label>` for more
information.

.. _powerconnections_label:

PowerConnection Overview
------------------------

PowerConnections can be used to represent non-material energy flow, like power
or heat. You can make use of generators, motors and buses.

Different use-cases for the implementation of :code:`PowerConnection` with the
respective power components can be:

- apply motor or generator efficiencies
- connect multiple turbomachines on a single shaft
- collect all electricity production and own consumption to calculate net
  power

The handling of the :code:`PowerConnection` and the respective components is
identical to standard components.

The :code:`PowerConnection` only holds a single parameter, namely the power
flow :code:`E` (:math:`\dot E`), which is measured in Watts. You can create a
:code:`PowerConnection` instance by connecting to a component that has a
respective inlet or outlet. For example, consider a turbine generating
electricity. First we can set up a system as we are used to do without any
:code:`PowerConnections`:

.. code-block:: python

    >>> from tespy.components import Source, Sink, Turbine
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(
    ... temperature="degC", pressure="bar", power="kW"
    ... )
    >>> so = Source("source")
    >>> turbine = Turbine("turbine")
    >>> si = Sink("sink")
    >>> c1 = Connection(so, "out1", turbine, "in1", label="c1")
    >>> c2 = Connection(turbine, "out1", si, "in1", label="c2")
    >>> nw.add_conns(c1, c2)

We can parametrize the model, e.g. consider the turbine part of a gas turbine,
which expands hot flue gases:

.. code-block:: python

    >>> c1.set_attr(fluid={"air": 1}, p=10, T=1000, m=1)
    >>> c2.set_attr(p=1)
    >>> turbine.set_attr(eta_s=0.9)
    >>> nw.solve("design")
    >>> round(turbine.P.val)
    -577

We can add a connection between the turbine and the grid. This will add one
extra variable to our problem (the energy flow :code:`E`) but also one extra
equation, namely the turbine energy balance. Therefore, after adding the new
connection, there is nothing to change to make the model solve.

.. code-block:: python

    >>> from tespy.connections import PowerConnection
    >>> from tespy.components import PowerSink
    >>> grid = PowerSink("grid")
    >>> e1 = PowerConnection(turbine, "power", grid, "power", label="e1")
    >>> nw.add_conns(e1)
    >>> nw.solve("design")
    >>> round(e1.E.val)
    577

.. note::

    Note that the value of the energy flow of a :code:`PowerConnection` will
    always be positive in the defined direction (from one component's outlet
    to another component's inlet).

To learn what power connections are available in each of the component classes
see the respective API documentation. There are a couple of example
applications available for the :code:`PowerConnection`
:ref:`in this section <powerconnection_examples_label>`.

Access from the :code:`Network` object
--------------------------------------

You may want to access the network's connections or powerconnections other than
using the variable names, for example in an imported network or connections
from a subsystem. It is possible to access these using the connection's label
similar as it is possible for components.
By default, the label is generated by this logic:

:code:`source:source_id_target:target_id`, where

- :code:`source` and :code:`target` are the  labels of the components that are
  connected.
- :code:`source_id` and :code:`target_id` are e.g. :code:`out1` and
  :code:`in2` respectively.

.. code-block:: python

    >>> conn = nw.get_conn('c1')
    >>> conn.label
    'c1'
    >>> powerconn = nw.get_conn('e1')
    >>> powerconn.label
    'e1'
