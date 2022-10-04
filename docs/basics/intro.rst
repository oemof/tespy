.. _tespy_basics_intro_label:

Introduction
============

In this section we provide you with the basics of the TESPy's modeling concept.
TESPy builds up on **components** that are connected using **connections** to
form a topological **network**. The figure below highlights these three core
components of the software at the example of a small heat pump.

.. figure:: /_static/images/basics/modeling_concept.svg
    :align: center
    :alt: TESPy's modeling concept

    Figure: TESPy's modeling concept.

Set up a plant
--------------

In order to simulate a plant we start by creating the network
(:py:class:`tespy.networks.network.Network`). The network is the main container
for the model. You need to specify a list of the fluids you require for the
calculation in your plant. For more information on the fluid properties go to
the :ref:`corresponding section <tespy_fluid_properties_label>` in the
documentation.

.. literalinclude:: /../tutorial/basics/heat_pump.py
    :language: python
    :start-after: [sec_1]
    :end-before: [sec_2]

On top of that, it is possible to specify a unit system and value ranges for
the network's variables. If you do not specify these, TESPy will use SI-units.
We will thus only specify the unit systems, in this case.

.. literalinclude:: /../tutorial/basics/heat_pump.py
    :language: python
    :start-after: [sec_2]
    :end-before: [sec_3]

Now you can start to create the components of the network.

Set up components
-----------------

The list of components available can be found
:ref:`here <tespy_modules_components_label>`. If you set up a component you
have to specify a (within one network) unique label. Moreover, it is possible
to specify parameters for the component, for example power :math:`P` for a pump
or upper terminal temperature difference :math:`ttd_\mathrm{u}` of a heat
exchanger. The full list of parameters for a specific component is stated in
the respective class documentation. This example uses a compressor, a control
valve two (simple) heat exchangers and a so called cycle closer.

.. note::

    Parameters for components are generally optional. Only the component's
    label is mandatory. If an optional parameter is not specified by the user,
    it will be a result of the plant's simulation. This way, the set of
    equations a component returns is determined by which parameters you
    specify. You can find all equations in each component's documentation as
    well.

.. literalinclude:: /../tutorial/basics/heat_pump.py
    :language: python
    :start-after: [sec_3]
    :end-before: [sec_4]


After creating the components the next step is to connect them in order to form
your topological network.

Establish connections
---------------------

Connections are used to link two components (outlet of component 1 to inlet of
component 2: source to target). If two components are connected with each other
the fluid properties at the source will be equal to the properties at the
target. It is possible to set the properties on each connection in a similar
way as parameters are set for components. The basic specification options are:

* mass flow (m)
* volumetric flow (v)
* pressure (p)
* enthalpy (h)
* temperature (T)
* a fluid vector (fluid)

.. seealso::

    There are more specification options available. Please refer to
    the :ref:`connections section <tespy_modules_connections_label>` in the
    TESPy modules documentation for detailed information. The specification
    options are stated in the connection class documentation, too:
    :py:class:`tespy.connections.connection.Connection`.

After creating the connections, we need to add them to the network. As the
connections hold the information, which components are connected in which way,
we do not need to pass the components to the network.

.. literalinclude:: /../tutorial/basics/heat_pump.py
    :language: python
    :start-after: [sec_4]
    :end-before: [sec_5]

.. note::

    The :code:`CycleCloser` is a necessary component when working with closed
    cycles, because a system would always be overdetermined, if, for example,
    a mass flow is specified at some point within the cycle. It would propagate
    through all of the components, since they have an equality constraint for
    the mass flow at their inlet and their outlet. With the example here, that
    would mean: **Without the cycle closer** specification of massflow at an
    connection would lead to the following set of equations for massflow, which
    is an overdetermination:

    .. math::

        \begin{split}
            \dot{m}_1 = &\;\text{5 kg/s}\\
            \dot{m}_1 = &\;\dot{m}_2\\
            \dot{m}_2 = &\;\dot{m}_3\\
            \dot{m}_3 = &\;\dot{m}_4\\
            \dot{m}_4 = &\;\dot{m}_1\\
        \end{split}

    Similarly, this applies to the fluid composition.

    The cycle closer will prompt a warning, if the mass flow and fluid
    composition are not equal at its inlet and outlet.

We can set the component and connection parameters. In this example, we specify
the pressure losses (by outlet to inlet pressure ratio :code:`pr`) in the
condenser and the evaporator as well as the efficiency :code:`eta_s` of the
compressor. On top of that, the heat production of the heat pump can be set
with :code:`Q` for the condenser. Since we are working in **subcritical**
regime in this tutorial, we set the state of the fluid at the evaporator's
outlet to fully saturated steam (:code:`x=1`) and at the condenser's outlet to
fully saturated liqud (:code:`x=0`). On top of that, we want to impose the
condensation and the evaporation temperature levels. Last, we have to specify
the fluid vector at one point in our network.

.. literalinclude:: /../tutorial/basics/heat_pump.py
    :language: python
    :start-after: [sec_5]
    :end-before: [sec_6]

.. note::

    The sign convention for energy transfer by components is always from the
    perspective of the component. Energy entering the component means positive
    sign, energy leaving the component's system boundary means negative sign.

Start your calculation
----------------------

After building your network, the components and the connections, add the
following line at the end of your script and run it. You can calculate the COP
with the respective component parameters.

.. literalinclude:: /../tutorial/basics/heat_pump.py
    :language: python
    :start-after: [sec_6]
    :end-before: [sec_7]

Next steps
----------

We highly recommend to check our other basic model examples on how to set up
different standard thermodynamic cycles in TESPy. The heat pump cycle in that
section builds on this heat pump. We will introduce couple of different inputs
and show, how to change the working fluid. The other tutorials show the usage
of more components, for example the combustion chamber and the turbine or a
condenser including the cooling water side of the system.

In the more advanced tutorials, you will learn, how to set up more complex
plants ste by step, make a design calculation of the plant as well as calculate
offdesign/partload performance.

In order to get a good overview of the TESPy functionalities, the sections on
the :ref:`TESPy modules <tespy_modules_label>` will guide you in detail.
