.. _tespy_basics_gas_turbine_label:

Gas Turbine
===========

.. figure:: /_static/images/basics/gas_turbine.svg
    :align: center
    :alt: Topology of the gas turbine
    :figclass: only-light

    Figure: Topology of the gas turbine

.. figure:: /_static/images/basics/gas_turbine_darkmode.svg
    :align: center
    :alt: Topology of the gas turbine
    :figclass: only-dark

    Figure: Topology of the gas turbine

This tutorial introduces a new component, the combustion chamber. You will
learn how to use the component and set up a simple open cycle gas turbine: It
compresses air and burns fuel in the combustion chamber. The hot and
pressurized flue gas expands in the turbine, which drives the compressor and
the generator. You will also learn, how to use the fluid compositon as a
variable in your simulation.

Setting up the Combustion Chamber
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We are setting up our system step by step. Especially for larger systems, it
is recommended you follow this approach, since TESPy highly relies on a set of
good starting values for good convergence. You can learn more about it in the
:ref:`advanced tutorial <tespy_tutorial_label>` section of the online
documentation.

.. note::

    There are two different types of combustion chambers available:

    - :py:class:`tespy.components.combustion.base.CombustionChamber` and
    - :py:class:`tespy.components.combustion.diabatic.DiabaticCombustionChamber`.

    Both can handle varying fluid compositions for the air and the fuel and
    calculate the fluid composition of the flue gas. Thus, it is possible to
    e.g. specify the oxygen mass fraction in the flue gas in a calculation.
    The difference between the components lies in the fact, that the
    :code:`CombustionChamber` does **not consider heat or pressure losses**,
    while :code:`DiabaticCombustionChamber` does so.

In this tutorial, we will use the
:py:class:`tespy.components.combustion.diabatic.DiabaticCombustionChamber`.
First, we set up a network and the components. The network's fluid list must
contain all fluid components used for the combustion chamber. **These are at**
**least the fuel, oxygen, carbon-dioxide and water**. For this example we
add Nitrogen, since it is the most important fresh air component.

.. literalinclude:: /../tutorial/basics/gas_turbine.py
    :language: python
    :start-after: [sec_1]
    :end-before: [sec_2]

In the first step, we do not connect the inlet of the combustion chamber with
the compressor but with the air source instead. Similarly, the outlet of the
combustion chamber is directly connected to the flue gas sink.

.. literalinclude:: /../tutorial/basics/gas_turbine.py
    :language: python
    :start-after: [sec_2]
    :end-before: [sec_3]

There are many different specifications possible. For the combustion chamber
we will specify its air to stoichiometric air ratio lamb and the thermal input
(:math:`LHV \cdot \dot{m}_{f}`). The ambient conditions as well as the fuel
gas inlet temperature are defined in the next step. The air and the fuel gas
composition are fully be stated, the component combustion chamber can not
handle "Air" as input fluid. Then, we can run the code.

.. literalinclude:: /../tutorial/basics/gas_turbine.py
    :language: python
    :start-after: [sec_3]
    :end-before: [sec_4]

Of course, you can change the parametrization in any desired way. For example
instead of stating the thermal input, you could choose any of the mass flows:

.. literalinclude:: /../tutorial/basics/gas_turbine.py
    :language: python
    :start-after: [sec_4]
    :end-before: [sec_5]

or instead of the air to stoichiometric air ratio you could specify the flue
gas temperature.

.. literalinclude:: /../tutorial/basics/gas_turbine.py
    :language: python
    :start-after: [sec_5]
    :end-before: [sec_6]

It is also possible to make modifications on the fluid
composition, for example stating the oxygen content in the flue gas or to
change the fuel composition. Make sure, all desired fuels of your fuel mixture
are also within the fluid_list of the network. For the example below we added
hydrogen to the fuel mixture.

.. literalinclude:: /../tutorial/basics/gas_turbine.py
    :language: python
    :start-after: [sec_6]
    :end-before: [sec_7]


.. note::

    A warning message is prompted at the end of the simulation, if the pressure
    of the inlet 2 is lower or equal to the pressure of inlet 1.

Furthermore, we specify the efficiency
:code:`eta` of the component, which determines the heat loss as ratio of the
thermal input. :code:`eta=1` means, no heat losses, thus adiabatic behavior.
On top of that, we set the pressure ratio :code:`pr`, which describes the
ratio of the pressure at the outlet to the pressure at **the inlet 1**. The
pressure value at the inlet 2 is detached from the other pressure values, it
must be a result of a different parameter specification. In this example, we
set it directly. To match the inputs of the first tutorial, we set
:code:`pr=1` and :code:`p=1` for connection :code:`sf_comb`.

Now, consider heat loss of the surface of the component. This is simply done by
specifying the value for :code:`eta`. We assume 4 % of thermal input as heat
loss and set that value accordingly. Furthermore, the pressure of the fuel is
set to 1.5 bar. The air inlet pressure will be the result of the specified
pressure ratio and the outlet pressure assuming 2 % pressure losses. All
other parameters stay untouched.

.. literalinclude:: /../tutorial/basics/gas_turbine.py
    :language: python
    :start-after: [sec_7]
    :end-before: [sec_8]

Setting up the Full System
^^^^^^^^^^^^^^^^^^^^^^^^^^
- Add remaining parts (del_conns/add_conns)

.. literalinclude:: /../tutorial/basics/gas_turbine.py
    :language: python
    :start-after: [sec_8]
    :end-before: [sec_9]

- Run full model

.. literalinclude:: /../tutorial/basics/gas_turbine.py
    :language: python
    :start-after: [sec_9]
    :end-before: [sec_10]

Fluid Composition Specifications
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- O2 in flue gas fraction

.. literalinclude:: /../tutorial/basics/gas_turbine.py
    :language: python
    :start-after: [sec_10]
    :end-before: [sec_11]

- CO2, CH4, H2 mixture, make it flexible

.. literalinclude:: /../tutorial/basics/gas_turbine.py
    :language: python
    :start-after: [sec_11]
    :end-before: [sec_12]
