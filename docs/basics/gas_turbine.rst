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
(:math:`LHV \cdot \dot{m}_{f}`).

Furthermore, we specify the efficiency :code:`eta` of the component, which
determines the heat loss as ratio of the thermal input. :code:`eta=1` means,
no heat losses, thus adiabatic behavior.

The pressure ratio :code:`pr` describes the ratio of the pressure at the
outlet to the pressure at **the inlet 1**. The pressure value at the inlet 2
is detached from the other pressure values, it must be a result of a different
parameter specification. In this example, we set it directly. Initially, we
assume adiabatic behavior :code:`eta=1` and no pressure losses :code:`pr=1`.

The ambient conditions as well as the fuel gas inlet temperature are defined
in the next step. The full vector for the air and the fuel gas composition
have to be defined. The component can not handle "Air" as input fluid. We can
run the code after the specifications.

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

It is also possible to make modifications on the fluid composition, for
example, we can add hydrogen to the fuel mixture.

.. literalinclude:: /../tutorial/basics/gas_turbine.py
    :language: python
    :start-after: [sec_6]
    :end-before: [sec_7]

The most convenient way to access the fluid composition is to access the
results dataframe for the connections.

.. literalinclude:: /../tutorial/basics/gas_turbine.py
    :language: python
    :start-after: [sec_7]
    :end-before: [sec_8]

.. note::

    All component and connection results are available in the :code:`results`
    dict of the :code:`Network` instance. The keys of the dictionary are the
    respective class names.

Setting up the Full System
^^^^^^^^^^^^^^^^^^^^^^^^^^
After learning more about the component, we are going to add the remaining
components: The turbine, the compressor and the generator. To do that, remove
the existing connections from the network, create the new connections and
add them to the network again. We also add a :code:`Bus` representing the
generator, assuming 98 % mechanical-electrical efficiency.

.. literalinclude:: /../tutorial/basics/gas_turbine.py
    :language: python
    :start-after: [sec_8]
    :end-before: [sec_9]

Since we deleted the connection 2 and 3, all specifications for those
connections have to be added again. The air fluid composition is specified on
connection 1 with ambient pressure and temperature. The compressor pressure
ratio is set to 15 bar, the turbine inlet temperature to 1200 °C. Finally, set
the gas turbine outlet pressure to ambient pressure as well as the
compressor's and turbine's efficiency.

.. literalinclude:: /../tutorial/basics/gas_turbine.py
    :language: python
    :start-after: [sec_9]
    :end-before: [sec_10]

Note, that the pressure of the fuel is lower than the pressure of the air at
the combustion chamber as we did not change the pressure of connection 5. A
respective warning is printed after the calculation. We can fix it like so:

.. literalinclude:: /../tutorial/basics/gas_turbine.py
    :language: python
    :start-after: [sec_10]
    :end-before: [sec_11]

We can investigate, how the turbine inlet temperature and the compressor
pressure ratio affect thermal efficiency and power generation. Also, we
assume 2 % heat losses and 3 % pressure losses in the combustion chamber.

.. dropdown:: Click to expand to code section

    .. literalinclude:: /../tutorial/basics/gas_turbine.py
        :language: python
        :start-after: [sec_11]
        :end-before: [sec_12]

.. figure:: /_static/images/basics/gas_turbine_parametric.svg
    :align: center
    :alt: Gas turbine performance at different compressor pressure ratios and turbine inlet temperatures

    Figure: Gas turbine performance at different compressor pressure ratios
    and turbine inlet temperatures.

Fluid Composition Specifications
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In this section you will learn how the fluid composition can be used as a
variable in such systems. To begin, we can impose the oxygen mass fraction on
the flue gas instead of the turbine inlet pressure, since it determines the
share of oxygen that is not required in the combustion. We can see, how the
turbine inlet temperature correlates with the oxygen mass fraction.

.. dropdown:: Click to expand to code section

    .. literalinclude:: /../tutorial/basics/gas_turbine.py
        :language: python
        :start-after: [sec_12]
        :end-before: [sec_13]

.. figure:: /_static/images/basics/gas_turbine_oxygen.svg
    :align: center
    :alt: Turbine inlet temperature at different levels of oxygen in the flue gas

    Figure: Turbine inlet temperature at different levels of oxygen in the
    flue gas.

Let us now assume, we do have an unknown shares of hydrogen and methane within
our fuel mixture. With the known mass flow of the fuel and an overall thermal
input, we can calculate both fractions by removing their respective values
from the input parameters and using the :code:`fluid_balance` keyword instead,
which automatically calculates the sum of all fluid mass fractions to be 1.

Investigate how changing the thermal input requires a different mixture of
hydrogen and methane.

.. attention::

    With this setup, a thermal input below the lower heating value of methane
    or above the lower heating value of hydrogen (each multiplied with the
    mass flow of 1 kg/s) does not make sense as input specification. This is
    individual of every fluid you use as fuel and you cannot easily abstract
    the values to any other combination.

.. dropdown:: Click to expand to code section

    .. literalinclude:: /../tutorial/basics/gas_turbine.py
        :language: python
        :start-after: [sec_13]
        :end-before: [sec_14]

.. figure:: /_static/images/basics/gas_turbine_fuel_composition.svg
    :align: center
    :alt: Mass fractions of H2 and CH4 in fuel mixture at different thermal input

    Figure: Mass fractions of H2 and CH4 in fuel mixture at varying thermal
    input and constant fuel mass flow.
