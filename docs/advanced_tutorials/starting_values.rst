.. _tutorial_starting_values_label:

How to Generate Stable Starting Values
--------------------------------------
Applying numerical algorithms and methods, the starting value of a variable
is the value used for the first iteration. With more complex TESPy models
it can happen that the simulation does not converge easily due to a
combination of "bad" starting values. The solver is especially vulnerable if
the specified parameters trigger complex equations with respect to the primary
variables.

The primary variables of TESPy are mass flow, pressure, enthalpy and fluid
composition. If such a value is directly specified by the user, the solver has
a solution for this value before starting the first iteration. Therefore,
specifying a set of parameters largely including primary variables will
improve the convergence significantly. Based on the converged solution of a
initial simulation, it is then possible to adjust the parameters, for example,
unsetting pressure values and specifying efficiencies instead.

Here we provide a short tutorial for you to better understand, how this
process could look like at the example of a subcritical heat pump with
different working fluids.

.. note::

    If the heat pump operates in trans- or supercritical range, some
    modifications have to be made on this setup. We plan to include respective
    examples here in the future.

You can download the full code of this example here:
:download:`starting_values.py </../tutorial/advanced/starting_values.py>`

Topology of the heat pump
^^^^^^^^^^^^^^^^^^^^^^^^^
Following the first tutorial a slightly different topology for a heat pump
with internal heat exchangers is considered instead of dumping the heat to the
ambient. You can see the plant topology in the figure below.

.. figure:: /_static/images/tutorials/heat_pump_starting_values/flowsheet.svg
    :align: center
    :alt: Topology of heat pump with internal heat exchanger
    :figclass: only-light

    Figure: Topology of heat pump with internal heat exchanger

.. figure:: /_static/images/tutorials/heat_pump_starting_values/flowsheet_darkmode.svg
    :align: center
    :alt: Topology of heat pump with internal heat exchanger
    :figclass: only-dark

    Figure: Topology of heat pump with internal heat exchanger

The system consists of a consumer system, a valve, an evaporator system, a
compressor and additionally an internal heat exchanger. In order to simulate
this heat pump, the TESPy model has to be built up. First, the network has to
be initialized, and the refrigerants used have to be specified. This example
shows how to make the heat pump model work with a variety of working fluids
with water on both the heat source and heat sink side of the system.

Running into errors
^^^^^^^^^^^^^^^^^^^
As always, we start by importing the necessary TESPy classes.

.. literalinclude:: /../tutorial/advanced/starting_values.py
    :language: python
    :start-after: [sec_1]
    :end-before: [sec_2]

Then, we can build the network by defining components and connections. The
working fluid will be set with the variable `wf`, `"NH3"` is used in the first
setup. This way, we will be able to change the working fluid in a flexible
way.

.. dropdown:: Click to expand to code section

    .. literalinclude:: /../tutorial/advanced/starting_values.py
        :language: python
        :start-after: [sec_2]
        :end-before: [sec_3]

After setting up the topology, the system's parameters should be set in the
following way:

- Heat sink temperature levels (`T` at 23 and 24)
- Heat source temperature levels (`T` at 11 and 13)
- Degree of overheating after the internal heat exchanger (`td_dew` at 2)
- Pinch point temperature difference at the evaporator (`ttd_l`) to derive
  evaporation pressure
- Temperature difference at the condenser (`ttd_u`) to derive condensation
  pressure
- Saturated gaseous state of the working fluid (`x=1`) after leaving the
  evaporator
- Efficiencies of pumps and the compressor (`eta_s`)
- Pressure losses in all heat exchangers (`pr1`, `pr2`, `pr`)
- Consumer heat demand (`Q`)

.. dropdown:: Click to expand to code section

    .. literalinclude:: /../tutorial/advanced/starting_values.py
        :language: python
        :start-after: [sec_3]
        :end-before: [sec_4]

The system should be well defined with the parameter settings, however no
solution can be found. We might run in some error, like

.. error::

    .. code-block:: bash

        ERROR:root:Singularity in jacobian matrix, calculation aborted! Make
        sure your network does not have any linear dependencies in the
        parametrisation. Other reasons might be

        -> given temperature with given pressure in two phase region, try
        setting enthalpy instead or provide accurate starting value for
        pressure.

        -> given logarithmic temperature differences or kA-values for heat
        exchangers,

        -> support better starting values.

        -> bad starting value for fuel mass flow of combustion chamber, provide
        small (near to zero, but not zero) starting value.

or simply not making progress in the convergence

.. error::

    .. code-block:: bash

        WARNING:root:The solver does not seem to make any progress, aborting
        calculation. Residual value is 7.43e+05. This frequently happens, if
        the solver pushes the fluid properties out of their feasible range.

Fixing the errors
^^^^^^^^^^^^^^^^^

To generate good starting values for the simulation, it is recommended to set
pressure and enthalpy values instead of temperature differences. In this
example, fixed points can be identified with the help of the logph diagram
which you can see in the figure below.

.. figure:: /_static/images/tutorials/heat_pump_starting_values/logph.svg
    :align: center

    Figure: Logph diagram of ammonia

A rough estimation of the evaporation and condensation pressure can be
obtained and will be used to replace the temperature differences at the
evaporator and the condenser for the starting value generator. After
condensation, the working fluid is in saturated liquid state. We can retrieve
the condensation pressure corresponding to a temperature slightly below the
heat sink temperature by using the CoolProp `PropsSI` interface with the
respective inputs. The same step can be carried out on the heat source side.
For the internal heat exchanger, an enthalpy value is specified instead of the
temperature difference to the boiling point as well. It is important to note
that the PropertySI function (PropsSI) is used with SI units, which differ
from the units defined in the network.

The temperature difference values are unset and pressure and enthalpy values
are set instead.

.. literalinclude:: /../tutorial/advanced/starting_values.py
    :language: python
    :start-after: [sec_4]
    :end-before: [sec_5]

The model was solved successfully and has stored the starting values for any
follow-up. Therefore, we can undo our recent changes and restart the
simulation. For example, the COP is then calculated.

.. literalinclude:: /../tutorial/advanced/starting_values.py
    :language: python
    :start-after: [sec_5]
    :end-before: [sec_6]

Expand fix to any working fluids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Finally, using this strategy, it is possible to build a generic function,
building a network, that works with a variety of working fluids.

.. dropdown:: Click to expand to code section

    .. literalinclude:: /../tutorial/advanced/starting_values.py
        :language: python
        :start-after: [sec_6]
        :end-before: [sec_7]

We can run that function for different working fluids and plot the results:

.. literalinclude:: /../tutorial/advanced/starting_values.py
    :language: python
    :start-after: [sec_7]
    :end-before: [sec_8]

.. figure:: /_static/images/tutorials/heat_pump_starting_values/COP_by_wf.svg
    :align: center
    :alt: Analysis of the COP using different working fluids
    :figclass: only-light

    Figure: Analysis of the COP using different working fluids

.. figure:: /_static/images/tutorials/heat_pump_starting_values/COP_by_wf_darkmode.svg
    :align: center
    :alt: Analysis of the COP using different working fluids
    :figclass: only-dark

    Figure: Analysis of the COP using different working fluids
