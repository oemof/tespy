.. _tespy_basics_heat_pump_label:

Heat Pump
=========

.. figure:: /_static/images/basics/heat_pump.svg
    :align: center
    :alt: Topology of the heat pump
    :figclass: only-light

    Figure: Topology of the heat pump

.. figure:: /_static/images/basics/heat_pump_darkmode.svg
    :align: center
    :alt: Topology of the heat pump
    :figclass: only-dark

    Figure: Topology of the heat pump

This tutorial is the continuation of the tutorial in the
:ref:`introduction <tespy_basics_intro_label>`. First, we have a look at
specification options and modeling flexibility in TESPy as well as some
typical errors, that might occur when using the software. Then we will make a
simple analysis of the COP of the heat pump based on several input parameters.

Download the full script here:
:download:`heat_pump.py </../tutorial/basics/heat_pump.py>`

Flexibility in Modeling
^^^^^^^^^^^^^^^^^^^^^^^
In TESPy the specifications for components and/or connections are
interchangable in every possible way, provided that the system of equations
representing the plant is well defined.

For example, instead of the heat provided by the condenser we could specify
the mass flow :code:`m` of the refrigerant. To unset a parameter you need to
set it to :code:`None`. To replace the specification, set the mass flow of
connection c1 to 5 kg/s:

.. literalinclude:: /../tutorial/basics/heat_pump.py
    :language: python
    :start-after: [sec_7]
    :end-before: [sec_8]

You can observe, that the heat transferred by the condenser now is a result of
the mass flow imposed. We could do similar things, for example with the heat
sink temperature. We imposed it in our initial set up. Now we want to insert
a compressor with a fixed output to input pressure ratio. In that case, we
cannot choose the condensation temperature but it will be a result of that
specification:

.. literalinclude:: /../tutorial/basics/heat_pump.py
    :language: python
    :start-after: [sec_8]
    :end-before: [sec_9]

Or, we have a plant running with data observation running. It tells us the
compressor outlet temperature and we want to know what the efficiency of the
compressor would be, in case we measure :code:`T=97.3` at connection 3.

.. literalinclude:: /../tutorial/basics/heat_pump.py
    :language: python
    :start-after: [sec_9]
    :end-before: [sec_10]

Typical Errors
^^^^^^^^^^^^^^
If you over- or underdetermine the system by specifying too few or too many
parameters, you will get an error message. We could set the heat demand and the
mass flow at the same time.

.. attention::

    The two code examples in this section are not included in the downloadable
    script!

.. code-block:: python

    co.set_attr(Q=-1e6)
    c1.set_attr(m=5)

    my_plant.solve('design')

.. error::

    .. code-block:: bash

        ERROR:root:You have provided too many parameters: 20 required, 21 supplied. Aborting calculation!
        Traceback (most recent call last):
        File ".\tutorial\basics\heat_pump.py", line 59, in <module>
            my_plant.solve('design')
        File "c:\users\user\documents\github\tespy\src\tespy\networks\network.py", line 1623, in solve
            self.solve_determination()
        File "c:\users\user\documents\github\tespy\src\tespy\networks\network.py", line 1749, in solve_determination
            raise hlp.TESPyNetworkError(msg)
        tespy.tools.helpers.TESPyNetworkError: You have provided too many parameters: 20 required, 21 supplied. Aborting calculation!

If you make a specification that leads to the correct amount of parameters but
causes a linear dependency in the system of equations, the error message cannot
be that clear to you. To make an easy example, we can set mass flow on the
connections 1 and 2 with the heat demand and the evaporation temperature unset.
In this case the number of equations will be correct, but the specification
obviously does not make any sense.

.. code-block:: python

    co.set_attr(Q=None)
    c1.set_attr(m=5)
    c2.set_attr(m=5, T=None)

    my_plant.solve('design')

.. error::

    .. code-block:: bash

        ERROR:root:Singularity in jacobian matrix, calculation aborted! Make sure your network does not have any linear dependencies in the parametrisation. Other reasons might be
        -> given temperature with given pressure in two phase region, try setting enthalpy instead or provide accurate starting value for pressure.
        -> given logarithmic temperature differences or kA-values for heat exchangers,
        -> support better starting values.
        -> bad starting value for fuel mass flow of combustion chamber, provide small (near to zero, but not zero) starting value.

.. seealso::

    For more detailed information about the number of variables involved and
    ways of parameter specifications, please go to the
    :ref:`TESPy modules section <tespy_modules_label>` inside the Documentation
    chapter.

    Another frequent reason for such errors are bad starting values. We have a
    tutorial specifically dedicated to this topic
    :ref:`here <tespy_tutorial_starting_values_label>`.

Parametric Analysis of COP
^^^^^^^^^^^^^^^^^^^^^^^^^^
For a constant amount of heat production, we will investigate the influence of

* the source temperature level
* the sink temperature level and
* the isentropic efficiency of the compressor.

To do this, we import the numpy and matplotlib package and define the ranges
and then iterate over a loop and restart the simulation for every input
parameter. After each loop, we set the respective parameter back to its
original value. We collect the results in lists and can finally make a scatter
plot using matplotlib.

.. dropdown:: Click to expand to code section

    .. literalinclude:: /../tutorial/basics/heat_pump.py
        :language: python
        :start-after: [sec_10]
        :end-before: [sec_11]

.. figure:: /_static/images/basics/heat_pump_parametric.svg
    :align: center
    :alt: Parametric analysis of the heat pump's COP
    :figclass: only-light

    Figure: Parametric analysis of the heat pump's COP

.. figure:: /_static/images/basics/heat_pump_parametric_darkmode.svg
    :align: center
    :alt: Parametric analysis of the heat pump's COP
    :figclass: only-dark

    Figure: Parametric analysis of the heat pump's COP

The figure shows the results of the COP analysis. The base case is at an
evaporation temperature of 20 °C, the condensation temperature at 80 °C and the
isentropic effficiency of the compressor at 85 %.
