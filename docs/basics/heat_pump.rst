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

Flexibility in Modeling
^^^^^^^^^^^^^^^^^^^^^^^
In TESPy the specifications for components and/or connections are
interchangable in every possible way, provided that the system of equations
representing the plant is well defined.

For example, instead of the heat provided by the condenser we could specify
the mass flow :code:`m` of the refrigerant. To unset a parameter you need to
set it to :code:`None`. To replace the specification, set the mass flow of
connection c1 to 5 kg/s:

.. code-block:: python

    co.set_attr(Q=None)
    c1.set_attr(m=5)

    my_plant.solve('design')
    my_plant.print_results()

You can observe, that the heat transferred by the condenser now is a result of
the mass flow imposed. We could do similar things, for example with the heat
sink temperature. We imposed it in our initial set up. Now we want to insert
a compressor with a fixed output to input pressure ratio. In that case, we
cannot choose the condensation temperature but it will be a result of that
specification:

.. code-block:: python

    cp.set_attr(pr=4)
    c4.set_attr(T=None)

    my_plant.solve('design')
    my_plant.print_results()

Or, we have a plant running with data observation running. It tells us the
compressor outlet temperature and we want to know what the efficiency of the
compressor would be, in case we measure :code:`T=97.3` at connection 3.

.. code-block:: python

    cp.set_attr(pr=None, eta_s=None)
    c3.set_attr(T=97.3)
    c4.set_attr(T=80)

    my_plant.solve('design')
    my_plant.print_results()

Typical Errors
^^^^^^^^^^^^^^
If you over- or underdetermine the system by specifying too few or too many
parameters, you will get an error message. We could set the heat demand and the
mass flow at the same time.

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

    .. code-block:: python

        import matplotlib.pyplot as plt
        import numpy as np


        data = {
            'T_source': np.linspace(0, 40, 11),
            'T_sink': np.linspace(60, 100, 11),
            'eta_s': np.linspace(0.75, 0.95, 11) * 100
        }
        COP = {
            'T_source': [],
            'T_sink': [],
            'eta_s': []
        }
        description = {
            'T_source': 'Evaporation temperature in °C',
            'T_sink': 'Condensation temperature in °C',
            'eta_s': 'Isentropic efficiency in %'
        }

        for T in data['T_source']:
            c2.set_attr(T=T)
            my_plant.solve('design')
            COP['T_source'] += [abs(co.Q.val) / cp.P.val]

        # reset to base temperature
        c2.set_attr(T=20)

        for T in data['T_sink']:
            c4.set_attr(T=T)
            my_plant.solve('design')
            COP['T_sink'] += [abs(co.Q.val) / cp.P.val]

        # reset to base temperature
        c4.set_attr(T=80)

        for eta_s in data['eta_s']:
            cp.set_attr(eta_s=eta_s / 100)
            my_plant.solve('design')
            COP['eta_s'] += [abs(co.Q.val) / cp.P.val]

        fig, ax = plt.subplots(1, 3, sharey=True, figsize=(16, 8))

        [a.grid() for a in ax]

        i = 0
        for key in data:
            ax[i].scatter(data[key], COP[key])
            ax[i].set_xlabel(description[key])
            i += 1

        ax[0].set_ylabel('COP of the heat pump')
        plt.tight_layout()
        fig.savefig('heat_pump_parametric.svg')

.. figure:: /_static/images/basics/heat_pump_parametric.svg
    :align: center
    :alt: Parametric analysis of the heat pump's COP

    Figure: Parametric analysis of the heat pump's COP

The figure shows the results of the COP analysis. The base case is at an
evaporation temperature of 20 °C, the condensation temperature at 80 °C and the
isentropic effficiency of the compressor at 85 %.
