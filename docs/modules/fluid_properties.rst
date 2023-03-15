.. _tespy_fluid_properties_label:

Fluid properties
================
The basic fluid properties are handled by
`CoolProp <http://www.coolprop.org/>`_. All available fluids can be found on
their homepage. Also see :cite:`Bell2014`.

CoolProp back ends
------------------
CoolProp provides multiple back ends for fluid property calculation. The
back ends vary in calculation speed and calculation accuracy. It is possible
to choose from the following back ends:

- :code:`HEOS`: Helmhotz Equation Of State with highest accuracy and lowest
  calculation speed. **This is the default back end!**
- :code:`REFPROP`: Highest accuracy and highest convergence stability.
  **This back end is not free**, a separate
  `REFPROP <https://www.nist.gov/srd/refprop>`__ license is required.
- :code:`BICUBIC`: Tabular back end with high accuracy and very high
  calculation speed.
- :code:`TTSE`: Tabular back end with lowest accuracy and very high calculation
  speed.
- :code:`INCOMP`: Back end for incompressible fluids.
- :code:`IF97`: Back end for the IAPWS-IF97 of water, very accurate and much
  higher calculation speed than :code:`HEOS`. Due to a bug in the CoolProp
  back end this option is available with a fix (not the original
  implementation), for more information see the
  `CoolProp issue #1918 <https://github.com/CoolProp/CoolProp/issues/1918/>`_.

For more information on the Back ends please visit the CoolProp online
documentation.

Pure and pseudo-pure fluids
---------------------------
If you use pure fluids, TESPy directly uses CoolProp functions to gather all
fluid properties. CoolProp covers the most important fluids such as water, air
as a pseudo-pure fluid as well as its components, several fuels and
refrigerants etc.. Look for the aliases in the list of
`fluids <http://www.coolprop.org/fluid_properties/PurePseudoPure.html>`__.
All fluids provided in this list cover liquid and gaseous state and the
two-phase region.

Incompressible fluids
---------------------
If you are looking for heat transfer fluids, the list of incompressible
`fluids <http://www.coolprop.org/fluid_properties/Incompressibles.html>`__
might be interesting for you. In contrast to the pure fluids, the properties
cover liquid state only.

Fluid mixtures
--------------
CoolProp provides fluid properties for two component mixtures. BUT: These are
NOT integrated in TESPy! Nevertheless, you can use fluid mixtures for gases.

Ideal mixtures of gaseous fluids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
TESPy can handle mixtures of gaseous fluids, by using the single fluid
properties from CoolProp together with corresponding equations for mixtures.
The equations can be found in the
:py:mod:`fluid_properties <tespy.tools.fluid_properties>` module and are
applied automatically to the fluid vector.

Other mixtures
^^^^^^^^^^^^^^
Apart from partially liquid water in flue gases it is **not possible** to use
mixtures of liquids and other liquids or gaseous fluids **at the moment**! If
you try to use a mixture of two liquid or gaseous fluids and liquid fluids,
e.g. water and methanol, the equations will still be applied, but obviously
return wrong values. If you have ideas for the implementation of new kinds of
mixtures we appreciate you contacting us.

.. _FluProDia_label:

Creating Fluid Property Diagrams
--------------------------------

.. figure:: /_static/images/modules/logph_diagram_states.svg
    :align: center
    :alt: logph diagram of NH3 with a simple heat pump cycle

    Figure: logph diagram of NH3 with a simple heat pump cycle

.. figure:: /_static/images/modules/Ts_diagram_states.svg
    :align: center
    :alt: Ts diagram of NH3 with a simple heat pump cycle

    Figure: Ts diagram of NH3 with a simple heat pump cycle

CoolProp has an inbuilt feature for creating fluid property diagrams.
Unfortunately, the handling is not very easy at the moment. We recommend using
fluprodia (Fluid Property Diagram) instead. You can create and customize
different types of diagrams for all pure and pseudo-pure fluids available in
CoolProp. In order to plot your process data into a diagram, you can use the
:code:`get_plotting_data` method of each component. The method returns a
dictionary, that can be passed as :code:`**kwargs` to the
:code:`calc_individual_isoline` method of a fluprodia
:code:`FluidPropertyDiagram` object. The fluprodia documentation provides
examples of how to plot a process into different diagrams, too. For more
information on fluprodia have a look at the
`online documentation <https://fluprodia.readthedocs.io/en/latest/>`_. You can
install the package with pip.

.. code-block:: bash

    pip install fluprodia

.. note::

    The plotting data a returned from the :code:`get_plotting_data` as a
    nested dictionary. The first level key contains the connection id of the
    state change (change state from incoming connection to outgoing
    connection). The table below shows the state change and the respective id.

    .. list-table:: State change and respective ids of dictionary
       :widths: 60 10 10 10
       :header-rows: 1

       * - component
         - state from
         - state to
         - id
       * - components with one inlet and one outlet only
         - :code:`in1`
         - :code:`out1`
         - :code:`1`
       * - class HeatExchanger and subclasses
         - :code:`in1`
         - :code:`out1`
         - :code:`1`
       * -
         - :code:`in2`
         - :code:`out2`
         - :code:`2`
       * - class ORCEvaporator
         - :code:`in1`
         - :code:`out1`
         - :code:`1`
       * -
         - :code:`in2`
         - :code:`out2`
         - :code:`2`
       * -
         - :code:`in3`
         - :code:`out3`
         - :code:`3`
       * - class Merge
         - :code:`in1`
         - :code:`out1`
         - :code:`1`
       * -
         - :code:`in2`
         - :code:`out1`
         - :code:`2`
       * -
         - ...
         - ...
         - ...
       * - class Drum
         - :code:`out1`
         - :code:`out2`
         - :code:`1`

    All other components do not return any information as either there is no
    change in state or the state change is accompanied by a change in fluid
    composition.
