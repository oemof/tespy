.. _fluprodia_label:

Creating Fluid Property Diagrams
================================

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
