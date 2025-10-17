.. _tespy_tutorial_heat_exchanger:

Overview of Heat exchanger models
---------------------------------

This tutorial shows an overview of the different heat exchanger models
available in tespy. The table below indicates the different types and how they
compare in general. Further down the page we have created a problem, which
models the heat exchanger with a variety of the mentioned types and we show
the differences in the results.

Overview table
++++++++++++++

+-----------------------------+-----+-----------+-----------+----------------+
|            type             |     |   speed   | accuracy* | internal pinch |
+=============================+=====+===========+===========+================+
| SimpleHeatExchanger         | 0D  | fastest   | lowest    | ❌              |
+-----------------------------+-----+-----------+-----------+----------------+
| HeatExchanger               | 0D  | very fast | lower     | ❌              |
+-----------------------------+-----+-----------+-----------+----------------+
| ParallelFlowHeatExchanger   | 0D  | very fast | lower     | ❌              |
+-----------------------------+-----+-----------+-----------+----------------+
| Desuperheater               | 0D  | very fast | lower     | ❌              |
+-----------------------------+-----+-----------+-----------+----------------+
| Condenser                   | 0D  | very fast | low       | (✅)            |
+-----------------------------+-----+-----------+-----------+----------------+
| MovingBoundaryHeatExchanger | 1D  | mid       | high      | ✅              |
+-----------------------------+-----+-----------+-----------+----------------+
| SectionedHeatExchanger      | 1D  | slow      | highest   | ✅              |
+-----------------------------+-----+-----------+-----------+----------------+

.. note::

    - The accuracy depends on the context. In many context, the standard 0D
      heat exchanger components can be just as accurate as the 1D models.
    - The calculation speed also depends on context for the
      SectionedHeatExchanger:

      1. The number of specified sections is influencing the speed.
      2. If you impose :code:`td_pinch` or :code:`UA` to your model it will
         be significantly slower compared to other models, especially with a
         high number of sections.
      3. If you do not specify :code:`td_pinch` or :code:`UA` the sectioning is
         only applied once in the postprocessing, this takes more time than
         other types of heat exchangers but is still relatively fast.

Model comparisons
+++++++++++++++++

For the model comparison we have selected a typical problem: the condensation
of a working fluid in a heat pump to heat up water. The boundary conditions
are listed in the table below.

+---------------+----------------------+-------+------+
|     label     |    Specification     | value | unit |
+===============+======================+=======+======+
| c1            | fluid                | R290  | 1    |
+---------------+----------------------+-------+------+
|               | mass flow            | 5     | kg/s |
+---------------+----------------------+-------+------+
|               | dew line temperature | 60    | °C   |
+---------------+----------------------+-------+------+
|               | superheating         | 50    | °C   |
+---------------+----------------------+-------+------+
| c2            | subcooling           | 5     | °C   |
+---------------+----------------------+-------+------+
| d1            | fluid                | water | 1    |
+---------------+----------------------+-------+------+
|               | temperature          | 45    | °C   |
+---------------+----------------------+-------+------+
|               | pressure             | 1     | bar  |
+---------------+----------------------+-------+------+
| d2            | temperature          | 55    | °C   |
+---------------+----------------------+-------+------+
| heatexchanger | pressure drops       | 0     | bar  |
+---------------+----------------------+-------+------+

Given these boundary conditions, the following results can be obtained:

+--------------------------------+-------------------+-----------+-----------+
| type                           | minimum pinch (K) | kA (kW/K) | UA (kW/K) |
+================================+===================+===========+===========+
| HeatExchanger                  | n/a (10.0)        |  75.5     | n/a       |
+--------------------------------+-------------------+-----------+-----------+
| Condenser                      | n/a (5.0)         | 276.1     | n/a       |
+--------------------------------+-------------------+-----------+-----------+
| MovingBoundaryHeatExchanger    | 8.08              |  75.5     | 149.4     |
+--------------------------------+-------------------+-----------+-----------+
| SectionedHeatExchanger (50)    | 8.20              |  75.5     | 150.0     |
+--------------------------------+-------------------+-----------+-----------+
| SectionedHeatExchanger (6)     | 8.33              |  75.5     | 150.7     |
+--------------------------------+-------------------+-----------+-----------+

.. attention::

    Keep in mind: under different boundary conditions (e.g. no phase change)
    results of this comparison may vary a lot. There might be conditions, where
    the 0D components yield very similar results, or where deviation is even
    higher.

For the calculation of the results the following equations apply:

.. math::

    \Delta \theta_\text{log} = \frac{\Delta T_{i} - \Delta T_{i+1}}{\ln \frac{\Delta T_{i}}{\Delta T_{i+1}}}\\
    kA=\frac{\dot Q}{\Delta \theta_\text{log}}\\
    UA_{i}=\frac{\dot Q_{i}}{\Delta \theta_{\text{log,}i}}\\
    UA=\sum UA_{i}

For the calculation of :math:`kA` the terminal temperature differences
:code:`ttd_u` and :code:`ttd_l` are considered as the :math:`\Delta T_0` and
:math:`\Delta T_1`. For the calculation of UA, the internal temperature
differences :math:`\Delta T_{i}` are employed as highlighted in the figures
below.

We have created QT graphs to visualize the reason for the differences.

HeatExchanger
^^^^^^^^^^^^^

.. figure:: /_static/images/tutorials/heat_exchangers/HeatExchanger.svg
    :align: center
    :alt: HeatExchanger
    :figclass: only-light

    QT diagram for HeatExchanger class

.. figure:: /_static/images/tutorials/heat_exchangers/HeatExchanger_darkmode.svg
    :align: center
    :alt: HeatExchanger
    :figclass: only-dark

    QT diagram for HeatExchanger class

Condenser
^^^^^^^^^

In the condenser the upper terminal temperature differences is assigned to the
temperature differences between the dew line temperature of the condensing
fluid and the outlet temperature of the cold fluid.

.. figure:: /_static/images/tutorials/heat_exchangers/HeatExchanger.svg
    :align: center
    :alt: HeatExchanger
    :figclass: only-light

    QT diagram for HeatExchanger class

.. figure:: /_static/images/tutorials/heat_exchangers/HeatExchanger_darkmode.svg
    :align: center
    :alt: HeatExchanger
    :figclass: only-dark

    QT diagram for HeatExchanger class

MovingBoundaryHeatExchanger
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The moving boundary model sections the heat exchange into three different
sections at the phase change points.

.. figure:: /_static/images/tutorials/heat_exchangers/MovingBoundaryHeatExchanger.svg
    :align: center
    :alt: MovingBoundaryHeatExchanger
    :figclass: only-light

    QT diagram for MovingBoundaryHeatExchanger class

.. figure:: /_static/images/tutorials/heat_exchangers/MovingBoundaryHeatExchanger_darkmode.svg
    :align: center
    :alt: MovingBoundaryHeatExchanger
    :figclass: only-dark

    QT diagram for MovingBoundaryHeatExchanger class

SectionedHeatExchanger
^^^^^^^^^^^^^^^^^^^^^^

The sectioned model sections the heat exchange into 50 sections by default.

.. figure:: /_static/images/tutorials/heat_exchangers/SectionedHeatExchanger.svg
    :align: center
    :alt: SectionedHeatExchanger
    :figclass: only-light

    QT diagram for SectionedHeatExchanger class

.. figure:: /_static/images/tutorials/heat_exchangers/SectionedHeatExchanger_darkmode.svg
    :align: center
    :alt: SectionedHeatExchanger
    :figclass: only-dark

    QT diagram for SectionedHeatExchanger class

MovingBoundary and Sectioned models
+++++++++++++++++++++++++++++++++++

Comparing these two models, we see almost identical results in the cases
shown above. However, this is not necessarily the case. There are situations,
where these models yield different results.

1. when specifying a different number of sections
2. with pressure losses along the flow
3. when there is curvature in the isobars (e.g. supercritical conditions near
   critical point)

Comparing different number of sections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For variant 1 the following graph shows the comparison of a 0D model to a
sectioned one.

.. figure:: /_static/images/tutorials/heat_exchangers/SectionedHeatExchanger_vs_HeatExchanger.svg
    :align: center
    :alt: SectionedHeatExchanger vs. HeatExchanger
    :figclass: only-light

    QT diagram comparison for SectionedHeatExchanger and HeatExchanger classes

.. figure:: /_static/images/tutorials/heat_exchangers/SectionedHeatExchanger_vs_HeatExchanger_darkmode.svg
    :align: center
    :alt: SectionedHeatExchanger vs. HeatExchanger
    :figclass: only-dark

    QT diagram comparison for SectionedHeatExchanger and HeatExchanger classes

And two sectioned models with different number of sections.

.. figure:: /_static/images/tutorials/heat_exchangers/SectionedHeatExchanger_sectionscompare.svg
    :align: center
    :alt: SectionedHeatExchanger: 50 vs. 6 sections
    :figclass: only-light

    QT diagram comparison for SectionedHeatExchanger with different numbers of
    sections

.. figure:: /_static/images/tutorials/heat_exchangers/SectionedHeatExchanger_sectionscompare_darkmode.svg
    :align: center
    :alt: SectionedHeatExchanger: 50 vs. 6 sections
    :figclass: only-dark

    QT diagram comparison for SectionedHeatExchanger with different numbers of
    sections

Considering pressure drop
^^^^^^^^^^^^^^^^^^^^^^^^^
For the variant 2 we can get differences in the results between the 1D models
when there is a pressure drop in the two-phase region.

.. figure:: /_static/images/tutorials/heat_exchangers/Sectioned_vs_Moving_pressure_drop.svg
    :align: center
    :alt: SectionedHeatExchanger vs. MovingBoundaryHeatExchanger with pressure drop
    :figclass: only-light

    QT diagram comparison for SectionedHeatExchanger with
    MovingBoundaryHeatExchanger considering pressure drop

.. figure:: /_static/images/tutorials/heat_exchangers/Sectioned_vs_Moving_pressure_drop_darkmode.svg
    :align: center
    :alt: SectionedHeatExchanger vs. MovingBoundaryHeatExchanger with pressure drop
    :figclass: only-dark

    QT diagram comparison for SectionedHeatExchanger with
    MovingBoundaryHeatExchanger considering pressure drop

Curvature of isobars
^^^^^^^^^^^^^^^^^^^^
For variant 3 we yield different results due to the curvature of the isobars
in the supercritical region.

.. figure:: /_static/images/tutorials/heat_exchangers/Sectioned_vs_Moving_near_critical.svg
    :align: center
    :alt: SectionedHeatExchanger vs. MovingBoundaryHeatExchanger near critical point
    :figclass: only-light

    QT diagram comparison for SectionedHeatExchanger with
    MovingBoundaryHeatExchanger when supercritical near critical point

.. figure:: /_static/images/tutorials/heat_exchangers/Sectioned_vs_Moving_near_critical_darkmode.svg
    :align: center
    :alt: SectionedHeatExchanger vs. MovingBoundaryHeatExchanger near critical point
    :figclass: only-dark

    QT diagram comparison for SectionedHeatExchanger with
    MovingBoundaryHeatExchanger when supercritical near critical point
