.. _tespy_examples_label:

~~~~~~~~~~~~~~~~~~~~
Example Applications
~~~~~~~~~~~~~~~~~~~~

On this page we collect example applications of TESPy. If you want to add your
example here, please open an issue on GitHub and let us know. The source code
of the application should be accessible freely, so other users can learn from
your project.


.. card::

    **Coupled Porous Media Storage and Power Plant Simulation**
    ^^^

    .. image:: /_static/images/examples/PM_CAES_graphical-abstract.svg
      :align: center
      :alt: Overall setup of a porous media compressed air energy storage in the energy system
      :class: only-light
      :target: https://github.com/fgasa/Integrated_CASE_Assessment

    .. image:: /_static/images/examples/PM_CAES_graphical-abstract_darkmode.svg
      :align: center
      :alt: Overall setup of a porous media compressed air energy storage in the energy system
      :class: only-dark
      :target: https://github.com/fgasa/Integrated_CASE_Assessment

    A porous media energy storage is coupled with a power plant simulation.
    Yearly operation of different power plant and storage setups is optimized
    within different energy system scenarios in a first step by transferring
    information from the complex system to a mixed-integer linear formulation.
    The optimized storage dispatch is then passed to the coupled simulation to
    find actual mass flow rates as well as pressure subsurface distribution.

    +++
    Title: Integration of geological compressed air energy storage into future energy supply systems dominated by renewable power sources

    Authors: Firdovsi Gasanzade, Francesco Witte, Ilja Tuschy, Sebastian Bauer

    Reference: :cite:`Gasanzade2023`


.. card::

    **Parametric Optimization of an Organic Rankine Cycle**
    ^^^

    .. image:: /_static/images/examples/ORC_parametric_flowsheet.svg
      :align: center
      :alt: Flowsheet and Ts-diagram of the ORC plant
      :class: only-light
      :target: https://github.com/fwitte/ORCSimulator

    .. image:: /_static/images/examples/ORC_parametric_flowsheet_darkmode.svg
      :align: center
      :alt: Flowsheet and Ts-diagram of the ORC plant
      :class: only-dark
      :target: https://github.com/fwitte/ORCSimulator

    An ORC power plant using two-phase geothermal sources is designed, and an
    optimization is carried out. The plant's performance is investigated for
    six different working fluids. Gross and net power output are optimized.
    The open source library pygmo :cite:`Biscani2020` is applied in combination
    with TESPy for the first time, setting the foundation for the optimization
    API of TESPy (:py:class:`tespy.tools.optimization.OptimizationProblem`).

    +++
    Title: Parametric optimization and comparative study of an organic Rankine
    cycle power plant for two-phase geothermal sources

    Authors: Chaofan Chen, Francesco Witte, Ilja Tuschy, Olaf Kolditz, Haibing
    Shao

    Reference: :cite:`Chen2022`


.. card::

    **Combined Heat and Power Organic Rankine Cycle**
    ^^^

    .. image:: /_static/images/examples/GRC_flowsheet.svg
      :align: center
      :alt: Development of the Electrical Power Output of the ORC for a District with 2 MW Peak Heat Load
      :class: only-light
      :target: https://github.com/fwitte/chp_orc

    .. image:: /_static/images/examples/GRC_flowsheet_darkmode.svg
      :align: center
      :alt: Development of the Electrical Power Output of the ORC for a District with 2 MW Peak Heat Load
      :class: only-dark
      :target: https://github.com/fwitte/chp_orc

    Starting from well production information for a geothermal energy reservoir
    over a lifetime of 40 years, the development of the electrical power output
    of an ORC is monitored within different designs of the plant. The
    geothermal heat source is exploited to provide heat to a district heating
    system and the residual heat is used to operate the orc cycle.

    +++
    Title: Computational Modeling of Organic Rankine Cycle Combined Heat and
    Power for Sedimentary Geothermal Exploitation

    Authors: Nicholas Fry, Jessica Eagle-Bluestone, Francesco Witte

    Reference: :cite:`Fry2022`
