.. _tespy_examples_label:

~~~~~~~~~~~~~~~~~~~~
Example Applications
~~~~~~~~~~~~~~~~~~~~

On this page we collect example applications of TESPy. If you want to add your
example here, please open an issue on GitHub and let us know. The source code
of the application should be accessible freely so other users can learn from
your project.

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

    An ORC power plant using two-phase geothermal sources is designed and an
    optimization is carried out. The plant's performance is investigated for
    six different working fluids. Gross and net power output are optimized
    The open source library pygmo :cite:`Biscani2020` is applied in combination
    with TESPy for the first time, setting the foundation for the optimization
    API of TESPy (see
    :py:class:`tespy.tools.optimization.OptimizationProblem`).

    +++
    Title: Parametric optimization and comparative study of an organic Rankine
    cycle power plant for two-phase geothermal sources

    Authors: Chaofan Chen, Francesco Witte, Ilja Tuschy, Olaf Kolditz, Haibing
    Shao

    Reference: :cite:`Chen2022`


.. card::

    **Combined Heat and Power Organic Rankine Cycle**
    ^^^

    .. grid:: 2

        .. grid-item::

            Starting from well production information for a geothermal energy
            reservoir over a lifetime of 40 years, the development of the
            electrical power output of an ORC is monitored within different
            designs of the plant. The geothermal heat source is exploted to
            provide heat to a district heating system and the residual heat is
            used to operate the orc cycle.

        .. grid-item::

            .. image:: /_static/images/examples/GRC_electrical_power_output.svg
              :align: center
              :alt: Development of the Electrical Power Output of the ORC for a District with 2 MW Peak Heat Load
              :class: only-light
              :target: https://github.com/fwitte/chp_orc

            .. image:: /_static/images/examples/GRC_electrical_power_output_darkmode.svg
              :align: center
              :alt: Development of the Electrical Power Output of the ORC for a District with 2 MW Peak Heat Load
              :class: only-dark
              :target: https://github.com/fwitte/chp_orc

    +++
    Title: Computational Modeling of Organic Rankine Cycle Combined Heat and
    Power for Sedimentary Geothermal Exploitation

    Authors: Nicholas Fry, Jessica Eagle-Bluestone, Francesco Witte

    Reference:
