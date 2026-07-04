.. _advanced_exergy_label:

~~~~~~~~~~~~~~~
Exergy analysis
~~~~~~~~~~~~~~~

Performing thermodynamic cycle analyses making use of the second law of
thermodynamics provides further process information and uncovers potentials for
improvement in power plant engineering. Therefore, this feature provides you
with an inbuilt and fully automatic exergy analysis.

We have published a paper with the features described in this section. The
publication is licensed under an open-access license, download the pdf
at https://doi.org/10.3390/en15114087, also see :cite:`Witte2022`.

.. note::

    The exergy analysis described on this page uses
    `exerpy <https://github.com/oemof/exerpy>`__, a dedicated external library
    for exergy analysis. The implementations are based on the original features
    developed and included in TESPy. These features have been removed from
    TESPy now. Exerpy is fully compatible with TESPy models. You can build your
    network as usual and pass the system boundary crossing streams to exerpy.
    Beyond physical exergy, exerpy also supports chemical exergy and
    exergoeconomic methods. You can also use the inbuilt integration of the
    :ref:`ModelTemplate <integration_model_class_template_label>` class. A
    complete tutorial on exergy analysis with the ModelTemplate is available in
    the heat pump
    :ref:`exergy analysis tutorial <tutorial_heat_pump_exergy_label>`. More
    examples are in the
    `exerpy documentation <https://exerpy.readthedocs.io/en/latest/examples.html>`__.

Fundamentals of exergy analysis
===============================
Energy is a concept of the first law of thermodynamics. It cannot be destroyed.
But regarding the design and analysis of thermal systems, the idea that
something can be destroyed is useful. According to the second law of
thermodynamics, the conversion of heat and internal energy into work is
limited. This constraint and the idea of destruction are applied to introduce a
new concept: "Exergy".

Exergy can be destroyed due to irreversibility and is able to describe the
quality of different energy forms. The difference in quality of different forms
of energy shall be illustrated by the following example. 1 kJ of electrical
energy is clearly more valuable than 1 kJ of energy in a glass of water at
ambient temperature :cite:`Bejan1996`.

In literature, exergy is defined as follows:

    *"An opportunity for doing useful work exists whenever two systems at
    different states are placed in communication, for in principle work can be
    developed as the two are allowed to come into equilibrium. When one of the
    two systems is a suitably idealized system called an environment and the
    other is some system of interest, exergy is the maximum theoretical useful
    work (shaft work or electrical work) obtainable as the systems interact to
    equilibrium, heat transfer occurring with the environment only."*
    :cite:`Bejan1996`

Terminology
===========
The definitions and nomenclature of the exergy analysis in TESPy are based on
:cite:`Tsatsaronis2007`. The exergy destruction ratios are described in more
detail in :cite:`Bejan1996`. Since the current version of the exergy analysis
in TESPy only focuses on physical exergy and does not include reaction
processes yet, chemical exergy is not considered. Changes in kinetic and
potential exergy are neglected and therefore not considered as well.

.. list-table:: Terminology
    :widths: 20 20 10 50
    :header-rows: 1
    :class: tight-table

    * - variable
      - name
      - symbol
      - description
    * - :code:`ex_physical`, :code:`Ex_physical`
      - (specific) physical exergy
      - :math:`e^\text{PH}`, :math:`E^\text{PH}`
      - due to the deviation of the temperature and pressure of the system from
        those of the environment
    * - :code:`ex_therm`, :code:`Ex_therm`
      - (specific) thermal exergy
      - :math:`e^\text{T}`, :math:`E^\text{T}`
      - associated with the system temperature
    * - :code:`ex_mech`, :code:`Ex_mech`
      - (specific) mechanical exergy
      - :math:`e^\text{M}`, :math:`E^\text{M}`
      - associated with the system pressure
    * - :code:`ex_chemical`, :code:`Ex_chemical`
      - (specific) chemical exergy
      - :math:`e^\text{CH}`, :math:`E^\text{CH}`
      - based on standard chemical exergy in ambient model, the `tespy.data`
        module provides three different datasets for standard exergy based on
        various sources, i.e. `Ahrendts`
        :cite:`Ahrendts1980,Ahrendts1977,Ahrendts1974`, `Szargut1988`
        :cite:`Szargut1988` and `Szargut2007` :cite:`Szargut2007,Bakshi2011`.
    * - :code:`E_P`
      - product exergy
      - :math:`\dot{E}_\text{P}`
      - represents the desired result(expressed in terms of exergy) generated
        by the system being considered represents the resources (expressed in
        terms of exergy)
    * - :code:`E_F`
      - fuel exergy
      - :math:`\dot{E}_\text{F}`
      - represents the resources (expressed in terms of exergy) expended to
        provide the product exergy
    * - :code:`E_D`
      - exergy destruction
      - :math:`\dot{E}_\text{D}`
      - thermodynamic inefficiencies associated with the irreversibility
        (entropy generation) within the system boundaries
    * - :code:`E_L`
      - exergy loss
      - :math:`\dot{E}_\text{L}`
      - thermodynamic inefficiencies associated with the transfer of exergy
        through material and energy streams to the surroundings
    * - :code:`epsilon`
      - exergetic efficiency
      - :math:`\varepsilon`
      - ratio between product exergy and fuel exergy
    * - :code:`y_D,k`
      - exergy destruction ratio
      - :math:`y_\text{D}`
      - rate of exergy destruction in a component compared to the exergy rate
        of the fuel provided to the overall system
    * - :code:`y*_D,k`
      - exergy destruction ratio
      - :math:`y^*_\text{D}`
      - rate of exergy destruction in a component compared to the total exergy
        destruction rate within the system
