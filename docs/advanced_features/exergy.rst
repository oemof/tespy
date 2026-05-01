.. _advanced_exergy_label:

~~~~~~~~~~~~~~~
Exergy analysis
~~~~~~~~~~~~~~~

.. note::

    The exergy analysis described on this page uses
    `exerpy <https://github.com/oemof/exerpy>`__, the dedicated external
    library for exergy analysis. Exerpy is fully compatible with TESPy models -
    build your network as usual and pass the system boundary crossing streams to
    exerpy. Beyond physical exergy, exerpy also supports chemical exergy and
    exergoeconomic methods. You can also use the inbuilt integration of the
    :ref:`ModelTemplate <integration_model_class_template_label>` class. A
    complete tutorial on exergy analysis with the ModelTemplate is available in
    the heat pump
    :ref:`exergy analysis tutorial <tutorial_heat_pump_exergy_label>`. More
    examples are in the
    `exerpy documentation <https://exerpy.readthedocs.io/en/latest/examples.html>`__.

Performing thermodynamic cycle analyses making use of the second law of
thermodynamics provides further process information and uncovers potentials for
improvement in power plant engineering. Therefore, TESPy's analyses module
provides you with an inbuilt and fully automatic exergy analysis.

We have published a paper with the features described in this section. The
publication is licensed under an open-access license, download the pdf
at https://doi.org/10.3390/en15114087, also see :cite:`Witte2022`.

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
      - :math:`e^\mathrm{PH}`, :math:`E^\mathrm{PH}`
      - due to the deviation of the temperature and pressure of the system from
        those of the environment
    * - :code:`ex_therm`, :code:`Ex_therm`
      - (specific) thermal exergy
      - :math:`e^\mathrm{T}`, :math:`E^\mathrm{T}`
      - associated with the system temperature
    * - :code:`ex_mech`, :code:`Ex_mech`
      - (specific) mechanical exergy
      - :math:`e^\mathrm{M}`, :math:`E^\mathrm{M}`
      - associated with the system pressure
    * - :code:`ex_chemical`, :code:`Ex_chemical`
      - (specific) chemical exergy
      - :math:`e^\mathrm{CH}`, :math:`E^\mathrm{CH}`
      - based on standard chemical exergy in ambient model, the `tespy.data`
        module provides three different datasets for standard exergy based on
        various sources, i.e. `Ahrendts`
        :cite:`Ahrendts1980,Ahrendts1977,Ahrendts1974`, `Szargut1988`
        :cite:`Szargut1988` and `Szargut2007` :cite:`Szargut2007,Bakshi2011`.
    * - :code:`E_P`
      - product exergy
      - :math:`\dot{E}_\mathrm{P}`
      - represents the desired result(expressed in terms of exergy) generated
        by the system being considered represents the resources (expressed in
        terms of exergy)
    * - :code:`E_F`
      - fuel exergy
      - :math:`\dot{E}_\mathrm{F}`
      - represents the resources (expressed in terms of exergy) expended to
        provide the product exergy
    * - :code:`E_D`
      - exergy destruction
      - :math:`\dot{E}_\mathrm{D}`
      - thermodynamic inefficiencies associated with the irreversibility
        (entropy generation) within the system boundaries
    * - :code:`E_L`
      - exergy loss
      - :math:`\dot{E}_\mathrm{L}`
      - thermodynamic inefficiencies associated with the transfer of exergy
        through material and energy streams to the surroundings
    * - :code:`epsilon`
      - exergetic efficiency
      - :math:`\varepsilon`
      - ratio between product exergy and fuel exergy
    * - :code:`y_D,k`
      - exergy destruction ratio
      - :math:`y_\mathrm{D}`
      - rate of exergy destruction in a component compared to the exergy rate
        of the fuel provided to the overall system
    * - :code:`y*_D,k`
      - exergy destruction ratio
      - :math:`y^*_\mathrm{D}`
      - rate of exergy destruction in a component compared to the total exergy
        destruction rate within the system

.. note::

    The generic exergy analysis balance equations have not yet been fully
    implemented and tested for the components `FuelCell`, `WaterElectrolzer`
    and `CombustionEngine`.

Tutorial
========
In this short tutorial, an exergy analysis is carried out for the so-called
"Solar Energy Generating System" (SEGS). The full python script is available on
GitHub in an individual repository: https://github.com/fwitte/SEGS_exergy.

.. tip::

  Two other full code examples are to be found at:

  - Supercritical CO\ :sub:`2` power cycle: https://github.com/fwitte/sCO2_exergy
  - Refrigeration machine: https://github.com/fwitte/refrigeration_cycle_exergy

SEGS consists of three main systems, the solar field, the steam cycle and the
cooling water system. In the solar field Therminol VP1 (TVP1) is used as heat
transfer fluid. In the steam generator and reheater the TVP1 is cooled down to
evaporate and overheat/reheat the water of the steam cycle. The turbine is
divided in a high pressure turbine and a low pressure turbine, which are
further subdivided in 2 parts (high pressure turbine) and 5 parts. In between
the stages steam is extracted for preheating. Finally, the main condenser of
the steam cycle is connected to an air cooling tower. The figure below shows
the topology of the model.

.. figure:: /_static/images/advanced/exergy/flowsheet.svg
    :align: center
    :alt: Topology of the Solar Energy Generating System (SEGS)
    :figclass: only-light

.. figure:: /_static/images/advanced/exergy/flowsheet_darkmode.svg
    :align: center
    :alt: Topology of the Solar Energy Generating System (SEGS)
    :figclass: only-dark

The input data are based on literature :cite:`Kearney1988`, which provides
measured data. Some parameters are however taken from a follow-up publication,
as the original data show some inconsistencies, e.g. higher enthalpy at the low
pressure turbine's last stage outlet than at its inlet :cite:`Lippke1995`. As
mentioned, you can find all data in the respective GitHub repository.

TESPy model
-----------
The TESPy model consists of 53 components. The feed water tank serves as mixing
preheater, thus can be modeled using a merge. All other components are modeled
highlighted in the flowsheet. The preheaters and the main condenser are modeled
as :code:`Condenser` instances, while all other heat exchangers are modeled
using :code:`HeatExchanger` instances. For the solar field a parabolic trough
is implemented, calculating the surface area required for the provision of the
heat input at optimal conditions.

All components are flagged with the :code:`fkt_group` parameter, which will
automatically create functional groups (component groups) for the exergy
analysis Grassmann diagram. The specification of this parameter is not required
for the exergy analysis itself, but helps to simplify the automatically
generated diagram. Components not assigned to any functional group will form
their respective group.

Regarding parameter specification, the following parameters are specified:

- isentropic efficiency values
- electrical conversion efficiencies of motors and generators
- terminal temperature difference values at preheaters
- pressure values of steam extraction
- pressure values in the preheating route
- pressure losses in the heat exchangers
- solar fluid temperature
- steam cycle live steam and reheat temperatures
- some temperature values in the cooling water system

The ambient state is defined as follows:

.. code-block:: python

    pamb = 1.013
    Tamb = 25

Pressure and temperature of the ambient air in the cooling tower are equal to
these values in the script provided.

For the exact values of the component parameters please see in the referenced
python script.

Due to the complexity of the plant, the solver sometimes struggles when given bad
starting values. Therefore, the TESPy model is built in two steps. After
solving the initial setup without both of the high pressure preheater
subcoolers, the missing connections and components are added in a second step
and the model is again solved.

Analysis setup
--------------
After the simulation of the plant, the exergy analysis can be carried out.
To perform it, all exergy streams leaving or entering the network's system
boundaries have to be defined by the user. These are:

- fuel exergy :code:`E_F`
- product exergy :code:`E_P`
- exergy loss streams :code:`E_L`
- internal exergy streams not bound to connections :code:`internal_busses`

In case of the solar thermal power plant, the fuel exergy is the heat input at
the parabolic trough. The product is the electricity produced by the plant,
i.e. the electricity generated by the turbine generators minus the electricity
consumed by the pumps and the fan. Lastly, exergy loss streams are the hot air
leaving the cooling tower and the cold air entering the cooling tower fan from
the ambient.

In exerpy, the fuel, product and loss are defined as dictionaries referencing
connection labels at the system boundary. The connection labels must be set
explicitly when building the network. For example, if the solar heat input
connection is labelled :code:`'heat_in'`, the grid power output connection
:code:`'power_out'`, and the cooling tower air connections :code:`'air_in'`
and :code:`'air_out'`:

.. code-block:: python

    from exerpy import ExergyAnalysis

    # exerpy requires SI units: Tamb in K, pamb in Pa
    ean = ExergyAnalysis.from_tespy(SEGSvi, Tamb + 273.15, pamb * 1e5)

    fuel = {"inputs": ["heat_in"], "outputs": []}
    product = {"inputs": ["power_out"], "outputs": []}
    loss = {"inputs": ["air_out"], "outputs": ["air_in"]}

    ean.analyse(E_F=fuel, E_P=product, E_L=loss)

The :code:`E_F`, :code:`E_P`, :code:`E_L` and :code:`epsilon` results for the
overall system are available as attributes after calling :code:`analyse`.
Component-level results are returned by :code:`exergy_results()`.

.. code-block:: python

    df_comp, df_material, df_power = ean.exergy_results()
    print(f"epsilon = {ean.epsilon:.3f}")

To repeat the analysis with a different ambient state, create a new instance
from the same (re-solved) network, e.g. for :code:`Tamb = 15 °C`:

.. code-block:: python

    ean2 = ExergyAnalysis.from_tespy(SEGSvi, 15 + 273.15, pamb * 1e5)
    ean2.analyse(E_F=fuel, E_P=product, E_L=loss)


.. note::

    If the network's topology changed a new :code:`ExergyAnalysis` instance
    must be created via :code:`from_tespy`.

Checking consistency
--------------------
An automatic check of consistency is performed by the analysis. The sum of all
exergy destruction values of the network's components is calculated. On top of
that, fuel, product exergy and exergy loss are determined. The total exergy
destruction must therefore be equal to the fuel exergy minus product exergy and
minus exergy loss. The deviation of that equation is then calculated and
checked versus a threshold value of :math:`10^{-3}` (to compensate for
rounding errors).

.. math::

    \dot{E}_\mathrm{D} = \dot{E}_\mathrm{F} - \dot{E}_\mathrm{P} -
    \dot{E}_\mathrm{L}

    \Delta \dot{E} = \dot{E}_\mathrm{F} - \dot{E}_\mathrm{P} -
    \dot{E}_\mathrm{L} - \dot{E}_\mathrm{D}

    \Delta \dot{E} \leq 10^{-3}

.. note::

    If the exergy analysis is carried out on a converged simulation and the
    analysis is set up correctly, this equation must be True. Otherwise, an
    error will be printed to the console, which means:

    - The simulation of your plant did not converge or
    - the exergy analysis has not been set up correctly. You should
      check, if the definition of the exergy streams :code:`E_F`, :code:`E_P`,
      :code:`E_L` and :code:`internal_busses` is correct.

    If you suspect a bug in the calculation, you are welcome to submit an issue
    on our GitHub page.

Printing and accessing the results is done with the
:py:meth:`exerpy.ExergyAnalysis.exergy_results` method, which returns three
pandas DataFrames and optionally prints tabular summaries to the console.

.. code-block:: python

    df_comp, df_material, df_power = ean.exergy_results()

The DataFrames contain:

- :code:`df_comp` - component-level results (:code:`E_F [kW]`, :code:`E_P [kW]`,
  :code:`E_D [kW]`, :code:`epsilon [%]`, :code:`y [%]`, :code:`y* [%]`), plus a
  :code:`"TOT"` row with overall system totals.
- :code:`df_material` - exergy values for material (fluid) connections.
- :code:`df_power` - exergy values for power connections.

A concise text summary of the system-level results is printed by:

.. code-block:: python

    ean.print_exergy_summary()

The overall system quantities are also available directly as attributes:

.. code-block:: python

    print(ean.E_F)      # total fuel exergy in W
    print(ean.E_P)      # total product exergy in W
    print(ean.E_D)      # total exergy destruction in W
    print(ean.E_L)      # total exergy loss in W
    print(ean.epsilon)  # overall exergetic efficiency

Plotting
--------
An exergy destruction waterfall diagram is generated with:

.. code-block:: python

    ean.plot_exergy_waterfall(title='SEGS Exergy Analysis')

.. figure:: /_static/images/advanced/exergy/sankey.svg
    :align: center
    :alt: Sankey diagram of the Solar Energy Generating System (SEGS)

The results can be exported to JSON for further processing or exchange with
other tools:

.. code-block:: python

    ean.export_to_json('segs_exergy_results.json')


Conclusion
==========
Further examples using TESPy with exerpy are available in the
`exerpy documentation <https://exerpy.readthedocs.io/en/latest/examples.html>`__.
Exerpy supports physical and chemical exergy as well as exergoeconomic methods.
If you are interested in contributing to exerpy, please file an issue at the
`exerpy GitHub page <https://github.com/oemof/exerpy>`__.
