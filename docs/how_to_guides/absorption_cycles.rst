.. _tutorial_absorption_cycles_label:

########################################
Model absorption chillers (experimental)
########################################

.. warning::

    The absorption cycle modeling is experimental and not thoroughly tested
    yet, results and the API may change. A :code:`UserWarning` reminds of
    this once per session. Feedback is very welcome, please open an issue on
    `GitHub <https://github.com/oemof/tespy/issues>`__.

Absorption chillers and heat pumps replace the compressor of a vapor
compression cycle by a solution circuit: the refrigerant vapor is absorbed in
a solution at low pressure, the solution is pumped to high pressure and the
refrigerant is driven out by heat. TESPy models the solution streams with the
:py:class:`SolutionConnection <tespy.connections.solutionconnection.SolutionConnection>`
and the sorption processes with the components

- :py:class:`Absorber <tespy.components.sorption.absorber.Absorber>` and
  :py:class:`CooledAbsorber <tespy.components.sorption.absorber.CooledAbsorber>`,
- :py:class:`Desorber <tespy.components.sorption.desorber.Desorber>` and
  :py:class:`HeatedDesorber <tespy.components.sorption.desorber.HeatedDesorber>`,
- :py:class:`Rectifier <tespy.components.sorption.rectifier.Rectifier>` for
  working pairs with a volatile absorbent.

Examples
--------

The example models live in the
`tespy-sorption-models <https://github.com/fwitte/tespy-sorption-models>`__
repository, together with instructions to install and run them:

.. list-table::
    :header-rows: 1
    :widths: 35 65

    * - Example
      - Content
    * - `libr_water_absorber.py <https://github.com/fwitte/tespy-sorption-models/blob/main/libr_water_absorber.py>`__
      - :code:`Absorber`: outlet enthalpy, mass flow, temperature and
        composition specifications, approach to equilibrium :code:`td_eq`
    * - `libr_water_absorber_cooled.py <https://github.com/fwitte/tespy-sorption-models/blob/main/libr_water_absorber_cooled.py>`__
      - :code:`CooledAbsorber` with cooling water
    * - `libr_water_desorber.py <https://github.com/fwitte/tespy-sorption-models/blob/main/libr_water_desorber.py>`__
      - :code:`Desorber`: the :code:`vapor_state` options and outlet
        specifications
    * - `libr_water_desorber_heated.py <https://github.com/fwitte/tespy-sorption-models/blob/main/libr_water_desorber_heated.py>`__
      - :code:`HeatedDesorber` with heating water
    * - `libr_water_chiller.py <https://github.com/fwitte/tespy-sorption-models/blob/main/libr_water_chiller.py>`__
      - Closed single-effect LiBr-water chiller with solution heat exchanger
    * - `nh3_water_absorber.py <https://github.com/fwitte/tespy-sorption-models/blob/main/nh3_water_absorber.py>`__
      - :code:`Absorber` with NH3-water
    * - `nh3_water_generator_rectifier.py <https://github.com/fwitte/tespy-sorption-models/blob/main/nh3_water_generator_rectifier.py>`__
      - :code:`Desorber` and :code:`Rectifier`: purification of the generator
        vapor
    * - `nh3_water_chiller.py <https://github.com/fwitte/tespy-sorption-models/blob/main/nh3_water_chiller.py>`__
      - Closed NH3-water chiller with rectifier, solution and refrigerant heat
        exchanger (CoolSolve Library CSL-0249)
    * - `custom_working_pair.py <https://github.com/fwitte/tespy-sorption-models/blob/main/custom_working_pair.py>`__
      - Custom working pair derived from LiBr-water
    * - `working_pair_properties.py <https://github.com/fwitte/tespy-sorption-models/blob/main/working_pair_properties.py>`__
      - Duehring chart of LiBr-water and bubble and dew lines of NH3-water

Working pairs
-------------

The properties of a solution come from a working pair, which TESPy treats as
mixing rule of the solution connections. Two pairs are available:

- :code:`"libr-water"`: LiBr-water after Pátek and Klomfar
  :cite:`Patek2006`, with the reference state of IAPWS-95. The vapor is pure
  water on a plain :code:`Connection`.
- :code:`"nh3-water"`: NH3-water after the IAPWS guideline G4-01
  :cite:`TillnerRoth1998`, evaluated with the compiled Helmholtz energy kernel
  of teqp :cite:`Bell2022`. teqp is an optional dependency, install it with
  :code:`pip install teqp`. Water is volatile, the vapor is a
  :code:`SolutionConnection` of the same pair.

The convenience parameter :code:`xi` of a solution connection is the mass
fraction of the key component of the pair, LiBr for LiBr-water and NH3 for
NH3-water. A solution connection takes its pair from the :code:`mixing_rule`
argument, from the fluid names or from the other connections of the same
fluid branch, it is sufficient to name it on one connection of a circuit.

Besides the usual results, a solution connection reports the mass fraction
:code:`xi`, the vapor mass fraction :code:`x` of a two-phase solution, the
bubble point temperature :code:`T_bubble` and the subcooling
:code:`td_bubble`. :code:`x` and :code:`td_bubble` can be specified like on a
pure fluid connection.

Modeling notes
--------------

- The solution leaving an absorber or desorber is at its bubble point, an
  approach to the bubble point can be specified with :code:`td_eq`.
- The vapor leaving a desorber is in equilibrium with the entering solution
  by default (:code:`vapor_state="inlet-equilibrium"`), as in a
  counter-current generator. :code:`"outlet-equilibrium"` assumes equilibrium
  with the leaving solution, :code:`"free"` leaves the vapor state to a
  specification on the vapor connection.
- A :code:`CycleCloser` does not close the variable composition of a solution
  circuit, specify :code:`fluid_balance=True` on the connection leaving it.
- With NH3-water the vapor leaving the generator contains water. A
  :code:`Rectifier` purifies it by partial condensation: the purified vapor
  leaves at its dew point at the composition specified on its outlet, the
  reflux is the liquid in equilibrium with the entering vapor. The rectifier
  does not impose its pressure on the reflux, merge the reflux with the
  solution entering the generator, e.g. with a :code:`Merge`, or specify its
  pressure.
- NH3-water refrigerant evaporates with a temperature glide, specify the
  evaporator outlet with the vapor mass fraction :code:`x` reached in
  evaporation.

Custom working pairs
--------------------

A working pair is a subclass of
:py:class:`SolutionPair <tespy.tools.fluid_properties.solutions.base.SolutionPair>`.
It implements the state functions in terms of pressure, temperature or
enthalpy and the mass fraction of the key component and is registered with
:py:func:`register_solution_pair <tespy.tools.fluid_properties.solutions.base.register_solution_pair>`,
which makes it available as mixing rule of the same name. A custom pair with
the same key component as a built-in pair must be named with
:code:`mixing_rule` on a solution connection of each circuit, it cannot be
inferred from the fluid names. See the
`custom_working_pair.py <https://github.com/fwitte/tespy-sorption-models/blob/main/custom_working_pair.py>`__
example.
