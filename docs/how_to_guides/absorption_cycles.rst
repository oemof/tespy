.. _tutorial_absorption_cycles_label:

########################################
Model absorption chillers (experimental)
########################################

.. warning::

    The absorption cycle modeling is experimental and not thoroughly tested.
    Results and the API may change in a future release. Feedback is very
    welcome, please open an issue on
    `GitHub <https://github.com/oemof/tespy/issues>`__.

Absorption chillers and heat pumps replace the compressor of a vapor
compression cycle by a solution cycle: the refrigerant vapor is absorbed in
a solution at low pressure, the solution is pumped to high pressure and the
refrigerant is driven out by heating the solution. TESPy models solutions with
:py:class:`SolutionConnection <tespy.connections.solutionconnection.SolutionConnection>`
and unit processes with

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

- :code:`"libr-water"`: LiBr-water after Pátek and Klomfar :cite:`Patek2006`,
  with the reference state of IAPWS-95.
- :code:`"nh3-water"`: NH3-water after the IAPWS guideline G4-01
  :cite:`TillnerRoth1998`, evaluated with the compiled Helmholtz energy kernel
  of teqp :cite:`Bell2022`.

The convenience parameter :code:`xi` of a solution connection is the mass
fraction LiBr for LiBr-water and NH3 for NH3-water. A solution connection
reports the mass fraction :code:`xi`, the vapor mass fraction :code:`x` of a
two-phase solution, the bubble point temperature :code:`T_bubble` and the
subcooling :code:`td_bubble`. :code:`x` and :code:`td_bubble` can be specified
like on a pure fluid connection.

Custom working pairs
--------------------

A working pair is a subclass of
:py:class:`SolutionPair <tespy.tools.fluid_properties.solutions.base.SolutionPair>`.
It implements the state functions in terms of pressure, temperature or enthalpy
and the mass fraction of the key component and is registered with
:py:func:`register_solution_pair <tespy.tools.fluid_properties.solutions.base.register_solution_pair>`,
which makes it available as mixing rule of the same name. A custom pair with
the same key component as a built-in pair must be named with
:code:`mixing_rule` on a solution connection of each circuit, it cannot be
inferred from the fluid names. See the
`custom_working_pair.py <https://github.com/fwitte/tespy-sorption-models/blob/main/custom_working_pair.py>`__
example.
