.. _tutorial_absorption_cycles_label:

#########################
Model absorption chillers
#########################

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

Single-effect LiBr-water chiller
--------------------------------

The cycle consists of the generator (desorber), condenser, evaporator,
absorber, solution pump, solution heat exchanger and the two valves. The
solution circuit is closed by a :code:`CycleCloser`.

.. literalinclude:: /../tutorial/advanced/absorption_chiller.py
    :language: python
    :start-after: [sec_1]
    :end-before: [sec_2]

The working pair is named on the connection leaving the absorber, all other
solution connections take it over.

.. literalinclude:: /../tutorial/advanced/absorption_chiller.py
    :language: python
    :start-after: [sec_2]
    :end-before: [sec_3]

The temperatures of the absorber and generator outlets, of the condenser and
of the evaporator define the pressures and the equilibrium compositions; the
cooling capacity defines the mass flows. The cycle closer does not close the
variable composition of the solution, therefore :code:`fluid_balance=True` is
specified on the connection leaving it.

.. literalinclude:: /../tutorial/advanced/absorption_chiller.py
    :language: python
    :start-after: [sec_3]
    :end-before: [sec_4]

The vapor leaving the generator is in equilibrium with the entering solution
by default (:code:`vapor_state="inlet-equilibrium"`), as in a counter-current
generator. :code:`"outlet-equilibrium"` assumes equilibrium with the leaving
solution, :code:`"free"` leaves the vapor state to a specification on the
vapor connection. An approach to the bubble point of the solution leaving
the absorber or desorber can be specified with :code:`td_eq`.

.. literalinclude:: /../tutorial/advanced/absorption_chiller.py
    :language: python
    :start-after: [sec_4]
    :end-before: [sec_5]

NH3-water chiller with rectifier
--------------------------------

With NH3-water the vapor leaving the generator contains water. A
:code:`Rectifier` purifies it by partial condensation: the purified vapor
leaves at its dew point at the composition specified on its outlet, the reflux
is the liquid in equilibrium with the entering vapor. The reflux is merged
with the solution entering the generator, e.g. with a :code:`Merge`; the
rectifier does not impose its pressure on the reflux, the merge does. All
streams of the cycle are solution connections, the condenser outlet is
specified with :code:`x=0` and the evaporator outlet with the vapor mass
fraction reached in evaporation, since the refrigerant evaporates with a
temperature glide. The test
:code:`tests/test_models/test_absorption_chillers.py` holds the complete
model of the cycle of the CoolSolve Library (CSL-0249).

Custom working pairs
--------------------

A working pair is a subclass of
:py:class:`SolutionPair <tespy.tools.fluid_properties.solutions.base.SolutionPair>`.
It implements the state functions in terms of pressure, temperature or
enthalpy and the mass fraction of the key component and is registered with
:py:func:`register_solution_pair <tespy.tools.fluid_properties.solutions.base.register_solution_pair>`,
which makes it available as mixing rule of the same name. The example derives
a pair from LiBr-water with a modified boiling point elevation.

.. literalinclude:: /../tutorial/advanced/absorption_chiller.py
    :language: python
    :start-after: [sec_5]
    :end-before: [sec_6]
