.. _roadmap_label:
.. _roadmap_timeline_label:

###################
Development roadmap
###################

TESPy development is heading for version 1.0. The page shows the most important
upcoming and planned changes and features. The release infos are split into two
subsections: **Framework & API** which refer to the simulation back-end and the
API of the software, and **New features & components** which are additions of
new features to extend modeling capabilities. Some things are not yet assigned
and collected in the **Under discussion** section.

.. tip::

   This roadmap is an intent and not a promise and written in stone. You can
   support the development or suggest own features or ideas. Open an issue or
   a pull request! Also see :ref:`how to contribute <developing_label>`.

.. rst-class:: roadmap-timeline

- .. rubric:: v0.12
     :class: roadmap-release

  .. rubric:: Framework & API

  - Backtracking line search solver strategy :issue:`986`
  - Custom per variable bounds :issue:`1072`
  - Refactoring of the data containers with clearly separated responsibilities
    and documented attributes for connections and components :issue:`1019`
    :issue:`383`
  - Version stamp and metadata block in all files tespy saves and exports
    :issue:`1069`
  - Some smaller fixes, e.g. silently ignored specifications of non-input
    properties :issue:`1070`
  - Deprecations of renamed parameters (:issue:`972`) and unit specification
    through the :code:`Network` :issue:`750`
  - Unified unit handling with a single source of truth per network
    :issue:`1071`
  - Thermopack as alternative property back end :issue:`624`
  - Connections through named ports, e.g.
    :code:`Connection(turbine.outlet, condenser.hot_inlet)` and
    :code:`component.hot_inlet` instead of :code:`inl[0]` :issue:`689`
  - Deprecate Python 3.11

  .. container:: roadmap-additions

     .. rubric:: New features & components

     - Zoned crossflow heat exchanger for fin-tube condensers and evaporators
       :issue:`1073`

- .. rubric:: v0.13
     :class: roadmap-release

  .. rubric:: Framework & API

  - Stable, versioned model schema for tools building on TESPy :issue:`1074`
  - Explicitly defined public namespace and a naming consistency deprecation
    :issue:`665`
  - Turn parts of a network on and off, e.g. solve a single subsystem or
    deactivate a branch :issue:`600`
  - Retrieve mass flow and fluid branches from the structure matrix
    :issue:`845`
  - Unit agnostic polynomial coefficients and reference states for the
    :code:`PolynomialCompressor` :issue:`1017`
  - Application of Dulmage-Mendelsohn deconmposition also for variable fluid
    problems :issue:`1075`
  - Incremental re-solve: only recalculate the equations affected by changed
    specifications :issue:`667` :issue:`1076`
  - Testable characteristic function evaluation limits, e.g. for :code:`Motor`
    and :code:`Generator` :issue:`700`
  - Well defined specification of fluid composition with partially unknown
    species :issue:`464`

  .. container:: roadmap-additions

     .. rubric:: New features & components

     - Humidity separation for humid air connections :issue:`1077`
     - Absorber and desorber components for absorption heat pumps and
       refrigeration machines :issue:`1016`
     - Variable frequency operation for the :code:`TurboCompressor` and
       :code:`PolynomialCompressor` :issue:`1078`
     - Allow entering the envelope space of the :code:`PolynomialCompressor`
       :issue:`784`

- .. rubric:: Under discussion
     :class: roadmap-release

  .. rubric:: Framework & API

  - Solving of isolated parallel networks :issue:`612`
  - Track which components and connections belong to which subsystem, including
    nested subsystems :issue:`1079`
  - Global and local user configuration files :issue:`1052`
  - Parallelization of component equations :issue:`611`
  - Revised sign convention of the vapor content property :issue:`709`

  .. container:: roadmap-additions

     .. rubric:: New features & components

     - Extended humid air capabilities to include freezing :issue:`1080`
     - Alpha scaling for the :code:`SectionedHeatExchanger` offdesign
       :issue:`1037`
     - Part load efficiency models for compressors :issue:`585`
     - Displacement and volumetric efficiency based compressor equation without
       frequency :issue:`981`

- .. rubric:: v1.0
     :class: roadmap-release

  .. rubric:: Goal: Long term stable API

  - Internal and public API is stable, all deprecations are executed and
    migration scripts/guides are available
  - The model serialization schema is declared stable
