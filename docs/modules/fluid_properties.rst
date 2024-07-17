.. _tespy_fluid_properties_label:

Fluid properties
================
The default fluid property engine `CoolProp <http://coolprop.org/>`_. All
available fluids can be found on their homepage. Also see :cite:`Bell2014`.
Since version 0.7 of TESPy it is possible to use other engines. TESPy comes with
two additional predefined engines, i.e.

- the `iapws <https://github.com/jjgomera/iapws/>`_ library and
- the `pyromat <http://pyromat.org>`_ library.

For each fluid you can specify, which library should be used, and you can easily
implement your own engine, for example, if your fluid is not available through
the predefined engines, or you need a different implementation of the fluid
properties of a specific fluid. Further below on this page we have added an
example implementation of the KKH polynomial formulation to illustrate, how you
can integrate your own engine into TESPy.

CoolProp
--------
CoolProp is the default fluid property back end of TESPy. It provides different
back ends for different applications and allows calling mixture functions.

Available back ends
+++++++++++++++++++
CoolProp provides multiple back ends for fluid property calculation. The
back ends vary in calculation speed and calculation accuracy. It is possible
to choose from the following back ends:

- :code:`HEOS`: Helmhotz Equation Of State with highest accuracy and lowest
  calculation speed. **This is the default back end!**
- :code:`REFPROP`: Highest accuracy and highest convergence stability.
  **This back end is not free**, a separate
  `REFPROP <https://www.nist.gov/srd/refprop>`__ license is required.
- :code:`BICUBIC`: Tabular back end with high accuracy and very high
  calculation speed.
- :code:`TTSE`: Tabular back end with lowest accuracy and very high calculation
  speed.
- :code:`INCOMP`: Back end for incompressible fluids.
- :code:`IF97`: Back end for the IAPWS-IF97 of water, very accurate and much
  higher calculation speed than :code:`HEOS`. Due to a bug in the CoolProp
  back end this option is instable, for more information see the
  `CoolProp issue #1918 <https://github.com/CoolProp/CoolProp/issues/1918/>`_.

For more information on the back ends please visit the CoolProp online
documentation.

Pure and pseudo-pure fluids
+++++++++++++++++++++++++++
CoolProp covers the most important fluids such as water, air as a pseudo-pure
fluid as well as its components, several fuels and refrigerants etc.. Look for
the aliases in the list of
`fluids <http://www.coolprop.org/fluid_properties/PurePseudoPure.html>`__.
All fluids provided in this list cover liquid and gaseous state and the
two-phase region.

Incompressible fluids
+++++++++++++++++++++
If you are looking for heat transfer fluids, the list of incompressible
`fluids <http://www.coolprop.org/fluid_properties/Incompressibles.html>`__
might be interesting for you. In contrast to the pure fluids, the properties
cover liquid state only.

Fluid mixtures
++++++++++++++
TESPy provides support for three types of mixtures:

- ideal: Mixtures for gases only.
- ideal-cond: Mixture for gases with condensation calculation for water share.
- incompressible: Mixtures for incompressible fluids.

Furthermore, CoolProp provides a back end for predefined mixtures, which is
rather instable using HEOS. Using the CoolProp mixture back-end is not tested,
reach out if you would like to support us in adopting the TESPy implementation.
In general, to use the mixture feature of CoolProp we recommend using the
REFPROP back end instead of HEOS.

Using other engines
-------------------
To use any of the other fluid property engines, you can do the following, e.g.
to use the `iapws` back end:

.. code-block:: python

    >>> from tespy.components import Sink
    >>> from tespy.components import Source
    >>> from tespy.components import Turbine
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> from tespy.tools.fluid_properties.wrappers import IAPWSWrapper

    >>> nwk = Network(iterinfo=False)

    >>> so = Source("Source")
    >>> tu = Turbine("Turbine")
    >>> si = Sink("Sink")

    >>> c1 = Connection(so, "out1", tu, "in1", label="1")
    >>> c2 = Connection(tu, "out1", si, "in1", label="2")

    >>> nwk.add_conns(c1, c2)

    >>> tu.set_attr(eta_s=0.9)

    >>> c1.set_attr(
    ...     v=1, p=1e5, T=500,
    ...     fluid={"IF97::H2O": 1}, fluid_engines={"H2O": IAPWSWrapper}
    ... )
    >>> c2.set_attr(p=1e4)

    >>> nwk.solve("design")
    >>> round(c2.x.val, 3)
    0.99

    >>> tu.set_attr(eta_s=None)
    >>> c2.set_attr(x=1)

    >>> nwk.solve("design")
    >>> round(tu.eta_s.val, 3)
    np.float64(0.841)


Implementing a custom engine
----------------------------
The fluid property calls to different engines have to be masqueraded with
respective wrappers. The implementation of the wrappers for `CoolProp`,
`iapws` and `pyromat` can be found in the
:py:mod:`fluid_properties.wrappers <tespy.tools.fluid_properties.wrappers>`
module, and serve as example implementations for your own wrappers:

The wrapper for your own engine (or an engine from a different library) has to
inherit from the
:py:class:`FluidPropertyWrapper <tespy.tools.fluid_properties.wrappers.FluidPropertyWrapper>`
class. Below we will use the polynomial formulation for **gaseous water** from
:cite:`Knacke1991` as an example. First we import the necessary dependencies.

.. code-block:: python

    >>> import numpy as np
    >>> from tespy.tools.fluid_properties.wrappers import FluidPropertyWrapper
    >>> from tespy.tools.global_vars import gas_constants

Then we set up a new class and implement the methods to calculate enthalpy and
entropy from (pressure and) temperature. The structure and names of the
functions have to match the pattern from the `FluidPropertyWrapper`, in this
case `h_pT`. On top of that, we add a backwards function `T_ph` and a function
to analytically calculate the heat capacity `_cp_pT`, the derivative of the
enthalpy to the temperature. Lastly, to make the calculation of isentropic
efficiencies possible, we can add the equation for change in enthalpy on
isentropic change of pressure for an ideal gas.


.. code-block:: python

    # coefficients       H+       S+       a      b       c    d        M
    >>> COEF = {
    ...    "H2O": [-253.871, -11.750, 34.376, 7.841, -0.423, 0.0, 18.0152],
    ... }

    >>> class KKHWrapper(FluidPropertyWrapper):
    ...
    ...     def __init__(self, fluid, back_end=None, reference_temperature=298.15) -> None:
    ...         super().__init__(fluid, back_end)
    ...
    ...         if self.fluid not in COEF:
    ...             msg = "Fluid not available in KKH database"
    ...             raise KeyError(msg)
    ...
    ...         self.coefficients = COEF[fluid]
    ...         self.h_ref = self._h_pT(None, reference_temperature)
    ...         self._molar_mass = self.coefficients[-1] * 1e-3
    ...         self._T_min = 100
    ...         self._T_max = 2000
    ...         self._p_min = 1000
    ...         self._p_max = 10000000
    ...
    ...     def cp_pT(self, p, T):
    ...         y = T * 1e-3
    ...         return 1e3 * (
    ...             self.coefficients[2]
    ...             + self.coefficients[3] * y
    ...             + self.coefficients[4] / (y ** 2)
    ...             + self.coefficients[5] * y ** 2
    ...         ) / self.coefficients[6]
    ...
    ...     def h_pT(self, p, T):
    ...         return self._h_pT(p, T) - self.h_ref
    ...
    ...     def _h_pT(self, p, T):
    ...         y = T * 1e-3
    ...         return 1e6 * (
    ...             self.coefficients[0]
    ...             + self.coefficients[2] * y
    ...             + self.coefficients[3] / 2 * y ** 2
    ...             - self.coefficients[4] / y
    ...             + self.coefficients[5] / 3 * y ** 3
    ...         ) / self.coefficients[6]
    ...
    ...     def T_ph(self, p, h):
    ...         return newton(self.h_pT, self.cp_pT, h, p)
    ...
    ...     def isentropic(self, p_1, h_1, p_2):
    ...         T_1 = self.T_ph(p_1, h_1)
    ...         cp = self.cp_pT(p_1, T_1)
    ...         kappa = cp / (cp - gas_constants["uni"] / self._molar_mass)
    ...         T_2 = T_1 * (p_2 / p_1) ** ((kappa - 1) / kappa)
    ...         return self.h_pT(p_2, T_2)

We can add a newton for the backwards function:

.. code-block:: python

    >>> def newton(func, deriv, h, p):
    ...     # default valaues
    ...     T = 300
    ...     valmin = 70
    ...     valmax = 3000
    ...     max_iter = 10
    ...     tol_rel = 1e-6
    ...
    ...     # start newton loop
    ...     expr = True
    ...     i = 0
    ...     while expr:
    ...         # calculate function residual and new value
    ...         res = h - func(p, T)
    ...         T += res / deriv(p, T)
    ...
    ...         # check for value ranges
    ...         if T < valmin:
    ...             T = valmin
    ...         if T > valmax:
    ...             T = valmax
    ...         i += 1
    ...
    ...         if i > max_iter:
    ...             break
    ...
    ...         expr = abs(res / h) >= tol_rel
    ...
    ...     return T

And then we can test a call to the interface and check our results, also compare
them to CoolProp.

.. code-block:: python

    >>> kkh_water = KKHWrapper("H2O", reference_temperature=298.15)  # same as in CoolProp
    >>> h = kkh_water.h_pT(1e5, 400)
    >>> T = kkh_water.T_ph(1e5, h)
    >>> round(T, 1)
    400.0
    >>> round(h)
    189769

    >>> from tespy.tools.fluid_properties import CoolPropWrapper

    >>> coolprop_water = CoolPropWrapper("H2O")
    >>> h_cp = coolprop_water.h_pT(1e5, 400)
    >>> T_cp = coolprop_water.T_ph(1e5, h_cp)
    >>> round(T_cp, 1)
    400.0
    >>> round(h)
    189769

To use this wrapper in a simple TESPy model, we can then proceed as we have in
the previous section:

.. code-block:: python


    >>> from tespy.components import Sink
    >>> from tespy.components import Source
    >>> from tespy.components import Turbine
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network

    >>> nwk = Network(T_unit="C", p_unit="MPa", iterinfo=False)

    >>> so = Source("Source")
    >>> tu = Turbine("Turbine")
    >>> si = Sink("Sink")

    >>> c1 = Connection(so, "out1", tu, "in1", label="1")
    >>> c2 = Connection(tu, "out1", si, "in1", label="2")

    >>> nwk.add_conns(c1, c2)

    >>> c1.set_attr(
    ...     m=1, p=10, T=600,
    ...     fluid={"H2O": 1}, fluid_engines={"H2O": KKHWrapper}
    ... )
    >>> c2.set_attr(p=1, T=400)

    >>> nwk.solve("design")

    >>> tu.set_attr(eta_s=0.9)
    >>> c2.set_attr(T=None)
    >>> nwk.solve("design")
    >>> round(c2.T.val, 1)
    306.3

Mixture routines in TESPy
-------------------------
Different types of mixture routines are implemented in TESPy. You can select,
which routine should be applied in each separated subnetwork of your system by
specifying a mixing rule. `ideal-cond` is the default mixing rule. The following
mixing rules are available at the moment:

- `ideal-cond`: gaseous fluids **with flash calculations for water**.
- `ideal`: gaseous fluids **without flash calculations**.
- `incompressible`: mass based mixtures of individual incompressible fluids.

The mixtures are calculated by using the pure fluid properties from the selected
fluid property engines and combining them through corresponding equations. The
equations are documented in the
:py:mod:`fluid_properties.mixtures <tespy.tools.fluid_properties.mixtures>`
module.

.. note::

    Similarly to the custom fluid property engine, you can implement your own
    mixture routines. If you are interested in doing so, you can get in contact
    via the :ref:`user meeting <tespy_community_label>` or the GitHub
    `discussion forum <https://github.com/oemof/tespy/discussions>`__.

.. _FluProDia_label:

Creating Fluid Property Diagrams
--------------------------------

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
       * - class ORCEvaporator
         - :code:`in1`
         - :code:`out1`
         - :code:`1`
       * -
         - :code:`in2`
         - :code:`out2`
         - :code:`2`
       * -
         - :code:`in3`
         - :code:`out3`
         - :code:`3`
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
