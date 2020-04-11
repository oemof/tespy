Fluid properties
================
The basic fluid properties are handled by
`CoolProp <http://www.coolprop.org/>`_. All available fluids can be found on
their homepage. Also see :cite:`Bell2014`.

CoolProp back-ends
------------------
CoolProp provides multiple back-ends for fluid property calculation. The
back-ends vary in calculation speed and calculation accuracy. It is possible
to choose from the following back-ends:

- :code:`HEOS`: Helmhotz Equation Of State with highest accuracy and lowest
  calculation speed. **This is the default back-end!**
- :code:`BICUBIC`: Tabular back-end with high accuracy and very high
  calculation speed.
- :code:`TTSE`: Tabular back-end with lowest accuracy and very high calculation
  speed.
- :code:`INCOMP`: Back-end for incompressible fluids.
- ~~:code:`IF97`: Back-end for the IAPWS-IF97 of water, very accurate and much
  higher calculation speed than :code:`HEOS`.~~ Due to a bug in the CoolProp
  back end this option is not available at the moment, for more information
  see this `github issue <https://github.com/CoolProp/CoolProp/issues/1918/>`_.

For more information on the back-ends please visit the CoolProp online
documentation.

Pure and pseudo-pure fluids
---------------------------
If you use pure fluids, TESPy directly uses CoolProp functions to gather all
fluid properties. CoolProp covers the most important fluids such as water, air
as a pseudo-pure fluid as well as its components, several fuels and
refrigerants etc.. Look for the aliases in the list of
`fluids <http://www.coolprop.org/fluid_properties/PurePseudoPure.html>`_.
All fluids provided in this list cover liquid and gaseous state and the
two-phase region.

Incompressible fluids
---------------------
If you are looking for heat transfer fluids, the list of incompressible
`fluids <http://www.coolprop.org/fluid_properties/Incompressibles.html>`_
might be interesting for you. In contrast to the pure fluids, the properties
cover liquid state only.

Fluid mixtures
--------------
CoolProp provides fluid properties for two component mixtures. BUT: These are
NOT integrated in TESPy! Nevertheless, you can use fluid mixtures for gases.

Ideal mixtures of gaseous fluids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
TESPy can handle mixtures of gaseous fluids, by using the single fluid
properties from CoolProp together with corresponding equations for mixtures.
The equations can be found in the
:py:mod:`fluid_properties <tespy.tools.fluid_properties>` module and are
applied automatically to the fluid vector.

It is also possible create lookup-tables for fluid mixtures with fixed mass
fractions of the components, as this reduces the amount of CoolProp fluid
property calls and speeds up your calculation. Look up the
:py:class:`tespy_fluids documentation <tespy.tools.fluid_properties.tespy_fluid>`
for more information.

Other mixtures
^^^^^^^^^^^^^^
It is **not possible** to use mixtures of liquid and other liquid or gaseous
fluids **at the moment**! If you try to use a mixture of two liquid or gaseous
fluids and liquid fluids, e. g. water and methanol or liquid water and air, the
equations will still be applied, but obviously return bad values. If you have
ideas for the implementation of new kinds of mixtures we appreciate you
contacting us.
