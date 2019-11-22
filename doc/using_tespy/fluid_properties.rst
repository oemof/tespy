.. _tespy_fluid_properties_label:

Fluid properties in TESPy
=========================

The basic fluid properties are handled by `CoolProp <http://www.coolprop.org/>`_. All available fluids can be found on their homepage.

Pure and pseudo-pure fluids
---------------------------

If you use pure fluids, TESPy directly uses CoolProp functions to gather all fluid properties.
CoolProp covers the most important fluids such as water, air as a pseudo-pure fluid as well as its components, several fuels and refrigerants etc..
Look for the aliases in the `list of fluids <http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids>`_. All fluids provided in this list cover liquid and gaseous state and the two-phase region.

Incompressible fluids
---------------------

If you are looking for heat transer fluids, the `list of incompressible fluids <http://www.coolprop.org/fluid_properties/Incompressibles.html>`_ might be interesting for you.
In contrast to the pure fluids, the properties cover liquid state only.

Fluid mixtures
--------------

CoolProp provides fluid properties for two component mixtures. BUT: These are NOT integrated in TESPy! Nevertheless, you can use fluid mixtures for gases:

Ideal mixtures of gaseous fluids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TESPy can handle mixtures of gaseous fluids, by using the single fluid properties from CoolProp together with corresponding equations for mixtures.
The equations can be found in the :py:mod:`tespy.tools.helpers module <tespy.tools.helpers>` and are applied automatically to the fluid vector.

It is also possible create lookup-tables for fluid mixtures with fixed mass fractions of the components, as this reduces the amount of CoolProp fluid property calls and speeds up your calculation. Look up the :py:class:`tespy_fluids documentation <tespy.tools.helpers.tespy_fluid>` for more information.

Other mixtures
^^^^^^^^^^^^^^

It is **not possible** to use mixtures of liquid and other liquid or gaseous fluids **at the moment**!
If you try to use a mixture of two liquid or gaseous fluids and liquid fluids, e. g. water and methanol or liquid water and air, the equations will still be applied, but obviously return bad values.
If you have ideas for the implementation of new kinds of mixtures we appreciate you contacting us.
