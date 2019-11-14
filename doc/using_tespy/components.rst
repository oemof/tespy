.. _using_tespy_components_label:

TESPy components
================

In this section we will introduce you into the details of component parametrisation and component characteristics. At the end of the section we show you, how to create custom components.

List of components
------------------

More information on the components can be gathered from the code documentation. We have linked the base class containing a figure and basic informations as well as the equations.

- :py:class:`Source <tespy.components.components.source>` (no equations)
- :py:class:`Sink <tespy.components.components.sink>` (no equations)
- Nodes (base class is node)
	- :py:class:`Node <tespy.components.components.node>` (:py:meth:`equations <tespy.components.components.node.equations>`)
	- :py:class:`Merge <tespy.components.components.merge>` (:py:meth:`equations <tespy.components.components.node.equations>`)
	- :py:class:`Splitter <tespy.components.components.splitter>` (:py:meth:`equations <tespy.components.components.node.equations>`)
	- :py:class:`Separator <tespy.components.components.separator>` (:py:meth:`equations <tespy.components.components.node.equations>`)
- :py:class:`Valve <tespy.components.components.valve>` (:py:meth:`equations <tespy.components.components.valve.equations>`)
- Turbomachines (base class is turbomachine)
	* :py:class:`Pump <tespy.components.components.pump>` (:py:meth:`equations <tespy.components.components.turbomachine.equations>`)
	* :py:class:`Compressor <tespy.components.components.compressor>` (:py:meth:`equations <tespy.components.components.turbomachine.equations>`)
	* :py:class:`Turbine <tespy.components.components.turbine>` (:py:meth:`equations <tespy.components.components.turbomachine.equations>`)
- Components with combustion (base class is combustion_chamber)
	* :py:class:`Combustion chamber <tespy.components.components.combustion_chamber>` (:py:meth:`equations <tespy.components.components.combustion_chamber.equations>`)
	* :py:class:`Combustion chamber stoichiometric <tespy.components.components.combustion_chamber_stoich>` (:py:meth:`equations <tespy.components.components.combustion_chamber_stoich.equations>`)
	* :py:class:`Cogeneration unit <tespy.components.components.cogeneration_unit>` (:py:meth:`equations <tespy.components.components.cogeneration_unit.equations>`)
- Heat exchangers (base class is heat_exchanger)
	* :py:class:`Heat exchanger <tespy.components.components.heat_exchanger>` (:py:meth:`equations <tespy.components.components.heat_exchanger.equations>`)
	* :py:class:`Condenser <tespy.components.components.condenser>` (:py:meth:`equations <tespy.components.components.heat_exchanger.equations>`)
	* :py:class:`Desuperheater <tespy.components.components.desuperheater>` (:py:meth:`equations <tespy.components.components.heat_exchanger.equations>`)
- Simplified heat exchangers (base class is heat_exchanger_simple)
	* :py:class:`Heat exchanger simple <tespy.components.components.heat_exchanger_simple>` (:py:meth:`equations <tespy.components.components.heat_exchanger_simple.equations>`)
	* :py:class:`Pipe <tespy.components.components.pipe>` (:py:meth:`equations <tespy.components.components.heat_exchanger_simple.equations>`)
	* :py:class:`Solar collector <tespy.components.components.solar_collector>` (:py:meth:`equations <tespy.components.components.heat_exchanger_simple.equations>`)
- :py:class:`Drum <tespy.components.components.drum>` (:py:meth:`equations <tespy.components.components.drum.equations>`)
- :py:class:`Water electrolyzer <tespy.components.components.water_electrolyzer>` (:py:meth:`equations <tespy.components.components.water_electrolyzer.equations>`)
- :py:class:`Subsystem interface <tespy.components.components.subsys_interface>` (:py:meth:`equations <tespy.components.components.subsys_interface.equations>`)

.. _using_tespy_components_parametrisation_label:

Component parametrisation
-------------------------

Component parameters can be set and accessed in various ways. All parameters of components are objects of a :code:`data_container` class. The data container for component parameters it is called :code:`dc_cp`, :code:`dc_cc` for component characteristics and :code:`dc_cm` for characteristic maps.
The main purpose of having a data container for the parameters (instead of pure numbers), is added flexibility for the user.

There are different ways for you to specify a component parameter, we use a heat exchanger as an example.

Parameters
^^^^^^^^^^

.. code-block:: python

	from tespy import cmp, hlp
	import numpy as np

	he = cmp.heat_exchanger('evaporator')

	# ways to specify (and set) value
	he.set_attr(kA=1e5)
	# specify data container (same result as above)
	he.set_attr(kA=hlp.dc_cp(val=1e5, is_set=True))

	# ways to unset value
	he.set_attr(kA=np.nan)
	he.kA.set_attr(is_set=False)

	# custom variables
	pipe = cmp.pipe('my pipe')

	# make diameter variable of system
	pipe.set_attr(D='var')
	# data container specification with identical result,
	# benefit: val is the starting value in this case
	pipe.set_attr(D=hlp.dc_cp(val=0.2, is_set=True, is_var=True))

	# data container specification with identical result,
	# benefit: specification of bounds will increase stability
	pipe.set_attr(D=hlp.dc_cp(val=0.2, is_set=True, is_var=True, min_val=0.1, max_val=0.3))


Characteristics
^^^^^^^^^^^^^^^

.. code-block:: python

	from tespy import cmp, hlp
	import numpy as np

	he = cmp.heat_exchanger('evaporator')

	# specify name of predefined method
	he.set_attr(kA_char1='EVA_HOT')
	he.set_attr(kA_char2='EVA_COLD')

	# specify data container (yields same result)
	he.set_attr(kA_char1=hlp.dc_cc(method='EVA_HOT', param='m'))

	# specify data container (custom interpolation points x and y)
	x = np.array([0, 0.5, 1, 2])
	y = np.array([0, 0.8, 1, 1.2])
	he.set_attr(kA_char1=hlp.dc_cc(param='m', x=x, y=y))

.. _component_characteristics_label:

Component characteristics
-------------------------

Characteristics are available for the following components and parameters:

- pump
	* :py:meth:`eta_s_char <tespy.components.components.pump.eta_s_char_func>`: isentropic efficiency vs. volumetric flow rate.
	* :py:meth:`flow_char <tespy.components.components.pump.flow_char_func>`: pressure rise vs. volumetric flow characteristic.
- compressor
	* :py:meth:`char_map <tespy.components.components.compressor.char_map_func>`: component map for isentropic efficiency and pressure rise.
	* :py:meth:`eta_s_char <tespy.components.components.compressor.eta_s_char_func>`: isentropic efficiency vs. pressure ratio.
- turbine
	* :py:meth:`eta_s_char <tespy.components.components.turbine.eta_s_char_func>`: isentropic efficiency vs. isentropic enthalpy difference/pressure ratio/volumetric flow/mass flow.
- heat exchangers:
	* :py:meth:`kA1_char, kA2_char <tespy.components.components.heat_exchanger.kA_func>`: heat transfer coefficient, various predefined types, mass flows as specification parameters.
- simple heat exchangers
	* :py:meth:`kA_char <tespy.components.components.heat_exchanger_simple.kA_func>`: e. g. pipe, see heat exchangers
- cogeneration unit
	* :py:meth:`tiP_char <tespy.components.components.cogeneration_unit.tiP_char_func>`: thermal input vs. power ratio.
	* :py:meth:`Q1_char <tespy.components.components.cogeneration_unit.Q1_char_func>`: heat output 1 vs. power ratio.
	* :py:meth:`Q2_char <tespy.components.components.cogeneration_unit.Q2_char_func>`: heat output 2 vs. power ratio.
	* :py:meth:`Qloss_char <tespy.components.components.cogeneration_unit.Qloss_char_func>`: heat loss vs. power ratio.

You can specify the name of a default characteristic line or you define the whole data container for this parameter. The default characteristic lines can be found in the :py:mod:`documentation <tespy.components.characteristics>`.

.. code-block:: python

	from tespy import cmp, hlp

	turb = cmp.turbine('turbine')
	# method specification (default characteristic line "TRAUPEL")
	turb.set_attr(eta_s_char='TRAUPEL')
	# data container specification
	turb.set_attr(eta_s_char=hlp.dc_cc(method='TRAUPEL', param='dh_s', x=None, y=None))

	# defining a custom line (this line overrides the default characteristic line, method does not need to be specified)
	x = np.array([0, 1, 2])
	y = np.array([0.95, 1, 0.95])
	turb.set_attr(eta_s_char=hlp.dc_cc(param='dh_s', x=x, y=y)

	# heat exchanger analogously
	he = cmp.heat_exchanger('evaporator')
	he.set_attr(kA_char1='EVA_HOT')
	he.set_attr(kA_char2='EVA_COLD')

Custom components
-----------------

If required, you can add custom components. These components should inherit from :py:class:`tespy.components.components.component class <tespy.components.components.component>` or its children.
In order to do that, create a python file in your working directory and import the :py:mod:`tespy.components.components module <tespy.components.components>`. The most important methods are

- :code:`attr(self)`,
- :code:`inlets(self)`,
- :code:`outlets(self)`,
- :code:`equations(self)`,
- :code:`derivatives(self, nw)` and
- :code:`calc_parameters(self, nw, mode)`,

where :code:`nw` is a :py:class:`tespy.networks.network object <tespy.networks.network>`.

The starting lines of your file would look like this:

.. code:: python

	from tespy import cmp


	class my_custom_component(cmp.component):


Attributes
^^^^^^^^^^

The attr method returns a dictionary with the attributes you are able to specify when you want to parametrize your component as keys. The values for each key are the type of data_container this parameter should hold.

.. code:: python

	def attr(self):
		return {'par1': dc_cp(), 'par2': dc_cc()}


Inlets and outlets
^^^^^^^^^^^^^^^^^^

:code:`inlets(self)` and :code:`outlets(self)` respectively must return a list of strings. The list may look like this:

.. code:: python

	def inlets(self):
		return ['in1', 'in2']

	def outlets(self):
		return ['out1', 'out2']

The number of inlets and outlets might even be generic, e. g. if you have added an attribute :code:`'num_in'` in :code:`attr(self)`:

.. code:: python

    def inlets(self):
        if self.num_in_set:
            return ['in' + str(i + 1) for i in range(self.num_in)]
        else:
            self.set_attr(num_in=2)
            return self.inlets()

Equations
^^^^^^^^^

The equations contain the information on the changes to the fluid properties within the component. Each equation must be defined in a way, that the correct result is zero, e. g.:

.. math::

	0 = \dot{m}_{in} - \dot{m}_{out}\\
	0 = \dot{p}_{in} - \dot{p}_{out} - \Delta p

The connections connected to your component are available as a list in :code:`self.inl` and :code:`self.outl` respectively.

.. code:: python

    def equations(self):

    	vec_res = []

		vec_res += [self.inl[0].m.val_SI - self.outl[0].m.val_SI]
		vec_res += [self.inl[0].p.val_SI - self.outl[0].p.val_SI - self.dp.val]

The equations are added to a list one after another, which will be returned at the end.

Derivatives
^^^^^^^^^^^

You need to calculate the partial derivatives of the equations to all variables of the network.
This means, that you have to calculate the partial derivatives to mass flow, pressure, enthalpy and all fluids in the fluid vector on each incomming or outgoing connection of the component.

Add all derivatives to a list (in the same order as the equations) and return the list as numpy array (:code:`np.asarray(list)`).
The derivatives can be calculated analytically or numerically by using the inbuilt function :code:`numeric_deriv(self, func, dx, pos, **kwargs)`.

- :code:`func` is the function you want to calculate the derivatives for,
- :code:`dx` is the variable you want to calculate the derivative to and
- :code:`pos` indicates the connection you want to calculate the derivative for, e. g. :code:`pos=1` means, that counting your inlets and outlets from low index to high index (first inlets, then outlets),
  the connection to be used is the second connection in that list.
- :code:`kwargs` are additional keyword arguments required for the function.

For a good start just look into the source code of the inbuilt components. If you have further questions feel free to contact us.
