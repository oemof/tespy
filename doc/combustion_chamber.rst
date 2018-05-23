~~~~~~~~~~~~~~~~~~~~~~~~~~~
Combustion Chamber Tutorial
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. contents::
    :depth: 1
    :local:
    :backlinks: top
	
The combustion chamber is an important component within thermal power plants, but unfortunately is the reason for many issues, as the solving algorithm is very sensitive to small changes e. g.
the fluid composition. We will demonstrate how to handle the combustion chamber in a very small, simple example. You can download the full code from the `github repository <https://github.com/oemof/tespy/tree/master/examples>`_.

First of all you need to define the network containing all fluid components used for the combustion chamber. **These are at least the fuel, oxygen, carbon-dioxide and water**. For this example we added Argon, and of course - as we are using Air for the combustion - Nitrogen.
On top, it is recommended to specify reasonable ranges for the fluid properties.

.. code-block:: python

	from tespy import con, cmp, nwk

	# %% network

	# define full fluid list for the network's variable space
	fluid_list = ['Ar', 'N2', 'O2', 'CO2', 'CH4', 'H2O']
	# define unit systems and fluid property ranges
	nw = nwk.network(fluids=fluid_list, p_unit='bar', T_unit='C',
					 p_range=[0.5, 10], T_range=[10, 1200])
					 
As components there are two sources required, one for the fresh air, one for the fuel, a sink for the flue gas and the combustion chamber. Connect the components and add the connections to your network afterwards.

.. code-block:: python

	# %% components

	# sinks & sources
	amb = cmp.source('ambient')
	sf = cmp.source('fuel')
	fg = cmp.sink('flue gas outlet')

	# combustion chamber
	comb=cmp.combustion_chamber(label='combustion chamber')

	# %% connections

	amb_comb = con.connection(amb, 'out1', comb, 'in1')
	sf_comb = con.connection(sf, 'out1', comb, 'in2')
	comb_fg = con.connection(comb, 'out1', fg, 'in1')

	nw.add_conns(sf_comb, amb_comb, comb_fg)
	
For the parametrisation specify the combustions chambers fuel (we chose methane for this example, for all available fuels check for the :py:class:`combustion chamber in the API <tespy.components.components.combustion_chamber>`.), the air to stoichmetric air ratio lamb and the thermal input (:math:`LHV \cdot \dot{m}_{f}`).

.. code-block:: python

	# %% component parameters

	# set combustion chamber fuel, air to stoichometric air ratio and thermal input
	comb.set_attr(fuel='CH4', lamb=3, ti=20000)
	
The ambient conditions as well as the fuel gas inlet temperature are defined in the next step. The air and the fuel gas composition must fully be stated, the component combustion chamber can not handle "Air" as input fluid!

.. code-block:: python

	# %% connection parameters
								 
	# air from abient (ambient pressure and temperature), air composition must be
	# stated component wise.
	amb_comb.set_attr(p=1, T=20,
					  fluid={'Ar': 0.0129, 'N2': 0.7553, 'H2O': 0,
							 'CH4': 0, 'CO2': 0.0004, 'O2': 0.2314})

	# fuel, pressure must not be stated, as pressure is the same at all inlets and
	# outlets of the combustion chamber
	sf_comb.set_attr(T=25,
					 fluid={'CO2': 0.04, 'Ar': 0, 'N2': 0,
							'O2': 0, 'H2O': 0, 'CH4': 0.96})
							
Finally run the code:

.. code-block:: python

	# %% solving

	mode = 'design'
	nw.solve(mode=mode)
	nw.print_results()
	nw.save('combustion')
	
Of course, you can change the parametrisation in any desired way. For example instead of stating the thermal input, you could choose any of the mass flows, or instead of the air to stoichometric air ratio you could specify the flue gas temperature.
It is also possible to make modifications on the fluid composition, for example stating the oxygen content of the flue gas.
