Combustion Chamber Tutorial
---------------------------

.. contents::
    :depth: 1
    :local:
    :backlinks: top
	
There are two different types of combustion chambers available. The combustion 
chamber can handle varying fluid compositions for the air and the fuel and
calculates the fluid composition of the flue gas. Thus, it is possible to e. g.
specify the oxygen mass fraction in the flue gas in a calculation. In
contrast, the stoichiometric combustion chamber uses fuel and air as pseudo
pure gases for the input at calculates a mixture of stoichiometric flue gas
and air at the outlet. The sacrifice of flexibility for parametrisation results
in a faster solution process. Thus, if the air composition and the fuel
composition are known prior to calculation, it is always recommended to use the
stoichiometric combustion chamber. We provide a tutorial for both components,
where you learn how they work, and what the differences are.
	
Combustion chamber
^^^^^^^^^^^^^^^^^^
	
The combustion chamber is an important component within thermal power plants,
but unfortunately is the reason for many issues, as the solving algorithm is
very sensitive to small changes e. g. the fluid composition. We will
demonstrate how to handle the combustion chamber in a very small, simple
example. You can download the full code from the `tespy_examples repository
<https://github.com/oemof/oemof-examples/tree/master/oemof_examples/tespy/combustion_chamber>`_.

First of all you need to define the network containing all fluid components
used for the combustion chamber. **These are at least the fuel, oxygen,
carbon-dioxide and water**. For this example we added Argon, and of course - as
we are using Air for the combustion - Nitrogen.
On top, it is highly recommended to specify reasonable ranges for the fluid
properties. If you have fluid mixtures within your network, the fluid property
ranges will keep your fluid properties within the specified ranges for the 
first three steps of the newton's algorithm in order to stabilize the 
calculation.

.. code-block:: python

    from tespy.networks import network

    # define full fluid list for the network's variable space
    fluid_list = ['Ar', 'N2', 'O2', 'CO2', 'CH4', 'H2O']

    # define unit systems and fluid property ranges
    nw = network(fluids=fluid_list, p_unit='bar', T_unit='C',
                 p_range=[0.5, 10], T_range=[10, 1200])
					 
As components there are two sources required, one for the fresh air, one for
the fuel, a sink for the flue gas and the combustion chamber. Connect the
components and add the connections to your network afterwards.

.. code-block:: python

    from tespy.components import sink, source, combustion_chamber

    # sinks & sources
    amb = source('ambient')
    sf = source('fuel')
    fg = sink('flue gas outlet')

    # combustion chamber
    comb=combustion_chamber(label='combustion chamber')

.. code-block:: python

    from tespy.connections import connection

    amb_comb = connection(amb, 'out1', comb, 'in1')
    sf_comb = connection(sf, 'out1', comb, 'in2')
    comb_fg = connection(comb, 'out1', fg, 'in1')

    nw.add_conns(sf_comb, amb_comb, comb_fg)
	
For the parametrisation we specify the combustions chambers the air to
stoichmetric air ratio lamb and the thermal input
(:math:`LHV \cdot \dot{m}_{f}`).

.. code-block:: python

    # set combustion chamber air to stoichometric air ratio and thermal input
    comb.set_attr(lamb=3, ti=2e6)
	
The ambient conditions as well as the fuel gas inlet temperature are defined in
the next step. The air and the fuel gas composition must fully be stated, the
component combustion chamber can not handle "Air" as input fluid!

.. code-block:: python

    # air from abient (ambient pressure and temperature), air composition must
    # be stated component wise.
    amb_comb.set_attr(p=1, T=20, fluid={'Ar': 0.0129, 'N2': 0.7553, 'H2O': 0,
                                        'CH4': 0, 'CO2': 0.0004, 'O2': 0.2314})

    # fuel, pressure must not be stated, as pressure is the same at all inlets
    # and outlets of the combustion chamber
    sf_comb.set_attr(T=25, fluid={'CO2': 0.04, 'Ar': 0, 'N2': 0, 'O2': 0,
                                  'H2O': 0, 'CH4': 0.96})
							
Finally run the code:

.. code-block:: python

    nw.solve('design')
    nw.print_results()
	
Of course, you can change the parametrisation in any desired way. For example
instead of stating the thermal input, you could choose any of the mass flows,
or instead of the air to stoichometric air ratio you could specify the flue
gas temperature. It is also possible to make modifications on the fluid
composition, for example stating the oxygen content in the flue gas. It is also
possible to change the fuel composition. Make sure, all desired fuels of your
fuel mixture are also within the fluid_list of the network. For the example
below we added some hydrogen to the fuel mixture.

.. code-block:: python

    from tespy.networks import network
    from tespy.components import sink, source, combustion_chamber
    from tespy.connections import connection

    # %% network

    fluid_list = ['Ar', 'N2', 'O2', 'CO2', 'CH4', 'H2O', 'H2']
    nw = network(fluids=fluid_list, p_unit='bar', T_unit='C',
                                 p_range=[0.5, 10], T_range=[10, 1200])

    # %% components

    # sinks & sources
    amb = source('ambient')
    sf = source('fuel')
    fg = sink('flue gas outlet')

    # combustion chamber
    comb=combustion_chamber(label='combustion chamber')

    # %% connections

    amb_comb = connection(amb, 'out1', comb, 'in1')
    sf_comb = connection(sf, 'out1', comb, 'in2')
    comb_fg = connection(comb, 'out1', fg, 'in1')

    nw.add_conns(sf_comb, amb_comb, comb_fg)

    # %% component parameters

    # set combustion chamber air to stoichometric air ratio and thermal input
    comb.set_attr(lamb=3, ti=2e6)

    # %% connection parameters

    amb_comb.set_attr(p=1, T=20,
                      fluid={'Ar': 0.0129, 'N2': 0.7553, 'H2O': 0, 'CH4': 0,
                             'CO2': 0.0004, 'O2': 0.2314, 'H2': 0.01})

    sf_comb.set_attr(T=25, fluid={'CO2': 0, 'Ar': 0, 'N2': 0,'O2': 0, 'H2O': 0,
                                  'CH4': 0.95, 'H2': 0.05})

    # %% solving

    nw.solve('design')
    nw.print_results()

Stoichiometric combustion chamber
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The example for the stoichiometric combustion chamber can as well be taken from
the `tespy_examples repository
<https://github.com/oemof/oemof-examples/tree/master/oemof_examples/tespy/combustion_chamber>`_.

Again, the network must have the information, which fluids will be part of the
fluid vector. In contrast to the normal combustion chamber, you will need the
following fluids: **Air, Fuel and Flue Gas**. For this tutorial we will call
them: **"TESPy::myAir", "TESPy::myFuel" and "TESPy::myFuel_fg"**, we will see,
why we chose these names for the fluids later. Do not forget to specify the
ranges for pressure and temperature. This is a very important stept for this
specific component, we will explain later, why it is.

.. code-block:: python

    from tespy.networks import network

    # define full fluid list for the network's variable space
    fluid_list = ['TESPy::myAir', 'TESPy::myFuel', 'TESPy::myFuel_fg']

    # define unit systems and fluid property ranges
    nw = network(fluids=fluid_list, p_unit='bar', T_unit='C',
                     p_range=[1, 10], T_range=[10, 2000])

The components required are then the same as in the first tutorial, the
stoichiometric combustion chamber's class is called
"combustion_chamber_stoich". As components there are two sources required, one
for the fresh air, one for the fuel, a sink for the flue gas and the combustion
chamber. Connect the components and add the connections to your network
afterwards.

.. code-block:: python

    from tespy.components import sink, source, combustion_chamber_stoich

    # sinks & sources
    amb = source('ambient')
    sf = source('fuel')
    fg = sink('flue gas outlet')

    # combustion chamber
    comb = combustion_chamber_stoich('stoichiometric combustion chamber')

.. code-block:: python

    from tespy.connections import connection

    amb_comb = connection(amb, 'out1', comb, 'in1')
    sf_comb = connection(sf, 'out1', comb, 'in2')
    comb_fg = connection(comb, 'out1', fg, 'in1')

    nw.add_conns(sf_comb, amb_comb, comb_fg)
	
The basic parametrisation of the stoichiometric combustion chamber is different
compared to the normal combustion chamber: We need to specify the air and the
fuel composition, and additionally, aliases for the these fluids. Since air and
fuel usually are mixtures of different gases, **TESPy will create lookup tables
for the fluid properties of the specified air and fuel composition and a third
lookup table for the flue gas**. TESPy will therefore calculate the
stoichiometric flue gas composition. The fluids will then be accessable with
the following aliases: **"TESPy::youraliasforair", "TESPy::youraliasforfuel"
and "TESPy::youraliasforfuel_fg"**. The creation of the lookup tables will use
your network's settings: **The fluid properties will be calculated within the
network's specified ranges for pressure and temperature.**

A folder called "LUT" will be created in your working directory containing all
fluid property lookup tables. As the creation of the lookup tables does take
some time, it is possible, to read the fluid properties from that folder: You
need to specify the path variable, like this: :code:`path='./LUT'`.

There are some important things to keep in mind, when reading the fluid
properties from path:

- **Do not specify the path in case**

	- you change the pressure range or the temperature range or
	- you change the air or the fuel composition.

- **For convergence stability choose large maximum temperatures**, much higher
than the highest temperature you are expecting at the combustion chambers 
outlet.
- **If you use more than one combustion chamber** do not use identical aliases,
if the fluid compositions are not identical.

As in the example above, we also specify thermal input and lambda, as well as
identical parameters for the connections. Thus the results should be exactly
the same.

.. code-block:: python

    # for the first calculation run
    comb.set_attr(fuel={'CH4': 0.96, 'CO2': 0.04},
                        air={'Ar': 0.0129, 'N2': 0.7553, 'H2O': 0,
                             'CH4': 0, 'CO2': 0.0004, 'O2': 0.2314},
                        fuel_alias='myFuel', air_alias='myAir',
                        lamb=3, ti=20000)

    # if there are existing lookup tables
    comb.set_attr(fuel={'CH4': 0.96, 'CO2': 0.04},
                        air={'Ar': 0.0129, 'N2': 0.7553, 'H2O': 0,
                        'CH4': 0, 'CO2': 0.0004, 'O2': 0.2314},
                         fuel_alias='myFuel', air_alias='myAir', path='./LUT',
                         lamb=3, ti=20000)
				  
.. code-block:: python

    # air from abient (ambient pressure and temperature), air composition must
    # be stated component wise.
    amb_comb.set_attr(T=20, p=1, fluid={'TESPy::myAir': 1, 'TESPy::myFuel': 0,
                                        'TESPy::myFuel_fg': 0})

    # fuel, pressure must not be stated, as pressure is the same at all inlets
    # and outlets of the combustion chamber
    sf_comb.set_attr(T=25, fluid={'TESPy::myAir': 0, 'TESPy::myFuel': 1,
                                  'TESPy::myFuel_fg': 0})
							
Finally run the code:

.. code-block:: python

	# %% solving

	mode = 'design'
	nw.solve(mode=mode)
	nw.print_results()
	nw.save('combustion')
