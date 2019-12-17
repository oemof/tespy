Heat pump tutorial
------------------

.. contents::
    :depth: 1
    :local:
    :backlinks: top
	
Task
^^^^

This tutorial introduces you in how to model a heat pump in TESPy. You can see 
the plants topology in figure 1. Also, you will find a fully working model in 
the last chapter of this tutorial.

.. figure:: api/_images/tutorial_heat_pump.svg
    :align: center
	
    Figure 1: Topology of the heat pump.
	
The main purpose of the heat pump is to deliver heat e. g. for the consumers of 
a heating system. Thus, the heat pump's parameters will be set in a way, which 
supports this target.
Generally, if systems are getting more complex, it is highly recommended to set 
up your plant in incremental steps. This tuturial divides the plant in three 
sections: The consumer part, the valve and the evaporator and the compressor as 
last element. Each new section will be appended to the existing ones.


Set up a Network
^^^^^^^^^^^^^^^^

In order to simulate our heat pump we have to create an instance of the 
tespy.network class. The network is the main container of the model and will be 
required in all following sections.
First, it is necessary to specify a list of the fluids used in the plant. In 
this example we will work with water (H\ :sub:`2`\O) and ammonia (NH\ :sub:`3`\). 
Water is used for the cold side of the heat exchanger, for the consumer and for 
the hot side of the environmental temperature. Ammonia is used as coolant within 
the heat pump circuit.
If you don’t specify the unit system, the variables are set to SI-Units.

.. code-block:: python

    from tespy.networks import network

    nw = network(fluids=['water', 'NH3'],
                 T_unit='C', p_unit='bar', h_unit='kJ / kg', m_unit='kg / s')
	
We will use °C, bar and kJ/kg as units for temperature and enthalpy.
	
Modeling the heat pump: Consumer system
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Components
++++++++++

We will start with the consumer as the plant will be designed to deliver a 
specific heat flow. From figure 1 you can determine the components of the 
consumer system: condenser, pump and the consumer (heat-exchanger-simple). 
Additionally we need a source and a sink for the consumer and the heat pump 
circuit respectively. We label the sink for the coolant "valve", as for our next 
calculation the valve (labeled "valve") will be attached there. In this way, the 
fluid properties can be initialised by csv at the interface-connection, too.

.. code-block:: python
	
    from tespy.components import sink, source, condenser, pump, 
                                 heat_excanger_simple
    
    # sources & sinks

    c_in = source('coolant in')
    cb = source('consumer back flow')
    cf = sink('consumer feed flow')

    va = sink('valve')

    # consumer system

    cd = condenser('condenser')
    rp = pump('recirculation pump')
    cons = heat_exchanger_simple('consumer'))
	
Connections
+++++++++++

In the next steps we will connect the components in order to form a network. 
Every connection requires the source, the source id, the target and the target 
id as arguments: the source is the component from which the connection 
originates, the source id is the outlet id of that component. This applies 
analogously to the target. To find all inlet and outlet ids of a component look 
up the class documentation.

.. code-block:: python

    from tespy.connections import connection
    
    # consumer system

    c_in_cd = connection(c_in, 'out1', cd, 'in1')

    cb_rp = connection(cb, 'out1', rp, 'in1')
    rp_cd = connection(rp, 'out1', cd, 'in2')
    cd_cons = connection(cd, 'out2', cons, 'in1')
    cons_cf = connection(cons, 'out1', cf, 'in1')

    nw.add_conns(c_in_cd, cb_rp, rp_cd, cd_cons, cons_cf)

    # connection condenser - evaporator system

    cd_va = connection(cd, 'out1', va, 'in1')

    nw.add_conns(cd_va)

Parametrization
+++++++++++++++

For the condenser we set pressure ratios on hot and cold side and additionally 
we set a value for the upper terminal temperature difference as design parameter 
and the heat transfer coefficient as offdesign parameter. The consumer will have 
a pressure ratio, too. Further we set the isentropic efficiency for the pump, 
the offdesign efficiency is calculated with a characteristic function. Thus, we 
set the efficiency as design parameter and the characteristic function as 
offdesign parameter. In offdesign calculation the consumer's pressure ratio will 
be a function of the mass flow, thus as offdesign parameter we select zeta. The 
most important parameter is the consumers heat demand. We marked this setting as 
key parameter.

.. code-block:: python

	cd.set_attr(pr1=0.99, pr2=0.99, ttd_u=5, design=['pr2', 'ttd_u'], 
                offdesign=['zeta2', 'kA'])
	rp.set_attr(eta_s=0.8, design=['eta_s'], offdesign=['eta_s_char'])
	cons.set_attr(pr=0.99, design=['pr'], offdesign=['zeta'])
	
.. note::

	In TESPy there are two different types of calculations: design point and 
    offdesign calculation. All parameters specified in the design attribute of a 
    component or connection, will be unset in a offdesign calculation, all 
    parameters specified in the offdesign attribute of a component or connection 
    will be set for the offdesign calculation. The value for these parameters is 
    the value deriven from the design-calculation.
	
	Generally, the design calculation is used for designing your system in the 
    way you want it to look like. This means, that you might want to specify a 
    design point isentropic efficiency, pressure loss or terminal temperature 
    difference.
	After you have designed your system, you are able to make offdesign 
    calculations with TESPy. The offdesign calculation is used to predict the 
    system's behaviour at different points of operation. For this case, this 
    might be different ambient temperature, different feed flow temperature, or 
    partial load.

In order to calculate this network further parametrization is necessary, as e. g. 
the fluids are not determined yet: At the hot inlet of the condensator we define 
the temperature and the fluid vector. In order to fully determine the fluid's 
state at this point, an information on the pressure is required. This is 
archieved by setting the terminal temperature difference (see above). The same 
needs to be done for the consumer cycle. We suggest to set the parameters at the 
pump's inlet. On top, we assume that the consumer requires a constant inlet 
temperature.

The last step is to define the fluid's state after the consumer, this is done 
with references to the pump's inlet, in order to grant that the fluid properties 
at the consumer's outlet are identical to those at the pump's inlet.

.. code-block:: python

    from tespy.connections import ref
    
    c_in_cd.set_attr(T=170, fluid={'water': 0, 'NH3': 1})
    cb_rp.set_attr(T=60, p=10, fluid={'water': 1, 'NH3': 0})
    cd_cons.set_attr(T=90)
    cons_cf.set_attr(h=ref(cb_rp, 1, 0), p=ref(cb_rp, 1, 0))

.. code-block:: python
	
    # %% key parameter
	
    cons.set_attr(Q=-230e3)	

Solve
+++++

After creating the system, we want to solve our network. First, we calculate the 
design case and directly after we can perform the offdesign calculation at a 
different value for our key parameter. For general information on the solving 
process in TESPy and available parameters check the corresponding section in 
:ref:`Using TESPy <using_tespy_networks_label>`.

.. code-block:: python

    nw.solve('design')
    nw.print_results()
    nw.save('condenser')

    cons.set_attr(Q=-200e3)

    nw.solve('offdesign', design_path='condenser')
    nw.print_results()


Valve and evaporator system
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Next we will add the valve and the evaporator system to our existing network.

Components
++++++++++

This part contains of a valve followed by a drum with evaporator in forced flow 
and a superheater. Do not forget to change the old sink labeled "valve" to an 
actual valve and the sink used in the previous calculation will represent the 
first compressor, labeled "compressor 1". Add the following components to the 
script. Additionally add the import of the components with heat_exchanger, drum 
and valve.

.. code-block:: python

    from tespy.components import heat_exchanger, drum, valve
	# sources & sinks
	
    amb_in = source('source ambient')
    amb_out = sink('sink ambient')

    cp1 = sink('compressor 1')


	# evaporator system

    va = valve('valve')
    dr = drum('drum')
    ev = heat_exchanger('evaporator')
    su = heat_exchanger('superheater')
    pu = pump('pump evaporator')

Connections
+++++++++++

As we already redefined our variable "va" to be a valve instead of a sink (see 
above), we do not need any adjustments to the connection between the condenser 
and the former sink "cd_va". The valve connects to the drum at the inlet 'in1'. 
The pump of the forced flow evaporation system connects to the drum's outlet 
'out1', the evaporator's cold side connects to the drum's inlet 'in2' and the 
superheater's cold side connects to the drum's outlet 'out2'. This will add the 
following connections to the model:

.. code-block:: python

    # evaporator system

    va_dr = connection(va, 'out1', dr, 'in1')
    dr_pu = connection(dr, 'out1', pu, 'in1')
    pu_ev = connection(pu, 'out1', ev, 'in2')
    ev_dr = connection(ev, 'out2', dr, 'in2')
    dr_su = connection(dr, 'out2', su, 'in2')

    nw.add_conns(va_dr, dr_pu, pu_ev, ev_dr, dr_su)

    amb_in_su = connection(amb_in, 'out1', su, 'in1')
    su_ev = connection(su, 'out1', ev, 'in1')
    ev_amb_out = connection(ev, 'out1', amb_out, 'in1')

    nw.add_conns(amb_in_su, su_ev, ev_amb_out)

    # connection evaporator system - compressor system

    su_cp1 = connection(su, 'out2', cp1, 'in1')

    nw.add_conns(su_cp1)

Parametrization
+++++++++++++++

Previous parametrization stays untouched. Regarding the evaporator, we specify 
pressure ratios on hot and cold side as well as the lower terminal temperature 
difference. We use the hot side pressure ratio and the lower terminal 
temperature difference as design parameteres and choose zeta as well as the area
independet heat transition coefficient as its offdesign parameters. On top of 
that, the characteristic function of the evaporator should follow the predefined
methods 'EVA_HOT' and 'EVA_COLD'. If you want to learn more about handling 
characteristic functions you should have a glance at the :ref:`TESPy components
section <using_tespy_components_label>`. The superheater will also use the 
pressure ratios on hot and cold side. Further we set a value for the upper 
terminal temperature difference. For the pump we set the isentropic efficiency. 
For offdesign and design parameter specification of these components the same 
logic as for the evaporator and the already existing part of the network is 
applied: The system designer has to answer the question: Which parameters are 
design point parameters and how does the component perform at a different 
operation point.

.. code-block:: python

    # evaporator system

    ev.set_attr(pr1=0.99, pr2=0.99, ttd_l=5,
                kA_char1='EVA_HOT', kA_char2='EVA_COLD',
    			design=['pr1', 'ttd_l'], offdesign=['zeta1', 'kA'])
    su.set_attr(pr1=0.99, pr2=0.99, ttd_u=2, design=['pr1', 'pr2', 'ttd_u'],
                offdesign=['zeta1', 'zeta2', 'kA'])
    pu.set_attr(eta_s=0.8, design=['eta_s'], offdesign=['eta_s_char'])
	
Next step is the connetion parametrization: The pressure in the drum and the 
enthalpy of the wet steam reentering the drum need to be determined. For the 
enthalpy we can specify a reference of the circulating mass flow to the main 
cycle mass flow. The pressure is archieved through the given lower terminal 
temperature difference of the evaporator and its hot side outlet temperature. As
we have specified a terminal temperature difference at the evaporator's cold 
side inlet (:code:`ttd_l`), it might be necessary to state a starting value for
the pressure, as we are near to the two-phase region. On the hot side inlet of 
the superheater we definde the temperature, pressure and the fluid. Since the 
pressure between superheater and first compressor will be a result of the 
pressure losses in the superheater and we set the terminal temperature 
difference there, bad starting values will lead to a linear dependency, as a 
temperature and a pressure are set while the fluid's state could be within the 
two phase region. Thus, we set starting values for pressure and for enthalpy at
this connection, to make sure the starting point is outside of the two phase 
region. At last we have to fully determine the state of the incoming fluid at 
the superheater's hot side. 

.. code-block:: python

    # evaporator system cold side

    pu_ev.set_attr(m=ref(va_dr, 4, 0), p0=5)
    su_cp1.set_attr(p0=5, h0=1700)

    # evaporator system hot side

    amb_in_su.set_attr(T=12, p=1, fluid={'water': 1, 'NH3': 0})
    ev_amb_out.set_attr(T=9))
	
Solve
+++++

Again, you should calculate your network after you added these parts. As we have
already calculated one part of our network, this time we can use the
:code:`init_path` for the design calculation and load the results from the 
previous network. This step is not required, but in larger, more complex 
networks, it might help, to archieve better convergence. For the offdesign 
calculation see part 3.1.4.

	
Compressor system
^^^^^^^^^^^^^^^^^

To complete the heat pump, we will add the compressor system to our existnetwork.

Components
++++++++++

This part contains two compressors with an intercooler between them. Again, add 
the component compressor to tespy.import. The cold side of the intercooler 
requires a source and a sink. Again, remember redefining the former sink "cp1" 
to a compressor and add a sink for the outlet of the coolant after the 
compressor system.

.. code-block:: python

    from tespy.components import compressor
    
    # sources & sinks

    ic_in = source('source intercool')
    ic_out = sink('sink intercool')

    c_out = sink('coolant out')

    # compressor-system

    cp1 = compressor('compressor 1')
    cp2 = compressor('compressor 2')
    he = heat_exchanger('intercooler')

Connections
+++++++++++

As done before, add the new connections to the script. After the second 
compressor we need to install a sink, because closing a circuit will always lead
to linear dependency. Just make sure, the fluid properties at the sink after 
the compressor are identical to the fluid properties at the source connected to
the condenser. Another way of doing this, is adding a merge and a splitter at 
some point of your network. Nevertheless, you will require a sink and a source.

.. code-block:: python

    # compressor-system

    cp1_he = connection(cp1, 'out1', he, 'in1')
    he_cp2 = connection(he, 'out1', cp2, 'in1')
    cp2_c_out = connection(cp2, 'out1', c_out, 'in1')

    ic_in_he = connection(ic_in, 'out1', he, 'in2')
    he_ic_out = connection(he, 'out2', ic_out, 'in1')

    nw.add_conns(cp1_he, he_cp2, ic_in_he, he_ic_out, cp2_c_out)

Parametrization
+++++++++++++++

For the two compressor we defined an isentropic efficency and for the offdesign
calculation a generic characteristic line for the isentropic efficiency will be
applied. The first compressor has a fixed pressure ratio, the seconds compressor
pressure ratio will result from the required pressure at the condenser. The heat
exchanger comes with pressure ratios on both sides. The parametrization of all
other components remains identical.

.. code-block:: python

    cp1.set_attr(eta_s=0.8, design=['eta_s'], offdesign=['eta_s_char'])
    cp2.set_attr(eta_s=0.8, pr=5, design=['eta_s'], offdesign=['eta_s_char'])
    he.set_attr(pr1=0.99, pr2=0.98, design=['pr1', 'pr2'],
                offdesign=['zeta1', 'zeta2', 'kA'])
	
Regarding the connections, on the hot side after the intercooler we set the 
temperature. For the cold side of the heat exchanger we set the temperature, the
pressure and the fluid on the inlet flow, at the outlet we specify the
temperature as a design parameter. In offdesign calculation, this will be a 
result from the given heat transfer coefficient (see parametrisation of 
intercooler, kA is an offdesign parameter). Last, make sure the fluid properties
 after the compressor outlet are identical to those at the condenser inlet using
 the references.

The last step leads to a necessary redefinition of the parametrization of the
existing model: As the enthalpy at the outlet of the second compressor is a 
result of the given pressure ratio and the isentropic efficiency, it is not 
allowed to set the temperature at the condenser's hot inlet anymore. This is 
due to forcing the fluid properties at the compressor's outlet and the 
condenser's hot side inlet to be identical with the references.

.. code-block:: python

    # condenser system

    c_in_cd.set_attr(fluid={'water': 0, 'NH3': 1})

    # compressor-system

    he_cp2.set_attr(T0=40, p0=10, design=['T'])
    ic_in_he.set_attr(p=1, T=20, fluid={'water': 1, 'NH3': 0})
    he_ic_out.set_attr(T=30, design=['T'])
    cp2_c_out.set_attr(p=ref(c_in_cd, 1, 0), h=ref(c_in_cd, 1, 0))

Solve
+++++

Here again, using the saved results from previous calculations is always 
favourable, but with the manually adjusted starting values, the calculation 
should still converge. Also see section 3.2.4. If you want to use the previous 
part to initialise start the solver with

.. code-block:: python

	nw.solve('design', init_path='condenser')


Further tasks
^^^^^^^^^^^^^

After successfully modeling the heat pump in design and offdesign cases, you can
now start using your model for further calculations. E. g., if you have a time 
series of required heat flow of your consumer, you can loop over the series and 
perform offdesign calculation adjusting the heat flow every time. Of course, 
this is possible with every offdesign parameter. We provide the scripts after 
each of the three steps of the tutorial:
:download:`Step 1 <../tutorial/step_1.py>`, 
:download:`Step 2 <../tutorial/step_2.py>`, 
:download:`Step 3 <../tutorial/step_3.py>`.

Have fun working with TESPy!
