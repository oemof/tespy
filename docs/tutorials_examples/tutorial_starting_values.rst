Starting value tutorial
-----------------------

.. contents::
    :depth: 1
    :local:
    :backlinks: top

Subjekt
^^^^^^^

In numerical and iterative methods, a start value is a certain value of a
variable with which a calculation is started. With more complex TESPy models
it can happen that the simulation does not convert. The numerics of the solver
is vulnerable if the specified variables are not primary variables. These are
variables that can be solved after one iteration. For example, pressures,
temperatures and enthalpies are among to them. Because of that it is
recommended to give these parameter for your simulation.

Task
^^^^

After learning to create simple heat pump with TESPy, the goal is to model
more complex systems. After the successful creation of the network,
components, connections and parameterization, it is often the problem that the
simulation does not convert.

.. error::
    Singularity in jacobian matrix, calculation aborted! Make sure
    your network does not have any linear dependencies in the parametrisation.
    Other reasons might be

    -> given temperature with given pressure in two phase region, try setting
    enthalpy instead or provide accurate starting value for pressure.

    -> given logarithmic temperature differences or kA-values for heat
    exchangers.

    -> support better starting values.

    -> bad starting value for fuel mass flow of combustion chamber, provide
    small (near to zero, but not zero) starting value.

To fix this error, it is recommended to create a stable calculation first.
This is used to have suitable starting values for the actual calculation.

Application
^^^^^^^^^^^

Following the first tutorial a heat pumps with internal heat exchangers is
considered. You can see the plant topology in the figure.

.. figure:: api/_images/tutorial_sv_heat_pump_intheatex.svg
    :align: center

    Figure: Topology of heat pump with internal heat exchanger

It consists of a consumer system, a valve, an evaporator system, a compressor
and additionally of an internal heat exchanger.
In order to simulate this heat pump, the TESPy model has to be built up.
First, the network has to be initialized and the refrigerants used have to be
specified. In this example the heat pump work with ammonia (NH\ :sub:`3`\) and
water (H\ :sub:`2`\O). Furthermore, it is recommended to define the unit
system not to work with variables set to SI-Units.

.. code-block:: python

    from tespy.networks import Network

    # network
    nw = Network(
        fluids=['water', 'NH3'],
        T_unit='C', p_unit='bar', h_unit='kJ / kg', m_unit='kg / s'
        )

After that the required components have to be created.

.. code-block:: python

    from tespy.components import (
        Condenser, Compressor, CycleCloser,  HeatExchanger,
        HeatExchangerSimple, Pump, Sink, Source, Valve
        )

    # components
    cycle_closer = CycleCloser('Refrigerant Cycle Closer')

    # heat source
    heatsource_feedflow = Source('Heat Source Feed Flow')
    heatsource_pump = Pump('Heat Source Recirculation Pump')
    heatsource_evaporator = HeatExchanger('Heat Source Evaporator')
    heatsource_backflow = Sink('Heat Source Back Flow')

    # compression
    compressor = Compressor('Compressor')

    # heat sink
    cons_pump = Pump('Heat Sink Recirculation Pump')
    condenser = Condenser('Heat Sink Condenser')
    cons_heatsink = HeatExchangerSimple('Heat Consumer')
    cons_cycle_closer = CycleCloser('Consumer Feed Flow')

    # internal heat exchange
    int_heatex = HeatExchanger('Internal Heat Exchanger')

    # expansion
    valve = Valve('Expansion Valve')

.. note::

    If the heat pump operates in a supercritical range, the condenser has to
    be replaced with a heat exchanger.

Now the connections according to the topology have to be linked. For a better
overview of the results it is recommended to label the connections.

.. code-block:: python

    from tespy.connections import Connection

    # connections
    # heat source
    cc2hs_eva = Connection(cycle_closer, 'out1', heatsource_evaporator, 'in2', label='cc2hs_eva')
    hs_feed2hs_pump = Connection(heatsource_feedflow, 'out1', heatsource_pump, 'in1', label='hs_feed2hs_pump')
    hs_pump2hs_eva = Connection(heatsource_pump, 'out1', heatsource_evaporator, 'in1', label='hs_pump2hs_eva')
    hs_eva2hs_back = Connection(
        heatsource_evaporator, 'out1', heatsource_backflow, 'in1', label='hs_eva2hs_back'
        )

    nw.add_conns(cc2hs_eva, hs_feed2hs_pump, hs_pump2hs_eva, hs_eva2hs_back)

    # internal heat exchange cold side
    hs_eva2int_heatex = Connection(
        heatsource_evaporator, 'out2', int_heatex, 'in2', label='hs_eva2int_heatex'
        )

    nw.add_conns(hs_eva2int_heatex)

    # compression
    int_heatex2comp = Connection(int_heatex, 'out2', compressor, 'in1', label='int_heatex2comp')
    comp2cond = Connection(compressor, 'out1', condenser, 'in1', label='comp2cond')

    nw.add_conns(int_heatex2comp, comp2cond)

    # heat sink
    cons_back2cons_pump = Connection(
    cons_cycle_closer, 'out1', cons_pump, 'in1', label='cons_back2cons_pump'
        )
    cons_pump2cond = Connection(cons_pump, 'out1', condenser, 'in2', label='cons_pump2cond')
    cond2cons_hs = Connection(condenser, 'out2', cons_heatsink, 'in1', label='cond2cons_hs')
    cons_hs2cons_feed = Connection(
        cons_heatsink, 'out1', cons_cycle_closer, 'in1', label='cons_hs2cons_feed'
        )

    nw.add_conns(cons_back2cons_pump, cons_pump2cond, cond2cons_hs, cons_hs2cons_feed)

    # internal heat exchange hot side
    cond2int_heatex = Connection(condenser, 'out1', int_heatex, 'in1', label='cond2int_heatex')

    nw.add_conns(cond2int_heatex)

    # expansion
    int_heatex2valve = Connection(int_heatex, 'out1', valve, 'in1', label='int_heatex2valve')
    valve2cc = Connection(valve, 'out1', cycle_closer, 'in1', label='valve2cc')

    nw.add_conns(int_heatex2valve, valve2cc)

After the initialization of the network and the creation of the components and
connections, a stable parameterization is built up to have suitable initial
values for the actual simulation.

.. note::

    To create a stable simulation, it is recommended to set pressure and
    enthalpie values instead of temperature values. In this example, fixed
    points can be identified with the help of the logph diagram which you can
    see in the figure.

    On the one hand the point behind the evaporator is fixed. At this point
    the vapor content of the ammonia is at 100% (x=1). Furthermore, it is
    recommended to specify the pressure in order to clearly determine the
    point. On the other hand the point behind the condenser is fixed, too.
    At these point the ammonia has a vapor content of 0% (x=0). As before, the
    pressure value has also to be set.

.. figure:: api/_images/tutorial_sv_logph.svg
    :align: center

    Figure: Logph diagram of ammonia

In addition to the fixed evaporation and condensation points, the fluids to be
used, the feedflow and backflow temperatures of the consumer and heat source
as well as the enthalpy between internal heat exchanger and valve have to be
defined.

To correctly determine the enthalpies and pressures, CoolProp is to be
imported. It is important to note that the PropertySI function (PropsSI) works
with SI unit. These may differ from the units defined in the network.

.. code-block:: python

    import CoolProp.CoolProp as CP

    # parametrization connections
    # set feedflow and backflow temperature of heat source and consumer
    T_hs_bf = 5
    T_hs_ff = 10
    T_cons_bf = 50
    T_cons_ff = 90

    # evaporation point
    h_eva = CP.PropsSI('H', 'Q', 1, 'T', T_hs_bf - 5 + 273, 'NH3') * 1e-3
    p_eva = CP.PropsSI('P', 'Q', 1, 'T', T_hs_bf - 5 + 273, 'NH3') * 1e-5
    hs_eva2int_heatex.set_attr(x=1, p=p_eva)

    # condensation point
    h_cond = CP.PropsSI('H', 'Q', 0, 'T', T_cons_ff + 5 + 273, 'NH3') * 1e-3
    p_cond = CP.PropsSI('P', 'Q', 0, 'T', T_cons_ff + 5 + 273, 'NH3') * 1e-5
    cond2int_heatex.set_attr(p=p_cond)

    # internal heat exchanger to valve
    int_heatex2valve.set_attr(h=h_cond * 0.99, fluid={'water': 0, 'NH3': 1})

    # consumer cycle
    cond2cons_hs.set_attr(T=T_cons_ff, p=10, fluid={'water': 1, 'NH3': 0})
    cons_hs2cons_feed.set_attr(T=T_cons_bf)

    # heat source cycle
    hs_feed2hs_pump.set_attr(T=T_hs_ff, p=1, fluid={'water': 1, 'NH3': 0})
    hs_eva2hs_back.set_attr(T=T_hs_bf, p=1)

Some components have to be parameterized. For the heat source and heat sink
recirculation pump as well as the conedenser the isentropic efficiency is to
be set. Further we set the pressure ratios on hot and cold side for the
condenser, evaporator and internal heat exchanger. The consumer will have
pressure losses, too.

.. code-block:: python

    # parametrization components
    # isentropic efficiency
    cons_pump.set_attr(eta_s=0.8)
    heatsource_pump.set_attr(eta_s=0.8)
    compressor.set_attr(eta_s=0.85)

    # pressure ratios
    condenser.set_attr(pr1=0.99, pr2=0.99)
    heatsource_evaporator.set_attr(pr1=0.98, pr2=0.98)
    cons_heatsink.set_attr(pr=0.99)
    int_heatex.set_attr(pr1=0.99, pr2=0.99)

The most important parameter is the consumers heat demand setting as
“key parameter”. After that the network can be solved and a stable simulation
can be used for further simulations.

.. code-block:: python

    # key parameter
    cons_heatsink.set_attr(P.val=-1e6)

    # solve the network
    nw.solve('design')
    nw.print_results()

    # calculate and print COP
    cop = abs(
        cons_heatsink.P.val
        / (cons_pump.P.val + heatsource_pump.P.val + compressor.P.val)
        )
    print(f'COP = {cop:.4}')

After that, enthalpies and pressures can be set as "None" and the desired
values for the upper or lower terminal temperature differences, references or
other unstable values can be set for the actual simulation.

.. code-block:: python

    # parametrization for the actual simulation
    hs_eva2int_heatex.set_attr(p=None)
    heatsource_evaporator.set_attr(ttd_l=5)

    cond2int_heatex.set_attr(p=None)
    condenser.set_attr(ttd_u=5)

    int_heatex2valve.set_attr(h=None)
    int_heatex2comp.set_attr(T=Ref(hs_eva2int_heatex, 1, 10))

    # solve the actual network
    nw.solve('design')
    nw.print_results()

    # calculate and print the actual COP
    cop = abs(
        cons_heatsink.P.val
        / (cons_pump.P.val + heatsource_pump.P.val + compressor.P.val)
        )
    print(f'COP = {cop:.4}')
