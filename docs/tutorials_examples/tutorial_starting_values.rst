Stable starting values for subcritical heat pumps
-------------------------------------------------

.. contents::
    :depth: 1
    :local:
    :backlinks: top

Applying numerical algorithms and methods, the starting value of a variable
is the value used for the first iteration. With more complex TESPy models
it can happen that the simulation does not converge easily due to a combination
of "bad" starting values. The solver is especially vulnerable if the specified
parameters trigger complex equations with respect to the primary variables.

The primary variables of TESPy are mass flow, pressure, enthalpy and fluid
composition. If such a value is directly specified by the user, the solver has
a solution for this value before starting the first iteration. Therefore,
specifying a set of parameters largely including primary variables will improve
the convergence significantly. Based on the converged solution of a initial
simulation, it is then possible to adjust the parameters, for example, unsetting
pressure values and specifying efficiencies instead.

Here we provide a short tutorial for you to better understand, how this process
could look like at the example of a subcritical heat pump with different working
fluids.

.. note::

    If the heat pump operates in trans- or supercritical range, some
    modifications have to be made on this setup. We plan to include respective
    examples here in the future.

Topology of the heat pump
^^^^^^^^^^^^^^^^^^^^^^^^^

Following the first tutorial a slightly different topology for a heat pump with
internal heat exchangers is considered instead of dumping the heat to the
ambient. You can see the plant topology in the figure below.

.. figure:: api/_images/tutorial_sv_heat_pump_intheatex.svg
    :align: center

    Figure: Topology of heat pump with internal heat exchanger

The system consists of a consumer system, a valve, an evaporator system, a
compressor and additionally an internal heat exchanger. In order to simulate
this heat pump, the TESPy model has to be built up. First, the network has to
be initialized and the refrigerants used have to be specified. This example
shows how to make the heat pump model work with a variety of working fluids with
water on both the heat source and heat sink side of the system.

Running into errors
^^^^^^^^^^^^^^^^^^^

As always, we start by importing the necessary TESPy classes.

.. code-block:: python

    from tespy.networks import Network

    from tespy.components import (
        Condenser, Compressor, CycleCloser,  HeatExchanger,
        HeatExchangerSimple, Pump, Sink, Source, Valve
        )

    from tespy.connections import Connection, Ref, Bus

Then, we can build the network by defining components and connections. The
working fluid will be set with the variable `wf`, `"NH3"` is used in the first
setup. This way, we will be able to change the working fluid in a flexible way.

.. code-block:: python

    wf = 'NH3'

    # network
    nw = Network(
        fluids=['water', wf],
        T_unit='C', p_unit='bar', h_unit='kJ / kg', m_unit='kg / s'
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

    # connections
    # main cycle
    cc2hs_eva = Connection(cycle_closer, 'out1', heatsource_evaporator, 'in2', label='0')
    hs_eva2int_heatex = Connection(heatsource_evaporator, 'out2', int_heatex, 'in2', label='1')
    int_heatex2comp = Connection(int_heatex, 'out2', compressor, 'in1', label='2')
    comp2cond = Connection(compressor, 'out1', condenser, 'in1', label='3')
    cond2int_heatex = Connection(condenser, 'out1', int_heatex, 'in1', label='4')
    int_heatex2valve = Connection(int_heatex, 'out1', valve, 'in1', label='5')
    valve2cc = Connection(valve, 'out1', cycle_closer, 'in1', label='6')

    nw.add_conns(
        cc2hs_eva, hs_eva2int_heatex, int_heatex2comp, comp2cond, cond2int_heatex,
        int_heatex2valve, valve2cc
        )

    # heat source
    hs_feed2hs_pump = Connection(heatsource_feedflow, 'out1', heatsource_pump, 'in1', label='11')
    hs_pump2hs_eva = Connection(heatsource_pump, 'out1', heatsource_evaporator, 'in1', label='12')
    hs_eva2hs_back = Connection(heatsource_evaporator, 'out1', heatsource_backflow, 'in1', label='13')

    nw.add_conns(hs_feed2hs_pump, hs_pump2hs_eva, hs_eva2hs_back)

    # heat sink
    cons_back2cons_pump = Connection(cons_cycle_closer, 'out1', cons_pump, 'in1', label='21')
    cons_pump2cond = Connection(cons_pump, 'out1', condenser, 'in2', label='22')
    cond2cons_hs = Connection(condenser, 'out2', cons_heatsink, 'in1', label='23')
    cons_hs2cons_feed = Connection(cons_heatsink, 'out1', cons_cycle_closer, 'in1', label='24')

    nw.add_conns(cons_back2cons_pump, cons_pump2cond, cond2cons_hs, cons_hs2cons_feed)

After setting up the topology, the system's parameters should be set in the
following way:

- Heat sink temperature levels (`T` at 23 and 24)
- Heat source temperature levels (`T` at 11 and 13)
- Degree of overheating after the internal heat exchanger (`Td_bp` at 2)
- Pinch point temperature difference at the evaporator (`ttd_l`) to derive
  evaporation pressure
- Temperature difference at the condenser (`ttd_u`) to derive condensation
  pressure
- Saturated gaseous state of the working fluid (`x=1`) after leaving the
  evaporator
- Efficiencies of pumps and the compressor (`eta_s`)
- Pressure losses in all heat exchangers (`pr1`, `pr2`, `pr`)
- Consumer heat demand (`Q`)

.. code-block:: python

    # parametrization connections
    # set feedflow and backflow temperature of heat source and consumer
    T_hs_bf = 10
    T_hs_ff = 15
    T_cons_bf = 50
    T_cons_ff = 90

    # consumer cycle
    cond2cons_hs.set_attr(T=T_cons_ff, p=10, fluid={'water': 1, wf: 0})
    cons_hs2cons_feed.set_attr(T=T_cons_bf)

    # heat source cycle
    hs_feed2hs_pump.set_attr(T=T_hs_ff, p=1, fluid={'water': 1, wf: 0})
    hs_eva2hs_back.set_attr(T=T_hs_bf, p=1)

    # evaporation to fully saturated gas
    hs_eva2int_heatex.set_attr(x=1, fluid={'water': 0, wf: 1})
    # degree of overheating after internal heat exchanger (evaporation side)
    int_heatex2comp.set_attr(Td_bp=10)

    # parametrization components
    # isentropic efficiency
    cons_pump.set_attr(eta_s=0.8)
    heatsource_pump.set_attr(eta_s=0.8)
    compressor.set_attr(eta_s=0.85)

    # pressure ratios
    condenser.set_attr(pr1=0.98, pr2=0.98)
    heatsource_evaporator.set_attr(pr1=0.98, pr2=0.98)
    cons_heatsink.set_attr(pr=0.99)
    int_heatex.set_attr(pr1=0.98, pr2=0.98)

    # temperature differences
    heatsource_evaporator.set_attr(ttd_l=5)
    condenser.set_attr(ttd_u=5)

    # consumer heat demand
    cons_heatsink.set_attr(Q=-1e6)

    nw.solve('design')

The system should be well defined with the parameter settings, however no
solution can be found. We might run in some error, like

.. error::

    .. code-block:: bash

        ERROR:root:Singularity in jacobian matrix, calculation aborted! Make
        sure your network does not have any linear dependencies in the
        parametrisation. Other reasons might be

        -> given temperature with given pressure in two phase region, try
        setting enthalpy instead or provide accurate starting value for
        pressure.

        -> given logarithmic temperature differences or kA-values for heat
        exchangers,

        -> support better starting values.

        -> bad starting value for fuel mass flow of combustion chamber, provide
        small (near to zero, but not zero) starting value.

or simply not making progress in the convergence

.. error::

    .. code-block:: bash

        WARNING:root:The solver does not seem to make any progress, aborting
        calculation. Residual value is 7.43e+05. This frequently happens, if
        the solver pushes the fluid properties out of their feasible range.

Fixing the errors
^^^^^^^^^^^^^^^^^

To generate good starting values for the simulation, it is recommended to set
pressure and enthalpy values instead of temperature differences. In this
example, fixed points can be identified with the help of the logph diagram
which you can see in the figure below.

.. figure:: api/_images/tutorial_sv_logph.svg
    :align: center

    Figure: Logph diagram of ammonia

A rough estimation of the evaporation and condensation pressure can be obtained
and will be used to replace the temperature differences at the evaporator and
the condenser for the starting value generator. After condensation, the working
fluid is in saturated liquid state. We can retrieve the condensation pressure
corresponding to a temperature slightly below the heat sink temperature by using
the CoolProp `PropsSI` interface with the respective inputs. The same step can
be carried out on the heat source side. For the internal heat exchanger, an
enthalpy value is specified instead of the temperature difference to the boiling
point as well. It is important to note that the PropertySI function (PropsSI) is
used with SI units, which differ from the units defined in the network.

The temperature difference values are unset and pressure and enthalpy values are
set instead.

.. code-block:: python

    import CoolProp.CoolProp as CP

    # evaporation point
    p_eva = CP.PropsSI('P', 'Q', 1, 'T', T_hs_bf - 5 + 273.15, wf) * 1e-5
    hs_eva2int_heatex.set_attr(p=p_eva)
    heatsource_evaporator.set_attr(ttd_l=None)

    # condensation point
    p_cond = CP.PropsSI('P', 'Q', 0, 'T', T_cons_ff + 5 + 273.15, wf) * 1e-5
    cond2int_heatex.set_attr(p=p_cond)
    condenser.set_attr(ttd_u=None)

    # internal heat exchanger to compressor enthalpy
    h_evap = CP.PropsSI('H', 'Q', 1, 'T', T_hs_bf - 5 + 273.15, wf) * 1e-3
    int_heatex2comp.set_attr(Td_bp=None, h=h_evap * 1.01)

    # solve the network again
    nw.solve('design')


The model was solved successfully and has stored the starting values for any
follow-up. Therefore, we can undo our recent changes and restart the
simulation. For example, the COP is then calculated.

.. code-block:: python

    # evaporation point
    hs_eva2int_heatex.set_attr(p=None)
    heatsource_evaporator.set_attr(ttd_l=5)

    # condensation point
    cond2int_heatex.set_attr(p=None)
    condenser.set_attr(ttd_u=5)

    # internal heat exchanger superheating
    int_heatex2comp.set_attr(Td_bp=5, h=None)

    # solve the network again
    nw.solve('design')

    # calculate the COP
    cop = abs(
        cons_heatsink.Q.val
        / (cons_pump.P.val + heatsource_pump.P.val + compressor.P.val)
    )

Expand fix to any working fluids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally, using this strategy, it is possible to build a generic function,
building a network, that works with a variety of working fluids.

.. code-block:: python

    import matplotlib.pyplot as plt
    import pandas as pd

    from tespy.networks import Network
    from tespy.components import (
        Condenser, Compressor, CycleCloser,  HeatExchanger,
        HeatExchangerSimple, Pump, Sink, Source, Valve
        )
    from tespy.connections import Connection, Ref, Bus
    import CoolProp.CoolProp as CP


    def generate_starting_values(wf):

        # network
        nw = Network(
            fluids=['water', wf],
            T_unit='C', p_unit='bar', h_unit='kJ / kg', m_unit='kg / s',
            iterinfo=False
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

        # connections
        # main cycle
        cc2hs_eva = Connection(cycle_closer, 'out1', heatsource_evaporator, 'in2', label='0')
        hs_eva2int_heatex = Connection(heatsource_evaporator, 'out2', int_heatex, 'in2', label='1')
        int_heatex2comp = Connection(int_heatex, 'out2', compressor, 'in1', label='2')
        comp2cond = Connection(compressor, 'out1', condenser, 'in1', label='3')
        cond2int_heatex = Connection(condenser, 'out1', int_heatex, 'in1', label='4')
        int_heatex2valve = Connection(int_heatex, 'out1', valve, 'in1', label='5')
        valve2cc = Connection(valve, 'out1', cycle_closer, 'in1', label='6')

        nw.add_conns(
            cc2hs_eva, hs_eva2int_heatex, int_heatex2comp, comp2cond, cond2int_heatex,
            int_heatex2valve, valve2cc
            )

        # heat source
        hs_feed2hs_pump = Connection(heatsource_feedflow, 'out1', heatsource_pump, 'in1', label='11')
        hs_pump2hs_eva = Connection(heatsource_pump, 'out1', heatsource_evaporator, 'in1', label='12')
        hs_eva2hs_back = Connection(heatsource_evaporator, 'out1', heatsource_backflow, 'in1', label='13')

        nw.add_conns(hs_feed2hs_pump, hs_pump2hs_eva, hs_eva2hs_back)

        # heat sink
        cons_back2cons_pump = Connection(cons_cycle_closer, 'out1', cons_pump, 'in1', label='20')
        cons_pump2cond = Connection(cons_pump, 'out1', condenser, 'in2', label='21')
        cond2cons_hs = Connection(condenser, 'out2', cons_heatsink, 'in1', label='22')
        cons_hs2cons_feed = Connection(cons_heatsink, 'out1', cons_cycle_closer, 'in1', label='23')

        nw.add_conns(cons_back2cons_pump, cons_pump2cond, cond2cons_hs, cons_hs2cons_feed)

        # set feedflow and backflow temperature of heat source and consumer
        T_hs_bf = 10
        T_hs_ff = 15
        T_cons_bf = 50
        T_cons_ff = 90

        # consumer cycle
        cond2cons_hs.set_attr(T=T_cons_ff, p=10, fluid={'water': 1, wf: 0})
        cons_hs2cons_feed.set_attr(T=T_cons_bf)

        # heat source cycle
        hs_feed2hs_pump.set_attr(T=T_hs_ff, p=1, fluid={'water': 1, wf: 0})
        hs_eva2hs_back.set_attr(T=T_hs_bf, p=1)

        # evaporation to fully saturated gas
        hs_eva2int_heatex.set_attr(x=1, fluid={'water': 0, wf: 1})

        # parametrization components
        # isentropic efficiency
        cons_pump.set_attr(eta_s=0.8)
        heatsource_pump.set_attr(eta_s=0.8)
        compressor.set_attr(eta_s=0.85)

        # pressure ratios
        condenser.set_attr(pr1=0.98, pr2=0.98)
        heatsource_evaporator.set_attr(pr1=0.98, pr2=0.98)
        cons_heatsink.set_attr(pr=0.99)
        int_heatex.set_attr(pr1=0.98, pr2=0.98)

        # evaporation point
        p_eva = CP.PropsSI('P', 'Q', 1, 'T', T_hs_bf - 5 + 273.15, wf) * 1e-5
        hs_eva2int_heatex.set_attr(p=p_eva)

        # condensation point
        p_cond = CP.PropsSI('P', 'Q', 0, 'T', T_cons_ff + 5 + 273.15, wf) * 1e-5
        cond2int_heatex.set_attr(p=p_cond)

        # internal heat exchanger to compressor enthalpy
        h_evap = CP.PropsSI('H', 'Q', 1, 'T', T_hs_bf - 5 + 273.15, wf) * 1e-3
        int_heatex2comp.set_attr(h=h_evap * 1.01)

        # consumer heat demand
        cons_heatsink.set_attr(Q=-1e6)

        power_bus = Bus('Total power input')
        heat_bus = Bus('Total heat production')
        power_bus.add_comps(
            {'comp': compressor, 'base': 'bus'},
            {'comp': cons_pump, 'base': 'bus'},
            {'comp': heatsource_pump, 'base': 'bus'},
        )
        heat_bus.add_comps({'comp': cons_heatsink})

        nw.add_busses(power_bus, heat_bus)

        nw.solve('design')

            # evaporation point
        hs_eva2int_heatex.set_attr(p=None)
        heatsource_evaporator.set_attr(ttd_l=5)

        # condensation point
        cond2int_heatex.set_attr(p=None)
        condenser.set_attr(ttd_u=5)

        # internal heat exchanger superheating
        int_heatex2comp.set_attr(Td_bp=5, h=None)

        # solve the network again
        nw.solve('design')

        return nw


    cop = pd.DataFrame(columns=["COP"])

    for wf in ['NH3', 'R22', 'R134a', 'R152a', 'R290', 'R718']:
        nw = generate_starting_values(wf)

        power = nw.busses['Total power input'].P.val
        heat = abs(nw.busses['Total heat production'].P.val)
        cop.loc[wf] = heat / power


    fig, ax = plt.subplots(1)

    cop.plot.bar(ax=ax, legend=False)

    ax.set_axisbelow(True)
    ax.yaxis.grid(linestyle='dashed')
    ax.set_xlabel('Name of working fluid')
    ax.set_ylabel('Coefficicent of performance')
    ax.set_title('Coefficicent of performance for different working fluids')
    plt.tight_layout()

    fig.savefig('tutorial_sv_COP_by_wf.svg')


.. figure:: api/_images/tutorial_sv_COP_by_wf.svg
    :align: center

    Figure: Topology of heat pump with internal heat exchanger

Of course, there are different strategies, which include building the plant
step by step and successively adding more and more components.
