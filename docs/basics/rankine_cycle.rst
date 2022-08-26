.. _tespy_basics_rankine_cycle_label:

Rankine Cycle
=============

.. figure:: /_static/images/basics/rankine_cycle.svg
    :align: center
    :alt: Topology of the rankine cycle
    :figclass: only-light

    Figure: Topology of the rankine cycle

.. figure:: /_static/images/basics/rankine_cycle_darkmode.svg
    :align: center
    :alt: Topology of the rankine cycle
    :figclass: only-dark

    Figure: Topology of the rankine cycle

Setting up the Cycle
^^^^^^^^^^^^^^^^^^^^
We will model the cycle including the cooling water of the condenser. For this
start with the :code:`Network` set up we already know.

.. code-block:: python

    from tespy.networks import Network

    # create a network object with R134a as fluid
    fluid_list = ['water']
    my_plant = Network(fluids=fluid_list)
    my_plant.set_attr(T_unit='C', p_unit='bar', h_unit='kJ / kg')

Following, we create the components and connect them. The :code:`Condenser` has
a hot side inlet and outlet as well as a cold side inlet and outlet. The hot
side is indicated by using the index 1 for the inlet and outlet :code:`in1` and
:code:`out1`, the cold side uses the index 2 (:code:`in2` and :code:`out2`).

Again, for the closed thermodynamic cycle we have to insert a cycle closer. The
cooling water inlet and the cooling water outlet of the condenser are directly
connected to a :code:`Source` and a :code:`Sink` respectively.

.. code-block:: python

    from tespy.components import (
        CycleCloser, Pump, Condenser, Turbine, HeatExchangerSimple, Source, Sink
    )

    cc = CycleCloser('cycle closer')
    sg = HeatExchangerSimple('steam generator')
    mc = Condenser('main condenser')
    tu = Turbine('steam turbine')
    fp = Pump('feed pump')

    cwso = Source('cooling water source')
    cwsi = Sink('cooling water sink')

    from tespy.connections import Connection

    c1 = Connection(cc, 'out1', tu, 'in1', label='1')
    c2 = Connection(tu, 'out1', mc, 'in1', label='2')
    c3 = Connection(mc, 'out1', fp, 'in1', label='3')
    c4 = Connection(fp, 'out1', sg, 'in1', label='4')
    c0 = Connection(sg, 'out1', cc, 'in1', label='0')

    my_plant.add_conns(c1, c2, c3, c4, c0)

    c11 = Connection(cwso, 'out1', mc, 'in2', label='11')
    c12 = Connection(mc, 'out2', cwsi, 'in1', label='12')

    my_plant.add_conns(c11, c12)

For the parameters, we predefine the pressure losses in the heat exchangers.
For the condenser, the hot side pressure losses are neglected :code:`pr1=1`,
for the cooling water side we assume pressure loss of 2 % :code:`pr2=0.98`. The
steam generator inflicts a pressure loss of 10 %.

The turbine and feed pump will have the isentropic efficiency specified. For
the connection parameters, the fluid has to be defined in both the main cycle
and the cooling water system. Furthermore, the live steam temperature, pressure
and mass flow are set. Lastly, we set the condensation pressure level and the
feed and return flow temperature of the cooling water as well as its feed
pressure.

.. code-block:: python

    mc.set_attr(pr1=1, pr2=0.98)
    sg.set_attr(pr=0.9)
    tu.set_attr(eta_s=0.9)
    fp.set_attr(eta_s=0.75)

    c11.set_attr(T=20, p=1.2, fluid={'water': 1})
    c12.set_attr(T=35)
    c1.set_attr(T=600, p=150, m=10, fluid={'water': 1})
    c2.set_attr(p=0.1)

    my_plant.solve(mode='design')
    my_plant.print_results()

After running the simulation, for example, we can observe the temperature
differences at the condenser. Instead of directly setting a pressure value for
condensation, we could also set the upper terminal temperature difference
:code:`ttd_u` instead. It is defined as the condensation temperature to cooling
water return flow temperature.

.. tip::

    You will find the documentation of each equation of the components in the
    respective seciton of the API documentation. For example, the condenser
    :py:class:`tespy.components.heat_exchangers.condenser.Condenser`.

.. code-block:: python

    mc.set_attr(pr1=1, pr2=0.98, ttd_u=4)
    c2.set_attr(p=None)

    my_plant.solve(mode='design')
    my_plant.print_results()

After rerunning, we will see that the condensation temperature and pressure
are both automatically calculated by the specified terminal temperature value.

Assess Electrical Power
^^^^^^^^^^^^^^^^^^^^^^^
To assess the electrical power output we want to consider the power generated
by the turbine as well as the power required to drive the feed pump. It is
possible to include both of the component's power values in a single electrical
:code:`Bus`. We can do this by importing the
:py:class:`tespy.connections.bus.Bus` class, creating an instance and adding
both components to the bus.

.. code-block:: python

    from tespy.connections import Bus

    powergen = Bus("electrical power output")

    powergen.add_comps(
        {"comp": tu, "char": 0.97, "base": "component"},
        {"comp": fp, "char": 0.97, "base": "bus"},
    )

    my_plant.add_busses(powergen)

    my_plant.solve(mode='design')
    my_plant.print_results()

.. note::

    The :code:`Bus` can take components which either produce or consume energy.
    Specifying :code:`'base': 'bus'` means, that the efficiency value is
    referenced to the electrical power

    .. math::

        \dot{W} = \dot{W}_\text{el} \cdot \eta

    while specifying :code:`'base': 'component'` (default) takes the component's
    power as base value.

    .. math::

        \dot{W}_\text{el} = \dot{W} \cdot \eta

The results for the bus are printed separately. Observe, that the steam
turbine's electrical power production (:code:`bus value`) is lower than the
:code:`component value`, while it is inverted for the feed pump.

You can also set the total desired power production of the system, for example
replacing the mass flow specification at connection 1:

.. code-block:: python

    powergen.set_attr(P=-10e6)
    c1.set_attr(m=None)

    my_plant.solve(mode='design')
    my_plant.print_results()

Analyze Efficiency and Powergeneration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In this section, we will analyze the power production and the efficiency
of the cycle, given constant steam mass flow and with varying values for the

- live steam pressure,
- live steam temperature and
- cooling water temperature level.

To do that, we are using a very similar setup as has been used in the
:ref:`heat pump tutorial <tespy_basics_heat_pump_label>`.

.. dropdown:: Click to expand to code section

    .. code-block:: python

        import matplotlib.pyplot as plt
        import numpy as np



        fig, ax = plt.subplots(2, 3, sharey=True, sharex=True, figsize=(16, 8))

        [a.grid() for a in ax]

        ax[0].set_ylabel('Efficiency of the rankine cycle in %')
        ax[3].set_ylabel('Power of the rankine cycle in %')
        plt.tight_layout()
        fig.savefig('eta_power_parametric.svg')

.. figure:: /_static/images/basics/eta_power_parametric.svg
    :align: center
    :alt: Parametric analysis of the efficiency and power output

    Figure: Parametric analysis of the efficiency and power output

Partload Simulation
^^^^^^^^^^^^^^^^^^^
In the partload simulation part, we are starting with a specific design of the
plant and calculate the partload perfomance with some assumptions on the
component's individual behavior. The table below summarizes the assumptions,
which we will keep as simple as possible in this moment. For more insights
have a look at the step by step
:ref:`heat pump tutorial <tespy_tutorial_heat_pump_label>` or at the
:ref:`Network documentation <tespy_modules_networks_label>`.

+-----------+---------------------------+-------------------------------+
| Component | Assumptions               | Settings                      |
+===========+===========================+===============================+
| Turbine   | cone law applies          | unset inlet pressure and      |
|           |                           | apply cone law                |
+-----------+---------------------------+-------------------------------+
| Condenser | constant heat transfer    | unset terminal temperature    |
|           | coefficient               | difference and set heat       |
|           |                           | transfer coefficient          |
+-----------+---------------------------+-------------------------------+
| Cooling   | constant volumetric flow  | unset return temperature      |
| water     |                           | value and set volumetric flow |
+-----------+---------------------------+-------------------------------+

With these specifications, the following physics are applied to the model:

- Due to the constant volumetric flow of water, the temperature of the cooling
  water returning from the condenser will react to the total heat transferred
  in the condensation: Increased heat transfer means incresing temperature,
  decreased heat transfer means decreased temperature.
- The constant heat transfer coefficient of the condenser will calculate the
  condensation temperature (and therefore pressure) based on the temperature
  regime in the cooling water side:

  - Increase in temperature for the cooling water leads to increased
    condensation temperature (at constant heat transfer).
  - Increase in heat transfer means increase in necessary temperature
    difference at the condenser (at constant cooling water inlet temperature).

- The cone law is a mathematical model to predict the pressure at the turbine's
  inlet based on the deviation from the design conditions. Generally,
  increased mass flow leads to higher inlet pressure (at constant inlet
  temperature and constant outlet pressure). However, this equation is more
  complex, since a lot more parameters are involved compared to the other
  equations applied.

.. code-block:: python

    some code
