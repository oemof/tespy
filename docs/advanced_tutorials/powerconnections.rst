.. _tutorial_powerconnection_label:

PowerConnection examples
^^^^^^^^^^^^^^^^^^^^^^^^
This page shows a couple of more elaborate examples involving
:code:`PowerConnection` instances.

Single shaft gas turbine
++++++++++++++++++++++++

We will implement an open gas turbine system using air as working fluid and a
heater. You can also model gas turbines with combustion, for this example, the
focus is on modeling the single shaft gas turbine system.

First, we import the necessary components and set up the material flow
system connecting the compressor to the heater and to the turbine.

.. code-block:: python

    >>> from tespy.connections import Connection, PowerConnection
    >>> from tespy.components import (
    ...     Turbine, Source, Sink, Compressor, SimpleHeatExchanger, PowerBus,
    ...     Generator, PowerSink
    ... )
    >>> from tespy.networks import Network
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(temperature="degC", pressure="bar")
    >>> so = Source("source")
    >>> heater = SimpleHeatExchanger("heater")
    >>> compressor = Compressor("compressor")
    >>> turbine = Turbine("turbine")
    >>> si = Sink("sink")
    >>> c1 = Connection(so, "out1", compressor, "in1", label="c1")
    >>> c2 = Connection(compressor, "out1", heater, "in1", label="c2")
    >>> c3 = Connection(heater, "out1", turbine, "in1", label="c3")
    >>> c4 = Connection(turbine, "out1", si, "in1", label="c4")

Next, we can set up the energy flows. Since the turbine and the compressor
are rotating on the same shaft, the turbine powers the compressor and the
generator at the same time. For this, we can use a PowerBus to represent
the shaft, which gets powered by the turbine. The turbine's power connector
is connected to the 'power_in1' connector of the shaft. The shaft
connects to the compressor's connector 'power' and to the generator's
connector 'power_in'. The generator then is connected at its outlet
'power_out' to the grid representation at the connector 'power'.

.. code-block:: python

    >>> shaft = PowerBus("shaft", num_in=1, num_out=2)
    >>> generator = Generator("generator")
    >>> grid = PowerSink("grid")
    >>> e1 = PowerConnection(turbine, "power", shaft, "power_in1", label="e1")
    >>> e2 = PowerConnection(shaft, "power_out1", compressor, "power", label="e2")
    >>> e3 = PowerConnection(shaft, "power_out2", generator, "power_in", label="e3")
    >>> e4 = PowerConnection(generator, "power_out", grid, "power", label="e4")
    >>> nw.add_conns(c1, c2, c3, c4, e1, e2, e3, e4)

We can parametrize the system, in this example, we fix the ambient air
temperature, pressure and mass flow, the turbine inlet temperature and the
turbine outlet pressure (to be equal to the ambient pressure).

.. code-block:: python

    >>> c1.set_attr(fluid={"air": 1}, m=1, p=1, T=25)
    >>> c3.set_attr(T=1000)
    >>> c4.set_attr(p=1)

In the components the turbine's and the compressor's efficiency are set as
well as the compressor's pressure ratio and the heater's pressure drop.

.. code-block:: python

    >>> turbine.set_attr(eta_s=0.9)
    >>> compressor.set_attr(eta_s=0.9, pr=15)
    >>> heater.set_attr(dp=0)

With the four power connections we have four additional variables in our
system. The compressor, the turbine and the shaft all deliver one equation
(their energy balance equation), meaning, one parameter is missing to fully
set up our problem. This could be the generator efficiency. With that, we
can solve the system, and check what amount of electricity is generated
through the generator.

.. code-block:: python

    >>> generator.set_attr(eta=0.95)
    >>> nw.solve("design")
    >>> round(e4.E.val_SI / 1e3, 1)
    246.9

Alternatively, we could also fix the electricity output to a specific
target value and unset the air mass flow. This will calculate the required
air mass flow to generate the desired amount of electricity.

.. code-block:: python

    >>> e4.set_attr(E=3e5)
    >>> c1.set_attr(m=None)
    >>> nw.solve("design")
    >>> round(c1.m.val, 3)
    1.215

Single shaft feed water pump powered by a turbine
+++++++++++++++++++++++++++++++++++++++++++++++++

Create a pump that is powered by a turbine. The turbine's :code:`turbine_fwp`
power output must therefore be equal to the pump's :code:`fwp` power
consumption.

.. code-block:: python

    >>> from tespy.networks import Network
    >>> from tespy.components import Pump, Turbine, Source, Sink
    >>> from tespy.connections import Connection, PowerConnection

    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(temperature="degC", pressure="bar")
    >>> cond = Source("condensate")
    >>> fwp = Pump("feed water pump")
    >>> feedwater = Sink("feedwater")
    >>> c1 = Connection(cond, "out1", fwp, "in1")
    >>> c2 = Connection(fwp, "out1", feedwater, "in1")
    >>> ls = Source("live steam")
    >>> turbine_fwp = Turbine("turbine fwp")
    >>> ws = Sink("Waste steam")
    >>> c11 = Connection(ls, "out1", turbine_fwp, "in1")
    >>> c12 = Connection(turbine_fwp, "out1", ws, "in1")
    >>> e1 = PowerConnection(turbine_fwp, "power", fwp, "power")
    >>> nw.add_conns(c1, c2, c11, c12, e1)

We can set up the system in a way, that calculates the required mass flow
of steam through the turbine to power the feed water pump and find the
power flow by accessing the respective attribute of the power connection.

.. code-block:: python

    >>> c1.set_attr(fluid={"water": 1}, p=0.5, x=0, m=10)
    >>> c2.set_attr(p=50)
    >>> fwp.set_attr(eta_s=0.75)
    >>> c11.set_attr(fluid={"water": 1}, p=40, T=500)
    >>> c12.set_attr(p=0.55)
    >>> turbine_fwp.set_attr(eta_s=0.9)
    >>> nw.solve("design")
    >>> nw.assert_convergence()
    >>> round(e1.E.val_SI / 1e3)
    68

Logic to force same power of two compressors
++++++++++++++++++++++++++++++++++++++++++++

In this example we combine the PowerConnection with a UserDefinedEquation.
Two air compressors should run in series and at identical power. For this the
intermediate pressure is variable.

    >>> from tespy.components import Source, Sink, Compressor, PowerSource
    >>> from tespy.connections import Connection, PowerConnection
    >>> from tespy.networks import Network
    >>> from tespy.tools import UserDefinedEquation
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(temperature="degC", pressure="bar")
    >>> so = Source("air source")
    >>> compressor1 = Compressor("compressor 1")
    >>> compressor2 = Compressor("compressor 2")
    >>> si = Sink("compressed air")
    >>> grid1 = PowerSource("grid compressor 1")
    >>> grid2 = PowerSource("grid compressor 2")
    >>> c1 = Connection(so, "out1", compressor1, "in1", label="c1")
    >>> c2 = Connection(compressor1, "out1", compressor2, "in1", label="c2")
    >>> c3 = Connection(compressor2, "out1", si, "in1", label="c3")
    >>> e1 = PowerConnection(grid1, "power", compressor1, "power")
    >>> e2 = PowerConnection(grid2, "power", compressor2, "power")
    >>> nw.add_conns(c1, c2, c3, e1, e2)
    >>> c1.set_attr(fluid={"air": 1}, m=1, p=1, T=25)
    >>> c3.set_attr(p=5)
    >>> compressor1.set_attr(eta_s=0.85)
    >>> compressor2.set_attr(eta_s=0.85)
    >>> def same_power_ude(ude):
    ...     e1, e2 = ude.conns
    ...     return e1.E.val_SI - e2.E.val_SI
    >>> def same_power_dependents(ude):
    ...     e1, e2 = ude.conns
    ...     return [c.E for c in ude.conns]
    >>> ude = UserDefinedEquation(
    ...     "power equality ude",
    ...     func=same_power_ude,
    ...     dependents=same_power_dependents,
    ...     conns=[e1, e2]
    ... )
    >>> nw.add_ude(ude)
    >>> nw.solve("design")
    >>> nw.assert_convergence()
    >>> round(e1.E.val / 1e3) == round(e2.E.val / 1e3)
    True
    >>> round(e1.E.val / 1e3)
    105

Including part load model for motor efficiency
++++++++++++++++++++++++++++++++++++++++++++++

This example how a partload efficiency curve can be applied to a motor. For
this, let's assume the motor powers a refrigeration compressor. We can set up
the model by connecting the compressor to the refrigerant flows as usual and
add the :code:`PowerConnection` to the electricity grid via a :code:`Motor`
instance.

.. code-block:: python

    >>> from tespy.components import Source, Sink, Compressor, PowerSource, Motor
    >>> from tespy.connections import Connection, PowerConnection
    >>> from tespy.networks import Network
    >>> from tespy.tools import CharLine

    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(temperature="degC", pressure="bar")
    >>> so = Source("evaporated refrigerant")
    >>> compressor = Compressor("compressor")
    >>> si = Sink("compressed refrigerant")
    >>> grid = PowerSource("grid")
    >>> motor = Motor("motor")
    >>> c1 = Connection(so, "out1", compressor, "in1", label="c1")
    >>> c2 = Connection(compressor, "out1", si, "in1", label="c2")
    >>> e1 = PowerConnection(grid, "power", motor, "power_in", label="e1")
    >>> e2 = PowerConnection(motor, "power_out", compressor, "power", label="e2")
    >>> nw.add_conns(c1, c2, e1, e2)

The design efficiency is 0.98, the compressor's design efficiency is 0.85. On
top we fix the inlet state and mass flow as well as the compressor's pressure
ratio. For the characteristics of the motor's efficiency we can pass data to a
:code:`CharLine` instance, which is set to be used for the :code:`eta_char`
method in the model of the motor.

.. code-block:: python

    >>> c1.set_attr(fluid={"R290": 1}, m=1, td_dew=10, T=10)
    >>> compressor.set_attr(pr=3, eta_s=0.85, design=["eta_s"], offdesign=["eta_s_char"])
    >>> motor.set_attr(eta_char=CharLine(x=[0.5, 0.75, 1, 1.25], y=[0.9, 0.975, 1, 0.975]))
    >>> motor.set_attr(eta=0.98, design=["eta"], offdesign=["eta_char"])
    >>> nw.solve("design")
    >>> nw.save("design.json")
    >>> nw.assert_convergence()

After performing the design simulation we can change the fluid mass flow and
observe the change in efficiency of the motor:

.. code-block:: python

    >>> c1.set_attr(m=0.8)
    >>> nw.solve("offdesign", design_path="design.json", init_path="design.json")
    >>> nw.assert_convergence()
    >>> round(motor.eta.val, 3)
    0.966

.. note::

    As mentioned in the component section: It is also possible to import your
    custom characteristics from the :code:`HOME/.tespy/data` folder. Read more
    about this :ref:`here <modules_characteristics_label>`.
