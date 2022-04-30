Distric heating system
----------------------

The district heating system is a great example for the usage of flexible
user-defined subsystems. The example system and data are based on the district
heating system Hamburg Wilhelmsburg :cite:`Lorenzen2014`. The source code for
this example can be found
`here <https://github.com/oemof/oemof-examples/tree/master/oemof_examples/tespy/district_heating>`__.
Although the structure of the system (see the Figure below) does not seem very
complex, it has more than 120 components. But we can easily determine repeating
structures for the consumers and this is, where the subsystems come in place.

.. figure:: api/_images/dhs.svg
    :align: center

    Figure: Topology of the heating system.

The single consumers are connected to the main grid with a control valve at
the outlet and each fork is connected with a pipe to the next fork. Also, the
main grid may have a dead end (e.g. in the housing areas, see subsystem
closed) or is open to connect to another part of the grid (industrial area,
see subsystem open). Additionally, each branch of the main grid is connected to
the upstream part with the fork subsystem (Ki, see subsystem fork).

.. figure:: api/_images/dhs_closed.svg
    :align: center

    Figure: Generic topology of the dead end subsystem.

.. figure:: api/_images/dhs_open.svg
    :align: center

    Figure: Generic topology of the open subsystem.

.. figure:: api/_images/dhs_forks.svg
    :align: center

    Figure: Generic topology of the forks (variable number of branches).

After the system has been set up, we designed the pipes' insulation in a way,
that the feed flow the temperature gradient is at 1 K / 100 m and the back flow
gradient is at 0.5 K / 100 m. Having designed the system, heat losses at
different ambient temperatures can be calculated, as the heat transfer
coefficient for the pipes has been calculated in the design case. By this way,
it is for example possible to apply load profiles for the consumers as well as
a profile for the ambient temperature to investigate the network heat losses
over a specific period of time.
