Thermal Power Plant Efficiency Optimization
-------------------------------------------

.. contents::
    :depth: 1
    :local:
    :backlinks: top

Task
^^^^

Designing a power plant meets multiple different tasks, such as finding the
optimal fresh steam temperature and pressure to reduce exhaust steam water
content, or the optimization of extraction pressures to maximize cycle
efficiency and many more.

In case of a rather simple power plant topologies the task of finding optimized
values for e.g. extraction pressures is still manageable without any
optimization tool. As the topology becomes more complex and boundary
conditions come into play the usage of additional tools is recommended. The
following tutorial is intended to show the usage of PyGMO in combination with
TESPy to **maximize the cycle efficiency of a power plant with two**
**extractions.**

The source code can be found at the `tespy_examples repository
<https://github.com/oemof/oemof-examples/tree/master/oemof_examples/tespy/efficiency_optimization>`_.

.. figure:: api/_images/power_plant_two_extractions.svg
    :align: center

    Figure: Topology of the power plant.

What is PyGMO?
^^^^^^^^^^^^^^

PyGMO (Python Parallel Global Multiobjective Optimizer, :cite:`Biscani2020`) is
a library that provides a large number of evolutionary optimization algorithms.
PyGMO can be used to solve constrained, unconstrained, single objective and
multi objective problems.

Evolutionary Algorithms
+++++++++++++++++++++++

Evolutionary Algorithms (EA) are optimization algorithms inspired by biological
evolution. In a given population the algorithm uses the so called fitness
function to determine the quality of the solutions to each individual (set of
decision variables) problem. The best possible solution of the population is
called champion. Via mutation, recombination and selection your population
evolves to find better solutions.

EA will never find an exact solution to your problem. They can only give an
approximation for the real optimum.

Install PyGMO
+++++++++++++


Conda
#####

With the `conda
<https://conda.io/en/latest/>`_ package manager PyGMO is 
available for Linux, OSX and Windows thanks to the infrastructure of `conda-forge
<https://conda-forge.org/>`_:

.. code-block:: bash

    conda install -c conda-forge pygmo

pip
###

On Linux you also have the option to use the 
`pip
<https://pip.pypa.io/en/stable/>`_ package installer:

.. code-block:: bash

    conda install -c conda-forge pygmo
    
Windows user can perform an installation from source as an alternative to anaconda. 
For further information on this process we recommend the `PyGMO installation
<https://esa.github.io/pygmo2/install.html>`_

Creating your TESPy-Model
^^^^^^^^^^^^^^^^^^^^^^^^^

It is necessary to use object oriented programming in PyGMO. Therefore we create
a class :code:`PowerPlant` which contains our TESPy-Model and a function to
return the cycle efficiency.

.. code-block:: python

    from tespy.networks import network
    from tespy.components import (
        turbine, splitter, merge, condenser, pump, sink, source,
        heat_exchanger_simple, desuperheater, cycle_closer
    )
    from tespy.connections import connection, bus
    from tespy.tools import logger
    import logging

    import numpy as np


    logger.define_logging(screen_level=logging.ERROR)


    class PowerPlant():

        def __init__(self):
            self.nw = network(
                fluids=['BICUBIC::water'],
                p_unit='bar', T_unit='C', h_unit='kJ / kg',
                iterinfo=False)
            # components
            # main cycle
            eco = heat_exchanger_simple('economizer')
            eva = heat_exchanger_simple('evaporator')
            sup = heat_exchanger_simple('superheater')
            cc = cycle_closer('cycle closer')
            hpt = turbine('high pressure turbine')
            sp1 = splitter('splitter 1', num_out=2)
            mpt = turbine('mid pressure turbine')
            sp2 = splitter('splitter 2', num_out=2)
            lpt = turbine('low pressure turbine')
            con = condenser('condenser')
            pu1 = pump('feed water pump')
            fwh1 = condenser('feed water preheater 1')
            fwh2 = condenser('feed water preheater 2')
            dsh = desuperheater('desuperheater')
            me2 = merge('merge2', num_in=2)
            pu2 = pump('feed water pump 2')
            pu3 = pump('feed water pump 3')
            me = merge('merge', num_in=2)

            # cooling water
            cwi = source('cooling water source')
            cwo = sink('cooling water sink')

            # connections
            # main cycle
            cc_hpt = connection(cc, 'out1', hpt, 'in1', label='feed steam')
            hpt_sp1 = connection(hpt, 'out1', sp1, 'in1', label='extraction1')
            sp1_mpt = connection(sp1, 'out1', mpt, 'in1', state='g')
            mpt_sp2 = connection(mpt, 'out1', sp2, 'in1', label='extraction2')
            sp2_lpt = connection(sp2, 'out1', lpt, 'in1')
            lpt_con = connection(lpt, 'out1', con, 'in1')
            con_pu1 = connection(con, 'out1', pu1, 'in1')
            pu1_fwh1 = connection(pu1, 'out1', fwh1, 'in2')
            fwh1_me = connection(fwh1, 'out2', me, 'in1', state='l')
            me_fwh2 = connection(me, 'out1', fwh2, 'in2', state='l')
            fwh2_dsh = connection(fwh2, 'out2', dsh, 'in2', state='l')
            dsh_me2 = connection(dsh, 'out2', me2, 'in1')
            me2_eco = connection(me2, 'out1', eco, 'in1', state='l')
            eco_eva = connection(eco, 'out1', eva, 'in1')
            eva_sup = connection(eva, 'out1', sup, 'in1')
            sup_cc = connection(sup, 'out1', cc, 'in1')

            self.nw.add_conns(cc_hpt, hpt_sp1, sp1_mpt, mpt_sp2, sp2_lpt,
                              lpt_con, con_pu1, pu1_fwh1, fwh1_me, me_fwh2,
                              fwh2_dsh, dsh_me2, me2_eco, eco_eva, eva_sup, sup_cc)

            # cooling water
            cwi_con = connection(cwi, 'out1', con, 'in2')
            con_cwo = connection(con, 'out2', cwo, 'in1')

            self.nw.add_conns(cwi_con, con_cwo)

            # preheating
            sp1_dsh = connection(sp1, 'out2', dsh, 'in1')
            dsh_fwh2 = connection(dsh, 'out1', fwh2, 'in1')
            fwh2_pu2 = connection(fwh2, 'out1', pu2, 'in1')
            pu2_me2 = connection(pu2, 'out1', me2, 'in2')

            sp2_fwh1 = connection(sp2, 'out2', fwh1, 'in1')
            fwh1_pu3 = connection(fwh1, 'out1', pu3, 'in1')
            pu3_me = connection(pu3, 'out1', me, 'in2')

            self.nw.add_conns(sp1_dsh, dsh_fwh2, fwh2_pu2, pu2_me2,
                              sp2_fwh1, fwh1_pu3, pu3_me)

            # busses
            # power bus
            self.power = bus('power')
            self.power.add_comps(
                {'comp': hpt, 'char': -1}, {'comp': mpt, 'char': -1},
                {'comp': lpt, 'char': -1}, {'comp': pu1, 'char': -1},
                {'comp': pu2, 'char': -1}, {'comp': pu3, 'char': -1})

            # heating bus
            self.heat = bus('heat')
            self.heat.add_comps(
                {'comp': eco, 'char': 1}, {'comp': eva, 'char': 1},
                {'comp': sup, 'char': 1})

            self.nw.add_busses(self.power, self.heat)

            # parametrization
            # components
            hpt.set_attr(eta_s=0.9)
            mpt.set_attr(eta_s=0.9)
            lpt.set_attr(eta_s=0.9)

            pu1.set_attr(eta_s=0.8)
            pu2.set_attr(eta_s=0.8)
            pu3.set_attr(eta_s=0.8)

            eco.set_attr(pr=0.99)
            eva.set_attr(pr=0.99)
            sup.set_attr(pr=0.99)

            con.set_attr(pr1=0.99, pr2=0.99, ttd_u=5)
            fwh1.set_attr(pr1=0.99, pr2=0.99, ttd_u=5)
            fwh2.set_attr(pr1=0.99, pr2=0.99, ttd_u=5)
            dsh.set_attr(pr1=0.99, pr2=0.99)

            # connections
            eco_eva.set_attr(x=0)
            eva_sup.set_attr(x=1)

            cc_hpt.set_attr(m=200, T=650, p=100, fluid={'water': 1})
            hpt_sp1.set_attr(p=20)
            mpt_sp2.set_attr(p=3)
            lpt_con.set_attr(p=0.05)

            cwi_con.set_attr(T=20, p=10, fluid={'water': 1})

        def calculate_efficiency(self, x):
            # set extraction pressure
            self.nw.connections['extraction1'].set_attr(p=x[0])
            self.nw.connections['extraction2'].set_attr(p=x[1])

            self.nw.solve('design')

            for cp in self.nw.components.values():
                if isinstance(cp, condenser) or isinstance(cp, desuperheater):
                    if cp.Q.val > 0:
                        return np.nan
                elif isinstance(cp, pump):
                    if cp.P.val < 0:
                        return np.nan
                elif isinstance(cp, turbine):
                    if cp.P.val > 0:
                        return np.nan

            if self.nw.res[-1] > 1e-3 or self.nw.lin_dep:
                return np.nan
            else:
                return self.nw.busses['power'].P.val / self.nw.busses['heat'].P.val

Note, that you have to label all busses and connections you want to access
later on with PyGMO. In :code:`calculate_efficiency(self, x)` the variable
:code:`x` is a list containing your decision variables. This function returns
the cycle efficiency for a specific set of decision variables. The efficiency
is defined by the ratio of total power transferred (including turbines and
pumps) to steam generator heat input.

Additionally, we have to make sure, only the result of physically feasible
solutions is returned. In case we have infeasible solutions, we can simply
return :code:`np.nan`. An infeasible solution is obtained in case the power
of a turbine is positive, the power of a pump is negative or the heat exchanged
in any of the preheaters is positive. We also check, if the calculation does
converge.

.. math::

    \eta_\mathrm{th}=\frac{|\sum P|}{\dot{Q}_{sg}}

Creating your PyGMO-Model
^^^^^^^^^^^^^^^^^^^^^^^^^

The optimization in PyGMO starts by defining the problem. You can set the
number of objectives your problem has in :code:`get_nobj()`. The number of
constraints is set in :code:`get_nec()` (equality constraints) and
:code:`get_nic()` (inequality constraints). In :code:`get_bounds()` you set the
bounds of your decision variables. Finally, you define your fitness function
and constraints in :code:`fitness(self, x)`:

.. code-block:: python

    import pygmo as pg


    class optimization_problem():

        def fitness(self, x):
            f1 = 1 / self.model.calculate_efficiency(x)
            ci1 = -x[0] + x[1]
            print(x)
            return [f1, ci1]

        def get_nobj(self):
            """Return number of objectives."""
            return 1

        # equality constraints
        def get_nec(self):
            return 0

        # inequality constraints
        def get_nic(self):
            return 1

        def get_bounds(self):
            """Return bounds of decision variables."""
            return ([1, 1], [40, 40])

By default PyGMO minimizes the fitness function. Therefore we set the fitness
function f1 to the reciprocal of the cycle efficiency. We set one inequality
constraint so that the pressure of the first extraction has to be bigger than
the second one:

.. math::

    p_{e,1} > p_{e,2}

In PyGMO your inequality constraint has to be in form of <0:

.. math::

    - p_{e,1} + p_{e,2} < 0

We expect that the extraction pressure won't be more than 40 bar and not less
1 bar. Therefore we set the bounds of our decision variables:

.. math::

    1 bar < p_{e,1} < 40 bar\\
    1 bar < p_{e,2} < 40 bar


Run PyGMO-Optimization
^^^^^^^^^^^^^^^^^^^^^^

The following code shows how to run the PyGMO optimization.

.. code-block:: python

    optimize = optimization_problem()
    optimize.model = PowerPlant()
    prob = pg.problem(optimize)
    num_gen = 15

    pop = pg.population(prob, size=10)
    algo = pg.algorithm(pg.ihs(gen=num_gen))


With optimize you tell PyGMO which problem you want to optimize. In the class
:code:`optimization_problem()` we defined our problem be setting fitness
function and inequality constraint. With :code:`optimize.model` we set the
model we want to optimize. In our case we want to optimize the extraction
pressures in our instance of class :code:`PowerPlant`. Finally, our problem is
set in :code:`prob = pg.problem(optimize)`.

With :code:`pop` we define the size of each population for the optimization,
:code:`algo` is used to set the algorithm you want to use. A list of available
algorithms can be found in
`List of algorithms <https://esa.github.io/pygmo2/overview.html#list-of-algorithms>`_.
The choice of your algorithm depends on the type of problem. Have you set
equality or inequality constraints? Do you perform a single- or multi-objective
optimization?

We choose a population size of 10 individuals and want to carry out 15
generations. We can evolve the population generation by generation, e.g. using
a for loop. At the end, we print out the information of the best individual.

.. code-block:: python

    for gen in range(num_gen):
        print('Evolution: {}'.format(gen))
        print('Efficiency: {} %'.format(round(100 / pop.champion_f[0], 4)))
        pop = algo.evolve(pop)

    print()
    print('Efficiency: {} %'.format(round(100 / pop.champion_f[0], 4)))
    print('Extraction 1: {} bar'.format(round(pop.champion_x[0], 4)))
    print('Extraction 2: {} bar'.format(round(pop.champion_x[1], 4)))

In our run, we got:

.. code:: bash

    Efficiency: 44.852 %
    Extraction 1: 26.62 bar
    Extraction 2: 2.825 bar


.. figure:: api/_images/scatterplot_efficiency_optimization.svg
    :align: center

    Figure: Scatter plot for all individuals during the optimization.
