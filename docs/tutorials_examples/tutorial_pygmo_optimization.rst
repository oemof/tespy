cycle optimization example using PyGMO
---------------------------------------

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
In case of a rather simple power plant topology the task of finding optimized 
values for e.g. extraction pressures is still managable without any optimization 
tool. As the topology becomes more complexe and boundary conditions come into play 
the usage of additional tools is recommended. 
The following example is intended to show the usage of PyGMO in combination 
with TESPy to maximize the cycle efficiency of a power plant with two extractions.

The source code can be found at the `tespy_examples repository
<https://github.com/oemof/oemof-examples/tree/master/oemof_examples/tespy/clausius_rankine>`_.  


PyGMO
^^^^^

PyGMO (Python Parallel Global Multiobjective Optimizer) is a library that provides 
a large number of evolutionary optimization algroithms. PyGMO can be used to 
solve constrained, unconstrained, single objective and multi objective problems.

Install PyGMO
+++++++++++++

We recommend installation with conda, which is included with the Anaconda distribution of Python. 
You can install pygmo by opening the anaconda prompt and executing the following::

    conda install -c conda-forge pygmo
    
The installation of all available optimization solvers is included in this process.


Creating your TESPy-Model
^^^^^^^^^^^^^^^^^^^^^^^^^

It is necessary to use object oriented programing in PyGMO. Therefore we creat 
a creat a class PowerPlant which contains our TESPy-Model and a function to return 
the cycle efficiency:

.. code-block:: python

    from tespy.networks import network
    from tespy.components import (
    turbine, splitter, merge, condenser, pump, sink, source,
    heat_exchanger_simple, desuperheater, cycle_closer
    )
    from tespy.connections import connection, bus
    
    class PowerPlant():
        def __init__(self):
            TESPy_Model
        
        def calculate_efficiency(self,x):
            # set extraction pressure
            self.nw.connections['extraction1'].set_attr(p=x[0])
            self.nw.connections['extraction2'].set_attr(p=x[1])
            
            self.nw.solve('design')
            self.nw.save('extraction')
                    
            return self.nw.busses['power'].P.val/self.nw.busses['heat'].P.val
        
In calculate_efficiency(self, x) the variable x is a list containing your 
decision variables. 




Defining the optimization problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The optimization in PyGMO starts be defining the problem at hand. You can set 
the number of objectives your problem has in get_nobj(). The number of constraints 
is set in get_nec() (equality constraints) and get_nic() (inequality constraints). 
In get_bounds() you set the bounds of your decision variables. Finally, you define 
your fitness function and constraints in fitness(self, x):

.. code-block:: python

    import pygmo as pg
    
    class ThermalPowerPlant():
        def fitness(self, x):
            f1 = 1/self.model.calculate_efficiency(x)
            ci1 = -x[0]+x[1]
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
    
        # # integer dimension
        # def get_nix(self):
        #     return 0
    
        def get_bounds(self):
            """Return bounds of decision variables."""
            return ([1,1], [40,40])
    
        def get_name(self):
            """Return function name."""
            return ""
            
In this case we want to maximize the efficiency of our powerplant. TESPy








