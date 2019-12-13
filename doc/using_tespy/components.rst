TESPy components
================

In this section we will introduce you into the details of component
parametrisation and component characteristics. At the end of the section we
show you, how to create custom components.

List of components
------------------

More information on the components can be gathered from the code documentation.
We have linked the base class containing a figure and basic informations as
well as the equations.

- Basics
    * :py:class:`Cycle closer<tespy.components.basics.cycle_closer>`
    * :py:class:`Sink <tespy.components.basics.sink>`
    * :py:class:`Source <tespy.components.basics.source>`
    * :py:class:`Subsystem interface <tespy.components.basics.subsystem_interface>`
- Combustion
    * :py:class:`Combustion chamber <tespy.components.combustion.combustion_chamber>`
    * :py:class:`Combustion chamber stoichiometric <tespy.components.combustion.combustion_chamber_stoich>`
    * :py:class:`Combustion engine <tespy.components.combustion.combustion_engine>`
- Heat exchangers
    * :py:class:`Heat exchanger simple <tespy.components.heat_exchangers.heat_exchanger_simple>`
    * :py:class:`Heat exchanger <tespy.components.heat_exchangers.heat_exchanger>`
    * :py:class:`Condenser <tespy.components.heat_exchangers.condenser>`
    * :py:class:`Desuperheater <tespy.components.heat_exchangers.desuperheater>`
- Nodes
    * :py:class:`Node <tespy.components.nodes.node>`
    * :py:class:`Drum <tespy.components.nodes.drum>`
    * :py:class:`Merge <tespy.components.nodes.merge>`
    * :py:class:`Separator <tespy.components.nodes.separator>`
    * :py:class:`Splitter <tespy.components.nodes.splitter>`
- Piping
    * :py:class:`Pipe <tespy.components.piping.pipe>`
    * :py:class:`Valve <tespy.components.piping.valve>`
- Reactors
    * :py:class:`Water electrolyzer <tespy.components.reactors.water_electrolyzer>`
- Turbomachinery
    * :py:class:`Compressor <tespy.components.turbomachinery.compressor>`
    * :py:class:`Pump <tespy.components.turbomachinery.pump>`
    * :py:class:`Turbine <tespy.components.turbomachinery.turbine>`

.. _using_tespy_components_parametrisation_label:

Component parametrisation
-------------------------

Component parameters can be set and accessed in various ways. All parameters of
components are objects of a :code:`data_container` class. The data container
for component parameters it is called :code:`dc_cp`, :code:`dc_cc` for
component characteristics and :code:`dc_cm` for characteristic maps. The main
purpose of having a data container for the parameters (instead of pure
numbers), is added flexibility for the user. There are different ways for you
to specify a component parameter.

Component parameters
^^^^^^^^^^^^^^^^^^^^

The example shows different ways to specify the heat transfer coefficient of an
evaporator and how to unset the parameter again.

.. code-block:: python

    from tespy.components import heat_exchanger
    from tespy.tools import dc_cp
    import numpy as np

    he = heat_exchanger('evaporator')

    # specify the value
    he.set_attr(kA=1e5)
    # create a data container
    he.set_attr(kA=dc_cp(val=1e5, is_set=True))
    # set data container parameters
    he.kA.set_attr(val=1e5, is_set=True)

    # unset value
    he.set_attr(kA=np.nan)
    he.kA.set_attr(is_set=False)


Grouped parameters
^^^^^^^^^^^^^^^^^^

Grouped parameters are used whenever a component property depends on multiple
parameters. For instance, the pressure loss calculation via Darcy-Weissbach
requires information about the length, diameter and roughness of the pipe.
The solver will prompt a warning, if you do not specify all parameters required
by a parameter group. If parameters of the group are missing, the equation will
not be implemented by the solver.

.. code-block:: python

    from tespy.components import pipe
    import numpy as np

    my_pipe = pipe('pipe')

    # specify grouped parameters
    my_pipe.set_attr(D=0.1, L=20, ks=0.00005)

    # the solver will not use the Darcy-Weissbach-equation in this case
    my_pipe.set_attr(D=0.1, ks=0.00005)

There are three components using parameter groups:

- heat_exchanger_simple and pipe
    * :code:`hydro_group` (:code:`D`, :code:`L`, :code:`ks`)
    * :code:`kA_group` (:code:`kA`, :code:`Tamb`)
- solar_collector
    * :code:`hydro_group` (:code:`D`, :code:`L`, :code:`ks`)
    * :code:`energy_group` (:code:`E`, :code:`eta_opt`, :code:`lkf_lin`,
      :code:`lkf_quad`, :code:`A`, :code:`Tamb`)

Custom variables
^^^^^^^^^^^^^^^^

It is possible to use component parameters as variables of your system of
equations. For example, give a pressure ratio :code:`pr`, length :code:`L` and
roughness :code:`ks` of a pipe you want to calculate the pipe's diameter
:code:`D` required to achieve the specified pressure ratio. In this case you
need to specify the diameter the following way.

.. code-block:: python

    from tespy.components import pipe
    import numpy as np

    # custom variables
    my_pipe = pipe('my pipe')

    # make diameter variable of system
    my_pipe.set_attr(pr=0.98, L=100, ks=0.00002, D='var')
    
    # a second way of specifying this is similar to the
    # way used in the component parameters section
    # the benefit is, that val will be the starting value
    my_pipe.set_attr(pr=0.98, L=100, ks=0.00002)
    my_pipe.set_attr(D=dc_cp(val=0.2, is_set=True, is_var=True))
    
It is also possible to set value boundaries for you custom variable. You can do
this, if you expect the result to be within a specific range. But beware: This
might result in a non converging simulation, if the actual value is out of your
specified range.

.. code-block:: python
    
    # data container specification with identical result,
    # benefit: specification of bounds will increase stability
    my_pipe.set_attr(D=dc_cp(val=0.2, is_set=True, is_var=True,
                             min_val=0.1, max_val=0.3))

REVIEW THE COMPONENT CHARACTERISTICS PART
      
Component characteristics
^^^^^^^^^^^^^^^^^^^^^^^^^

Component characteristics are implemented for parameters in 
:ref:`several components <component_characteristics_label>`.

.. code-block:: python

    from tespy.components import heat_exchanger
    from tespy.tools import dc_cc
    import numpy as np

    he = heat_exchanger('evaporator')

    # specify name of predefined method
    he.set_attr(kA_char1='EVA_HOT')
    he.set_attr(kA_char2='EVA_COLD')

    # specify data container (yields same result)
    he.set_attr(kA_char1=dc_cc(method='EVA_HOT', param='m'))

    # specify data container (custom interpolation points x and y)
    x = np.array([0, 0.5, 1, 2])
    y = np.array([0, 0.8, 1, 1.2])
    he.set_attr(kA_char1=dc_cc(param='m', x=x, y=y))

.. _component_characteristics_label:

Component characteristics
-------------------------

Characteristics are available for the following components and parameters:

- pump
    * :py:meth:`eta_s_char <tespy.components.components.pump.eta_s_char_func>`: isentropic efficiency vs. volumetric flow rate.
    * :py:meth:`flow_char <tespy.components.components.pump.flow_char_func>`: pressure rise vs. volumetric flow characteristic.
- compressor
    * :py:meth:`char_map <tespy.components.components.compressor.char_map_func>`: component map for isentropic efficiency and pressure rise.
    * :py:meth:`eta_s_char <tespy.components.components.compressor.eta_s_char_func>`: isentropic efficiency vs. pressure ratio.
- turbine
    * :py:meth:`eta_s_char <tespy.components.components.turbine.eta_s_char_func>`: isentropic efficiency vs. isentropic enthalpy difference/pressure ratio/volumetric flow/mass flow.
- heat exchangers:
    * :py:meth:`kA1_char, kA2_char <tespy.components.components.heat_exchanger.kA_func>`: heat transfer coefficient, various predefined types, mass flows as specification parameters.
- simple heat exchangers
    * :py:meth:`kA_char <tespy.components.components.heat_exchanger_simple.kA_func>`: e. g. pipe, see heat exchangers
- cogeneration unit
    * :py:meth:`tiP_char <tespy.components.components.cogeneration_unit.tiP_char_func>`: thermal input vs. power ratio.
    * :py:meth:`Q1_char <tespy.components.components.cogeneration_unit.Q1_char_func>`: heat output 1 vs. power ratio.
    * :py:meth:`Q2_char <tespy.components.components.cogeneration_unit.Q2_char_func>`: heat output 2 vs. power ratio.
    * :py:meth:`Qloss_char <tespy.components.components.cogeneration_unit.Qloss_char_func>`: heat loss vs. power ratio.

You can specify the name of a default characteristic line or you define the whole data container for this parameter. The default characteristic lines can be found in the :py:mod:`documentation <tespy.components.characteristics>`.

.. code-block:: python

    from tespy.components import turbine, heat_exchanger
    from tespy.tools import dc_cc

    turb = turbine('turbine')
    # method specification (default characteristic line "TRAUPEL")
    turb.set_attr(eta_s_char='TRAUPEL')
    # data container specification
    turb.set_attr(eta_s_char=dc_cc(method='TRAUPEL', param='dh_s', x=None, y=None))

    # defining a custom line (this line overrides the default characteristic line, method does not need to be specified)
    x = np.array([0, 1, 2])
    y = np.array([0.95, 1, 0.95])
    turb.set_attr(eta_s_char=dc_cc(param='dh_s', x=x, y=y)

    # heat exchanger analogously
    he = heat_exchanger('evaporator')
    he.set_attr(kA_char1='EVA_HOT')
    he.set_attr(kA_char2='EVA_COLD')

Custom components
-----------------

You can add own components. The class should inherit from
:py:class:`tespy.components.components.component class <tespy.components.components.component>`
or its children. In order to do that, create a python file in your working
directory and import the base class for your custom component, e. g. the
:py:mod:`tespy.components.components.component <tespy.components.components>`
class. Now add create a class for your component and at least add the following
methods.

- :code:`attr(self)`,
- :code:`inlets(self)`,
- :code:`outlets(self)`,
- :code:`equations(self)`,
- :code:`derivatives(self)` as well as
- :code:`calc_parameters(self)`.

The starting lines of your file would look like this:

.. code:: python

    from tespy.components.components import component
    from tespy.tools import dc_cc, dc_cp

    class my_custom_component(component):
        """
        This is a custom component.
        
        You can add your documentation here. From this part, it should be clear
        for the user, which parameters are available, which mandatory equations
        are applied and which optional equations can be applied using the
        component parameters.
        """


Attributes
^^^^^^^^^^

The attr method must return a dictionary with the attributes you want to
specify. The keys represent the attributes and the respective values the type
of data_container used. Using the data_container attributes it is possible to
add defaults. Defaults for characteristic lines or characteristic maps are
loaded automatically by the component initialisation method of class
:code:`component`.

ADD DEFAULT CHAR DESCRIPTION

.. code:: python

    def attr(self):
        return {'par1': dc_cp(min_val=0, max_val=1),
                'par2': dc_cc(param='m')}


Inlets and outlets
^^^^^^^^^^^^^^^^^^

:code:`inlets(self)` and :code:`outlets(self)` respectively must return a list
of strings. The list may look like this:

.. code:: python

    def inlets(self):
        return ['in1', 'in2']

    def outlets(self):
        return ['out1', 'out2']

The number of inlets and outlets might even be generic, e. g. if you have added
an attribute :code:`'num_in'` in :code:`attr(self)`:

.. code:: python

    def inlets(self):
        if self.num_in_set:
            return ['in' + str(i + 1) for i in range(self.num_in)]
        else:
            self.set_attr(num_in=2)
            return self.inlets()

Equations
^^^^^^^^^

The equations contain the information on the changes to the fluid properties within the component. Each equation must be defined in a way, that the correct result is zero, e. g.:

.. math::

    0 = \dot{m}_{in} - \dot{m}_{out}\\
    0 = \dot{p}_{in} - \dot{p}_{out} - \Delta p

The connections connected to your component are available as a list in :code:`self.inl` and :code:`self.outl` respectively.

.. code:: python

    def equations(self):

        vec_res = []

        vec_res += [self.inl[0].m.val_SI - self.outl[0].m.val_SI]
        vec_res += [self.inl[0].p.val_SI - self.outl[0].p.val_SI - self.dp.val]

The equations are added to a list one after another, which will be returned at the end.

Derivatives
^^^^^^^^^^^

You need to calculate the partial derivatives of the equations to all variables of the network.
This means, that you have to calculate the partial derivatives to mass flow, pressure, enthalpy and all fluids in the fluid vector on each incomming or outgoing connection of the component.

Add all derivatives to a list (in the same order as the equations) and return the list as numpy array (:code:`np.asarray(list)`).
The derivatives can be calculated analytically or numerically by using the inbuilt function :code:`numeric_deriv(self, func, dx, pos, **kwargs)`.

- :code:`func` is the function you want to calculate the derivatives for,
- :code:`dx` is the variable you want to calculate the derivative to and
- :code:`pos` indicates the connection you want to calculate the derivative for, e. g. :code:`pos=1` means, that counting your inlets and outlets from low index to high index (first inlets, then outlets),
  the connection to be used is the second connection in that list.
- :code:`kwargs` are additional keyword arguments required for the function.

For a good start just look into the source code of the inbuilt components. If you have further questions feel free to contact us.
