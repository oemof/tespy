.. _faq_label:

Frequently Asked Questions
~~~~~~~~~~~~~~~~~~~~~~~~~~

This section contains answers to some frequently arising questions, which are
grouped by theme and can be navigated with the sidebar to the right. If you can
not find a question that relates to your specific problem, please consult the
full API documentation or search for existing questions in the
`discussion tab <https://github.com/oemof/tespy/discussions>`__ on GitHub. You
can also submit a new question there, if you feel like your specific problem
could affect more users than yourself. Furthermore, if you encounter a bug,
please submit an `issue <https://github.com/oemof/tespy/issues>`__ on GitHub.

.. _faq_best_practices_label:

Modeling best practices
-----------------------

.. container:: accordion-group

    .. dropdown:: Which of the available HeatExchanger classes should I use?

        There are multiple different types of :code:`HeatExchanger` available
        in TESPy. These are

        - Energy balance heat exchangers, which only consider one side of the
          heat exchange
        - 0D models which transfer heat from one fluid to the other and
        - discretized 1D models, which allow the specification of internal
          pinch.

        We have created an :ref:`overview <tutorial_heat_exchanger>` on the
        different 0D and 1D types available and when to use which model.

.. _faq_errors_label:

Instability, Errors and Debugging
---------------------------------

.. container:: accordion-group

    .. dropdown:: Why does my simulation model crash with a CoolProp related error or does not converge towards a solution?

        TESPy uses
        `Newton's Method <https://en.wikipedia.org/wiki/Newton%27s_method>`__
        in order to numerically solve the system of equations created according
        to your model design and parametrization. The algorithm iteratively
        approaches a solution for the unknown variables within your model with
        the help of the derivatives of the equations. Since the algorithm
        itself is not aware of the physical boundaries limiting inputs to
        CoolProp such errors may occur, most of the time due to the following
        two reasons:

        - The model inputs or the system design are not reasonable.
        - The initial guess for the variables is far away from the actual
          solution.

        While there is little you can do to prevent the first one from the
        software side (we recommend you reach out in the forum or get into
        exchange with your friends or colleagues), the initial value issue can
        be tackled:

        If you set up a model without providing any start values, these are
        automatically generated. If those do not fit your specific model,
        instability in the solving process can be the result and lead to
        CoolProp errors.

        The
        :ref:`tutorial about starting values <tutorial_starting_values_label>`
        is a good introduction into a systematic approach to prevent
        non-converging modelling. Furthermore, the network's
        :py:meth:`.solve() <tespy.networks.Network.solve>` method automatically
        starts at the previous solution. If your new inputs are far away from
        that, you can iteratively change the inputs to arrive at the desired
        destination. Furthermore, it is possible to provide starting values
        from a saved simulation through the :code:`init_path`. Find a more
        thorough description about how the solving process is implemented and
        configured in TESPy in the respective chapter in the
        :ref:`Network documentation <networks_solving_label>`.
        It also contains a comprehensive section about
        :ref:`convergence stability <module_convergence_label>`.

    .. dropdown:: How do I know, which equations are applied in my model?

        You can use the :ref:`debugging <tutorial_debugging_label>` functions
        to identify the variables and the equations in your model. The
        equations available are displayed in a compact overview in separate
        tables for all :ref:`component classes <modules_components_label>` and
        for :ref:`connections <modules_connections_label>`.

    .. dropdown:: How do I handle errors I am seeing about my variable specifications in the presolve phase?

        If you are seeing an error similar to this:

        .. error::

            .. code-block:: bash

                TESPyNetworkError: You specified more than one variable of the linear dependent variables: (e1: E), (e2: E).

        there is something wrong with the variables you specified in your
        model. The specific error above is caused by setting two individual
        variables, which are dependent on each other through some linear
        relation and therefore can not both be variables. Other reasons for
        similar errors are circular dependency of a number of variables as well
        as over-determination of a single variable via multiple equations.
        These errors are explained and handled for an example model in the
        :ref:`debugging tutorial about presolve errors <tutorial_debugging_error_presolve_label>`.

    .. dropdown:: How do I deal with a singularity in the Jacobian matrix?

        If you are seeing an error similar to this:

        .. error::

            .. code-block:: bash

                Detected singularity in Jacobian matrix. This singularity is most likely caused by the parametrization of your problem and NOT a numerical issue. Double check your setup.
                The following variables of your problem are not in connection with any equation: (2, 'p')

        there is something wrong with the variables you specified in your
        model or a numerical problem with the partial derivatives of the
        underlying equations. The
        :ref:`debugging tutorial about linear dependency <tutorial_debugging_linear_dependency_label>`
        illustrates how you can inspect your model for reasons that cause
        linear dependencies resulting in a singularity in the Jacobian matrix.

.. _faq_customization_label:

Customizing model behavior
--------------------------

.. container:: accordion-group

    .. dropdown:: How do I create a relation between two unknown variables?

        You can use the :code:`Ref` specification on connection properties for
        that. Consider a splitter, where you do not know the absolute mass flow
        at the inflow, but you want to fix a certain split ratio.

        .. code-block:: python

            >>> from tespy.components import Source, Sink, SimpleHeatExchanger, Splitter
            >>> from tespy.connections import Connection, Ref
            >>> from tespy.networks import Network

            >>> nw = Network(iterinfo=False)
            >>> nw.units.set_defaults(
            ...     temperature="°C",
            ...     pressure="bar",
            ...     heat="kW"
            ... )

            >>> water = Source("water inflow")
            >>> boiler = SimpleHeatExchanger("boiler")
            >>> splitter = Splitter("splitter")
            >>> steam1 = Sink("steam 1")
            >>> steam2 = Sink("steam 2")

            >>> c1 = Connection(water, "out1", boiler, "in1", label="c1")
            >>> c2 = Connection(boiler, "out1", splitter, "in1", label="c2")
            >>> c3 = Connection(splitter, "out1", steam1, "in1", label="c3")
            >>> c4 = Connection(splitter , "out2", steam2, "in1", label="c4")

            >>> nw.add_conns(c1, c2, c3, c4)

            >>> c1.set_attr(fluid={"water": 1}, p=50, T=25)
            >>> c2.set_attr(x=1)
            >>> c3.set_attr(m=Ref(c2, 0.9, 0))  # c3 mass flow is 90 % of c2 mass flow
            >>> boiler.set_attr(Q=100, dp=0)

            >>> nw.solve("design")
            >>> round(c3.m.val / c2.m.val, 1)
            0.9

    .. dropdown:: How can I add custom constraints or equations to my model?

        You can make use of the :py:class:`UserDefinedEquation <tespy.tools.helpers.UserDefinedEquation>`.
        With this class you can create a function that returns the residual
        value of your equation and define what variables the function depends
        on. Everything else is cared for by the solver. For a couple of
        examples you can also check :ref:`this page <ude_label>` in docs. Using
        it is the simple entry point before customizing existing or developing
        new components.

    .. dropdown:: I have two compressors in a model, both should have the same isentropic efficiency as a result.

        This is possible with the :ref:`UserDefinedEquation <ude_label>`. In a
        :code:`UserDefinedEquation` you can build arbitrary equations and apply
        them to your model by the variables registered in the solver. The
        variables are mostly representing connection properties,
        :ref:`some component properties <modules_components_label>` can also be
        variables, the isentropic efficiency of a compressor is not one of
        these. But, you can still implement an equation for this problem, by
        writing it as function of the pressure and enthalpy values involved:

        .. code-block::

            >>> from tespy.tools.fluid_properties import isentropic


            >>> def equal_efficiency_ude(ude):
            ...     cp1, cp2 = ude.comps
            ...
            ...     cp1_in, cp1_out, cp2_in, cp2_out = cp1.inl[0], cp1.outl[0], cp2.inl[0], cp2.outl[0]
            ...     return (
            ...         (
            ...             isentropic(
            ...                 cp1_in.p.val_SI,
            ...                 cp1_in.h.val_SI,
            ...                 cp1_out.p.val_SI,
            ...                 cp1_in.fluid_data,
            ...                 cp1_in.mixing_rule
            ...             ) - cp1_in.h.val_SI
            ...         ) * (cp2_out.h.val_SI - cp2_in.h.val_SI)
            ...         - (
            ...             isentropic(
            ...                 cp2_in.p.val_SI,
            ...                 cp2_in.h.val_SI,
            ...                 cp2_out.p.val_SI,
            ...                 cp2_in.fluid_data,
            ...                 cp2_in.mixing_rule
            ...             ) - cp2_in.h.val_SI
            ...         ) * (cp1_out.h.val_SI - cp1_in.h.val_SI)
            ...     )

            >>> def equal_efficiency_ude_dependents(ude):
            ...     cp1, cp2 = ude.comps
            ...     cp1_in, cp1_out, cp2_in, cp2_out = cp1.inl[0], cp1.outl[0], cp2.inl[0], cp2.outl[0]
            ...     return [
            ...         var for c in [cp1_in, cp1_out, cp2_in, cp2_out]
            ...         for var in [c.p, c.h]
            ...     ]

    .. dropdown:: I have custom fluid property data, how can I integrate them into my model?

        In the documentation section
        :ref:`on fluid property engines <fluid_properties_label>` you will find
        a lot of helpful information. Specifically for liquids/incompressibles
        there is already an interface, that takes your datatables and then fits
        equations to them and integrates them  into your tespy model. Check out
        the information on that topic in
        :ref:`this section <incompressible_wrapper_label>`. You can also
        implement your own class, that handles your fluid property equations.
        Follow the examples given in the sections mentioned to learn, how that
        can be accomplished.

.. _faq_postprocessing_label:

Visualization, post-processing and cycle analysis
-------------------------------------------------

.. container:: accordion-group

    .. dropdown:: How can I plot log(p)-h and T-s diagrams from my simulation?

        Fluid property diagrams can help you analyze and improve your
        thermodynamic processes. CoolProp offers a somewhat unwieldly
        `inbuilt feature <https://coolprop.org/coolprop/python-plotting.html>`__
        for creating those, but we recommend using the officially supported
        companion module fluprodia. It allows for the creation and
        customization of different types of diagrams for all pure and
        pseudo-pure fluids available in CoolProp. There is an introductory
        :ref:`FluProDia Example <fluprodia_label>` that shows its usage as well
        as TESPy's API to extract the necessary plotting data. A thorough
        explanation can be found in
        `fluprodia's online documentation <https://fluprodia.readthedocs.io/en/latest/>`__.

    .. dropdown:: How can I plot Q-T diagrams of heat exchangers?

        If you want to analyze the heat transfer within a heat exchanger, e.g.
        to investigate the phase change of a working fluid, creating a Q-T
        diagram can be very helpful. While the heat transfer of simple heat
        exchanger models like the
        :py:class:`HeatExchanger <tespy.components.heat_exchangers.base.HeatExchanger>`
        or :py:class:`Condenser <tespy.components.heat_exchangers.condenser.Condenser>`
        is easily calculated by hand, this can not be said about more complex
        variants like the
        :py:class:`SectionedHeatExchanger <tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger>`
        or :py:class:`MovingBoundaryHeatExchanger <tespy.components.heat_exchangers.movingboundary.MovingBoundaryHeatExchanger>`.
        Fortunately, the latter possess the
        :py:meth:`.calc_sections() <tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger.calc_sections>`
        method that returns information about their heat transfer, including
        the hot and cold side temperatures at different steps. Below is an
        example snippet of how to generate a Q-T diagram from the results of a
        :py:class:`SectionedHeatExchanger <tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger>`
        /:py:class:`MovingBoundaryHeatExchanger <tespy.components.heat_exchangers.movingboundary.MovingBoundaryHeatExchanger>`.

        .. code-block:: python

            from matplotlib import pyplot as plt

            ...

            heatex = SectionedHeatExchanger("heat exchanger")  # or MovingBoundaryHeatExchanger()

            ...

            nw.solve("design")

            heat, T_hot, T_cold, heat_per_section, td_log_per_section = heatex.calc_sections()
            heat /= 1e6

            fig, ax = plt.subplots(figsize=(10, 6))

            ax.plot(heat, T_hot, "o-", color="red")
            ax.plot(heat, T_cold, "o-", color="blue")

            ax.set_ylabel("temperature in K")
            ax.set_xlabel("heat transferred in MW")

            plt.show()

    .. dropdown:: How can I apply Second Law of thermodynamics analysis methods?

        TESPy has an integrated exergy analysis feature. This was developed
        originally in two publications: :cite:`Witte2022,Hofmann2022`. The
        feature automatically calculates the exergy of all streams in your
        model and applies the exergy balance equations for all components to
        make the analysis of the complete model. Recently, a new project to
        extend this work has started. For this, the feature is extracted from
        tespy and implemented in an own library:
        `exerpy <https://exerpy.readthedocs.io/en/latest/>`__. exerpy is
        compatible with TESPy and you can seamlessly pass your TESPy model to
        the exergy analysis methods.

    .. dropdown:: Is it possible to create an optimization problem with my model?

        You can couple your simulation model with any kind of optimization
        library for non-linear optimization e.g. scipy, pygmo, pymoo etc.. For
        the coupling with pymoo there is a dedicated API available. To use it
        you have to provide your model in the form of a
        :ref:`model class <integration_model_class_template_label>`. This will
        let you use the pymoo integration as described in the
        :ref:`Optimization Example <tutorial_optimization_label>`.
        Furthermore, see the publication by Chen et al. :cite:`Chen2022` for a
        thorough application example in the scientific context.
