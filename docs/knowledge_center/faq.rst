.. _faq_label:

Frequently Asked Questions
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. container:: accordion-group

    .. dropdown:: Why does my simulation model crash with a CoolProp related error or does not converge towards a solution?

        - Short introduction into sources of numerical instability
        - Point towards :ref:`Starting Value tutorial <tutorial_starting_values_label>`
        - Discuss things like :code:`init_path`, :code:`design_path` etc. for off-design simulations
        - Best practices for automatically making small steps towards a changed input

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

        - Introduce :py:class:`tespy.tools.helpers.UserDefinedEquation`
        - Maybe small example here
        - Point towards full :ref:`UDE Example <ude_label>` in docs
        - One step further is to adjust existing or create new components (:ref:`Custom component tutorial <components_custom_components_label>` and :ref:`Documentation on how to contribute <developing_label>`)

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
        `fluprodia's online documentation <https://https://fluprodia.readthedocs.io/en/latest/>`__.

    .. dropdown:: How can I plot Q-T diagrams of heat exchangers?

        If you want to analyze the heat transfer within a heat exchangerk, e.g.
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

        - Introduce exergy analysis and how it can help
        - `exerpy <https://exerpy.readthedocs.io/en/latest/>`__
        - Reference other methods like entropy analysis etc.

    .. dropdown:: Is it possible to create an optimization problem with my model?

        - Introduce the optimization API :py:class:`tespy.tools.optimization.OptimizationProblem`
        - Point towards full :ref:`Optimization Example <tutorial_pygmo_optimization_label>`

    .. dropdown:: Which of the available HeatExchanger classes should I use?

        - Point towards :ref:`Heat Exchanger Overview <tutorial_heat_exchanger>`
