.. _faq_label:

Frequently Asked Questions
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. dropdown:: Why is my model unstable or not convering towards a solution?

    - Short introduction into sources of numerical instability
    - Point towards :ref:`Starting Value tutorial <tutorial_starting_values_label>`
    - Discuss things like :code:`init_path`, :code:`design_path` etc. for off-design simulations

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

.. dropdown:: How can I define custom rules or behavior for my model?

    - Introduce :py:class:`tespy.tools.helpers.UserDefinedEquation`
    - Maybe small example here
    - Point towards full :ref:`UDE Example <ude_label>` in docs
    - One step further is to adjust existing or create new components (:ref:`Custom component tutorial <components_custom_components_label>` and :ref:`Documentation on how to contribute <developing_label>`)

.. dropdown:: How can I visualize my results?

    - Introduce and link to :ref:`FluProDia Example <fluprodia_label>`

.. dropdown:: How can I analyze potentials and improve upon my model?

    - Introduce exergy analysis and how it can help
    - Point towards `ExerPy <https://exerpy.readthedocs.io/en/latest/>`__
    - Reference other methods like entropy analysis etc.

.. dropdown:: Can I optimize certain properties of my model?

    - Introduce the optimization API :py:class:`tespy.tools.optimization.OptimizationProblem`
    - Point towards full :ref:`Optimization Example <tutorial_pygmo_optimization_label>`

.. dropdown:: Which HeatExchanger model is best-suited for me?

    - Point towards :ref:`Heat Exchanger Overview <tutorial_heat_exchanger>`
