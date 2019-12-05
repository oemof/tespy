Combustion engine
-----------------

We have added a cogeneration unit in version 0.0.5 of TESPy. The cogeneration unit is based on the combustion chamber, and additionally provides a power and two heating outlets,
heat losses can be taken into account as well. Power output, heat production and losses are all linked to the thermal input of the combustion engine by characteristic lines, which usually are provided by the manufacturer.
TESPy provides a set of predefined characteristics (documented in the :py:class:`characteristics module <tespy.components.characteristics.characteristics.default>`).

.. figure:: api/_images/cogeneration_unit.svg
    :align: center
	
    Figure: Topology of a cogeneration unit.
	
The characteristics take the power ratio (:math:`\frac{P}{P_{ref}}`) as argument. For a design case simulation the power ratio is always assumed to be equal to 1.
For offdesign calculation TESPy will automatically take the rated power from the design case and use it to determine the power ratio. Still it is possible to specify the rated power manually, if you like.

In contrast to other components, the cogeneration unit has several busses, which are accessible by specifying the corresponding bus parameter:

- TI (thermal input),
- Q (total heat output),
- Q1 and Q2 (heat output 1 and 2),
- QLOSS (heat losses) and
- P (power output).

If you want to add a bus to the example from the tespy_examples repository, your code could look like this:

.. code-block:: python

	chp = cmp.cogeneration_unit('chp')

	bus = con.bus('some label')
	# for thermal input
	bus.add_comps({'c': chp, 'p': 'TI'})
	# for power output
	bus.add_comps({'c': chp, 'p': 'P'})
	
Enjoy fiddling around with the source code of the `cogeneration unit <https://github.com/oemof/oemof-examples/tree/master/oemof_examples/tespy/cogeneration_unit>`_ in the examples repository!
