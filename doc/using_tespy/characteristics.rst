TESPy characteristics
=====================

In this part we give an introduction on the TESPy characteritics
implementation. There two different types of characteristics available in
TESPy: lines (:py:class:`char_line <tespy.tools.characteristics.char_line>`)
and maps (:py:class:`char_line <tespy.tools.characteristics.char_map>`).

Characteristic lines
--------------------

The characteristic lines use linear interpolation in order to determine the
y-value given a x-value as input, where x is between the lower boundary x0 and
the upper boundary x1. y0 and y1 are the corresponding y-values.

.. math::

    y = y0 + \frac{x-x0}{x1-x0} \cdot \left(y1-y0 \right)
	
If the x value is above the maximum x value or below the minium x value of the
characteristic line the y value corresponding to the the maximum/minimum value
is returned instead.

Characteristic maps
-------------------

.. _import_custom_characteristics_label:

Import your own characteristics
-------------------------------

It is possible to import your own characteristic lines or maps instead of
writing the x, y (z1 and z2) data into your python script. Go to your .tespy
folder in your HOME directory. Create a folder named :code:`data`, if it does
not exist. In this folder you can place two json-files for your
characteristics.

- :code:`char_lines.json`
- :code:`char_maps.json`

Your custom definitions of characteristic lines go into the
:code:`char_lines.json` and your characteristic map definitions go into the
:code:`char_maps.json` document.

The :code:`char_lines.json` should have names for identification of the
characteristic lines on the first level. On the second level the x and y data
are assigned to the name of the characteristic line. The x and y data must be
stated as list.

.. code-block:: json

    {
        "name_of_char_line": {
            "x": [0, 0.5, 1, 1.5, 2],
            "y": [0.8, 0.9, 1, 1.1, 1.2]
        },
        "name_of_2nd_char_line": {
            "x": [0, 0.5, 1, 1.5, 2],
            "y": [2, 1.1, 1, 1.2, 1.7]
        },
        ...
        "name_of_last_char_line": {
            "x": [0, 0.5, 1, 1.5, 2],
            "y": [0.8, 0.95, 1, 0.95, 0.8]
        }
    }

The :code:`char_maps.json` should also have names for identification of the
characteristic lines on the first level. On the second level we aditionally
need z1 and z2 data. The x data are a list of values, the y, z1 and z2 data
are arrays with a list of values for each dimension of the x data. The example
below has 3 x values, thus the y, z1 and z2 data must contain 3 sets of values.

.. code-block:: json

    {
        "name_of_char_map": {
                "x": [0.971, 1, 1.029],
                "y": [[0.93, 0.943, 0.953, 0.961, 0.962, 0.963],
                      [0.987, 0.995, 1.0, 1.002, 1.005, 1.005],
                      [1.02, 1.023, 1.026,1.028, 1.03, 1.032]],
                "z1": [[0.982, 0.939, 0.895, 0.851, 0.806, 0.762],
                       [1.102, 1.052, 1.0, 0.951, 0.9, 0.85],
                       [1.213, 1.149, 1.085, 1.022, 0.958, 0.894]],
                "z2": [[0.981, 0.995, 1.007, 1.002, 0.981, 0.961],
                       [0.969, 0.984, 1.0, 0.985, 0.967, 0.95],
                       [0.962, 0.949, 0.935, 0.922, 0.908, 0.895]]
            }
    }
