Characteristics
===============

In this part we give an introduction on the TESPy characteristics
implementation. There two different types of characteristics available in
TESPy: lines (:py:class:`char_line <tespy.tools.characteristics.char_line>`)
and maps (:py:class:`char_map <tespy.tools.characteristics.char_map>`).
The default characteristics available are to be found in the
:py:mod:`tespy.data <tespy.data>` module documentation.

Characteristic lines
--------------------

The characteristic lines use linear interpolation in order to determine the
y-value given a x-value as input, where x is between the lower boundary
:math:`x_0` and the upper boundary :math:`x_1`. :math:`y_0` and :math:`y_1` are
the corresponding y-values.

.. math::

    y = y_0 + \frac{x-x_0}{x_1-x_0} \cdot \left(y_1-y_0 \right)

It is possible to specify an :code:`extrapolate` parameter. If the value is
:code:`False` (default state) and the x-value is above the maximum or below the
minimum value of the characteristic line the y-value corresponding to the the
maximum/minimum value is returned instead. If the :code:`extrapolate` is
:code:`True` linear extrapolation is performed using the two lower most or
upper moste value pairs respectively.

Characteristic maps
-------------------

The characteristic maps use linear interpolation as well. First step is
interpolation on the x-dimension similar to the characteristic line
functionality. As the y, z1 and z2 data are two-dimensional, **each row of**
**the data corresponds to one x value**. Thus, the result of the first step is
a vector for each dimesion (y, z1 and z2).

.. math::

    \vec{y} = \vec{y_0} + \frac{x-x_0}{x_1-x_0} \cdot \left(\vec{y_1}-
    \vec{y_0} \right)

    \vec{z1} = \vec{z1_0} + \frac{x-x_0}{x_1-x_0} \cdot \left(\vec{z1_1}-
    \vec{z1_0} \right)

    \vec{z2} = \vec{z2_0} + \frac{x-x_0}{x1-x_0} \cdot \left(\vec{z2_1}-
    \vec{z2_0}\right)

Using the y value as second input dimension the corresponding z1 and z2 values
are calculated, again using linear interpolation.

.. math::

    z1 = z1_0 + \frac{y-y_0}{y_1-y_0} \cdot \left(z1_1-z1_0 \right)

    z2 = z2_0 + \frac{y-y_0}{y_1-y_0} \cdot \left(z2_1-z2_0 \right)

.. _import_custom_characteristics_label:

Import your own characteristics
-------------------------------

It is possible to import your own characteristic lines or maps instead of
writing the x, y (z1 and z2) data into your python script, for example:

.. code-block:: python

    from tespy.tools.characteristics import load_custom_char, char_line

    gen_char = load_custom_char('generator', char_line)


For the imports to work in the way shown, navigate to your .tespy folder in
your HOME directory :code:`HOME/.tespy`. Create a folder named :code:`data`, if
it does not exist. In this folder, you can place two json-files for your
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
        "name_of_last_char_line": {
            "x": [0, 0.5, 1, 1.5, 2],
            "y": [0.8, 0.95, 1, 0.95, 0.8]
        }
    }

The :code:`char_maps.json` should also have names for identification of the
characteristic lines on the first level. On the second level we additionally
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
