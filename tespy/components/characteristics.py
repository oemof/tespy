# -*- coding: utf-8

"""
.. module:: components.characteristics
    :synopsis:

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>
"""

import numpy as np

import logging

# %%


class characteristics:
    r"""
    Class characteristics for components.

    Parameters
    ----------
    x : ndarray
        An array for the x-values of the lookup table. Number of x and y values must be identical.

    y : ndarray
        The corresponding y-values for the lookup table. Number of x and y values must be identical.

    method : str
        Specify a method to choose from the default characteristic lines. If you specify custom x and y values, this parameter will be ignored.

    comp : str
        Component base name, see :func:`tespy.components.components.component.comp` method.

    Note
    ----
    This class generates a lookup table from the given input data x and y, then performs linear interpolation.
    The x and y values may be specified by the user. There are some default characteristic lines for different
    components, see the :func:`tespy.components.characteristics.characteristics.default` method.
    If you neither specify the method to use from the defaults nor specify x and y values,
    the characteristic line generated will be :code:`x: [1, 2, 3, 4], y: [1, 1, 1, 1]`.
    """

    def __init__(self, **kwargs):

        for key in kwargs:
            if key not in self.attr():
                msg = ('Invalid keyword ' + key + '. Available keywords for kwargs are: ' + str(self.attr()) + '.')
                logging.error(msg)
                raise KeyError(msg)

        # method to read from default characteristic lines
        method = kwargs.get('method', 'default')

        self.x = kwargs.get('x', None)
        self.y = kwargs.get('y', None)
        self.comp = kwargs.get('comp', None)

        if self.x is None:
            self.x = self.default(method)[0]
        if self.y is None:
            self.y = self.default(method)[1]

        if isinstance(self.x, list):
            self.x = np.array(self.x)
        if isinstance(self.y, list):
            self.y = np.array(self.y)

        if len(self.x) != len(self.y):
            msg = ('Please provide the same amount of x-values and y-values. Number of x-values: ' +
                   str(len(self.x)) + ', number of y-values: ' + str(len(self.y)) + '.')
            logging.error(msg)
            raise ValueError(msg)

        msg = 'Created characteristic function for component of type ' + str(self.comp) + ' with default method ' + method +'.'
        logging.debug(msg)

    def default(self, key):
        r"""

        **default characteristic lines for turbines**

        .. math::

            \frac{\eta_\mathrm{s,t}}{\eta_\mathrm{s,t,ref}}=f\left(X \right)

        .. math::

            \text{choose calculation method for X}

            X = \begin{cases}
            \frac{\dot{m}}{\dot{m}_{ref}} & \text{mass flow}\\
            \frac{\dot{V}}{\dot{V}_{ref}} & \text{volumetric flow}\\
            \frac{p_1 \cdot p_{2,ref}}{p_{1,ref} \cdot p_2} &
            \text{pressure ratio}\\
            \sqrt{\frac{\Delta h_\mathrm{s,ref}}{\Delta h_\mathrm{s}}} &
            \text{isentropic enthalpy difference}
            \end{cases}

        **GENERIC**

        .. image:: _images/turbine_GENERIC.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

        **TRAUPEL**

        .. image:: _images/turbine_TRAUPEL.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

        **literature**

        - Walter Traupel (2001): Thermische Turbomaschinen Band 2. Berlin:
          Spinger.
          -> TRAUPEL

        **default characteristic lines for compressors**

        .. math::

            \frac{\eta_\mathrm{s,t}}{\eta_\mathrm{s,t,ref}}=f\left(X \right)

        .. math::

            \text{choose calculation method for X}

            X = \begin{cases}
            \frac{\dot{m}}{\dot{m}_{ref}} & \text{mass flow}\\
            \frac{p_1 \cdot p_{2,ref}}{p_{1,ref} \cdot p_2} &
            \text{pressure ratio}\\
            \end{cases}

        **GENERIC**

        .. image:: _images/compressor_GENERIC.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

        **default characteristic lines for pumps**

        .. math::

            \frac{\eta_\mathrm{s,t}}{\eta_\mathrm{s,t,ref}}=
            f\left(\frac{\dot{V}}{\dot{V}_{ref}} \right)

        **GENERIC**

        .. image:: _images/pump_GENERIC.svg
           :scale: 100 %
           :alt: alternative text
           :align: center


        **default characteristic lines for cogeneration units**

        .. math::

            \frac{X}{P}=f\left(\frac{P}{P_{ref}} \right)

        **thermal input** (TI)

        .. math::
            X = TI = \dot{m}_f \cdot LHV

        .. image:: _images/TI.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

        **heat production** (Q1)

        .. math::
            X = \dot{Q}_1 = \dot{m}_1 \cdot \left( h_{out,1} - h_{in,1} \right)

        .. image:: _images/Q1.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

        **heat production** (Q2)

        .. math::
            X = \dot{Q}_2 = \dot{m}_2 \cdot \left( h_{out,2} - h_{in,2} \right)

        .. image:: _images/Q2.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

        **heat loss** (QLOSS)

        .. math::
            X = \dot{Q}_{loss}

        .. image:: _images/QLOSS.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

        **default characteristic lines for heat exchangers**

        .. math::

            \frac{kA}{kA_\mathrm{ref}}=f_1\left(x_1 \right)
            \cdot f_2\left(x_2 \right)

        available characteristic lines:

        **condensing fluid** (COND)

        .. image:: _images/COND_HOT.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

        .. image:: _images/COND_COLD.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

        **economiser, evaporator, superheater** (EVA)

        .. image:: _images/EVA_HOT.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

        .. image:: _images/EVA_COLD.svg
           :scale: 100 %
           :alt: alternative text
           :align: center


        **heat exchanger without phase change** (HE)

        .. image:: _images/HE_HOT.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

        .. image:: _images/HE_COLD.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

        """

        if key == 'default':
            return np.array([0, 1, 2, 3]), np.array([1, 1, 1, 1])

        x = {}
        y = {}

        if self.comp == 'turbine':

            x['GENERIC'] = np.array(
                    [0.000, 0.300, 0.600, 0.700, 0.800, 0.900, 1.000, 1.100,
                     1.200, 1.300, 1.400, 1.500])
            y['GENERIC'] = np.array(
                    [0.950, 0.980, 0.993, 0.996, 0.998, 0.9995, 1.000, 0.999,
                     0.996, 0.990, 0.980, 0.960])

            x['TRAUPEL'] = np.array(
                    [0.0000, 0.1905, 0.3810, 0.5714, 0.7619, 0.9524, 1.0000,
                     1.1429, 1.3333, 1.5238, 1.7143, 1.9048])
            y['TRAUPEL'] = np.array(
                    [0.0000, 0.3975, 0.6772, 0.8581, 0.9593, 0.9985, 1.0000,
                     0.9875, 0.9357, 0.8464, 0.7219, 0.5643])

        elif self.comp == 'compressor':

            x['GENERIC'] = np.array(
                    [0.000, 0.400, 1.000, 1.200])
            y['GENERIC'] = np.array(
                    [0.500, 0.900, 1.000, 1.100])

        elif self.comp == 'pump':

            x['GENERIC'] = np.array(
                    [0.071, 0.282, 0.635, 0.776, 0.917, 1.000, 1.128, 1.270,
                     1.410, 1.763, 2.115, 2.500])
            y['GENERIC'] = np.array(
                    [0.250, 0.547, 0.900, 0.965, 0.995, 1.000, 0.990, 0.959,
                     0.911, 0.737, 0.519, 0.250])

        elif self.comp == 'cogeneration unit':

            x['TI'] = np.array([0.50, 0.75, 0.90, 1.00])
            y['TI'] = np.array([2.50, 2.33, 2.27, 2.25])

            x['Q1'] = np.array([0.660, 0.770, 0.880, 0.990, 1.100])
            y['Q1'] = np.array([0.215, 0.197, 0.185, 0.175, 0.168])

            x['Q2'] = np.array([0.660, 0.770, 0.880, 0.990, 1.100])
            y['Q2'] = np.array([0.215, 0.197, 0.185, 0.175, 0.168])

            x['QLOSS'] = np.array([0.50, 0.7500, 0.90, 1.000])
            y['QLOSS'] = np.array([0.32, 0.3067, 0.30, 0.295])

        elif (self.comp == 'heat exchanger' or self.comp == 'desuperheater' or
              self.comp == 'pipe' or self.comp == 'heat exchanger simple' or
              self.comp == 'condenser'):

            x['EVA_HOT'] = np.array(
                    [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
                     0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1,
                     1.5, 2])
            y['EVA_HOT'] = np.array(
                    [0.030, 0.158, 0.245, 0.313, 0.373, 0.427, 0.477, 0.524,
                     0.569, 0.611, 0.652, 0.692, 0.730, 0.767, 0.803, 0.838,
                     0.872, 0.905, 0.937, 0.969, 1.000, 1.281, 1.523])

            x['EVA_COLD'] = np.array(
                    [0.01, 0.04, 0.07, 0.11, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
                     0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 1, 2])
            y['EVA_COLD'] = np.array(
                    [0.018, 0.075, 0.134, 0.215, 0.300, 0.412, 0.531, 0.658,
                     0.794, 0.934, 0.988, 0.991, 0.994, 0.995, 0.997, 0.998,
                     0.999, 1.000, 1.001])

            x['HE_HOT'] = np.array([
                    0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
                    0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1,
                    1.5, 2])
            y['HE_HOT'] = np.array(
                    [0.030, 0.158, 0.344, 0.469, 0.535, 0.590, 0.638, 0.680,
                     0.718, 0.752, 0.783, 0.812, 0.839, 0.864, 0.887, 0.909,
                     0.929, 0.948, 0.966, 0.984, 1.000, 1.128, 1.216])

            x['HE_COLD'] = np.array([
                    0.01, 0.04, 0.07, 0.11, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
                    0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9,
                    0.95, 1, 1.5, 2])
            y['HE_COLD'] = np.array(
                    [0.018, 0.075, 0.134, 0.215, 0.300, 0.412, 0.507, 0.564,
                     0.614, 0.660, 0.701, 0.739, 0.774, 0.806, 0.836, 0.864,
                     0.890, 0.915, 0.938, 0.960, 0.981, 1.000, 1.151, 1.253])

            x['COND_HOT'] = np.array(
                    [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
                     0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1,
                     1.5, 2])
            y['COND_HOT'] = np.array(
                    [0.030, 0.158, 0.344, 0.567, 0.838, 0.888, 0.906, 0.921,
                     0.933, 0.943, 0.952, 0.959, 0.966, 0.972, 0.977, 0.982,
                     0.986, 0.990, 0.994, 0.997, 1.000, 1.021, 1.033])

            x['COND_COLD'] = np.array([
                    0.01, 0.04, 0.07, 0.11, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
                    0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9,
                    0.95, 1, 1.5, 2])
            y['COND_COLD'] = np.array(
                    [0.019, 0.075, 0.134, 0.192, 0.243, 0.303, 0.359, 0.412,
                     0.463, 0.512, 0.559, 0.604, 0.648, 0.691, 0.733, 0.774,
                     0.813, 0.852, 0.890, 0.928, 0.964, 1.000, 1.327, 1.612])

        else:
            return np.array([0, 1, 2, 3]), np.array([1, 1, 1, 1])

        return x[key], y[key]

    def attr(self):
        return ['x', 'y', 'method', 'comp']

    def f_x(self, x):
        r"""
        Returns characteristic line evaluation at x.

        Parameters
        ----------
        x : float
            Input value for lookup table.

        Returns
        -------
        y : float
            Evaluation of characteristic line at x.

        Note
        ----
        This methods checks for the value range first. If the x-value is outside of the specified range,
        the function will return the values at the corresponding boundary.
        """
        xpos = np.searchsorted(self.x, x)
        if xpos == len(self.x):
            y = self.y[xpos - 1]
        elif xpos == 0:
            y = self.y[0]
        else:
            yfrac = (x - self.x[xpos - 1]) / (self.x[xpos] - self.x[xpos - 1])
            y = self.y[xpos - 1] + yfrac * (self.y[xpos] - self.y[xpos - 1])
        return y

    def get_bound_errors(self, x):
        r"""
        Returns error messages, if operation is out of bounds of characteristc line.

        Parameters
        ----------
        x : float
            Input value for lookup table.

        Returns
        -------
        msg : str
            Error message.
        """
        if x > self.x[-1]:
            msg = ('Operating point above characteristic line range: '
                   'X=' + str(round(x, 3)) + ' with maximum of ' + str(self.x[-1]))
            logging.warning(msg)
        elif x < self.x[0]:
            msg = ('Operating point below characteristic line range: '
                   'X=' + str(round(x, 3)) + ' with minimum of ' + str(self.x[0]))
            logging.warning(msg)

    def get_attr(self, key):
        r"""
        get the value of a characteristics attribute

        :param key: attribute to return its value
        :type key: str
        :returns:
            - :code:`self.__dict__[key]` if object has attribute key
            - :code:`None` if object has no attribute key
        """
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            msg = 'Characteristics has no attribute \"' + key + '\".'
            logging.error(msg)
            raise KeyError(msg)


# %%


class char_map(characteristics):
    r"""
    Class for characteristic maps.

    Parameters
    ----------
    x : ndarray
        An array for the x-values of the map.

    y : ndarray
        An array of the y-values of the map.

    z1 : ndarray
        An array of the z1-values of the map.

    z2 : ndarray
        An array of the z2-values of the map.

    method : str
        Specify a method to choose from the default characteristic maps. If you specify custom x, y, z1 and z2 values, this parameter will be ignored.

    Note
    ----
    This class generates a lookup table from the given input data x, y, z1 and z2, then performs linear interpolation.
    The output parameters are z1 and z2 to be calculated as functions from x and y. The x, y, z1 and z2 values may be specified by the user.
    There is a default characteristic map for axial compressors (GENERIC), see :func:`tespy.components.characteristics.char_map.default` method.

    If you want to use your own map, see the :func:`tespy.components.characteristics.char_map.default` method for more information.
    """

    def __init__(self, **kwargs):

        for key in kwargs:
            if key not in self.attr():
                msg = ('Invalid keyword ' + key + '. Available keywords for kwargs are: ' + str(self.attr()) + '.')
                logging.error(msg)
                raise KeyError(msg)

        # in case of various default characteristics
        method = kwargs.get('method', 'default')

        self.x = kwargs.get('x', None)
        self.y = kwargs.get('y', None)
        self.z1 = kwargs.get('z1', None)
        self.z2 = kwargs.get('z2', None)
        self.comp = kwargs.get('comp', None)

        if self.x is None:
            self.x = self.default(method)[0]
        if self.y is None:
            self.y = self.default(method)[1]

        if self.z1 is None:
            self.z1 = self.default(method)[2]
        if self.z2 is None:
            self.z2 = self.default(method)[3]

        if isinstance(self.x, list):
            self.x = np.array(self.x)
        if isinstance(self.y, list):
            self.y = np.array(self.y)
        if isinstance(self.z1, list):
            self.z1 = np.array(self.z1)
        if isinstance(self.z2, list):
            self.z2 = np.array(self.z2)

        if self.x.shape[0] != self.y.shape[0]:
            msg = ('The number of x-values determines the number of dimension for the characteristic map. You have provided ' +
                   str(len(self.x)) + 'x-values. Thus, the y-, z1- and z2-arrays must have ' + str(len(self.x)) +' number of dimensions.')
            logging.error(msg)
            raise ValueError(msg)
        elif self.y.shape != self.z1.shape or self.y.shape != self.z2.shape:
            msg = 'Make sure that the number of dimensions and the number of values in the y-, z1- and z2-arrays are identical!'
            logging.error(msg)
            raise ValueError(msg)

        msg = 'Created characteristic map for component of type ' + str(self.comp) + ' with default method ' + method + '.'
        logging.debug(msg)

    def attr(self):
        return ['x', 'y', 'z1', 'z2', 'method', 'comp']

    def default(self, key):
        r"""

        Generic characteristic map for axial compressors.

        .. image:: _images/CMAP_GENERIC_PR.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

        .. image:: _images/CMAP_GENERIC_ETA.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

        **literature**

        compressor map:

        - Marcin Plis, Henryk Rusinowski (2016): Mathematical modeling of an
          axial compressor in a gas turbine cycle. Journal of Power
          Technologies 96 (3), pp. 194-199.

        vigv:

        - GasTurb GmbH (2015): GasTurb 12.

        .. note::

            The x-values represent the speedlines, y-values represent the corrected mass flow for each speedline.
            The z1-values are the pressure ratio to nominal pressure ratio and z2-values are the isentropic efficiency to nominal isentropic efficiency.
            Thus, the number of points in the x-array equals the number of dimensions of the y, z1 and z2-array.
            E. g., if you specify 5 speedlines in the x-array, you will need to have an y-array with 5 dimensions,
            where the points of each dimension represent one of the speedlines. The same logic applies for the z1 and z2 arrays!

            For calculation/equations see :func:`tespy.components.characteristics.char_map.get_pr_eta`.
        """

        if key == 'default':
            x = np.array([0, 1, 2])
            y = np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]])
            z1 = y
            z2 = y
            return x, y, z1, z2

        x = {}
        y = {}
        z1 = {}
        z2 = {}

        if self.comp == 'compressor':

            x['GENERIC'] = np.array([0.810, 0.870, 0.946, 0.971, 1, 1.029, 1.062])
            y['GENERIC'] = np.array([[0.460, 0.481, 0.502, 0.523, 0.543,
                                      0.562, 0.583, 0.598, 0.606, 0.612],
                                     [0.590, 0.605, 0.620, 0.640, 0.660,
                                      0.685, 0.703, 0.710, 0.711, 0.713],
                                     [0.767, 0.805, 0.838, 0.859, 0.87,
                                      0.876, 0.878, 0.878, 0.879, 0.88],
                                     [0.874, 0.908, 0.93, 0.943, 0.953,
                                      0.961, 0.962, 0.963, 0.963, 0.964],
                                     [0.948, 0.974, 0.987, 0.995, 1.0,
                                      1.002, 1.005, 1.005, 1.006, 1.006],
                                     [1.014, 1.017, 1.02, 1.023, 1.026,
                                      1.028, 1.03, 1.032, 1.034, 1.036],
                                     [1.045, 1.047, 1.049, 1.051, 1.052,
                                      1.053, 1.054, 1.054, 1.055, 1.056]])

            z1['GENERIC'] = np.array([[0.502, 0.493, 0.485, 0.467, 0.442,
                                       0.411, 0.378, 0.344, 0.31, 0.276],
                                      [0.65, 0.637, 0.617, 0.589, 0.556,
                                       0.519, 0.482, 0.445, 0.407, 0.37],
                                      [0.931, 0.917, 0.893, 0.859, 0.82,
                                       0.779, 0.738, 0.698, 0.657, 0.616],
                                      [1.05, 1.02, 0.982, 0.939, 0.895,
                                       0.851, 0.806, 0.762, 0.717, 0.672],
                                      [1.195, 1.151, 1.102, 1.052, 1.0,
                                       0.951, 0.9, 0.85, 0.799, 0.748],
                                      [1.34, 1.276, 1.213, 1.149, 1.085,
                                       1.022, 0.958, 0.894, 0.831, 0.767],
                                      [1.441, 1.37, 1.3, 1.229, 1.158,
                                       1.088, 1.017, 0.946, 0.876, 0.805]])

            z2['GENERIC'] = np.array([[0.872, 0.885, 0.898, 0.911, 0.925,
                                       0.94, 0.945, 0.926, 0.903, 0.879],
                                      [0.887, 0.909, 0.93, 0.947, 0.963,
                                       0.971, 0.965, 0.939, 0.913, 0.887],
                                      [0.891, 0.918, 0.946, 0.973, 1.001,
                                       1.014, 1.015, 0.986, 0.955, 0.925],
                                      [0.977, 0.977, 0.981, 0.995, 1.007,
                                       1.002, 0.981, 0.961, 0.94, 0.92],
                                      [0.956, 0.959, 0.969, 0.984, 1.0,
                                       0.985, 0.967, 0.95, 0.932, 0.914],
                                      [0.948, 0.959, 0.962, 0.949, 0.935,
                                       0.922, 0.908, 0.895, 0.881, 0.868],
                                      [0.879, 0.888, 0.898, 0.907, 0.916,
                                       0.924, 0.915, 0.906, 0.896, 0.887]])
        else:
            x = np.array([0, 1, 2])
            y = np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]])
            z1 = y
            z2 = y
            return x, y, z1, z2

        return x[key], y[key], z1[key], z2[key]

    def get_pr_eta(self, x, y, igva):
        r"""
        Calculates pressure ratio and isentropic efficiency at given speedline
        and corrected mass flow.

        Parameters
        ----------
        x : float
            Speedline.

        y : float
            Corrected mass flow.

        igva : float
            Inlet guide vane angle.

        Returns
        -------
        pr : float
            Pressure ratio to nominal pressure ratio (Z1).

            .. math::

                Z1 = \frac{p_2 \cdot p_\mathrm{1,ref}}{p_1 \cdot p_\mathrm{2,ref}}=
                f\left(X, Y \right)

        eta : float
            Isentropic efficiency to nominal isentropic efficiency ratio (Z2).

            .. math::

                Z2 = \frac{\eta_\mathrm{s,c}}{\eta_\mathrm{s,c,ref}}=
                f\left(X, Y \right)

        Note
        ----
        .. math::

            X = \sqrt{\frac{T_\mathrm{1,ref}}{T_\mathrm{1}}}

            Y = \frac{\dot{m}_\mathrm{1} \cdot p_\mathrm{1,ref}}
            {\dot{m}_\mathrm{1,ref} \cdot p_\mathrm{1} \cdot X}
        """
        xpos = np.searchsorted(self.x, x)
        if xpos == len(self.x):
            yarr = self.y[xpos - 1]
            z1 = self.z1[xpos - 1]
            z2 = self.z2[xpos - 1]
        elif xpos == 0:
            yarr = self.y[0]
            z1 = self.z1[0]
            z2 = self.z2[0]
        else:
            yfrac = (x - self.x[xpos - 1]) / (self.x[xpos] - self.x[xpos - 1])
            yarr = self.y[xpos - 1] + yfrac * (self.y[xpos] - self.y[xpos - 1])
            z1 = self.z1[xpos - 1] + yfrac * (
                    self.z1[xpos] - self.z1[xpos - 1])
            z2 = self.z2[xpos - 1] + yfrac * (
                    self.z2[xpos] - self.z2[xpos - 1])

        yarr *= (1 - igva / 100)
        z1 *= (1 - igva / 100)
        z2 *= (1 - igva ** 2 / 10000)

        ypos = np.searchsorted(yarr, y)
        if ypos == len(yarr):
            return z1[ypos - 1], z2[ypos - 1]
        elif ypos == 0:
            return z1[0], z2[0]
        else:
            zfrac = (y - yarr[ypos - 1]) / (yarr[ypos] - yarr[ypos - 1])
            pr = z1[ypos - 1] + zfrac * (z1[ypos] - z1[ypos - 1])
            eta = z2[ypos - 1] + zfrac * (z2[ypos] - z2[ypos - 1])
            return pr, eta

    def get_bound_errors(self, x, y, igva):
        r"""
        Returns error messages, if operation is out of bounds of compressor map.

        Parameters
        ----------
        x : float
            Speedline.

        y : float
            Corrected mass flow.

        igva : float
            Inlet guide vane angle.

        Returns
        -------
        msg : str
            Error message.
        """
        xpos = np.searchsorted(self.x, x)
        if xpos == len(self.x) and x != self.x[-1]:
            yarr = self.y[xpos - 1]
            msg = ('Operating point above compressor map range: '
                   'X=' + str(round(x, 3)) + ' with maximum of ' + str(self.x[-1]) + '.')
            logging.warning(msg)
        elif xpos == 0 and y != self.x[0]:
            yarr = self.y[0]
            msg = ('Operating point below compressor map range: '
                   'X=' + str(round(x, 3)) + ' with minimum of ' + str(self.x[0]) + '.')
            logging.warning(msg)
        else:
            yfrac = (x - self.x[xpos - 1]) / (self.x[xpos] - self.x[xpos - 1])
            yarr = self.y[xpos - 1] + yfrac * (self.y[xpos] - self.y[xpos - 1])

        yarr *= (1 - igva / 100)

        ypos = np.searchsorted(yarr, y)
        if ypos == len(yarr) and y != yarr[-1]:
            msg = ('Operating point above compressor map range: '
                   'Y=' + str(round(y, 3)) + ' with maximum of ' + str(yarr[-1]) + '.')
            logging.warning(msg)
            return msg
        elif ypos == 0 and y != yarr[0]:
            msg = ('Operating point below compressor map range: '
                   'Y=' + str(round(y, 3)) + ' with minimum of ' + str(yarr[0]) + '.')
            logging.warning(msg)
