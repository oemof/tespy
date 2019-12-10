# -*- coding: utf-8

"""Module for characteristic functions.

The characteristics module provides the integration of characteristic lines
and characteristic maps. The user can create custom characteristic lines or
maps with individual data.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/tools/characteristics.py

SPDX-License-Identifier: MIT
"""

import numpy as np

import logging

# %%


class char_line:
    r"""
    Class characteristics for components.

    Parameters
    ----------
    x : ndarray
        An array for the x-values of the lookup table. Number of x and y
        values must be identical.

    y : ndarray
        The corresponding y-values for the lookup table. Number of x and y
        values must be identical.

    Note
    ----
    This class generates a lookup table from the given input data x and y,
    then performs linear interpolation. The x and y values may be specified by
    the user. There are some default characteristic lines for different
    components, see the :py:mod:`tespy.data` module. If you neither specify the
    method to use from the defaults nor specify x and y values, the
    characteristic line generated will be
    :code:`x = [0, 1], y = [1, 1]`.
    """

    def __init__(self, x=np.array([0, 1]), y=np.array([1, 1])):

        self.x = x
        self.y = y

        if isinstance(self.x, list):
            self.x = np.array(self.x)
        if isinstance(self.y, list):
            self.y = np.array(self.y)

        if len(self.x) != len(self.y):
            msg = ('Please provide the same amount of x-values and y-values. '
                   'Number of x-values is ' + str(len(self.x)) + ', number of '
                   'y-values is ' + str(len(self.y)) + ' for char_line.')
            logging.error(msg)
            raise ValueError(msg)

        msg = ('Created characteristic function.')
        logging.debug(msg)

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
        This methods checks for the value range first. If the x-value is
        outside of the specified range, the function will return the values at
        the corresponding boundary.
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

    def get_bound_errors(self, x, c):
        r"""
        Returns error messages, if operation is out of bounds of characteristc
        line.

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
                   'X=' + str(round(x, 3)) + ' with maximum of ' +
                   str(self.x[-1]) + ' at component ' + c + '.')
            logging.warning(msg)
        elif x < self.x[0]:
            msg = ('Operating point below characteristic line range: '
                   'X=' + str(round(x, 3)) + ' with minimum of ' +
                   str(self.x[0]) + ' at component ' + c + '.')
            logging.warning(msg)

# %%


class char_map:
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
        Specify a method to choose from the default characteristic maps. If you
        specify custom x, y, z1 and z2 values, this parameter will be ignored.

    Note
    ----
    This class generates a lookup table from the given input data x, y, z1 and
    z2, then performs linear interpolation. The output parameters are z1 and z2
    to be calculated as functions from x and y. The x, y, z1 and z2 values may
    be specified by the user. There is a default characteristic map for axial
    compressors (GENERIC),
    see :func:`tespy.components.characteristics.char_map.default` method.

    If you want to use your own map, see the
    :func:`tespy.components.characteristics.char_map.default` method for more
    information.

    TODO: The part on default characteristics will be adjusted for version
    0.2.0. Some parts of the docstrings will have to be reworked.
    """

    def __init__(self, **kwargs):

        for key in kwargs:
            if key not in self.attr():
                msg = ('Invalid keyword ' + key +
                       '. Available keywords are: ' + str(self.attr()) + '.')
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
            msg = ('The number of x-values determines the number of dimension '
                   'for the characteristic map. You have provided ' +
                   str(len(self.x)) + 'x-values. Thus, the y-, z1- and '
                   'z2-arrays must have ' + str(len(self.x)) +
                   ' number of dimensions.')
            logging.error(msg)
            raise ValueError(msg)
        elif self.y.shape != self.z1.shape or self.y.shape != self.z2.shape:
            msg = ('Make sure that the number of dimensions and the number of '
                   'values in the y-, z1- and z2-arrays are identical!')
            logging.error(msg)
            raise ValueError(msg)

        msg = ('Created characteristic map for component of type ' +
               str(self.comp) + ' with default method ' + method + '.')
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

            The x-values represent the speedlines, y-values represent the
            corrected mass flow for each speedline. The z1-values are the
            pressure ratio to nominal pressure ratio and z2-values are the
            isentropic efficiency to nominal isentropic efficiency. Thus, the
            number of points in the x-array equals the number of dimensions of
            the y, z1 and z2-array. E. g., if you specify 5 speedlines in the
            x-array, you will need to have an y-array with 5 dimensions, where
            the points of each dimension represent one of the speedlines.
            The same logic applies for the z1 and z2 arrays!

            For calculation/equations see
            :func:`tespy.components.characteristics.char_map.get_pr_eta`.
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

            x['GENERIC'] = np.array([0.810, 0.870, 0.946, 0.971, 1, 1.029,
                                     1.062])

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

                Z1 = \frac{p_2 \cdot p_\mathrm{1,ref}}{p_1 \cdot
                p_\mathrm{2,ref}}= f\left(X, Y \right)

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

    def get_bound_errors(self, x, y, igva, c):
        r"""
        Returns error message, if operation is out of bounds of compressor map.

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
                   'X=' + str(round(x, 3)) + ' with maximum of ' +
                   str(self.x[-1]) + ' at component ' + c + '.')
            logging.warning(msg)
        elif xpos == 0 and y != self.x[0]:
            yarr = self.y[0]
            msg = ('Operating point below compressor map range: '
                   'X=' + str(round(x, 3)) + ' with minimum of ' +
                   str(self.x[0]) + ' at component ' + c + '.')
            logging.warning(msg)
        else:
            yfrac = (x - self.x[xpos - 1]) / (self.x[xpos] - self.x[xpos - 1])
            yarr = self.y[xpos - 1] + yfrac * (self.y[xpos] - self.y[xpos - 1])

        yarr *= (1 - igva / 100)

        ypos = np.searchsorted(yarr, y)
        if ypos == len(yarr) and y != yarr[-1]:
            msg = ('Operating point above compressor map range: '
                   'Y=' + str(round(y, 3)) + ' with maximum of ' +
                   str(yarr[-1]) + ' at component ' + c + '.')
            logging.warning(msg)
            return msg
        elif ypos == 0 and y != yarr[0]:
            msg = ('Operating point below compressor map range: '
                   'Y=' + str(round(y, 3)) + ' with minimum of ' +
                   str(yarr[0]) + ' at component ' + c + '.')
            logging.warning(msg)
