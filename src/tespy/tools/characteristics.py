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

import json
import logging
import os

import numpy as np
from pkg_resources import resource_filename

from tespy.tools.helpers import extend_basic_path

# %%


class char_line:
    r"""
    Class for characteristc lines.

    Parameters
    ----------
    x : ndarray
        An array for the x-values of the lookup table. Number of x and y
        values must be identical.

    y : ndarray
        The corresponding y-values for the lookup table. Number of x and y
        values must be identical.

    extrapolate : boolean
        If :code:`True` linear extrapolation is performed when the x value is
        out of the defined value range.

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

    def __init__(
            self, x=np.array([0, 1]), y=np.array([1, 1]), extrapolate=False):

        self.x = x
        self.y = y
        self.extrapolate = extrapolate

        if isinstance(self.x, list):
            self.x = np.array(self.x)
        if isinstance(self.y, list):
            self.y = np.array(self.y)

        self.x = self.x.astype(float)
        self.y = self.y.astype(float)

        if len(self.x) != len(self.y):
            msg = ('Please provide the same amount of x-values and y-values. '
                   'Number of x-values is ' + str(len(self.x)) + ', number of '
                   'y-values is ' + str(len(self.y)) + ' for char_line.')
            logging.error(msg)
            raise ValueError(msg)

        msg = ('Created characteristic line function.')
        logging.debug(msg)

    def evaluate(self, x):
        r"""
        Return characteristic line evaluation at x.

        Parameters
        ----------
        x : float
            Input value for linear interpolation.

        Returns
        -------
        y : float
            Evaluation of characteristic line at x.

        Note
        ----
        This methods checks for the value range first. If :code:`extrapolate`
        is :code:`False` (default) and the x-value is outside of the specified
        range, the function will return the values at the corresponding
        boundary. If :code:`extrapolate` is :code:`True` the y-value is
        calculated by linear extrapolation.

        .. math::

            y = y_0 + \frac{x-x_0}{x_1-x_0} \cdot \left(y_1-y_0 \right)

        where the index :math:`x_0` represents the lower and :math:`x_1` the
        upper adjacent x-value. :math:`y_0` and :math:`y_1` are the
        corresponding y-values. On extrapolation the two smallest or the two
        largest value pairs are used respectively.
        """
        xpos = np.searchsorted(self.x, x)
        if xpos == len(self.x):
            if self.extrapolate is True:
                xpos = -1
            else:
                return self.y[-1]
        elif xpos == 0:
            if self.extrapolate is True:
                xpos = 1
            else:
                return self.y[0]

        yfrac = (x - self.x[xpos - 1]) / (self.x[xpos] - self.x[xpos - 1])
        return self.y[xpos - 1] + yfrac * (self.y[xpos] - self.y[xpos - 1])

    def get_bound_errors(self, x, c):
        r"""
        Prompt error messages, if x value is out of bounds.

        Parameters
        ----------
        x : float
            Input value for linear interpolation.

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

    def get_attr(self, key):
        r"""
        Get the value of an attribute.

        Parameters
        ----------
        key : str
            Object attribute to get value of.

        Returns
        -------
        value : object
            Value of object attribute key.
        """
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            msg = 'Char_map has no attribute \"' + key + '\".'
            logging.error(msg)
            raise KeyError(msg)

# %%


class char_map:
    r"""
    Class for characteristic maps.

    Parameters
    ----------
    x : ndarray
        An array for the first dimension input of the map.

    y : ndarray
        A two-dimensional array of the second dimension input of the map.

    z1 : ndarray
        A two-dimensional array of the first dimension output of the map.

    z2 : ndarray
        A two-dimensional array of the second dimension output of the map.

    Note
    ----
    This class generates a lookup table from the given input data x, y, z1 and
    z2, then performs linear interpolation. The output parameters are z1 and z2
    to be calculated as functions from x and y. The x, y, z1 and z2 values may
    be specified by the user.
    """

    def __init__(self, x=np.array([0, 1]), y=np.array([[1, 1], [1, 1]]),
                 z1=np.array([[1, 1], [1, 1]]), z2=np.array([[1, 1], [1, 1]])):

        self.x = x
        self.y = y
        self.z1 = z1
        self.z2 = z2

        if isinstance(self.x, list):
            self.x = np.array(self.x)
        if isinstance(self.y, list):
            self.y = np.array(self.y)
        if isinstance(self.z1, list):
            self.z1 = np.array(self.z1)
        if isinstance(self.z2, list):
            self.z2 = np.array(self.z2)

        self.x = self.x.astype(float)
        self.y = self.y.astype(float)
        self.z1 = self.z1.astype(float)
        self.z2 = self.z2.astype(float)

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

        msg = ('Created characteristic map function.')
        logging.debug(msg)

    def evaluate_x(self, x):
        r"""
        Evaluate char_map for x inputs.

        Parameters
        ----------
        x : float
            Input for first dimension of char_map.

        Returns
        -------
        yarr : ndarray
            Second dimension input array of char_map calculated from first
            dimension input.

        z1arr : ndarray
            First dimension output array of char_map calculated from first
            dimension input.

        z2arr : ndarray
            Second dimension output array of char_map calculated from first
            dimension input.
        """
        xpos = np.searchsorted(self.x, x)
        if xpos == len(self.x):
            yarr = self.y[xpos - 1]
            z1arr = self.z1[xpos - 1]
            z2arr = self.z2[xpos - 1]
        elif xpos == 0:
            yarr = self.y[0]
            z1arr = self.z1[0]
            z2arr = self.z2[0]
        else:
            yfrac = (x - self.x[xpos - 1]) / (self.x[xpos] - self.x[xpos - 1])
            yarr = self.y[xpos - 1] + yfrac * (self.y[xpos] - self.y[xpos - 1])
            z1arr = self.z1[xpos - 1] + yfrac * (self.z1[xpos] -
                                                 self.z1[xpos - 1])
            z2arr = self.z2[xpos - 1] + yfrac * (self.z2[xpos] -
                                                 self.z2[xpos - 1])

        return yarr, z1arr, z2arr

    def evaluate_y(self, y, yarr, z1arr, z2arr):
        r"""
        Evaluate char_map for y inputs.

        Parameters
        ----------
        y : float
            Input for second dimension of char_map.

        yarr : ndarray
            Second dimension array of char_map calculated from first dimension
            input.

        z1arr : ndarray
            First dimension output array of char_map calculated from first
            dimension input.

        z2arr : ndarray
            Second dimension output array of char_map calculated from first
            dimension input.
        """
        ypos = np.searchsorted(yarr, y)
        if ypos == len(yarr):
            return z1arr[ypos - 1], z2arr[ypos - 1]
        elif ypos == 0:
            return z1arr[0], z2arr[0]
        else:
            zfrac = (y - yarr[ypos - 1]) / (yarr[ypos] - yarr[ypos - 1])
            z1 = z1arr[ypos - 1] + zfrac * (z1arr[ypos] - z1arr[ypos - 1])
            z2 = z2arr[ypos - 1] + zfrac * (z2arr[ypos] - z2arr[ypos - 1])
            return z1, z2

    def evaluate(self, x, y, **kwargs):
        r"""
        Evaluate char_map for x and y inputs.

        Parameters
        ----------
        x : float
            Input for first dimension of char_map.

        y : float
            Input for second dimension of char_map.

        Returns
        -------
        z1 : float
            Resulting z1 value.

        z2 : float
            Resulting z2 value.

        Note
        ----
        This methods checks for the value range first. If the x-value is
        outside of the specified range, the function will return the arrays
        for :math:`y`, :math:`z1` and :math:`z2`.

        .. math::

            \vec{y} = \vec{y_0} + \frac{x-x_0}{x_1-x_0} \cdot
            \left(\vec{y_1}-\vec{y_0} \right)\\
            \vec{z1} = \vec{z1_0} + \frac{x-x_0}{x_1-x_0} \cdot
            \left(\vec{z1_1}-\vec{z1_0} \right)\\
            \vec{z2} = \vec{z2_0} + \frac{x-x_0}{x_1-x_0} \cdot
            \left(\vec{z2_1}-\vec{z2_0}\right)

        The index :math:`x_0` represents the lower and :math:`x_1` the
        upper adjacent x-value. Using the y-value as second input dimension
        the corresponding z1- and z2-values are calculated, again using linear
        interpolation.

        .. math::

            z1 = z1_0 + \frac{y-y_0}{y_1-y_0} \cdot \left(z1_1-z1_0 \right)\\
            z2 = z2_0 + \frac{y-y_0}{y_1-y_0} \cdot \left(z2_1-z2_0 \right)
        """
        yarr, z1arr, z2arr = self.evaluate_x(x)
        z1, z2 = self.evaluate_y(y, yarr, z1arr, z2arr)

        return z1, z2

    def get_bound_errors_x(self, x, c):
        r"""
        Prompt error message, if operation is out bounds in first dimension.

        Parameters
        ----------
        x : float
            Input for first dimension of char_map.

        c : str
            Label of the component, the char_map is applied on.

        Returns
        -------
        yarr : ndarray
            Second dimension input array of char_map calculated from first
            dimension input.
        """
        xpos = np.searchsorted(self.x, x)
        if xpos == len(self.x) and x != self.x[-1]:
            yarr = self.y[xpos - 1]
            msg = ('Operating point above compressor map range: '
                   'X=' + str(round(x, 3)) + ' with maximum of ' +
                   str(self.x[-1]) + ' at component ' + c + '.')
            logging.warning(msg)
        elif xpos == 0 and x != self.x[0]:
            yarr = self.y[0]
            msg = ('Operating point below compressor map range: '
                   'X=' + str(round(x, 3)) + ' with minimum of ' +
                   str(self.x[0]) + ' at component ' + c + '.')
            logging.warning(msg)
        else:
            yfrac = (x - self.x[xpos - 1]) / (self.x[xpos] - self.x[xpos - 1])
            yarr = self.y[xpos - 1] + yfrac * (self.y[xpos] - self.y[xpos - 1])

        return yarr

    def get_bound_errors_y(self, y, yarr, c):
        r"""
        Prompt error message, if operation is out bounds in second dimension.

        Parameters
        ----------
        y : float
            Input for second dimension of char_map.

        yarr : ndarray
            Second dimension input array of char_map calculated from first
            dimension input.

        c : str
            Label of the component, the char_map is applied on.
        """
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

    def get_bound_errors(self, x, y, c):
        r"""
        Check the char_map for bound violations.

        Parameters
        ----------
        x : float
            Input for first dimension of char_map.

        y : float
            Input for second dimension of char_map.
        """
        yarr = self.get_bound_errors_x(x, c)
        self.get_bound_errors_y(y, yarr, c)

    def get_attr(self, key):
        r"""
        Get the value of an attribute.

        Parameters
        ----------
        key : str
            Object attribute to get value of.

        Returns
        -------
        value : object
            Value of object attribute key.
        """
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            msg = 'Char_map has no attribute \"' + key + '\".'
            logging.error(msg)
            raise KeyError(msg)


class compressor_map(char_map):
    r"""
    Class for compressor maps.

    Parameters
    ----------
    x : ndarray
        An array for the first dimension input of the map.

    y : ndarray
        A two-dimensional array of the second dimension input of the map.

    z1 : ndarray
        A two-dimensional array of the first dimension output of the map.

    z2 : ndarray
        A two-dimensional array of the second dimension output of the map.

    Note
    ----
    Ìn contrast to the :py:class:`tespy.tools.characteristics.char_map` this
    classes output is manipulated with the inlet guide vane angle
    :code:`igva`.
    """

    def evaluate(self, x, y, **kwargs):
        r"""
        Evaluate compressor_map for x and y inputs.

        This method is different from the base method. The second dimension
        array is manipulated by the inlet guide vane angle igva.

        Parameters
        ----------
        x : float
            Input for first dimension of char_map.

        y : float
            Input for second dimension of char_map.

        igva : float
            Inlet guide vane angle of the compressor.

        Returns
        -------
        z1 : float
            Pressure ratio of compressor.

        z2 : float
            Isentropic efficiency of compressor.

        Note
        ----
        In contrast to the :py:class:`tespy.tools.characteristics.char_map`
        the values are manipulated by the inlet guide vane angle (igva):

        .. math::

            \vec{y} = \vec{y} \cdot \left( 1 - \frac{igva}{100} \right)\\
            \vec{z1} = \vec{z1} \cdot \left( 1 - \frac{igva}{100} \right)\\
            \vec{z2} = \vec{z2} \cdot \left(
            1 - \left(\frac{igva}{100}\right)^2 \right)

        Reference: :cite:`GasTurb2018`.
        """
        yarr, z1arr, z2arr = self.evaluate_x(x)
        igva = kwargs.get('igva', 0)

        yarr *= (1 - igva / 100)
        z1arr *= (1 - igva / 100)
        z2arr *= (1 - igva ** 2 / 10000)

        z1, z2 = self.evaluate_y(y, yarr, z1arr, z2arr)

        return z1, z2

    def get_bound_errors(self, x, y, igva, c):
        r"""
        Check the compressor_map for bound violations.

        This method is different from the base method. The second dimension
        array is manipulated by the inlet guide vane angle igva.

        Parameters
        ----------
        x : float
            Input for first dimension of char_map.

        y : float
            Input for second dimension of char_map.

        igva : float
            Inlet guide vane angle of the compressor.
        """
        yarr = self.get_bound_errors_x(x, c)
        yarr *= (1 - igva / 100)
        self.get_bound_errors_y(y, yarr, c)


def load_default_char(component, parameter, function_name, char_type):
    r"""
    Load a characteristic line of map.

    Parameters
    ----------
    component : str
        Type of component.

    parameter : str
        Component parameter using the characteristics.

    function_name : str
        Name of the characteristics.

    char_type : class
        Class to generate an instance of.

    Returns
    -------
    obj : object
        The characteristics (char_line, char_map, compressor_map) object.
    """
    if char_type == char_line:
        path = resource_filename('tespy.data', 'char_lines.json')
    else:
        path = resource_filename('tespy.data', 'char_maps.json')

    with open(path) as f:
        data = json.loads(f.read())

    if char_type == char_line:
        x = data[component][parameter][function_name]['x']
        y = data[component][parameter][function_name]['y']
        obj = char_type(x, y)

    else:
        x = data[component][parameter][function_name]['x']
        y = data[component][parameter][function_name]['y']
        z1 = data[component][parameter][function_name]['z1']
        z2 = data[component][parameter][function_name]['z2']
        obj = char_type(x, y, z1, z2)

    return obj


def load_custom_char(name, char_type):
    r"""
    Load a characteristic line of map.

    Parameters
    ----------
    name : str
        Name of the characteristics.

    char_type : class
        Class to be generate the object of.

    Returns
    -------
    obj : object
        The characteristics (char_line, char_map, compressor_map) object.
    """
    path = extend_basic_path('data')

    if char_type == char_line:
        path = os.path.join(path, 'char_lines.json')
    else:
        path = os.path.join(path, 'char_maps.json')

    if os.path.isfile(path):

        with open(path) as f:
            data = json.loads(f.read())

        if char_type == char_line:
            x = data[name]['x']
            y = data[name]['y']
            obj = char_type(x, y)

        else:
            x = data[name]['x']
            y = data[name]['y']
            z1 = data[name]['z1']
            z2 = data[name]['z2']
            obj = char_type(x, y, z1, z2)

        return obj

    else:
        msg = ('The file containing your custom charactersitics could not be '
               'found on your system. The path should be ' + path + '. Please '
               'make sure the .tespy/data path exists in your home directory.')
        logging.error(msg)
        raise FileNotFoundError(msg)
