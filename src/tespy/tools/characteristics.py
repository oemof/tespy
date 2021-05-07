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
from matplotlib import pyplot as plt
from pkg_resources import resource_filename

from tespy.tools.helpers import extend_basic_path

# %%


class CharLine:
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
            self, x=np.array([0, 1]), y=np.ones((2)), extrapolate=False):

        self.x = x
        self.y = y
        self.extrapolate = extrapolate

        if isinstance(self.x, list):
            self.x = np.asarray(self.x)
        if isinstance(self.y, list):
            self.y = np.asarray(self.y)

        self.x = self.x.astype(float)
        self.y = self.y.astype(float)

        if len(self.x) != len(self.y):
            msg = ('Please provide the same amount of x-values and y-values. '
                   'Number of x-values is ' + str(len(self.x)) + ', number of '
                   'y-values is ' + str(len(self.y)) + ' for CharLine.')
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
            if self.extrapolate:
                xpos = -1
            else:
                return self.y[-1]
        elif xpos == 0:
            if self.extrapolate:
                xpos = 1
            else:
                return self.y[0]

        yfrac = (x - self.x[xpos - 1]) / (self.x[xpos] - self.x[xpos - 1])
        return self.y[xpos - 1] + yfrac * (self.y[xpos] - self.y[xpos - 1])

    def get_domain_errors(self, x, c):
        r"""
        Prompt error messages, if x value is out of bounds.

        Parameters
        ----------
        x : float
            Input value for linear interpolation.
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

    def plot(self, path, title, xlabel, ylabel):

        # plotting
        fig = plt.figure()
        ax = plt.subplot()
        ax.plot(self.x, self.y, 'x', mew=2)
        plt.grid(linestyle='dotted')
        # formatting
        plt.tight_layout()
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        # export
        fig.savefig(path, bbox_inches='tight')
        plt.close(fig)

# %%


class CharMap:
    r"""
    Class for characteristic maps.

    Parameters
    ----------
    x : ndarray
        An array for the first dimension input of the map.

    y : ndarray
        A two-dimensional array of the second dimension input of the map.

    z : ndarray
        A two-dimensional array of the output of the map.

    Note
    ----
    This class generates a lookup table from the given input data x, y and z,
    then performs linear interpolation. The output parameter is z to be
    calculated as functions from x and y.
    """

    def __init__(self, x=np.array([0, 1]), y=np.ones((2, 2)),
                 z=np.ones((2, 2))):

        self.x = x
        self.y = y
        self.z = z

        if isinstance(self.x, list):
            self.x = np.array(self.x)
        if isinstance(self.y, list):
            self.y = np.array(self.y, dtype=object)
        if isinstance(self.z, list):
            self.z = np.array(self.z, dtype=object)

        self.x = self.x.astype(float)
        self.y = self.y.astype(float)
        self.z = self.z.astype(float)

        if self.x.shape[0] != self.y.shape[0]:
            msg = (
                'The number of x-values determines the number of dimension '
                'for the characteristic map. You have provided ' +
                str(len(self.x)) + 'x-values. Thus, the y- and z-arrays must '
                'have ' + str(len(self.x)) + ' number of dimensions.')
            logging.error(msg)
            raise ValueError(msg)
        elif self.y.shape != self.z.shape:
            msg = (
                'Make sure that the number of dimensions and the number of '
                'values in the y-, z-arrays are identical!')
            logging.error(msg)
            raise ValueError(msg)

        msg = ('Created characteristic map function.')
        logging.debug(msg)

    def evaluate_x(self, x):
        r"""
        Evaluate CharMap for x inputs.

        Parameters
        ----------
        x : float
            Input for first dimension of CharMap.

        Returns
        -------
        yarr : ndarray
            Second dimension input array of CharMap calculated from first
            dimension input.

        zarr : ndarray
            Output array of CharMap calculated from first dimension input.
        """
        xpos = np.searchsorted(self.x, x)
        if xpos == len(self.x):
            yarr = self.y[xpos - 1]
            zarr = self.z[xpos - 1]
        elif xpos == 0:
            yarr = self.y[0]
            zarr = self.z[0]
        else:
            yfrac = (x - self.x[xpos - 1]) / (self.x[xpos] - self.x[xpos - 1])
            yarr = self.y[xpos - 1] + yfrac * (self.y[xpos] - self.y[xpos - 1])
            zarr = self.z[xpos - 1] + yfrac * (self.z[xpos] - self.z[xpos - 1])

        return yarr, zarr

    def evaluate_y(self, y, yarr, zarr):
        r"""
        Evaluate CharMap for y inputs.

        Parameters
        ----------
        y : float
            Input for second dimension of CharMap.

        yarr : ndarray
            Second dimension array of CharMap calculated from first dimension
            input.

        zarr : ndarray
            Output array of CharMap calculated from first dimension input.
        """
        ypos = np.searchsorted(yarr, y)
        if ypos == len(yarr):
            return zarr[ypos - 1]
        elif ypos == 0:
            return zarr[0]
        else:
            zfrac = (y - yarr[ypos - 1]) / (yarr[ypos] - yarr[ypos - 1])
            z = zarr[ypos - 1] + zfrac * (zarr[ypos] - zarr[ypos - 1])
            return z

    def evaluate(self, x, y):
        r"""
        Evaluate CharMap for x and y inputs.

        Parameters
        ----------
        x : float
            Input for first dimension of CharMap.

        y : float
            Input for second dimension of CharMap.

        Returns
        -------

        z : float
            Resulting z value.

        Note
        ----
        .. math::

            \vec{y} = \vec{y_0} + \frac{x-x_0}{x_1-x_0} \cdot
            \left(\vec{y_1}-\vec{y_0} \right)\\
            \vec{z} = \vec{z1_0} + \frac{x-x_0}{x_1-x_0} \cdot
            \left(\vec{z_1}-\vec{z_0} \right)

        The index :code:`0` represents the lower and :code:`1` the
        upper adjacent x-value. Using the y-value as second input dimension
        the corresponding z-values are calculated, again using linear
        interpolation.

        .. math::

            z = z_0 + \frac{y-y_0}{y_1-y_0} \cdot \left(z_1-z_0 \right)
        """
        return self.evaluate_y(y, *self.evaluate_x(x))

    def get_domain_errors_x(self, x, c):
        r"""
        Prompt error message, if operation is out bounds in first dimension.

        Parameters
        ----------
        x : float
            Input for first dimension of CharMap.

        c : str
            Label of the component, the CharMap is applied on.

        Returns
        -------
        yarr : ndarray
            Second dimension input array of CharMap calculated from first
            dimension input.
        """
        xpos = np.searchsorted(self.x, x)
        if xpos == len(self.x) and x != self.x[-1]:
            yarr = self.y[xpos - 1]
            msg = ('Operating point above CharMap range: '
                   'X=' + str(round(x, 3)) + ' with maximum of ' +
                   str(self.x[-1]) + ' at component ' + c + '.')
            logging.warning(msg)
        elif xpos == 0 and x != self.x[0]:
            yarr = self.y[0]
            msg = ('Operating point below CharMap range: '
                   'X=' + str(round(x, 3)) + ' with minimum of ' +
                   str(self.x[0]) + ' at component ' + c + '.')
            logging.warning(msg)
        else:
            yfrac = (x - self.x[xpos - 1]) / (self.x[xpos] - self.x[xpos - 1])
            yarr = self.y[xpos - 1] + yfrac * (self.y[xpos] - self.y[xpos - 1])

        return yarr

    def get_domain_errors_y(self, y, yarr, c):
        r"""
        Prompt error message, if operation is out bounds in second dimension.

        Parameters
        ----------
        y : float
            Input for second dimension of CharMap.

        yarr : ndarray
            Second dimension input array of CharMap calculated from first
            dimension input.

        c : str
            Label of the component, the CharMap is applied on.
        """
        ypos = np.searchsorted(yarr, y)
        if ypos == len(yarr) and y != yarr[-1]:
            msg = (
                'Operating point above compressor map range: Y=' +
                str(round(y, 3)) + ' with maximum of ' + str(yarr[-1]) +
                ' at component ' + c + '.')
            logging.warning(msg)
        elif ypos == 0 and y != yarr[0]:
            msg = (
                'Operating point below compressor map range: Y=' +
                str(round(y, 3)) + ' with minimum of ' + str(yarr[0]) +
                ' at component ' + c + '.')
            logging.warning(msg)

    def get_domain_errors(self, x, y, c):
        r"""
        Check the CharMap for bound violations.

        Parameters
        ----------
        x : float
            Input for first dimension of CharMap.

        y : float
            Input for second dimension of CharMap.
        """
        yarr = self.get_domain_errors_x(x, c)
        self.get_domain_errors_y(y, yarr, c)

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

    def plot(self, path, title, xlabel, ylabel):

        # plotting
        fig = plt.figure()
        ax = plt.subplot()
        for datapoint in range(len(self.x)):
            ax.plot(self.y[datapoint], self.z[datapoint], 'x', mew=2)
        plt.grid(linestyle='dotted')
        # formatting
        plt.tight_layout()
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        # export
        fig.savefig(path, bbox_inches='tight')
        plt.close(fig)


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
        The characteristics (CharLine, CharMap) object.
    """
    if char_type == CharLine:
        path = resource_filename('tespy.data', 'char_lines.json')
    else:
        path = resource_filename('tespy.data', 'char_maps.json')

    with open(path) as f:
        data = json.loads(f.read())

    if char_type == CharLine:
        x = data[component][parameter][function_name]['x']
        y = data[component][parameter][function_name]['y']
        obj = CharLine(x, y)

    else:
        x = data[component][parameter][function_name]['x']
        y = data[component][parameter][function_name]['y']
        z = data[component][parameter][function_name]['z']
        obj = CharMap(x, y, z)

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
        The characteristics (CharLine, CharMap) object.
    """
    path = extend_basic_path('data')

    if char_type == CharLine:
        path = os.path.join(path, 'char_lines.json')
    else:
        path = os.path.join(path, 'char_maps.json')

    if os.path.isfile(path):

        with open(path) as f:
            data = json.loads(f.read())

        if char_type == CharLine:
            x = data[name]['x']
            y = data[name]['y']
            obj = CharLine(x, y)

        else:
            x = data[name]['x']
            y = data[name]['y']
            z = data[name]['z']
            obj = CharMap(x, y, z)

        return obj

    else:
        msg = ('The file containing your custom charactersitics could not be '
               'found on your system. The path should be ' + path + '. Please '
               'make sure the .tespy/data path exists in your home directory.')
        logging.error(msg)
        raise FileNotFoundError(msg)
