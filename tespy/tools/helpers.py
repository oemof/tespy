# -*- coding: utf-8

"""
.. module:: helpers
    :synopsis: helpers for frequently used functionalities

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>
"""

import CoolProp as CP
from CoolProp.CoolProp import PropsSI as CPPSI

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import math
import numpy as np
import sys
from scipy import interpolate
import pandas as pd
import os
import collections

import logging

import warnings
warnings.simplefilter("ignore", RuntimeWarning)

global err
err = 1e-6
global molar_masses
molar_masses = {}
global gas_constants
gas_constants = {}
gas_constants['uni'] = 8.3144598

# %%


class data_container:
    r"""
    Class data_container is the base class for dc_cc, dc_cp, dc_flu, dc_prop.

    Parameters
    ----------
    **kwargs :
        See the class documentation of desired data_container for available keywords.

    Note
    ----
    The initialisation method (__init__), setter method (set_attr) and getter method (get_attr)
    are used for instances of class data_container and its children.

    Example
    -------
    TESPy uses different data_containers for specific tasks:
    Component characteristics (dc_cc), component maps (dc_cm), component properties (dc_cp),
    grouped component properites (dc_gcp), fluid composition (dc_flu),
    fluid properties (dc_prop). Grouped component properties are used, if more than one
    component propertie has to be specified in order to apply one equation, e. g. pressure drop in pipes by
    specified length, diameter and roughness. If you specify all three of these properties,
    the data_container for the group will be created automatically!

    For the full list of available parameters for each data container, see its documentation.

    >>> from tespy import hlp, cmp
    >>> type(hlp.dc_cm(is_set=True))
    <class 'tespy.tools.helpers.dc_cm'>
    >>> type(hlp.dc_cc(x=[1, 2, 3, 4], y=[1, 4, 9, 16], is_set=True))
    <class 'tespy.tools.helpers.dc_cc'>
    >>> type(hlp.dc_cp(val=100, is_set=True, is_var=True, printout=True,
    ...     max_val=1000, min_val=1))
    <class 'tespy.tools.helpers.dc_cp'>
    >>> pipe = cmp.pipe('testpipe', L=100, D=0.5, ks=5e-5)
    >>> type(hlp.dc_gcp(is_set=True, elements=[pipe.L, pipe.D, pipe.ks],
    ...     method='default'))
    <class 'tespy.tools.helpers.dc_gcp'>
    >>> type(hlp.dc_flu(val={'CO2': 0.1, 'H2O': 0.11, 'N2': 0.75, 'O2': 0.03},
    ...     val_set={'CO2': False, 'H2O': False, 'N2': False, 'O2': True},
    ...     balance=False))
    <class 'tespy.tools.helpers.dc_flu'>
    >>> type(hlp.dc_prop(val=5, val_SI=500000, val_set=True, unit='bar',
    ...     unit_set=False, ref=None, ref_set=False))
    <class 'tespy.tools.helpers.dc_prop'>
    >>> type(hlp.dc_simple(val=5, val_set=False))
    <class 'tespy.tools.helpers.dc_simple'>
    """

    def __init__(self, **kwargs):

        var = self.attr()

        # default values
        for key in var.keys():
            self.__dict__.update({key: var[key]})

        self.set_attr(**kwargs)

    def set_attr(self, **kwargs):
        r"""
        Sets, resets or unsets attributes of a data_container type object.

        Parameters
        ----------
        **kwargs :
            See the class documentation of desired data_container for available keywords.
        """
        var = self.attr()
        # specify values
        for key in kwargs:
            if key in var.keys():
                self.__dict__.update({key: kwargs[key]})

    def get_attr(self, key):
        r"""
        Get the value of a data_container's attribute.

        Parameters
        ----------
        key : str
            The attribute you want to retrieve.

        Returns
        -------
        out :
            Specified attribute.
        """
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            msg = 'Datacontainer of type ' + self.__class__.__name__ + ' has no attribute \"' + str(key) + '\".'
            logging.error(msg)
            raise KeyError(msg)

    def attr(self):
        r"""
        Return the available attributes for a data_container type object.

        Returns
        -------
        out : dict
            Dictionary of available attributes (dictionary keys) with default values.
        """
        return {}


class dc_prop(data_container):
    r"""
    Parameters
    ----------
    val : float
        Value in user specified unit (or network unit) if unit is unspecified,
        default: val=np.nan.

    val0 : float
        Starting value in user specified unit (or network unit) if unit is unspecified,
        default: val0=np.nan.

    val_SI : float
        Value in SI_unit, default: val_SI=0.

    val_set : boolean
        Has the value for this property been set?, default: val_set=False.

    ref : tespy.connections.ref
        Reference object, default: ref=None.

    ref_set : boolean
        Has a value for this property been referenced to another connection?, default: ref_set=False.

    unit : str
        Unit for this property, default: ref=None.

    unit : boolean
        Has the unit for this property been specified manually by the user?, default: unit_set=False.
    """
    def attr(self):
        return {'val': np.nan, 'val0': np.nan, 'val_SI': 0, 'val_set': False,
                'ref': None, 'ref_set': False,
                'unit': None, 'unit_set': False, 'design': np.nan}


class dc_simple(data_container):
    r"""
    Simple data container without data type restrictions to val field.

    Parameters
    ----------
    val : no specific datatype
        Value for the property, no predefined datatype. Unset this property by
        stating val=np.nan.

    val_set : boolean
        Has the value for this property been set? default: val_set=False.
    """
    def attr(self):
        return {'val': np.nan, 'val_set': False}


class dc_flu(data_container):
    r"""
    Parameters
    ----------

    val : dict
        Mass fractions of the fluids in a mixture, default: val={}.
        Pattern for dictionary: keys are fluid name, values are mass fractions.

    val0 : dict
        Starting values for mass fractions of the fluids in a mixture, default: val0={}.
        Pattern for dictionary: keys are fluid name, values are mass fractions.

    val_set : dict
        Which fluid mass fractions have been set, default val_set={}.
        Pattern for dictionary: keys are fluid name, values are True or False.

    balance : boolean
        Should the fluid balance equation be applied for this mixture? default: False.
    """
    def attr(self):
        r"""
        Return the available attributes for a data_container type object.

        Returns
        -------
        out : dict
            Dictionary of available attributes (dictionary keys) with default values.
        """
        return {'val': {}, 'val0': {}, 'val_set': {}, 'design': collections.OrderedDict(), 'balance': False}


class dc_cp(data_container):
    r"""
    Parameters
    ----------

    val : float
        Value for this component attribute, default: val=1.

    val_SI : float
        Value in SI_unit (available for temperatures only, unit transformation according to network's temperature unit),
        default: val_SI=0.

    is_set : boolean
        Has the value for this attribute been set?, default: is_set=False.

    is_var : boolean
        Is this attribute part of the system variables?, default: is_var=False.

    d : float
        Interval width for numerical calculation of partial derivative towards this attribute,
        it is part of the system variables, default d=1e-4.

    min_val : float
        Minimum value for this attribute, used if attribute is part of the system variables,
        default: min_val=1.1e-4.

    max_val : float
        Maximum value for this attribute, used if attribute is part of the system variables,
        default: max_val=1e12.

    printout : boolean
        Should the value of this attribute be printed in the results overview?
    """
    def attr(self):
        r"""
        Return the available attributes for a data_container type object.

        Returns
        -------
        out : dict
            Dictionary of available attributes (dictionary keys) with default values.
        """
        return {'val': 1, 'val_SI': 0, 'is_set': False, 'printout': True,
                'd': 1e-4, 'min_val': 0, 'max_val': 1e12, 'is_var': False,
                'val_ref': 1, 'design': np.nan}


class dc_cc(data_container):
    r"""
    Parameters
    ----------

    func : tespy.components.characteristics.characteristics
        Function to be applied for this characteristics, default: None.

    is_set : boolean
        Should this equation be applied?, default: is_set=False.

    method : str
        Which default method for this characteristic function should be used?, default: method='default'.

    param : str
        Which parameter should be applied as the x value?, default: method='default'.

    x : numpy.array
        Array for the x-values of the characteristic line, default x=None.

    y : numpy.array
        Array for the y-values of the characteristic line, default y=None.

    Note
    ----
    If you do not specify x-values or y-values, default values according to the
    specified method will be used. If you specify a method as well as x-values and/or
    y-values, these will override the defaults values of the chosen method.
    """
    def attr(self):
        r"""
        Return the available attributes for a data_container type object.

        Returns
        -------
        out : dict
            Dictionary of available attributes (dictionary keys) with default values.
        """
        return {'func': None, 'is_set': False,
                'method': 'default', 'param': None,
                'x': None, 'y': None}


class dc_cm(data_container):
    r"""
    Parameters
    ----------

    func : tespy.components.characteristics.characteristics
        Function to be applied for this characteristic map, default: None.

    is_set : boolean
        Should this equation be applied?, default: is_set=False.

    method : str
        Which default method for this characteristic function should be used?, default: method='default'.

    param : str
        Which parameter should be applied as the x value?, default: method='default'.

    x : numpy.array
        Array for the x-values of the characteristic line, default x=None.

    y : numpy.array
        Array for the y-values of the characteristic line, default y=None.

    z1 : numpy.array
        Array for the y-values of the characteristic line, default y=None.

    z2 : numpy.array
        Array for the y-values of the characteristic line, default y=None.

    Note
    ----
    If you do not specify any interpolation points (x, y, z1, z2), default values according to the
    specified method will be used. If you specify a method as well as interpolation points,
    these will override the defaults values of the chosen method.
    """
    def attr(self):
        r"""
        Return the available attributes for a data_container type object.

        Returns
        -------
        out : dict
            Dictionary of available attributes (dictionary keys) with default values.
        """
        return {'func': None, 'is_set': False,
                'method': 'default', 'param': None,
                'x': None, 'y': None, 'z1': None, 'z2': None}


class dc_gcp(data_container):
    r"""
    Parameters
    ----------
    is_set : boolean
        Should the equation for this parameter group be applied? default: is_set=False.

    method : str
        Which calculation method for this parameter group should be used?, default: method='default'.

    elements : list
        Which component properties are part of this component group? default elements=[].
    """
    def attr(self):
        r"""
        Return the available attributes for a data_container type object.

        Returns
        -------
        out : dict
            Dictionary of available attributes (dictionary keys) with default values.
        """
        return {'is_set': False, 'method': 'default', 'elements': []}

# %%


class TESPyNetworkError(Exception):
    pass


class TESPyConnectionError(Exception):
    pass


class TESPyComponentError(Exception):
    pass

#
#def query_yes_no(question, default='yes'):
#    r"""
#    Parameters
#    ----------
#    question : str
#        Question to be asked.
#
#    default : str
#        Default answer: default='yes'.
#
#    elements : list
#        Which component properties are part of this component group? default elements=[].
#
#    Returns
#    -------
#    answer : boolean
#        Answer.
#    """
#    valid = {'yes': True,
#             'y': True,
#             'ye': True,
#             'no': False,
#             'n': False}
#    if default is None:
#        prompt = '[y / n]'
#    elif default == 'yes':
#        prompt = '[Y / n]'
#    elif default == 'no':
#        prompt = '[y / N]'
#
#    while True:
#        sys.stdout.write(question + prompt)
#        choice = input().lower()
#        if default is not None and choice == '':
#            return valid[default]
#        elif choice in valid:
#            return valid[choice]
#        else:
#            sys.stdout.write('Please respond with \'yes\' or \'no\' (or \'y\' or \'n\').\n')
#
# %%


class tespy_fluid:
    r"""
    The tespy_fluid class allows the creation of custom fluid properies for a
    specified mixture of fluids. The created fluid properties adress an ideal
    mixture of real fluids.

    Parameters
    ----------
    alias : str
        Alias for the fluid. Please note: The alias of a tespy_fluid class object
        will always start with :code:`TESPy::`. See the example for more information!

    fluid : dict
        Fluid vector specifying the fluid composition of the TESPy fluid.

    p_range : list/ndarray
        Pressure range for the new fluid lookup table.

    T_range : list/ndarray
        Temperature range for the new fluid lookup table.

    path : str
        Path to importing tespy fluid from.

    plot : boolean
        Plot the lookup tables after creation?

    Note
    ----
    Creates lookup tables for

    - enthalpy (h),
    - entropy (s),
    - density (d) and
    - viscoity (visc)

    from pressure and temperature. Additionally molar mass and gas constant
    will be calculated. Inverse functions, e. g. entropy from pressure and
    enthalpy are calculated via newton algorithm from these tables.

    Example
    -------
    >>> from tespy import con, cmp, hlp, nwk
    >>> import shutil
    >>> fluidvec = {'CO2': 0.05, 'H2O': 0.06, 'O2': 0.10, 'N2': 0.76, 'Ar': 0.01}
    >>> p = np.array([0.1, 10]) * 1e5
    >>> T = np.array([280, 1280])
    >>> myfluid = hlp.tespy_fluid('flue gas', fluid=fluidvec, p_range=p, T_range=T)
    >>> type(myfluid)
    <class 'tespy.tools.helpers.tespy_fluid'>
    >>> nw = nwk.network(fluids=[myfluid.alias], h_unit='kJ / kg', T_unit='C', p_unit='bar')
    >>> nw.set_printoptions(print_level='none')
    >>> source = cmp.source('source')
    >>> sink = cmp.sink('sink')
    >>> c = con.connection(source, 'out1', sink, 'in1')
    >>> nw.add_conns(c)
    >>> c.set_attr(m=1, T=500, p=5, fluid={'TESPy::flue gas': 1})
    >>> nw.solve('design')
    >>> round(hlp.h_mix_pT(c.to_flow(), c.T.val_SI), 0)
    957564.0
    >>> loadfluid = hlp.tespy_fluid('flue gas', fluid=fluidvec, p_range=p, T_range=T, path='./LUT')
    >>> shutil.rmtree('./LUT', ignore_errors=True)
    """

    def __init__(self, alias, fluid, p_range, T_range, path=None, plot=False):

        if not isinstance(alias, str):
            msg = 'Alias must be of type String.'
            logging.error(msg)
            raise TypeError(msg)

        if 'IDGAS::' in alias:
            msg = 'You are not allowed to use "IDGAS::" within your alias.'
            logging.error(msg)
            raise ValueError(msg)

        # process parameters
        if 'TESPy::' in alias:
            self.alias = alias
        else:
            self.alias = 'TESPy::' + alias
        self.fluid = fluid

        # adjust value ranges according to specified unit system
        self.p_range = np.array(p_range)
        self.T_range = np.array(T_range)

        # set up grid
        self.p = np.linspace(self.p_range[0], self.p_range[1])
        self.T = np.linspace(self.T_range[0], self.T_range[1])

        # plotting
        self.plot = plot

        # path for loading
        self.path = path

        # calculate molar mass and gas constant
        for f in self.fluid:
            molar_masses[f] = CPPSI('M', f)
            gas_constants[f] = CPPSI('GAS_CONSTANT', f)

        molar_masses[self.alias] = 1 / molar_mass_flow(self.fluid)
        gas_constants[self.alias] = (gas_constants['uni'] / molar_masses[self.alias])

        # create look up tables
        tespy_fluid.fluids[self.alias] = {}
        memorise.add_fluids(self.fluid.keys())

        params = {}

        params['h_pT'] = h_mix_pT
        params['s_pT'] = s_mix_pT
        params['d_pT'] = d_mix_pT
        params['visc_pT'] = visc_mix_pT

        self.funcs = {}

        if self.path is None:
            # generate fluid properties
            msg = 'Generating lookup-tables from CoolProp fluid properties.'
            logging.debug(msg)
            for key in params.keys():
                self.funcs[key] = self.generate_lookup(key, params[key])
                msg = 'Loading function values for function ' + key + '.'
                logging.debug(msg)

        else:
            # load fluid properties from specified path
            msg = 'Generating lookup-tables from base path ' + self.path + '/' + self.alias + '/.'
            logging.debug(msg)
            for key in params.keys():
                self.funcs[key] = self.load_lookup(key)
                msg = 'Loading function values for function ' + key + '.'
                logging.debug(msg)

        tespy_fluid.fluids[self.alias] = self

        msg = 'Successfully created look-up-tables for custom fluid ' + self.alias + '.'
        logging.debug(msg)


    def generate_lookup(self, name, func):
        r"""
        Create lookup table from CoolProp-database

        Parameters
        ----------
        name : str
            Name of the lookup table.

        func : function
            Function to create lookup table from.
        """
        x1 = self.p
        x2 = self.T

        y = np.empty((0, x1.shape[0]), float)

        # iterate
        for p in x1:
            row = []
            for T in x2:
                row += [func([0, p, 0, self.fluid], T)]

            y = np.append(y, [np.array(row)], axis=0)
#
#        # plot table after creation?
#        if self.plot:
#            fig = plt.figure()
#            ax = fig.add_subplot(111, projection='3d')
#            ax.plot_wireframe(np.meshgrid(x2, x1)[0],
#                              np.meshgrid(x2, x1)[1], y)
#            ax.set_xlabel('temperature')
#            ax.set_ylabel('pressure')
#            ax.set_zlabel(name)
#            ax.view_init(10, 225)
#            plt.show()

        self.save_lookup(name, x1, x2, y)

        func = interpolate.RectBivariateSpline(x1, x2, y)
        return func

    def save_lookup(self, name, x1, x2, y):
        r"""
        Save lookup table to working dir in new folder :code:`./LUT/fluid_alias/`.

        Parameters
        ----------
        name : str
            Name of the lookup table.

        x1 : ndarray
            Pressure.

        x2 : ndarray
            Temperature.

        y : ndarray
            Lookup value (enthalpy, entropy, density or viscosity)
        """
        df = pd.DataFrame(y, columns=x2, index=x1)
        alias = self.alias.replace('::','_')
        path = './LUT/' + alias + '/'
        if not os.path.exists(path):
            os.makedirs(path)
        df.to_csv(path + name + '.csv')

    def load_lookup(self, name):
        r"""
        Load lookup table from specified path base path and alias.

        Parameters
        ----------
        name : str
            Name of the lookup table.
        """
        alias = self.alias.replace('::','_')
        path = self.path + '/' + alias + '/' + name + '.csv'
        df = pd.read_csv(path, index_col=0)

        x1 = df.index.get_values()
        x2 = np.array(list(map(float, list(df))))
        y = df.values
#
#        # plot table after creation?
#        if self.plot:
#            fig = plt.figure()
#            ax = fig.add_subplot(111, projection='3d')
#            ax.plot_wireframe(np.meshgrid(x2, x1)[0], np.meshgrid(x2, x1)[1], y)
#            ax.set_xlabel('temperature')
#            ax.set_ylabel('pressure')
#            ax.set_zlabel(name)
#            ax.view_init(10, 225)
#            plt.show()

        func = interpolate.RectBivariateSpline(x1, x2, y)
        return func

# create dict for tespy fluids
tespy_fluid.fluids = {}


def reverse_2d(params, y):
    r"""
    Reverse function for lookup table.

    Parameters
    ----------
    params : list
        Variable function parameters.

    y : float
        Function value of function :math:`y = f \left( x_1, x_2 \right)`.

    Returns
    ------
    deriv : float
        Residual value of inverse function :math:`x_2 - f\left(x_1, y \right)`.
    """
    func, x1, x2 = params[0], params[1], params[2]
    return x2 - func.ev(x1, y)


def reverse_2d_deriv(params, y):
    r"""
    Derivative of the reverse function for a lookup table.

    params : list
        Variable function parameters.

    y : float
        Function value of function :math:`y = f \left( x_1, x_2 \right)`,
        so that :math:`x_2 - f\left(x_1, y \right) = 0`

    Returns
    ------
    deriv : float
        Partial derivative :math:`\frac{\partial f}{\partial y}`.
    """
    func, x1 = params[0], params[1]
    return - func.ev(x1, y, dy=1)

# %%


class memorise:
    r"""
    Memorization of fluid properties.
    """

    def add_fluids(fluids):
        r"""
        Add list of fluids to fluid memorisation class.

        - Generate arrays for fluid property lookup.
        - Calculate/set fluid property value ranges for convergence checks.

        Parameters
        ----------
        fluids : list
            List of fluid for fluid property memorization.

        Note
        ----
        The memorise class creates globally accessible variables for different fluid
        property calls as dictionaries:

            - T(p,h)
            - T(p,s)
            - v(p,h)
            - visc(p,h)
            - s(p,h)

        Each dictionary uses the list of fluids passed to the memorise class as
        identifier for the fluid property memorisation. The fluid properties are
        stored as numpy array, where each column represents the mass fraction of the
        respective fluid and the additional columns are the values for the fluid
        properties. The fluid property function will then look for identical fluid
        property inputs (p, h, (s), fluid mass fraction). If the inputs are in the
        array, the first column of that row is returned, see example.

        Example
        -------
        T(p,h) for set of fluids ('water', 'air'):

            - row 1: [282.64527752319697, 10000, 40000, 1, 0]
            - row 2: [284.3140698256616, 10000, 47000, 1, 0]
        """
        # number of fluids
        num_fl = len(fluids)
        if num_fl > 0:
            fl = tuple(fluids)
            # fluid property tables
            memorise.T_ph[fl] = np.empty((0, num_fl + 3), float)
            memorise.T_ps[fl] = np.empty((0, num_fl + 4), float)
            memorise.v_ph[fl] = np.empty((0, num_fl + 3), float)
            memorise.visc_ph[fl] = np.empty((0, num_fl + 3), float)
            memorise.s_ph[fl] = np.empty((0, num_fl + 3), float)
            # lists for memory cache, values not in these lists will be deleted
            # from the table after every tespy.networks.network.solve call.
            memorise.T_ph_f[fl] = []
            memorise.T_ps_f[fl] = []
            memorise.v_ph_f[fl] = []
            memorise.visc_ph_f[fl] = []
            memorise.s_ph_f[fl] = []
            memorise.count = 0

            msg = 'Added fluids ' + str(fl) +' to memorise lookup tables.'
            logging.debug(msg)

        for f in fluids:
            if 'TESPy::' in f:
                if f in tespy_fluid.fluids.keys():
                    pmin, pmax = tespy_fluid.fluids[f].p_range[0], tespy_fluid.fluids[f].p_range[1]
                    Tmin, Tmax = tespy_fluid.fluids[f].T_range[0], tespy_fluid.fluids[f].T_range[1]
                    msg = 'Loading fluid property ranges for TESPy-fluid ' + f + '.'
                    logging.debug(msg)
                    # value range for fluid properties
                    memorise.vrange[f] = [pmin, pmax, Tmin, Tmax]
                    msg = 'Specifying fluid property ranges for pressure and temperature for convergence check.'
                    logging.debug(msg)
                else:
                    memorise.vrange[f] = [2000, 2000000, 300, 2000]

            elif 'INCOMP::' in f:
                # temperature range available only for incompressibles
                Tmin, Tmax = CPPSI('TMIN', f), CPPSI('TMAX', f)
                memorise.vrange[f] = [2000, 2000000, Tmin, Tmax]

            else:
                if f not in memorise.heos.keys():
                    # abstractstate object
                    memorise.heos[f] = CP.AbstractState('HEOS', f)
                    msg = 'Created CoolProp.AbstractState object for fluid ' + f + ' in memorise class.'
                    logging.debug(msg)
                    # pressure range
                    pmin, pmax = CPPSI('PMIN', f), CPPSI('PMAX', f)
                    # temperature range
                    Tmin, Tmax = CPPSI('TMIN', f), CPPSI('TMAX', f)
                    # value range for fluid properties
                    memorise.vrange[f] = [pmin, pmax, Tmin, Tmax]
                    msg = 'Specifying fluid property ranges for pressure and temperature for convergence check of fluid ' + f + '.'
                    logging.debug(msg)


    def del_memory(fluids):
        r"""
        Deletes non frequently used fluid property values from memorise class.

        Parameters
        ----------
        fluids : list
            List of fluid for fluid property memorization.
        """

        fl = tuple(fluids)

        # delete memory
        mask = np.isin(memorise.T_ph[fl][:, -1], memorise.T_ph_f[fl])
        memorise.T_ph[fl] = (memorise.T_ph[fl][mask])

        mask = np.isin(memorise.T_ps[fl][:, -1], memorise.T_ps_f[fl])
        memorise.T_ps[fl] = (memorise.T_ps[fl][mask])

        mask = np.isin(memorise.v_ph[fl][:, -1], memorise.v_ph_f[fl])
        memorise.v_ph[fl] = (memorise.v_ph[fl][mask])

        mask = np.isin(memorise.visc_ph[fl][:, -1], memorise.visc_ph_f[fl])
        memorise.visc_ph[fl] = (memorise.visc_ph[fl][mask])

        mask = np.isin(memorise.s_ph[fl][:, -1], memorise.s_ph_f[fl])
        memorise.s_ph[fl] = (memorise.s_ph[fl][mask])

        # refresh cache
        memorise.T_ph_f[fl] = []
        memorise.T_ps_f[fl] = []
        memorise.v_ph_f[fl] = []
        memorise.visc_ph_f[fl] = []
        memorise.s_ph_f[fl] = []

        msg = 'Dropping not frequently used fluid property values from memorise class for fluids ' + str(fl) + '.'
        logging.debug(msg)


# create memorise dictionaries
memorise.heos = {}
memorise.T_ph = {}
memorise.T_ph_f = {}
memorise.T_ps = {}
memorise.T_ps_f = {}
memorise.v_ph = {}
memorise.v_ph_f = {}
memorise.visc_ph = {}
memorise.visc_ph_f = {}
memorise.s_ph = {}
memorise.s_ph_f = {}
memorise.vrange = {}

# %%


def newton(func, deriv, params, y, **kwargs):
    r"""
    1-D newton algorithm to find zero crossings of function func with its derivative
    deriv.

    Parameters
    ----------
    func : function
        Function to find zero crossing in, :math:`0=y-func\left(x,\text{params}\right)`.

    deriv : function
        First derivative of the function.

    params : list
        Additional parameters for function, optional.

    y : float
        Target function value.

    val0 : float
        Starting value, default: val0=300.

    valmin : float
        Lower value boundary, default: valmin=70.

    valmax : float
        Upper value boundary, default: valmax=3000.

    max_iter : int
        Maximum number of iterations, default: max_iter=10.

    Returns
    -------
    val : float
        x-value of zero crossing.

    Note
    ----
    Algorithm

    .. math::

        x_{i+1} = x_{i} - \frac{f(x_{i})}{\frac{df}{dx}(x_{i})}\\
        f(x_{i}) \leq \epsilon
    """
    # default valaues
    x = kwargs.get('val0', 300)
    valmin = kwargs.get('valmin', 70)
    valmax = kwargs.get('valmax', 3000)
    max_iter = kwargs.get('max_iter', 10)

    # start newton loop
    res = 1
    i = 0
    while abs(res) >= err:
        # calculate function residual and new value
        res = y - func(params, x)
        x += res / deriv(params, x)

        # check for value ranges
        if x < valmin:
            x = valmin
        if x > valmax:
            x = valmax
        i += 1

        if i > max_iter:
            msg = ('Newton algorithm was not able to find a feasible value for function ' + str(func) + '. '
                   'Current value with x=' + str(x) + ' is ' + str(func(params, x)) + ', target value is ' + str(y) + '.')
            logging.debug(msg)

            break

    return x

# %%


def T_mix_ph(flow):
    r"""
    Calculates the temperature from pressure and enthalpy.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and fluid composition.

    Returns
    -------
    T : float
        Temperature T / K.

    Note
    ----
    First, check if fluid property has been memorised already.
    If this is the case, return stored value, otherwise calculate value and
    store it in the memorisation class.

    Uses CoolProp interface for pure fluids, newton algorithm for mixtures:

    .. math::

        T_{mix}\left(p,h\right) = T_{i}\left(p,h_{i}\right)\;
        \forall i \in \text{fluid components}\\

        h_{i} = h \left(pp_{i}, T_{mix} \right)\\
        pp: \text{partial pressure}

    """
    # check if fluid properties have been calculated before
    fl = tuple(flow[3].keys())
    a = memorise.T_ph[fl][:, 0:-1]
    b = np.array([flow[1], flow[2]] + list(flow[3].values()))
    ix = np.where(np.all(abs(a - b) <= err, axis=1))[0]

    if ix.size == 1:
        # known fluid properties
        T = memorise.T_ph[fl][ix, -1][0]
        memorise.T_ph_f[fl] += [T]
        return T
    else:
        # unknown fluid properties
        if num_fluids(flow[3]) > 1:
            # calculate the fluid properties for fluid mixtures
            val = newton(h_mix_pT, dh_mix_pdT, flow, flow[2], val0=300, valmin=70, valmax=3000, imax=10)
            new = np.array([[flow[1], flow[2]] + list(flow[3].values()) + [val]])
            # memorise the newly calculated value
            memorise.T_ph[fl] = np.append(memorise.T_ph[fl], new, axis=0)
            return val
        else:
            # calculate fluid property for pure fluids
            for fluid, x in flow[3].items():
                if x > err:
                    val = T_ph(flow[1], flow[2], fluid)
                    new = np.array([[flow[1], flow[2]] + list(flow[3].values()) + [val]])
                    # memorise the newly calculated value
                    memorise.T_ph[fl] = np.append(memorise.T_ph[fl], new, axis=0)
                    return val


def T_ph(p, h, fluid):
    r"""
    Calculates the temperature from pressure and enthalpy for a pure fluid.

    Parameters
    ----------
    p : float
        Pressure p / Pa.

    h : float
        Specific enthalpy h / (J/kg).

    fluid : str
        Fluid name.

    Returns
    -------
    T : float
        Temperature T / K.
    """
    if 'IDGAS::' in fluid:
        print('Ideal gas calculation not available by now.')
    elif 'TESPy::' in fluid:
        db = tespy_fluid.fluids[fluid].funcs['h_pT']
        return newton(reverse_2d, reverse_2d_deriv, [db, p, h], 0)
    elif 'INCOMP::' in fluid:
        return CPPSI('T', 'P', p, 'H', h, fluid)
    else:
        memorise.heos[fluid].update(CP.HmassP_INPUTS, h, p)
        return memorise.heos[fluid].T()


def dT_mix_dph(flow):
    r"""
    Calculate partial derivate of temperature to pressure at constant enthalpy and fluid composition.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and fluid composition.

    Returns
    -------
    dT / dp : float
        Partial derivative of temperature to pressure dT /dp / (K/Pa).

        .. math::

            \frac{\partial T_{mix}}{\partial p} = \frac{T_{mix}(p+d,h)-
            T_{mix}(p-d,h)}{2 \cdot d}
    """
    d = 1
    u = flow.copy()
    l = flow.copy()
    u[1] += d
    l[1] -= d
    return (T_mix_ph(u) - T_mix_ph(l)) / (2 * d)


def dT_mix_pdh(flow):
    r"""
    Calculate partial derivate of temperature to enthalpy at constant pressure and fluid composition.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and fluid composition.

    Returns
    -------
    dT / dh : float
        Partial derivative of temperature to enthalpy dT /dh / ((kgK)/J).

        .. math::

            \frac{\partial T_{mix}}{\partial h} = \frac{T_{mix}(p,h+d)-
            T_{mix}(p,h-d)}{2 \cdot d}
    """
    d = 1
    u = flow.copy()
    l = flow.copy()
    u[2] += d
    l[2] -= d
    return (T_mix_ph(u) - T_mix_ph(l)) / (2 * d)


def dT_mix_ph_dfluid(flow):
    r"""
    Calculate partial derivate of temperature to fluid composition at constant pressure and enthalpy.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and fluid composition.

    Returns
    -------
    dT / dfluid : ndarray
        Partial derivatives of temperature to fluid composition dT / dfluid / K.

        .. math::

            \frac{\partial T_{mix}}{\partial fluid_{i}} =
            \frac{T_{mix}(p,h,fluid_{i}+d)-
            T_{mix}(p,h,fluid_{i}-d)}{2 \cdot d}
    """
    d = 1e-5
    u = flow.copy()
    l = flow.copy()
    vec_deriv = []
    for fluid, x in flow[3].items():
        if x > err:
            u[3][fluid] += d
            l[3][fluid] -= d
            vec_deriv += [(T_mix_ph(u) - T_mix_ph(l)) / (2 * d)]
            u[3][fluid] -= d
            l[3][fluid] += d
        else:
            vec_deriv += [0]

    return np.asarray(vec_deriv)

# %%


def T_mix_ps(flow, s):
    r"""
    Calculates the temperature from pressure and entropy.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and fluid composition.

    s : float
        Entropy of flow in J / (kgK).

    Returns
    -------
    T : float
        Temperature T / K.

    Note
    ----
    First, check if fluid property has been memorised already.
    If this is the case, return stored value, otherwise calculate value and
    store it in the memorisation class.

    Uses CoolProp interface for pure fluids, newton algorithm for mixtures:

    .. math::

        T_{mix}\left(p,s\right) = T_{i}\left(p,s_{i}\right)\;
        \forall i \in \text{fluid components}\\

        s_{i} = s \left(pp_{i}, T_{mix} \right)\\
        pp: \text{partial pressure}

    """
    # check if fluid properties have been calculated before
    fl = tuple(flow[3].keys())
    a = memorise.T_ps[fl][:, 0:-1]
    b = np.array([flow[1], flow[2]] + list(flow[3].values()) + [s])
    ix = np.where(np.all(abs(a - b) <= err, axis=1))[0]
    if ix.size == 1:
        # known fluid properties
        T = memorise.T_ps[fl][ix, -1][0]
        memorise.T_ps_f[fl] += [T]
        return T
    else:
        # unknown fluid properties
        if num_fluids(flow[3]) > 1:
            # calculate the fluid properties for fluid mixtures
            val = newton(s_mix_pT, ds_mix_pdT, flow, s, val0=300, valmin=70, valmax=3000, imax=10)
            new = np.array([[flow[1], flow[2]] + list(flow[3].values()) + [s, val]])
            # memorise the newly calculated value
            memorise.T_ps[fl] = np.append(memorise.T_ps[fl], new, axis=0)
            return val
        else:
            # calculate fluid property for pure fluids
            for fluid, x in flow[3].items():
                if x > err:
                    val = T_ps(flow[1], s, fluid)
                    new = np.array([[flow[1], flow[2]] + list(flow[3].values()) + [s, val]])
                    # memorise the newly calculated value
                    memorise.T_ps[fl] = np.append(memorise.T_ps[fl], new, axis=0)
                    return val


def T_ps(p, s, fluid):
    r"""
    Calculates the temperature from pressure and entropy for a pure fluid.

    Parameters
    ----------
    p : float
        Pressure p / Pa.

    s : float
        Specific entropy h / (J/(kgK)).

    fluid : str
        Fluid name.

    Returns
    -------
    T : float
        Temperature T / K.
    """
    if 'IDGAS::' in fluid:
        print('Ideal gas calculation not available by now.')
    elif 'TESPy::' in fluid:
        db = tespy_fluid.fluids[fluid].funcs['s_pT']
        return newton(reverse_2d, reverse_2d_deriv, [db, p, s], 0)
    elif 'INCOMP::' in fluid:
        return CPPSI('T', 'P', p, 'H', s, fluid)
    else:
        memorise.heos[fluid].update(CP.PSmass_INPUTS, p, s)
        return memorise.heos[fluid].T()

# %%


def h_mix_pT(flow, T):
    r"""
    Calculates the enthalpy from pressure and Temperature.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and fluid composition.

    T : float
        Temperature of flow T / K.

    Returns
    -------
    h : float
        Enthalpy h / (J/kg).

    Note
    ----
    Calculation for fluid mixtures.

    .. math::

        h_{mix}(p,T)=\sum_{i} h(pp_{i},T,fluid_{i})\;
        \forall i \in \text{fluid components}\\
        pp: \text{partial pressure}
    """
    n = molar_mass_flow(flow[3])

    h = 0
    for fluid, x in flow[3].items():
        if x > err:
            ni = x / molar_masses[fluid]
            h += h_pT(flow[1] * ni / n, T, fluid) * x

    return h


def h_pT(p, T, fluid):
    r"""
    Calculates the enthalpy from pressure and temperature for a pure fluid.

    Parameters
    ----------
    p : float
        Pressure p / Pa.

    T : float
        Temperature T / K.

    fluid : str
        Fluid name.

    Returns
    -------
    h : float
        Specific enthalpy h / (J/kg).
    """
    if 'IDGAS::' in fluid:
        print('Ideal gas calculation not available by now.')
    elif 'TESPy::' in fluid:
        return tespy_fluid.fluids[fluid].funcs['h_pT'].ev(p, T)
    elif 'INCOMP::' in fluid:
        return CPPSI('H', 'P', p, 'T', T, fluid)
    else:
        memorise.heos[fluid].update(CP.PT_INPUTS, p, T)
        return memorise.heos[fluid].hmass()


def dh_mix_pdT(flow, T):
    r"""
    Calculate partial derivate of enthalpy to temperature at constant pressure and fluid composition.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and fluid composition.

    T : float
        Temperature T / K.

    Returns
    -------
    dh / dT : float
        Partial derivative of enthalpy to temperature dh / dT / (J/(kgK)).

        .. math::

            \frac{\partial h_{mix}}{\partial T} =
            \frac{h_{mix}(p,T+d)-h_{mix}(p,T-d)}{2 \cdot d}
    """
    d = 2
    return (h_mix_pT(flow, T + d) - h_mix_pT(flow, T - d)) / (2 * d)

# %%


def h_mix_ps(flow, s):
    r"""
    Calculates the enthalpy from pressure and temperature.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and fluid composition.

    s : float
        Specific entropy of flow s / (J/(kgK)).

    Returns
    -------
    h : float
        Specific enthalpy h / (J/kg).

    Note
    ----
    Calculation for fluid mixtures.

    .. math::

        h_{mix}\left(p,s\right)=h\left(p, T_{mix}\left(p,s\right)\right)
    """
    return h_mix_pT(flow, T_mix_ps(flow, s))


def h_ps(p, s, fluid):
    r"""
    Calculates the enthalpy from pressure and entropy for a pure fluid.

    Parameters
    ----------
    p : float
        Pressure p / Pa.

    s : float
        Specific entropy h / (J/(kgK)).

    fluid : str
        Fluid name.

    Returns
    -------
    h : float
        Specific enthalpy h / (J/kg).
    """
    if 'IDGAS::' in fluid:
        print('Ideal gas calculation not available by now.')
    elif 'TESPy::' in fluid:
        db = tespy_fluid.fluids[fluid].funcs['s_pT']
        T = newton(reverse_2d, reverse_2d_deriv, [db, p, s], 0)
        return tespy_fluid.fluids[fluid].funcs['h_pT'].ev(p, T)
    elif 'INCOMP::' in fluid:
        return CPPSI('H', 'P', p, 'S', s, fluid)
    else:
        memorise.heos[fluid].update(CP.PSmass_INPUTS, p, s)
        return memorise.heos[fluid].hmass()

# %%


def h_mix_pQ(flow, Q):
    r"""
    Calculates the enthalpy from pressure and vapour mass fraction.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and fluid composition.

    Q : float
        Vapour mass fraction Q / 1.

    Returns
    -------
    h : float
        Specific enthalpy h / (J/kg).

    Note
    ----
    This function works for pure fluids only!
    """
    n = molar_mass_flow(flow[3])

    h = 0
    for fluid, x in flow[3].items():
        if x > err:
            pp = flow[1] * x / (molar_masses[fluid] * n)
            pcrit = CPPSI('Pcrit', fluid)
            if pp > pcrit:
                pp = pcrit * 0.95

            memorise.heos[fluid].update(CP.PQ_INPUTS, pp, Q)
            h += memorise.heos[fluid].hmass() * x

    return h


def dh_mix_dpQ(flow, Q):
    r"""
    Calculate partial derivate of enthalpy to vapour mass fraction at constant pressure.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and fluid composition.

    Q : float
        Vapour mass fraction Q / 1.

    Returns
    -------
    dh / dQ : float
        Partial derivative of enthalpy to vapour mass fraction dh / dQ / (J/kg).

        .. math::

            \frac{\partial h_{mix}}{\partial p} =
            \frac{h_{mix}(p+d,Q)-h_{mix}(p-d,Q)}{2 \cdot d}\\
            Q: \text{vapour mass fraction}

    Note
    ----
    This works for pure fluids only!
    """
    d = 1
    u = flow.copy()
    l = flow.copy()
    u[1] += d
    l[1] -= d
    return (h_mix_pQ(u, Q) - h_mix_pQ(l, Q)) / (2 * d)

# %%


def T_bp_p(flow):
    r"""
    Calculates temperature from boiling point pressure.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and fluid composition.

    Returns
    -------
    T : float
        Temperature at boiling point.

    Note
    ----
    This function works for pure fluids only!
    """
    n = molar_mass_flow(flow[3])

    for fluid, x in flow[3].items():
        if x > err:
            pp = flow[1] * x / (molar_masses[fluid] * n)

            memorise.heos[fluid].update(CP.PQ_INPUTS, pp, 1)
            return memorise.heos[fluid].T()


def dT_bp_dp(flow):
    r"""
    Calculate partial derivate of temperature to boiling point pressure.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and fluid composition.

    Returns
    -------
    dT / dp : float
        Partial derivative of temperature to boiling point pressure in K / Pa.

        .. math::

            \frac{\partial h_{mix}}{\partial p} =
            \frac{T_{bp}(p+d)-T_{bp}(p-d)}{2 \cdot d}\\
            Q: \text{vapour mass fraction}

    Note
    ----
    This works for pure fluids only!
    """
    d = 1
    u = flow.copy()
    l = flow.copy()
    u[1] += d
    l[1] -= d
    return (T_bp_p(u) - T_bp_p(l)) / (2 * d)

# %%


def v_mix_ph(flow):
    r"""
    Calculates the specific volume from pressure and enthalpy.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and fluid composition.

    Returns
    -------
    v : float
        Specific volume v / (:math:`\mathrm{m}^3`/kg).

    Note
    ----
    First, check if fluid property has been memorised already.
    If this is the case, return stored value, otherwise calculate value and
    store it in the memorisation class.

    Uses CoolProp interface for pure fluids, newton algorithm for mixtures:

    .. math::

        v_{mix}\left(p,h\right) = v\left(p,T_{mix}(p,h)\right)
    """
    # check if fluid properties have been calculated before
    fl = tuple(flow[3].keys())
    a = memorise.v_ph[fl][:, 0:-1]
    b = np.array([flow[1], flow[2]] + list(flow[3].values()))
    ix = np.where(np.all(abs(a - b) <= err, axis=1))[0]
    if ix.size == 1:
        # known fluid properties
        v = memorise.v_ph[fl][ix, -1][0]
        memorise.v_ph_f[fl] += [v]
        return v
    else:
        # unknown fluid properties
        if num_fluids(flow[3]) > 1:
            # calculate the fluid properties for fluid mixtures
            val = v_mix_pT(flow, T_mix_ph(flow))
            new = np.array([[flow[1], flow[2]] + list(flow[3].values()) + [val]])
            # memorise the newly calculated value
            memorise.v_ph[fl] = np.append(memorise.v_ph[fl], new, axis=0)
            return val
        else:
            # calculate fluid property for pure fluids
            for fluid, x in flow[3].items():
                if x > err:
                    val = 1 / d_ph(flow[1], flow[2], fluid)
                    new = np.array([[flow[1], flow[2]] + list(flow[3].values()) + [val]])
                    # memorise the newly calculated value
                    memorise.v_ph[fl] = np.append(memorise.v_ph[fl], new, axis=0)
                    return val


def d_ph(p, h, fluid):
    r"""
    Calculates the density from pressure and enthalpy for a pure fluid.

    Parameters
    ----------
    p : float
        Pressure p / Pa.

    h : float
        Specific enthalpy h / (J/kg).

    fluid : str
        Fluid name.

    Returns
    -------
    d : float
        Density d / (kg/:math:`\mathrm{m}^3`).
    """
    if 'IDGAS::' in fluid:
        print('Ideal gas calculation not available by now.')
    elif 'TESPy::' in fluid:
        db = tespy_fluid.fluids[fluid].funcs['h_pT']
        T = newton(reverse_2d, reverse_2d_deriv, [db, p, h], 0)
        return tespy_fluid.fluids[fluid].funcs['d_pT'].ev(p, T)
    elif 'INCOMP::' in fluid:
        return CPPSI('D', 'P', p, 'H', h, fluid)
    else:
        memorise.heos[fluid].update(CP.HmassP_INPUTS, h, p)
        return memorise.heos[fluid].rhomass()


def dv_mix_dph(flow):
    r"""
    Calculate partial derivate of specific volume to pressure at constant enthalpy and fluid composition.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and fluid composition.

    Returns
    -------
    dv / dp : float
        Partial derivative of specific volume to pressure dv /dp / (:math:`\mathrm{m}^3`/(Pa kg)).

        .. math::

            \frac{\partial v_{mix}}{\partial p} = \frac{v_{mix}(p+d,h)-
            v_{mix}(p-d,h)}{2 \cdot d}
    """
    d = 1
    u = flow.copy()
    l = flow.copy()
    u[1] += d
    l[1] -= d
    return (v_mix_ph(u) - v_mix_ph(l)) / (2 * d)


def dv_mix_pdh(flow):
    r"""
    Calculate partial derivate of specific volume to enthalpy at constant pressure and fluid composition.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and fluid composition.

    Returns
    -------
    dv / dh : float
        Partial derivative of specific volume to enthalpy dv /dh / (:math:`\mathrm{m}^3`/J).

        .. math::

            \frac{\partial v_{mix}}{\partial h} = \frac{v_{mix}(p,h+d)-
            v_{mix}(p,h-d)}{2 \cdot d}
    """
    d = 1
    u = flow.copy()
    l = flow.copy()
    u[2] += d
    l[2] -= d
    return (v_mix_ph(u) - v_mix_ph(l)) / (2 * d)

# %%


def v_mix_pT(flow, T):
    r"""
    Calculates the specific volume from pressure and temperature.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and fluid composition.

    T : float
        Temperature T / K.

    Returns
    -------
    v : float
        Specific volume v / (:math:`\mathrm{m}^3`/kg).

    Note
    ----
    Calculation for fluid mixtures.

    .. math::

        v_{mix}(p,T)=\sum_{i} \frac{x_i}{\rho(p, T, fluid_{i})}
    """
    v = 0
    for fluid, x in flow[3].items():
        if x > err:
            v += x / d_pT(flow[1], T, fluid)

    return v


def d_mix_pT(flow, T):
    r"""
    Calculates the density from pressure and temperature.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and fluid composition.

    T : float
        Temperature T / K.

    Returns
    -------
    d : float
        Density d / (kg/:math:`\mathrm{m}^3`).

    Note
    ----
    Calculation for fluid mixtures.

    .. math::

        \rho_{mix}\left(p,T\right)=\frac{1}{v_{mix}\left(p,T\right)}
    """
    return 1 / v_mix_pT(flow, T)


def d_pT(p, T, fluid):
    r"""
    Calculates the density from pressure and temperature for a pure fluid.

    Parameters
    ----------
    p : float
        Pressure p / Pa.

    T : float
        Temperature T / K.

    fluid : str
        Fluid name.

    Returns
    -------
    d : float
        Density d / (kg/:math:`\mathrm{m}^3`).
    """
    if 'IDGAS::' in fluid:
        print('Ideal gas calculation not available by now.')
    elif 'TESPy::' in fluid:
        return tespy_fluid.fluids[fluid].funcs['d_pT'].ev(p, T)
    elif 'INCOMP::' in fluid:
        return CPPSI('D', 'P', p, 'T', T, fluid)
    else:
        memorise.heos[fluid].update(CP.PT_INPUTS, p, T)
        return memorise.heos[fluid].rhomass()

# %%


def visc_mix_ph(flow):
    r"""
    Calculates the dynamic viscorsity from pressure and enthalpy.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and fluid composition.

    Returns
    -------
    visc : float
        Dynamic viscosity visc / Pa s.

    Note
    ----
    First, check if fluid property has been memorised already.
    If this is the case, return stored value, otherwise calculate value and
    store it in the memorisation class.

    Uses CoolProp interface for pure fluids, newton algorithm for mixtures:

    .. math::

        \eta_{mix}\left(p,h\right) = \eta\left(p,T_{mix}(p,h)\right)
    """
    # check if fluid properties have been calculated before
    fl = tuple(flow[3].keys())
    a = memorise.visc_ph[fl][:, 0:-1]
    b = np.array([flow[1], flow[2]] + list(flow[3].values()))
    ix = np.where(np.all(abs(a - b) <= err, axis=1))[0]
    if ix.size == 1:
        # known fluid properties
        visc = memorise.visc_ph[fl][ix, -1][0]
        memorise.visc_ph_f[fl] += [visc]
        return visc
    else:
        # unknown fluid properties
        if num_fluids(flow[3]) > 1:
            # calculate the fluid properties for fluid mixtures
            val = visc_mix_pT(flow, T_mix_ph(flow))
            new = np.array([[flow[1], flow[2]] + list(flow[3].values()) + [val]])
            # memorise the newly calculated value
            memorise.visc_ph[fl] = np.append(memorise.visc_ph[fl], new, axis=0)
            return val
        else:
            # calculate fluid property for pure fluids
            for fluid, x in flow[3].items():
                if x > err:
                    val = visc_ph(flow[1], flow[2], fluid)
                    new = np.array([[flow[1], flow[2]] + list(flow[3].values()) + [val]])
                    # memorise the newly calculated value
                    memorise.visc_ph[fl] = np.append(memorise.visc_ph[fl], new, axis=0)
                    return val


def visc_ph(p, h, fluid):
    r"""
    Calculates the dynamic viscosity from pressure and enthalpy for a pure fluid.

    Parameters
    ----------
    p : float
        Pressure p / Pa.

    h : float
        Specific enthalpy h / (J/kg).

    fluid : str
        Fluid name.

    Returns
    -------
    visc : float
        Viscosity visc / Pa s.
    """
    if 'IDGAS::' in fluid:
        print('Ideal gas calculation not available by now.')
    elif 'TESPy::' in fluid:
        db = tespy_fluid.fluids[fluid].funcs['h_pT']
        T = newton(reverse_2d, reverse_2d_deriv, [db, p, h], 0)
        return tespy_fluid.fluids[fluid].funcs['visc_pT'].ev(p, T)
    elif 'INCOMP::' in fluid:
        return CPPSI('V', 'P', p, 'H', h, fluid)
    else:
        memorise.heos[fluid].update(CP.HmassP_INPUTS, h, p)
        return memorise.heos[fluid].viscosity()

# %%


def visc_mix_pT(flow, T):
    r"""
    Calculates the dynamic viscosity from pressure and temperature.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and fluid composition.

    T : float
        Temperature T / K.

    Returns
    -------
    visc : float
        Dynamic viscosity visc / Pa s.

    Note
    ----
    Calculation for fluid mixtures.

    .. math::

        \eta_{mix}(p,T)=\frac{\sum_{i} \left( \eta(p,T,fluid_{i}) \cdot y_{i}
        \cdot \sqrt{M_{i}} \right)}
        {\sum_{i} \left(y_{i} \cdot \sqrt{M_{i}} \right)}\;
        \forall i \in \text{fluid components}\\
        y: \text{volume fraction}\\
        M: \text{molar mass}
    """
    n = molar_mass_flow(flow[3])

    a = 0
    b = 0
    for fluid, x in flow[3].items():
        if x > err:
            bi = x * math.sqrt(molar_masses[fluid]) / (molar_masses[fluid] * n)
            b += bi
            a += bi * visc_pT(flow[1], T, fluid)

    return a / b


def visc_pT(p, T, fluid):
    r"""
    Calculates the dynamic viscosity from pressure and temperature for a pure fluid.

    Parameters
    ----------
    p : float
        Pressure p / Pa.

    T : float
        Temperature T / K.

    fluid : str
        Fluid name.

    Returns
    -------
    visc : float
        Viscosity visc / Pa s.
    """
    if 'IDGAS::' in fluid:
        print('Ideal gas calculation not available by now.')
    elif 'TESPy::' in fluid:
        return tespy_fluid.fluids[fluid].funcs['visc_pT'].ev(p, T)
    elif 'INCOMP::' in fluid:
        return CPPSI('V', 'P', p, 'T', T, fluid)
    else:
        memorise.heos[fluid].update(CP.PT_INPUTS, p, T)
        return memorise.heos[fluid].viscosity()

# %%


def s_mix_ph(flow):
    r"""
    Calculates the entropy from pressure and enthalpy.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and fluid composition.

    Returns
    -------
    s : float
        Specific entropy s / (J/(kgK)).

    Note
    ----
    First, check if fluid property has been memorised already.
    If this is the case, return stored value, otherwise calculate value and
    store it in the memorisation class.

    Uses CoolProp interface for pure fluids, newton algorithm for mixtures:

    .. math::

        s_{mix}\left(p,h\right) = s\left(p,T_{mix}(p,h)\right)
    """
    # check if fluid properties have been calculated before
    fl = tuple(flow[3].keys())
    a = memorise.s_ph[fl][:, 0:-1]
    b = np.array([flow[1], flow[2]] + list(flow[3].values()))
    ix = np.where(np.all(abs(a - b) <= err, axis=1))[0]
    if ix.size == 1:
        # known fluid properties
        s = memorise.s_ph[fl][ix, -1][0]
        memorise.s_ph_f[fl] += [s]
        return s
    else:
        # unknown fluid properties
        if num_fluids(flow[3]) > 1:
            # calculate the fluid properties for fluid mixtures
            val = s_mix_pT(flow, T_mix_ph(flow))
            new = np.array([[flow[1], flow[2]] + list(flow[3].values()) + [val]])
            # memorise the newly calculated value
            memorise.s_ph[fl] = np.append(memorise.s_ph[fl], new, axis=0)
            return val
        else:
            # calculate fluid property for pure fluids
            for fluid, x in flow[3].items():
                if x > err:
                    val = s_ph(flow[1], flow[2], fluid)
                    new = np.array([[flow[1], flow[2]] + list(flow[3].values()) + [val]])
                    # memorise the newly calculated value
                    memorise.s_ph[fl] = np.append(memorise.s_ph[fl], new, axis=0)
                    return val


def s_ph(p, h, fluid):
    r"""
    Calculates the entropy from pressure and enthalpy for a pure fluid.

    Parameters
    ----------
    p : float
        Pressure p / Pa.

    h : float
        Specific enthalpy h / (J/kg).

    fluid : str
        Fluid name.

    Returns
    -------
    s : float
        Specific entropy s / (J/(kgK)).
    """
    if 'IDGAS::' in fluid:
        print('Ideal gas calculation not available by now.')
    elif 'TESPy::' in fluid:
        db = tespy_fluid.fluids[fluid].funcs['h_pT']
        T = newton(reverse_2d, reverse_2d_deriv, [db, p, h], 0)
        return tespy_fluid.fluids[fluid].funcs['s_pT'].ev(p, T)
    elif 'INCOMP::' in fluid:
        return CPPSI('S', 'P', p, 'H', h, fluid)
    else:
        memorise.heos[fluid].update(CP.HmassP_INPUTS, h, p)
        return memorise.heos[fluid].smass()

# %%


def s_mix_pT(flow, T):
    r"""
    Calculates the entropy from pressure and temperature.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and fluid composition.

    T : float
        Temperature T / K.

    Returns
    -------
    s : float
        Specific entropy s / (J/(kgK)).

    Note
    ----
    Calculation for fluid mixtures.

    .. math::

        s_{mix}(p,T)=\sum_{i} x_{i} \cdot s(pp_{i},T,fluid_{i})-
        \sum_{i} x_{i} \cdot R_{i} \cdot \ln \frac{pp_{i}}{p}\;
        \forall i \in \text{fluid components}\\
        pp: \text{partial pressure}\\
        R: \text{gas constant}
    """
    n = molar_mass_flow(flow[3])

    fluid_comps = {}

    tespy_fluids = [s for s in flow[3].keys() if "TESPy::" in s]
    for f in tespy_fluids:
        for it, x in tespy_fluid.fluids[f].fluid.items():
            if it in fluid_comps.keys():
                fluid_comps[it] += x * flow[3][f]
            else:
                fluid_comps[it] = x * flow[3][f]

    fluids = [s for s in flow[3].keys() if "TESPy::" not in s]
    for f in fluids:
        if f in fluid_comps.keys():
            fluid_comps[f] += flow[3][f]
        else:
            fluid_comps[f] = flow[3][f]

    s = 0
    for fluid, x in fluid_comps.items():
        if x > err:
            pp = flow[1] * x / (molar_masses[fluid] * n)
            s += s_pT(pp, T, fluid) * x
            s -= (x * gas_constants[fluid] / molar_masses[fluid] * math.log(pp / flow[1]))

    return s


def s_pT(p, T, fluid):
    r"""
    Calculates the entropy from pressure and temperature for a pure fluid.

    Parameters
    ----------
    p : float
        Pressure p / Pa.

    T : float
        Temperature T / K.

    fluid : str
        Fluid name.

    Returns
    -------
    s : float
        Specific entropy s / (J/(kgK)).
    """
    if 'IDGAS::' in fluid:
        print('Ideal gas calculation not available by now.')
    elif 'TESPy::' in fluid:
        return tespy_fluid.fluids[fluid].funcs['s_pT'].ev(p, T)
    elif 'INCOMP::' in fluid:
        return CPPSI('S', 'P', p, 'T', T, fluid)
    else:
        memorise.heos[fluid].update(CP.PT_INPUTS, p, T)
        return memorise.heos[fluid].smass()


def ds_mix_pdT(flow, T):
    r"""
    Calculate partial derivate of entropy to temperature at constant pressure and fluid composition.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and fluid composition.

    T : float
        Temperature T / K.

    Returns
    -------
    ds / dT : float
        Partial derivative of specific entropy to temperature ds / dT / (J/(kg :math:`\mathrm{K}^2`)).

        .. math::

            \frac{\partial s_{mix}}{\partial T} =
            \frac{s_{mix}(p,T+d)-s_{mix}(p,T-d)}{2 \cdot d}
    """
    d = 2
    return (s_mix_pT(flow, T + d) - s_mix_pT(flow, T - d)) / (2 * d)

# %%


def molar_mass_flow(flow):
    r"""
    Calculates molar mass flow.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and fluid composition.

    Returns
    -------
    mm : float
        Molar mass flow mm / (mol/s).

        .. math::

            mm = \sum_{i} \left( \frac{x_{i}}{M_{i}} \right)
    """
    mm = 0
    for fluid, x in flow.items():
        if x > err:
            mm += x / molar_masses[fluid]
    return mm

# %%


def num_fluids(fluids):
    r"""
    Returns number of fluids in fluid mixture.

    Parameters
    ----------
    fluids : dict
        Fluid mass fractions.

    Returns
    -------
    n : int
        Number of fluids in fluid mixture n / 1.

        .. math::

            n = \sum_{i} \left( \begin{cases}
            0 & x_{i} < \epsilon \\
            1 & x_{i} \geq \epsilon
            \end{cases} \right)\;
            \forall i \in \text{network fluids}
    """
    n = 0
    for fluid, x in fluids.items():
        if x > err:
            n += 1

    return n

# %%


def single_fluid(fluids):
    r"""
    Returns the name of the pure fluid in a fluid vector.

    Parameters
    ----------
    fluids : dict
        Fluid mass fractions.

    Returns
    -------
    fluid : str
        Name of the single fluid.
    """
    if num_fluids(fluids) == 1:
        for fluid, x in fluids.items():
            if x > err:
                return fluid
    else:
        return []

# %%


def fluid_structure(fluid):
    r"""
    Returns the checmical formula of fluid.

    Parameters
    ----------
    fluid : str
        Name of the fluid.

    Returns
    -------
    parts : dict
        Dictionary of the chemical base elements as keys and the number of atoms in a molecule as values.
    """
    parts = {}
    for element in CP.CoolProp.get_fluid_param_string(fluid, 'formula').split('}'):
        if element != '':
            el = element.split('_{')
            parts[el[0]] = int(el[1])

    return parts

# %%


def lamb(re, ks, d):
    r"""
    Calculates the darcy friction factor from the moody diagram.

    Parameters
    ----------
    re : float
        Reynolds number re / 1.

    ks : float
        Pipe roughness ks / m.

    d : float
        Pipe diameter/characteristic lenght d / m.

    Returns
    -------
    lamb : float
        Darcy friction factor lamb / 1

    Note
    ----
    **Laminar flow** (:math:`re \leq 2320`)

    .. math::
        \lambda = \frac{64}{re}

    **turbulent flow** (:math:`re > 2320`)

    *hydraulically smooth:* :math:`\frac{re \cdot k_{s}}{d} < 65`

    .. math::
        \lambda = \begin{cases}
        0.03164 \cdot re^{-0.25} & re \leq 10^5\\
        0.0032 + 0.221 \cdot re^{-0.237} & 10^5 < re < 5 \cdot 10^6\\
        solve \left(0 = 2 \cdot \log\left(re \cdot \sqrt{\lambda} \right) -0.8
        - \frac{1}{\sqrt{\lambda}}\right) & re \geq 5 \cdot 10^6 \\
        \end{cases}

    *transition zone:* :math:`65 \leq \frac{re \cdot k_{s}}{d} \leq 1300`

    .. math::
        \lambda = solve \left( 0 = 2 \cdot \log \left( \frac{2.51}{re \cdot
        \sqrt{\lambda}} + \frac{k_{s}}{d} \cdot 0.269 \right) -
        \frac{1}{\sqrt{\lambda}} \right)

    *hydraulically rough:* :math:`\frac{re \cdot k_{s}}{d} > 1300`

    .. math::
        \lambda = \frac{1}{\left( 2\cdot \log \left( \frac{3.71 \cdot d}{k_{s}}
        \right) \right)}

    Example
    -------
    >>> from tespy import hlp
    >>> ks = 5e-5
    >>> d = 0.05
    >>> re_laminar = 2000
    >>> re_turb_smooth = 20000
    >>> re_turb_trans = 70000
    >>> re_high = 1000000
    >>> d_high = 0.8
    >>> re_very_high = 6000000
    >>> d_very_high = 1
    >>> ks_low = 1e-5
    >>> ks_rough = 1e-3
    >>> hlp.lamb(re_laminar, ks, d)
    0.032
    >>> round(hlp.lamb(re_turb_smooth, ks, d), 3)
    0.027
    >>> round(hlp.lamb(re_turb_trans, ks, d), 3)
    0.023
    >>> round(hlp.lamb(re_turb_trans, ks_rough, d), 3)
    0.049
    >>> round(hlp.lamb(re_high, ks, d_high), 3)
    0.012
    >>> round(hlp.lamb(re_very_high, ks_low, d_very_high), 3)
    0.009
    """
    if re <= 2320:
        return 64 / re
    else:
        if re * ks / d < 65:
            if re <= 1e5:
                return 0.3164 * re ** (-0.25)
            elif re > 1e5 and re < 5e6:
                return 0.0032 + 0.221 * re ** (-0.237)
            else:
                l0 = 0.0001
                return newton(lamb_smooth, dlamb_smooth_dl, [re], 0, val0=l0, valmin=0.00001, valmax=0.2)

        elif re * ks / d > 1300:
            return 1 / (2 * math.log(3.71 * d / ks, 10)) ** 2

        else:
            l0 = 0.002
            return newton(lamb_trans, dlamb_trans_dl, [re, ks, d], 0, val0=l0, valmin=0.0001, valmax=0.2)


def lamb_smooth(params, l):
    re = params[0]
    return 2 * math.log(re * math.sqrt(l), 10) - 0.8 - 1 / math.sqrt(l)


def dlamb_smooth_dl(params, l):
    return 1 / (l * math.log(10)) + 1 / 2 * l ** (-1.5)


def lamb_trans(params, l):
    re, ks, d = params[0], params[1], params[2]
    return (2 * math.log(2.51 / (re * math.sqrt(l)) + ks / d * 0.269, 10) + 1 / math.sqrt(l))


def dlamb_trans_dl(params, l):
    d = 0.001
    return (lamb_trans(params, l + d) - lamb_trans(params, l - d)) / (2 * d)


def modify_path_os(path):
    """
    Modify the path according the os. Also detects weather the path
    specification is absolute or relative and adjusts the path respectively.

    Parameters
    ----------
    path : str
        Path to modify.

    Returns
    -------
    path : str
        Modified path.
    """

    if os.name == 'nt':
        # windows
        path = path.replace('/', '\\')
        if path[0] != '\\' and path[1:2] != ':' and path[0] != '.':
            # relative path
            path = '.\\' + path
    elif os.name == 'posix':
        # linux, max
        path = path.replace('\\', '/')
        if path[0] != '/' and path[0] != '.':
            # absolute path
            path = './' + path
    else:
        # unkown os
        msg = 'Unknown operating system, using posix pathing logic.'
        logging.warning(msg)

    return path
