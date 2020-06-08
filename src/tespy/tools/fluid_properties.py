# -*- coding: utf-8

"""Module for fluid property integration.

TESPy uses the CoolProp python interface for all fluid property functions. The
tespy_fluid class allows the creation of lookup table for custom fluid
mixtures.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/tools/fluid_properties.py

SPDX-License-Identifier: MIT
"""

import logging
import os
from collections import OrderedDict

import CoolProp as CP
import numpy as np
import pandas as pd
from CoolProp.CoolProp import PropsSI as CPPSI
from scipy import interpolate

from tespy.tools.global_vars import err
from tespy.tools.global_vars import gas_constants
from tespy.tools.global_vars import molar_masses
from tespy.tools.helpers import molar_mass_flow
from tespy.tools.helpers import newton
from tespy.tools.helpers import reverse_2d
from tespy.tools.helpers import reverse_2d_deriv
from tespy.tools.helpers import single_fluid

# %%


class tespy_fluid:
    r"""
    Allow the creation of lookup tables for specified mixtures.

    The created fluid property lookup table adress an ideal mixture of real
    fluids.

    Parameters
    ----------
    alias : str
        Alias for the fluid. Please note: The alias of a tespy_fluid class
        object will always start with :code:`TESPy::`. See the example for more
        information!

    fluid : dict
        Fluid vector specifying the fluid composition of the TESPy fluid.

    p_range : list/ndarray
        Pressure range for the new fluid lookup table.

    T_range : list/ndarray
        Temperature range for the new fluid lookup table (optional).

    path : str
        Path to importing tespy fluid from.

    Note
    ----
    Creates lookup tables for

    - enthalpy (h),
    - entropy (s),
    - density (d) and
    - viscoity (visc)

    from pressure and temperature. Additionally molar mass and gas constant
    will be calculated. Inverse functions, e.g. entropy from pressure and
    enthalpy are calculated via newton algorithm from these tables.

    Example
    -------
    Create a custom fluid from specified composition within defined pressure
    and temperature limits. We define dry air component wise as our custom
    fluid.

    >>> from tespy.connections import connection
    >>> from tespy.components.basics import sink, source
    >>> from tespy.tools.fluid_properties import (tespy_fluid, h_mix_pT,
    ... s_mix_pT, v_mix_pT, visc_mix_pT)
    >>> from tespy.networks.networks import network
    >>> import shutil
    >>> fluidvec = {'N2': 0.7552, 'O2': 0.2314, 'CO2': 0.0005, 'Ar': 0.0129}
    >>> p_arr = np.array([0.1, 10]) * 1e5
    >>> myfluid = tespy_fluid('dry air', fluid=fluidvec, p_range=p_arr,
    ... T_range=[300, 1200])

    Check if the fluid creation was successful and compare some fluid
    properties to the CoolProp air implementation. We have to add the CoolProp
    air implementation to the memorise class first. The relative deviation
    should be very small (< 0.01), we check for enthalpy, volume, entropy and
    viscosity. Specific volume and viscosity are absolute values, thus no
    difference is calculated.

    >>> memorise.add_fluids({'air': 'HEOS'})

    >>> type(myfluid)
    <class 'tespy.tools.fluid_properties.tespy_fluid'>
    >>> p = 3e5
    >>> fluid_props = [0, p, 0, {myfluid.alias: 1}]
    >>> T1 = 400
    >>> T2 = 1000
    >>> delta_h_tespy = h_mix_pT(fluid_props, T2) - h_mix_pT(fluid_props, T1)
    >>> fluid_props_CP = [0, p, 0, {'air': 1}]
    >>> delta_h = h_mix_pT(fluid_props_CP, T2) - h_mix_pT(fluid_props_CP, T1)
    >>> round(abs(delta_h_tespy - delta_h) / delta_h, 2)
    0.0

    >>> v_tespy = v_mix_pT(fluid_props, T2)
    >>> v = v_mix_pT(fluid_props_CP, T2)
    >>> round(abs(v_tespy - v) / v, 2)
    0.0

    >>> s_tespy = s_mix_pT(fluid_props, T2) - s_mix_pT(fluid_props, T1)
    >>> s = s_mix_pT(fluid_props_CP, T2) - s_mix_pT(fluid_props_CP, T1)
    >>> round(abs(s_tespy - s) / s, 2)
    0.0

    >>> visc_tespy = visc_mix_pT(fluid_props, T2)
    >>> visc = visc_mix_pT(fluid_props_CP, T2)
    >>> round(abs(visc_tespy - visc) / visc, 2)
    0.0

    The fluid had been saved automatically, load it now.

    >>> loadfluid = tespy_fluid('dry air', fluid=fluidvec, p_range=p_arr,
    ... path='./LUT')
    >>> type(loadfluid)
    <class 'tespy.tools.fluid_properties.tespy_fluid'>
    >>> shutil.rmtree('./LUT', ignore_errors=True)
    """

    def __init__(self, alias, fluid, p_range, T_range=None, path=None):

        if not isinstance(alias, str):
            msg = 'Alias must be of type String.'
            logging.error(msg)
            raise TypeError(msg)

        # process parameters
        if '::' in alias:
            msg = 'The sequence "::" is not allowed in your fluid alias.'
            logging.error(msg)
            raise ValueError(msg)
        else:
            self.alias = alias

        self.fluid = fluid

        # path for loading
        self.path = path

        # create look up tables
        tespy_fluid.fluids[self.alias] = self

        self.fluids_back_ends = OrderedDict()

        for fluid in sorted(list(self.fluid.keys()) + [self.alias]):
            try:
                data = fluid.split('::')
                back_end = data[0]
            except IndexError:
                back_end = 'HEOS'
            if fluid != self.alias:
                self.fluids_back_ends[fluid] = back_end

        memorise.add_fluids(self.fluids_back_ends)

        # set up grid
        self.p = np.geomspace(p_range[0], p_range[1], 100)

        if T_range is None:
            T_range = [max([memorise.value_range[f][2] for f
                            in self.fluids_back_ends.keys()]) + 1, 2000]

        self.T = np.geomspace(T_range[0], T_range[1], 100)

        memorise.add_fluids({self.alias: 'TESPy'})

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
            msg = ('Generating lookup-tables from base path ' + self.path +
                   '/' + self.alias + '/.')
            logging.debug(msg)
            for key in params.keys():
                self.funcs[key] = self.load_lookup(key)
                msg = 'Loading function values for function ' + key + '.'
                logging.debug(msg)

        msg = ('Successfully created look-up-tables for custom fluid ' +
               self.alias + '.')
        logging.debug(msg)

    def generate_lookup(self, name, func):
        r"""
        Create lookup table from CoolProp-database.

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

        self.save_lookup(name, x1, x2, y)

        func = interpolate.RectBivariateSpline(x1, x2, y)
        return func

    def save_lookup(self, name, x1, x2, y):
        r"""
        Save lookup table to working dir :code:`./LUT/fluid_alias/`.

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
        alias = self.alias.replace('::', '_')
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
        alias = self.alias.replace('::', '_')
        path = self.path + '/' + alias + '/' + name + '.csv'
        df = pd.read_csv(path, index_col=0)

        x1 = df.index.values
        x2 = np.array(list(map(float, list(df))))
        y = df.values

        func = interpolate.RectBivariateSpline(x1, x2, y)
        return func


# create dict for tespy fluids
tespy_fluid.fluids = {}

# %%


class memorise:
    r"""Memorization of fluid properties."""

    @staticmethod
    def add_fluids(fluids):
        r"""
        Add list of fluids to fluid memorisation class.

        - Generate arrays for fluid property lookup.
        - Calculate/set fluid property value ranges for convergence checks.

        Parameters
        ----------
        fluids : dict
            Dict of fluid and corresponding CoolProp back end for fluid
            property memorization.

        Note
        ----
        The memorise class creates globally accessible variables for different
        fluid property calls as dictionaries:

        - T(p,h)
        - T(p,s)
        - v(p,h)
        - visc(p,h)
        - s(p,h)

        Each dictionary uses the list of fluids passed to the memorise class as
        identifier for the fluid property memorisation. The fluid properties
        are stored as numpy array, where each column represents the mass
        fraction of the respective fluid and the additional columns are the
        values for the fluid properties. The fluid property function will then
        look for identical fluid property inputs (p, h, (s), fluid mass
        fraction). If the inputs are in the array, the first column of that row
        is returned, see example.

        Example
        -------
        T(p,h) for set of fluids ('water', 'air'):

        - row 1: [282.64527752319697, 10000, 40000, 1, 0]
        - row 2: [284.3140698256616, 10000, 47000, 1, 0]
        """
        # number of fluids
        num_fl = len(fluids)
        if num_fl > 0:
            fl = tuple(fluids.keys())
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

            msg = 'Added fluids ' + str(fl) + ' to memorise lookup tables.'
            logging.debug(msg)

        for f, back_end in fluids.items():
            if back_end == 'IF97':
                msg = (
                    'Due to a bug in the IF97 CoolProp back end, it is not '
                    'possible to use this back end at the moment. For more '
                    'information see '
                    'https://github.com/CoolProp/CoolProp/issues/1918.')
                logging.error(msg)
                raise ValueError(msg)
            elif f not in memorise.state.keys() and back_end != 'TESPy':
                # create CoolProp.AbstractState object
                try:
                    memorise.state[f] = CP.AbstractState(back_end, f)
                except ValueError:
                    msg = (
                        'Could not find the fluid "' + f + '" in the fluid '
                        'property database. If you are using a stoichimetric '
                        'combustion chamber, this warning can be ignored in '
                        'case the fluid is your fuel, your flue gas or your '
                        'air.')
                    logging.warning(msg)
                    continue

                msg = (
                    'Created CoolProp.AbstractState object for fluid ' +
                    f + ' with back end ' + back_end + '.')
                logging.debug(msg)
                # pressure range
                try:
                    pmin = memorise.state[f].trivial_keyed_output(CP.iP_min)
                    pmax = memorise.state[f].trivial_keyed_output(CP.iP_max)
                except ValueError:
                    pmin = 1e4
                    pmax = 1e8
                    msg = (
                        'Could not find values for maximum and minimum '
                        'pressure.')
                    logging.warning(msg)

                # temperature range
                Tmin = memorise.state[f].trivial_keyed_output(CP.iT_min)
                Tmax = memorise.state[f].trivial_keyed_output(CP.iT_max)

                # value range for fluid properties
                memorise.value_range[f] = [pmin, pmax, Tmin, Tmax]

                try:
                    molar_masses[f] = memorise.state[f].molar_mass()
                    gas_constants[f] = memorise.state[f].gas_constant()
                except ValueError:
                    try:
                        molar_masses[f] = CPPSI('M', f)
                        gas_constants[f] = CPPSI('GAS_CONSTANT', f)
                    except ValueError:
                        molar_masses[f] = 1
                        gas_constants[f] = 1
                        msg = (
                            'Could not find values for molar mass and gas '
                            'constant.')
                        logging.warning(msg)

                msg = (
                    'Specifying fluid property ranges for pressure and '
                    'temperature for convergence check of fluid ' + f + '.')
                logging.debug(msg)

            elif back_end == 'TESPy':
                pmin = tespy_fluid.fluids[f].p[0]
                pmax = tespy_fluid.fluids[f].p[-1]
                Tmin = tespy_fluid.fluids[f].T[0]
                Tmax = tespy_fluid.fluids[f].T[-1]
                msg = (
                    'Loading fluid property ranges for TESPy-fluid ' + f + '.')
                logging.debug(msg)
                # value range for fluid properties
                memorise.value_range[f] = [pmin, pmax, Tmin, Tmax]
                molar_masses[f] = 1 / molar_mass_flow(
                    tespy_fluid.fluids[f].fluid)
                gas_constants[f] = (
                    gas_constants['uni'] / molar_masses[f])
                msg = ('Specifying fluid property ranges for pressure and '
                       'temperature for convergence check.')
                logging.debug(msg)

    @staticmethod
    def del_memory(fluids):
        r"""
        Delete non frequently used fluid property values from memorise class.

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

        msg = ('Dropping not frequently used fluid property values from '
               'memorise class for fluids ' + str(fl) + '.')
        logging.debug(msg)


# create memorise dictionaries
memorise.state = {}
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
memorise.value_range = {}

# %%


def T_mix_ph(flow, T0=300):
    r"""
    Calculate the temperature from pressure and enthalpy.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

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
    memorisation = fl in memorise.T_ph.keys()
    if memorisation is True:
        a = memorise.T_ph[fl][:, :-1]
        b = np.array([flow[1], flow[2]] + list(flow[3].values()))
        ix = np.where(np.all(abs(a - b) <= err, axis=1))[0]

        if ix.size == 1:
            # known fluid properties
            T = memorise.T_ph[fl][ix, -1][0]
            memorise.T_ph_f[fl] += [T]
            return T

    # unknown fluid properties
    fluid = single_fluid(flow[3])
    if fluid is None:
        # calculate the fluid properties for fluid mixtures
        if memorisation is True:
            valmin = max(
                [memorise.value_range[f][2] for f in fl if flow[3][f] > err]
            ) + 0.1
            if T0 < valmin or np.isnan(T0):
                T0 = valmin * 1.1
        else:
            valmin = 70
        val = newton(h_mix_pT, dh_mix_pdT, flow, flow[2], val0=T0,
                     valmin=valmin, valmax=3000, imax=10)
    else:
        # calculate fluid property for pure fluids
        val = T_ph(flow[1], flow[2], fluid)

    if memorisation is True:
        # memorise the newly calculated value
        new = np.asarray([[flow[1], flow[2]] + list(flow[3].values()) + [val]])
        memorise.T_ph[fl] = np.append(memorise.T_ph[fl], new, axis=0)

    return val


def T_ph(p, h, fluid):
    r"""
    Calculate the temperature from pressure and enthalpy for a pure fluid.

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
    if fluid in tespy_fluid.fluids.keys():
        db = tespy_fluid.fluids[fluid].funcs['h_pT']
        return newton(reverse_2d, reverse_2d_deriv, [db, p, h], 0)
    else:
        memorise.state[fluid].update(CP.HmassP_INPUTS, h, p)
        return memorise.state[fluid].T()


def dT_mix_dph(flow, T0=300):
    r"""
    Calculate partial derivate of temperature to pressure.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

    Returns
    -------
    dT / dp : float
        Partial derivative of temperature to pressure dT /dp / (K/Pa).

        .. math::

            \frac{\partial T_{mix}}{\partial p} = \frac{T_{mix}(p+d,h)-
            T_{mix}(p-d,h)}{2 \cdot d}
    """
    d = 1
    up = flow.copy()
    lo = flow.copy()
    up[1] += d
    lo[1] -= d
    return (T_mix_ph(up, T0=T0) - T_mix_ph(lo, T0=T0)) / (2 * d)


def dT_mix_pdh(flow, T0=300):
    r"""
    Calculate partial derivate of temperature to enthalpy.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

    Returns
    -------
    dT / dh : float
        Partial derivative of temperature to enthalpy dT /dh / ((kgK)/J).

        .. math::

            \frac{\partial T_{mix}}{\partial h} = \frac{T_{mix}(p,h+d)-
            T_{mix}(p,h-d)}{2 \cdot d}
    """
    d = 1
    up = flow.copy()
    lo = flow.copy()
    up[2] += d
    lo[2] -= d
    return (T_mix_ph(up, T0=T0) - T_mix_ph(lo, T0=T0)) / (2 * d)


def dT_mix_ph_dfluid(flow, T0=300):
    r"""
    Calculate partial derivate of temperature to fluid composition.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

    Returns
    -------
    dT / dfluid : ndarray
        Partial derivatives of temperature to fluid composition
        dT / dfluid / K.

        .. math::

            \frac{\partial T_{mix}}{\partial fluid_{i}} =
            \frac{T_{mix}(p,h,fluid_{i}+d)-
            T_{mix}(p,h,fluid_{i}-d)}{2 \cdot d}
    """
    d = 1e-5
    up = flow.copy()
    lo = flow.copy()
    vec_deriv = []
    for fluid, x in flow[3].items():
        if x > err:
            up[3][fluid] += d
            lo[3][fluid] -= d
            vec_deriv += [
                (T_mix_ph(up, T0=T0) - T_mix_ph(lo, T0=T0)) / (2 * d)]
            up[3][fluid] -= d
            lo[3][fluid] += d
        else:
            vec_deriv += [0]

    return np.asarray(vec_deriv)

# %%


def T_mix_ps(flow, s, T0=300):
    r"""
    Calculate the temperature from pressure and entropy.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

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
    memorisation = fl in memorise.T_ps.keys()
    if memorisation is True:
        a = memorise.T_ps[fl][:, :-1]
        b = np.asarray([flow[1], flow[2]] + list(flow[3].values()) + [s])
        ix = np.where(np.all(abs(a - b) <= err, axis=1))[0]
        if ix.size == 1:
            # known fluid properties
            T = memorise.T_ps[fl][ix, -1][0]
            memorise.T_ps_f[fl] += [T]
            return T

    # unknown fluid properties
    fluid = single_fluid(flow[3])
    if fluid is None:
        # calculate the fluid properties for fluid mixtures
        if memorisation is True:
            valmin = max(
                [memorise.value_range[f][2] for f in fl if flow[3][f] > err]
            ) + 0.1
            if T0 < valmin or np.isnan(T0):
                T0 = valmin * 1.1
        else:
            valmin = 70

        val = newton(s_mix_pT, ds_mix_pdT, flow, s, val0=T0,
                     valmin=valmin, valmax=3000, imax=10)

    else:
        # calculate fluid property for pure fluids
        val = T_ps(flow[1], s, fluid)

    if memorisation is True:
        new = np.asarray(
            [[flow[1], flow[2]] + list(flow[3].values()) + [s, val]])
        # memorise the newly calculated value
        memorise.T_ps[fl] = np.append(memorise.T_ps[fl], new, axis=0)

    return val


def T_ps(p, s, fluid):
    r"""
    Calculate the temperature from pressure and entropy for a pure fluid.

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
    if fluid in tespy_fluid.fluids.keys():
        db = tespy_fluid.fluids[fluid].funcs['s_pT']
        return newton(reverse_2d, reverse_2d_deriv, [db, p, s], 0)
    else:
        memorise.state[fluid].update(CP.PSmass_INPUTS, p, s)
        return memorise.state[fluid].T()

# %%


def h_mix_pT(flow, T):
    r"""
    Calculate the enthalpy from pressure and Temperature.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

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
    Calculate the enthalpy from pressure and temperature for a pure fluid.

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
    if fluid in tespy_fluid.fluids.keys():
        return tespy_fluid.fluids[fluid].funcs['h_pT'].ev(p, T)
    else:
        memorise.state[fluid].update(CP.PT_INPUTS, p, T)
        return memorise.state[fluid].hmass()


def dh_mix_pdT(flow, T):
    r"""
    Calculate partial derivate of enthalpy to temperature.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

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
    d = 0.1
    return (h_mix_pT(flow, T + d) - h_mix_pT(flow, T - d)) / (2 * d)

# %%


def h_mix_ps(flow, s, T0=300):
    r"""
    Calculate the enthalpy from pressure and temperature.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

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
    return h_mix_pT(flow, T_mix_ps(flow, s, T0=T0))


def h_ps(p, s, fluid):
    r"""
    Calculate the enthalpy from pressure and entropy for a pure fluid.

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
    if fluid in tespy_fluid.fluids.keys():
        db = tespy_fluid.fluids[fluid].funcs['s_pT']
        T = newton(reverse_2d, reverse_2d_deriv, [db, p, s], 0)
        return tespy_fluid.fluids[fluid].funcs['h_pT'].ev(p, T)
    else:
        memorise.state[fluid].update(CP.PSmass_INPUTS, p, s)
        return memorise.state[fluid].hmass()

# %%


def h_mix_pQ(flow, Q):
    r"""
    Calculate the enthalpy from pressure and vapour mass fraction.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

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
    fluid = single_fluid(flow[3])
    if fluid is None:
        msg = 'The function h_mix_pQ can only be used for pure fluids.'
        logging.error(msg)
        raise ValueError(msg)

    try:
        memorise.state[fluid].update(CP.PQ_INPUTS, flow[1], Q)
    except ValueError:
        pcrit = memorise.state[fluid].trivial_keyed_output(CP.iP_critical)
        memorise.state[fluid].update(CP.PQ_INPUTS, pcrit * 0.99, Q)

    return memorise.state[fluid].hmass()


def dh_mix_dpQ(flow, Q):
    r"""
    Calculate partial derivate of enthalpy to vapour mass fraction.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

    Q : float
        Vapour mass fraction Q / 1.

    Returns
    -------
    dh / dQ : float
        Partial derivative of enthalpy to vapour mass fraction
        dh / dQ / (J/kg).

        .. math::

            \frac{\partial h_{mix}}{\partial p} =
            \frac{h_{mix}(p+d,Q)-h_{mix}(p-d,Q)}{2 \cdot d}\\
            Q: \text{vapour mass fraction}

    Note
    ----
    This works for pure fluids only!
    """
    d = 1
    up = flow.copy()
    lo = flow.copy()
    up[1] += d
    lo[1] -= d
    return (h_mix_pQ(up, Q) - h_mix_pQ(lo, Q)) / (2 * d)

# %%


def T_bp_p(flow):
    r"""
    Calculate temperature from boiling point pressure.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

    Returns
    -------
    T : float
        Temperature at boiling point.

    Note
    ----
    This function works for pure fluids only!
    """
    fluid = single_fluid(flow[3])
    pcrit = memorise.state[fluid].trivial_keyed_output(CP.iP_critical)
    if flow[1] > pcrit:
        memorise.state[fluid].update(CP.PQ_INPUTS, pcrit * 0.99, 1)
    else:
        memorise.state[fluid].update(CP.PQ_INPUTS, flow[1], 1)
    return memorise.state[fluid].T()


def dT_bp_dp(flow):
    r"""
    Calculate partial derivate of temperature to boiling point pressure.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

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
    up = flow.copy()
    lo = flow.copy()
    up[1] += d
    lo[1] -= d
    return (T_bp_p(up) - T_bp_p(lo)) / (2 * d)

# %%


def v_mix_ph(flow, T0=300):
    r"""
    Calculate the specific volume from pressure and enthalpy.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

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
    memorisation = fl in memorise.v_ph.keys()
    if memorisation is True:
        a = memorise.v_ph[fl][:, :-1]
        b = np.asarray([flow[1], flow[2]] + list(flow[3].values()))
        ix = np.where(np.all(abs(a - b) <= err, axis=1))[0]
        if ix.size == 1:
            # known fluid properties
            v = memorise.v_ph[fl][ix, -1][0]
            memorise.v_ph_f[fl] += [v]
            return v

    # unknown fluid properties
    fluid = single_fluid(flow[3])
    if fluid is None:
        # calculate the fluid properties for fluid mixtures
        val = v_mix_pT(flow, T_mix_ph(flow, T0=T0))
    else:
        # calculate fluid property for pure fluids
        val = 1 / d_ph(flow[1], flow[2], fluid)

    if memorisation is True:
        # memorise the newly calculated value
        new = np.asarray([[flow[1], flow[2]] + list(flow[3].values()) + [val]])
        memorise.v_ph[fl] = np.append(memorise.v_ph[fl], new, axis=0)

    return val


def d_ph(p, h, fluid):
    r"""
    Calculate the density from pressure and enthalpy for a pure fluid.

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
    if fluid in tespy_fluid.fluids.keys():
        db = tespy_fluid.fluids[fluid].funcs['h_pT']
        T = newton(reverse_2d, reverse_2d_deriv, [db, p, h], 0)
        return tespy_fluid.fluids[fluid].funcs['d_pT'].ev(p, T)
    else:
        memorise.state[fluid].update(CP.HmassP_INPUTS, h, p)
        return memorise.state[fluid].rhomass()

# %%


def Q_ph(p, h, fluid):
    r"""
    Calculate vapor mass fraction from pressure and enthalpy for a pure fluid.

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
    x : float
        Vapor mass fraction.
    """
    try:
        memorise.state[fluid].update(CP.HmassP_INPUTS, h, p)
        return memorise.state[fluid].Q()
    except (KeyError, ValueError):
        return np.nan

# %%


def dv_mix_dph(flow, T0=300):
    r"""
    Calculate partial derivate of specific volume to pressure.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

    Returns
    -------
    dv / dp : float
        Partial derivative of specific volume to pressure
        dv /dp / (:math:`\mathrm{m}^3`/(Pa kg)).

        .. math::

            \frac{\partial v_{mix}}{\partial p} = \frac{v_{mix}(p+d,h)-
            v_{mix}(p-d,h)}{2 \cdot d}
    """
    d = 1
    up = flow.copy()
    lo = flow.copy()
    up[1] += d
    lo[1] -= d
    return (v_mix_ph(up, T0=T0) - v_mix_ph(lo, T0=T0)) / (2 * d)


def dv_mix_pdh(flow, T0=300):
    r"""
    Calculate partial derivate of specific volume to enthalpy.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

    Returns
    -------
    dv / dh : float
        Partial derivative of specific volume to enthalpy
        dv /dh / (:math:`\mathrm{m}^3`/J).

        .. math::

            \frac{\partial v_{mix}}{\partial h} = \frac{v_{mix}(p,h+d)-
            v_{mix}(p,h-d)}{2 \cdot d}
    """
    d = 1
    up = flow.copy()
    lo = flow.copy()
    up[2] += d
    lo[2] -= d
    return (v_mix_ph(up, T0=T0) - v_mix_ph(lo, T0=T0)) / (2 * d)

# %%


def v_mix_pT(flow, T):
    r"""
    Calculate the specific volume from pressure and temperature.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

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

        v_{mix}(p,T)=\frac{1}{\sum_{i} \rho(pp_{i}, T, fluid_{i})}\;
        \forall i \in \text{fluid components}\\
        pp: \text{partial pressure}
    """
    n = molar_mass_flow(flow[3])

    d = 0
    for fluid, x in flow[3].items():
        if x > err:
            ni = x / molar_masses[fluid]
            d += d_pT(flow[1] * ni / n, T, fluid)

    return 1 / d


def d_mix_pT(flow, T):
    r"""
    Calculate the density from pressure and temperature.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

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
    Calculate the density from pressure and temperature for a pure fluid.

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
    if fluid in tespy_fluid.fluids.keys():
        return tespy_fluid.fluids[fluid].funcs['d_pT'].ev(p, T)
    else:
        memorise.state[fluid].update(CP.PT_INPUTS, p, T)
        return memorise.state[fluid].rhomass()

# %%


def visc_mix_ph(flow, T0=300):
    r"""
    Calculate the dynamic viscorsity from pressure and enthalpy.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

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
    memorisation = fl in memorise.visc_ph.keys()
    if memorisation is True:
        a = memorise.visc_ph[fl][:, :-1]
        b = np.asarray([flow[1], flow[2]] + list(flow[3].values()))
        ix = np.where(np.all(abs(a - b) <= err, axis=1))[0]
        if ix.size == 1:
            # known fluid properties
            visc = memorise.visc_ph[fl][ix, -1][0]
            memorise.visc_ph_f[fl] += [visc]
            return visc

    # unknown fluid properties
    fluid = single_fluid(flow[3])
    if fluid is None:
        # calculate the fluid properties for fluid mixtures
        val = visc_mix_pT(flow, T_mix_ph(flow, T0=T0))
    else:
        # calculate the fluid properties for pure fluids
        val = visc_ph(flow[1], flow[2], fluid)

    if memorisation is True:
        # memorise the newly calculated value
        new = np.asarray([[flow[1], flow[2]] + list(flow[3].values()) + [val]])
        memorise.visc_ph[fl] = np.append(memorise.visc_ph[fl], new, axis=0)
    return val


def visc_ph(p, h, fluid):
    r"""
    Calculate dynamic viscosity from pressure and enthalpy for a pure fluid.

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
    if fluid in tespy_fluid.fluids.keys():
        db = tespy_fluid.fluids[fluid].funcs['h_pT']
        T = newton(reverse_2d, reverse_2d_deriv, [db, p, h], 0)
        return tespy_fluid.fluids[fluid].funcs['visc_pT'].ev(p, T)
    else:
        memorise.state[fluid].update(CP.HmassP_INPUTS, h, p)
        return memorise.state[fluid].viscosity()

# %%


def visc_mix_pT(flow, T):
    r"""
    Calculate dynamic viscosity from pressure and temperature.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

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

    Reference: :cite:`Herning1936`.
    """
    n = molar_mass_flow(flow[3])

    a = 0
    b = 0
    for fluid, x in flow[3].items():
        if x > err:
            bi = x * np.sqrt(molar_masses[fluid]) / (molar_masses[fluid] * n)
            b += bi
            a += bi * visc_pT(flow[1], T, fluid)

    return a / b


def visc_pT(p, T, fluid):
    r"""
    Calculate dynamic viscosity from pressure and temperature for a pure fluid.

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
    if fluid in tespy_fluid.fluids.keys():
        return tespy_fluid.fluids[fluid].funcs['visc_pT'].ev(p, T)
    else:
        memorise.state[fluid].update(CP.PT_INPUTS, p, T)
        return memorise.state[fluid].viscosity()

# %%


def s_mix_ph(flow, T0=300):
    r"""
    Calculate the entropy from pressure and enthalpy.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

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
    memorisation = fl in memorise.s_ph.keys()
    if memorisation is True:
        a = memorise.s_ph[fl][:, :-1]
        b = np.asarray([flow[1], flow[2]] + list(flow[3].values()))
        ix = np.where(np.all(abs(a - b) <= err, axis=1))[0]
        if ix.size == 1:
            # known fluid properties
            s = memorise.s_ph[fl][ix, -1][0]
            memorise.s_ph_f[fl] += [s]
            return s

    # unknown fluid properties
    fluid = single_fluid(flow[3])
    if fluid is None:
        # calculate the fluid properties for fluid mixtures
        val = s_mix_pT(flow, T_mix_ph(flow, T0=T0))
    else:
        # calculate fluid property for pure fluids
        val = s_ph(flow[1], flow[2], fluid)

    if memorisation is True:
        # memorise the newly calculated value
        new = np.asarray([[flow[1], flow[2]] + list(flow[3].values()) + [val]])
        memorise.s_ph[fl] = np.append(memorise.s_ph[fl], new, axis=0)

    return val


def s_ph(p, h, fluid):
    r"""
    Calculate the entropy from pressure and enthalpy for a pure fluid.

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
    if fluid in tespy_fluid.fluids.keys():
        db = tespy_fluid.fluids[fluid].funcs['h_pT']
        T = newton(reverse_2d, reverse_2d_deriv, [db, p, h], 0)
        return tespy_fluid.fluids[fluid].funcs['s_pT'].ev(p, T)
    else:
        memorise.state[fluid].update(CP.HmassP_INPUTS, h, p)
        return memorise.state[fluid].smass()

# %%


def s_mix_pT(flow, T):
    r"""
    Calculate the entropy from pressure and temperature.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

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

    tespy_fluids = [
        s for s in flow[3].keys() if s in tespy_fluid.fluids.keys()]
    for f in tespy_fluids:
        for it, x in tespy_fluid.fluids[f].fluid.items():
            if it in fluid_comps.keys():
                fluid_comps[it] += x * flow[3][f]
            else:
                fluid_comps[it] = x * flow[3][f]

    fluids = [s for s in flow[3].keys() if s not in tespy_fluid.fluids.keys()]
    for f in fluids:
        if f in fluid_comps.keys():
            fluid_comps[f] += flow[3][f]
        else:
            fluid_comps[f] = flow[3][f]

    # print(fluid_comps)

    s = 0
    for fluid, x in fluid_comps.items():
        if x > err:
            pp = flow[1] * x / (molar_masses[fluid] * n)
            s += s_pT(pp, T, fluid) * x
            s -= (x * gas_constants[fluid] / molar_masses[fluid] *
                  np.log(pp / flow[1]))

    return s


def s_pT(p, T, fluid):
    r"""
    Calculate the entropy from pressure and temperature for a pure fluid.

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
    if fluid in tespy_fluid.fluids.keys():
        return tespy_fluid.fluids[fluid].funcs['s_pT'].ev(p, T)
    else:
        memorise.state[fluid].update(CP.PT_INPUTS, p, T)
        return memorise.state[fluid].smass()


def ds_mix_pdT(flow, T):
    r"""
    Calculate partial derivate of entropy to temperature.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

    T : float
        Temperature T / K.

    Returns
    -------
    ds / dT : float
        Partial derivative of specific entropy to temperature
        ds / dT / (J/(kg :math:`\mathrm{K}^2`)).

        .. math::

            \frac{\partial s_{mix}}{\partial T} =
            \frac{s_{mix}(p,T+d)-s_{mix}(p,T-d)}{2 \cdot d}
    """
    d = 0.1
    return (s_mix_pT(flow, T + d) - s_mix_pT(flow, T - d)) / (2 * d)
