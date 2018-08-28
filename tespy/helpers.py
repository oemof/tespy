"""
.. module:: helpers
    :synopsis: helpers for frequently used functionalities

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>
"""

import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI as CPPSI

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import math
import numpy as np
import sys
from scipy import interpolate
import pandas as pd
import os

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
    """r

    The data container stores data on components and connections attributes.
    There are subclasses for the following applications:

    - mass flow, pressure, enthalpy and temperature
    - fluid
    - component parameters
    - component characteristics

    **allowed keywords** in kwargs:

    - see data_container.attr()
    """

    def __init__(self, **kwargs):

        invalid = []
        var = self.attr()

        # default values
        for key in var.keys():
            self.__dict__.update({key: var[key]})

        # specify values
        for key in kwargs:
            if key not in var.keys():
                invalid += []
            self.__dict__.update({key: kwargs[key]})

        # print invalid keywords
        if len(invalid) > 0:
            print('The following keys are not available: ' + str(invalid))

    def set_attr(self, **kwargs):

        invalid = []
        var = self.attr()

        # specify values
        for key in kwargs:
            if key not in var.keys():
                invalid += []
            self.__dict__.update({key: kwargs[key]})

        # print invalid keywords
        if len(invalid) > 0:
            print('The following keys are not available: ' + str(invalid))

    def get_attr(self, key):
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            print('No attribute \"', key, '\" available!')
            return None

    def attr(self):
        return {}


class dc_prop(data_container):
    """r

    data container for fluid properties

    **value specification**

    - val (*numeric*) - user specified value
    - val0 (*numeric*) - user specified starting value
    - val_SI (*numeric*) - value in SI unit
    - val_set (*bool*) - is the specified value a parameter?

    **reference specification**

    - ref (*numeric*) - referenced connection
    - ref_set (*bool*) - is the reference a parameter?

    **units**

    - unit (*str*) - unit
    - unit_set (*bool*) - is the unit set for the corresponding value? if not,
      network unit will be used in calculation (default)
    """
    def attr(self):
        return {'val': np.nan, 'val0': np.nan, 'val_SI': 0, 'val_set': False,
                'ref': None, 'ref_set': False,
                'unit': None, 'unit_set': False}


class dc_flu(data_container):
    """r

    data container for fluid vector

    - val (*dict*) - user specified values
    - val0 (*dict*) - user specified starting values
    - val_set (*dict*) - which components of the fluid vector are set?
    - balance (*bool*) - apply fluid balance equation?
    """
    def attr(self):
        return {'val': {}, 'val0': {}, 'val_set': {}, 'balance': False}


class dc_cp(data_container):
    """r

    data container for component properties

    - val (*numeric*) - user specified value
    - val_SI (*numeric*) - value in SI units
    - val_set (*bool*) - is the specified value set?
    - is_var (*bool*) - make this parameter a variable of the system? if so,
      val will be used as starting value
    """
    def attr(self):
        return {'val': 0, 'val_SI': 0, 'is_set': False, 'is_var': False}


class dc_cc(data_container):
    """r

    data container for component characteristics

    - func (*tespy.components.characteristics.characteristics object*) -
      characteristic function to be applied
    - func_set (*bool*) - is the characteristic function set?

    **using default characteristics**

    see tespy.components.characteristics module for default methods and
    parameters, also see tespy.components.components module for available
    parameters.

    - method (*str*) - which method of the characteristic function should be
      applied?
    - param (*str*) - to which parameter should the characteristic function be
      applied?

    **using custom characteristics**

    linear interpolation will be applied, it is possible to use default
    characteristics and overwrite x-values or y-values

    - x (*np.array*) - array for the x-values of the characteristic line
    - y (*np.array*) - array for the y-values of the characteristic line

    """
    def attr(self):
        return {'func': None, 'is_set': False,
                'method': 'default', 'param': None,
                'x': None, 'y': None}


class dc_gcp(data_container):
    """r

    data container for grouped component properties

    - is_set (*bool*) - is the group set
    - method (*str*) - calculation method for identical property groups
    - elements (*list*) - list of elements for this group, if you want to make
      use of a group, every element's is_set value must be True.

    """
    def attr(self):
        return {'is_set': False, 'method': 'default', 'elements': []}

# %%


class MyNetworkError(Exception):
    pass


class MyConnectionError(Exception):
    pass


class MyComponentError(Exception):
    pass


class MyConvergenceError(Exception):
    pass


def query_yes_no(question, default='yes'):
    """
    in prompt query

    :param question: question to ask in prompt
    :type question: str
    :param default: default answer
    :type default: str
    :returns: bool
    """
    valid = {'yes': True,
             'y': True,
             'ye': True,
             'no': False,
             'n': False}
    if default is None:
        prompt = '[y / n]'
    elif default == 'yes':
        prompt = '[Y / n]'
    elif default == 'no':
        prompt = '[y / N]'

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write('Please respond with \'yes\' or \'no\' '
                             '(or \'y\' or \'n\').\n')

# %%


class tespy_fluid:
    """r

    The tespy_fluid class allows the creation of custom fluid properies for a
    specified mixture of fluids. The created fluid properties adress an ideal
    mixture of real fluids.

    Creates lookup-tables for

    - enthalpy,
    - entropy,
    - density and
    - viscoity

    from pressure and temperature. Additionally molar mass and gas constant
    will be calculated. Inverse functions, e. g. entropy from pressure and
    enthalpy are calculated via newton algorithm from these tables.

    :param alias: name of the fluid mixture will be "TESPy::alias"
    :type alias: str
    :param fluid: fluid vector for composition {fluid_i: mass fraction, ...}
    :type fluid: dict
    :param p_range: range of feasible pressures for newly created fluid
                    (provide in SI units)
    :type p_range: list
    :param T_range: range of feasible temperatures for newly created fluid
                    (provide in SI units)
    :type T_range: list
    :returns: no return value
    :raises: - :code:`TypeError`, if alias is not of type string
             - :code:`ValueError`, if the alias contains "IDGAS::"

    **allowed keywords** in kwargs:

    - plot (*bool*), plot the lookup table after creation
    """

    def __init__(self, alias, fluid, p_range, T_range, **kwargs):

        if not hasattr(tespy_fluid, 'fluids'):
            tespy_fluid.fluids = {}

        if not isinstance(alias, str):
            msg = 'Alias must be of type String.'
            raise TypeError(msg)

        if 'IDGAS::' in alias:
            msg = 'You are not allowed to use "IDGAS::" within your alias.'
            raise ValueError(msg)

        # process parameters
        if 'TESPy::' in alias:
            self.alias = alias
        else:
            self.alias = 'TESPy::' + alias
        self.fluid = fluid

        # load LUT from this path
        self.path = kwargs.get('path', dc_cp())

        # adjust value ranges according to specified unit system
        self.p_range = np.array(p_range)
        self.T_range = np.array(T_range)

        # set up grid
        self.p = np.linspace(self.p_range[0], self.p_range[1])
        self.T = np.linspace(self.T_range[0], self.T_range[1])

        # plotting
        self.plot = kwargs.get('plot', False)

        # calculate molar mass and gas constant
        for f in self.fluid:
            molar_masses[f] = CPPSI('M', f)
            gas_constants[f] = CPPSI('GAS_CONSTANT', f)

        molar_masses[self.alias] = 1 / molar_massflow(self.fluid)
        gas_constants[self.alias] = (gas_constants['uni'] /
                                     molar_masses[self.alias])

        # create look up tables
        tespy_fluid.fluids[self.alias] = {}

        params = {}

        params['h_pT'] = h_mix_pT
        params['s_pT'] = s_mix_pT
        params['d_pT'] = d_mix_pT
        params['visc_pT'] = visc_mix_pT

        self.funcs = {}

        if not self.path.is_set:

            for key in params.keys():

                self.funcs[key] = self.generate_lookup(key, params[key])

        else:

            for key in params.keys():
                self.funcs[key] = self.load_lookup(key)

        tespy_fluid.fluids[self.alias] = self

        print('Successfully created LUTs for custom fluid ' + self.alias)

    def generate_lookup(self, name, func):
        """
        create lookup table

        .. math::

        :param func: function to create lookup from
        :type func: callable function
        :returns: y (*scipy.interpolate.RectBivariateSpline*) - lookup table
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

        # plot table after creation?
        if self.plot:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.plot_wireframe(np.meshgrid(x2, x1)[0],
                              np.meshgrid(x2, x1)[1], y)
            ax.set_xlabel('temperature')
            ax.set_ylabel('pressure')
            ax.set_zlabel(name)
            ax.view_init(10, 225)
            plt.show()

        self.save_lookup(name, x1, x2, y)

        func = interpolate.RectBivariateSpline(x1, x2, y)
        return func

    def save_lookup(self, name, x1, x2, y):

        df = pd.DataFrame(y, columns=x2, index=x1)
        path = './LUT/' + self.alias + '/'
        if not os.path.exists(path):
            os.makedirs(path)
        df.to_csv(path + name + '.csv')

    def load_lookup(self, name):

        path = self.path.val + '/' + self.alias + '/' + name + '.csv'
        df = pd.read_csv(path, index_col=0)

        x1 = df.index.get_values()
        x2 = np.array(list(map(float, list(df))))
        y = df.as_matrix()

        # plot table after creation?
        if self.plot:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.plot_wireframe(np.meshgrid(x2, x1)[0],
                              np.meshgrid(x2, x1)[1], y)
            ax.set_xlabel('temperature')
            ax.set_ylabel('pressure')
            ax.set_zlabel(name)
            ax.view_init(10, 225)
            plt.show()

        func = interpolate.RectBivariateSpline(x1, x2, y)
        return func


def reverse_2d(params, y):
    r"""
    reverse function for lookup table

    :param params: variable function parameters
    :type params: list
    :param y: functional value, so that :math:`x_2 -
              f\left(x_1, y \right) = 0`
    :type y: float
    :returns: residual value of the function :math:`x_2 -
              f\left(x_1, y \right)`
    """
    func, x1, x2 = params[0], params[1], params[2]
    return x2 - func.ev(x1, y)


def reverse_2d_deriv(params, y):
    r"""
    derivative of the reverse function for a lookup table

    :param params: variable function parameters
    :type params: list
    :param y: functional value, so that :math:`x_2 -
              f\left(x_1, y \right) = 0`
    :type y: float
    :returns: partial derivative :math:`\frac{\partial f}{\partial y}`
    """
    func, x1 = params[0], params[1]
    return - func.ev(x1, y, dy=1)

# %%


class memorise:

    def __init__(self, fluids):

        num_fl = len(fluids)
        if num_fl > 0:
            fl = tuple(fluids)
            memorise.T_ph[fl] = np.empty((0, num_fl + 3), float)
            memorise.T_ph_f[fl] = []
            memorise.T_ps[fl] = np.empty((0, num_fl + 4), float)
            memorise.T_ps_f[fl] = []
            memorise.v_ph[fl] = np.empty((0, num_fl + 3), float)
            memorise.v_ph_f[fl] = []
            memorise.visc_ph[fl] = np.empty((0, num_fl + 3), float)
            memorise.visc_ph_f[fl] = []
            memorise.s_ph[fl] = np.empty((0, num_fl + 3), float)
            memorise.s_ph_f[fl] = []
            memorise.count = 0

    def del_memory(fluids):

        fl = tuple(fluids)

        mask = np.isin(memorise.T_ph[fl][:, -1],
                       memorise.T_ph_f[fl])
        memorise.T_ph[fl] = (memorise.T_ph[fl][mask])
        memorise.T_ph_f[fl] = []

        mask = np.isin(memorise.T_ps[fl][:, -1],
                       memorise.T_ps_f[fl])
        memorise.T_ps[fl] = (memorise.T_ps[fl][mask])
        memorise.T_ps_f[fl] = []

        mask = np.isin(memorise.v_ph[fl][:, -1],
                       memorise.v_ph_f[fl])
        memorise.v_ph[fl] = (memorise.v_ph[fl][mask])
        memorise.v_ph_f[fl] = []

        mask = np.isin(memorise.visc_ph[fl][:, -1],
                       memorise.visc_ph_f[fl])
        memorise.visc_ph[fl] = (memorise.visc_ph[fl][mask])
        memorise.visc_ph_f[fl] = []

        mask = np.isin(memorise.s_ph[fl][:, -1],
                       memorise.s_ph_f[fl])
        memorise.s_ph[fl] = (memorise.s_ph[fl][mask])
        memorise.s_ph_f[fl] = []


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

# %%


def newton(func, deriv, params, k, **kwargs):
    r"""
    find zero crossings of function func with 1-D newton algorithm,
    required for reverse functions of fluid mixtures

    :param func: function to find zero crossing in
    :type func: function
    :param deriv: derivative of the function
    :type deriv: function
    :param params: vector containing parameters for func
    :type params: list
    :param k: target value for function func
    :type k: numeric
    :returns: val (float) - val, so that func(params, val) = k

    **allowed keywords** in kwargs:

    - val0 (*numeric*) - starting value
    - valmin (*numeric*) - minimum value
    - valmax (*numeric*) - maximum value
    - imax (*numeric*) - maximum number of iterations

    .. math::

        x_{i+1} = x_{i} - \frac{f(x_{i})}{\frac{df}{dx}(x_{i})}\\
        f(x_{n}) \leq \epsilon, \; n < 10\\
        n: \text{number of iterations}
    """

    # default valaues
    val = kwargs.get('val0', 300)
    valmin = kwargs.get('valmin', 70)
    valmax = kwargs.get('valmax', 3000)
    imax = kwargs.get('imax', 10)

    # start newton loop
    res = 1
    i = 0
    while abs(res) >= err:
        # calculate function residual
        res = k - func(params, val)
        # calculate new value
        val += res / deriv(params, val)

        # check for value ranges
        if val < valmin:
            val = valmin
        if val > valmax:
            val = valmax
        i += 1

        if i > imax:
#            print('Newton algorithm was not able to find a feasible '
#                  'value for function ' + str(func) + '.')

            break

    return val

# %%


def T_mix_ph(flow):
    r"""
    calculates the temperature from pressure and enthalpy,
    uses CoolProp reverse functions for pure fluids, newton for mixtures

    - check if property has already been memorised
    - calculate property otherwise

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :returns: T (float) - temperature in K

    **fluid mixtures**

    .. math::

        T_{mix}\left(p,h\right) = T_{i}\left(pp_{i},h_{i}\right)\;
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
            val = newton(h_mix_pT, dh_mix_pdT, flow, flow[2],
                         val0=300, valmin=70, valmax=3000, imax=10)
            new = np.array([[flow[1], flow[2]] + list(flow[3].values()) +
                            [val]])
            # memorise the newly calculated value
            memorise.T_ph[fl] = np.append(memorise.T_ph[fl], new, axis=0)
            return val
        else:
            # calculate fluid property for pure fluids
            for fluid, x in flow[3].items():
                if x > err:
                    val = T_ph(flow[1], flow[2], fluid)
                    new = np.array([[flow[1], flow[2]] +
                                    list(flow[3].values()) + [val]])
                    # memorise the newly calculated value
                    memorise.T_ph[fl] = np.append(memorise.T_ph[fl],
                                                  new, axis=0)
                    return val


def T_ph(p, h, fluid):
    r"""
    returns the temperature of a pure fluid given pressure and enthalpy

    :param p: pressure
    :type p: float
    :param h: enthalpy
    :type h: float
    :param fluid: fluid alias
    :type fluid: str
    :returns: T (float) - temperature in K
    """
    if 'IDGAS::' in fluid:
        print('Ideal gas calculation not available by now.')
    elif 'TESPy::' in fluid:
        db = tespy_fluid.fluids[fluid].funcs['h_pT']
        return newton(reverse_2d, reverse_2d_deriv, [db, p, h], 0)
    else:
        return CPPSI('T', 'P', p, 'H', h, fluid)


def dT_mix_dph(flow):
    r"""
    calculates partial derivate of temperature to pressure at
    constant enthalpy and fluid composition

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :returns: dT / dp (float) - derivative in K / Pa

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
    method to calculate partial derivate of temperature to enthalpy at
    constant pressure and fluid composition

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :returns: dT / dh (float) - derivative in (K * kg) / J

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
    calculates partial derivates of temperature to fluid composition at
    constant pressure and enthalpy

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :returns: dT / dfluid (np.array of floats) - derivatives in K

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
    calculates the temperature from pressure and entropy,
    uses CoolProp reverse functions for pure fluids, newton for mixtures

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :param s: entropy in J / (kg * K)
    :type s: numeric
    :returns: T (float) - temperature in K

    **fluid mixtures**

    .. math::

        T_{mix}\left(p,s\right) = T_{i}\left(pp_{i},s_{i}\right)\;
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
            val = newton(s_mix_pT, ds_mix_pdT, flow, s,
                         val0=300, valmin=70, valmax=3000, imax=10)
            new = np.array([[flow[1], flow[2]] + list(flow[3].values()) +
                            [s, val]])
            # memorise the newly calculated value
            memorise.T_ps[fl] = np.append(memorise.T_ps[fl], new, axis=0)
            return val
        else:
            # calculate fluid property for pure fluids
            for fluid, x in flow[3].items():
                if x > err:
                    val = T_ps(flow[1], s, fluid)
                    new = np.array([[flow[1], flow[2]] +
                                    list(flow[3].values()) + [s, val]])
                    # memorise the newly calculated value
                    memorise.T_ps[fl] = np.append(memorise.T_ps[fl],
                                                  new, axis=0)
                    return val


def T_ps(p, s, fluid):
    r"""
    returns the temperature of a pure fluid given pressure and entropy

    :param p: pressure
    :type p: float
    :param s: entropy
    :type s: float
    :param fluid: fluid alias
    :type fluid: str
    :returns: T (float) - temperature in K
    """
    if 'IDGAS::' in fluid:
        print('Ideal gas calculation not available by now.')
    elif 'TESPy::' in fluid:
        db = tespy_fluid.fluids[fluid].funcs['s_pT']
        return newton(reverse_2d, reverse_2d_deriv, [db, p, s], 0)
    else:
        return CPPSI('T', 'P', p, 'S', s, fluid)

# %%


def h_mix_pT(flow, T):
    r"""
    calculates enthalpy from pressure and temperature

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :param T: temperature in K
    :type T: numeric
    :returns: h (float) - enthalpy in J / kg

    .. math::
        h_{mix}(p,T)=\sum_{i} h(pp_{i},T,fluid_{i})\;
        \forall i \in \text{fluid components}\\
        pp: \text{partial pressure}

    """

    n = molar_massflow(flow[3])

    h = 0
    for fluid, x in flow[3].items():
        if x > err:
            ni = x / molar_masses[fluid]
            h += h_pT(flow[1] * ni / n, T, fluid) * x

    return h


def h_pT(p, T, fluid):
    r"""
    returns the enthalpy of a pure fluid given pressure and temperature

    :param p: pressure
    :type p: float
    :param T: temperature
    :type T: float
    :param fluid: fluid alias
    :type fluid: str
    :returns: h (float) - enthalpy in J / kg
    """
    if 'IDGAS::' in fluid:
        print('Ideal gas calculation not available by now.')
    elif 'TESPy::' in fluid:
        return tespy_fluid.fluids[fluid].funcs['h_pT'].ev(p, T)
    else:
        return CPPSI('H', 'P', p, 'T', T, fluid)


def dh_mix_pdT(flow, T):
    r"""
    calculates partial derivate of enthalpy to temperature at constant pressure

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :param T: temperature in K
    :type T: numeric
    :returns: dh / dT (float) - derivative in J / (kg * K)

    .. math::

        \frac{\partial h_{mix}}{\partial T} =
        \frac{h_{mix}(p,T+d)-h_{mix}(p,T-d)}{2 \cdot d}
    """
    d = 2
    return (h_mix_pT(flow, T + d) - h_mix_pT(flow, T - d)) / (2 * d)

# %%


def h_mix_ps(flow, s):
    r"""
    calculates enthalpy from pressure and entropy

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :param s: entropy in J / (kg * K)
    :type s: numeric
    :returns: h (float) - enthalpy in J / kg

    .. math::
        h_{mix}(p,s)=\sum_{i} h(pp_{i},T,fluid_{i})\;
        \forall i \in \text{fluid components}\\
        pp: \text{partial pressure}

    """

    return h_mix_pT(flow, T_mix_ps(flow, s))


def h_ps(p, s, fluid):
    r"""
    returns the enthalpy of a pure fluid given pressure and entropy

    :param p: pressure
    :type p: float
    :param s: entropy
    :type s: float
    :param fluid: fluid alias
    :type fluid: str
    :returns: h (float) - enthalpy in J / kg
    """
    if 'IDGAS::' in fluid:
        print('Ideal gas calculation not available by now.')
    elif 'TESPy::' in fluid:
        db = tespy_fluid.fluids[fluid].funcs['s_pT']
        T = newton(reverse_2d, reverse_2d_deriv, [db, p, s], 0)
        return tespy_fluid.fluids[fluid].funcs['h_pT'].ev(p, T)
    else:
        return CPPSI('H', 'P', p, 'S', s, fluid)

# %%


def h_mix_pQ(flow, Q):
    """
    calculates enthalpy from pressure and quality

    .. note::

       This function works for pure fluids only!

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :param Q: fraction of vapour mass to total mass in 1
    :type Q: numeric
    :returns: h (float) - enthalpy in J / kg
    """
    n = molar_massflow(flow[3])

    h = 0
    for fluid, x in flow[3].items():
        if x > err:
            pp = flow[1] * x / (molar_masses[fluid] * n)
            pcrit = CPPSI('Pcrit', fluid)
            if pp > pcrit:
                pp = pcrit * 0.95

            h += CPPSI('H', 'P', pp, 'Q', Q, fluid) * x

    return h


def dh_mix_dpQ(flow, Q):
    r"""
    calculates partial derivative of enthalpy to pressure at constant quality

    .. note::

       This function works for pure fluids only!

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :param Q: fraction of vapour mass to total mass in 1
    :type Q: numeric
    :returns: dh / dp (float) - derivative in J / (kg * Pa)

    .. math::

        \frac{\partial h_{mix}}{\partial p} =
        \frac{h_{mix}(p+d,Q)-h_{mix}(p-d,Q)}{2 \cdot d}\\
        Q: \text{vapour mass fraction}
    """
    d = 1
    u = flow.copy()
    l = flow.copy()
    u[1] += d
    l[1] -= d
    return (h_mix_pQ(u, Q) - h_mix_pQ(l, Q)) / (2 * d)

# %%


def v_mix_ph(flow):
    r"""
    calculates specific volume from pressure and enthalpy
    uses CoolProp reverse functions for pure fluids, newton for mixtures

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :returns: v (float) - specific volume in kg / m :sup:`3`

    **fluid mixtures**

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
            new = np.array([[flow[1], flow[2]] + list(flow[3].values()) +
                            [val]])
            # memorise the newly calculated value
            memorise.v_ph[fl] = np.append(memorise.v_ph[fl], new, axis=0)
            return val
        else:
            # calculate fluid property for pure fluids
            for fluid, x in flow[3].items():
                if x > err:
                    val = 1 / d_ph(flow[1], flow[2], fluid)
                    new = np.array([[flow[1], flow[2]] +
                                    list(flow[3].values()) + [val]])
                    # memorise the newly calculated value
                    memorise.v_ph[fl] = np.append(memorise.v_ph[fl],
                                                  new, axis=0)
                    return val


def d_ph(p, h, fluid):
    r"""
    returns the density of a pure fluid given pressure and enthalpy

    :param p: pressure
    :type p: float
    :param h: enthalpy
    :type h: float
    :param fluid: fluid alias
    :type fluid: str
    :returns: d (float) - density in kg / m^3
    """
    if 'IDGAS::' in fluid:
        print('Ideal gas calculation not available by now.')
    elif 'TESPy::' in fluid:
        db = tespy_fluid.fluids[fluid].funcs['h_pT']
        T = newton(reverse_2d, reverse_2d_deriv, [db, p, h], 0)
        return tespy_fluid.fluids[fluid].funcs['d_pT'].ev(p, T)
    else:
        return CPPSI('D', 'P', p, 'H', h, fluid)

# %%


def v_mix_pT(flow, T):
    r"""
    calculates specific volume from pressure and temperature

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :param T: temperature in K
    :type T: numeric
    :returns: v (float) - specific volume in kg / m :sup:`3`

    .. math::
        v_{mix}(p,T)=\sum_{i} v(pp_{i},T,fluid_{i})\;
        \forall i \in \text{fluid components}\\
        pp: \text{partial pressure}
    """
    n = molar_massflow(flow[3])

    d = 0
    for fluid, x in flow[3].items():
        if x > err:
            pp = flow[1] * x / (molar_masses[fluid] * n)
            d += d_pT(pp, T, fluid) * x

    return 1 / d


def d_mix_pT(flow, T):
    r"""
    calculates specific volume from pressure and temperature

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :param T: temperature in K
    :type T: numeric
    :returns: d (float) - density in m :sup:`3` / kg

    .. math::
        \rho_{mix}(p,T)=\sum_{i} \rho(pp_{i},T,fluid_{i})\;
        \forall i \in \text{fluid components}\\
        pp: \text{partial pressure}
    """
    return 1 / v_mix_pT(flow, T)


def d_pT(p, T, fluid):
    r"""
    returns the density of a pure fluid given pressure and temperature

    :param p: pressure
    :type p: float
    :param T: temperature
    :type T: float
    :param fluid: fluid alias
    :type fluid: str
    :returns: d (float) - density in kg / m^3
    """
    if 'IDGAS::' in fluid:
        print('Ideal gas calculation not available by now.')
    elif 'TESPy::' in fluid:
        return tespy_fluid.fluids[fluid].funcs['d_pT'].ev(p, T)
    else:
        return CPPSI('D', 'P', p, 'T', T, fluid)

# %%


def visc_mix_ph(flow):
    r"""
    calculates dynamic viscosity from pressure and enthalpy,
    uses CoolProp reverse functions for pure fluids, newton for mixtures

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :returns: v (float) - specific volume in Pa s

    **fluid mixtures**

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
            new = np.array([[flow[1], flow[2]] + list(flow[3].values()) +
                            [val]])
            # memorise the newly calculated value
            memorise.visc_ph[fl] = np.append(memorise.visc_ph[fl], new, axis=0)
            return val
        else:
            # calculate fluid property for pure fluids
            for fluid, x in flow[3].items():
                if x > err:
                    val = visc_ph(flow[1], flow[2], fluid)
                    new = np.array([[flow[1], flow[2]] +
                                    list(flow[3].values()) + [val]])
                    # memorise the newly calculated value
                    memorise.visc_ph[fl] = np.append(memorise.visc_ph[fl],
                                                     new, axis=0)
                    return val


def visc_ph(p, h, fluid):
    r"""
    returns the dynamic viscosity of a pure fluid given pressure and enthalpy

    :param p: pressure
    :type p: float
    :param h: enthalpy
    :type h: float
    :param fluid: fluid alias
    :type fluid: str
    :returns: visc (float) - dynamic viscosity in Pa s
    """
    if 'IDGAS::' in fluid:
        print('Ideal gas calculation not available by now.')
    elif 'TESPy::' in fluid:
        db = tespy_fluid.fluids[fluid].funcs['h_pT']
        T = newton(reverse_2d, reverse_2d_deriv, [db, p, h], 0)
        return tespy_fluid.fluids[fluid].funcs['visc_pT'].ev(p, T)
    else:
        return CPPSI('V', 'P', p, 'H', h, fluid)

# %%


def visc_mix_pT(flow, T):
    r"""
    calculates dynamic viscosity from pressure and temperature

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :param T: temperature in K
    :type T: numeric
    :returns: v (float) - specific volume in kg / m :sup:`3`

    .. math::
        \eta_{mix}(p,T)=\frac{\sum_{i} \left( \eta(p,T,fluid_{i}) \cdot y_{i}
        \cdot \sqrt{M_{i}} \right)}
        {\sum_{i} \left(y_{i} \cdot \sqrt{M_{i}} \right)}\;
        \forall i \in \text{fluid components}\\
        y: \text{volume fraction}\\
        M: \text{molar mass}
    """
    n = molar_massflow(flow[3])

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
    returns the dynamic viscosity of a pure fluid given pressure and
    temperature

    :param p: pressure
    :type p: float
    :param T: temperature
    :type T: float
    :param fluid: fluid alias
    :type fluid: str
    :returns: visc (float) - dynamic viscosity in Pa s
    """
    if 'IDGAS::' in fluid:
        print('Ideal gas calculation not available by now.')
    elif 'TESPy::' in fluid:
        return tespy_fluid.fluids[fluid].funcs['visc_pT'].ev(p, T)
    else:
        return CPPSI('V', 'P', p, 'T', T, fluid)

# %%


def s_mix_ph(flow):
    r"""
    calculates entropy from pressure and enthalpy
    uses CoolProp reverse functions for pure fluids, newton for mixtures

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :returns: s (float) - entropy in J / (kg * K)

    **fluid mixtures**

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
            new = np.array([[flow[1], flow[2]] + list(flow[3].values()) +
                            [val]])
            # memorise the newly calculated value
            memorise.s_ph[fl] = np.append(memorise.s_ph[fl], new, axis=0)
            return val
        else:
            # calculate fluid property for pure fluids
            for fluid, x in flow[3].items():
                if x > err:
                    val = s_ph(flow[1], flow[2], fluid)
                    new = np.array([[flow[1], flow[2]] +
                                    list(flow[3].values()) + [val]])
                    # memorise the newly calculated value
                    memorise.s_ph[fl] = np.append(memorise.s_ph[fl],
                                                  new, axis=0)
                    return val


def s_ph(p, h, fluid):
    r"""
    returns the entropy of a pure fluid given pressure and enthalpy

    :param p: pressure
    :type p: float
    :param h: enthalpy
    :type h: float
    :param fluid: fluid alias
    :type fluid: str
    :returns: s (float) - entropy in J / (kg * K)
    """
    if 'IDGAS::' in fluid:
        print('Ideal gas calculation not available by now.')
    elif 'TESPy::' in fluid:
        db = tespy_fluid.fluids[fluid].funcs['h_pT']
        T = newton(reverse_2d, reverse_2d_deriv, [db, p, h], 0)
        return tespy_fluid.fluids[fluid].funcs['s_pT'].ev(p, T)
    else:
        return CPPSI('S', 'P', p, 'H', h, fluid)

# %%


def s_mix_pT(flow, T):
    r"""
    calculates entropy from pressure and temperature
    uses CoolProp reverse functions for pure fluids, newton for mixtures

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :param T: temperature in K
    :type T: numeric
    :returns: s (float) - entropy in J / (kg * K)

    .. math::
        s_{mix}(p,T)=\sum_{i} x_{i} \cdot s(pp_{i},T,fluid_{i})-
        \sum_{i} x_{i} \cdot R_{i} \cdot \ln \frac{pp_{i}}{p}\;
        \forall i \in \text{fluid components}\\
        pp: \text{partial pressure}\\
        R: \text{gas constant}
    """

# second method seems to be faster, maybe speed can further be improved
#    n = molar_massflow(flow[3])
#
#    s = 0
#    for fluid, x in flow[3].items():
#        if x > err:
#            if 'TESPy::' in fluid:
#                for f, xi in tespy_fluid.fluids[fluid].fluid.items():
#                    if xi > err:
#                        pp = flow[1] * xi * x / (molar_masses[f] * n)
#                        s += s_pT(pp, T, f) * xi * x
#                        s -= (xi * x * gas_constants[f] / molar_masses[f] *
#                              math.log(pp / flow[1]))
#
#            else:
#                pp = flow[1] * x / (molar_masses[fluid] * n)
#                s += s_pT(pp, T, fluid) * x
#                s -= (x * gas_constants[fluid] / molar_masses[fluid] *
#                      math.log(pp / flow[1]))
#
#    return s

    n = molar_massflow(flow[3])

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
            s -= (x * gas_constants[fluid] / molar_masses[fluid] *
                  math.log(pp / flow[1]))

    return s


def s_pT(p, T, fluid):
    r"""
    returns the entropy of a pure fluid given pressure and temperature

    :param p: pressure
    :type p: float
    :param T: temperature
    :type T: float
    :param fluid: fluid alias
    :type fluid: str
    :returns: s (float) - entropy in J / (kg * K)
    """
    if 'IDGAS::' in fluid:
        print('Ideal gas calculation not available by now.')
    elif 'TESPy::' in fluid:
        return tespy_fluid.fluids[fluid].funcs['s_pT'].ev(p, T)
    else:
        return CPPSI('S', 'P', p, 'T', T, fluid)


def ds_mix_pdT(flow, T):
    r"""
    calculates partial derivate of entropy to temperature at constant pressure

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :param T: temperature in K
    :type T: numeric
    :returns: ds / dT (float) - derivative in J / (kg * K :sup:`2`)

    .. math::

        \frac{\partial s_{mix}}{\partial T} =
        \frac{s_{mix}(p,T+d)-s_{mix}(p,T-d)}{2 \cdot d}
    """

    d = 2
    return (s_mix_pT(flow, T + d) - s_mix_pT(flow, T - d)) / (2 * d)

# %%


def molar_massflow(flow):
    r"""
    calculates molar massflow

    :param flow: vector containing [mass flow, pressure, enthalpy, fluid]
    :type flow: list
    :returns: mm (float) - molar massflow in mol / s

    .. math::
        mm = \sum_{i} \left( \frac{x_{i}}{M_{i}} \right)
    """
    mm = 0
    for fluid, x in flow.items():
        if x > err:
            try:
                mm += x / molar_masses[fluid]
            except:
                mm += x / CPPSI('molar_mass', fluid)

    return mm

# %%


def num_fluids(fluids):
    r"""
    calculates number of fluids in fluid vector

    :param fluids: fluid vector {fluid: mass fraction}
    :type fluids: dict
    :returns: n (int) - number of fluids in fluid vector in 1

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
    returns the name of the single fluid (x=1) in a fluid vector

    :param fluids: fluid vector {fluid: mass fraction}
    :type fluids: dict
    :returns: fluid (str) - name of the fluid
    """

    if num_fluids(fluids) == 1:
        for fluid, x in fluids.items():
            if x > err:
                return fluid
    else:
        return []

# %%


def fluid_structure(fluid):
    """
    gets the chemical formular of a fluid

    :param fluid: alias of the fluid
    :type fluid: str
    :returns: parts (dict) - returns the elements of the fluid {'element': n}
    """
    parts = {}
    for element in CP.get_fluid_param_string(fluid, 'formula').split('}'):
        if element != '':
            el = element.split('_{')
            parts[el[0]] = int(el[1])

    return parts

# %%


def lamb(re, ks, d):
    r"""
    calculates darcy friction factor

    :param re: reynolds number in 1
    :type re: numeric
    :param ks: roughness in m
    :type ks: numeric
    :param d: pipe diameter in m
    :type d: numeric
    :returns: lambda (float) - darcy friction factor in 1

    **laminar flow:** :math:`re \leq 2320`

    .. math::
        \lambda = \frac{64}{re}

    **turbulent flow:** :math:`re > 2320`

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
                return newton(lamb_smooth, dlamb_smooth_dl, [re], 0,
                              val0=l0, valmin=0.00001, valmax=0.2)

        elif re * ks / d > 1300:
            return 1 / (2 * math.log(3.71 * d / ks, 10)) ** 2

        else:
            l0 = 0.002
            return newton(lamb_trans, dlamb_trans_dl, [re, ks, d], 0,
                          val0=l0, valmin=0.0001, valmax=0.2)


def lamb_smooth(params, l):
    re = params[0]
    return 2 * math.log(re * math.sqrt(l), 10) - 0.8 - 1 / math.sqrt(l)


def dlamb_smooth_dl(params, l):
    return 1 / (l * math.log(10)) + 1 / 2 * l ** (-1.5)


def lamb_trans(params, l):
    re, ks, d = params[0], params[1], params[2]
    return (2 * math.log(2.51 / (re * math.sqrt(l)) + ks / d * 0.269, 10) +
            1 / math.sqrt(l))


def dlamb_trans_dl(params, l):
    d = 0.001
    return (lamb_trans(params, l + d) - lamb_trans(params, l - d)) / (2 * d)
