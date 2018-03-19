"""
.. module:: components.characteristics
    :synopsis:

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>

**TODO**

- add documentation

**Improvements**

- build generic architecture for characteristics
"""

from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import copy
import math


def find_length(xs, ys):
    dx = np.diff(xs)
    dy = np.diff(ys)
    return np.sum(np.sqrt(dx ** 2 + dy ** 2))


class compressor:
    r"""

    generic characteristic map for axial compressors

    - links mass flow to pressure rise and isentropic efficiency
    - includes a vigv angle

    the map can be plotted using :code:`map.plot()`

    **literature**

    compressor map:

    - Marcin Plis, Henryk Rusinowski (2016): Mathematical modeling of an
      axial compressor in a gas turbine cycle. Journal of Power
      Technologies 96 (3), pp. 194-199.

    vigv:

    - GasTurb GmbH (2015): GasTurb 12.
    """

    def __init__(self,  ** kwargs):

        self.raw_data = None
        self.beta = None

        for key in kwargs:
            if key in self.keywordargs():
                self.__dict__.update({key: kwargs[key]})

        for key in self.keywordargs():
            if self.__dict__[key] is None:
                self.__dict__.update({key: self.default(key)})

        self.pr = copy.deepcopy(self.raw_data)
        self.eta = copy.deepcopy(self.raw_data)

# create lut for each speedline
        self.pr_lut = {}
        self.eta_lut = {}
        for key in sorted(self.raw_data):
# pressure ratio
            x = [it[0] for it in self.raw_data[key]]
            y = [it[1] for it in self.raw_data[key]]
# sorted values of x and y
            self.pr[key] = interp1d(x, y, 'linear')
            y = self.pr[key].y
            x = self.pr[key].x
            l = np.insert(np.cumsum(np.sqrt(np.diff(x) ** 2 +
                                            np.diff(y) ** 2)), 0, 0)
            self.pr_lut[key] = pd.DataFrame({'x': x, 'y': y}, index=l)

# add values for beta lines
            for i in np.linspace(0, self.pr_lut[key].index[-1], self.beta):
                if i not in self.pr_lut[key].index:
                    self.pr_lut[key].loc[i] = [np.nan, np.nan]

            self.pr_lut[key] = self.pr_lut[key].sort_index()
            self.pr_lut[key] = self.pr_lut[key].interpolate(method='index')

# isentropic efficiency
            x = [it[0] for it in self.raw_data[key]]
            y = [it[2] for it in self.raw_data[key]]
# sorted values of x and y
            self.eta[key] = interp1d(x, y, 'linear')
            y = self.eta[key].y
            x = self.eta[key].x
            l = np.insert(np.cumsum(np.sqrt(np.diff(x) ** 2 +
                                            np.diff(y) ** 2)), 0, 0)
            self.eta_lut[key] = pd.DataFrame({'x': x, 'y': y}, index=l)

# add values for beta lines
            for i in np.linspace(0, self.eta_lut[key].index[-1], self.beta):
                if i not in self.eta_lut[key].index:
                    self.eta_lut[key].loc[i] = [np.nan, np.nan]

            self.eta_lut[key] = self.eta_lut[key].sort_index()
            self.eta_lut[key] = self.eta_lut[key].interpolate(method='index')

        self.pr_beta = {}
        self.eta_beta = {}
        for i in range(self.beta):
            x = []
            y = []
            for key in sorted(self.raw_data):
                x += [self.pr_lut[key].loc[np.linspace(
                    0, self.pr_lut[key].index[-1], self.beta
                )[i]].x]
                y += [self.pr_lut[key].loc[np.linspace(
                    0, self.pr_lut[key].index[-1], self.beta
                )[i]].y]
            self.pr_beta[i] = pd.DataFrame({'x': x, 'y': y},
                                           index=sorted(self.raw_data))
            x = []
            y = []
            for key in sorted(self.raw_data):
                x += [self.eta_lut[key].loc[np.linspace(
                    0, self.eta_lut[key].index[-1], self.beta
                )[i]].x]
                y += [self.eta_lut[key].loc[np.linspace(
                    0, self.eta_lut[key].index[-1], self.beta
                )[i]].y]
            self.eta_beta[i] = pd.DataFrame({'x': x, 'y': y},
                                            index=sorted(self.raw_data))

# allowed keyword arguments
    def keywordargs(self):
        return ['raw_data', 'beta']

# default values
    def default(self, key):
        """
        source of map:
            Marcin Plis, Henryk Rusinowski (2016): Mathematical modeling of an
            axial compressor in a gas turbine cycle. Journal of Power
            Technologies 96 (3), pp. 194-199.
        """

        # default map
        default_map = {}
        default_map[1.062] = [[1.045, 1.441, 0.879], [1.056, 0.805, 0.887],
                              [1.052, 1.176, 0.925]]

        default_map[1.029] = [[1.014, 1.340, 0.948], [1.026, 1.082, 0.967],
                              [1.036, 0.767, 0.868]]

        default_map[1.000] = [[0.948, 1.195, 0.956], [0.961, 1.176, 0.958],
                              [0.974, 1.151, 0.962], [0.987, 1.101, 0.981],
                              [0.994, 1.057, 0.992], [1.000, 1.000, 1.000],
                              [1.005, 0.893, 0.9909], [1.006, 0.748, 0.914]]

        default_map[0.971] = [[0.874, 1.050, 0.977], [0.909, 1.019, 0.977],
                              [0.922, 1.000, 0.987], [0.935, 0.969, 1.004],
                              [0.948, 0.918, 1.008], [0.961, 0.861, 1.004],
                              [0.964, 0.672, 0.920]]

        default_map[0.946] = [[0.767, 0.931, 0.891], [0.819, 0.912, 1.002],
                              [0.838, 0.893, 1.013], [0.858, 0.861, 1.017],
                              [0.871, 0.817, 1.017], [0.877, 0.767, 1.013],
                              [0.880, 0.616, 0.925]]

        default_map[0.870] = [[0.59, 0.65, 0.887], [0.618, 0.641, 0.929],
                              [0.657, 0.616, 0.962], [0.676, 0.597, 0.969],
                              [0.7, 0.553, 0.975], [0.702, 0.528, 0.967],
                              [0.709, 0.477, 0.948], [0.713, 0.37, 0.887]]

        default_map[0.810] = [[0.46, 0.502, 0.872], [0.534, 0.483, 0.918],
                              [0.573, 0.452, 0.948], [0.586, 0.424, 0.944],
                              [0.599, 0.38, 0.925], [0.605, 0.348, 0.906],
                              [0.612, 0.276, 0.879]]

        default_val = {
            'raw_data': default_map,
            'beta': 10
        }

        return default_val[key]

# adding speedlines to the compressor map
    def add_speedline(self, n):
        n = round(n, 3)
        if n not in self.pr.keys():
            if n > max(self.pr.keys()) or n < min(self.pr.keys()):
                return

            x = []
            y = []
            for i in range(self.beta):
                self.pr_beta[i].loc[n] = [np.nan, np.nan]
                self.pr_beta[i] = self.pr_beta[i].sort_index()
                self.pr_beta[i] = self.pr_beta[i].interpolate(method='index')
                x += [self.pr_beta[i].loc[n].x]
                y += [self.pr_beta[i].loc[n].y]

            self.pr[n] = interp1d(x, y, 'linear')

            x = []
            y = []
            for i in range(self.beta):
                self.eta_beta[i].loc[n] = [np.nan, np.nan]
                self.eta_beta[i] = self.eta_beta[i].sort_index()
                self.eta_beta[i] = self.eta_beta[i].interpolate(method='index')
                x += [self.eta_beta[i].loc[n].x]
                y += [self.eta_beta[i].loc[n].y]

            self.eta[n] = interp1d(x, y, 'linear')

# get speedline as interpolation object
    def get_speedline(self, n, vigv):
        """
        source of speedline adaption by igv:
            GasTurb GmbH (2015): GasTurb 12.
        """
        n = round(n, 3)
        self.add_speedline(n)
        pr = interp1d(self.pr[n].x * (1 - vigv / 100),
                      self.pr[n].y * (1 - vigv / 100), 'linear')
        eta = interp1d(self.eta[n].x * (1 - vigv / 100),
                       self.eta[n].y * (1 - vigv ** 2 / 10000), 'linear')
        return pr, eta

# calculate feasible vigv angles for given speedline and non dimensional mass flow
    def get_vigv_range(self, n, m):
        n = round(n, 3)
        if n not in self.pr.keys():
            self.add_speedline(n)

        vigv_min = 100 * (1 - m / self.eta_beta[0].loc[n].x)
        vigv_max = 100 * (1 - m / self.eta_beta[self.beta - 1].loc[n].x)

        return vigv_min, vigv_max

# get pressure ratio
    def get_pr(self, n, m, vigv):
        return self.get_speedline(n, vigv)[0](m)

# get isentropic efficiency
    def get_eta(self, n, m, vigv):
        return self.get_speedline(n, vigv)[1](m)

# calculate vigv angle for given speedline, non dimensional mass flow and pressure ratio
    def get_vigv(self, n, m, p):
        n = round(n, 3)
        if n not in self.pr.keys():
            self.add_speedline(n)

        tolerance = 1e-12
        d = 1e-3
        res = 1
        deriv = 1

        vigv_range = self.get_vigv_range(n, m)
        vigv = (vigv_range[0] + vigv_range[1]) / 2
        vigv_hist = [vigv]
        z_hist = [res / deriv]

        while abs(res) >= tolerance:
            try:
                res = p - self.get_pr(n, m, vigv)
                deriv = ((self.get_pr(n, m, vigv + d) -
                          self.get_pr(n, m, vigv - d)) / (2 * d))
                vigv += res / deriv
                z_hist += [res / deriv]
            except:
                if vigv < vigv_range[0]:
                    vigv = vigv_range[0] + 1e-2
                if vigv > vigv_range[1]:
                    vigv = vigv_range[1] - 1e-2

            vigv_hist += [vigv]

            if ((len(vigv_hist) > 10 and
                vigv_hist[(len(vigv_hist) - 10):] == 10 * [vigv_hist[-1]]) or
                (len(z_hist) > 5 and
                z_hist[(len(z_hist) - 5):] == 5 * [z_hist[-1]])):
                raise ValueError('Given pressure ratio can not be archieved'
                                 ' with given speedline.')

#        print(time.time() - tmp)
        return vigv

# plot the compressor map
    def plot(self):
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        for i in range(self.beta):
            ax1.plot(self.pr_beta[i].x, self.pr_beta[i].y, 'xk', ms=3)

        for i in range(self.beta):
            ax2.plot(self.eta_beta[i].x, self.eta_beta[i].y, 'xr', ms=3)

        ax1.set_ylabel('$p$ / $p_\mathrm{ref}$')
        ax2.set_ylabel('$\eta$ / $\eta_\mathrm{ref}$')
        ax1.set_xlabel('$m$ / $m_\mathrm{ref}$')
        ax2.set_ylim([0, 1.1])
        ax1.set_ylim([0, 2])
        plt.sca(ax1)
        plt.show()


class pump:
    r"""

    generic characteristic for pumps

    - links isentropic efficiency :math:`\eta_{s,p}` to volumetric flow
      :math:`\dot{V}`
    - uses a distorted quadratic function:

    .. math::

        \eta_{s,p} = \left( a \cdot \dot{V}^2 + b \cdot \dot{V} \right) \cdot
        e^{k \cdot \dot{V}}

    - function parameters are deriven from design status
    - specific rotational speed :math:`n_q`:

    .. math::

        n_q = \frac{n \cdot \sqrt{\dot{V}}}{H^{0,75}}

    .. note::

        The calculation with pump characteristics is unstable, better use
        constant value for isentropic efficiency!

    **literature**

    - Wolfgang Wesche (2012): Radiale Kreiselpumpen - Berechnung und
      Konstruktion der hydrodynamischen Komponenten. Berlin: Springer.
    """

    def __init__(self, v_opt, eta_s, H_opt):

        n_q = 3000 * math.sqrt(v_opt) / ((H_opt) ** 0.75)
        v_0 = v_opt * 3.1 * n_q ** (-0.15)
        self.k = (v_0 - 2 * v_opt) / (v_opt ** 2 - v_0 * v_opt)
        self.a = eta_s / ((v_opt ** 2 - v_0 * v_opt) *
                          math.exp(self.k * v_opt))
        self.b = -self.a * v_0

    def eta(self, v):
        return (self.a * v ** 2 + self.b * v) * math.exp(self.k * v)


class turbine:
    r"""

    generic characteristics for turbine isentropic efficiency

    - links isentropic efficiency :math:`\eta_\mathrm{s,t}` to keyfigure
      :math:`\nu`

    .. math::

        \eta_\mathrm{s,t}=f\left(\frac{\nu}{\nu_\mathrm{ref}} \right)

        \frac{\nu}{\nu_\mathrm{ref}}=\frac{\sqrt{\Delta h_\mathrm{s,ref}}}
        {\sqrt{\Delta h_\mathrm{s}}}

    - values from Traupel (see literature)
    - maximum value of isentropic efficiency in characteristic is assigned to
      isentropic efficiency in reference state with
      :math:`\frac{\nu}{\nu_{ref}}=1`.

    **literature**

    - Walter Traupel (2001): Thermische Turbomaschinen Band 2. Berlin: Spinger.
    """

    def __init__(self, eta_s0):

        self.nu = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        self.eta_s = np.array([0.0, 0.6, 0.85, 0.875, 0.75, 0.5])
        self.char = interp1d(self.nu, self.eta_s, kind='cubic')
        self.nu = np.linspace(self.nu[0], self.nu[-1], 20)
        self.eta_s = self.char(self.nu)
        self.char = interp1d(self.nu, self.eta_s /
                             self.eta_s[np.argmax(self.char(self.nu))],
                             kind='linear')

        print(self.char.x, self.char.y)

    def eta(self, nu):
        return self.char(nu)


class characteristics:

    def __init__(self, **kwargs):

        for key in kwargs:
            if key not in self.attr():
                msg = ('Invalid keyword ' + key + '. Available keywords for '
                       'kwargs are: ' + str(self.attr()) + '.')
                raise KeyError(msg)

        x_default, y_default = self.default

        self.x = kwargs.get('x', x_default)
        self.y = kwargs.get('y', y_default)
        self.char = interp1d(x, y, kind='linear', bounds_error=True)

    def default (self):
        x = np.array([0, 1, 2])
        y = np.array([1, 1, 1])
        return x, y

    def attr(self):
        return ['x', 'y']


    def f_x(self, x):
        r"""

        """
        return self.char(x)


class turbine(characteristics):

    def default(self, key):

        x = {}
        y = {}

        x['EBS_ST'] = np.array([0, 0.5, 0.6, 0.7, 0.8,
                                0.9, 0.95, 1, 1.05, 1.1])
        y['EBS_ST'] = np.array([0.98, 0.991, 0.993, 0.995, 0.9975,
                                0.999, 0.9998, 1, 0.9995, 0.995])

        x['EBS_GT'] = np.array([0, 0.4, 0.7, 1, 1.2])
        y['EBS_GT'] = np.array([0.85, 0.9, 0.95, 1, 1.1])

        y['TRAUPEL'] =

        return x[key], y[key]
