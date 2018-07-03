"""
.. module:: components.characteristics
    :synopsis:

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>

**TODO**

- add documentation for compressor maps

**Improvements**

- add pump and compressor characteristics to characteristics class
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


class characteristics:
    r"""

    characteristics for components performs linear interpolation on given pairs
    of values. Value pairs may be user-specified, default values are used
    instead.

    :param method: what
    :type method: tespy.components.components.component
    :returns: no return value

    **allowed keywords** in kwargs (also see characteristics.attr()):

    - x, y (*numeric*) - values for function parameters x and function values y
    - method (*str*) - keyword method is necessary, if you do not provide any
      x or y data. TESPy will use the characteristics as stated in the
      subclasses
    """

    def __init__(self, **kwargs):

        for key in kwargs:
            if key not in self.attr():
                msg = ('Invalid keyword ' + key + '. Available keywords for '
                       'kwargs are: ' + str(self.attr()) + '.')
                raise KeyError(msg)

        # in case of various default characteristics
        method = kwargs.get('method', 'default')

        self.x = kwargs.get('x', None)
        self.y = kwargs.get('y', None)

        if self.x is None:
            self.x = self.default(method)[0]
        if self.y is None:
            self.y = self.default(method)[1]

        self.char = interp1d(self.x, self.y, kind='linear', bounds_error=True)

    def default(self, key):

        x = {}
        y = {}

        x['default'] = np.array([0, 1, 2])
        y['default'] = np.array([1, 1, 1])

        return x[key], y[key]

    def attr(self):
        return ['x', 'y', 'method']

    def f_x(self, x):
        return self.char(x)

    def get_attr(self, key):
        """
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
            return None


class turbine(characteristics):
    r"""

    generic characteristics for turbine isentropic efficiency

    - links isentropic efficiencay :math:`\eta_\mathrm{s,t}` to keyfigure
    - three default characteristic lines available (see literature and method
      turbine.default())

    **literature**

    - Walter Traupel (2001): Thermische Turbomaschinen Band 2. Berlin: Spinger.
      -> TRAUPEL
    """

    def default(self, key):
        r"""

        default **characteristic lines** for turbines **are designed for the
        following cases**:

            \frac{\eta_\mathrm{s,t}}{\eta_\mathrm{s,t,ref}}=f\left(X \right)

        available lines characteristics:

        **GENERIC**

        .. math::

            \frac{\eta_\mathrm{s,t}}{\eta_\mathrm{s,t,ref}}=f\left(X \right)

            \text{choose calculation method for X}

            X = \begin{cases}
            \frac{\dot{m}}{\dot{m}_{ref}} & \text{mass flow}\\
            \frac{\dot{V}}{\dot{V}_{ref}} & \text{volumetric flow}\\
            \frac{p_1 \cdot p_{2,ref}}{p_{1,ref} \cdot p_2} &
            \text{pressure ratio}
            \end{cases}

        .. image:: _images/GENERIC.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

        **TRAUPEL**

        .. math::

           X=\frac{
            \sqrt{\Delta h_\mathrm{s,ref}}}{\sqrt{\Delta h_\mathrm{s}}}

        .. image:: _images/TRAUPEL.svg
           :scale: 100 %
           :alt: alternative text
           :align: center
        """

        if key == 'default':
            return np.array([0, 1, 2]), np.array([1, 1, 1])

        x = {}
        y = {}

        x['GENERIC'] = np.array([0, 0.5, 0.8, 0.95, 1, 1.05, 1.2])
        y['GENERIC'] = np.array([0.975, 0.985, 0.994, 0.999, 1, 0.999, 0.99])
        x['TRAUPEL'] = np.array([0.0, 0.1905, 0.3810, 0.5714, 0.7619, 0.9524,
                                 1.0, 1.1429, 1.3333, 1.5238, 1.7143, 1.9048])
        y['TRAUPEL'] = np.array([0.0, 0.3975, 0.6772, 0.8581, 0.9593, 0.9985,
                                 1.0, 0.9875, 0.9357, 0.8464, 0.7219, 0.5643])

        return x[key], y[key]


class heat_ex(characteristics):
    r"""

    generic characteristics for heat exchanger heat transfer coefficient

    - links heat transfer coefficient :math:`kA` to keyfigures
    - different default characteristic lines available for different types of
      heat exchangers (see method heat_exchanger.default())
    """

    def default(self, key):
        r"""

        default **characteristic lines** for heat exchangers **are designed for
        the following cases**:

        .. math::

            \frac{kA_\mathrm{s,t}}{kA_\mathrm{s,t,ref}}=f_1\left(x_1 \right)
            \cdot f_2\left(x_2 \right)

        available lines characteristics:

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
            return np.array([0, 1, 2]), np.array([1, 1, 1])

        x = {}
        y = {}

        x['COND_HOT'] = np.array([0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35,
                                  0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,
                                  0.8, 0.85, 0.9, 0.95, 1, 1.5, 2])
        y['COND_HOT'] = np.array(
                [0.030, 0.158, 0.344, 0.567, 0.838, 0.888, 0.906, 0.921, 0.933,
                 0.943, 0.952, 0.959, 0.966, 0.972, 0.977, 0.982, 0.986, 0.990,
                 0.994, 0.997, 1.000, 1.021, 1.033])

        x['COND_COLD'] = np.array([0.01, 0.04, 0.07, 0.11, 0.15, 0.2, 0.25,
                                   0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65,
                                   0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.5, 2])
        y['COND_COLD'] = np.array(
                [0.019, 0.075, 0.134, 0.192, 0.243, 0.303, 0.359, 0.412, 0.463,
                 0.512, 0.559, 0.604, 0.648, 0.691, 0.733, 0.774, 0.813, 0.852,
                 0.890, 0.928, 0.964, 1.000, 1.327, 1.612])

        x['EVA_HOT'] = np.array([0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35,
                                 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,
                                 0.8, 0.85, 0.9, 0.95, 1, 1.5, 2])
        y['EVA_HOT'] = np.array(
                [0.030, 0.158, 0.245, 0.313, 0.373, 0.427, 0.477, 0.524, 0.569,
                 0.611, 0.652, 0.692, 0.730, 0.767, 0.803, 0.838, 0.872, 0.905,
                 0.937, 0.969, 1.000, 1.281, 1.523])

        x['EVA_COLD'] = np.array([0.01, 0.04, 0.07, 0.11, 0.15, 0.2, 0.25,
                                  0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65,
                                  0.7, 0.8, 1, 2])
        y['EVA_COLD'] = np.array(
                [0.018, 0.075, 0.134, 0.215, 0.300, 0.412, 0.531, 0.658, 0.794,
                 0.934, 0.988, 0.991, 0.994, 0.995, 0.997, 0.998, 0.999, 1.000,
                 1.001])

        x['HE_HOT'] = np.array([0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35,
                                0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,
                                0.8, 0.85, 0.9, 0.95, 1, 1.5, 2])
        y['HE_HOT'] = np.array(
                [0.030, 0.158, 0.344, 0.469, 0.535, 0.590, 0.638, 0.680, 0.718,
                 0.752, 0.783, 0.812, 0.839, 0.864, 0.887, 0.909, 0.929, 0.948,
                 0.966, 0.984, 1.000, 1.128, 1.216])

        x['HE_COLD'] = np.array([0.01, 0.04, 0.07, 0.11, 0.15, 0.2, 0.25,
                                 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65,
                                 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.5, 2])
        y['HE_COLD'] = np.array(
                [0.018, 0.075, 0.134, 0.215, 0.300, 0.412, 0.507, 0.564, 0.614,
                 0.660, 0.701, 0.739, 0.774, 0.806, 0.836, 0.864, 0.890, 0.915,
                 0.938, 0.960, 0.981, 1.000, 1.151, 1.253])

        return x[key], y[key]


class pump(characteristics):
    r"""

    generic characteristic for pumps

    - links isentropic efficiency :math:`\eta_{s,p}` to volumetric flow
      :math:`\dot{V}`
    - uses a distorted quadratic function:

    .. math::

        \eta_{s,p} = \left( a \cdot \dot{V}^2 + b \cdot \dot{V}+c \right) \cdot
        e^{k \cdot \dot{V}}

    - constraints

    .. math::

        \eta_{s,p}\left(0 \right) = 0

        \eta_{s,p}\left(\dot{V}_0 \right) = 0

        \eta_{s,p}\left(\dot{V}_{ref} \right) = 1

        \frac{\partial \eta_{s,p}}
        {\partial \dot{V}}\left(\dot{V}_{ref} \right) = 0

    - function parameters

    .. math::

        k = \frac{\dot{V}_0 - 2 \cdot \dot{V}_{ref}}
        {\dot{V}_{ref}^2-\dot{V}_0 \cdot \dot{V}_{ref}}

        a = \frac{1}{\left( \dot{V}_{ref}^2 - \dot{V}_0 \cdot \dot{V}_{ref}
        \right) \cdot e^{k \cdot \dot{V}_{ref}}}

        b = -a \cdot \dot{V}_0

        c = 0

    - volume flow without pressure rise :math:`\dot{V}_0`:

    .. math::

        \dot{V}_0 = \dot{V}_{ref} \cdot 3.1 \cdot n_q ^{-0.15}

    - specific rotational speed :math:`n_q`:

    .. math::

        n_q = \frac{333 \cdot n \cdot \sqrt{\dot{V}_{ref}}}
        {\left(g \cdot H_{ref}\right)^{0,75}}

    .. note::

        In order to stabilize the calculation with pump characteristics
        the minimum value for isentropic efficiency is limited to a quater
        of reference state efficiency!

    **literature**

    - Wolfgang Wesche (2012): Radiale Kreiselpumpen - Berechnung und
      Konstruktion der hydrodynamischen Komponenten. Berlin: Springer.

    - KSB SE & Co. KGaA (2018): Kreiselpumpenlexikon - Spezifische Drehzahl.
      Available at:
      https://www.ksb.com/kreiselpumpenlexikon/spezifische-drehzahl/186490,
      accessed on 04.04.2018.
    """

    def __init__(self, v_opt, H_opt):

        n_q = 333 * 50 * math.sqrt(v_opt) / ((9.81 * H_opt) ** 0.75)
        v_0 = v_opt * 3.1 * n_q ** (-0.15)
        self.k = (v_0 - 2 * v_opt) / (v_opt ** 2 - v_0 * v_opt)
        self.a = 1 / ((v_opt ** 2 - v_0 * v_opt) * math.exp(self.k * v_opt))
        self.b = -self.a * v_0

    def char(self, v):
        eta = (self.a * v ** 2 + self.b * v) * math.exp(self.k * v)
        if eta < 0.25:
            eta = 0.25
        return eta

    def f_x(self, x):
        return self.char(x)
