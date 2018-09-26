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
import math


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

        self.char = interp1d(self.x, self.y, kind='cubic', bounds_error=True)

    def default(self, key):

        x = {}
        y = {}

        x['default'] = np.array([0, 1, 2, 3])
        y['default'] = np.array([1, 1, 1, 1])

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

        default characteristic lines for turbines:

        .. math::

            \frac{\eta_\mathrm{s,t}}{\eta_\mathrm{s,t,ref}}=f\left(X \right)

        **GENERIC**

        .. math::

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
        {\left(g \cdot H_{ref}\right)^{0,75}}\\
        \text{assuming n=}\frac{50}{s}

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


class compressor(characteristics):
    r"""

    generic characteristic map for axial compressors

    - links mass flow to pressure rise and isentropic efficiency

    the map can be plotted using :code:`map.plot()`

    **literature**

    compressor map:

    - Marcin Plis, Henryk Rusinowski (2016): Mathematical modeling of an
      axial compressor in a gas turbine cycle. Journal of Power
      Technologies 96 (3), pp. 194-199.

    vigv:

    - GasTurb GmbH (2015): GasTurb 12.
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
        self.z1 = kwargs.get('z1', None)
        self.z2 = kwargs.get('z2', None)

        if self.x is None:
            self.x = self.default(method)[0]
        if self.y is None:
            self.y = self.default(method)[1]

        if self.z1 is None:
            self.z1 = self.default(method)[2]
        if self.z2 is None:
            self.z2 = self.default(method)[3]

    def default(self, key):

        r"""

        default characteristic map for compressor

        .. math::

            X = \sqrt{\frac{T_\mathrm{1,ref}}{T_\mathrm{1}}}

            Y = \frac{\dot{m}_\mathrm{1} \cdot p_\mathrm{1,ref}}
            {\dot{m}_\mathrm{1,ref} \cdot p_\mathrm{1} \cdot X}

            Z1 = \frac{p_2 \cdot p_\mathrm{1,ref}}{p_1 \cdot p_\mathrm{2,ref}}=
            f\left(X, Y \right)

            Z2 = \frac{\eta_\mathrm{s,c}}{\eta_\mathrm{s,c,ref}}=
            f\left(X, Y \right)

        .. image:: _images/CMAP_GENERIC_PR.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

        .. image:: _images/CMAP_GENERIC_ETA.svg
           :scale: 100 %
           :alt: alternative text
           :align: center
        """

        if key == 'default':
            return np.array([0, 1, 2]), np.array([1, 1, 1])

        x = {}
        y = {}
        z1 = {}
        z2 = {}

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

        return x[key], y[key], z1[key], z2[key]

    def get_pr_eta(self, x, y, igva):
        """
        returns the pressure ratio and isentropic efficiency at given speedline
        and correxted mass flow

        :param x: speedline
        :type x: float
        :param y: corrected mass flow
        :type y: float
        :returns: - pr (*float*) - pressure ratio
                  - eta (*float*) - isentropic efficiency
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
        """
        returns error messages for operation out of the maps bounds

        :param x: speedline
        :type x: float
        :param y: corrected mass flow
        :type y: float
        :returns: - msg (*float*) - errormessage
                  - ypos (*float*) - position of corrected mass flow
        """
        xpos = np.searchsorted(self.x, x)
        if xpos == len(self.x):
            yarr = self.y[xpos - 1]
            msg = ('##### WARNING #####\n'
                   'Operating point above compressor map range: '
                   'X=' + str(round(x, 3)) + ' with maximum of ' +
                   str(self.x[-1]))
            return msg
        elif xpos == 0:
            yarr = self.y[0]
            msg = ('##### WARNING #####\n'
                   'Operating point below compressor map range: '
                   'X=' + str(round(x, 3)) + ' with minimum of ' +
                   str(self.x[0]))
            return msg
        else:
            yfrac = (x - self.x[xpos - 1]) / (self.x[xpos] - self.x[xpos - 1])
            yarr = self.y[xpos - 1] + yfrac * (self.y[xpos] - self.y[xpos - 1])

        yarr *= (1 - igva / 100)

        ypos = np.searchsorted(yarr, y)
        if ypos == len(yarr):
            msg = ('##### WARNING #####\n'
                   'Operating point above compressor map range: '
                   'Y=' + str(round(y, 3)) + ' with maximum of ' +
                   str(yarr[-1]))
            return msg
        elif ypos == 0:
            msg = ('##### WARNING #####\n'
                   'Operating point below compressor map range: '
                   'Y=' + str(round(y, 3)) + ' with minimum of ' +
                   str(yarr[0]))
            return msg
        else:
            return None
