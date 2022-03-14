# -*- coding: utf-8 -*-
"""
@author: Mathias Hofmann

Technische Universitaet Berlin
Fachgebiet Energietechnik und Umweltschutz

@author: Francesco Witte

German Areospace Center (DLR), Institute of Networked Energy Systems

Stoffwertpolynome aus Formelsammlung

O. Knacke, O. Kubschewski, K. Hesselmann, Thermochemical Properties of
Inorganic Substances, 2nd ed., Springer, Berlin, 1991.

updated values:
A. Bejan, G. Tsatsaronis, M. Moran: Thermal Design and Optimization,
Wiley, New York, 1996.
"""

import numpy as np

# coefficients    H+        S+     a        b       c        d       M
c = {'C':      [-2.101,  -6.540,  0.109,  38.940, -0.146, -17.385, 12.011],
     'S':      [-5.242, -59.014, 14.795,  24.075,  0.071,   0.0,   32.06],
    #  'N2':     [-7.069,  51.539, 24.229,  10.521,  0.180,  -2.315, 28.0134],
     'N2':     [-9.982,  16.203, 30.418,  2.544,  -0.238,  0.0, 28.0134],
     'O2':     [-9.589,  36.116, 29.154,   6.477, -0.184,  -1.017, 31.9988],
     'H2':     [-7.823, -22.966, 26.882,   3.586,  0.105,   0.0,    2.0158],
     'CO':   [-120.809,  18.937, 30.962,   2.439, -0.280,   0.0,   28.0104],
     'CO2':  [-413.886, -87.078, 51.128,   4.368, -1.469,   0.0,   44.0098],
     'H2O': [-253.871, -11.750, 34.376,   7.841, -0.423,   0.0,   18.0152],
     'H2Ol': [-289.932, -67.147, 20.355, 109.198,  2.033,   0.0,   18.0152],
     'CH4':   [-81.242,  96.731, 11.933,  77.647,  0.142, -18.414, 16.0426],
     'SO2':  [-315.422, -43.725, 49.936,   4.766, -1.046,   0.0,   64.0588],
     'H2S':   [-32.887,   1.142, 34.911,  10.686, -0.448,   0.0,   34.0758],
     'NH3':   [-60.244, -29.402, 37.321,  18.661, -0.649,   0.0,   17.0304]}


def specheat(b, t):
    y = t*1E-3
    cp = c[b][2]+c[b][3]*y+c[b][4]*y**-2+c[b][5]*y**2
    return cp

def specheat_mass(b, t):
    y = t*1E-3
    cpm = (c[b][2]+c[b][3]*y+c[b][4]*y**-2+c[b][5]*y**2)/c[b][6]
    return cpm

def enthalpy(b, t):
    y = t*1E-3
    h = 1E3*(c[b][0]+c[b][2]*y+c[b][3]/2*y**2-c[b][4]/y+c[b][5]/3*y**3)
    return h

def enthalpy_mass(b, t):
    y = t*1E-3
    hm = 1E3 * (1E3*(c[b][0]+c[b][2]*y+c[b][3]/2*y**2-c[b][4]/y+c[b][5]/3*y**3))/c[b][6]
    return hm

def T_enthalpy_mass(params, t):
    b, h = params
    return h - enthalpy_mass(b, t)

def T_enthalpy_mass_deriv(params, t):
    d = 0.001
    return (T_enthalpy_mass(params, t + d) -
            T_enthalpy_mass(params, t - d)) / (2 * d)

def entropy(b, t):
    y = t*1E-3
    s = c[b][1]+c[b][2]*np.log(t)+c[b][3]*y-c[b][4]/2*y**-2+c[b][5]/2*y**2
    return s

def entropy_mass(b, t):
    y = t*1E-3
    sm = (c[b][1]+c[b][2]*np.log(t)+c[b][3]*y-c[b][4]/2*y**-2+c[b][5]/2*y**2)/c[b][6]
    return sm
