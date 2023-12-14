# -*- coding: utf-8

"""Module for global variables used by other modules of the tespy package.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/tools/global_vars.py

SPDX-License-Identifier: MIT
"""
from collections import OrderedDict

ERR = 1e-6
molar_masses = {}
gas_constants = {}
gas_constants['uni'] = 8.314462618

component_property_data = OrderedDict({
    'Q': {
        'text': 'heat flow',
        'SI_unit': 'W',
        'units': {
            'W': 1, 'kW': 1000, 'MW': 1e6, 'GW': 1e9, 'TW': 1e12,
            'J / s': 1 ,     'J / h': 1 / 3.6e3,     'J / y': 1 / (3.6e3*24*365),
            'kJ / s': 1e3,  'kJ / h': 1e3 / 3.6e3,  'kJ / y': 1e3 / (3.6e3*24*365),
            'MJ / s': 1e6,  'MJ / h': 1e6 / 3.6e3,  'GJ / y': 1e6 / (3.6e3*24*365),
            'GJ / s': 1e9,  'GJ / h': 1e9 / 3.6e3,  'MJ / y': 1e9 / (3.6e3*24*365),
            'TJ / s': 1e12, 'TJ / h': 1e12 / 3.6e3, 'TJ / y': 1e12 / (3.6e3*24*365),
            'Wh / s': 3.6e3 ,   'Wh / h': 3.6e3 / 3.6e3,    'Wh / y': 3.6e3 / (3.6e3*24*365),
            'kWh / s': 3.6e6,  'kWh / h': 3.6e6 / 3.6e3,   'kWh / y': 3.6e6 / (3.6e3*24*365),
            'MWh / s': 3.6e9,  'MWh / h': 3.6e9 / 3.6e3,   'GWh / y': 3.6e9 / (3.6e3*24*365),
            'GWh / s': 3.6e12, 'GWh / h': 3.6e12 / 3.6e3,  'MWh / y': 3.6e12 / (3.6e3*24*365),
            'TWh / s': 3.6e15, 'TWh / h': 3.6e15 / 3.6e3,  'TWh / y': 3.6e15 / (3.6e3*24*365),
        },
        #'latex_eq': r'0 = \dot{m} - \dot{m}_\mathrm{spec}',
        #'documentation': {'float_fmt': '{:,.3f}'}
    },
    'kA': {
        'text': 'heat transfer conductance',
        'SI_unit': 'W / K',
        'units': {
            'W / K': 1, 'kW / K': 1000, 'MW / K': 1e6, 'GW / K': 1e9, 'TW / K': 1e12
        },
        #'latex_eq': r'0 = \dot{m} - \dot{m}_\mathrm{spec}',
        #'documentation': {'float_fmt': '{:,.3f}'}
    },
    'KPI': {
        'text': 'KPI scaling with Q',
        'SI_unit': 'Wx',
        'units': {
            'J / kg': 1 ,         'J / t': 1 / 1e3,   
            'kJ / kg': 1e3 ,     'kJ / t': 1e3 / 1e3,   
            'MJ / kg': 1e6 ,     'MJ / t': 1e6 / 1e3,   
            'GJ / kg': 1e9 ,     'GJ / t': 1e9 / 1e3,   
            'TJ / kg': 1e12 ,     'TJ / t': 1e12 / 1e3,   
           
            'Wh / kg': 3.6e3 ,      'Wh / t': 3.6e3 / 1e3, 
            'kWh / kg': 3.6e6 ,    'kWh / t': 3.6e6 / 1e3, 
            'MWh / kg': 3.6e9 ,    'MWh / t': 3.6e9 / 1e3, 
            'GWh / kg': 3.6e12 ,   'GWh / t': 3.6e12 / 1e3, 
            'TWh / kg': 3.6e15 ,   'TWh / t': 3.6e15 / 1e3, 
        },
        #'latex_eq': r'0 = \dot{m} - \dot{m}_\mathrm{spec}',
        #'documentation': {'float_fmt': '{:,.3f}'}
    },   
    'SF': {
        'text': 'species split is a mass flow',
        'SI_unit': 'kg / s',
        'units': {
            'kg / s': 1, 'kg / min': 1 / 60, 'kg / h': 1 / 3.6e3, 
            't / h': 1 / 3.6, 'g / s': 1 / 1e3, 't / y': 1e3 / (3600*24*365), 't / s': 1e3 / 1
        },
        #'latex_eq': r'0 = \dot{m} - \dot{m}_\mathrm{spec}',
        #'documentation': {'float_fmt': '{:,.3f}'}
    }         

    })

component_property_data['Q_loss'] = component_property_data['Q']
component_property_data['Q_total'] = component_property_data['Q']



fluid_property_data = OrderedDict({
    'm': {
        'text': 'mass flow',
        'SI_unit': 'kg / s',
        'units': {
            'kg / s': 1, 'kg / min': 1 / 60, 'kg / h': 1 / 3.6e3,
            't / h': 1 / 3.6, 'g / s': 1 / 1e3, 't / y': 1e3 / (3600*24*365), 't / s': 1e3 / 1
        },
        'latex_eq': r'0 = \dot{m} - \dot{m}_\mathrm{spec}',
        'documentation': {'float_fmt': '{:,.3f}'}
    },
    'v': {
        'text': 'volumetric flow',
        'SI_unit': 'm3 / s',
        'units': {
            'm3 / s': 1, 'm3 / min': 1 / 60, 'm3 / h': 1 / 3.6e3,
            'l / s': 1 / 1e3, 'l / min': 1 / 60e3, 'l / h': 1 / 3.6e6
        },
        'latex_eq': (
            r'0 = \dot{m} \cdot v \left(p,h\right)- \dot{V}_\mathrm{spec}'),
        'documentation': {'float_fmt': '{:,.3f}'}
    },
    'p': {
        'text': 'pressure',
        'SI_unit': 'Pa',
        'units': {
            'Pa': 1, 'kPa': 1e3, 'psi': 6.8948e3,
            'bar': 1e5, 'atm': 1.01325e5, 'MPa': 1e6
        },
        'latex_eq': r'0 = p - p_\mathrm{spec}',
        'documentation': {'float_fmt': '{:,.3f}'}
    },
    'h': {
        'text': 'enthalpy',
        'SI_unit': 'J / kg',
        'units': {
            'J / kg': 1, 'kJ / kg': 1e3, 'MJ / kg': 1e6,
            'cal / kg': 4.184, 'kcal / kg': 4.184e3,
            'Wh / kg': 3.6e3, 'kWh / kg': 3.6e6
        },
        'latex_eq': r'0 = h - h_\mathrm{spec}',
        'documentation': {'float_fmt': '{:,.3f}'}
    },
    'T': {
        'text': 'temperature',
        'SI_unit': 'K',
        'units': {
            'K': [0, 1], 'R': [0, 5 / 9],
            'C': [273.15, 1], 'F': [459.67, 5 / 9]
        },
        'latex_eq': r'0 = T \left(p, h \right) - T_\mathrm{spec}',
        'documentation': {'float_fmt': '{:,.1f}'}
    },
    'Td_bp': {
        'text': 'temperature difference to boiling point',
        'SI_unit': 'K',
        'units': {
            'K': 1, 'R': 5 / 9, 'C': 1, 'F': 5 / 9
        },
        'latex_eq': r'0 = \Delta T_\mathrm{spec}- T_\mathrm{sat}\left(p\right)',
        'documentation': {'float_fmt': '{:,.1f}'}
    },
    'vol': {
        'text': 'specific volume',
        'SI_unit': 'm3 / kg',
        'units': {'m3 / kg': 1, 'l / kg': 1e-3},
        'latex_eq': (
            r'0 = v\left(p,h\right) \cdot \dot{m} - \dot{V}_\mathrm{spec}'),
        'documentation': {'float_fmt': '{:,.3f}'}
    },
    'x': {
        'text': 'vapor mass fraction',
        'SI_unit': '-',
        'units': {'-': 1, '%': 1e-2, 'ppm': 1e-6},
        'latex_eq': r'0 = h - h\left(p, x_\mathrm{spec}\right)',
        'documentation': {'float_fmt': '{:,.2f}'}
    },
    's': {
        'text': 'entropy',
        'SI_unit': 'J / kgK',
        'units': {'J / kgK': 1, 'kJ / kgK': 1e3, 'MJ / kgK': 1e6},
        'latex_eq': r'0 = s_\mathrm{spec} - s\left(p, h \right)',
        'documentation': {'float_fmt': '{:,.2f}'}
    }
})

combustion_gases = ['methane', 'ethane', 'propane', 'butane', 'hydrogen', 'nDodecane']
