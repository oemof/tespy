#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 19:31:51 2022

@author: karim
"""

# -*- coding: utf-8 -*-
"""
@author: Mathias Hofmann, TU Berlin
"""

# Values of standard chemical exergies of selected substances, adopted from
# Ahrendts, J.: Die Exergie chemisch reaktionsfähiger Systeme, Ruhr-Universität-Bochum, Dissertation, 1974
# Ahrendts, J.: Die Exergie chemisch reaktionsfähiger Systeme, VDI-Forschungsheft 579, VDI-Verlag, Düsseldorf, 1977
# Ahrendts, J.: Reference states, Energy (5), 1980, pp. 667-677
# [kJ/mol]
#
Chem_Ex = {                         # (0) CAS-No.     (1) solid     (2) liquid    (3) gaseous   (4) at T0,p0
    '1BUTENE':                      ['106-98-9',      'NaN',        'NaN',        'NaN',         'NaN'],
    'ACETONE':                      ['67-64-1',       'NaN',        'NaN',        'NaN',         'NaN'],
    'AMMONIA':                      ['7664-41-7',     'NaN',        'NaN',        336.684,       3],  
    'ARGON':                        ['7440-37-1',     'NaN',        'NaN',        11.627,        3],
    'BENZENE':                      ['71-43-2',       'NaN',        'NaN',        'NaN',         'NaN'],
    'BENZENE(l)':                   ['71-43-2',       'NaN',        'NaN',        'NaN',         'NaN'],
    'CARBONDIOXIDE':                ['124-38-9',      'NaN',        'NaN',        14.176,        3],
    'CARBONMONOXIDE':               ['630-08-0',      'NaN',        'NaN',        269.412,       3],
    'CARBONYLSULFIDE':              ['463-58-1',      'NaN',        'NaN',        'NaN',         'NaN'],
    'CYCLOHEXANE':                  ['110-82-7',      'NaN',        'NaN',        'NaN',         'NaN'],
    'CYCLOHEXANE(l)':               ['110-82-7',      'NaN',        'NaN',        'NaN',         'NaN'],
    'CYCLOPROPANE':                 ['75-19-4',       'NaN',        'NaN',        'NaN',         'NaN'],
    'CYCLOPENTANE':                 ['287-92-3',      'NaN',        'NaN',        'NaN',         'NaN'],
    'OCTAMETHYLCYCLOTETRASILOXANE': ['556-67-2',      'NaN',        'NaN',        'NaN',         'NaN'],
    'DECAMETHYLCYCLOPENTASILOXANE': ['541-02-6',      'NaN',        'NaN',        'NaN',         'NaN'],
    'DODECAMETHYLCYCLOHEXASILOXANE':['540-97-6',      'NaN',        'NaN',        'NaN',         'NaN'],
    'DEUTERIUM':                    ['7782-39-0',     'NaN',        'NaN',        3.951,         3],
    '1,2-DICHLOROETHANE':           ['107-06-2',      'NaN',        'NaN',        'NaN',         'NaN'],
    'DEE':                          ['60-29-7',       'NaN',        'NaN',        'NaN',         'NaN'],
    'DIMETHYLCARBONATE':            ['616-38-6',      'NaN',        'NaN',        'NaN',         'NaN'],
    'DIMETHYLETHER':                ['15-10-6',       'NaN',        'NaN',        'NaN',         'NaN'],
    'ETHANE':                       ['74-84-0',       'NaN',        'NaN',        1482.033,      3],
    'ETHANOL':                      ['64-17-5',       'NaN',        1342.086,     1348.328,      2],
    'ETHYLBENZENE':                 ['100-41-4',      'NaN',        'NaN',        'NaN',         'NaN'],
    'ETHYLBENZENE(l)':              ['100-41-4',      'NaN',        'NaN',        'NaN',         'NaN'],
    'ETHYLENE':                     ['74-85-1',       'NaN',        'NaN',        'NaN',         'NaN'],
    'ETHYLENEOXIDE':                ['75-21-8',       'NaN',        'NaN',        'NaN',         'NaN'],
    'FLUORINE':                     ['7782-41-4',     'NaN',        'NaN',        'NaN',         'NaN'],
    'HFE143M':                      ['421-14-7',      'NaN',        'NaN',        'NaN',         'NaN'],
    'HEAVYWATER':                   ['7789-20-0',     'NaN',        'NaN',        'NaN',         'NaN'],
    'HEAVYWATER(l)':                ['7789-20-0',     'NaN',        'NaN',        'NaN',         'NaN'],
    'HELIUM':                       ['7440-59-7',     'NaN',        'NaN',        'NaN',         'NaN'],
    'HYDROGEN':                     ['1333-74-0',     'NaN',        'NaN',        235.249,        3],
    'HYDROGENCHLORIDE':             ['7647-01-0',     'NaN',        'NaN',        79.759,         3],
    'HYDROGENSULFIDE':              ['7783-06-4',     'NaN',        'NaN',        799.890,        3],
    'ISOBUTAN':                     ['75-28-5',       'NaN',        'NaN',        'NaN',         'NaN'],
    'ISOBUTENE':                    ['115-11-7',      'NaN',        'NaN',        'NaN',         'NaN'],
    'ISOHEXANE':                    ['107-83-5',      'NaN',        'NaN',        'NaN',         'NaN'],
    'ISOPENTANE':                   ['78-78-4',       'NaN',        'NaN',        'NaN',         'NaN'],
    'KRYPTON':                      ['7439-90-9',     'NaN',        'NaN',        'NaN',         'NaN'],
    'DECAMETHYLTETRASILOXANE':      ['141-62-8',      'NaN',        'NaN',        'NaN',         'NaN'],
    'DODECAMETHYLPENTASILOXANE':    ['141-63-9',      'NaN',        'NaN',        'NaN',         'NaN'],
    'TETRADECAMETHYLHEXASILOXANE':  ['107-52-8',      'NaN',        'NaN',        'NaN',         'NaN'],
    'OCTAMETHYLTRISILOXANE':        ['107-51-7',      'NaN',        'NaN',        'NaN',         'NaN'],
    'HEXAMETHYLDISILOXANE':         ['107-46-0',      'NaN',        'NaN',        'NaN',         'NaN'],
    'METHANE':                      ['74-82-8',       'NaN',        'NaN',        824.348,        3],
    'METHANOL':                     ['67-56-1',       'NaN',        715.069,      710.747,        2],
    'METHYLLINOLEATE':              ['112-63-0',      'NaN',        'NaN',        'NaN',         'NaN'],
    'METHYLLINOLENATE':             ['301-00-8',      'NaN',        'NaN',        'NaN',         'NaN'],
    'METHYLOLEATE':                 ['112-62-9',      'NaN',        'NaN',        'NaN',         'NaN'],
    'METHYLPALMITATE':              ['112-39-0',      'NaN',        'NaN',        'NaN',         'NaN'],
    'METHYLSTEARATE':               ['112-61-8',      'NaN',        'NaN',        'NaN',         'NaN'],
    'NEON':                         ['7440-01-9',     'NaN',        'NaN',        'NaN',         'NaN'],
    'NEOPENTANE':                   ['463-82-1',      'NaN',        'NaN',        'NaN',         'NaN'],
    'NITROGEN':                     ['7727-37-9',     'NaN',        'NaN',        0.639,         3],
    'NITROUSOXIDE':                 ['10024-97-2',    'NaN',        'NaN',        106.807,       3],
    'NOVEC649':                     ['756-13-8',      'NaN',        'NaN',        'NaN',         'NaN'],
    'ORTHODEUTERIUM':               ['7782-39-0o',    'NaN',        'NaN',        'NaN',         'NaN'],
    'ORTHOHYDROGEN':                ['1333-74-0o',    'NaN',        'NaN',        'NaN',         'NaN'],
    'OXYGEN':                       ['7782-44-7',     'NaN',        'NaN',        3.951,         3],
    'PARADEUTERIUM':                ['7782-39-0p',    'NaN',        'NaN',        'NaN',         'NaN'],
    'PARAHYDROGEN':                 ['1333-74-0p',    'NaN',        'NaN',        'NaN',         'NaN'],
    'PROPYLENE':                    ['115-07-1',      'NaN',        'NaN',        'NaN',         'NaN'],
    'PROPYNE':                      ['74-99-7',       'NaN',        'NaN',        'NaN',         'NaN'],
    'R11':                          ['75-69-4',       'NaN',        'NaN',        'NaN',         'NaN'],
# R ...
    'SULFURDIOXIDE':                ['7446-09-5',     'NaN',        'NaN',        301.939,       3],
    'SULFURHEXAFLUORIDE':           ['2551-62-4',     'NaN',        'NaN',        'NaN',         'NaN'],
    'TOLUENE':                      ['108-88-3',      'NaN',        'NaN',        'NaN',         'NaN'],
    'TOLUENE(l)':                   ['108-88-3',      'NaN',        'NaN',        'NaN',         'NaN'],
    'WATER':                        ['7732-18-5',     'NaN',        0.045,        8.636,         2],
    'XENON':                        ['7440-63-3',     'NaN',        'NaN',        'NaN',         'NaN'],
    'C2BUTENE':                     ['590-18-1',      'NaN',        'NaN',        'NaN',         'NaN'],
    'M-XYLENE':                     ['108-38-3',      'NaN',        'NaN',        'NaN',         'NaN'],
    'N-BUTANE':                     ['106-97-8',      'NaN',        'NaN',        'NaN',         'NaN'],
    'N-DECANE':                     ['124-18-5',      'NaN',        'NaN',        'NaN',         'NaN'],
    'N-DODECANE':                   ['112-40-3',      'NaN',        'NaN',        'NaN',         'NaN'],
    'N-HEPTANE':                    ['142-82-5',      'NaN',        'NaN',        'NaN',         'NaN'],
    'N-HEXANE':                     ['110-54-3',      'NaN',        'NaN',        'NaN',         'NaN'],
    'N-HEXANE(l)':                  ['110-54-3',      'NaN',        'NaN',        'NaN',         'NaN'],
    'N-NONANE':                     ['111-84-2',      'NaN',        'NaN',        'NaN',         'NaN'],
    'N-OCTANE':                     ['111-65-9',      'NaN',        'NaN',        'NaN',         'NaN'],
    'N-PENTANE':                    ['109-66-0',      'NaN',        'NaN',        'NaN',         'NaN'],
    'N-PENTANE(l)':                 ['109-66-0',      'NaN',        'NaN',        'NaN',         'NaN'],
    'N-PROPANE':                    ['74-98-6',       'NaN',        'NaN',        'NaN',         'NaN'],
    'N-UNDECANE':                   ['1120-21-4',     'NaN',        'NaN',        'NaN',         'NaN'],
    'O-XYLENE':                     ['95-47-6',       'NaN',        'NaN',        'NaN',         'NaN'],
    'P-XYLENE':                     ['106-42-3',      'NaN',        'NaN',        'NaN',         'NaN'],
    'T2BUTENE':                     ['624-64-6',      'NaN',        'NaN',        'NaN',         'NaN'],
}