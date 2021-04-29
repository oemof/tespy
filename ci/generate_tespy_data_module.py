# -*- coding: utf-8 -*-

import json

from matplotlib import pyplot as plt
from pkg_resources import resource_filename

import tespy


def get_char_data(filename):
    path = resource_filename('tespy.data', filename + '.json')

    with open(path) as f:
        data = json.loads(f.read())

    return data


def plot_line(component, parameter, name, data):

    char = tespy.tools.characteristics.CharLine(x=data['x'], y=data['y'])

    title = ('Characteristic line "' + name + '" for parameter "' +
             parameter + '".')
    xlabel = '$X$'
    ylabel = r'$f\left(X\right)$'

    path = component + '_' + parameter + '_' + name + '.svg'
    char.plot(path.replace(' ', '_'), title, xlabel, ylabel)


def plot_map(component, parameter, name, data):

    char = tespy.tools.characteristics.CharMap(
        x=data['x'], y=data['y'], z=data['z'])

    title = ('Characteristic line "' + name + '" for parameter "' +
             parameter + '".')
    xlabel = '$Y$'
    ylabel = r'$f\left(X,Y\right)$'

    path = component + '_' + parameter + '_' + name + '.svg'
    char.plot(path.replace(' ', '_'), title, xlabel, ylabel)


def generate_api_doc(component, parameter, name, char_type, ref):
    path = '_images/' + component + '_' + parameter + '_' + name + '.svg'
    rst = (
        '.. figure:: ' + path.replace(' ', '_') + '\n'
        '    :alt: Characteristic ' + char_type + ' "' + name +
        '" for parameter "' + parameter + '".\n'
        '    :align: center\n\n'
    )
    if ref:
        rst += '    Reference: :cite:`' + ref + '`.\n\n'
    else:
        rst += '    Reference: Generic data.\n\n'

    return rst


rst = (
    'tespy.data module\n'
    '=================\n\n'
)

rst += (
    'Module contents\n'
    '---------------\n\n'
    '.. automodule:: tespy.data\n'
    '    :members:\n'
    '    :undoc-members:\n'
    '    :show-inheritance:\n\n'
    'Default characteristics\n'
    '-----------------------\n\n'
)

rst += (
    'Characteristic lines\n'
    '^^^^^^^^^^^^^^^^^^^^\n'
)

for component, params in get_char_data('char_lines').items():
    rst += '**' + component + '**\n\n'
    for param, lines in params.items():
        for line, data in lines.items():
            plot_line(component, param, line, data)
            rst += generate_api_doc(
                component, param, line, 'line', data['ref'])

rst += (
    'Characteristic maps\n'
    '^^^^^^^^^^^^^^^^^^^\n\n'
)

for component, params in get_char_data('char_maps').items():
    rst += '**' + component + '**\n\n'
    for param, chars in params.items():
        for char, data in chars.items():
            plot_map(component, param, char, data)
            rst += generate_api_doc(component, param, char, 'map', data['ref'])

with open('tespy.data.rst', 'w') as f:
    f.write(rst)
    f.close()
