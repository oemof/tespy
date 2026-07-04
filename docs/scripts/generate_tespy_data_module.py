# -*- coding: utf-8 -*-

import json
import os

from matplotlib import pyplot as plt

import tespy

DOCS_DIR = os.path.join(os.path.dirname(__file__), "..")
IMAGES_DIR = os.path.join(DOCS_DIR, "api", "_images", "characteristics")
API_DIR = os.path.join(DOCS_DIR, "api")

os.makedirs(IMAGES_DIR, exist_ok=True)


def get_char_data(filename):
    path = os.path.join(tespy.__datapath__, f'{filename}.json')

    with open(path) as f:
        data = json.load(f)

    return data


def plot_line(component, parameter, name, data, suffix=""):

    char = tespy.tools.characteristics.CharLine(x=data['x'], y=data['y'])

    title = f'Characteristic line "{name}" for parameter "{parameter}".'
    xlabel = '$X$'
    ylabel = r'$f\left(X\right)$'

    filename = f'{component}_{parameter}_{name}{suffix}.svg'.replace(' ', '_')
    char.plot(os.path.join(IMAGES_DIR, filename), title, xlabel, ylabel)


def plot_map(component, parameter, name, data, suffix=""):

    char = tespy.tools.characteristics.CharMap(
        x=data['x'], y=data['y'], z=data['z'])

    title = f'Characteristic line "{name}" for parameter "{parameter}".'
    xlabel = '$Y$'
    ylabel = r'$f\left(X,Y\right)$'

    filename = f'{component}_{parameter}_{name}{suffix}.svg'.replace(' ', '_')
    char.plot(os.path.join(IMAGES_DIR, filename), title, xlabel, ylabel)


def generate_api_doc(component, parameter, name, char_type, ref):
    base = f'{component}_{parameter}_{name}'.replace(' ', '_')
    alt = f'Characteristic {char_type} "{name}" for parameter "{parameter}".'
    caption = f'    Reference: :cite:`{ref}`.\n\n' if ref else '    Reference: Generic data.\n\n'

    def figure(filename, extra_class):
        return (
            f'.. figure:: /api/_images/characteristics/{filename}.svg\n'
            f'    :alt: {alt}\n'
            f'    :align: center\n'
            f'    :figclass: {extra_class}\n\n'
            + caption
        )

    return figure(base, 'only-light') + figure(f'{base}_darkmode', 'only-dark')


def generate_plots(suffix=""):
    for component, params in get_char_data('char_lines').items():
        for param, lines in params.items():
            for line, data in lines.items():
                plot_line(component, param, line, data, suffix)

    for component, params in get_char_data('char_maps').items():
        for param, chars in params.items():
            for char, data in chars.items():
                plot_map(component, param, char, data, suffix)


generate_plots()

with plt.style.context("dark_background"):
    generate_plots(suffix="_darkmode")


rst = (
    '.. _data_label:\n\n'
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
    rst += f'**{component}**\n\n'
    for param, lines in params.items():
        for line, data in lines.items():
            rst += generate_api_doc(component, param, line, 'line', data['ref'])

rst += (
    'Characteristic maps\n'
    '^^^^^^^^^^^^^^^^^^^\n\n'
)

for component, params in get_char_data('char_maps').items():
    rst += f'**{component}**\n\n'
    for param, chars in params.items():
        for char, data in chars.items():
            rst += generate_api_doc(component, param, char, 'map', data['ref'])

with open(os.path.join(API_DIR, 'data.rst'), 'w') as f:
    f.write(rst)
