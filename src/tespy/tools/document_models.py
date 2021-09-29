# -*- coding: utf-8

"""Module for helper functions used by several other modules.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/tools/document_models.py

SPDX-License-Identifier: MIT
"""
import os
import sys
from collections.abc import Mapping
from copy import deepcopy
from datetime import date

import CoolProp as CP
import numpy as np
import pandas as pd

from tespy.tools import helpers as hlp
from tespy.tools.data_containers import ComponentCharacteristicMaps as dc_cm
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.global_vars import fluid_property_data as fpd
from tespy.tools.logger import check_git_branch
from tespy.tools.logger import check_version


def document_model(nw, path='report', filename='report.tex', fmt={}):
    """Generate LaTeX documentation for a TESPy model.

    - The documentation is stored at path/filename
    - Generated figures are stored at path/figures/

    Parameters
    ----------
    nw : tespy.networks.network.Network
        Network instance to document.

    path : str
        Folder for the documentation, default :code:`report`.

    filename : str
        Desired filename for the LaTeX document, default :code:`report.tex`.

    fmt : dict
        Dictionary for formatting the report, for sample see respective
        section in online documentation.
    """
    # prepare filestructure
    if path[-1] != '/' and path[-1] != '\\':
        path += '/'

    path = hlp.modify_path_os(path)
    fig_path = hlp.modify_path_os(path + 'figures/')

    # create paths, if non existent
    if not os.path.exists(path):
        os.makedirs(path)
    if not os.path.exists(fig_path):
        os.makedirs(fig_path)

    rpt = set_defaults(nw)
    rpt = merge_dicts(rpt, fmt)

    rpt['path'] = path

    latex = document_software_info(rpt)
    latex += document_connections(nw, rpt)
    latex += document_ude(nw, rpt['path'])
    latex += document_components(nw, rpt)
    latex += document_busses(nw, rpt)
    if rpt['latex_body']:
        latex += r'\end{document}'

    with open(path + filename, 'w') as f:
        f.write(latex)
        f.close()


def merge_dicts(dict1, dict2):
    """Return a new dictionary by merging two dictionaries recursively."""

    result = deepcopy(dict1)

    for key, value in dict2.items():
        if isinstance(value, Mapping):
            result[key] = merge_dicts(result.get(key, {}), value)
        else:
            result[key] = deepcopy(dict2[key])

    return result


def set_defaults(nw):
    """
    Set up defaults for report formatting.

    Parameters
    ----------
    nw : tespy.networks.network.Network
        TESPy Network instance.

    Returns
    -------
    rpt : dict
        Dictionary containting the default formatting data.
    """
    rpt = {
        'draft': True,
        'latex_body': True,
        'include_results': True,
        'Bus': {'float_fmt': '{:,.2f}'},
        'Connection': {
            key: data['documentation'] for key, data in fpd.items()}
    }

    classes = [
        nw.comps[nw.comps['comp_type'] == cp]['object'][0]
        for cp in nw.comps['comp_type'].unique()]

    for c in classes:
        rpt[c.__class__.__name__] = {'params': []}
        rpt[c.__class__.__name__].update({
            param: {'float_fmt': '{:,.2f}'}
            for param, data in c.variables.items()
            if isinstance(data, dc_cp)
        })

    rpt['Connection']['fluid'] = {
        'float_fmt': '{:.3f}', 'include_results': True}
    rpt['Connection']['params'] = ['m', 'p', 'h', 'T', 's']
    return rpt


def document_software_info(rpt):
    """Get software information.

    Parameters
    ----------
    rpt : dict
        Formatting data for the report.

    Returns
    -------
    latex : str
        LaTeX code for software information.
    """
    latex = ''
    if rpt['latex_body']:
        latex += (
            r'\documentclass[]{article}' + '\n'
            r'\usepackage{geometry}' + '\n'
            r'\geometry{a4paper, left=20mm, top=20mm,}' + '\n'
            r'\usepackage{graphicx}' + '\n'
            r'\usepackage{float}' + '\n'
            r'\usepackage{hyperref}' + '\n'
            r'\usepackage{booktabs}' + '\n'
            r'\usepackage{amsmath}' + '\n'
            r'\usepackage{units}' + '\n'
            r'\usepackage{cleveref}' + '\n\n'
            r'\usepackage{longtable}' + '\n\n'
            r'\newcommand{\iftab}{\fontshape{sl}\selectfont}' + '\n\n'
            r'\newcommand{\bftab}{\fontseries{b}\selectfont}' + '\n\n'
            r'\begin{document}' + '\n\n')

    latex += r'\section*{Software Information}' + '\n\n'

    if rpt['draft']:
        latex += r'\begin{itemize}' + '\n'
        latex += (
            r'\item Please check, whether your inputs, the equations '
            'applied and the charactersitics are displayed correctly.\n')
        latex += (
            r'\item You are welcome to send your feedback via '
            r'\url{https://github.com/oemof/tespy/issues}.' + '\n')
        latex += r'\item \LaTeX packages required are:' + '\n'
        latex += r'\begin{itemize}' + '\n'
        latex += r'\item graphicx' + '\n'
        latex += r'\item float' + '\n'
        latex += r'\item hyperref' + '\n'
        latex += r'\item booktabs' + '\n'
        latex += r'\item amsmath' + '\n'
        latex += r'\item units' + '\n'
        latex += r'\item cleveref' + '\n'
        latex += r'\item longtable' + '\n'
        latex += r'\end{itemize}' + '\n'
        latex += (
            'Additionally, you will need to make the following '
            'definitions:\n')
        latex += r'\begin{itemize}' + '\n'
        latex += r'\item \textbackslash newcommand\{\textbackslash iftab\}'
        latex += r'\{\textbackslash fontshape\{sl\}\textbackslash selectfont\}'
        latex += '\n'
        latex += r'\item \textbackslash newcommand\{\textbackslash bftab\}'
        latex += r'\{\textbackslash fontseries\{b\}\textbackslash selectfont\}'
        latex += '\n'
        latex += r'\end{itemize}' + '\n'
        latex += (
            r'\item To suppress these messages, call the model '
            'documentation with the keyword draft=False in the formatting '
            'dict.\n')
        latex += r'\end{itemize}' + '\n\n'

    latex += r'\begin{table}[H]' + '\n'
    latex += r'\begin{tabular}{ll}' + '\n'
    version = check_version().replace('_', r'\_')
    latex += r'\bftab General information&\\' + '\n'
    latex += r'& \\' + '\n'
    latex += 'TESPy Version:&' + version + r'\\' + '\n'
    try:
        git = check_git_branch().replace('_', r'\_')
    except FileNotFoundError:
        git = 'Installation from git not found'
    latex += 'Commit:&' + git + r'\\' + '\n'
    latex += 'CoolProp version:&' + CP.__version__ + r'\\' + '\n'
    latex += 'Python version:&' + sys.version + r'\\' + '\n'
    timestamp = date.today().strftime('%B %d, %Y')
    latex += 'Documentation generated:&' + timestamp + r'\\' + '\n'
    latex += r'& \\' + '\n'
    latex += r'\bftab Parameter highlighting&\\' + '\n'
    latex += r'& \\' + '\n'
    latex += r'Variable component parameters:& \iftab italic\\' + '\n'
    if rpt['include_results']:
        latex += r'Specified input parameter:& \bftab bold\\' + '\n'
        latex += r'Results of simulation:& normalfont \\' + '\n'
        latex += r'& \\' + '\n'
        latex += (
            r'\multicolumn{2}{l}{\iftab Equations are displayed for input '
            r'parameters only.}\\' + '\n')
    else:
        latex += r'Specified input parameter:& normalfont \\' + '\n'
    latex += r'\end{tabular}' + '\n'
    latex += r'\end{table}' + '\n'

    latex += r'\newpage'
    return latex


def document_connections(nw, rpt):
    """Document connection specifications.

    Parameters
    ----------
    nw : tespy.networks.network.Network
        TESPy model.

    rpt : dict
        Formatting data for the report.

    Returns
    -------
    latex : str
        LaTeX code for all connections.
    """

    ref_data = {'m': [], 'p': [], 'h': [], 'T': []}

    cols = nw.results['Connection'].columns
    conn_data = nw.results['Connection'].copy().loc[:, ~cols.isin(nw.fluids)]
    fluid_data = nw.results['Connection'].copy().loc[:, nw.fluids]

    specs = nw.specifications['Connection'].copy()
    if not rpt['include_results']:
        conn_data = conn_data[specs]
        fluid_data = fluid_data[specs]
    # it is possible to exclude fluid results
    elif not rpt['Connection']['fluid']['include_results']:
        fluid_data = fluid_data[specs]

    ref_spec = nw.specifications['Ref'].dropna(
        how='all').dropna(how='all', axis=1)

    # get some Connection object for equation generator
    c = nw.get_conn(specs.index[0])

    for c in nw.get_conn(ref_spec.index):
        for param in ref_data.keys():
            if c.get_attr(param).ref_set:
                ref_dict = {'label': c.label.replace('_', r'\_')}
                ref_dict.update(
                    {'reference':
                     c.get_attr(param).ref.obj.label.replace('_', r'\_'),
                     'factor in -': c.get_attr(param).ref.factor,
                     'delta in ' + hlp.latex_unit(
                         nw.get_attr(param + '_unit')):
                     c.get_attr(param).ref.delta})

                ref_data[param] += [ref_dict]

    latex = r'\section{Connections in ' + nw.mode + ' mode}' + '\n\n'

    # if list is empty, all parameters will be included
    if len(rpt['Connection']['params']) > 0:
        for col in conn_data.columns:
            if col not in rpt['Connection']['params'] and not any(specs[col]):
                conn_data[col] = np.nan

    df = data_to_df(conn_data)
    if len(df) > 0:
        eqs = df[specs].dropna(how='all').dropna(how='all', axis=1).columns
        latex += document_connection_params(nw, df, specs, eqs, c, rpt)

    df = data_to_df(fluid_data)
    if len(df) > 0:
        eqs = df[specs].dropna(how='all').dropna(how='all', axis=1).columns
        latex += document_connection_fluids(df, specs, eqs, c, rpt)

    for property, data in ref_data.items():
        df = data_to_df(data)
        if len(df) > 0:
            latex += document_connection_ref(df, property, c)

    return latex


def document_connection_params(nw, df, specs, eqs, c, rpt):
    """Document parameter specification of connections.

    Parameters
    ----------
    nw : tespy.networks.network.Network
        Network object for unit information.

    df : pandas.core.frame.DataFrame
        DataFrame containing the connection parameter data.

    specs : pandas.core.frame.DataFrame
        DataFrame containing information on model input specifications.

    eqs : list
        List of parameters to generate equations for.

    c : tespy.connections.connection.Connection
        Connection object, required for LaTeX equation generation.

    rpt : dict
        Formatting data for the report.

    Returns
    -------
    latex : str
        LaTeX code for all connections.
    """
    if rpt['include_results']:
        label = 'Connection specifications and results'
    else:
        label = 'Specified connection parameters'
    latex = r'\subsection{' + label + '}' + '\n\n'

    equations = ''
    for col in df.columns:
        unit = col + '_unit'
        if col == 'Td_bp':
            unit = 'T_unit'
        col_header = (
            col.replace('_', r'\_') + ' in ' +
            hlp.latex_unit(nw.get_attr(unit)))
        if col in eqs:
            col_header += (
                r' (\ref{eq:Connection_' + fpd[col]['text'] + '})')
            equations += generate_latex_eq(
                c, fpd[col]['latex_eq'], fpd[col]['text']) + '\n\n'

        for row in df.index:
            fmt = rpt['Connection'][col]['float_fmt']
            if specs.loc[row, col] and rpt['include_results']:
                df.loc[row, col] = r'\bftab ' + fmt.format(df.loc[row, col])
            else:
                df.loc[row, col] = fmt.format(df.loc[row, col])

        df.rename(columns={col: col_header}, inplace=True)

    num_col = len(df.columns)

    latex += create_latex_table(df, label, col_fmt='l' + num_col * 'r')
    latex += r'\subsection{Equations applied}' + '\n\n'
    latex += equations
    return latex


def document_connection_fluids(df, specs, eqs, c, rpt):
    """Document fluid specifications of connections.

    Parameters
    ----------
    df : pandas.core.frame.DataFrame
        DataFrame containing the connection fluid data.

    specs : pandas.core.frame.DataFrame
        DataFrame containing information on model input specifications.

    eqs : list
        List of parameters to generate equations for.

    c : tespy.connections.connection.Connection
        Connection object, required for LaTeX equation generation.

    rpt : dict
        Formatting data for the report.

    Returns
    -------
    latex : str
        LaTeX code for all connections.
    """
    label = 'Specified fluids'
    latex = r'\subsection{' + label + '}' + '\n\n'

    equations = ''
    fmt = rpt['Connection']['fluid']['float_fmt']
    for col in eqs:
        if col == 'balance':
            eq = r'0=1-\sum x_{fl}\;\forall fl\in\text{network fluids}'
            equations += generate_latex_eq(c, eq, col) + '\n\n'
        else:
            eq = (
                r'0 = x_\mathrm{' + col + r'} - x_\mathrm{' +
                col + ',spec}')
            equations += generate_latex_eq(c, eq, col) + '\n\n'

            for row in df.index:
                if specs.loc[row, col] and rpt['include_results']:
                    df.loc[row, col] = r'\bftab ' + fmt.format(
                        df.loc[row, col])
                else:
                    df.loc[row, col] = fmt.format(df.loc[row, col])

        col_header = (
            col.replace('_', r'\_') + ' ('
            r'\ref{eq:Connection_' + col + '})')
        df.rename(columns={col: col_header}, inplace=True)

    num_col = len(df.columns)

    latex += create_latex_table(df, label, col_fmt='l' + num_col * 'r')
    latex += r'\subsection{Equations applied}' + '\n\n'
    latex += equations
    return latex


def document_connection_ref(df, property, c):
    """Document referenced connection properties

    Parameters
    ----------
    df : pandas.core.frame.DataFrame
        DataFrame containing the referenced connection data.

    property : str
        Short name of specified property (:code:`'m', 'p', ...`).

    c : tespy.connections.connection.Connection
        Connection object, required for LaTeX equation generation.

    Returns
    -------
    latex : str
        LaTeX code for all connections.
    """
    label = fpd[property]['text']
    caption = 'Specified reference values for ' + label
    latex = r'\subsection{Referenced ' + label + '}' + '\n\n'
    latex += create_latex_table(df, caption, col_fmt='llrr')

    latex += r'\subsection{Equation applied}' + '\n\n'
    eq = (
        r'0 = \text{value} - \text{value}_\mathrm{ref} '
        r'\cdot \mathrm{factor} + \text{delta}')
    latex += generate_latex_eq(c, eq, 'ref') + '\n\n'
    return latex


def document_ude(nw, path):
    """Document UserDefinedEquation specifications.

    Parameters
    ----------
    nw : tespy.networks.network.Network
        TESPy model.

    path : str
        Folder for the documentation, default :code:`report`.

    Returns
    -------
    latex : str
        LaTeX code for all UserDefinedEquations.
    """
    if len(nw.user_defined_eq) == 0:
        return ''

    latex = (
        r'\section{User defined equations in ' + nw.mode + ' mode}' + '\n\n')
    for label, ude_data in nw.user_defined_eq.items():

        eq_label = (
            r'(\ref{eq:UserDefinedEquation_' + label.replace(' ', '_') + '})')
        latex += (
            r'\subsection{Equation for ``' + label + '\'\'' + eq_label +
            r'}' + '\n\n')

        latex += generate_latex_eq(
            ude_data, ude_data.latex['equation'], label.replace(' ', '_'))
        figures = []
        i = 1
        for line in ude_data.latex['lines']:
            local_path = (
                'figures/UDE_CharLine_' +
                ude_data.label.replace(' ', '_') + '_' + str(i) + '.pdf')
            figname = path + local_path
            label = 'UDE_CharLine_' + ude_data.label + '_' + str(i)
            xlabel = '$X$'
            ylabel = r'$f\left(X\right)$'
            line.plot(figname, '', xlabel, ylabel)
            figures += [create_latex_figure(
                local_path, 'CharLine ' + str(i) + ' of ' + ude_data.label +
                ' ' + eq_label, label)]
            i += 1

        i = 1
        for map in ude_data.latex['maps']:
            local_path = (
                'figures/UDE_CharMap_' +
                ude_data.label.replace(' ', '_') + '_' + str(i) + '.pdf')
            figname = path + local_path
            label = 'UDE_CharLine_' + ude_data.label + '_' + str(i)
            xlabel = '$Y$'
            ylabel = r'$f\left(Y,\vec{Y},\vec{Z}\right)$'
            map.plot(figname, '', xlabel, ylabel)
            figures += [create_latex_figure(
                local_path, 'CharMap ' + str(i) + ' of ' + ude_data.label +
                ' ' + eq_label, label)]
            i += 1

        latex += place_figures(figures)
    return latex


def document_components(nw, rpt):
    """Document component specifications.

    Parameters
    ----------
    nw : tespy.networks.network.Network
        TESPy model.

    rpt : dict
        Formatting data for the report.

    Returns
    -------
    latex : str
        LaTeX code for all components.
    """
    latex = ''
    for cp in nw.comps['comp_type'].unique():

        component_list = nw.comps[nw.comps['comp_type'] == cp]['object']
        latex += get_component_mandatory_constraints(
            cp, component_list, rpt['path'])
        latex += get_component_specifications(nw, cp, rpt)

    if latex != '':
        latex = (
            r'\section{Components in ' + nw.mode + ' mode}' + '\n\n' + latex)

    return latex


def get_component_mandatory_constraints(cp, component_list, path):
    """Get latex code for mandatory constraints of component type cp.

    Parameters
    ----------
    cp : str
        Classname of the current class.

    component_list : pandas.core.frame.DataFrame
        DataFrame of the components of Class cp.

    path : str
        Folder for the documentation, default :code:`report`.

    Returns
    -------
    latex : str
        LaTeX code for mandatory component constraints.
    """
    latex = ''

    num_mandatory_eq = 0
    mandatory_eq = ''
    figures = []
    for label, data in component_list[0].constraints.items():
        if 'char' in data:
            for component in component_list:
                local_path = (
                    'figures/' + cp + '_CharLine_' + label + '_' +
                    component.label.replace(' ', '_') + '.pdf')
                figname = path + local_path
                xlabel = r'$X$'
                ylabel = r'$f\left(X\right)$'
                component.get_attr(data['char']).char_func.plot(
                    figname, '', xlabel, ylabel)
                figures += [create_latex_figure(
                    local_path,
                    'Characteristics of ' +
                    component.label.replace('_', r'\_') +
                    r' (eq. \ref{eq:' + cp + '_' + label + '})',
                    'CharLine_' + label + '_' + component.label)]
        mandatory_eq += data['latex'](label) + '\n\n'
        num_mandatory_eq += 1

    if num_mandatory_eq > 0:
        latex += r'\subsection{Components of type ' + cp + '}\n\n'
        latex += r'\subsubsection{Mandatory constraints}' + '\n\n'
        latex += mandatory_eq
        latex += place_figures(figures)

    return latex


def get_component_specifications(nw, cp, rpt):
    """Get latex code for component specifications of component type cp.

    Parameters
    ----------
    cp : str
        Classname of the current class.

    component_list : pandas.core.frame.DataFrame
        DataFrame of the components of Class cp.

    rpt : dict
        Formatting data for the report.

    Returns
    -------
    latex : str
        LaTeX code for component parameter specification.
    """
    figures = []
    col_headers = {}
    equations = ''

    result = nw.results[cp].copy()
    specs = nw.specifications[cp]

    if not rpt['include_results']:
        result = result[specs['properties'] | specs['variables']]
    elif len(rpt[cp]['params']) > 0:
        for col in result.columns:
            if (col not in rpt[cp]['params']
                    and not any(specs['properties'][col])
                    and not any(specs['variables'][col])):
                result[col] = np.nan

    result = result.dropna(how='all', axis=1)
    cols = result.columns.tolist()

    for col in cols:
        fmt = rpt[cp][col]['float_fmt']
        for row in result.index:
            if specs['variables'].loc[row, col]:
                result.loc[row, col] = (
                    r'\iftab ' + fmt.format(result.loc[row, col]))
            elif specs['properties'].loc[row, col] and rpt['include_results']:
                result.loc[row, col] = (
                    r'\bftab ' + fmt.format(result.loc[row, col]))
            else:
                result.loc[row, col] = fmt.format(result.loc[row, col])

    group_data = specs['groups'][specs['groups']].dropna(how='all', axis=1)
    char_data = specs['chars'][specs['chars']].dropna(how='all', axis=1)

    specs = pd.concat(
        [specs['properties'] | specs['variables'],
         specs['groups'], specs['chars']], axis=1)

    df_data = pd.concat([result, group_data, char_data], axis=1)

    for col in char_data.columns:
        for row in char_data.index:
            component = nw.get_comp(row)
            if char_data.loc[row, col]:
                data = component.get_attr(col)
                figures += [get_char_specification(
                    component, col, data, rpt['path'])]

    data_dict_gcp = {}
    group_elements = []
    for col in group_data.columns:
        for row in group_data.index:
            component = nw.get_comp(row)
            if group_data.loc[row, col]:
                data = component.get_attr(col)
                for element in data.elements:
                    element_data = component.get_attr(element)
                    figures += [get_char_specification(
                        component, element, element_data, rpt['path'],
                        group=col)]

        elements = [el for el in data.elements if el in df_data.columns]
        data_dict_gcp[col] = df_data[elements]
        group_elements += data.elements

    # remove gouped parameters from main parameter list
    df_data = df_data[
        [col for col in df_data.columns if col not in group_elements]]

    if len(df_data.index) == 0:
        return ''

    # replace column headers
    for col in df_data.columns:
        if any(specs[col]) and col not in group_elements:
            data = nw.get_comp(row).get_attr(col)
            if data.latex is None:
                df_data[col] = np.nan
            else:
                col_headers[col] = (
                    col.replace('_', r'\_') +
                    r' (\ref{eq:' + cp + '_' + col + '})')
                equations += data.latex(col, **data.func_params) + '\n\n'
        else:
            col_headers[col] = col.replace('_', r'\_')

    df_data.dropna(how='all', axis=1, inplace=True)
    df_data.rename(columns=col_headers, inplace=True)

    if rpt['include_results']:
        latex = r'\subsubsection{Specifications and results}' + '\n\n'
    else:
        latex = r'\subsubsection{Inputs specified}' + '\n\n'

    caption = 'Parameters of components of type ' + cp
    num_col = len(df_data.columns)
    latex += create_latex_table(df_data, caption, col_fmt='l' + num_col * 'r')

    # # get parameter groups tables
    for param, data in data_dict_gcp.items():
        df_data_gcp = pd.DataFrame(data, dtype='object')
        if df_data_gcp.size > 0:
            for col in df_data_gcp.columns:
                col_headers[col] = col.replace('_', r'\_')

            df_data_gcp.rename(columns=col_headers, inplace=True)
            caption = 'Parametergroup ' + param.replace('_', r'\_')
            latex += create_latex_table(df_data_gcp, caption)

    # write equations and figures of characteristics applied
    if equations != '':
        latex += r'\subsubsection{Equations applied}' + '\n\n'
        latex += equations
        latex += place_figures([fig for fig in figures if fig is not None])
    return latex


def document_busses(nw, rpt):
    """Document bus specifications.

    Parameters
    ----------
    nw : tespy.networks.network.Network
        TESPy model.

    rpt : dict
        Formatting data for the report.

    Returns
    -------
    latex : str
        LaTeX code for all busses.
    """
    if len(nw.busses) > 0:
        latex = r'\section{Busses in ' + nw.mode + ' mode}' + '\n\n'
    else:
        return ''

    chars_plotted = {}
    fmt = rpt['Bus']['float_fmt']
    for label, b in nw.busses.items():
        if rpt['include_results']:
            df = nw.results[label][
                ['component value', 'bus value', 'efficiency']].copy()
            df.loc['total'] = df.sum()
            df.loc['total', 'efficiency'] = np.nan
            df.loc['total', 'component value'] = (
                fmt.format(df.loc['total', 'component value']))
            if b.P.is_set:
                df.loc['total', 'bus value'] = (
                    r'\bftab' + fmt.format(
                        df.loc['total', 'bus value']))
            else:
                df.loc['total', 'bus value'] = (
                    fmt.format(df.loc['total', 'bus value']))
        else:
            df = pd.DataFrame(
                columns=['comp eq', 'bus eq', 'eta ref'], dtype='object')

        figures = []
        for cp in b.comps.index:
            if rpt['include_results']:
                # format cols
                df.loc[cp.label, 'bus value'] = (
                    fmt.format(df.loc[cp.label, 'bus value']))
                df.loc[cp.label, 'component value'] = (
                    fmt.format(df.loc[cp.label, 'component value']))
                df.loc[cp.label, 'efficiency'] = (
                    fmt.format(df.loc[cp.label, 'efficiency']))

            cp_data = b.comps.loc[cp]
            char = cp_data['char']
            if np.all(char.y == char.y[0]):
                if rpt['include_results']:
                    eta = np.nan
                else:
                    eta = fmt.format(char.y[0])
            else:
                key = (char, cp_data['base'])
                if key in chars_plotted:
                    eta = (
                        r'$f\left(X\right)$ (\ref{fig:' +
                        chars_plotted[key]['label'] + '})')
                else:
                    chars_plotted[key] = {
                        'path':
                            'figures/Bus_CharLine_' +
                            cp.label.replace(' ', '_') + nw.mode + '.pdf',
                        'label':
                            'Bus_CharLine_' + cp.label + nw.mode
                    }
                    figname = rpt['path'] + chars_plotted[key]['path']
                    if nw.mode == 'design':
                        xlabel = (
                            r'Energy flow ratio $X$ ($X=1$ in design mode)')
                    elif cp_data['base'] == 'bus':
                        xlabel = (
                            r'Energy flow ratio $X=\frac{\dot{E}_'
                            r'\mathrm{bus}}{\dot{E}_\mathrm{bus,design}}$')
                    else:
                        xlabel = (
                            r'Energy flow ratio $X=\frac{\dot{E}_\mathrm{'
                            r'comp}}{\dot{E}_\mathrm{comp,design}}$')
                    ylabel = r'Efficiency $\eta$'
                    char.plot(figname, '', xlabel, ylabel)
                    figures += [create_latex_figure(
                        chars_plotted[key]['path'],
                        'Bus efficiency characteristic',
                        chars_plotted[key]['label'])]
                    eta = (
                        r'$f\left(X\right)$ (\ref{fig:' +
                        chars_plotted[key]['label'] + '})')

            comp_eq = cp.bus_func_doc(cp_data)
            if comp_eq is None:
                df.loc[cp.label, 'comp eq'] = np.nan
                df.loc[cp.label, 'bus eq'] = np.nan
                df.loc[cp.label, 'eta ref'] = np.nan
                continue

            df.loc[cp.label, 'comp eq'] = '$' + comp_eq + '$'

            if cp_data['base'] == 'bus':
                eq = r'$\frac{\dot{E}_\mathrm{comp}}{\eta}$'
            else:
                eq = r'$\dot{E}_\mathrm{comp} \cdot \eta$'

            df.loc[cp.label, 'bus eq'] = eq
            df.loc[cp.label, 'eta ref'] = eta

        if all(df['comp eq'].isnull()):
            continue

        # reorder and rename columns
        if rpt['include_results']:
            col_order = [
                'comp eq', 'component value', 'bus eq', 'bus value',
                'eta ref', 'efficiency']
            df = df[col_order]

        rename_dict = {
            'component value': r'$\dot{E}_\mathrm{comp,result}$',
            'bus value': r'$\dot{E}_\mathrm{bus,result}$',
            'efficiency': r'$\eta_\mathrm{result}$',
            'comp eq': r'$\dot{E}_\mathrm{comp}$',
            'bus eq': r'$\dot{E}_\mathrm{bus}$',
            'eta ref': r'$\eta$'
        }
        df.rename(columns=rename_dict, inplace=True)

        df = data_to_df(df)

        latex += r'\subsection{Bus ``' + label + '\'\'}\n\n'
        if b.P.is_set:
            latex += (
                r'Specified total value of energy flow:'
                r' $\dot{E}_\mathrm{bus} = \unit[' +
                fmt.format(b.P.val) + ']{W}$\n\n')

            eq = r'0=\dot{E}_\mathrm{bus} -\sum_i \dot{E}_{\mathrm{bus,}i}'
            latex += generate_latex_eq(b, eq, 'energy_flow_sum') + '\n\n'
        else:
            latex += 'This bus is used for postprocessing only.\n\n'
        num_col = len(df.columns)
        latex += create_latex_table(
            df, 'Results overview for bus ' + label,
            col_fmt='l' + num_col * 'r') + '\n\n'
        latex += place_figures(figures)

    return latex


def data_to_df(data):
    """Create pandas DataFrame from list of dictionaries, remove nan columns.

    Parameters
    ----------
    data : list
        Rows for the DataFrame.

    Returns
    -------
    df : pandas.core.frame.DataFrame
        Polished DataFrame.
    """
    if not isinstance(data, pd.DataFrame):
        df = pd.DataFrame(data, dtype='object')
    else:
        df = data
    to_drop = [n for n in df.columns if n != 'label']
    df.dropna(subset=to_drop, how='all', axis=0, inplace=True)
    df.dropna(how='all', axis=1, inplace=True)
    return df


def create_latex_table(df, caption, col_fmt=None):
    """Create LaTeX table environment from DataFrame df.

    Parameters
    ----------
    df : pandas.core.frame.DataFrame
        DataFrame to export.

    caption : str
        Caption for the table.

    Returns
    -------
    latex : str
        LaTeX code for table.
    """
    df['label'] = df.index.astype('str')
    df['label'] = df['label'].str.replace('_', r'\_')
    df.set_index('label', inplace=True)
    try:
        df.replace({'nan': '-'}, inplace=True)
    except TypeError:
        # dataframes with bool data only
        pass
    longtable = False
    if len(df.index) > 60:
        longtable = True
    with pd.option_context('max_colwidth', 2000):
        latex = df.to_latex(
            index=True, escape=False, na_rep='-', column_format=col_fmt,
            longtable=longtable, caption=caption, position='H')
    return latex


def create_latex_figure(path, caption, label):
    """Create LaTeX figure environment.

    Parameters
    ----------
    path : str
        Path to the figure.

    caption : str
        Caption of the figure.

    label : str
        LaTeX label for the figure.

    Returns
    -------
    latex : str
        LaTeX code for figure.
    """
    latex = ''
    latex += r'\begin{figure}[H]\begin{center}' + '\n'
    latex += r'\includegraphics[width=\textwidth]{' + path + '}' + '\n'
    latex += r'\caption{' + caption + '}' + '\n'
    latex += r'\label{fig:' + label + '}' + '\n'
    latex += r'\end{center}\end{figure}' + '\n\n'
    return latex


def generate_latex_eq(obj, eqn, label):
    """Generate LaTeX code for equations.

    Parameters
    ----------
    obj : object
        Object equation is applied for.

    eqn : str
        LaTeX code of the equation core.

    label : str
        LaTeX label for the equation.

    Returns
    -------
    latex : str
        LaTeX code for equation.
    """
    latex = (
        r'\begin{equation}' + '\n' + r'\label{eq:' +
        obj.__class__.__name__ + '_' + label + r'}' + '\n'
    )
    latex += eqn + '\n'
    latex += r'\end{equation}'
    return latex


def create_latex_CharLine(component, param, data, path, group=None):
    """Generate image and create LaTeX code for CharLine documentation.

    Parameters
    ----------
    component : object
        Component or Bus object the characteristics are applied on.

    param : str
        Name of the parameter holding the CharLine information.

    data : tespy.tools.data_containers.ComponentCharacteristics
        DataContainer holding the CharLine information.

    path : str
        Basepath of the report.

    group : str
        Name of the group if the parameter is part of a group, else None.

    Returns
    -------
    latex : str
        LaTeX code for figure.
    """
    cp = component.__class__.__name__
    if group is None:
        group = param
    local_path = (
        'figures/' + cp + '_CharLine_' + param + '_' +
        component.label.replace(' ', '_') + '.pdf')
    figname = path + local_path
    xlabel = (
        r'$X=' + component.get_char_expr_doc(
            data.param, **data.char_params) + '$')
    ylabel = r'$f\left(X\right)$'
    data.char_func.plot(figname, '', xlabel, ylabel)
    return create_latex_figure(
        local_path,
        'Characteristics of ' + component.label.replace('_', r'\_') +
        r' (eq. \ref{eq:' + cp + '_' + group + '})',
        'CharLine_' + param + '_' + component.label)


def create_latex_CharMap(component, param, data, path, group=None):
    """Generate image and create LaTeX code for CharMap documentation.

    Parameters
    ----------
    component : object
        Component or Bus object the characteristics are applied on.

    param : str
        Name of the parameter holding the CharLine information.

    data : tespy.tools.data_containers.ComponentCharacteristicMaps
        DataContainer holding the CharMap information.

    path : str
        Basepath of the report.

    group : str
        Name of the group if the parameter is part of a group, else None.

    Returns
    -------
    latex : str
        LaTeX code for figure.
    """
    cp = component.__class__.__name__
    if group is None:
        group = param
    local_path = (
        'figures/' + cp + '_CharMap_' + param + '_' +
        component.label.replace(' ', '_') + '.pdf')
    figname = path + local_path
    xlabel = ('$Y$')
    ylabel = r'$f\left(Y,\vec{Y},\vec{Z}\right)$'
    data.char_func.plot(figname, '', xlabel, ylabel)
    return create_latex_figure(
        local_path,
        'Characteristics of ' + component.label.replace('_', r'\_') +
        r' (eq. \ref{eq:' + cp + '_' + group + '})',
        'CharMap_' + param + '_' + component.label)


def place_figures(figures):
    """Generate LaTeX code for figure placement.

    Parameters
    ----------
    figures : list
        List holding LaTeX code of individual figures to be placed in document.

    Returns
    -------
    latex : str
        LaTeX code for figure alignment.
    """
    latex = ''
    num_per_row = 2
    num_rows = len(figures) // num_per_row + 1
    for row in range(num_rows):
        for figure in figures[num_per_row * row:num_per_row * (row + 1)]:
            latex += (
                r'\begin{minipage}{' + str(1 / num_per_row) +
                r'\textwidth}' + '\n')
            latex += figure
            latex += r'\end{minipage}' + '\n'
        latex += '\n'
    return latex


def get_char_specification(component, param, data, path, group=None):
    """Get CharLine or CharMap plotting latex code.

    Parameters
    ----------
    component : object
        Component or Bus object the characteristics are applied on.

    param : str
        Name of the parameter holding the CharLine information.

    data : tespy.tools.data_containers.DataContainer
        DataContainer holding the CharMap or CharLine information.

    path : str
        Basepath of the report.

    group : str
        Name of the group if the parameter is part of a group, else None.

    Returns
    -------
    latex : str
        LaTeX code for characteristic figures.
    """
    if isinstance(data, dc_cc):
        return create_latex_CharLine(component, param, data, path, group=group)
    elif isinstance(data, dc_cm):
        return create_latex_CharMap(component, param, data, path, group=group)
