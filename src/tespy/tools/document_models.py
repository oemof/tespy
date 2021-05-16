# -*- coding: utf-8

"""Module for helper functions used by several other modules.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/tools/document_models.py

SPDX-License-Identifier: MIT
"""
import os
import sys
from datetime import date

import CoolProp as CP
import numpy as np
import pandas as pd

from tespy.tools import helpers as hlp
from tespy.tools.data_containers import ComponentCharacteristicMaps as dc_cm
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentCharacteristics as dc_gcc
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.global_vars import fluid_property_data as fpd
from tespy.tools.logger import check_git_branch
from tespy.tools.logger import check_version


def document_model(nw, path='report', filename='report.tex', draft=True,
                   latex_body=False, include_results=True):
    """Generate LaTeX documentation for a TESPy model.

    - The documentation is stored at path/filename
    - Generated figures are stored at path/figures/
    - Disable draft mode to skip on general info at the beginning of the report

    Parameters
    ----------
    nw : tespy.networks.network.Network
        Network instance to document.

    path : str
        Folder for the documentation, default :code:`report`.

    filename : str
        Desired filename for the LaTeX document, default :code:`report.tex`.

    draft : boolean
        Add general usage information at beginning of report,
        default :code:`True`.

    latex_body : boolean
        Add a LaTeX body to enable compilation out of the box,
        default :code:`False`.

    include_results : boolean
        Include the results in the report, default :code:`True`.
    """
    # prepare filestructure
    if path[-1] != '/' and path[-1] != '\\':
        path += '/'

    path = hlp.modify_path_os(path)
    fig_path = hlp.modify_path_os(path + 'figures/')

    # creat paths, if non existent
    if not os.path.exists(path):
        os.makedirs(path)
    if not os.path.exists(fig_path):
        os.makedirs(fig_path)

    latex = document_software_info(draft, latex_body, include_results)
    latex += document_connections(nw, include_results)
    latex += document_ude(nw, path)
    latex += document_components(nw, path, include_results)
    latex += document_busses(nw, path, include_results)
    if latex_body:
        latex += r'\end{document}'

    with open(path + filename, 'w') as f:
        f.write(latex)
        f.close()


def document_software_info(draft, latex_body, include_results):
    """Get software information.

    Parameters
    ----------
    draft : boolean
        Add general usage information at beginning of report.

    latex_body : boolean
        Add a LaTeX body to enable compilation out of the box.

    Returns
    -------
    latex : str
        LaTeX code for software information.
    """
    latex = ''
    if latex_body:
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
            r'\newcommand{\bftab}{\fontseries{b}\selectfont}'
            r'\begin{document}')

    latex += r'\section*{Software Information}' + '\n\n'

    if draft:
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
            r'\item To supress these messages, call the model '
            'documentation with the keyword draft=False.\n')
        latex += r'\end{itemize}' + '\n\n'

    latex += r'\begin{table}[H]' + '\n'
    latex += r'\begin{tabular}{ll}' + '\n'
    version = check_version().replace('_', r'\_')
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
    latex += r'\textbf{Parameter highlighting}&\\' + '\n'
    latex += r'Component variables:& \textit{italic}\\' + '\n'
    if include_results:
        latex += r'Input parameter:& \textbf{bold}\\' + '\n'
        latex += r'Results:& normalfont \\' + '\n'
    latex += r'\end{tabular}' + '\n'
    latex += r'\end{table}' + '\n'

    latex += r'\newpage'
    return latex


def document_connections(nw, include_results):
    """Document connection specifications.

    Parameters
    ----------
    nw : tespy.networks.network.Network
        TESPy model.

    include_results : boolean
        Include the results in the report.

    Returns
    -------
    latex : str
        LaTeX code for all connections.
    """

    ref_data = {'m': [], 'p': [], 'h': [], 'T': []}

    cols = nw.results['Connection'].columns
    conn_data = nw.results['Connection'].copy().loc[:, ~cols.isin(nw.fluids)]
    fluid_data = nw.results['Connection'].copy().loc[:, nw.fluids]

    lables = [label.replace('_', r'\_') for label in conn_data.index]

    value_spec = nw.specifications['Connection'].copy()
    if not include_results:
        conn_data = conn_data[value_spec]
        fluid_data = fluid_data[value_spec]

    ref_spec = nw.specifications['Ref'].dropna(
        how='all').dropna(how='all', axis=1)

    # get some Connection object for equation generator
    c = nw.get_conn(value_spec.index[0])

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
    df = data_to_df(conn_data)
    if len(df) > 0:
        eqs = df[value_spec].dropna(
            how='all').dropna(how='all', axis=1).columns
        latex += document_connection_params(
            nw, df, value_spec, eqs, c, include_results)

    df = data_to_df(fluid_data)
    if len(df) > 0:
        eqs = df[value_spec].dropna(
            how='all').dropna(how='all', axis=1).columns
        latex += document_connection_fluids(
            df, value_spec, eqs, c, include_results)

    for property, data in ref_data.items():
        df = data_to_df(data)
        if len(df) > 0:
            latex += document_connection_ref(df, property, c)

    return latex


def document_connection_params(nw, df, value_spec, eqs, c, include_results, format_string='%.2f'):
    """Document parameter specification of connections.

    Parameters
    ----------
    nw : tespy.networks.network.Network
        Network object for unit information.

    df : pandas.core.frame.DataFrame
        DataFrame containing the connection parameter data.

    c : tespy.connections.connection.Connection
        Connection object, required for LaTeX equation generation.

    Returns
    -------
    latex : str
        LaTeX code for all connections.
    """
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
            if value_spec.loc[row, col] and include_results:
                df.loc[row, col] = (
                    r'\bftab %s' % format_string % df.loc[row, col])
            else:
                df.loc[row, col] = '%s' % format_string % df.loc[row, col]

        df.rename(columns={col: col_header}, inplace=True)

    num_col = len(df.columns)

    latex += create_latex_table(df, label, col_fmt='l' + num_col * 'r')
    latex += r'\subsection{Equations applied}' + '\n\n'
    latex += equations
    return latex


def document_connection_fluids(df, value_spec, eqs, c, include_results, format_string='%.2f'):
    """Document fluid specifications of connections.

    Parameters
    ----------
    df : pandas.core.frame.DataFrame
        DataFrame containing the connection fluid data.

    c : tespy.connections.connection.Connection
        Connection object, required for LaTeX equation generation.

    Returns
    -------
    latex : str
        LaTeX code for all connections.
    """
    label = 'Specified fluids'
    latex = r'\subsection{' + label + '}' + '\n\n'

    equations = ''
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
                if value_spec.loc[row, col] and include_results:
                    df.loc[row, col] = (
                        r'\bftab %s' % format_string % df.loc[row, col])
                else:
                    df.loc[row, col] = '%s' % format_string % df.loc[row, col]

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
    caption = 'Referenced values for ' + label
    latex = r'\subsection{Referenced values for ' + label + '}' + '\n\n'
    latex += create_latex_table(df, caption)

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
        Basepath of the report.

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


def document_components(nw, path, include_results):
    """Document component specifications.

    Parameters
    ----------
    nw : tespy.networks.network.Network
        TESPy model.

    path : str
        Basepath of the report.

    Returns
    -------
    latex : str
        LaTeX code for all components.
    """
    latex = ''
    for cp in nw.comps['comp_type'].unique():

        component_list = nw.comps[nw.comps['comp_type'] == cp]['object']
        latex += get_component_mandatory_constraints(cp, component_list, path)

        latex += get_component_specifications(nw, cp, path, include_results)

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


def get_component_specifications(nw, cp, path, include_results, format_string='%.2f'):
    """Get latex code for component specifications of component type cp.

    Parameters
    ----------
    cp : str
        Classname of the current class.

    component_list : pandas.core.frame.DataFrame
        DataFrame of the components of Class cp.

    path : str
        Basepath of the report.

    Returns
    -------
    latex : str
        LaTeX code for component parameter specification.
    """
    figures = []
    col_headers = {}
    equations = ''

    cp_result = nw.results[cp].copy()
    cp_spec = nw.specifications[cp]
    lables = [label.replace('_', r'\_') for label in cp_result.index]

    if not include_results:
        cp_result = cp_result[cp_spec['properties'] | cp_spec['variables']]

    cp_result = cp_result.dropna(how='all', axis=1)
    cols = cp_result.columns.tolist()

    for row in cp_result.index:
        for col in cols:
            if cp_spec['variables'].loc[row, col]:
                cp_result.loc[row, col] = (
                    r'\textit{%s}' % format_string %
                    cp_result.loc[row, col])
            elif cp_spec['properties'].loc[row, col] and include_results:
                cp_result.loc[row, col] = (
                    r'\bftab %s' % format_string %
                    cp_result.loc[row, col])
            else:
                cp_result.loc[row, col] = (
                    '%s' % format_string % cp_result.loc[row, col])

    group_data = cp_spec['groups'][cp_spec['groups']].dropna(how='all', axis=1)
    char_data = cp_spec['chars'][cp_spec['chars']].dropna(how='all', axis=1)

    specs = pd.concat(
        [cp_spec['properties'] | cp_spec['variables'],
         cp_spec['groups'], cp_spec['chars']], axis=1)

    df_data = pd.concat([cp_result, group_data, char_data], axis=1)

    for col in char_data.columns:
        for row in char_data.index:
            component = nw.get_comp(row)
            if char_data.loc[row, col]:
                data = component.get_attr(col)
                figures += [get_char_specification(component, col, data, path)]

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
                        component, element, element_data, path,
                        group=col)]

        elements = [el for el in data.elements if el in df_data.columns]
        data_dict_gcp[col] = df_data[elements]
        group_elements += data.elements

    # remove gouped parameters from main parameter list
    df_data = df_data[
        [col for col in df_data.columns if col not in group_elements]]

    if df_data.size == 0:
        return ''

    # replace column headers
    for col in df_data.columns:
        if any(specs[col]) and col not in group_elements:
            col_headers[col] = (
                col.replace('_', r'\_') +
                r' (\ref{eq:' + cp + '_' + col + '})')
            data = nw.get_comp(row).get_attr(col)
            equations += data.latex(col, **data.func_params) + '\n\n'
        else:
            col_headers[col] = col.replace('_', r'\_')

        df_data.rename(columns=col_headers, inplace=True)

    if include_results:
        latex = r'\subsubsection{Component parameters}' + '\n\n'
    else:
        latex = r'\subsubsection{Inputs specified}' + '\n\n'

    caption = 'Parameters of components of type ' + cp
    num_col = len(df_data.columns)
    latex += create_latex_table(df_data, caption, col_fmt='l' + num_col * 'r')

    # # get parameter groups tables
    for param, data in data_dict_gcp.items():
        df_data_gcp = pd.DataFrame(data)
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


def document_busses(nw, path, include_results, format_string='{:,.2f}'):
    """Document bus specifications.

    Parameters
    ----------
    nw : tespy.networks.network.Network
        TESPy model.

    path : str
        Basepath of the report.

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
    for label, b in nw.busses.items():
        if include_results:
            df = nw.results[label][
                ['component value', 'bus value', 'efficiency']].copy()
            df.loc['total'] = df.sum()
            df.loc['total', 'efficiency'] = np.nan
            df.loc['total', 'component value'] = (
                format_string.format(df.loc['total', 'component value']))
            if b.P.is_set:
                df.loc['total', 'bus value'] = (
                    r'\bftab' + format_string.format(
                        df.loc['total', 'bus value']))
            else:
                df.loc['total', 'bus value'] = (
                    format_string.format(df.loc['total', 'bus value']))
        else:
            df = pd.DataFrame(columns=['comp eq', 'bus eq', 'eta ref'])

        figures = []
        for cp in b.comps.index:
            if include_results:
                # format cols
                df.loc[cp.label, 'bus value'] = (
                    format_string.format(df.loc[cp.label, 'bus value']))
                df.loc[cp.label, 'component value'] = (
                    format_string.format(df.loc[cp.label, 'component value']))
                df.loc[cp.label, 'efficiency'] = (
                    format_string.format(df.loc[cp.label, 'efficiency']))

            cp_data = b.comps.loc[cp]
            char = cp_data['char']
            if np.all(char.y == char.y[0]):
                if include_results:
                    eta = np.nan
                else:
                    eta = format_string.format(char.y[0])
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
                    figname = path + chars_plotted[key]['path']
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
        if include_results:
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
                '{:.3f}'.format(b.P.val) + ']{W}$\n\n')

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
        df = pd.DataFrame(data)
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
    df['label'] = df.index
    df['label'] = df['label'].str.replace('_', r'\_')
    df.set_index('label', inplace=True)
    try:
        df.replace({'nan': '-'}, inplace=True)
    except TypeError:
        # dataframes with bool data only
        pass
    longtable = False
    if df.size > 60:
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


def get_parameter_specification(data, include_results):
    """Short summary.

    Parameters
    ----------
    data : tespy.tools.data_containers.DataContainer
        Datacontainer holding the parameter's specification information.

    Returns
    -------
    value : float, boolean
        Value for specified component properties, True for characteristics,
        :code:`np.nan` if parmater is not set.
    """
    if data.is_set or include_results:
        if isinstance(data, dc_cp):
            return data.val
        else:
            return True
    else:
        return np.nan


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
    if isinstance(data, dc_cc) and data.is_set:
        return create_latex_CharLine(component, param, data, path, group=group)
    elif isinstance(data, dc_cm) and data.is_set:
        return create_latex_CharMap(component, param, data, path, group=group)
