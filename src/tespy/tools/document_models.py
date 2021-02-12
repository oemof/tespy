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


def document_model(nw, path='report', filename='report.tex', draft=True):
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
        default :code:`True`
    """
    # prepare filestructure
    if path[-1] != '/' and path[-1] != '\\':
        path += '/'
    path = hlp.modify_path_os(path)
    # creat path, if non existent
    if not os.path.exists(path):
        os.makedirs(path)

    fig_path = path + 'figures'
    if fig_path[-1] != '/' and fig_path[-1] != '\\':
        fig_path += '/'
    fig_path = hlp.modify_path_os(fig_path)
    # creat path, if non existent
    if not os.path.exists(fig_path):
        os.makedirs(fig_path)

    latex = document_software_info(draft)
    latex += document_connections(nw)
    latex += document_ude(nw, path)
    latex += document_components(nw, path)
    latex += document_busses(nw, path)

    with open(path + filename, 'w') as f:
        f.write(latex)
        f.close()


def document_software_info(draft):
    """Get software information.

    Parameters
    ----------
    draft : boolean
        Add general usage information at beginning of report.

    Returns
    -------
    latex : str
        LaTeX code for software information.
    """
    latex = r'\section*{Software Information}' + '\n\n'

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
    latex += r'\end{tabular}' + '\n'
    latex += r'\end{table}' + '\n'

    latex += r'\newpage'
    return latex


def document_connections(nw):
    """Document connection specifications.

    Parameters
    ----------
    nw : tespy.networks.network.Network
        TESPy model.

    Returns
    -------
    latex : str
        LaTeX code for all connections.
    """
    conn_data = []

    ref_data = {'m': [], 'p': [], 'h': [], 'T': []}
    ref_params = ['m', 'p', 'h', 'T']

    fluid_data = []

    for c in nw.conns['object']:
        data_dict = {'label': c.label.replace('_', r'\_')}
        fluid_dict = data_dict.copy()

        data_dict.update(
            {param: c.get_attr(param).val for param in fpd.keys()
             if c.get_attr(param).val_set})

        conn_data += [data_dict]

        for param in ref_params:
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

        fluid_dict.update(
            {fluid: c.fluid.val[fluid] for fluid in nw.fluids
             if c.fluid.val_set[fluid]})
        if c.fluid.balance:
            fluid_dict.update({'balance': c.fluid.balance})

        fluid_data += [fluid_dict]

    latex = r'\section{Connections in ' + nw.mode + ' mode}' + '\n\n'
    df = data_to_df(conn_data)
    if len(df) > 0:
        latex += document_connection_params(nw, df, c)

    df = data_to_df(fluid_data)
    if len(df) > 0:
        latex += document_connection_fluids(df, c)

    for property, data in ref_data.items():
        df = data_to_df(data)
        if len(df) > 0:
            latex += document_connection_ref(df, property, c)

    return latex


def document_connection_params(nw, df, c):
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
        if col == 'label':
            continue
        equations += generate_latex_eq(
            c, fpd[col]['latex_eq'], fpd[col]['text']) + '\n\n'
        col_header = (
            col.replace('_', r'\_') + ' in ' +
            hlp.latex_unit(nw.get_attr(col + '_unit')) + ' ('
            r'\ref{eq:Connection_' + fpd[col]['text'] + '})')
        df.rename(columns={col: col_header}, inplace=True)

    latex += create_latex_table(df, label)
    latex += r'\subsection{Equations applied}' + '\n\n'
    latex += equations
    return latex


def document_connection_fluids(df, c):
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
    for col in df.columns:
        if col == 'label':
            continue
        elif col == 'balance':
            eq = r'0=1-\sum x_{fl}\;\forall fl\in\text{network fluids}'
            equations += generate_latex_eq(c, eq, col) + '\n\n'
        else:
            eq = (
                r'0 = x_\mathrm{' + col + r'} - x_\mathrm{' +
                col + ',spec}')
            equations += generate_latex_eq(c, eq, col) + '\n\n'
        col_header = (
            col.replace('_', r'\_') + ' ('
            r'\ref{eq:Connection_' + col + '})')
        df.rename(columns={col: col_header}, inplace=True)

    latex += create_latex_table(df, label)
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


def document_components(nw, path):
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

        component_list = nw.comps[nw.comps['comp_type'] == cp]

        latex += get_component_mandatory_constraints(cp, component_list, path)
        latex += get_component_specifications(cp, component_list, path)

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
    for label, data in component_list['object'][0].constraints.items():
        if 'char' in data.keys():
            for component in component_list['object']:
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


def get_component_specifications(cp, component_list, path):
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
    cp_data = []
    data_dict_gcp = {}
    figures = []

    # loop through all components of type cp in component_list
    for component in component_list['object']:
        data_dict = {'label': component.label.replace('_', r'\_')}
        for param, data in component.variables.items():
            if data.latex is not None:
                data_dict[param] = get_parameter_specification(data)

                figures += [get_char_specification(
                    component, param, data, path)]

                # get group parameter specifications
                if isinstance(data, dc_gcp) and data.is_set:
                    gcp_data = {'label': component.label.replace('_', r'\_')}

                    for element in data.elements:
                        element_data = component.get_attr(element)
                        figures += [get_char_specification(
                            component, element, element_data, path,
                            group=param)]

                        gcp_data.update(
                            {element.replace('_', r'\_'):
                             get_parameter_specification(element_data)})

                    if param in data_dict_gcp.keys():
                        data_dict_gcp[param] += [gcp_data]
                    else:
                        data_dict_gcp[param] = [gcp_data]

                # get group characteristics figures
                elif isinstance(data, dc_gcc) and data.is_set:
                    for element in data.elements:
                        element_data = component.get_attr(element)
                        figures += [create_latex_CharLine(
                            component, element, element_data, path,
                            group=param)]

        cp_data += [data_dict]

    df_data = pd.DataFrame(cp_data)
    # drop empty columns
    df_data.dropna(how='all', axis=1, inplace=True)
    # if ramaining column is only the label jump to next component class
    if len(df_data.columns) < 2:
        return ''

    # get applied equations
    equations = ''
    for col in df_data.columns:
        if col == 'label':
            continue
        # can use last instance of component for this
        equations += component.get_attr(col).latex(
            col, **component.get_attr(col).func_params) + '\n\n'
        col_header = (
            col.replace('_', r'\_') + ' ('
            r'\ref{eq:' + cp + '_' + col + '})')
        df_data.rename(columns={col: col_header}, inplace=True)

    latex = r'\subsubsection{Inputs specified}' + '\n\n'

    caption = 'Parameters of components of type ' + cp
    latex += create_latex_table(df_data, caption)

    # get parameter groups tables
    for param, data in data_dict_gcp.items():
        df_data_gcp = pd.DataFrame(data)
        caption = 'Parametergroup ' + param.replace('_', r'\_')
        latex += create_latex_table(df_data_gcp, caption)

    # write equations and figures of characteristics applied
    latex += r'\subsubsection{Equations applied}' + '\n\n'
    latex += equations
    latex += place_figures([fig for fig in figures if fig is not None])
    return latex


def document_busses(nw, path):
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
        bus_data = []
        figures = []
        for cp in b.comps.index:
            bus_dict = {'label': cp.label.replace('_', r'\_')}

            cp_data = b.comps.loc[cp]
            char = cp_data['char']
            if np.all(char.y == char.y[0]):
                eta = '{:.3f}'.format(char.y[0])
            else:

                key = (char, cp_data['base'])
                if key in chars_plotted.keys():
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
                continue

            bus_dict[r'$\dot{E}_\mathrm{comp}$'] = '$' + comp_eq + '$'

            if cp_data['base'] == 'bus':
                bus_dict[r'$\dot{E}_\mathrm{bus}$'] = (
                    r'$\frac{\dot{E}_\mathrm{comp}}{\eta}$')
            else:
                bus_dict[r'$\dot{E}_\mathrm{bus}$'] = (
                    r'$\dot{E}_\mathrm{comp} \cdot \eta$')

            bus_dict[r'$\eta$'] = eta

            bus_data += [bus_dict]

        df = data_to_df(bus_data)

        if len(df) > 0:
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
            latex += create_latex_table(df, label) + '\n\n'
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
    df = pd.DataFrame(data)
    to_drop = [n for n in df.columns if n != 'label']
    df.dropna(subset=to_drop, how='all', axis=0, inplace=True)
    return df


def create_latex_table(df, caption):
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
    latex = ''
    latex += r'\begin{table}[H]\begin{center}' + '\n'
    with pd.option_context('max_colwidth', 2000):
        latex += df.to_latex(
            index=False, escape=False, na_rep='-', float_format='%.3f')
    latex += r'\caption{' + caption + '}' + '\n'
    latex += r'\end{center}\end{table}' + '\n\n'
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


def get_parameter_specification(data):
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
    if data.is_set:
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
