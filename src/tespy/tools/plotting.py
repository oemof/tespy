# -*- coding: utf-8

"""Module for cycle plotting data extraction.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/tools/plotting.py

SPDX-License-Identifier: MIT
"""
import numpy as np

from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.fluid_properties import s_mix_ph


def get_plotting_data(nw, connection_label):
    """
    Retrieve the process information and all state points of the fluid in a
    given part of the network identified by a connection label.

    Parameters
    ----------
    nw : tespy.networks.network.Network
        Network object.
    connection_label : str
        Label of any connection in the desired part of the network.

    Returns
    -------
    tuple
        Tuple with dictionary of process line data that are passed to the
        :code:`calc_individual_isoline` method of :code:`fluprodia` and process
        points dictionary with the state information.
    """
    # these need to be imported here as there will be circular imports
    # otherwise
    from tespy.components.basics.source import Source
    from tespy.components.heat_exchangers.base import HeatExchanger
    from tespy.components.nodes.drum import Drum


    branch_label = _get_wrapper_by_connection_label(
        nw.fluid_wrapper_branches, connection_label
    )
    connections = nw.fluid_wrapper_branches[branch_label]["connections"]
    components = nw.fluid_wrapper_branches[branch_label]["components"]

    points = {}
    processes = {}

    for component in components:
        data = component.get_plotting_data()
        if data is None:
            if isinstance(component, Source):
                processes[component.label] = None
            continue

        if isinstance(component, HeatExchanger):
            if component.inl[0] in connections and component.inl[1] in connections:
                processes[f"{component.label}_hot"] = data[1]
                processes[f"{component.label}_cold"] = data[2]

            elif component.inl[0] in connections:
                processes[component.label] = data[1]

            else:
                processes[component.label] = data[2]

        elif component.num_i == 1 and component.num_o == 1:
            processes[component.label] = data[1]

        elif isinstance(component, Drum):
            processes[f"{component.label}_1"] = data[1]
            processes[f"{component.label}_2"] = data[2]

        else:
            for key, value in data.items():
                processes[f"{component.label}_{key}"] = value

    points = {
        c.label: {
            key.replace("vol", "v"): c.get_attr(key).val
            for key in ["p", "h", "T", "s", "vol"]
        }
        for c in connections
    }

    return processes, points


def _get_wrapper_by_connection(wrappers, connection):
    """Get the wrapper branch which holds the specified connection object

    Parameters
    ----------
    wrappers : dict
        Dictionary of network fluid property wrapper branches
    connection : tespy.connections.connection.Connection
        Connection object to find wrapper branch for.

    Returns
    -------
    str
        Key of the wrapper branch
    """
    for wrapper, data in wrappers.items():
        if connection in data["connections"]:
            return wrapper


def _get_wrapper_by_connection_label(wrappers, label):
    """Get the wrapper branch which holds the specified connection label

    Parameters
    ----------
    wrappers : dict
        Dictionary of network fluid property wrapper branches
    label : str
        Connection label to find wrapper branch for.

    Returns
    -------
    str
        Key of the wrapper branch
    """
    for wrapper, data in wrappers.items():
        connection_labels = [c.label for c in data["connections"]]
        if label in connection_labels:
            return wrapper


def get_heatexchanger_secondary_Ts(nw, connection_label):
    """
    Get the **fake** process lines and points to visualize heat exchange to
    secondary fluids of a cycle.

    Temperature values are taken from the secondary side and they are matched
    to the entropy values of the primary side (the branch the connection label)
    was specified for.

    Parameters
    ----------
    nw : tespy.networks.network.Network
        Network object.
    connection_label : str
        Label of any connection in the desired part of the network.

    Returns
    -------
    tuple
        Tuple with dictionary of process line data of the secondary sides of
        heat exchangers and process points dictionary with the state
        information to plot the secondary sides in a Ts diagram
    """
    # these need to be imported here as there will be circular imports
    # otherwise
    from tespy.components.heat_exchangers.base import HeatExchanger
    from tespy.components.heat_exchangers.parallel import ParallelFlowHeatExchanger

    branch_label = _get_wrapper_by_connection_label(
        nw.fluid_wrapper_branches, connection_label
    )
    connections = nw.fluid_wrapper_branches[branch_label]["connections"]
    components = nw.fluid_wrapper_branches[branch_label]["components"]

    other_processes = {}
    other_points = {}

    for component in components:
        if not isinstance(component, HeatExchanger):
            continue

        elif component.inl[0] in connections and component.inl[1] in connections:
            continue

        elif component.inl[0] in connections:
            connection_idx = 1
        else:
            connection_idx = 0

        connection_in = component.inl[connection_idx]
        connection_out = component.outl[connection_idx]

        wrapper = _get_wrapper_by_connection(
            nw.fluid_wrapper_branches, connection_in
        )

        if wrapper not in other_processes:
            other_processes[wrapper] = {}
            other_points[wrapper] = {}

        p_range = np.linspace(connection_in.p.val_SI, connection_out.p.val_SI)
        h_range = np.linspace(connection_in.h.val_SI, connection_out.h.val_SI)
        c = connection_in

        ureg = nw.units._ureg
        T_range = ureg.Quantity(
            np.array([
                T_mix_ph(p, h, c.fluid_data, c.mixing_rule)
                for p, h in zip(p_range, h_range)
            ]),
            "K"
        ).to(c.T.unit).magnitude
        label_inflow = connection_in.label
        label_outflow = connection_out.label
        other_points[wrapper][label_inflow] = {"T": T_range[0]}
        other_points[wrapper][label_outflow] = {"T": T_range[-1]}

        # is comes from other side
        connection_in = component.inl[1 - connection_idx]
        connection_out = component.outl[1 - connection_idx]
        p_range = np.linspace(connection_in.p.val_SI, connection_out.p.val_SI)
        h_range = np.linspace(connection_in.h.val_SI, connection_out.h.val_SI)
        c = connection_in

        s_range = ureg.Quantity(
            np.array([
                s_mix_ph(p, h, c.fluid_data, c.mixing_rule)
                for p, h in zip(p_range, h_range)
            ]),
            "J/(kgK)"
        ).to(c.s.unit).magnitude

        if isinstance(component, ParallelFlowHeatExchanger):
            other_points[wrapper][label_inflow]["s"] = s_range[0]
            other_points[wrapper][label_outflow]["s"] = s_range[-1]

        else:
            other_points[wrapper][label_inflow]["s"] = s_range[-1]
            other_points[wrapper][label_outflow]["s"] = s_range[0]
            s_range = s_range[::-1]

        other_processes[wrapper][component.label] = {
            "T": T_range,
            "s": s_range
        }

    return other_processes, other_points
