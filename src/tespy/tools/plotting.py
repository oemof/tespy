def get_plotting_data(nw, wrapper_branch, include_heatexchanger_secondary=False):
    """

    Parameters
    ----------
    wrapper_branch : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """
    from tespy.components.basics.source import Source
    from tespy.components.heat_exchangers.base import HeatExchanger
    from tespy.components.nodes.drum import Drum


    connections = nw.fluid_wrapper_branches[wrapper_branch]["connections"]
    components = nw.fluid_wrapper_branches[wrapper_branch]["components"]

    points = {}
    processes = {}

    for component in components:
        data = component.get_plotting_data()
        if data is None:
            if isinstance(component, Source):
                processes[component.label] = None
                points[component.label] = component.outl[0].label
            continue

        if isinstance(component, HeatExchanger):
            if component.inl[0] in connections and component.inl[1] in connections:
                processes[f"{component.label}_hot"] = data[1]
                processes[f"{component.label}_cold"] = data[2]
                points[f"{component.label}_hot"] = component.outl[0].label
                points[f"{component.label}_cold"] = component.outl[1].label
            elif component.inl[0] in connections:
                processes[component.label] = data[1]
                points[component.label] = component.outl[0].label
            else:
                processes[component.label] = data[2]
                points[component.label] = component.outl[1].label

        elif component.num_i == 1 and component.num_o == 1:
            processes[component.label] = data[1]
            points[component.label] = component.outl[0].label

        elif isinstance(component, Drum):
            processes[f"{component.label}_1"] = data[1]
            points[f"{component.label}_1"] = component.outl[0].label
            processes[f"{component.label}_2"] = data[2]
            points[f"{component.label}_2"] = component.outl[1].label

        else:
            for key, value in data.items():
                processes[f"{component.label}_{key}"] = value
                if component.num_o == 1:
                    points[f"{component.label}_{key}"] = component.outl[0].label
                else:
                    points[f"{component.label}_{key}"] = component.outl[key - 1].label

    return processes, points
