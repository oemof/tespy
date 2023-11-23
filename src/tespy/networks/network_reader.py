# -*- coding: utf-8

"""Module for loading a tespy network from saved state.

Use the method :func:`tespy.networks.network_reader.load_network` for importing
a network from a saved state.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/networks/network_reader.py

SPDX-License-Identifier: MIT
"""
import json
import os

import pandas as pd

from tespy.components import CombustionChamber
from tespy.components import CombustionEngine
from tespy.components import Compressor
from tespy.components import Condenser
from tespy.components import CycleCloser
from tespy.components import Desuperheater
from tespy.components import DropletSeparator
from tespy.components import Drum
from tespy.components import HeatExchanger
from tespy.components import HeatExchangerSimple
from tespy.components import Merge
from tespy.components import ParabolicTrough
from tespy.components import Pipe
from tespy.components import Pump
from tespy.components import Separator
from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import SolarCollector
from tespy.components import Source
from tespy.components import Splitter
from tespy.components import SubsystemInterface
from tespy.components import Turbine
from tespy.components import Valve
from tespy.components import WaterElectrolyzer
from tespy.components.newcomponents import *
from tespy.connections import Bus
from tespy.connections import Connection
from tespy.connections import Ref
from tespy.networks.network import Network
from tespy.tools import logger
from tespy.tools.characteristics import CharLine
from tespy.tools.characteristics import CharMap
from tespy.tools.data_containers import ComponentCharacteristicMaps as dc_cm
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import DataContainer as dc
from tespy.tools.data_containers import FluidProperties as dc_prop
from tespy.tools.fluid_properties.wrappers import CoolPropWrapper
from tespy.tools.fluid_properties.wrappers import IAPWSWrapper
from tespy.tools.fluid_properties.wrappers import PyromatWrapper
from tespy.tools.helpers import modify_path_os

COMP_TARGET_CLASSES = {
    'CycleCloser': CycleCloser,
    'Sink': Sink,
    'Source': Source,
    'SubsystemInterface': SubsystemInterface,
    'CombustionChamber': CombustionChamber,
    'CombustionEngine': CombustionEngine,
    'Condenser': Condenser,
    'Desuperheater': Desuperheater,
    'HeatExchanger': HeatExchanger,
    'HeatExchangerSimple': HeatExchangerSimple,
    'SimpleHeatExchanger': SimpleHeatExchanger,
    'SolarCollector': SolarCollector,
    'ParabolicTrough': ParabolicTrough,
    'DropletSeparator': DropletSeparator,
    'Drum': Drum,
    'Merge': Merge,
    'Separator': Separator,
    'Splitter': Splitter,
    'Pipe': Pipe,
    'Valve': Valve,
    'WaterElectrolyzer': WaterElectrolyzer,
    'Compressor': Compressor,
    'Pump': Pump,
    'Turbine': Turbine,
    'MergeDeltaP' :MergeDeltaP,
    'SeparatorWithSpeciesSplits' : SeparatorWithSpeciesSplits,
    'SeparatorWithSpeciesSplitsDeltaT' : SeparatorWithSpeciesSplitsDeltaT,
    'SeparatorWithSpeciesSplitsDeltaTDeltaP' : SeparatorWithSpeciesSplitsDeltaTDeltaP,
    'SeparatorWithSpeciesSplitsDeltaP' : SeparatorWithSpeciesSplitsDeltaP,
    'SimpleHeatExchangerDeltaPLossFactor' : SimpleHeatExchangerDeltaPLossFactor,
    'SimpleHeatExchangerDeltaP' : SimpleHeatExchangerDeltaP,
    'SimpleHeatExchangerDeltaPLfKpi' : SimpleHeatExchangerDeltaPLfKpi,
    'DrierWithAir' : DrierWithAir,
}

ENGINE_TARGET_CLASSES = {
    "CoolPropWrapper": CoolPropWrapper,
    "IAPWSWrapper": IAPWSWrapper,
    "PyromatWrapper": PyromatWrapper,
}

def load_network(path):
    r"""
    Load a network from a base path.

    Parameters
    ----------
    path : str
        The path to the network data.

    Returns
    -------
    nw : tespy.networks.network.Network
        TESPy networks object.

    Note
    ----
    If you export the network structure of an existing TESPy network, it will be
    saved to the path you specified. The structure of the saved data in that
    path is the structure you need to provide in the path for loading the
    network.

    The structure of the path must be as follows:

    - Folder: path (e.g. 'mynetwork')
    - Subfolder: components ('mynetwork/components') containing
      {component_class_name}.json (e.g. HeatExchanger.json)

    - connections.json
    - busses.json
    - network.json

    Example
    -------
    Create a network and export it. This is followed by loading the network
    with the network_reader module. All network information stored will be
    passed to a new network object. Components, connections and busses will
    be accessible by label. The following example setup is simple gas turbine
    setup with compressor, combustion chamber and turbine. The fuel is fed
    from a pipeline and throttled to the required pressure while keeping the
    temperature at a constant value.

    >>> import numpy as np
    >>> from tespy.components import (Sink, Source, CombustionChamber,
    ... Compressor, Turbine, SimpleHeatExchanger)
    >>> from tespy.connections import Connection, Ref, Bus
    >>> from tespy.networks import load_network, Network
    >>> import shutil
    >>> nw = Network(p_unit='bar', T_unit='C', h_unit='kJ / kg', iterinfo=False)
    >>> air = Source('air')
    >>> f = Source('fuel')
    >>> c = Compressor('compressor')
    >>> comb = CombustionChamber('combustion')
    >>> t = Turbine('turbine')
    >>> p = SimpleHeatExchanger('fuel preheater')
    >>> si = Sink('sink')
    >>> inc = Connection(air, 'out1', c, 'in1', label='ambient air')
    >>> cc = Connection(c, 'out1', comb, 'in1')
    >>> fp = Connection(f, 'out1', p, 'in1')
    >>> pc = Connection(p, 'out1', comb, 'in2')
    >>> ct = Connection(comb, 'out1', t, 'in1')
    >>> outg = Connection(t, 'out1', si, 'in1')
    >>> nw.add_conns(inc, cc, fp, pc, ct, outg)

    Specify component and connection properties. The intlet pressure at the
    compressor and the outlet pressure after the turbine are identical. For the
    compressor, the pressure ratio and isentropic efficiency are design
    parameters. A compressor map (efficiency vs. mass flow and pressure rise
    vs. mass flow) is selected for the compressor. Fuel is Methane.

    >>> c.set_attr(pr=10, eta_s=0.88, design=['eta_s', 'pr'],
    ... offdesign=['char_map_eta_s', 'char_map_pr'])
    >>> t.set_attr(eta_s=0.9, design=['eta_s'],
    ... offdesign=['eta_s_char', 'cone'])
    >>> comb.set_attr(lamb=2)
    >>> inc.set_attr(fluid={'N2': 0.7556, 'O2': 0.2315, 'Ar': 0.0129}, T=25, p=1)
    >>> fp.set_attr(fluid={'CH4': 0.96, 'CO2': 0.04}, T=25, p=40)
    >>> pc.set_attr(T=25)
    >>> outg.set_attr(p=Ref(inc, 1, 0))
    >>> power = Bus('total power output')
    >>> power.add_comps({"comp": c, "base": "bus"}, {"comp": t})
    >>> nw.add_busses(power)

    For a stable start, we specify the fresh air mass flow.

    >>> inc.set_attr(m=3)
    >>> nw.solve('design')

    The total power output is set to 1 MW, electrical or mechanical
    efficiencies are not considered in this example. The documentation
    example in class :py:class:`tespy.connections.bus.Bus` provides more
    information on efficiencies of generators, for instance.

    >>> comb.set_attr(lamb=None)
    >>> ct.set_attr(T=1100)
    >>> inc.set_attr(m=None)
    >>> power.set_attr(P=-1e6)
    >>> nw.solve('design')
    >>> nw.lin_dep
    False
    >>> nw.save('design_state')
    >>> nw.export('exported_nwk')
    >>> mass_flow = round(nw.get_conn('ambient air').m.val_SI, 1)
    >>> c.set_attr(igva='var')
    >>> nw.solve('offdesign', design_path='design_state')
    >>> round(t.eta_s.val, 1)
    0.9
    >>> power.set_attr(P=-0.75e6)
    >>> nw.solve('offdesign', design_path='design_state')
    >>> nw.lin_dep
    False
    >>> eta_s_t = round(t.eta_s.val, 3)
    >>> igva = round(c.igva.val, 3)
    >>> eta_s_t
    0.898
    >>> igva
    20.138

    The designed network is exported to the path 'exported_nwk'. Now import the
    network and recalculate. Check if the results match with the previous
    calculation in design and offdesign case.

    >>> imported_nwk = load_network('exported_nwk')
    >>> imported_nwk.set_attr(iterinfo=False)
    >>> imported_nwk.solve('design')
    >>> imported_nwk.lin_dep
    False
    >>> round(imported_nwk.get_conn('ambient air').m.val_SI, 1) == mass_flow
    True
    >>> round(imported_nwk.get_comp('turbine').eta_s.val, 3)
    0.9
    >>> imported_nwk.get_comp('compressor').set_attr(igva='var')
    >>> imported_nwk.solve('offdesign', design_path='design_state')
    >>> round(imported_nwk.get_comp('turbine').eta_s.val, 3)
    0.9
    >>> imported_nwk.busses['total power output'].set_attr(P=-0.75e6)
    >>> imported_nwk.solve('offdesign', design_path='design_state')
    >>> round(imported_nwk.get_comp('turbine').eta_s.val, 3) == eta_s_t
    True
    >>> round(imported_nwk.get_comp('compressor').igva.val, 3) == igva
    True
    >>> shutil.rmtree('./exported_nwk', ignore_errors=True)
    >>> shutil.rmtree('./design_state', ignore_errors=True)
    """
    if path[-1] != '/' and path[-1] != '\\':
        path += '/'

    path_comps = modify_path_os(path + 'components/')
    path = modify_path_os(path)

    msg = 'Reading network data from base path ' + path + '.'
    logger.info(msg)

    # load components
    comps = {}

    files = os.listdir(path_comps)
    for f in files:
        fn = path_comps + f
        component = f.replace(".json", "")

        msg = f"Reading component data ({component}) from {fn}."
        logger.debug(msg)

        with open(path_comps + f, "r", encoding="utf-8") as c:
            data = json.loads(c.read())

        comps.update(construct_components(component, data))

    msg = 'Created network components.'
    logger.info(msg)

    # create network
    nw = construct_network(path)

    # load connections
    fn = path + 'connections.json'
    msg = f"Reading connection data from {fn}."
    logger.debug(msg)

    with open(fn, "r", encoding="utf-8") as c:
        data = json.loads(c.read())

    conns = construct_connections(data, comps)

    # add connections to network
    for c in conns.values():
        nw.add_conns(c)

    msg = 'Created connections.'
    logger.info(msg)

    # load busses
    fn = path + 'busses.json'
    if os.path.isfile(fn):

        msg = f"Reading bus data from {fn}."
        logger.debug(msg)

        with open(fn, "r", encoding="utf-8") as c:
            data = json.loads(c.read())

        busses = construct_busses(data, comps)
        # add busses to network
        for b in busses.values():
            nw.add_busses(b)

        msg = 'Created busses.'
        logger.info(msg)

    else:
        msg = 'No bus data found!'
        logger.debug(msg)

    msg = 'Created network.'
    logger.info(msg)

    nw.check_network()

    return nw


def construct_components(component, data):
    r"""
    Create TESPy component from class name and set parameters.

    Parameters
    ----------
    component : str
        Name of the component class to be constructed.

    data : dict
        Dictionary with component information.

    Returns
    -------
    dict
        Dictionary of all components of the specified type.
    """
    target_class = COMP_TARGET_CLASSES[component]
    instances = {}
    for cp, cp_data in data.items():
        instances[cp] = target_class(cp)
        for param, param_data in cp_data.items():
            container = instances[cp].get_attr(param)
            if isinstance(container, dc):
                if isinstance(container, dc_cc):
                    if 'char_func' in param_data.keys():
                        param_data["char_func"] = CharLine(**param_data["char_func"])
                elif isinstance(container, dc_cm):
                    if 'char_func' in param_data.keys():
                        param_data["char_func"] = CharMap(**param_data["char_func"])
                if isinstance(container, dc_prop):
                    param_data["val0"] = param_data["val"]
                container.set_attr(**param_data)
            else:
                instances[cp].set_attr(**{param: param_data})

    return instances


def construct_network(path):
    r"""
    Create TESPy network from the path provided by the user.

    Parameters
    ----------
    path : str
        Base-path to stored network data.

    Returns
    -------
    nw : tespy.networks.network.Network
        TESPy network object.
    """
    # read network .json-file
    with open(path + 'network.json', 'r') as f:
        data = json.loads(f.read())

    # create network object with its properties
    return Network(**data)


def construct_connections(data, comps):
    r"""
    Create TESPy connection from data in the .json-file and its parameters.

    Parameters
    ----------
    data : dict
        Dictionary with connection data.

    comps : dict
        Dictionary of constructed components.

    Returns
    -------
    dict
        Dictionary of TESPy connection objects.
    """
    conns = {}

    arglist = [
        _ for _ in data[list(data.keys())[0]]
        if _ not in ["source", "source_id", "target", "target_id", "label", "fluid"]
        and "ref" not in _
    ]
    arglist_ref = [_ for _ in data[list(data.keys())[0]] if "ref" in _]
    for label, conn in data.items():
        conns[label] = Connection(
            comps[conn["source"]], conn["source_id"],
            comps[conn["target"]], conn["target_id"],
            label=label
        )
        for arg in arglist:
            container = conns[label].get_attr(arg)
            if isinstance(container, dc):
                container.set_attr(**conn[arg])
            else:
                conns[label].set_attr(**{arg: conn[arg]})

        for f, engine in conn["fluid"]["engine"].items():
            conn["fluid"]["engine"][f] = ENGINE_TARGET_CLASSES[engine]

        conns[label].fluid.set_attr(**conn["fluid"])
        conns[label]._create_fluid_wrapper()

    for label, conn in data.items():
        for arg in arglist_ref:
            if len(conn[arg]) > 0:
                param = arg.replace("_ref", "")
                ref = Ref(conns[conn[arg]["conn"]], conn[arg]["factor"], conn[arg]["delta"])
                conns[label].set_attr(**{param: ref})

    return conns


def construct_busses(data, comps):
    r"""
    Create busses of the network.

    Parameters
    ----------
    data : dict
        Bus information from .json file.

    comps : dict
        TESPy components dictionary.

    Returns
    -------
    dict
        Dict with TESPy bus objects.
    """
    busses = {}

    for label, bus_data in data.items():
        busses[label] = Bus(label)
        busses[label].P.set_attr(**bus_data["P"])

        components = [_ for _ in bus_data if _ != "P"]
        for cp in components:
            char = CharLine(**bus_data[cp]["char"])
            component_data = {
                "comp": comps[cp], "param": bus_data[cp]["param"],
                "base": bus_data[cp]["base"], "char": char
            }
            busses[label].add_comps(component_data)

    return busses
