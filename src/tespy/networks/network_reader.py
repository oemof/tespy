# -*- coding: utf-8

"""Module for loading a tespy network from saved state.

Use the method :func:`tespy.networks.network_reader.load_network` for importing
a network from a saved state.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/networks/network_reader.py

SPDX-License-Identifier: MIT
"""
import warnings

from .network import Network


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
      - Component.json
      - Connection.json
      - Bus.json
      - Network.json

    Example
    -------
    Create a network and export it. This is followed by loading the network
    with the network_reader module. All network information stored will be
    passed to a new network object. Components, connections and busses will
    be accessible by label. The following example setup is simple gas turbine
    setup with compressor, combustion chamber and turbine. The fuel is fed
    from a pipeline and throttled to the required pressure while keeping the
    temperature at a constant value.

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
    >>> nw.save('design_state.json')
    >>> _ = nw.export('exported_nwk.json')
    >>> mass_flow = round(nw.get_conn('ambient air').m.val_SI, 1)
    >>> c.set_attr(igva='var')
    >>> nw.solve('offdesign', design_path='design_state.json')
    >>> round(t.eta_s.val, 1)
    0.9
    >>> power.set_attr(P=-0.75e6)
    >>> nw.solve('offdesign', design_path='design_state.json')
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

    >>> imported_nwk = load_network('exported_nwk.json')
    >>> imported_nwk.set_attr(iterinfo=False)
    >>> imported_nwk.solve('design')
    >>> imported_nwk.lin_dep
    False
    >>> round(imported_nwk.get_conn('ambient air').m.val_SI, 1) == mass_flow
    True
    >>> round(imported_nwk.get_comp('turbine').eta_s.val, 3)
    0.9
    >>> imported_nwk.get_comp('compressor').set_attr(igva='var')
    >>> imported_nwk.solve('offdesign', design_path='design_state.json')
    >>> round(imported_nwk.get_comp('turbine').eta_s.val, 3)
    0.9
    >>> imported_nwk.busses['total power output'].set_attr(P=-0.75e6)
    >>> imported_nwk.solve('offdesign', design_path='design_state.json')
    >>> round(imported_nwk.get_comp('turbine').eta_s.val, 3) == eta_s_t
    True
    >>> round(imported_nwk.get_comp('compressor').igva.val, 3) == igva
    True
    >>> shutil.rmtree('./exported_nwk', ignore_errors=True)
    >>> shutil.rmtree('./design_state', ignore_errors=True)
    """
    msg = (
        f'The load_network method is deprecated and will be removed in the '
        'next major release of TESPy. Please use the Network class function '
        'Network.from_json() instead.'
    )
    warnings.warn(msg, FutureWarning)
    return Network.from_json(path)