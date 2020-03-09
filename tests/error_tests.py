# -*- coding: utf-8

"""Module for testing program errors.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/error_tests.py

SPDX-License-Identifier: MIT
"""
import pytest

from tespy.connections import connection, bus, ref
from tespy.components import (basics, combustion, components, heat_exchangers,
                              nodes, piping, reactors, subsystems,
                              turbomachinery)
from tespy.networks.network_reader import load_network
from tespy.networks.networks import network
from tespy.tools.helpers import (TESPyComponentError, TESPyConnectionError,
                                 TESPyNetworkError)
from tespy.tools.data_containers import data_container, dc_cc, dc_cp, dc_flu
from tespy.tools.fluid_properties import tespy_fluid
from tespy.tools.characteristics import char_map, char_line, load_custom_char

import shutil
import csv


# %% bulk tests


class specification_error_tests:

    def setup(self):
        self.nw = network(['water', 'air'])
        self.comp = combustion.combustion_engine('combustion engine')
        self.pipe = piping.pipe('pipe')
        self.conn = connection(self.comp, 'out1', self.pipe, 'in1')
        self.bus = bus('mybus')
        self.sub = subsystems.subsystem('MySub')

    def cmp_instanciation_ValueError(self, label):
        with pytest.raises(ValueError):
            combustion.combustion_engine(label)

    def set_attr_ValueError(self, instance, **kwargs):
        with pytest.raises(ValueError):
            instance.set_attr(**kwargs)

    def create_connection_ValueError(self, pos):
        with pytest.raises(ValueError):
            if pos == 'source':
                connection(self.comp, 'out6', self.pipe, 'in1')
            elif pos == 'target':
                connection(self.comp, 'out1', self.pipe, 'in5')

    def set_attr_TypeError(self, instance, **kwargs):
        with pytest.raises(TypeError):
            instance.set_attr(**kwargs)

    def bus_add_comps_TypeError(self, c):
        with pytest.raises(TypeError):
            self.bus.add_comps(c)

    def test_create_conn_TypeError(self):
        with pytest.raises(TypeError):
            connection(self.comp, 'out1', self.conn, 'in1')

    def create_ref_TypeError(self, params):
        with pytest.raises(TypeError):
            ref(params[0], params[1], params[2])

    def test_network_bus_type(self):
        with pytest.raises(TypeError):
            self.nw.add_busses(self.conn)

    def set_attr_KeyError(self, instance, **kwargs):
        with pytest.raises(KeyError):
            instance.set_attr(**kwargs)

    def get_attr_KeyError(self, instance, key):
        with pytest.raises(KeyError):
            instance.get_attr(key)

    def test_con_TESPyConnectionError(self):
        with pytest.raises(TESPyConnectionError):
            connection(self.comp, 'out1', self.comp, 'in1')

    def test_bus_TESPyConnectionError(self):
        with pytest.raises(TESPyConnectionError):
            self.bus.add_comps(self.comp)

    def test_network_bus_duplicate(self):
        with pytest.raises(TESPyNetworkError):
            self.nw.add_busses(self.bus, self.bus)

    def test_network_buslabel_duplicate(self):
        with pytest.raises(TESPyNetworkError):
            b = bus('mybus')
            self.nw.add_busses(self.bus)
            self.nw.add_busses(b)

    def test_set_attr_errors(self):
        #
        labels = [5, 'Label,', 'Labe;l', 'Label.']
        for l in labels:
            self.cmp_instanciation_ValueError(l)

        # ValueErrors
        self.set_attr_ValueError(self.comp, offdesign=['Q'])

        self.set_attr_ValueError(self.conn, offdesign=['f'])
        self.set_attr_ValueError(self.conn, state='f')

        self.set_attr_ValueError(self.nw, m_unit='kg')
        self.set_attr_ValueError(self.nw, h_unit='kg')
        self.set_attr_ValueError(self.nw, p_unit='kg')
        self.set_attr_ValueError(self.nw, T_unit='kg')
        self.set_attr_ValueError(self.nw, v_unit='kg')

        self.create_connection_ValueError('source')
        self.create_connection_ValueError('target')

        # TypeErrors
        self.set_attr_TypeError(self.comp, P=[5])
        self.set_attr_TypeError(self.comp, tiP_char=None)
        self.set_attr_TypeError(self.comp, design='f')
        self.set_attr_TypeError(self.comp, lamb=dc_cc())
        self.set_attr_TypeError(self.comp, design_path=7)
        self.set_attr_TypeError(self.comp, local_design=5)
        self.set_attr_TypeError(self.comp, local_offdesign=5)
        self.set_attr_TypeError(self.pipe, hydro_group=5)
        self.set_attr_TypeError(self.comp, printout=5)

        self.set_attr_TypeError(self.conn, design='h')
        self.set_attr_TypeError(self.conn, fluid_balance=1)
        self.set_attr_TypeError(self.conn, h0=[4])
        self.set_attr_TypeError(self.conn, fluid=5)
        self.set_attr_TypeError(self.conn, state=5)
        self.set_attr_TypeError(self.conn, design_path=5)
        self.set_attr_TypeError(self.conn, local_design=5)
        self.set_attr_TypeError(self.conn, local_offdesign=5)
        self.set_attr_TypeError(self.conn, printout=5)

        self.set_attr_TypeError(self.nw, m_range=5)
        self.set_attr_TypeError(self.nw, p_range=5)
        self.set_attr_TypeError(self.nw, h_range=5)
        self.set_attr_TypeError(self.nw, T_range=5)
        self.set_attr_TypeError(self.nw, iterinfo=5)
        self.set_attr_TypeError(self.bus, P='some value')
        self.set_attr_TypeError(self.bus, printout=5)

        self.bus_add_comps_TypeError({'c': self.conn})
        self.bus_add_comps_TypeError({'f': self.comp})
        self.bus_add_comps_TypeError({'c': self.comp, 'char': 'Hi'})
        self.bus_add_comps_TypeError({'c': self.comp, 'p': 5})
        self.bus_add_comps_TypeError({'c': self.comp, 'P_ref': 'what'})

        self.create_ref_TypeError([self.conn, 7, 'hi'])
        self.create_ref_TypeError([self.conn, 'hi', 0])
        self.create_ref_TypeError([self.comp, 1, 0])

        # KeyErrors
        self.set_attr_KeyError(dc_cc(), x=7)
        self.set_attr_KeyError(self.comp, wow=5)
        self.set_attr_KeyError(self.conn, jey=5)
        self.set_attr_KeyError(self.bus, power_output=100000)

        self.get_attr_KeyError(self.comp, 'wow')
        self.get_attr_KeyError(self.conn, 'key')
        self.get_attr_KeyError(self.bus, 'components')
        self.get_attr_KeyError(self.nw, 'missing')
        self.get_attr_KeyError(ref(self.conn, 1, 0), 'comp')
        self.get_attr_KeyError(self.sub, 'test')
        self.get_attr_KeyError(char_line(), 'test')
        self.get_attr_KeyError(data_container(), 'somekey')
        self.get_attr_KeyError(char_map(), 'Stuff')


##############################################################################
# components


class combustion_chamber_error_tests:

    def test_combustion_chamber_missing_fuel(self):
        """
        Test no fuel in network.
        """
        nw = network(['H2O', 'N2', 'O2', 'Ar', 'CO2'])
        instance = combustion.combustion_chamber('combustion chamber')
        c1 = connection(basics.source('air'), 'out1', instance, 'in1')
        c2 = connection(basics.source('fuel'), 'out1', instance, 'in2')
        c3 = connection(instance, 'out1', basics.sink('flue gas'), 'in1')
        nw.add_conns(c1, c2, c3)
        with pytest.raises(TESPyComponentError):
            nw.solve('design', init_only=True)


class combustion_chamber_stoich_error_tests:

    def setup(self):
        self.nw = network(['TESPy::fuel', 'TESPy::fuel_fg', 'Air'])
        label = 'combustion chamber'
        self.instance = combustion.combustion_chamber_stoich(label)
        c1 = connection(basics.source('air'), 'out1', self.instance, 'in1')
        c2 = connection(basics.source('fuel'), 'out1', self.instance, 'in2')
        c3 = connection(self.instance, 'out1', basics.sink('flue gas'), 'in1')
        self.nw.add_conns(c1, c2, c3)

    def test_cc_stoich_missing_fuel(self):
        """
        Test missing fuel composition.
        """
        with pytest.raises(TESPyComponentError):
            self.nw.solve('design', init_only=True)

    def test_cc_stoich_missing_fuel_alias(self):
        """
        Test missing fuel alias.
        """
        self.instance.set_attr(fuel={'CH4': 1})
        with pytest.raises(TESPyComponentError):
            self.nw.solve('design', init_only=True)

    def test_cc_stoich_bad_fuel_alias(self):
        """
        Test bad name for fuel alias.
        """
        self.instance.set_attr(fuel={'CH4': 1}, fuel_alias='TESPy::fuel')
        with pytest.raises(TESPyComponentError):
            self.nw.solve('design', init_only=True)

    def test_cc_stoich_missing_air(self):
        """
        Test missing air composition.
        """
        self.instance.set_attr(fuel={'CH4': 1}, fuel_alias='fuel')
        with pytest.raises(TESPyComponentError):
            self.nw.solve('design', init_only=True)

    def test_cc_stoich_missing_air_alias(self):
        """
        Test missing air alias.
        """
        self.instance.set_attr(fuel={'CH4': 1}, fuel_alias='fuel',
                               air={'N2': 0.76, 'O2': 0.24})
        with pytest.raises(TESPyComponentError):
            self.nw.solve('design', init_only=True)

    def test_cc_stoich_bad_air_alias(self):
        """
        Test bad name for air alias.
        """
        self.instance.set_attr(fuel={'CH4': 1}, fuel_alias='fuel',
                               air={'N2': 0.76, 'O2': 0.24},
                               air_alias='TESPy::air')
        with pytest.raises(TESPyComponentError):
            self.nw.solve('design', init_only=True)

    def test_cc_stoich_missing_oxygen(self):
        """
        Test bad name for air alias.
        """
        self.instance.set_attr(fuel={'CH4': 1}, fuel_alias='fuel',
                               air={'N2': 1}, air_alias='myair')
        with pytest.raises(TESPyComponentError):
            self.nw.solve('design', init_only=True)


class combustion_engine_bus_error_tests:

    def setup(self):
        self.nw = network(['water', 'air'])
        self.instance = combustion.combustion_engine('combustion engine')
        self.bus = bus('power')
        self.bus.add_comps({'c': self.instance, 'p': 'Param'})

    def test_missing_bus_param_func(self):
        """
        Test missing bus parameter in bus function for combustion engine.
        """
        with pytest.raises(ValueError):
            self.instance.bus_func(self.bus.comps.loc[self.instance])

    def test_missing_bus_param_deriv(self):
        """
        Test missing bus parameter in bus derivatives for combustion engine.
        """
        # both values do not matter, but are required for the test
        self.instance.num_nw_vars = 1
        self.instance.num_vars = 1
        self.instance.inl = [connection(self.instance, 'out1',
                                        basics.sink('sink'), 'in1')]
        self.instance.inl[0].fluid = dc_flu(val={'water': 1})
        with pytest.raises(ValueError):
            self.instance.bus_deriv(self.bus.comps.loc[self.instance])


class water_electrolyzer_error_tests:

    def setup_electrolyzer_network(self):
        """
        Set up network for electrolyzer tests.
        """
        self.instance = reactors.water_electrolyzer('electrolyzer')

        fw = basics.source('feed water')
        cw_in = basics.source('cooling water')
        o2 = basics.sink('oxygen sink')
        h2 = basics.sink('hydrogen sink')
        cw_out = basics.sink('cooling water sink')

        cw_el = connection(cw_in, 'out1', self.instance, 'in1')
        el_cw = connection(self.instance, 'out1', cw_out, 'in1')

        self.nw.add_conns(cw_el, el_cw)

        fw_el = connection(fw, 'out1', self.instance, 'in2')
        el_o2 = connection(self.instance, 'out2', o2, 'in1')
        el_h2 = connection(self.instance, 'out3', h2, 'in1')

        self.nw.add_conns(fw_el, el_o2, el_h2)

    def test_missing_hydrogen_in_network(self):
        """
        Test missing hydrogen in network fluids with water electrolyzer.
        """
        self.nw = network(['H2O', 'O2'])
        self.setup_electrolyzer_network()
        with pytest.raises(TESPyComponentError):
            self.nw.solve('design')

    def test_missing_oxygen_in_network(self):
        """
        Test missing hydrogen in network fluids with water electrolyzer.
        """
        self.nw = network(['H2O', 'H2'])
        self.setup_electrolyzer_network()
        with pytest.raises(TESPyComponentError):
            self.nw.solve('design')

    def test_missing_water_in_network(self):
        """
        Test missing hydrogen in network fluids with water electrolyzer.
        """
        self.nw = network(['O2', 'H2'])
        self.setup_electrolyzer_network()
        with pytest.raises(TESPyComponentError):
            self.nw.solve('design')

    def test_wrong_bus_param_func(self):
        """
        Test missing/wrong bus parameter specification for water electrolyzer.
        """
        # this test does not need setup, since the function is called without
        # network initialisation
        self.nw = network(['H2O', 'O2', 'H2'])
        self.instance = reactors.water_electrolyzer('electrolyzer')
        some_bus = bus('some_bus')
        param = 'G'
        some_bus.add_comps({'c': self.instance, 'p': param})
        with pytest.raises(ValueError):
            self.instance.bus_func(some_bus.comps.loc[self.instance])

    def test_wrong_bus_param_func(self):
        """
        Test wrong bus parameter specification for water electrolyzer.
        """
        # this test does not need setup, since the function is called without
        # network initialisation
        self.nw = network(['H2O', 'O2', 'H2'])
        self.instance = reactors.water_electrolyzer('electrolyzer')
        # required for calling bus_deriv method without network initialisation
        self.instance.num_vars = 1
        self.instance.num_nw_fluids = 1
        self.instance.num_nw_vars = 1
        some_bus = bus('some_bus')
        param = 'G'
        some_bus.add_comps({'c': self.instance, 'p': param})
        with pytest.raises(ValueError):
            self.instance.bus_deriv(some_bus.comps.loc[self.instance])

##############################################################################
# networks


def test_network_instanciation_no_fluids():
    nw = network([])
    so = basics.source('source')
    si = basics.sink('sink')
    conn = connection(so, 'out1', si, 'in1')
    nw.add_conns(conn)
    with pytest.raises(TESPyNetworkError):
        nw.solve('design', init_only=True)


def test_network_instanciation_single_fluid():
    with pytest.raises(TypeError):
        network('water')


def test_network_add_conns():
    with pytest.raises(TypeError):
        network(['water']).add_conns(components.component('test'))


def test_network_no_connections_error():
    nw = network(['water'])
    with pytest.raises(TESPyNetworkError):
        nw.solve('design')


def test_network_connection_error_source():
    nw = network(['water'])
    source = basics.source('source')
    sink1 = basics.sink('sink1')
    sink2 = basics.sink('sink2')
    a = connection(source, 'out1', sink1, 'in1')
    b = connection(source, 'out1', sink2, 'in1')
    nw.add_conns(a, b)
    with pytest.raises(TESPyNetworkError):
        nw.check_network()


def test_network_connection_error_target():
    nw = network(['water'])
    source1 = basics.source('source1')
    source2 = basics.source('source2')
    sink = basics.sink('sink')
    a = connection(source1, 'out1', sink, 'in1')
    b = connection(source2, 'out1', sink, 'in1')
    nw.add_conns(a, b)
    with pytest.raises(TESPyNetworkError):
        nw.check_network()


def test_network_network_consistency_inlets():
    nw = network(['water'])
    merge = nodes.merge('merge')
    sink = basics.sink('label')
    a = connection(merge, 'out1', sink, 'in1')
    nw.add_conns(a)
    with pytest.raises(TESPyNetworkError):
        nw.check_network()


def test_network_network_consistency_outlets():
    nw = network(['water', 'air'])
    source = basics.source('source')
    splitter = nodes.splitter('splitter')
    a = connection(source, 'out1', splitter, 'in1')
    nw.add_conns(a)
    with pytest.raises(TESPyNetworkError):
        nw.check_network()


def test_network_component_labels():
    nw = network(['water'])
    source = basics.source('label')
    sink = basics.sink('label')
    a = connection(source, 'out1', sink, 'in1')
    nw.add_conns(a)
    with pytest.raises(TESPyNetworkError):
        nw.check_network()


def test_network_offdesign_path():
    nw = network(['water'])
    source = basics.source('source')
    sink = basics.sink('sink')
    a = connection(source, 'out1', sink, 'in1')
    nw.add_conns(a)
    with pytest.raises(TESPyNetworkError):
        nw.solve('offdesign')


def test_network_mode():
    nw = network(['water'])
    source = basics.source('source')
    sink = basics.sink('sink')
    a = connection(source, 'out1', sink, 'in1')
    nw.add_conns(a)
    with pytest.raises(ValueError):
        nw.solve('ofdesign')


def test_network_underdetermination():
    nw = network(['water'])
    source = basics.source('source')
    sink = basics.sink('sink')
    a = connection(source, 'out1', sink, 'in1', m=1)
    nw.add_conns(a)
    with pytest.raises(TESPyNetworkError):
        nw.solve('design')


def test_network_overdetermination():
    nw = network(['water'])
    source = basics.source('source')
    sink = basics.sink('sink')
    a = connection(source, 'out1', sink, 'in1', m=1, p=1e5, x=1, h=1e6,
                   fluid={'water': 1}, fluid_balance=True)
    nw.add_conns(a)
    with pytest.raises(TESPyNetworkError):
        nw.solve('design')


def test_network_linear_dependency():
    nw = network(['water'])
    source = basics.source('source')
    sink = basics.sink('sink')
    a = connection(source, 'out1', sink, 'in1', m=1, p=1e5, h=1e6, x=1)
    nw.add_conns(a)
    nw.solve('design')
    msg = ('This test must result in a linear dependency of the jacobian '
           'matrix.')
    assert nw.lin_dep is True, msg


def test_network_no_progress():
    nw = network(['water'])
    source = basics.source('source')
    pipe = piping.pipe('pipe', pr=1, Q=-100e3)
    sink = basics.sink('sink')
    a = connection(source, 'out1', pipe, 'in1', m=1, p=1e5, T=280,
                   fluid={'water': 1})
    b = connection(pipe, 'out1', sink, 'in1')
    nw.add_conns(a, b)
    nw.solve('design')
    msg = ('This test must result in a calculation making no progress, as the '
           'pipe\'s outlet enthalpy is below fluid property range.')
    assert nw.progress is False, msg


def test_network_max_iter():
    nw = network(['water'])
    source = basics.source('source')
    pipe = piping.pipe('pipe', pr=1, Q=100e3)
    sink = basics.sink('sink')
    a = connection(source, 'out1', pipe, 'in1', m=1, p=1e5, T=280,
                   fluid={'water': 1})
    b = connection(pipe, 'out1', sink, 'in1')
    nw.add_conns(a, b)
    nw.solve('design', max_iter=2)
    msg = ('This test must result in the itercount being equal to the max '
           'iter statement.')
    assert nw.max_iter == nw.iter + 1, msg

##############################################################################
# subsystems


def test_subsys_label_str():
    with pytest.raises(ValueError):
        subsystems.subsystem(5)


def test_subsys_label_forbidden():
    with pytest.raises(ValueError):
        subsystems.subsystem('label;')

##############################################################################
# characteristics


def test_char_number_of_points():
    with pytest.raises(ValueError):
        char_line(x=[0, 1, 2], y=[1, 2, 3, 4])


def test_char_map_number_of_points():
    with pytest.raises(ValueError):
        char_map(x=[0, 1, 2], y=[[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3]])


def test_char_map_number_of_dimensions():
    with pytest.raises(ValueError):
        char_map(x=[0, 1, 2], y=[[1, 2, 3, 4], [1, 2, 3, 4]])


def test_char_map_y_z_dimension_mismatch():
    with pytest.raises(ValueError):
        char_map(x=[0, 1], y=[[1, 2, 3, 4], [1, 2, 3, 4]],
                 z1=[[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4]],
                 z2=[[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4]])


def test_missing_char_line_files():
    """
    This test might not work locally.

    If you are using custom char_lines in your HOME/.tespy/data folder, this
    test will throw a KeyError instead of a FileNotFoundError.
    """
    with pytest.raises(FileNotFoundError):
        load_custom_char('stuff', char_line)


def test_missing_char_map_files():
    """
    This test might not work locally.

    If you are using custom char_maps in your HOME/.tespy/data folder, this
    test will throw a KeyError instead of a FileNotFoundError.
    """
    with pytest.raises(FileNotFoundError):
        load_custom_char('some other stuff', char_map)

##############################################################################
# tespy fluid


def test_tespy_fluid_alias_type():
    with pytest.raises(TypeError):
        tespy_fluid(5, {'water': 1}, [0, 1], [0, 1])


def test_tespy_fluid_alias_value():
    with pytest.raises(ValueError):
        tespy_fluid('IDGAS::water', {'water': 1}, [0, 1], [0, 1])


def test_compressor_missing_char_parameter():
    """
    Compressor with invalid parameter for eta_s_char function.
    """
    nw = network(['CH4'])
    so = basics.source('source')
    si = basics.sink('sink')
    instance = turbomachinery.compressor('compressor')
    c1 = connection(so, 'out1', instance, 'in1')
    c2 = connection(instance, 'out1', si, 'in1')
    nw.add_conns(c1, c2)
    instance.set_attr(eta_s_char=dc_cc(func=char_line([0, 1], [1, 2]),
                                       is_set=True, param=None))
    nw.solve('design', init_only=True)
    with pytest.raises(ValueError):
        instance.eta_s_char_func()


def test_turbine_missing_char_parameter():
    """
    Turbine with invalid parameter for eta_s_char function.
    """
    nw = network(['CH4'])
    so = basics.source('source')
    si = basics.sink('sink')
    instance = turbomachinery.turbine('turbine')
    c1 = connection(so, 'out1', instance, 'in1')
    c2 = connection(instance, 'out1', si, 'in1')
    nw.add_conns(c1, c2)
    instance.set_attr(eta_s_char=dc_cc(func=char_line([0, 1], [1, 2]),
                                       is_set=True, param=None))
    nw.solve('design', init_only=True)
    with pytest.raises(ValueError):
        instance.eta_s_char_func()
