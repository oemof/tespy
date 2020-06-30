# -*- coding: utf-8

"""Module for testing program errors.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_errors.py

SPDX-License-Identifier: MIT
"""
import os
import shutil

from pytest import raises

from tespy.components import basics
from tespy.components import combustion
from tespy.components import components
from tespy.components import nodes
from tespy.components import piping
from tespy.components import reactors
from tespy.components import subsystems
from tespy.components import turbomachinery
from tespy.connections import bus
from tespy.connections import connection
from tespy.connections import ref
from tespy.networks.networks import network
from tespy.tools.characteristics import char_line
from tespy.tools.characteristics import char_map
from tespy.tools.characteristics import load_custom_char
from tespy.tools.data_containers import data_container
from tespy.tools.data_containers import dc_cc
from tespy.tools.data_containers import dc_flu
from tespy.tools.fluid_properties import h_mix_pQ
from tespy.tools.fluid_properties import memorise
from tespy.tools.fluid_properties import tespy_fluid
from tespy.tools.helpers import TESPyComponentError
from tespy.tools.helpers import TESPyConnectionError
from tespy.tools.helpers import TESPyNetworkError
from tespy.tools.helpers import extend_basic_path

##############################################################################
# test errors of set_attr and get_attr methods


def get_attr_KeyError(instance, key):
    with raises(KeyError):
        instance.get_attr(key)


def set_attr_KeyError(instance, **kwargs):
    with raises(KeyError):
        instance.set_attr(**kwargs)


def set_attr_NotImplementedError(instance, **kwargs):
    with raises(NotImplementedError):
        instance.set_attr(**kwargs)


def set_attr_TypeError(instance, **kwargs):
    with raises(TypeError):
        instance.set_attr(**kwargs)


def set_attr_ValueError(instance, **kwargs):
    with raises(ValueError):
        instance.set_attr(**kwargs)


def test_set_attr_errors():
    """Test errors of set_attr methods."""
    nw = network(['water', 'air'])
    comb = combustion.combustion_engine('combustion engine')
    pipeline = piping.pipe('pipeline')
    conn = connection(comb, 'out1', pipeline, 'in1')
    mybus = bus('mybus')

    # ValueErrors
    set_attr_ValueError(comb, offdesign=['Q'])

    set_attr_ValueError(conn, offdesign=['f'])

    set_attr_ValueError(nw, m_unit='kg')
    set_attr_ValueError(nw, h_unit='kg')
    set_attr_ValueError(nw, p_unit='kg')
    set_attr_ValueError(nw, T_unit='kg')
    set_attr_ValueError(nw, v_unit='kg')
    set_attr_ValueError(conn, state=5)

    # TypeErrors
    set_attr_TypeError(comb, P=[5])
    set_attr_TypeError(comb, tiP_char=7)
    set_attr_TypeError(comb, design='f')
    set_attr_TypeError(comb, lamb=dc_cc())
    set_attr_TypeError(comb, design_path=7)
    set_attr_TypeError(comb, local_design=5)
    set_attr_TypeError(comb, local_offdesign=5)
    set_attr_TypeError(pipeline, hydro_group=5)
    set_attr_TypeError(comb, printout=5)

    set_attr_TypeError(conn, design='h')
    set_attr_TypeError(conn, fluid_balance=1)
    set_attr_TypeError(conn, h0=[4])
    set_attr_TypeError(conn, fluid=5)
    set_attr_TypeError(conn, design_path=5)
    set_attr_TypeError(conn, local_design=5)
    set_attr_TypeError(conn, local_offdesign=5)
    set_attr_TypeError(conn, printout=5)
    set_attr_TypeError(conn, label=5)
    set_attr_TypeError(conn, state='f')

    set_attr_TypeError(nw, m_range=5)
    set_attr_TypeError(nw, p_range=5)
    set_attr_TypeError(nw, h_range=5)
    set_attr_TypeError(nw, iterinfo=5)
    set_attr_TypeError(mybus, P='some value')
    set_attr_TypeError(mybus, printout=5)

    # KeyErrors
    set_attr_KeyError(dc_cc(), x=7)
    set_attr_KeyError(comb, wow=5)
    set_attr_KeyError(conn, jey=5)
    set_attr_KeyError(mybus, power_output=100000)

    # NotImplementedError
    set_attr_NotImplementedError(conn, Td_bp=ref(conn, 1, 0))


def test_get_attr_errors():
    """Test errors of get_attr methods."""
    nw = network(['water', 'air'])
    comb = combustion.combustion_engine('combustion engine')
    pipeline = piping.pipe('pipeline')
    conn = connection(comb, 'out1', pipeline, 'in1')
    mybus = bus('mybus')
    sub = subsystems.subsystem('MySub')

    get_attr_KeyError(comb, 'wow')
    get_attr_KeyError(conn, 'key')
    get_attr_KeyError(mybus, 'components')
    get_attr_KeyError(nw, 'missing')
    get_attr_KeyError(ref(conn, 1, 0), 'comp')
    get_attr_KeyError(sub, 'test')
    get_attr_KeyError(char_line(), 'test')
    get_attr_KeyError(data_container(), 'somekey')
    get_attr_KeyError(char_map(), 'Stuff')

##############################################################################
# test error in component label


def test_cmp_instanciation_ValueError():
    """Test bad label specification for component."""
    labels = [5, 'Label,', 'Labe;l', 'Label.']
    for label in labels:
        with raises(ValueError):
            combustion.combustion_engine(label)

##############################################################################
# test errors in connection classes
##############################################################################
# connection


def test_connection_creation_ValueError():
    """Test ValueErrors creating connections."""
    comb = combustion.combustion_engine('combustion engine')
    pipeline = piping.pipe('pipeline')

    with raises(ValueError):
        connection(comb, 'out6', pipeline, 'in1')

    with raises(ValueError):
        connection(comb, 'out1', pipeline, 'in5')


def test_connection_creation_TypeError():
    """Test TypeErrors creating connections."""
    comb = combustion.combustion_engine('combustion engine')
    with raises(TypeError):
        connection(comb, 'out1', 7, 'in1')


def test_connection_creation_TESPyConnectionError():
    comb = combustion.combustion_engine('combustion engine')
    with raises(TESPyConnectionError):
        connection(comb, 'out1', comb, 'in1')

##############################################################################
# ref


def create_ref_TypeError(params):
    with raises(TypeError):
        ref(params[0], params[1], params[2])


def test_ref_creation_error():
    """Test errors creating reference objects."""
    comb = combustion.combustion_engine('combustion engine')
    pipeline = piping.pipe('pipeline')
    conn = connection(comb, 'out1', pipeline, 'in1')
    create_ref_TypeError([conn, 7, 'hi'])
    create_ref_TypeError([conn, 'hi', 0])
    create_ref_TypeError([comb, 1, 0])

##############################################################################
# bus


def bus_add_comps_TypeError(b, c):
    with raises(TypeError):
        b.add_comps(c)


def test_bus_add_comps_errors():
    """Test errors adding components to busses."""
    mybus = bus('mybus')
    comb = combustion.combustion_engine('combustion engine')
    pipeline = piping.pipe('pipeline')
    conn = connection(comb, 'out1', pipeline, 'in1')
    bus_add_comps_TypeError(mybus, {'comp': conn})
    bus_add_comps_TypeError(mybus, {'f': comb})
    bus_add_comps_TypeError(mybus, {'comp': comb, 'char': 'Hi'})
    bus_add_comps_TypeError(mybus, {'comp': comb, 'param': 5})
    bus_add_comps_TypeError(mybus, {'comp': comb, 'P_ref': 'what'})
    bus_add_comps_TypeError(mybus, comb)

    with raises(ValueError):
        mybus.add_comps({'comp': comb, 'base': 5})

##############################################################################
# test errors of component classes

##############################################################################
# combustion_chamber


def test_combustion_chamber_missing_fuel():
    """Test no fuel in network."""
    nw = network(['H2O', 'N2', 'O2', 'Ar', 'CO2'])
    instance = combustion.combustion_chamber('combustion chamber')
    c1 = connection(basics.source('air'), 'out1', instance, 'in1')
    c2 = connection(basics.source('fuel'), 'out1', instance, 'in2')
    c3 = connection(instance, 'out1', basics.sink('flue gas'), 'in1')
    nw.add_conns(c1, c2, c3)
    with raises(TESPyComponentError):
        nw.solve('design', init_only=True)

##############################################################################
# combustion_chamber_stoich


class TestCombustionChamberStoichErrors:

    def setup_combustion_chamber_stoich_error_tests(self):
        self.nw = network(['fuel', 'fuel_fg', 'Air'], p_range=[1e4, 1e6])
        label = 'combustion chamber'
        self.instance = combustion.combustion_chamber_stoich(label)
        c1 = connection(basics.source('air'), 'out1', self.instance, 'in1')
        c2 = connection(basics.source('fuel'), 'out1', self.instance, 'in2')
        c3 = connection(self.instance, 'out1', basics.sink('flue gas'), 'in1')
        self.nw.add_conns(c1, c2, c3)

    def test_cc_stoich_unset_alias(self):
        """This test unsets the alias."""
        self.setup_combustion_chamber_stoich_error_tests()
        self.instance.set_attr(air_alias='some alias')
        msg = 'The air_alias has been set, is_set should be True.'
        assert self.instance.air_alias.is_set is True, msg

        self.instance.set_attr(air_alias=None)
        msg = 'The air_alias has been unset, is_set should be False.'
        assert self.instance.air_alias.is_set is False, msg

    def test_cc_stoich_missing_fuel(self):
        """Test missing fuel composition."""
        self.setup_combustion_chamber_stoich_error_tests()
        with raises(TESPyComponentError):
            self.nw.solve('design', init_only=True)

    def test_cc_stoich_missing_fuel_alias(self):
        """Test missing fuel alias."""
        self.setup_combustion_chamber_stoich_error_tests()
        self.instance.set_attr(fuel={'CH4': 1})
        with raises(TESPyComponentError):
            self.nw.solve('design', init_only=True)

    def test_cc_stoich_bad_fuel_alias(self):
        """Test bad name for fuel alias."""
        self.setup_combustion_chamber_stoich_error_tests()
        self.instance.set_attr(fuel={'CH4': 1},
                               air={'N2': 0.76, 'O2': 0.24},
                               fuel_alias='TESPy::fuel',
                               air_alias='myair')
        with raises(ValueError):
            self.nw.solve('design', init_only=True)

    def test_cc_stoich_missing_air(self):
        """Test missing air composition."""
        self.setup_combustion_chamber_stoich_error_tests()
        self.instance.set_attr(fuel={'CH4': 1}, fuel_alias='fuel')
        with raises(TESPyComponentError):
            self.nw.solve('design', init_only=True)

    def test_cc_stoich_missing_air_alias(self):
        """Test missing air alias."""
        self.setup_combustion_chamber_stoich_error_tests()
        self.instance.set_attr(fuel={'CH4': 1}, fuel_alias='fuel',
                               air={'N2': 0.76, 'O2': 0.24})
        with raises(TESPyComponentError):
            self.nw.solve('design', init_only=True)

    def test_cc_stoich_bad_air_alias(self):
        """Test bad name for air alias."""
        self.setup_combustion_chamber_stoich_error_tests()
        self.instance.set_attr(fuel={'CH4': 1}, fuel_alias='fuel',
                               air={'N2': 0.76, 'O2': 0.24},
                               air_alias='TESPy::air')
        with raises(ValueError):
            self.nw.solve('design', init_only=True)

    def test_cc_stoich_missing_oxygen(self):
        """Test bad name for air alias."""
        self.setup_combustion_chamber_stoich_error_tests()
        self.instance.set_attr(fuel={'CH4': 1}, fuel_alias='fuel',
                               air={'N2': 1}, air_alias='myair')
        with raises(TESPyComponentError):
            self.nw.solve('design', init_only=True)

        shutil.rmtree('LUT', ignore_errors=True)

##############################################################################
# combustion_engine


class TestCombustionEngineBusErrors:

    def setup(self):
        self.nw = network(['water', 'air'])
        self.instance = combustion.combustion_engine('combustion engine')
        self.bus = bus('power')
        self.bus.add_comps({'comp': self.instance, 'param': 'Param'})

    def test_missing_bus_param_func(self):
        """Test wrong/missing bus parameter in bus function."""
        with raises(ValueError):
            self.instance.bus_func(self.bus.comps.loc[self.instance])

    def test_missing_bus_param_deriv(self):
        """Test wrong/missing bus parameter in bus derivatives."""
        # both values do not matter, but are required for the test
        self.instance.num_nw_vars = 1
        self.instance.num_vars = 1
        self.instance.inl = [connection(self.instance, 'out1',
                                        basics.sink('sink'), 'in1')]
        self.instance.inl[0].fluid = dc_flu(val={'water': 1})
        with raises(ValueError):
            self.instance.bus_deriv(self.bus)

##############################################################################
# compressor


def test_compressor_missing_char_parameter():
    """Compressor with invalid parameter for eta_s_char function."""
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
    with raises(ValueError):
        instance.eta_s_char_func()

##############################################################################
# subsystems


def test_subsys_label_str():
    with raises(ValueError):
        subsystems.subsystem(5)


def test_subsys_label_forbidden():
    with raises(ValueError):
        subsystems.subsystem('label;')

##############################################################################
# turbine


def test_turbine_missing_char_parameter():
    """Turbine with invalid parameter for eta_s_char function."""
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
    with raises(ValueError):
        instance.eta_s_char_func()

##############################################################################
# water_electrolyzer


class TestWaterElectrolyzerErrors:

    def setup_electrolyzer_network(self):
        """Set up network for electrolyzer tests."""
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
        """Test missing hydrogen in network fluids with water electrolyzer."""
        self.nw = network(['H2O', 'O2'])
        self.setup_electrolyzer_network()
        with raises(TESPyComponentError):
            self.nw.solve('design')

    def test_missing_oxygen_in_network(self):
        """Test missing oxygen in network fluids with water electrolyzer."""
        self.nw = network(['H2O', 'H2'])
        self.setup_electrolyzer_network()
        with raises(TESPyComponentError):
            self.nw.solve('design')

    def test_missing_water_in_network(self):
        """Test missing water in network fluids with water electrolyzer."""
        self.nw = network(['O2', 'H2'])
        self.setup_electrolyzer_network()
        with raises(TESPyComponentError):
            self.nw.solve('design')


def test_wrong_bus_param_func():
    """Test missing/wrong bus parameter specification in equations."""
    # this test does not need setup, since the function is called without
    # network initialisation
    instance = reactors.water_electrolyzer('electrolyzer')
    some_bus = bus('some_bus')
    param = 'G'
    some_bus.add_comps({'comp': instance, 'param': param})
    with raises(ValueError):
        instance.bus_func(some_bus.comps.loc[instance])


def test_wrong_bus_param_deriv():
    """Test missing/wrong bus parameter specification in derivatives."""
    # this test does not need setup, since the function is called without
    # network initialisation
    instance = reactors.water_electrolyzer('electrolyzer')
    # required for calling bus_deriv method without network initialisation
    instance.num_vars = 1
    instance.num_nw_fluids = 1
    instance.num_nw_vars = 1
    some_bus = bus('some_bus')
    param = 'G'
    some_bus.add_comps({'comp': instance, 'param': param})
    with raises(ValueError):
        instance.bus_deriv(some_bus)


##############################################################################
# test errors of network class

class TestNetworkErrors:

    def setup(self):
        self.nw = network(['water'])

    def test_add_conns_TypeError(self):
        with raises(TypeError):
            self.nw.add_conns(components.component('test'))

    def test_no_connections_error(self):
        with raises(TESPyNetworkError):
            self.nw.solve('design')

    def test_duplicate_connection_labels(self):
        source1 = basics.source('source1')
        source2 = basics.source('source2')
        sink1 = basics.sink('sink1')
        sink2 = basics.sink('sink2')
        a = connection(source1, 'out1', sink1, 'in1', label='myconn')
        b = connection(source2, 'out1', sink2, 'in1', label='myconn')
        with raises(ValueError):
            self.nw.add_conns(a, b)

    def test_connection_error_source(self):
        source = basics.source('source')
        sink1 = basics.sink('sink1')
        sink2 = basics.sink('sink2')
        a = connection(source, 'out1', sink1, 'in1')
        b = connection(source, 'out1', sink2, 'in1')
        self.nw.add_conns(a, b)
        with raises(TESPyNetworkError):
            self.nw.check_network()

    def test_connection_error_target(self):
        source1 = basics.source('source1')
        source2 = basics.source('source2')
        sink = basics.sink('sink')
        a = connection(source1, 'out1', sink, 'in1')
        b = connection(source2, 'out1', sink, 'in1')
        self.nw.add_conns(a, b)
        with raises(TESPyNetworkError):
            self.nw.check_network()

    def test_consistency_inlets(self):
        merge = nodes.merge('merge')
        sink = basics.sink('label')
        a = connection(merge, 'out1', sink, 'in1')
        self.nw.add_conns(a)
        with raises(TESPyNetworkError):
            self.nw.check_network()

    def test_consistency_outlets(self):
        source = basics.source('source')
        splitter = nodes.splitter('splitter')
        a = connection(source, 'out1', splitter, 'in1')
        self.nw.add_conns(a)
        with raises(TESPyNetworkError):
            self.nw.check_network()

    def test_component_label_duplicates(self):
        source = basics.source('label')
        sink = basics.sink('label')
        a = connection(source, 'out1', sink, 'in1')
        self.nw.add_conns(a)
        with raises(TESPyNetworkError):
            self.nw.check_network()

    def test_missing_offdesign_path(self):
        source = basics.source('source')
        sink = basics.sink('sink')
        a = connection(source, 'out1', sink, 'in1')
        self.nw.add_conns(a)
        with raises(TESPyNetworkError):
            self.nw.solve('offdesign')

    def test_bad_mode_specification(self):
        source = basics.source('source')
        sink = basics.sink('sink')
        a = connection(source, 'out1', sink, 'in1')
        self.nw.add_conns(a)
        with raises(ValueError):
            self.nw.solve('ofdesign')

    def test_underdetermination(self):
        source = basics.source('source')
        sink = basics.sink('sink')
        a = connection(source, 'out1', sink, 'in1', m=1)
        self.nw.add_conns(a)
        with raises(TESPyNetworkError):
            self.nw.solve('design')

    def test_overdetermination(self):
        source = basics.source('source')
        sink = basics.sink('sink')
        a = connection(source, 'out1', sink, 'in1', m=1, p=1e5, x=1, h=1e6,
                       fluid={'water': 1}, fluid_balance=True)
        self.nw.add_conns(a)
        with raises(TESPyNetworkError):
            self.nw.solve('design')

    def test_add_bus_TypeError(self):
        source = basics.source('label')
        sink = basics.sink('label')
        a = connection(source, 'out1', sink, 'in1')
        with raises(TypeError):
            self.nw.add_busses(a)

    def test_bus_duplicate(self):
        with raises(TESPyNetworkError):
            b = bus('mybus')
            self.nw.add_busses(b, b)

    def test_buslabel_duplicate(self):
        with raises(TESPyNetworkError):
            a = bus('mybus')
            b = bus('mybus')
            self.nw.add_busses(a, b)


def test_network_instanciation_no_fluids():
    nw = network([])
    so = basics.source('source')
    si = basics.sink('sink')
    conn = connection(so, 'out1', si, 'in1')
    nw.add_conns(conn)
    with raises(TESPyNetworkError):
        nw.solve('design', init_only=True)


def test_network_instanciation_single_fluid():
    with raises(TypeError):
        network('water')

##############################################################################
# test errors of characteristics classes


def test_char_number_of_points():
    with raises(ValueError):
        char_line(x=[0, 1, 2], y=[1, 2, 3, 4])


def test_char_map_number_of_points():
    with raises(ValueError):
        char_map(x=[0, 1, 2], y=[[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3]])


def test_char_map_number_of_dimensions():
    with raises(ValueError):
        char_map(x=[0, 1, 2], y=[[1, 2, 3, 4], [1, 2, 3, 4]])


def test_char_map_y_z_dimension_mismatch():
    with raises(ValueError):
        char_map(x=[0, 1], y=[[1, 2, 3, 4], [1, 2, 3, 4]],
                 z1=[[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4]],
                 z2=[[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4]])


def test_missing_char_line_files():
    """Test missing files."""
    path = extend_basic_path('data')
    tmp_path = extend_basic_path('tmp_dir_for_testing')

    if os.path.exists(path):
        for f in os.listdir(path):
            shutil.copy(src=path + '/' + f, dst=tmp_path)

        shutil.rmtree(path, ignore_errors=True)

    with raises(FileNotFoundError):
        load_custom_char('stuff', char_line)

    if os.path.exists(tmp_path):
        for f in os.listdir(tmp_path):
            shutil.copy(src=tmp_path + '/' + f, dst=path)

        shutil.rmtree(tmp_path, ignore_errors=True)


def test_missing_char_map_files():
    """Test missing files."""
    path = extend_basic_path('data')
    tmp_path = extend_basic_path('tmp_dir_for_testing')

    if os.path.exists(path):
        for f in os.listdir(path):
            shutil.copy(src=path + '/' + f, dst=tmp_path)

        shutil.rmtree(path, ignore_errors=True)

    with raises(FileNotFoundError):
        load_custom_char('some other stuff', char_map)

    if os.path.exists(tmp_path):
        for f in os.listdir(tmp_path):
            shutil.copy(src=tmp_path + '/' + f, dst=path)

        shutil.rmtree(tmp_path, ignore_errors=True)

##############################################################################
# test errors of tespy fluid class


def test_tespy_fluid_alias_type():
    with raises(TypeError):
        tespy_fluid(5, {'water': 1}, [0, 1], [0, 1])


def test_tespy_fluid_alias_value():
    with raises(ValueError):
        tespy_fluid('IDGAS::water', {'water': 1}, [0, 1], [0, 1])


##############################################################################
# test errors in fluid porperty functions


def test_IF97_back_end():
    if 'water' in memorise.state.keys():
        del memorise.state['water']
    with raises(ValueError):
        memorise.add_fluids({'water': 'IF97'})


def test_h_mix_pQ_on_mixtures():
    with raises(ValueError):
        h_mix_pQ([0, 0, 0, {'O2': 0.24, 'N2': 0.76}], 0.75)
