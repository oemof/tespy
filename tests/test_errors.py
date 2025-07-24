# -*- coding: utf-8

"""Module for testing program errors.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_errors.py

SPDX-License-Identifier: MIT
"""
import logging
import os
import shutil
import warnings

from pytest import raises
from pytest import warns

from tespy.components import CombustionChamber
from tespy.components import CombustionEngine
from tespy.components import Compressor
from tespy.components import Merge
from tespy.components import Pipe
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Splitter
from tespy.components import Subsystem
from tespy.components import Turbine
from tespy.components import WaterElectrolyzer
from tespy.components.component import Component
from tespy.connections import Bus
from tespy.connections import Connection
from tespy.connections import Ref
from tespy.networks import Network
from tespy.tools.characteristics import CharLine
from tespy.tools.characteristics import CharMap
from tespy.tools.characteristics import load_custom_char
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import DataContainer
from tespy.tools.data_containers import FluidComposition as dc_flu
from tespy.tools.fluid_properties import h_mix_pQ
from tespy.tools.helpers import TESPyComponentError
from tespy.tools.helpers import TESPyConnectionError
from tespy.tools.helpers import TESPyNetworkError
from tespy.tools.helpers import UserDefinedEquation
from tespy.tools.helpers import extend_basic_path
from tespy.tools.logger import FutureWarningHandler

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
    nw = Network()
    comb = CombustionEngine('combustion engine')
    pipeline = Pipe('pipeline')
    conn = Connection(comb, 'out1', pipeline, 'in1')
    mybus = Bus('mybus')

    # ValueErrors
    set_attr_ValueError(comb, offdesign=['Q'])

    set_attr_ValueError(conn, offdesign=['f'])

    set_attr_ValueError(nw, m_unit='kg')
    set_attr_ValueError(nw, h_unit='kg')
    set_attr_ValueError(nw, p_unit='kg')
    set_attr_ValueError(nw, T_unit='kg')
    set_attr_ValueError(nw, v_unit='kg')

    # TypeErrors
    set_attr_TypeError(comb, P=[5])
    set_attr_TypeError(comb, P=[5])
    set_attr_TypeError(comb, tiP_char=7)
    set_attr_TypeError(comb, design='f')
    set_attr_TypeError(comb, lamb=dc_cc())
    set_attr_TypeError(comb, local_design=5)
    set_attr_TypeError(comb, local_offdesign=5)
    set_attr_TypeError(comb, printout=5)

    set_attr_TypeError(conn, design='h')
    set_attr_TypeError(conn, fluid_balance=1)
    set_attr_TypeError(conn, h0=[4])
    set_attr_TypeError(conn, fluid=5)
    set_attr_TypeError(conn, local_design=5)
    set_attr_TypeError(conn, local_offdesign=5)
    set_attr_TypeError(conn, printout=5)
    set_attr_TypeError(conn, state=5)

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
    set_attr_NotImplementedError(conn, Td_bp=Ref(conn, 1, 0))
    set_attr_NotImplementedError(conn, x=Ref(conn, 1, 0))


def test_get_attr_errors():
    """Test errors of get_attr methods."""
    nw = Network()
    comb = CombustionEngine('combustion engine')
    pipeline = Pipe('pipeline')
    conn = Connection(comb, 'out1', pipeline, 'in1')
    mybus = Bus('mybus')

    get_attr_KeyError(comb, 'wow')
    get_attr_KeyError(conn, 'key')
    get_attr_KeyError(mybus, 'components')
    get_attr_KeyError(nw, 'missing')
    get_attr_KeyError(Ref(conn, 1, 0), 'comp')
    get_attr_KeyError(CharLine(), 'test')
    get_attr_KeyError(DataContainer(), 'somekey')
    get_attr_KeyError(CharMap(), 'Stuff')

##############################################################################
# test error in component label


def test_cmp_instanciation_ValueError():
    """Test bad label specification for component."""
    labels = [5, 'Label,', 'Labe;l', 'Label.']
    for label in labels:
        with raises(ValueError):
            CombustionEngine(label)

##############################################################################
# test errors in connection classes
##############################################################################
# connection


def test_Connection_creation_ValueError():
    """Test ValueErrors creating connections."""
    comb = CombustionEngine('combustion engine')
    pipeline = Pipe('pipeline')

    with raises(ValueError):
        Connection(comb, 'out6', pipeline, 'in1')

    with raises(ValueError):
        Connection(comb, 'out1', pipeline, 'in5')


def test_Connection_creation_TypeError():
    """Test TypeErrors creating connections."""
    comb = CombustionEngine('combustion engine')
    pipeline = Pipe('pipeline')
    with raises(TypeError):
        Connection(comb, 'out1', 7, 'in1')

    with raises(TypeError):
        Connection(comb, 'out1', pipeline, 'in1', label=5)


def test_Connection_creation_TESPyConnectionError():
    comb = CombustionEngine('combustion engine')
    with raises(TESPyConnectionError):
        Connection(comb, 'out1', comb, 'in1')

##############################################################################
# ref


def create_ref_TypeError(params):
    with raises(TypeError):
        Ref(params[0], params[1], params[2])


def test_ref_creation_error():
    """Test errors creating reference objects."""
    comb = CombustionEngine('combustion engine')
    pipeline = Pipe('pipeline')
    conn = Connection(comb, 'out1', pipeline, 'in1')
    create_ref_TypeError([conn, 7, 'hi'])
    create_ref_TypeError([conn, 'hi', 0])
    create_ref_TypeError([comb, 1, 0])

##############################################################################
# bus


def bus_add_comps_TypeError(b, c):
    with raises(TypeError):
        b.add_comps(c)


def test_Bus_add_comps_errors():
    """Test errors adding components to busses."""
    mybus = Bus('mybus')
    comb = CombustionEngine('combustion engine')
    pipeline = Pipe('pipeline')
    conn = Connection(comb, 'out1', pipeline, 'in1')
    bus_add_comps_TypeError(mybus, {'comp': conn})
    bus_add_comps_TypeError(mybus, {'f': comb})
    bus_add_comps_TypeError(mybus, {'comp': comb, 'char': 'Hi'})
    bus_add_comps_TypeError(mybus, {'comp': comb, 'param': 5})
    bus_add_comps_TypeError(mybus, {'comp': comb, 'P_ref': 'what'})
    bus_add_comps_TypeError(mybus, comb)

    with raises(ValueError):
        mybus.add_comps({'comp': comb, 'base': 5})

##############################################################################
# test errors of UserDefinedEquation


def udf_dummy():
    return


def test_UserDefinedEquation_errors():
    with raises(TypeError):
        UserDefinedEquation(7, udf_dummy, udf_dummy, [])

##############################################################################
# test errors of component classes

##############################################################################
# CombustionChamber


def test_CombustionChamber_missing_fuel():
    """Test no fuel in network."""
    nw = Network()
    instance = CombustionChamber('combustion chamber')
    c1 = Connection(Source('air'), 'out1', instance, 'in1')
    c2 = Connection(Source('fuel'), 'out1', instance, 'in2')
    c3 = Connection(instance, 'out1', Sink('flue gas'), 'in1')
    nw.add_conns(c1, c2, c3)
    with raises(TESPyComponentError):
        nw.solve('design', init_only=True)


def test_CombustionChamber_missing_oxygen():
    """Test no oxygen in network."""
    nw = Network()
    instance = CombustionChamber('combustion chamber')
    c1 = Connection(Source('air'), 'out1', instance, 'in1')
    c2 = Connection(Source('fuel'), 'out1', instance, 'in2')
    c3 = Connection(instance, 'out1', Sink('flue gas'), 'in1')
    nw.add_conns(c1, c2, c3)
    with raises(TESPyComponentError):
        nw.solve('design', init_only=True)

##############################################################################
# combustion_engine


class TestCombustionEngineBusErrors:

    def setup_method(self):
        self.nw = Network()
        self.instance = CombustionEngine('combustion engine')
        self.bus = Bus('power')
        self.bus.add_comps({'comp': self.instance, 'param': 'Param'})

    def test_missing_Bus_param_func(self):
        """Test wrong/missing bus parameter in bus function."""
        with raises(ValueError):
            self.instance.bus_func(self.bus.comps.loc[self.instance])

##############################################################################
# compressor


def test_compressor_missing_char_parameter():
    """Compressor with invalid parameter for eta_s_char function."""
    nw = Network()
    so = Source('source')
    si = Sink('sink')
    instance = Compressor('compressor')
    c1 = Connection(so, 'out1', instance, 'in1')
    c2 = Connection(instance, 'out1', si, 'in1')
    nw.add_conns(c1, c2)
    c1.set_attr(fluid={"CH4": 1})
    instance.set_attr(eta_s_char={
        'func': CharLine([0, 1], [1, 2]), 'is_set': True, 'param': None})
    nw.solve('design', init_only=True)
    with raises(ValueError):
        instance.eta_s_char_func()

##############################################################################
# subsystems


def test_subsys_label_str():
    with raises(ValueError):
        Subsystem(5)


def test_subsys_label_forbidden():
    with raises(ValueError):
        Subsystem('label;')


def test_subsys_no_create_network():
    with raises(NotImplementedError):
        Subsystem("bare subsystem")

##############################################################################
# turbine


def test_Turbine_missing_char_parameter():
    """Turbine with invalid parameter for eta_s_char function."""
    nw = Network()
    so = Source('source')
    si = Sink('sink')
    instance = Turbine('turbine')
    c1 = Connection(so, 'out1', instance, 'in1')
    c2 = Connection(instance, 'out1', si, 'in1')
    nw.add_conns(c1, c2)
    c1.set_attr(fluid={"CH4": 1})
    instance.set_attr(eta_s_char={
        'char_func': CharLine([0, 1], [1, 2]), 'is_set': True, 'param': None})
    nw.solve('design', init_only=True)
    with raises(ValueError):
        instance.eta_s_char_func()


def test_wrong_Bus_param_func():
    """Test missing/wrong bus parameter specification in equations."""
    # this test does not need setup, since the function is called without
    # network initialisation
    instance = WaterElectrolyzer('electrolyzer')
    some_bus = Bus('some_bus')
    param = 'G'
    some_bus.add_comps({'comp': instance, 'param': param})
    with raises(ValueError):
        instance.bus_func(some_bus.comps.loc[instance])


##############################################################################
# test errors of Network class

class TestNetworkErrors:

    def setup_method(self):
        self.nw = Network()

    def test_add_conns_TypeError(self):
        with raises(TypeError):
            self.nw.add_conns(Component('test'))

    def test_no_connections_error(self):
        with raises(TESPyNetworkError):
            self.nw.solve('design')

    def test_bad_fluids_in_fluid_vector(self):
        source1 = Source('source1')
        sink1 = Sink('sink1')
        a = Connection(source1, 'out1', sink1, 'in1', fluid={'air': 1})
        self.nw.add_conns(a)
        with raises(TESPyNetworkError):
            self.nw.solve('design')

    def test_duplicate_Connection_labels(self):
        source1 = Source('source1')
        source2 = Source('source2')
        sink1 = Sink('sink1')
        sink2 = Sink('sink2')
        a = Connection(source1, 'out1', sink1, 'in1', label='myconn')
        b = Connection(source2, 'out1', sink2, 'in1', label='myconn')
        with raises(ValueError):
            self.nw.add_conns(a, b)

    def test_Connection_error_source(self):
        source = Source('source')
        sink1 = Sink('sink1')
        sink2 = Sink('sink2')
        a = Connection(source, 'out1', sink1, 'in1')
        b = Connection(source, 'out1', sink2, 'in1')
        self.nw.add_conns(a, b)
        with raises(TESPyNetworkError):
            self.nw.check_topology()

    def test_Connection_error_target(self):
        source1 = Source('source1')
        source2 = Source('source2')
        sink = Sink('sink')
        a = Connection(source1, 'out1', sink, 'in1')
        b = Connection(source2, 'out1', sink, 'in1')
        self.nw.add_conns(a, b)
        with raises(TESPyNetworkError):
            self.nw.check_topology()

    def test_consistency_inlets(self):
        merge = Merge('merge')
        sink = Sink('label')
        a = Connection(merge, 'out1', sink, 'in1')
        self.nw.add_conns(a)
        with raises(TESPyNetworkError):
            self.nw.check_topology()

    def test_consistency_outlets(self):
        source = Source('source')
        splitter = Splitter('splitter')
        a = Connection(source, 'out1', splitter, 'in1')
        self.nw.add_conns(a)
        with raises(TESPyNetworkError):
            self.nw.check_topology()

    def test_component_label_duplicates(self):
        source = Source('label')
        sink = Sink('label')
        a = Connection(source, 'out1', sink, 'in1')
        with raises(TESPyNetworkError):
            self.nw.add_conns(a)

    def test_missing_offdesign_path(self):
        source = Source('source')
        sink = Sink('sink')
        a = Connection(source, 'out1', sink, 'in1')
        self.nw.add_conns(a)
        with raises(TESPyNetworkError):
            self.nw.solve('offdesign')

    def test_bad_mode_specification(self):
        source = Source('source')
        sink = Sink('sink')
        a = Connection(source, 'out1', sink, 'in1')
        self.nw.add_conns(a)
        with raises(ValueError):
            self.nw.solve('ofdesign')

    def test_underdetermination(self):
        source = Source('source')
        sink = Sink('sink')
        a = Connection(source, 'out1', sink, 'in1', m=1)
        self.nw.add_conns(a)
        with raises(TESPyNetworkError):
            self.nw.solve('design')

    def test_overdetermination(self):
        source = Source('source')
        sink = Sink('sink')
        a = Connection(source, 'out1', sink, 'in1', m=1, p=1e5, x=1, h=1e6,
                       fluid={'water': 1}, fluid_balance=True)
        self.nw.add_conns(a)
        with raises(TESPyNetworkError):
            self.nw.solve('design')

    def test_add_Bus_TypeError(self):
        source = Source('label')
        sink = Sink('label')
        a = Connection(source, 'out1', sink, 'in1')
        with raises(TypeError):
            self.nw.add_busses(a)

    def test_Bus_duplicate(self):
        with raises(TESPyNetworkError):
            b = Bus('mybus')
            self.nw.add_busses(b, b)

    def test_buslabel_duplicate(self):
        with raises(TESPyNetworkError):
            a = Bus('mybus')
            b = Bus('mybus')
            self.nw.add_busses(a, b)

##############################################################################
# test errors of characteristics classes


def test_char_number_of_points():
    with raises(ValueError):
        CharLine(x=[0, 1, 2], y=[1, 2, 3, 4])


def test_CharMap_number_of_points():
    with raises(ValueError):
        CharMap(x=[0, 1, 2], y=[[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3]])


def test_CharMap_number_of_dimensions():
    with raises(ValueError):
        CharMap(x=[0, 1, 2], y=[[1, 2, 3, 4], [1, 2, 3, 4]])


def test_CharMap_y_z_dimension_mismatch():
    with raises(ValueError):
        CharMap(x=[0, 1], y=[[1, 2, 3, 4], [1, 2, 3, 4]],
                z=[[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4]])


def test_missing_CharLine_files():
    """Test missing files."""
    path = extend_basic_path('data')
    tmp_path = extend_basic_path('tmp_dir_for_testing')

    if os.path.exists(path):
        for f in os.listdir(path):
            shutil.copy(src=path + '/' + f, dst=tmp_path)

        shutil.rmtree(path, ignore_errors=True)

    with raises(FileNotFoundError):
        load_custom_char('stuff', CharLine)

    if os.path.exists(tmp_path):
        for f in os.listdir(tmp_path):
            shutil.copy(src=tmp_path + '/' + f, dst=path)

        shutil.rmtree(tmp_path, ignore_errors=True)


def test_missing_CharMap_files():
    """Test missing files."""
    path = extend_basic_path('data')
    tmp_path = extend_basic_path('tmp_dir_for_testing')

    if os.path.exists(path):
        for f in os.listdir(path):
            shutil.copy(src=path + '/' + f, dst=tmp_path)

        shutil.rmtree(path, ignore_errors=True)

    with raises(FileNotFoundError):
        load_custom_char('some other stuff', CharMap)

    if os.path.exists(tmp_path):
        for f in os.listdir(tmp_path):
            shutil.copy(src=tmp_path + '/' + f, dst=path)

        shutil.rmtree(tmp_path, ignore_errors=True)


##############################################################################
# test errors in fluid porperty functions


def test_h_mix_pQ_on_mixtures():
    c = Connection(Source("test"), "out1", Sink("test2"), "in1")
    c.set_attr(fluid={"O2": 0.24, "N2": 0.76})
    c._create_fluid_wrapper()
    with raises(ValueError):
        h_mix_pQ(1e5, 0.5, c.fluid_data, c.mixing_rule)


def test_warning_logged_with_correct_category(caplog):
    logger_name = "TESPyLogger"
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.WARNING)

    with caplog.at_level(logging.WARNING, logger=logger_name):
        warnings.showwarning = FutureWarningHandler(logger)
        warnings.warn("Custom test message", category=UserWarning)

        assert any([
            "UserWarning: Custom test message" in msg
            for msg in caplog.messages
        ])
