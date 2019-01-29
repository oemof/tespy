# -*- coding: utf-8

from nose.tools import ok_, eq_, raises, with_setup

from tespy import nwk, cmp, con, hlp
from CoolProp.CoolProp import PropsSI as CP
import numpy as np


# %% bulk tests


class specification_error_tests:

    def setup(self):
        self.nw = nwk.network(['water', 'air'])
        self.comp = cmp.cogeneration_unit('cogeneration unit')
        self.pipe = cmp.pipe('pipe')
        self.conn = con.connection(self.comp, 'out1', self.pipe, 'in1')
        self.bus = con.bus('mybus')

    @raises(ValueError)
    def cmp_instanciation_ValueError(self, label):
        cmp.cogeneration_unit(label)

    @raises(ValueError)
    def set_attr_ValueError(self, instance, **kwargs):
        instance.set_attr(**kwargs)

    @raises(ValueError)
    def create_connection_ValueError(self, pos):
        if pos == 'source':
            con.connection(self.comp, 'out6', self.pipe, 'in1')
        elif pos == 'target':
            con.connection(self.comp, 'out1', self.pipe, 'in5')

    @raises(TypeError)
    def set_attr_TypeError(self, instance, **kwargs):
        instance.set_attr(**kwargs)

    @raises(TypeError)
    def bus_add_comps_TypeError(self, c):
        self.bus.add_comps(c)

    @raises(TypeError)
    def test_create_conn_TypeError(self):
        con.connection(self.comp, 'out1', self.conn, 'in1')

    @raises(TypeError)
    def create_ref_TypeError(self, params):
        con.ref(params[0], params[1], params[2])

    @raises(KeyError)
    def set_attr_KeyError(self, instance, **kwargs):
        instance.set_attr(**kwargs)

    @raises(KeyError)
    def get_attr_KeyError(self, instance, key):
        instance.get_attr(key)

    @raises(hlp.TESPyComponentError)
    def test_TESPyComponentError(self):
        self.comp.set_attr(interface=True)

    @raises(hlp.TESPyConnectionError)
    def test_con_TESPyConnectionError(self):
        con.connection(self.comp, 'out1', self.comp, 'in1')

    @raises(hlp.TESPyConnectionError)
    def test_bus_TESPyConnectionError(self):
        self.bus.add_comps(self.comp)

    def test_set_attr_errors(self):
        #
        labels = [5, 'Label,', 'Labe;l', 'Label.']
        for l in labels:
            self.cmp_instanciation_ValueError(l)

        # ValueErrors
        self.set_attr_ValueError(self.comp, mode=5)
        self.set_attr_ValueError(self.comp, offdesign=['Q'])

        self.set_attr_ValueError(self.conn, offdesign=['f'])

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

        self.set_attr_TypeError(self.conn, design='h')
        self.set_attr_TypeError(self.conn, fluid_balance=1)
        self.set_attr_TypeError(self.conn, h0=[4])
        self.set_attr_TypeError(self.conn, fluid=5)

        self.set_attr_TypeError(self.nw, p_range=5)
        self.set_attr_TypeError(self.nw, h_range=5)
        self.set_attr_TypeError(self.nw, T_range=5)

        self.bus_add_comps_TypeError({'c': self.conn})
        self.bus_add_comps_TypeError({'f': self.comp})
        self.bus_add_comps_TypeError({'c': self.comp, 'char': 'Hi'})
        self.bus_add_comps_TypeError({'c': self.comp, 'p': 5})
        self.bus_add_comps_TypeError({'c': self.comp, 'P_ref': 'what'})

        self.create_ref_TypeError([self.conn, 7, 'hi'])
        self.create_ref_TypeError([self.conn, 'hi', 0])
        self.create_ref_TypeError([self.comp, 1, 0])

        # KeyErrors
        self.set_attr_KeyError(self.comp, wow=5)
        self.set_attr_KeyError(self.conn, jey=5)

        self.get_attr_KeyError(self.comp, 'wow')
        self.get_attr_KeyError(self.conn, 'key')
        self.get_attr_KeyError(self.bus, 'components')
        self.get_attr_KeyError(con.ref(self.conn, 1, 0), 'comp')
        self.get_attr_KeyError(self.nw, 'test')

# %% Single tests

@raises(ValueError)
def test_interface_ValueError():
    # interface specification
    cmp.sink('sink', interface=5)

@raises(ValueError)
def test_network_print_level():
    nwk.network(['INCOMP::DowQ']).set_printoptions(print_level='error')

@raises(TypeError)
def test_network_instanciation_single_fluid():
    nwk.network('water')

@raises(TypeError)
def test_network_instanciation_no_fluids():
    nwk.network(['water']).add_conns(cmp.component('test'))

@raises(hlp.TESPyNetworkError)
def test_network_connection_error_source():
    nw = nwk.network(['water'])
    source = cmp.source('source')
    sink1 = cmp.sink('sink1')
    sink2 = cmp.sink('sink2')
    a = con.connection(source, 'out1', sink1, 'in1')
    b = con.connection(source, 'out1', sink2, 'in1')
    nw.add_conns(a, b)
    nw.check_network()

@raises(hlp.TESPyNetworkError)
def test_network_connection_error_target():
    nw = nwk.network(['water'])
    source1 = cmp.source('source1')
    source2 = cmp.source('source2')
    sink = cmp.sink('sink')
    a = con.connection(source1, 'out1', sink, 'in1')
    b = con.connection(source2, 'out1', sink, 'in1')
    nw.add_conns(a, b)
    nw.check_network()

@raises(hlp.TESPyNetworkError)
def test_network_component_labels():
    nw = nwk.network(['water'])
    source = cmp.source('label')
    sink = cmp.sink('label')
    a = con.connection(source, 'out1', sink, 'in1')
    nw.add_conns(a)
    nw.check_network()

@raises(hlp.TESPyNetworkError)
def test_network_offdesign_path():
    nw = nwk.network(['water'])
    source = cmp.source('label')
    sink = cmp.sink('label')
    a = con.connection(source, 'out1', sink, 'in1')
    nw.add_conns(a)
    nw.solve('offdesign')

@raises(hlp.TESPyNetworkError)
def test_network_underdetermination():
    nw = nwk.network(['water'])
    source = cmp.source('label')
    sink = cmp.sink('label')
    a = con.connection(source, 'out1', sink, 'in1', m=1)
    nw.add_conns(a)
    nw.solve('design')

@raises(hlp.TESPyNetworkError)
def test_network_overdetermination():
    nw = nwk.network(['water'])
    source = cmp.source('label')
    sink = cmp.sink('label')
    a = con.connection(source, 'out1', sink, 'in1', m=1, p=1e5, x=1, h=1e6, fluid={'water': 1}, fluid_balance=True)
    nw.add_conns(a)
    nw.solve('design')

@raises(ValueError)
def test_network_mode():
    nw = nwk.network(['water'])
    source = cmp.source('label')
    sink = cmp.sink('label')
    a = con.connection(source, 'out1', sink, 'in1')
    nw.add_conns(a)
    nw.solve('ofdesign')

@raises(hlp.TESPyNetworkError)
def test_network_instanciation_no_fluids():
    nwk.network([])