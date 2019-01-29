# -*- coding: utf-8

from nose.tools import ok_, eq_, raises, with_setup

from tespy import nwk, cmp, con, hlp
from CoolProp.CoolProp import PropsSI as CP
import numpy as np


# %% bulk tests


class specification_error_tests:

    def setup(self):
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

        self.create_connection_ValueError('source')
        self.create_connection_ValueError('target')

        # TypeErrors
        self.set_attr_TypeError(self.comp, P=[5])
        self.set_attr_TypeError(self.comp, tiP_char=None)

        self.set_attr_TypeError(self.conn, design='h')
        self.set_attr_TypeError(self.conn, fluid_balance=1)
        self.set_attr_TypeError(self.conn, h0=[4])
        self.set_attr_TypeError(self.conn, fluid=5)

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

# %% Single tests

@raises(ValueError)
def test_interface_ValueError():
    # interface specification
    cmp.sink('sink', interface=5)