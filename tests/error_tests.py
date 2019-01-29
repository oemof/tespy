# -*- coding: utf-8

from nose.tools import ok_, eq_, raises, with_setup

from tespy import nwk, cmp, con, hlp
from CoolProp.CoolProp import PropsSI as CP
import numpy as np


class specification_error_tests:

    def setup(self):
        self.comp = cmp.cogeneration_unit('cogeneration unit')
        self.pipe = cmp.pipe('pipe')
        self.conn = con.connection(self.comp, 'out1', self.pipe, 'in1')
        self.bus = con.bus('mybus')

    @raises(ValueError)
    def test_cmp_instanciation(self):
        labels = [5, 'Label,', 'Labe;l', 'Label.']
        for l in labels:
            cmp.cogeneration_unit(l)

        # interface specification
        cmp.sink('sink', interface=5)

    @raises(ValueError)
    def set_attr_ValueError(self, instance, **kwargs):
        instance.set_attr(**kwargs)

    @raises(TypeError)
    def set_attr_TypeError(self, instance, **kwargs):
        instance.set_attr(**kwargs)

    @raises(KeyError)
    def set_attr_KeyError(self, instance, **kwargs):
        instance.set_attr(**kwargs)

    @raises(KeyError)
    def set_attr_KeyError(self, instance, key):
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
        # ValueErrors
        self.set_attr_ValueError(self.comp, mode=5)
        self.set_attr_ValueError(self.comp, offdesign=['Q'])
        self.set_attr_ValueError(self.conn, offdesign=['f'])


#        con.connection(self.comp, 'out1', self.pipe, 'in5')
#        con.connection(self.comp, 'out6', self.pipe, 'in1')

        # TypeErrors
        self.set_attr_TypeError(self.comp, P=[5])
        self.set_attr_TypeError(self.comp, tiP_char=None)

        self.set_attr_TypeError(self.conn, design='h')
        self.set_attr_TypeError(self.conn, fluid_balance=1)
        self.set_attr_TypeError(self.conn, h0=[4])
        self.set_attr_TypeError(self.conn, fluid=5)

#        con.connection(self.comp, 'out1', self.conn, 'in5')
#        self.bus.add_comps({'c': self.conn})
#        self.bus.add_comps({'f': self.comp})
#        self.bus.add_comps({'c': self.comp, 'char': 'Hi'})
#        self.bus.add_comps({'c': self.comp, 'p': 5})
#        self.bus.add_comps({'c': self.comp, 'P_ref': 'what'})
#        con.connection(self.comp, 'out1', self.pipe, 'in1', m=con.ref(self.conn, 7, 'hi'))
#        con.connection(self.comp, 'out1', self.pipe, 'in1', m=con.ref(self.conn, 'Hi', 0))
#        con.connection(self.comp, 'out1', self.pipe, 'in1', m=con.ref(self.comp, 1, 0))

        # KeyErrors
        self.set_attr_KeyError(self.comp, wow=5)
        self.set_attr_KeyError(self.conn, jey=5)

        self.get_attr_KeyError(self.comp, 'wow')
        self.get_attr_KeyError(self.conn, 'key')
        self.get_attr_KeyError(self.bus, 'components')
        self.get_attr_KeyError(con.ref(self.conn, 1, 0), 'comp')
