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
    def test_set_attr_ValueError(self):
        self.comp.set_attr(mode=5)
        self.comp.set_attr(offdesign=['Q'])
        con.connection(self.comp, 'out1', self.pipe, 'in5')
        con.connection(self.comp, 'out6', self.pipe, 'in1')
        self.conn.set_attr(offdesign=['f'])

    @raises(TypeError)
    def test_set_attr_TypeError(self):
        self.comp.set_attr(P=[5], tiP_char=None)
        con.connection(self.comp, 'out1', self.conn, 'in5')
        self.conn.set_attr(design='h', fluid_balance=1, fluid=5, h0=[4])
        self.bus.add_comps({'c': self.conn})
        self.bus.add_comps({'f': self.comp})
        self.bus.add_comps({'c': self.comp, 'char': 'Hi'})
        self.bus.add_comps({'c': self.comp, 'p': 5})
        self.bus.add_comps({'c': self.comp, 'P_ref': 'what'})
        con.connection(self.comp, 'out1', self.pipe, 'in1', m=con.ref(self.conn, 7, 'hi'))
        con.connection(self.comp, 'out1', self.pipe, 'in1', m=con.ref(self.conn, 'Hi', 0))
        con.connection(self.comp, 'out1', self.pipe, 'in1', m=con.ref(self.comp, 1, 0))

    @raises(KeyError)
    def test_set_attr_KeyError(self):
        self.comp.set_attr(wow=5)
        self.conn.set_attr(key=7)
        self.comp.get_attr('wow')
        self.conn.get_attr('key')
        self.bus.get_attr('components')
        con.ref(self.conn, 1, 0).get_attr('comp')

    @raises(hlp.TESPyComponentError)
    def test_set_attr_TESPyError(self):
        self.comp.set_attr(interface=True)

    @raises(hlp.TESPyConnectionError)
    def test_con_instanciation(self):
        con.connection(self.comp, 'out1', self.comp, 'in1')
        self.bus.add_comps(self.comp)

