# -*- coding: utf-8

from nose.tools import ok_, eq_, raises, with_setup

from tespy import nwk, cmp, con, hlp
from CoolProp.CoolProp import PropsSI as CP
import numpy as np


class specification_error_tests:

    def setup(self):
        self.comp = cmp.cogeneration_unit('cogeneration unit')

    @raises(ValueError)
    def test_cmp_instanciation(self):
        labels = [5, 'Label,', 'Labe;l', 'Label.']
        for l in labels:
            cmp.cogeneration_unit(l)

        # interface specification
        cmp.sink('sink', interface=5)

    @raises(ValueError)
    def cmp_set_attr_ValueError(self, **kwargs):
        self.comp.set_attr(**kwargs)

    @raises(TypeError)
    def cmp_set_attr_TypeError(self, **kwargs):
        self.comp.set_attr(**kwargs)

    @raises(KeyError)
    def cmp_set_attr_KeyError(self, **kwargs):
        self.comp.set_attr(**kwargs)

    @raises(hlp.TESPyComponentError)
    def cmp_set_attr_TESPyError(self, **kwargs):
        self.comp.set_attr(**kwargs)

    @with_setup(setup)
    def test_cmp_parameterisation(self):
        # value errors
        self.cmp_set_attr_ValueError(mode=5)
        self.cmp_set_attr_ValueError(design=['P', 'Q1'], offdesign=['Q'])
        # key errors
        self.cmp_set_attr_KeyError(wow=5)
        # type errors
        self.cmp_set_attr_TypeError(P=[5])
        self.cmp_set_attr_TypeError(tiP_char=None)
        # TESPy errors
        self.cmp_set_attr_TESPyError(interface=True)

        try:
            self.comp.set_attr(interface=True)
        except hlp.TESPyComponentError:
            pass
