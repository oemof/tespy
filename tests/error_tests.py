# -*- coding: utf-8

from nose.tools import ok_, eq_, raises, with_setup

from tespy import nwk, cmp, con, hlp
from CoolProp.CoolProp import PropsSI as CP
import numpy as np


class specification_error_tests:

    def setup(self):
        self.comp = cmp.cogeneration_unit('cogeneration unit')

    def test_component_instanciation(self):
        labels = [5, 'Label,', 'Labe;l', 'Label.']
        for l in labels:
            try:
                cmp.cogeneration_unit(l)
            except ValueError:
                continue
        pass

    @with_setup(setup)
    def test_component_parameterisation(self):
        # component properties
        try:
            self.comp.set_attr(P=[5])
        except TypeError:
            pass

        # component characteristics
        try:
            self.comp.set_attr(tiP_char=None)
        except TypeError:
            pass

        # non existent property
        try:
            self.comp.set_attr(wow=5)
        except ValueError:
            pass

        # mode specification
        try:
            self.comp.set_attr(mode=5)
        except ValueError:
            pass

        # design/offdesign parameter specification
        try:
            self.comp.set_attr(design=['P', 'Q1'], offdesign=['Q'])
        except ValueError:
            pass

        # interface specification
        try:
            self.comp.set_attr(interface=True)
        except hlp.TESPyComponentError:
            pass

        # interface specification
        try:
            cmp.sink('sink', interface=5)
        except ValueError:
            pass
