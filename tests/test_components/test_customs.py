# -*- coding: utf-8

"""Module for testing components of type orc evaporator.
This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_customs.py
SPDX-License-Identifier: MIT
"""
import shutil

import numpy as np

from tespy.components import ORCEvaporator
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Bus
from tespy.connections import Connection
from tespy.networks import Network
from tespy.tools.fluid_properties import T_bp_p


def convergence_check(lin_dep):
    """Check convergence status of a simulation."""
    msg = 'Calculation did not converge!'
    assert lin_dep is False, msg


class TestOrcEvaporator:

    def setup(self):
        self.nw = Network(['water', 'Isopentane'], T_unit='C', p_unit='bar',
                          h_unit='kJ / kg')
        self.inl1 = Source('inlet 1')
        self.outl1 = Sink('outlet 1')

        self.inl2 = Source('inlet 2')
        self.outl2 = Sink('outlet 2')

        self.inl3 = Source('inlet 3')
        self.outl3 = Sink('outlet 3')

        self.instance = ORCEvaporator('orc evaporator')

        self.c1 = Connection(self.inl1, 'out1', self.instance, 'in1')
        self.c2 = Connection(self.instance, 'out1', self.outl1, 'in1')
        self.c3 = Connection(self.inl2, 'out1', self.instance, 'in2')
        self.c4 = Connection(self.instance, 'out2', self.outl2, 'in1')
        self.c5 = Connection(self.inl3, 'out1', self.instance, 'in3')
        self.c6 = Connection(self.instance, 'out3', self.outl3, 'in1')

        self.nw.add_conns(self.c1, self.c2, self.c3, self.c4, self.c5, self.c6)

    def test_ORCEvaporator(self):
        """Test component properties of orc evaporator."""
        # design specification
        self.instance.set_attr(pr1=0.95, pr2=0.975, pr3=0.975,
                               design=['pr1', 'pr2', 'pr3'],
                               offdesign=['zeta1', 'zeta2', 'zeta3'])
        self.c1.set_attr(T=146.6, p=4.34, m=20.4, state='g',
                         fluid={'water': 1, 'Isopentane': 0})
        self.c3.set_attr(T=146.6, p=10.2,
                         fluid={'water': 1, 'Isopentane': 0})
        self.c4.set_attr(T=118.6)
        self.c5.set_attr(T=111.6, p=10.8,
                         fluid={'water': 0, 'Isopentane': 1})

        # test heat transfer
        Q = -6.64e+07
        self.instance.set_attr(Q=Q)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        Q_is = -self.c5.m.val_SI * (self.c6.h.val_SI - self.c5.h.val_SI)
        msg = ('Value of heat flow must be ' + str(round(Q, 0)) +
               ', is ' + str(round(Q_is, 0)) + '.')
        assert round(Q, 0) == round(Q_is, 0), msg

        # test bus
        self.instance.set_attr(Q=np.nan)
        P = -6.64e+07
        b = Bus('heat transfer', P=P)
        b.add_comps({'comp': self.instance})
        self.nw.add_busses(b)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        self.nw.save('tmp')

        Q_is = -self.c5.m.val_SI * (self.c6.h.val_SI - self.c5.h.val_SI)
        msg = ('Value of heat flow must be ' + str(round(P, 0)) +
               ', is ' + str(round(Q_is, 0)) + '.')
        assert round(P, 0) == round(Q_is, 0), msg

        # Check the state of the steam and working fluid outlet:
        x_outl1_calc = self.c2.x.val
        x_outl3_calc = self.c6.x.val
        zeta1 = self.instance.zeta1.val
        zeta2 = self.instance.zeta2.val
        zeta3 = self.instance.zeta3.val

        msg = ('Vapor mass fraction of steam outlet must be 0.0, is ' +
               str(round(x_outl1_calc, 1)) + '.')
        assert round(x_outl1_calc, 1) == 0.0, msg

        msg = ('Vapor mass fraction of working fluid outlet must be 1.0, is ' +
               str(round(x_outl3_calc, 1)) + '.')
        assert round(x_outl3_calc, 1) == 1.0, msg

        # Check offdesign by zeta values
        # geometry independent friction coefficient
        self.nw.solve('offdesign', design_path='tmp')
        convergence_check(self.nw.lin_dep)

        msg = ('Geometry independent friction coefficient '
               'at hot side 1 (steam) '
               'must be ' + str(round(zeta1, 1)) + ', is ' +
               str(round(self.instance.zeta1.val, 1)) + '.')
        assert round(self.instance.zeta1.val, 1) == round(zeta1, 1), msg
        msg = ('Geometry independent friction coefficient at '
               'hot side 2 (brine) '
               'must be ' + str(round(zeta2, 1)) + ', is ' +
               str(round(self.instance.zeta2.val, 1)) + '.')
        assert round(self.instance.zeta2.val, 1) == round(zeta2, 1), msg
        msg = ('Geometry independent friction coefficient at cold side '
               '(Isopentane) must be ' + str(round(zeta3, 1)) + ', is ' +
               str(round(self.instance.zeta3.val, 1)) + '.')
        assert round(self.instance.zeta3.val, 1) == round(zeta3, 1), msg

        # test parameters of 'subcooling' and 'overheating'
        self.instance.set_attr(subcooling=True, overheating=True)
        dT = 0.5
        self.c2.set_attr(Td_bp=-dT)
        self.c6.set_attr(Td_bp=dT)
        self.nw.solve('offdesign', design_path='tmp')
        convergence_check(self.nw.lin_dep)

        T_steam = T_bp_p(self.c2.get_flow()) - dT
        T_isop = T_bp_p(self.c6.get_flow()) + dT

        msg = ('Temperature of working fluid outlet must be ' +
               str(round(T_isop, 1)) + ', is ' +
               str(round(self.c6.T.val_SI, 1)) + '.')
        assert round(T_isop, 1) == round(self.c6.T.val_SI, 1), msg

        msg = ('Temperature of steam outlet must be ' +
               str(round(T_steam, 1)) + ', is ' +
               str(round(self.c2.T.val_SI, 1)) + '.')
        assert round(T_steam, 1) == round(self.c2.T.val_SI, 1), msg

        shutil.rmtree('./tmp', ignore_errors=True)
