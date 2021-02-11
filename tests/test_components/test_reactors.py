# -*- coding: utf-8

"""Module for testing components of type reactor.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_reactors.py

SPDX-License-Identifier: MIT
"""

import shutil

import numpy as np

from tespy.components import Sink
from tespy.components import Source
from tespy.components import WaterElectrolyzer
from tespy.connections import Bus
from tespy.connections import Connection
from tespy.networks import Network


def convergence_check(lin_dep):
    """Check convergence status of a simulation."""
    msg = 'Calculation did not converge!'
    assert lin_dep is False, msg


class TestReactors:

    def setup(self):
        """Set up network for electrolyzer tests."""
        self.nw = Network(['O2', 'H2', 'H2O'], T_unit='C', p_unit='bar')
        self.instance = WaterElectrolyzer('electrolyzer')

        fw = Source('feed water')
        cw_in = Source('cooling water')
        o2 = Sink('oxygen sink')
        h2 = Sink('hydrogen sink')
        cw_out = Sink('cooling water sink')

        self.instance.set_attr(pr=0.99, eta=1)

        cw_el = Connection(cw_in, 'out1', self.instance, 'in1',
                           fluid={'H2O': 1, 'H2': 0, 'O2': 0}, T=20, p=1)
        el_cw = Connection(self.instance, 'out1', cw_out, 'in1', T=45)

        self.nw.add_conns(cw_el, el_cw)

        fw_el = Connection(fw, 'out1', self.instance, 'in2', label='h2o')
        el_o2 = Connection(self.instance, 'out2', o2, 'in1')
        el_h2 = Connection(self.instance, 'out3', h2, 'in1', label='h2')

        self.nw.add_conns(fw_el, el_o2, el_h2)

    def test_WaterElectrolyzer(self):
        """Test component properties of water electrolyzer."""
        # check bus function:
        # power output on component and bus must be indentical
        self.nw.get_conn('h2o').set_attr(T=25, p=1)
        self.nw.get_conn('h2').set_attr(T=25)
        power = Bus('power')
        power.add_comps({'comp': self.instance, 'param': 'P', 'base': 'bus'})
        power.set_attr(P=2.5e6)
        self.nw.add_busses(power)

        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        msg = ('Value of power must be ' + str(power.P.val) + ', is ' +
               str(self.instance.P.val) + '.')
        assert round(power.P.val, 1) == round(self.instance.P.val), msg

        # effieciency was set to 100 % with inlet and outlet states of the
        # reaction educts and products beeing identical to reference state
        # therefore Q must be equal to 0
        msg = ('Value of heat output must be 0.0, is ' +
               str(self.instance.Q.val) + '.')
        assert round(self.instance.Q.val, 4) == 0.0, msg

        # reset power, change efficiency value and specify heat bus value
        power.set_attr(P=np.nan)
        self.nw.get_conn('h2o').set_attr(T=25, p=1)
        self.nw.get_conn('h2').set_attr(T=50)
        self.instance.set_attr(eta=0.8)
        # check bus function:
        # heat output on component and bus must be indentical
        heat = Bus('heat')
        heat.add_comps({'comp': self.instance, 'param': 'Q'})
        heat.set_attr(P=-8e5)
        self.nw.add_busses(heat)

        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        msg = ('Value of heat flow must be ' + str(heat.P.val) +
               ', is ' + str(self.instance.Q.val) + '.')
        assert round(heat.P.val, 1) == round(self.instance.Q.val), msg
        self.nw.save('tmp')

        # check bus function:
        # heat output on component and bus must identical (offdesign test)
        Q = heat.P.val * 0.9
        heat.set_attr(P=Q)
        self.nw.solve('offdesign', design_path='tmp')
        convergence_check(self.nw.lin_dep)
        msg = ('Value of heat flow must be ' + str(Q) +
               ', is ' + str(self.instance.Q.val) + '.')
        assert round(Q, 1) == round(self.instance.Q.val), msg

        # delete both busses again
        self.nw.del_busses(heat, power)

        # test efficiency vs. specific energy consumption
        self.nw.get_conn('h2').set_attr(m=0.1)
        self.instance.set_attr(eta=0.9, e='var')
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        msg = ('Value of efficiency must be ' + str(self.instance.eta.val) +
               ', is ' + str(self.instance.e0 / self.instance.e.val) + '.')
        eta = round(self.instance.eta.val, 2)
        eta_calc = round(self.instance.e0 / self.instance.e.val, 2)
        assert eta == eta_calc, msg

        # test efficiency value > 1, Q must be larger than 0
        e = 130e6
        self.instance.set_attr(e=np.nan, eta=np.nan)
        self.instance.set_attr(e=e)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        # test efficiency
        msg = ('Value of efficiency must be ' + str(self.instance.e0 / e) +
               ', is ' + str(self.instance.eta.val) + '.')
        eta = round(self.instance.e0 / e, 2)
        eta_calc = round(self.instance.eta.val, 2)
        assert eta == eta_calc, msg
        # test Q
        msg = ('Value of heat must be larger than zero, is ' +
               str(self.instance.Q.val) + '.')
        assert self.instance.Q.val > 0, msg

        # test specific energy consumption
        e = 150e6
        self.instance.set_attr(e=np.nan, eta=np.nan)
        self.instance.set_attr(e=e)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        msg = ('Value of specific energy consumption e must be ' + str(e) +
               ', is ' + str(self.instance.e.val) + '.')
        assert round(e, 1) == round(self.instance.e.val, 1), msg

        # test cooling loop pressure ratio, zeta as variable value
        pr = 0.95
        self.instance.set_attr(
            pr=pr, e=None, eta=None, zeta='var', P=2e7, design=['pr'])
        self.nw.solve('design')
        shutil.rmtree('./tmp', ignore_errors=True)
        self.nw.save('tmp')
        convergence_check(self.nw.lin_dep)
        msg = ('Value of pressure ratio must be ' + str(pr) + ', is ' +
               str(self.instance.pr.val) + '.')
        assert round(pr, 2) == round(self.instance.pr.val, 2), msg

        # use zeta as offdesign parameter, at design point pressure
        # ratio must not change
        self.instance.set_attr(zeta=np.nan, offdesign=['zeta'])
        self.nw.solve('offdesign', design_path='tmp')
        convergence_check(self.nw.lin_dep)
        msg = ('Value of pressure ratio must be ' + str(pr) + ', is ' +
               str(self.instance.pr.val) + '.')
        assert round(pr, 2) == round(self.instance.pr.val, 2), msg

        # test heat output specification in offdesign mode
        Q = self.instance.Q.val * 0.9
        self.instance.set_attr(Q=Q, P=np.nan)
        self.nw.solve('offdesign', design_path='tmp')
        convergence_check(self.nw.lin_dep)
        msg = ('Value of heat must be ' + str(Q) + ', is ' +
               str(self.instance.Q.val) + '.')
        assert round(Q, 0) == round(self.instance.Q.val, 0), msg
        shutil.rmtree('./tmp', ignore_errors=True)
