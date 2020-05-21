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

from tespy.components.basics import sink
from tespy.components.basics import source
from tespy.components.reactors import water_electrolyzer
from tespy.connections import bus
from tespy.connections import connection
from tespy.networks.networks import network


def convergence_check(lin_dep):
    """Check convergence status of a simulation."""
    msg = 'Calculation did not converge!'
    assert lin_dep is False, msg


class TestReactors:

    def setup(self):
        """Set up network for electrolyzer tests."""
        self.nw = network(['O2', 'H2', 'H2O'], T_unit='C', p_unit='bar')
        self.instance = water_electrolyzer('electrolyzer')

        fw = source('feed water')
        cw_in = source('cooling water')
        o2 = sink('oxygen sink')
        h2 = sink('hydrogen sink')
        cw_out = sink('cooling water sink')

        self.instance.set_attr(pr_c=0.99)

        cw_el = connection(cw_in, 'out1', self.instance, 'in1',
                           fluid={'H2O': 1, 'H2': 0, 'O2': 0}, T=20, p=1)
        el_cw = connection(self.instance, 'out1', cw_out, 'in1', T=45)

        self.nw.add_conns(cw_el, el_cw)

        fw_el = connection(fw, 'out1', self.instance, 'in2', m=0.1, T=20, p=10)
        el_o2 = connection(self.instance, 'out2', o2, 'in1')
        el_h2 = connection(self.instance, 'out3', h2, 'in1', T=50)

        self.nw.add_conns(fw_el, el_o2, el_h2)

    def test_water_electrolyzer(self):
        """Test component properties of water electrolyzer."""
        # check bus function:
        # power output on component and bus must be indentical
        power = bus('power')
        power.add_comps({'comp': self.instance, 'param': 'P'})
        power.set_attr(P=2.5e6)
        self.nw.add_busses(power)

        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        msg = ('Value of power must be ' + str(power.P.val) + ', is ' +
               str(self.instance.P.val) + '.')
        assert round(power.P.val, 1) == round(self.instance.P.val), msg
        power.set_attr(P=np.nan)

        # check bus function:
        # heat output on component and bus must be indentical
        heat = bus('heat')
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
        self.instance.set_attr(eta=0.9, e='var')
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        msg = ('Value of efficiency must be ' + str(self.instance.eta.val) +
               ', is ' + str(self.instance.e0 / self.instance.e.val) + '.')
        eta = round(self.instance.eta.val, 2)
        eta_calc = round(self.instance.e0 / self.instance.e.val, 2)
        assert eta == eta_calc, msg

        # test efficiency value > 1
        e = 130e6
        self.instance.set_attr(e=np.nan, eta=np.nan)
        self.instance.set_attr(e=e)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        msg = ('Value of efficiency must be ' + str(self.instance.e0 / e) +
               ', is ' + str(self.instance.eta.val) + '.')
        eta = round(self.instance.e0 / e, 2)
        eta_calc = round(self.instance.eta.val, 2)
        assert eta == eta_calc, msg

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
        self.instance.set_attr(pr_c=pr, e=np.nan, zeta='var',
                               P=2.5e6, design=['pr_c'])
        self.nw.solve('design')
        shutil.rmtree('./tmp', ignore_errors=True)
        self.nw.save('tmp')
        convergence_check(self.nw.lin_dep)
        msg = ('Value of pressure ratio must be ' + str(pr) + ', is ' +
               str(self.instance.pr_c.val) + '.')
        assert round(pr, 2) == round(self.instance.pr_c.val, 2), msg

        # use zeta as offdesign parameter, at design point pressure
        # ratio must not change
        self.instance.set_attr(zeta=np.nan, offdesign=['zeta'])
        self.nw.solve('offdesign', design_path='tmp')
        convergence_check(self.nw.lin_dep)
        msg = ('Value of pressure ratio must be ' + str(pr) + ', is ' +
               str(self.instance.pr_c.val) + '.')
        assert round(pr, 2) == round(self.instance.pr_c.val, 2), msg

        # test heat output specification in offdesign mode
        Q = self.instance.Q.val * 0.9
        self.instance.set_attr(Q=Q, P=np.nan)
        self.nw.solve('offdesign', design_path='tmp')
        convergence_check(self.nw.lin_dep)
        msg = ('Value of heat must be ' + str(Q) + ', is ' +
               str(self.instance.Q.val) + '.')
        assert round(Q, 0) == round(self.instance.Q.val, 0), msg
        shutil.rmtree('./tmp', ignore_errors=True)
