import sys
import math
import logging
import unittest
sys.path.append('.\\.\\')
from tespy import nwk, con, cmp, subsys
from tespy.tools import logger
import shutil

logger.define_logging(log_path=True, log_version=True,
                      screen_level=logging.INFO, file_level=logging.DEBUG)


class dh_subsystem_tests(unittest.TestCase):
    def setUp(self):
        print("\n\n\nRunning Test...\n\n\n")
        self.nw = nwk.network(fluids=['water'], T_unit='C', p_unit='bar', h_unit='J / kg')
        self.so = cmp.source('source')
        self.si = cmp.sink('sink')
        self.heat_losses = con.bus('network losses')
        self.heat_consumer = con.bus('network consumer')

    def tearDown(self):
        shutil.rmtree('./tmp', ignore_errors=True)

    def set_busses(self):
        self.nw.check_network()
        for comp in self.nw.comps.index:
            if isinstance(comp, cmp.district_heating_pipe):
                self.heat_losses.add_comps({'c': comp})
            if (isinstance(comp, cmp.heat_exchanger_simple) and
                    not isinstance(comp, cmp.pipe)):
                self.heat_consumer.add_comps({'c': comp})
        self.nw.add_busses(self.heat_losses, self.heat_consumer)

    def test_1_consumer(self):
        self.setUp()

        # Add components
        heat_source = cmp.heat_exchanger_simple('heat_source', Q=4e4, pr=1)
        block1 = subsys.dh_consumer_block('consumer block 1', 1)
        block1.set_attr(Tamb=0,
                        DN_main=40,
                        ks_main=1e-5,
                        lambda_ins=0.03,
                        lambda_soil=1.2,
                        depth=.6,
                        dist=0.2)
        block1.set_attr(pr0=.99,
                        Q0=-3e4,
                        DN_type0=20,
                        L0=20,
                        ks0=1e-5,
                        T0=40,
                        root_dist0=30)
        self.nw.add_subsys(block1)

        # Add connections
        source_con = con.connection(self.so, 'out1', heat_source, 'in1', p=2.9,
                                    # v=0.8/3600,
                                    fluid={'water': 1})
        feed_con = con.connection(heat_source, 'out1', block1.inlet, 'in1')
        back_con = con.connection(block1.outlet, 'out1', self.si, 'in1', h=con.ref(source_con, 1, 0))
        self.nw.add_conns(source_con, feed_con, back_con)

        self.set_busses()

        self.nw.solve('design')
        self.nw.save('tmp')

