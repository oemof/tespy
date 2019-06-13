import sys
import math
import logging
import unittest
sys.path.append('.\\.\\')
from tespy import nwk, con, cmp
from tespy.tools import logger
import shutil

logger.define_logging(log_path=True, log_version=True,
                      screen_level=logging.INFO, file_level=logging.DEBUG)


class DistrictHeatingPipeTests(unittest.TestCase):
    def setUp(self):
        print("\n\n\nRunning Test...\n\n\n")
        self.nw = nwk.network(fluids=['water'], T_unit='C', p_unit='bar', h_unit='J / kg')
        self.so = cmp.source('source')
        self.si = cmp.sink('sink')
        self.house = cmp.heat_exchanger_simple('house')
        self.heat_losses = con.bus('network losses')
        self.heat_consumer = con.bus('network consumer')

    def tearDown(self):
        shutil.rmtree('./tmp', ignore_errors=True)

    def set_busses(self):
        self.nw.check_network()
        for comp in self.nw.comps.index:
            if isinstance(comp, cmp.district_heating_pipe):
                comp.set_attr(Tamb=0)
                self.heat_losses.add_comps({'c': comp})
            if (isinstance(comp, cmp.heat_exchanger_simple) and
                    not isinstance(comp, cmp.pipe)):
                self.heat_consumer.add_comps({'c': comp})
        self.nw.add_busses(self.heat_losses, self.heat_consumer)

    def test_different_t_ambs(self):
        self.house.set_attr(Q=-5e5, pr=.98)

        pipe_diameter = .168
        t_ground = 10
        pipe = cmp.district_heating_pipe('pipe',
                                         ks=1e-5,
                                         L=60,
                                         D=pipe_diameter,
                                         lambda_ins=0.035,
                                         lambda_soil=1.6,
                                         depth=.6,
                                         Dout=.35,
                                         dist=0.15,
                                         Tamb=t_ground)

        # %% connections
        feed_con = con.connection(self.so, 'out1', pipe, 'in1', T=80, p=5,
                                  v=(pipe_diameter / 2) ** 2 * math.pi * 1,
                                  fluid={'water': 1})
        feed_house_con = con.connection(pipe, 'out1', self.house, 'in1',
                                        )
        back_house_con = con.connection(self.house, 'out1', pipe, 'in2')
        back_con = con.connection(pipe, 'out2', self.si, 'in1')
        # design case: 0 °C ambient temperature
        self.nw.add_conns(feed_con, feed_house_con, back_con, back_house_con)

        self.set_busses()

        self.nw.solve('design')
        self.nw.print_results()
        self.nw.save('tmp')

        print('Heat demand consumer:', self.heat_consumer.P.val)
        print('network losses at 0 °C outside temperature (design):', self.heat_losses.P.val)
        print('relative network losses at 0 °C outside temperature:', self.heat_losses.P.val / self.heat_consumer.P.val * 100,
              '%')

        print('Pipe losses: ', )

        # # offdesign case: 10 °C ambient temperature

        for comp in self.nw.comps.index:
            if isinstance(comp, cmp.district_heating_pipe):
                comp.set_attr(Tamb=10)

        self.nw.solve('offdesign', design_path='tmp')
        self.nw.print_results()
        print('Heat demand consumer:', self.heat_consumer.P.val)
        print('network losses at 10 °C outside temperature:', self.heat_losses.P.val)
        print('relative network losses at 10 °C outside temperature:', self.heat_losses.P.val / self.heat_consumer.P.val * 100,
              '%')

        for comp in self.nw.comps.index:
            if isinstance(comp, cmp.district_heating_pipe):
                comp.set_attr(Tamb=20)

        self.nw.solve('offdesign', design_path='tmp')
        self.nw.print_results()
        print('Heat demand consumer:', self.heat_consumer.P.val)
        print('network losses at 20 °C outside temperature:', self.heat_losses.P.val)
        print('relative network losses at 20 °C outside temperature:', self.heat_losses.P.val / self.heat_consumer.P.val * 100,
              '%')

        # # offdesign case: -10 °C ambient temperature

        for comp in self.nw.comps.index:
            if isinstance(comp, cmp.district_heating_pipe):
                comp.set_attr(Tamb=-10)

        self.nw.solve('offdesign', design_path='tmp')
        self.nw.print_results()
        print('Heat demand consumer:', self.heat_consumer.P.val)
        print('network losses at -10 °C outside temperature:', self.heat_losses.P.val)
        print('relative network losses at -10 °C outside temperature:', self.heat_losses.P.val / self.heat_consumer.P.val * 100,
              '%')

        # # offdesign case: -20 °C ambient temperature

        for comp in self.nw.comps.index:
            if isinstance(comp, cmp.district_heating_pipe):
                comp.set_attr(Tamb=-20)

        self.nw.solve('offdesign', design_path='tmp')
        self.nw.print_results()
        print('Heat demand consumer:', self.heat_consumer.P.val)
        print('network losses at -20 °C outside temperature:', self.heat_losses.P.val)
        print('relative network losses at -20 °C outside temperature:', self.heat_losses.P.val / self.heat_consumer.P.val * 100,
              '%')

    def test_in80out50c(self):
        self.house.set_attr(Q=-9.5e5, pr=.98)
        pipe_diameter = .1
        t_ground = 10  # °C
        pipe = cmp.district_heating_pipe('pipe',
                                         ks=1e-5,
                                         L=1,
                                         D=pipe_diameter,
                                         lambda_ins=0.035,
                                         lambda_soil=1.6,
                                         depth=.6,
                                         Dout=.5,
                                         dist=0.15,
                                         Tamb=t_ground)
        feed_con = con.connection(self.so, 'out1', pipe, 'in1', T=80, p=5,
                                  v=(pipe_diameter / 2) ** 2 * math.pi * 1,
                                  fluid={'water': 1})
        feed_house_con = con.connection(pipe, 'out1', self.house, 'in1',
                                        )
        back_house_con = con.connection(self.house, 'out1', pipe, 'in2')
        back_con = con.connection(pipe, 'out2', self.si, 'in1')
        self.nw.add_conns(feed_con, feed_house_con, back_con, back_house_con)

        self.set_busses()

        self.nw.solve('design')
        self.nw.save('tmp')

        self.assertGreater(pipe.Q_u.val, -15, "Heat losses in feed pipe expected to be between 10 and 15 W/m for this setting.")
        self.assertLess(pipe.Q_u.val, -10, "Heat losses in feed pipe expected to be between 10 and 15 W/m for this setting.")
        self.nw.print_results()
