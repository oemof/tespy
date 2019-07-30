import sys
sys.path.append('.\\.\\')
from tespy import nwk, con, cmp
from tespy.tools.timeSeries import TimeSeries as tS


class TestEnvironment:
    def __init__(self):
        self.nw = nwk.network(fluids=['water'], T_unit='C', p_unit='bar', h_unit='J / kg')
        self.so = cmp.source('source')
        self.si = cmp.sink('sink')
        self.house = cmp.heat_exchanger_simple('house')
        self.heat_losses = con.bus('network losses')
        self.heat_consumer = con.bus('network consumer')
        self.heat_demand = -9.5e5
        self.house.set_attr(Q=self.heat_demand, pr=.98)
        self.source_temp = 80
        self.tamb = 10  # °C

        # Components
        main_pipe = cmp.district_heating_pipe('pipe',
                                              ks=1e-5,
                                              L=100,
                                              DN_type=20,
                                              lambda_ins=.03,
                                              lambda_soil=1.2,
                                              depth=.6,
                                              dist=.2,
                                              Tamb=self.tamb)
        heat_source = cmp.heat_exchanger_simple('heat_source', Q=4e4, pr=1)

        # Connections
        source_con = con.connection(self.so, 'out1', heat_source, 'in1', p=2.9, T=self.source_temp, #TODO
                                    v=0.8 / 3600,
                                    fluid={'water': 1})
        feed_con = con.connection(heat_source, 'out1', main_pipe, 'in1')
        feed_house_con = con.connection(main_pipe, 'out1', self.house, 'in1')
        back_house_con = con.connection(self.house, 'out1', main_pipe, 'in2')
        back_con = con.connection(main_pipe, 'out2', self.si, 'in1')
        self.nw.add_conns(source_con, feed_con, feed_house_con, back_con, back_house_con)

        self.set_busses()

    def set_busses(self):
        self.nw.check_network()
        for comp in self.nw.comps.index:
            if 'house' in comp.component():
                self.heat_losses.add_comps({'c': comp})
            if (isinstance(comp, cmp.heat_exchanger_simple) and
                    not isinstance(comp, cmp.pipe)):
                self.heat_consumer.add_comps({'c': comp})
        self.nw.add_busses(self.heat_losses, self.heat_consumer)

    def set_attributes(self, demand, tamb):
        self.heat_demand = demand
        self.tamb = tamb
        for comp in self.nw.comps.index:
            if isinstance(comp, cmp.district_heating_pipe):
                comp.set_attr(Tamb=self.tamb)
            if isinstance(comp, cmp.heat_exchanger_simple):
                comp.set_attr(Q=self.heat_demand)


timeSeries = tS()
timeSeries.run_time_series(TestEnvironment())
