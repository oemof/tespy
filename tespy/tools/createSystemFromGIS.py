# -*- coding: utf-8

""" Class for importing GIS exported district heating networks from a .CSV file. The file has to contain the columns
'Shape_Length', 'pipe_ID', 'prior_ID', 'passage', 'capacity' and 'DN'. """

import pandas as pd
import os
import logging
from tespy import nwk, cmp, con
from tespy.tools import logger, helpers

logger.define_logging(log_path=True, log_version=True,
                      screen_level=logging.INFO, file_level=logging.DEBUG)


class CreateSystemFromGIS:
    def __init__(self, file_path, file_name):
        self.pathToFile = file_path
        self.filename = file_name
        self.so = cmp.source('source')
        self.si = cmp.sink('sink')
        self.heat_losses = con.bus('network losses')
        self.heat_consumer = con.bus('network consumer')
        self.pipes = pd.DataFrame()
        self.nw = nwk.network(fluids=['water'], T_unit='C', p_unit='Pa', h_unit='J / kg', v_unit='m3 / h')
        self.load_file()

        # Definition of a single heat source for the system. It is directly connected to the network's source
        self.heat_source = cmp.heat_exchanger_simple(label='Heat source',
                                                     Q=self.pipes['capacity'][0]*1.3,
                                                     pr=1)
        self.add_pipe_to_network([], [], 0)
        self.nw.set_printoptions(print_level='info', iterinfo=True)
        # Network's back temperature.TODO Shall be overwritten in a loop with last frame's temperature at sink
        self.back_pipe_temp = 80
        self.source_con = con.connection(self.so, 'out1', self.heat_source, 'in1',
                                         state='l',
                                         p=5e5,
                                         v=2,
                                         fluid={'water': 1})
        self.source_con.set_attr(T=self.back_pipe_temp)
        self.nw.add_conns(self.source_con)
        self.set_busses()
        print(200 * '#' + '\nDesign calculation with ' + str(self.back_pipe_temp) + 'C heat source:\n' + 200 * '#')
        self.nw.solve(mode='design')   # , init_path='../../../init/' + self.filename)
        self.nw.print_results()
        print('\nTotal heat demand consumers: ', -round(self.heat_consumer.P.val), 'W')
        print('Total pipe heat losses: ', -round(self.heat_losses.P.val), 'W\n')
        for comp in self.nw.comps.index:
            # if isinstance(comp, cmp.valve) and 'output' in comp.label:
            #     comp.set_attr(zeta=(0.2 if 'Straight' in comp.label else 2))
            # elif isinstance(comp, cmp.valve) and 'input' in comp.label:
            #     comp.set_attr(zeta=(1 if 'Straight' in comp.label else 1.5))
            # if isinstance(comp, cmp.valve):
            #     comp.get_attr('zeta').set_attr(is_set=True)
            #     comp.get_attr('pr').set_attr(is_set=False)

            # Find lowest consumer Temperature
            lowest_consumer_temp = 1000
            if isinstance(comp, cmp.heat_exchanger_simple) and 'house' in comp.label:
                temp = helpers.T_mix_ph(comp.inl[0].to_flow()) - 273.15
                print(comp.label + ' temperature \t{:.2f}°C'.format(temp))
                if temp < lowest_consumer_temp:
                    lowest_consumer_temp = temp
                    lowest_temp_house = comp.label
        print('Lowest consumer temperature: ' + '{:.2f}'.format(lowest_consumer_temp) +
              '°C in ' + lowest_temp_house)
        self.nw.save('../../../init/' + self.filename)
        ########################
        # Offdesign Calculations
        ########################
        # self.back_pipe_temp = 75
        # self.source_con.set_attr(T=self.back_pipe_temp)
        # print(200 * '#' + '\nOffdesign calculation with ' + str(self.back_pipe_temp) + 'C heat source:\n' + 200 * '#')
        # self.nw.solve('offdesign', design_path='../../../init/' + self.filename)
        # self.nw.print_results()
        # print('Total heat demand consumers: ', -round(self.heat_consumer.P.val), 'W')
        # print('Total pipe heat losses: ', -round(self.heat_losses.P.val), 'W')
        # self.nw.save('../../../init/' + self.filename + '2')
        #
        # self.back_pipe_temp = 95
        # self.source_con.set_attr(T=self.back_pipe_temp)
        # print(200 * '#' + '\nOffdesign calculation with ' + (self.back_pipe_temp) + 'C heat source:\n' + 200 * '#')
        # self.nw.solve('offdesign', design_path='../../../init/' + self.filename)
        # self.nw.print_results()
        # print('Total heat demand consumers: ', -round(self.heat_consumer.P.val), 'W')
        # print('Total pipe heat losses: ', -round(self.heat_losses.P.val), 'W')

    def load_file(self):
        self.pipes = pd.read_csv(os.path.join(self.pathToFile, self.filename), sep=';')
        self.pipes['capacity'] *= 1e3
        logging.debug('Loaded Data: \n\n' + str(self.pipes))

    # Recursively build a Tespy network from the data imported
    def add_pipe_to_network(self, parentout, parentin, cmp_id):
        children = pd.DataFrame(self.pipes.loc[self.pipes['prior_ID'] == cmp_id])
        num_children = children.shape[0]
        str_id = '{:02d}'.format(cmp_id)
        if cmp_id == 0:
            self.add_pipe_to_network([self.heat_source, 'out1'], [self.si, 'in1'], children['pipe_ID'][0])
        else:
            pipe = self.pipes.loc[self.pipes['pipe_ID'] == cmp_id].iloc[0]

            # Create pipe component
            p = cmp.district_heating_pipe(label='pipe ID ' + str_id,
                                          ks=1e-5,
                                          L=pipe['Shape_Length'],
                                          DN_type=pipe['DN'],
                                          lambda_ins=.03,
                                          lambda_soil=1.2,
                                          depth=.6,
                                          dist=.2,
                                          Tamb=10)
            logging.debug('Created component ' + p.get_attr('label'))

            feed_in = con.connection(parentout[0], parentout[1], p, 'in1', state='l')
            feed_out = con.connection(p, 'out2', parentin[0], parentin[1], state='l')
            self.nw.add_conns(feed_in, feed_out)

            # Create consumer if pipe is a house connector
            if num_children == 0:
                consumer = cmp.heat_exchanger_simple(label='house ID ' + str_id,
                                                     Q=pipe['capacity'] * -1,  # has to be overwritten by demand data
                                                     pr=.995)      # TODO: get realistic pressure drop value or define zeta
                pipe_valve = con.connection(p, 'out1', consumer, 'in1', state='l')
                consumer_out = con.connection(consumer, 'out1', p, 'in2', state='l')
                self.nw.add_conns(pipe_valve, consumer_out)
                logging.debug('Created consumer ' + consumer.get_attr('label'))

            # If there is only one connection to another pipe
            elif num_children == 1:
                self.add_pipe_to_network([p, 'out1'], [p, 'in2'], children['pipe_ID'].iloc[0])
                # TODO: add valve with 90 deg zeta

            # Create splitters and mergers if multiple connecting pipes are present
            else:
                splitter = cmp.splitter(label='Splitter behind ID ' + str_id,
                                        num_out=num_children)
                splitter_in = con.connection(p, 'out1', splitter, 'in1', state='l')
                merger = cmp.merge(label='Merger at ID: ' + str_id,
                                   num_in=num_children)
                merge_out = con.connection(merger, 'out1', p, 'in2', state='l')
                self.nw.add_conns(splitter_in, merge_out)
                for i in range(num_children):
                    child = children.astype(int).iloc[i]
                    is_straight = child['passage'] == 0
                    valve_out = cmp.valve(label=('Straight' if is_straight else '') + 'Valve ' + str(i) +
                                          ' at splitter output ID' + str_id,
                                          pr=1)
                                          # zeta=(0.2 if is_straight else 2)/(p.D.val**4))
                    valve_in = cmp.valve(label=('Straight' if is_straight else '') + 'Valve ' + str(i) +
                                         ' at merger input ID' + str_id,
                                         pr=1)
                                         # zeta=(1 if is_straight else 1.5)/(p.D.val**4))
                    con_valve_out = con.connection(splitter, 'out' + str(i+1), valve_out, 'in1', state='l')
                    con_valve_in = con.connection(valve_in, 'out1', merger, 'in' + str(i+1), state='l')
                    self.nw.add_conns(con_valve_out, con_valve_in)
                    self.add_pipe_to_network([valve_out, 'out1'],
                                             [valve_in, 'in1'],
                                             child['pipe_ID'])

    def set_busses(self):
        self.nw.check_network()
        for comp in self.nw.comps.index:
            if isinstance(comp, cmp.district_heating_pipe):
                self.heat_losses.add_comps({'c': comp})
            if isinstance(comp, cmp.heat_exchanger_simple) and 'source' not in comp.get_attr('label'):
                self.heat_consumer.add_comps({'c': comp})
        self.nw.add_busses(self.heat_losses, self.heat_consumer)


CreateSystemFromGIS(file_path='../../../TestFiles', file_name='191016_kleinesGebiet-4haus-nurInt.csv')
