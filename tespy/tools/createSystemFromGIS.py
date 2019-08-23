import pandas as pd
import os
import logging
from pathlib import Path
from tespy import nwk, cmp, con
from tespy.tools import logger

logger.define_logging(log_path=True, log_version=True,
                      screen_level=logging.INFO, file_level=logging.DEBUG)


class CreateSystemFromGIS:
    def __init__(self):
        self.pathToFile = ""
        self.filename = ""
        self.so = cmp.source('source')
        self.si = cmp.sink('sink')
        self.heat_losses = con.bus('network losses')
        self.heat_consumer = con.bus('network consumer')
        self.pipes = pd.DataFrame()
        self.nw = nwk.network(fluids=['water'], T_unit='C', p_unit='bar', h_unit='J / kg',
                              T_range=[0, 130], p_range=[0.05, 15])
        self.load_file()
        self.heat_source = cmp.heat_exchanger_simple(label='Heat source',
                                                     Q=151e3,
                                                     pr=1)
        self.source_con = con.connection(self.so, 'out1', self.heat_source, 'in1',
                                         state='l',
                                         p=8,
                                         T=80,
                                         v=7/3600,
                                         fluid={'water': 1})
        # sink_con = con.connection(self.heat_source, 'out2', self.si, 'in1', h=con.ref(self.source_con, 1, 0))
        self.nw.add_conns(self.source_con)
        self.add_pipe_to_network([], [], 0)
        self.nw.check_network()
        self.nw.set_printoptions(print_level='info', iterinfo=True)
        self.nw.solve('design')
        self.nw.print_results()

    def load_file(self):
        self.filename = 'testGrid.csv'
        self.pathToFile = str(Path.home())
        self.pipes = pd.read_csv(os.path.join(self.pathToFile, self.filename), sep=';')
        logging.debug('Loaded Data: \n\n' + str(self.pipes))

    def add_pipe_to_network(self, parentout, parentin, cmp_id):
        children = pd.DataFrame(self.pipes.loc[self.pipes['Vorgaenger_ID'] == cmp_id])
        num_children = children.shape[0]

        if cmp_id == 0:
            self.add_pipe_to_network([self.heat_source, 'out1'], [self.si, 'in1'], children['FWL_ID'][0])
        else:
            pipe = self.pipes.loc[self.pipes['FWL_ID'] == cmp_id].iloc[0]

            # Create pipe component
            p = cmp.district_heating_pipe(label='pipe ID ' + str(cmp_id),
                                          ks=1e-5,
                                          L=pipe['SHAPE_Length'],
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
            if pipe['HA'] == 1:
                consumer = cmp.heat_exchanger_simple(label='house ID ' + str(cmp_id),
                                                     Q=pipe['Leistung '] * -1e3,  # has to be overwritten by demand data
                                                     pr=1)      # TODO: get realistic pressure drop value or define zeta
                valve = cmp.valve(label='valve at house ID ' + str(cmp_id), pr=1)
                pipe_valve = con.connection(p, 'out1', valve, 'in1', state='l')
                valve_consumer_in = con.connection(valve, 'out1', consumer, 'in1', state='l')
                consumer_out = con.connection(consumer, 'out1', p, 'in2', state='l')
                self.nw.add_conns(pipe_valve, valve_consumer_in, consumer_out)
                logging.debug('Created consumer ' + consumer.get_attr('label'))

            # If there is only one connection to another pipe
            elif num_children == 1:
                self.add_pipe_to_network([p, 'out1'], [p, 'in2'], children['FWL_ID'].iloc[0])

            # Create splitters and mergers if multiple connecting pipes are present
            else:
                splitter = cmp.splitter(label='Splitter behind ID ' + str(cmp_id),
                                        num_out=num_children)
                splitter_in = con.connection(p, 'out1', splitter, 'in1', state='l')
                merger = cmp.merge(label='Merger at ID: ' + str(cmp_id),
                                   num_in=num_children)
                merge_out = con.connection(merger, 'out1', p, 'in2', state='l')
                self.nw.add_conns(splitter_in, merge_out)
                for i in range(num_children):
                    child = children.iloc[i]
                    is_straight = child['Uebergang_ID'] == 0
                    valve_out = cmp.valve(label='Valve ' + str(i) + ' at splitter output ID' + str(cmp_id),
                                          pr=1)
                                          # zeta=0.2 if is_straight else 2)
                    valve_in = cmp.valve(label='Valve ' + str(i) + ' at merger input ID' + str(cmp_id),
                                         pr=1)
                                         # zeta=1 if is_straight else 1.5)
                    con_valve_out = con.connection(splitter, 'out' + str(i+1), valve_out, 'in1', state='l')
                    con_valve_in = con.connection(valve_in, 'out1', merger, 'in' + str(i+1), state='l')
                    self.nw.add_conns(con_valve_out, con_valve_in)
                    self.add_pipe_to_network([valve_out, 'out1'],
                                             [valve_in, 'in1'],
                                             child['FWL_ID'])
# Uebergang_ID == 0: Durchgang


test = CreateSystemFromGIS()
