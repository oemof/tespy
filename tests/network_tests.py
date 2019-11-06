# -*- coding: utf-8

from nose.tools import eq_

from tespy import nwk, cmp, con, hlp
import numpy as np
import shutil


class network_tests:

    def setup(self):
        self.nw = nwk.network(['CH4'],
                              T_unit='C', p_unit='bar', v_unit='m3 / s')
        self.source = cmp.source('source')
        self.sink = cmp.sink('sink')

    def test_network_delete_conns(self):
        """
        Test deleting a network's connection.
        """
        a = con.connection(self.source, 'out1', self.sink, 'in1')
        self.nw.add_conns(a)
        self.nw.check_network()
        msg = ('After the network check, the .checked-property must be True.')
        eq_(self.nw.checked, True, msg)

        self.nw.del_conns(a)
        msg = ('A connection has been deleted, the network consistency check '
               'must be repeated (.checked-property must be False).')
        eq_(self.nw.checked, False, msg)

    def test_network_missing_connection_in_init_path(self):
        """
        Test debug message for missing connection in init_path.
        """
        IF = cmp.subsys_interface('IF')
        a = con.connection(self.source, 'out1', self.sink, 'in1')
        self.nw.add_conns(a)
        self.nw.solve('design', init_only=True)
        self.nw.save('tmp')
        msg = ('After the network check, the .checked-property must be True.')
        eq_(self.nw.checked, True, msg)

        self.nw.del_conns(a)
        a = con.connection(self.source, 'out1', IF, 'in1')
        b = con.connection(IF, 'out1', self.sink, 'in1')
        self.nw.add_conns(a, b)
        self.nw.solve('design', init_path='tmp', init_only=True)
        msg = ('After the network check, the .checked-property must be True.')
        eq_(self.nw.checked, True, msg)

        shutil.rmtree('./tmp', ignore_errors=True)


def test_network_missing_data_in_design_case_files():
    nw = nwk.network(['water'])
    source = cmp.source('source')
    pipe = cmp.pipe('pipe', Q=0, pr=0.95, design=['pr'], offdesign=['zeta'])
    sink = cmp.sink('sink')
    a = con.connection(source, 'out1', pipe, 'in1', m=1, p=1e5, T=293.15,
                       fluid={'water': 1})
    b = con.connection(pipe, 'out1', sink, 'in1')
    nw.add_conns(a, b)
    nw.solve('design')
    nw.save('tmp')
    nw.save('tmp2')

    inputs = open('./tmp/conn.csv')
    all_lines = inputs.readlines()
    all_lines.pop(len(all_lines) - 1)
    inputs.close()

    with open('./tmp2/conn.csv', 'w') as out:
        for line in all_lines:
            out.write(line.strip() + '\n')

    try:
        nw.solve('offdesign', design_path='tmp2')
    except hlp.TESPyNetworkError:
        pass

    shutil.rmtree('./tmp', ignore_errors=True)
    shutil.rmtree('./tmp2', ignore_errors=True)


def test_network_missing_data_in_individual_design_case_file():
    nw = nwk.network(['water'])
    source = cmp.source('source')
    pipe = cmp.pipe('pipe', Q=0, pr=0.95, design=['pr'], offdesign=['zeta'])
    sink = cmp.sink('sink')
    a = con.connection(source, 'out1', pipe, 'in1', m=1, p=1e5, T=293.15,
                       fluid={'water': 1})
    b = con.connection(pipe, 'out1', sink, 'in1', design_path='tmp2')
    nw.add_conns(a, b)
    nw.solve('design')
    nw.save('tmp')
    nw.save('tmp2')

    inputs = open('./tmp/conn.csv')
    all_lines = inputs.readlines()
    all_lines.pop(len(all_lines) - 1)
    inputs.close()

    with open('./tmp2/conn.csv', 'w') as out:
        for line in all_lines:
            out.write(line.strip() + '\n')
    try:
        nw.solve('offdesign', design_path='tmp')
    except hlp.TESPyNetworkError:
        pass

    shutil.rmtree('./tmp', ignore_errors=True)
    shutil.rmtree('./tmp2', ignore_errors=True)
