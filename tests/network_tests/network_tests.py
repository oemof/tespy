# -*- coding: utf-8

from nose.tools import eq_, raises

from tespy.connections import connection
from tespy.components import (basics, heat_exchangers, nodes, piping,
                              turbomachinery)

from tespy.networks.networks import network
from tespy.tools.helpers import TESPyNetworkError

import numpy as np
import shutil


class network_tests:

    def setup(self):
        self.nw = network(['CH4'],
                              T_unit='C', p_unit='bar', v_unit='m3 / s')
        self.source = basics.source('source')
        self.sink = basics.sink('sink')

    def test_network_delete_conns(self):
        """
        Test deleting a network's connection.
        """
        a = connection(self.source, 'out1', self.sink, 'in1')
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
        IF = basics.subsystem_interface('IF')
        a = connection(self.source, 'out1', self.sink, 'in1')
        self.nw.add_conns(a)
        self.nw.solve('design', init_only=True)
        self.nw.save('tmp')
        msg = ('After the network check, the .checked-property must be True.')
        eq_(self.nw.checked, True, msg)

        self.nw.del_conns(a)
        a = connection(self.source, 'out1', IF, 'in1')
        b = connection(IF, 'out1', self.sink, 'in1')
        self.nw.add_conns(a, b)
        self.nw.solve('design', init_path='tmp', init_only=True)
        msg = ('After the network check, the .checked-property must be True.')
        eq_(self.nw.checked, True, msg)

        shutil.rmtree('./tmp', ignore_errors=True)


@raises(TESPyNetworkError)
def offdesign_TESPyNetworkError(nw, **kwargs):
    nw.solve('offdesign', kwargs)


def test_network_missing_data_in_design_case_files():
    """
    Test for missing data in design case files.
    """
    nw = network(['water'])
    source = basics.source('source')
    pipe = piping.pipe('pipe', Q=0, pr=0.95, design=['pr'], offdesign=['zeta'])
    sink = basics.sink('sink')
    a = connection(source, 'out1', pipe, 'in1', m=1, p=1e5, T=293.15,
                       fluid={'water': 1})
    b = connection(pipe, 'out1', sink, 'in1')
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

    offdesign_TESPyNetworkError(nw, design_path='tmp2', init_only=True)

    shutil.rmtree('./tmp', ignore_errors=True)
    shutil.rmtree('./tmp2', ignore_errors=True)


def test_network_missing_data_in_individual_design_case_file():
    """
    Test for missing data in individual design case files.
    """
    nw = network(['water'])
    source = basics.source('source')
    pipe = piping.pipe('pipe', Q=0, pr=0.95, design=['pr'], offdesign=['zeta'])
    sink = basics.sink('sink')
    a = connection(source, 'out1', pipe, 'in1', m=1, p=1e5, T=293.15,
                       fluid={'water': 1})
    b = connection(pipe, 'out1', sink, 'in1', design_path='tmp2')
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

    offdesign_TESPyNetworkError(nw, design_path='tmp', init_only=True)

    shutil.rmtree('./tmp', ignore_errors=True)
    shutil.rmtree('./tmp2', ignore_errors=True)


class test_network_individual_offdesign:

    def setup_network_individual_offdesign(self):
        """
        Set up network for individual offdesign tests.
        """
        self.nw = network(['H2O'], T_unit='C', p_unit='bar', v_unit='m3 / s')

        so = basics.source('source')
        sp = nodes.splitter('splitter', num_out=2)
        self.pump1 = turbomachinery.pump('pump 1')
        self.sc1 = heat_exchangers.solar_collector('collector field 1')
        v1 = piping.valve('valve1')
        self.pump2 = turbomachinery.pump('pump 2')
        self.sc2 = heat_exchangers.solar_collector('collector field 2')
        v2 = piping.valve('valve2')
        me = nodes.merge('merge', num_in=2)
        si = basics.sink('sink')

        self.pump1.set_attr(eta_s=0.8, design=['eta_s'],
                            offdesign=['eta_s_char'])
        self.pump2.set_attr(eta_s=0.8, design=['eta_s'],
                            offdesign=['eta_s_char'])
        self.sc1.set_attr(pr=0.95, lkf_lin=3.33, lkf_quad=0.011, A=1252, E=700,
                          Tamb=20, eta_opt=0.92, design=['pr'],
                          offdesign=['zeta'])
        self.sc2.set_attr(pr=0.95, lkf_lin=3.5, lkf_quad=0.011, A=700, E=800,
                          Tamb=20, eta_opt=0.92, design=['pr'],
                          offdesign=['zeta'])

        fl = {'H2O': 1}
        inlet = connection(so, 'out1', sp, 'in1', T=50, p=3, fluid=fl)
        outlet = connection(me, 'out1', si, 'in1', p=3)

        sp_p1 = connection(sp, 'out1', self.pump1, 'in1')
        p1_sc1 = connection(self.pump1, 'out1', self.sc1, 'in1')
        self.sc1_v1 = connection(self.sc1, 'out1', v1, 'in1', p=3.1, T=90)
        v1_me = connection(v1, 'out1', me, 'in1')

        self.sp_p2 = connection(sp, 'out2', self.pump2, 'in1')
        self.p2_sc2 = connection(self.pump2, 'out1', self.sc2, 'in1')
        self.sc2_v2 = connection(self.sc2, 'out1', v2, 'in1', p=3.1, m=0.1)
        v2_me = connection(v2, 'out1', me, 'in2')

        self.nw.add_conns(inlet, outlet, sp_p1, p1_sc1, self.sc1_v1, v1_me,
                          self.sp_p2, self.p2_sc2, self.sc2_v2, v2_me)

    def test_individual_design_path_on_connections_and_components(self):
        """
        Test individual design path specification.
        """
        self.setup_network_individual_offdesign()
        self.nw.solve('design')
        self.sc2_v2.set_attr(m=0)
        self.nw.solve('design')
        self.nw.save('design1')
        v1_design = self.sc1_v1.v.val_SI
        zeta_sc1_design = self.sc1.zeta.val

        self.sc2_v2.set_attr(T=95, m=np.nan)
        self.sc1_v1.set_attr(T=np.nan, m=0.001)
        self.nw.solve('design')
        self.nw.save('design2')
        v2_design = self.sc2_v2.v.val_SI
        zeta_sc2_design = self.sc2.zeta.val

        self.sc1_v1.set_attr(m=np.nan)
        self.sc1_v1.set_attr(design=['T'], offdesign=['v'], state='l')
        self.sc2_v2.set_attr(design=['T'], offdesign=['v'], state='l')

        self.sc2.set_attr(design_path='design2')
        self.pump2.set_attr(design_path='design2')
        self.sp_p2.set_attr(design_path='design2')
        self.p2_sc2.set_attr(design_path='design2')
        self.sc2_v2.set_attr(design_path='design2')
        self.nw.solve('offdesign', design_path='design1')

        self.sc1.set_attr(E=500)
        self.sc2.set_attr(E=950)

        self.nw.solve('offdesign', design_path='design1')
        self.sc2_v2.set_attr(design_path=np.nan)

        # volumetric flow comparison
        msg = ('Design path was set to None, is ' +
               str(self.sc2_v2.design_path) + '.')
        eq_(None, self.sc2_v2.design_path, msg)

        # volumetric flow comparison
        msg = ('Value of volumetric flow must be ' + str(v1_design) + ', is ' +
               str(self.sc1_v1.v.val_SI) + '.')
        eq_(round(v1_design, 5), round(self.sc1_v1.v.val_SI, 5), msg)

        msg = ('Value of volumetric flow must be ' + str(v2_design) + ', is ' +
               str(self.sc2_v2.v.val_SI) + '.')
        eq_(round(v2_design, 5), round(self.sc2_v2.v.val_SI, 5), msg)

        # zeta value of solar collector comparison
        msg = ('Value of zeta must be ' + str(zeta_sc1_design) + ', is ' +
               str(self.sc1.zeta.val) + '.')
        eq_(round(zeta_sc1_design, 0), round(self.sc1.zeta.val, 0), msg)

        msg = ('Value of zeta must be ' + str(zeta_sc2_design) + ', is ' +
               str(self.sc2.zeta.val) + '.')
        eq_(round(zeta_sc2_design, 0), round(self.sc2.zeta.val, 0), msg)

        shutil.rmtree('./design1', ignore_errors=True)
        shutil.rmtree('./design2', ignore_errors=True)
