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



def test_network_stuff()

    def setup(self):
        self.nw = nwk.network(['INCOMP::DowQ', 'H2O', 'NH3', 'N2',
                               'O2', 'Ar', 'CO2', 'CH4'],
                              T_unit='C', p_unit='bar', v_unit='m3 / s')
        self.source = cmp.source('source')
        self.sink = cmp.sink('sink')
        c1 = con.connection(self.source, 'out1', instance, 'in1')
        c2 = con.connection(cmp.source('fuel'), 'out1', instance, 'in2')
        c3 = con.connection(instance, 'out1', self.sink, 'in1')
        self.nw.add_conns(c1, c2, c3)
        return c1, c2, c3

    def setup_network_individual_offdesign(self):
        """
        Set up network for individual offdesign tests.
        """
        so = cmp.source('source')
        sp = cmp.splitter('splitter', num_out=2)
        self.pump1 = cmp.pump('pump 1')
        self.sc1 = cmp.solar_collector('collector field 1')
        v1 = cmp.valve('valve1')
        self.pump2 = cmp.pump('pump 2')
        self.sc2 = cmp.solar_collector('collector field 2')
        v2 = cmp.valve('valve2')
        me = cmp.merge('merge', num_in=2)
        si = cmp.sink('sink')

        self.pump1.set_attr(eta_s=0.8, design=['eta_s'],
                            offdesign=['eta_s_char'])
        self.pump2.set_attr(eta_s=0.8, design=['eta_s'],
                            offdesign=['eta_s_char'])
        self.sc1.set_attr(pr=0.95, lkf_lin=3.33, lkf_quad=0.011, A=1252, E=700,
                          Tamb=20, design=['pr'], offdesign=['zeta'])
        self.sc2.set_attr(pr=0.95, lkf_lin=3.5, lkf_quad=0.011, A=700, E=800,
                          Tamb=20, design=['pr'], offdesign=['zeta'])

        fl = {'N2': 0, 'O2': 0, 'Ar': 0, 'INCOMP::DowQ': 0,
              'H2O': 1, 'NH3': 0, 'CO2': 0, 'CH4': 0}

        inlet = con.connection(so, 'out1', sp, 'in1', T=50, p=3, fluid=fl)
        outlet = con.connection(me, 'out1', si, 'in1', p=3)

        sp_p1 = con.connection(sp, 'out1', self.pump1, 'in1')
        p1_sc1 = con.connection(self.pump1, 'out1', self.sc1, 'in1')
        self.sc1_v1 = con.connection(self.sc1, 'out1', v1, 'in1', p=3.1, T=90)
        v1_me = con.connection(v1, 'out1', me, 'in1')

        self.sp_p2 = con.connection(sp, 'out2', self.pump2, 'in1')
        self.p2_sc2 = con.connection(self.pump2, 'out1', self.sc2, 'in1')
        self.sc2_v2 = con.connection(self.sc2, 'out1', v2, 'in1', p=3.1, m=0.1)
        v2_me = con.connection(v2, 'out1', me, 'in2')

        self.nw.add_conns(inlet, outlet, sp_p1, p1_sc1, self.sc1_v1, v1_me,
                          self.sp_p2, self.p2_sc2, self.sc2_v2, v2_me)

    def test_individual_design_path_on_connections_and_components(self):

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
