# -*- coding: utf-8

"""Module for testing network properties.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_network.py

SPDX-License-Identifier: MIT
"""

import os
import shutil

import numpy as np
from pytest import raises

from tespy.components import basics
from tespy.components import heat_exchangers
from tespy.components import nodes
from tespy.components import piping
from tespy.components import turbomachinery
from tespy.connections import connection
from tespy.networks.network_reader import load_network
from tespy.networks.networks import network
from tespy.tools.helpers import TESPyNetworkError


class TestNetworks:
    def setup_network_tests(self):
        self.nw = network(['water'], p_unit='bar', v_unit='m3 / s')
        self.source = basics.source('source')
        self.sink = basics.sink('sink')

    def offdesign_TESPyNetworkError(self, **kwargs):
        with raises(TESPyNetworkError):
            self.nw.solve('offdesign', **kwargs)

    def test_network_linear_dependency(self):
        """Test network linear dependency."""
        self.setup_network_tests()
        a = connection(self.source, 'out1', self.sink, 'in1', m=1, p=1, x=1,
                       T=280)
        self.nw.add_conns(a)
        self.nw.solve('design')
        msg = ('This test must result in a linear dependency of the jacobian '
               'matrix.')
        assert self.nw.lin_dep is True, msg

    def test_network_no_progress(self):
        """Test no convergence progress."""
        self.setup_network_tests()
        pipe = piping.pipe('pipe', pr=1, Q=-100e3)
        a = connection(self.source, 'out1', pipe, 'in1', m=1, p=1, T=280,
                       fluid={'water': 1})
        b = connection(pipe, 'out1', self.sink, 'in1')
        self.nw.add_conns(a, b)
        self.nw.solve('design')
        msg = ('This test must result in a calculation making no progress, as '
               'the pipe\'s outlet enthalpy is below fluid property range.')
        assert self.nw.progress is False, msg

    def test_network_max_iter(self):
        """Test reaching maximum iteration count."""
        self.setup_network_tests()
        pipe = piping.pipe('pipe', pr=1, Q=100e3)
        a = connection(self.source, 'out1', pipe, 'in1', m=1, p=1, T=280,
                       fluid={'water': 1})
        b = connection(pipe, 'out1', self.sink, 'in1')
        self.nw.add_conns(a, b)
        self.nw.solve('design', max_iter=2)
        msg = ('This test must result in the itercount being equal to the max '
               'iter statement.')
        assert self.nw.max_iter == self.nw.iter + 1, msg

    def test_network_delete_conns(self):
        """Test deleting a network's connection."""
        self.setup_network_tests()
        a = connection(self.source, 'out1', self.sink, 'in1')
        self.nw.add_conns(a)
        self.nw.check_network()
        msg = ('After the network check, the .checked-property must be True.')
        assert self.nw.checked is True, msg

        self.nw.del_conns(a)
        msg = ('A connection has been deleted, the network consistency check '
               'must be repeated (.checked-property must be False).')
        assert self.nw.checked is False, msg

    def test_network_missing_connection_in_init_path(self):
        """Test debug message for missing connection in init_path."""
        self.setup_network_tests()
        IF = basics.subsystem_interface('IF')
        a = connection(self.source, 'out1', self.sink, 'in1')
        self.nw.add_conns(a)
        self.nw.solve('design', init_only=True)
        self.nw.save('tmp')
        msg = ('After the network check, the .checked-property must be True.')
        assert self.nw.checked is True, msg

        self.nw.del_conns(a)
        a = connection(self.source, 'out1', IF, 'in1')
        b = connection(IF, 'out1', self.sink, 'in1')
        self.nw.add_conns(a, b)
        self.nw.solve('design', init_path='tmp', init_only=True)
        msg = ('After the network check, the .checked-property must be True.')
        assert self.nw.checked is True, msg

        shutil.rmtree('./tmp', ignore_errors=True)

    def test_network_export_no_chars_busses(self):
        """Test export of network without characteristics or busses."""
        self.setup_network_tests()
        a = connection(self.source, 'out1', self.sink, 'in1')
        self.nw.add_conns(a)
        self.nw.solve('design', init_only=True)
        self.nw.save('tmp')
        msg = ('The exported network does not contain any char_line, there '
               'must be no file char_line.csv!')
        assert os.path.isfile('tmp/components/char_line.csv') is False, msg

        msg = ('The exported network does not contain any char_map, there '
               'must be no file char_map.csv!')
        assert os.path.isfile('tmp/components/char_map.csv') is False, msg

        msg = ('The exported network does not contain any busses, there '
               'must be no file bus.csv!')
        assert os.path.isfile('tmp/components/bus.csv') is False, msg
        shutil.rmtree('./tmp', ignore_errors=True)

    def test_network_reader_no_chars_busses(self):
        """Test import of network without characteristics or busses."""
        self.setup_network_tests()
        a = connection(self.source, 'out1', self.sink, 'in1')
        self.nw.add_conns(a)
        self.nw.solve('design', init_only=True)
        self.nw.save('tmp')

        imported_nwk = load_network('tmp')
        imported_nwk.solve('design', init_only=True)
        msg = ('If the network import was successful the network check '
               'should have been successful, too, but it is not.')
        assert imported_nwk.checked is True, msg
        shutil.rmtree('./tmp', ignore_errors=True)

    def test_network_reader_deleted_chars(self):
        """Test import of network with missing characteristics."""
        self.setup_network_tests()
        comp = turbomachinery.compressor('compressor')
        a = connection(self.source, 'out1', comp, 'in1')
        b = connection(comp, 'out1', self.sink, 'in1')
        self.nw.add_conns(a, b)
        self.nw.solve('design', init_only=True)
        self.nw.save('tmp')

        # # remove char_line and char_map folders
        os.unlink('tmp/components/char_line.csv')
        os.unlink('tmp/components/char_map.csv')

        # import network with missing files
        imported_nwk = load_network('tmp')
        imported_nwk.solve('design', init_only=True)
        msg = ('If the network import was successful the network check '
               'should have been successful, too, but it is not.')
        assert imported_nwk.checked is True, msg
        shutil.rmtree('./tmp', ignore_errors=True)

    def test_network_missing_data_in_design_case_files(self):
        """Test for missing data in design case files."""
        self.setup_network_tests()
        pipe = piping.pipe('pipe', Q=0, pr=0.95, design=['pr'],
                           offdesign=['zeta'])
        a = connection(self.source, 'out1', pipe, 'in1', m=1, p=1, T=293.15,
                       fluid={'water': 1})
        b = connection(pipe, 'out1', self.sink, 'in1')
        self.nw.add_conns(a, b)
        self.nw.solve('design')
        self.nw.save('tmp')
        self.nw.save('tmp2')

        inputs = open('./tmp/connections.csv')
        all_lines = inputs.readlines()
        all_lines.pop(len(all_lines) - 1)
        inputs.close()

        with open('./tmp2/connections.csv', 'w') as out:
            for line in all_lines:
                out.write(line.strip() + '\n')

        self.offdesign_TESPyNetworkError(design_path='tmp2', init_only=True)

        shutil.rmtree('./tmp', ignore_errors=True)
        shutil.rmtree('./tmp2', ignore_errors=True)

    def test_network_missing_data_in_individual_design_case_file(self):
        """Test for missing data in individual design case files."""
        self.setup_network_tests()
        pipe = piping.pipe('pipe', Q=0, pr=0.95, design=['pr'],
                           offdesign=['zeta'])
        a = connection(self.source, 'out1', pipe, 'in1', m=1, p=1, T=293.15,
                       fluid={'water': 1})
        b = connection(pipe, 'out1', self.sink, 'in1', design_path='tmp2')
        self.nw.add_conns(a, b)
        self.nw.solve('design')
        self.nw.save('tmp')
        self.nw.save('tmp2')

        inputs = open('./tmp/connections.csv')
        all_lines = inputs.readlines()
        all_lines.pop(len(all_lines) - 1)
        inputs.close()

        with open('./tmp2/connections.csv', 'w') as out:
            for line in all_lines:
                out.write(line.strip() + '\n')

        self.offdesign_TESPyNetworkError(design_path='tmp', init_only=True)

        shutil.rmtree('./tmp', ignore_errors=True)
        shutil.rmtree('./tmp2', ignore_errors=True)

    def test_network_missing_connection_in_design_path(self):
        """Test for missing connection data in design case files."""
        self.setup_network_tests()
        pipe = piping.pipe('pipe', Q=0, pr=0.95, design=['pr'],
                           offdesign=['zeta'])
        a = connection(self.source, 'out1', pipe, 'in1', m=1, p=1, T=293.15,
                       fluid={'water': 1})
        b = connection(pipe, 'out1', self.sink, 'in1')
        self.nw.add_conns(a, b)
        self.nw.solve('design')
        self.nw.save('tmp')

        inputs = open('./tmp/connections.csv')
        all_lines = inputs.readlines()
        all_lines.pop(len(all_lines) - 1)
        inputs.close()

        with open('./tmp/connections.csv', 'w') as out:
            for line in all_lines:
                out.write(line.strip() + '\n')

        self.offdesign_TESPyNetworkError(design_path='tmp')

        shutil.rmtree('./tmp', ignore_errors=True)


class TestNetworkIndividualOffdesign:

    def setup_network_individual_offdesign(self):
        """Set up network for individual offdesign tests."""
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

        self.sp_p1 = connection(sp, 'out1', self.pump1, 'in1')
        self.p1_sc1 = connection(self.pump1, 'out1', self.sc1, 'in1')
        self.sc1_v1 = connection(self.sc1, 'out1', v1, 'in1', p=3.1, T=90)
        v1_me = connection(v1, 'out1', me, 'in1')

        self.sp_p2 = connection(sp, 'out2', self.pump2, 'in1')
        self.p2_sc2 = connection(self.pump2, 'out1', self.sc2, 'in1')
        self.sc2_v2 = connection(self.sc2, 'out1', v2, 'in1', p=3.1, m=0.1)
        v2_me = connection(v2, 'out1', me, 'in2')

        self.nw.add_conns(inlet, outlet, self.sp_p1, self.p1_sc1, self.sc1_v1,
                          v1_me, self.sp_p2, self.p2_sc2, self.sc2_v2, v2_me)

    def test_individual_design_path_on_connections_and_components(self):
        """Test individual design path specification."""
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
        assert self.sc2_v2.design_path is None, msg

        # volumetric flow comparison
        msg = ('Value of volumetric flow must be ' + str(v1_design) + ', is ' +
               str(self.sc1_v1.v.val_SI) + '.')
        assert round(v1_design, 5) == round(self.sc1_v1.v.val_SI, 5), msg

        msg = ('Value of volumetric flow must be ' + str(v2_design) + ', is ' +
               str(self.sc2_v2.v.val_SI) + '.')
        assert round(v2_design, 5) == round(self.sc2_v2.v.val_SI, 5), msg

        # zeta value of solar collector comparison
        msg = ('Value of zeta must be ' + str(zeta_sc1_design) + ', is ' +
               str(self.sc1.zeta.val) + '.')
        assert round(zeta_sc1_design, 0) == round(self.sc1.zeta.val, 0), msg

        msg = ('Value of zeta must be ' + str(zeta_sc2_design) + ', is ' +
               str(self.sc2.zeta.val) + '.')
        assert round(zeta_sc2_design, 0) == round(self.sc2.zeta.val, 0), msg

        shutil.rmtree('./design1', ignore_errors=True)
        shutil.rmtree('./design2', ignore_errors=True)

    def test_local_offdesign_on_connections_and_components(self):
        """Test local offdesign feature."""
        self.setup_network_individual_offdesign()
        self.nw.solve('design')
        self.sc2_v2.set_attr(m=0)
        self.nw.solve('design')
        self.nw.save('design1')

        self.sc1_v1.set_attr(design=['T'], offdesign=['v'], state='l')
        self.sc2_v2.set_attr(design=['T'], offdesign=['v'], state='l')

        self.sc1.set_attr(local_offdesign=True, design_path='design1')
        self.pump1.set_attr(local_offdesign=True, design_path='design1')
        self.sp_p1.set_attr(local_offdesign=True, design_path='design1')
        self.p1_sc1.set_attr(local_offdesign=True, design_path='design1')
        self.sc1_v1.set_attr(local_offdesign=True, design_path='design1')
        self.sc1.set_attr(E=500)

        self.sc2_v2.set_attr(T=95, m=np.nan)
        self.nw.solve('design')
        self.nw.save('design2')

        # connections and components on side 1 must have switched to offdesign

        msg = ('Solar collector outlet temperature must be different from ' +
               'design value ' + str(round(self.sc1_v1.T.design - 273.15, 1)) +
               ', is ' + str(round(self.sc1_v1.T.val, 1)) + '.')
        assert self.sc1_v1.T.design > self.sc1_v1.T.val, msg

        msg = ('Parameter eta_s_char must be set for pump one.')
        assert self.pump1.eta_s_char.is_set is True, msg

        msg = ('Parameter v must be set for connection from solar collector1 '
               'to pump1.')
        assert self.sc1_v1.v.val_set is True, msg

        shutil.rmtree('./design1', ignore_errors=True)
        shutil.rmtree('./design2', ignore_errors=True)

    def test_missing_design_path_local_offdesign_on_connections(self):
        """Test missing design path on connections in local offdesign mode."""
        self.setup_network_individual_offdesign()
        self.nw.solve('design')
        self.sc2_v2.set_attr(m=0)
        self.nw.solve('design')
        self.nw.save('design1')

        self.sc1_v1.set_attr(design=['T'], offdesign=['v'], state='l')
        self.sc2_v2.set_attr(design=['T'], offdesign=['v'], state='l')

        self.sc1.set_attr(local_offdesign=True, design_path='design1')
        self.pump1.set_attr(local_offdesign=True, design_path='design1')
        self.sp_p1.set_attr(local_offdesign=True, design_path='design1')
        self.p1_sc1.set_attr(local_offdesign=True, design_path='design1')
        self.sc1_v1.set_attr(local_offdesign=True)
        self.sc1.set_attr(E=500)

        self.sc2_v2.set_attr(T=95, m=np.nan)
        try:
            self.nw.solve('design', init_only=True)
        except TESPyNetworkError:
            pass

        shutil.rmtree('./design1', ignore_errors=True)
