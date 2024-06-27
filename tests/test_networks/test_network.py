# -*- coding: utf-8

"""Module for testing network properties.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_network.py

SPDX-License-Identifier: MIT
"""

import os

from pytest import mark
from pytest import raises

from tespy.components import Compressor
from tespy.components import Merge
from tespy.components import Pipe
from tespy.components import Pump
from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import SolarCollector
from tespy.components import Source
from tespy.components import Splitter
from tespy.components import SubsystemInterface
from tespy.components import Turbine
from tespy.components import Valve
from tespy.connections import Connection
from tespy.connections import Ref
from tespy.networks import Network
from tespy.networks import load_network
from tespy.tools.helpers import TESPyNetworkError


class TestNetworks:
    def setup_method(self):
        self.nw = Network(p_unit='bar', v_unit='m3 / s', T_unit='C')
        self.source = Source('source')
        self.sink = Sink('sink')

    def offdesign_TESPyNetworkError(self, **kwargs):
        with raises(TESPyNetworkError):
            self.nw.solve('offdesign', **kwargs)

    def test_Network_linear_dependency(self):
        """Test network linear dependency."""
        a = Connection(
            self.source, 'out1', self.sink, 'in1', p=5, x=1, T=7, fluid={"H2": 1}
        )
        self.nw.add_conns(a)
        self.nw.solve('design')
        msg = 'This test must result in a linear dependency of the jacobian matrix.'
        assert self.nw.lin_dep, msg

    def test_Network_no_progress(self):
        """Test no convergence progress."""
        pi = Pipe('pipe', pr=1, Q=-100e3)
        a = Connection(
            self.source, 'out1', pi, 'in1', m=1, p=1, T=7, fluid={'water': 1}
        )
        b = Connection(pi, 'out1', self.sink, 'in1')
        self.nw.add_conns(a, b)
        self.nw.solve('design')
        msg = ('This test must result in a calculation making no progress, as '
               'the pipe\'s outlet enthalpy is below fluid property range.')
        assert not self.nw.progress, msg

    def test_Network_max_iter(self):
        """Test reaching maximum iteration count."""
        pi = Pipe('pipe', pr=1, Q=100e3)
        a = Connection(
            self.source, 'out1', pi, 'in1', m=1, p=1, T=7, fluid={'water': 1}
        )
        b = Connection(pi, 'out1', self.sink, 'in1')
        self.nw.add_conns(a, b)
        self.nw.solve('design', max_iter=2)
        msg = ('This test must result in the itercount being equal to the max '
               'iter statement.')
        assert self.nw.max_iter == self.nw.iter + 1, msg

    def test_Network_delete_conns(self):
        """Test deleting a network's connection."""
        a = Connection(self.source, 'out1', self.sink, 'in1')
        self.nw.add_conns(a)
        self.nw.check_network()
        msg = ('After the network check, the .checked-property must be True.')
        assert self.nw.checked, msg

        self.nw.del_conns(a)
        msg = ('A connection has been deleted, the network consistency check '
               'must be repeated (.checked-property must be False).')
        assert not self.nw.checked, msg

    def test_Network_delete_comps(self):
        """Test deleting components by deleting connections."""
        p = Pipe("Dummy")
        a = Connection(self.source, 'out1', self.sink, 'in1')
        self.nw.add_conns(a)
        a.set_attr(fluid={"water": 1}, m=1, p=1, T=25)
        self.nw.solve("design")
        msg = ('After the network check, the .checked-property must be True.')
        assert self.nw.checked, msg

        # if a connection is deleted, the respective components must be deleted
        # as well, if there is no other connection in the system containing
        # the respective component
        self.nw.del_conns(a)
        new_sink = Sink("A different sink")
        a = Connection(self.source, 'out1', p, 'in1')
        b = Connection(p, 'out1', new_sink, 'in1')
        self.nw.add_conns(a, b)
        a.set_attr(fluid={"water": 1}, m=1, p=1, T=25)
        b.set_attr(p=1, T=25)
        self.nw.solve("design")
        self.nw._convergence_check()

    def test_Network_missing_connection_in_init_path(self, tmp_path):
        """Test debug message for missing connection in init_path."""
        IF = SubsystemInterface('IF')
        a = Connection(self.source, 'out1', self.sink, 'in1')
        a.set_attr(fluid={"Air": 1})
        self.nw.add_conns(a)
        self.nw.solve('design', init_only=True)
        self.nw.save(tmp_path)
        msg = ('After the network check, the .checked-property must be True.')
        assert self.nw.checked, msg

        self.nw.del_conns(a)
        a = Connection(self.source, 'out1', IF, 'in1')
        b = Connection(IF, 'out1', self.sink, 'in1')
        a.set_attr(fluid={"Air": 1})
        self.nw.add_conns(a, b)
        self.nw.solve('design', init_path=tmp_path, init_only=True)
        msg = ('After the network check, the .checked-property must be True.')
        assert self.nw.checked, msg

    def test_Network_reader_checked(self, tmp_path):
        """Test state of network if loaded successfully from export."""
        a = Connection(self.source, 'out1', self.sink, 'in1')
        self.nw.add_conns(a)
        a.set_attr(fluid={"H2O": 1})
        self.nw.solve('design', init_only=True)
        self.nw.export(tmp_path)
        imported_nwk = load_network(tmp_path)
        imported_nwk.solve('design', init_only=True)
        msg = ('If the network import was successful the network check '
               'should have been successful, too, but it is not.')
        assert imported_nwk.checked, msg

    def test_Network_missing_data_in_design_case_files(self, tmp_path_factory):
        """Test for missing data in design case files."""
        tmp_path = tmp_path_factory.mktemp("tmp")
        tmp_path2 = tmp_path_factory.mktemp("tmp2")
        pi = Pipe('pipe', Q=0, pr=0.95, design=['pr'], offdesign=['zeta'])
        a = Connection(self.source, 'out1', pi, 'in1')
        a.set_attr(m=1, p=1, T=20, fluid={'water': 1})
        b = Connection(pi, 'out1', self.sink, 'in1')
        self.nw.add_conns(a, b)
        self.nw.solve('design')
        self.nw.save(tmp_path)
        self.nw.save(tmp_path2)

        path = os.path.join(tmp_path, "connections.csv")
        inputs = open(path)
        all_lines = inputs.readlines()
        all_lines.pop(len(all_lines) - 1)
        inputs.close()

        with open(os.path.join(tmp_path2, 'connections.csv'), 'w') as out:
            for line in all_lines:
                out.write(line.strip() + '\n')

        self.offdesign_TESPyNetworkError(design_path=tmp_path2, init_only=True)

    def test_Network_missing_data_in_individual_design_case_file(self, tmp_path_factory):
        """Test for missing data in individual design case files."""
        tmp_path = tmp_path_factory.mktemp("tmp")
        tmp_path2 = tmp_path_factory.mktemp("tmp2")
        pi = Pipe('pipe', Q=0, pr=0.95, design=['pr'], offdesign=['zeta'])
        a = Connection(self.source, 'out1', pi, 'in1')
        a.set_attr(m=1, p=1, T=293.15, fluid={'water': 1})
        b = Connection(pi, 'out1', self.sink, 'in1')
        b.set_attr(design_path=tmp_path2)
        self.nw.add_conns(a, b)
        self.nw.solve('design')
        self.nw.save(tmp_path)
        self.nw.save(tmp_path2)

        path = os.path.join(tmp_path, "connections.csv")
        inputs = open(path)
        all_lines = inputs.readlines()
        all_lines.pop(len(all_lines) - 1)
        inputs.close()

        with open(os.path.join(tmp_path2, 'connections.csv'), 'w') as out:
            for line in all_lines:
                out.write(line.strip() + '\n')

        self.offdesign_TESPyNetworkError(design_path=tmp_path, init_only=True)

    def test_Network_missing_connection_in_design_path(self, tmp_path):
        """Test for missing connection data in design case files."""
        pi = Pipe('pipe', Q=0, pr=0.95, design=['pr'], offdesign=['zeta'])
        a = Connection(
            self.source, 'out1', pi, 'in1', m=1, p=1, T=293.15,
            fluid={'water': 1}
        )
        b = Connection(pi, 'out1', self.sink, 'in1')
        self.nw.add_conns(a, b)
        self.nw.solve('design')
        self.nw.save(tmp_path)

        path = os.path.join(tmp_path, "connections.csv")
        inputs = open(path)
        all_lines = inputs.readlines()
        all_lines.pop(len(all_lines) - 1)
        inputs.close()

        with open(path, 'w') as out:
            for line in all_lines:
                out.write(line.strip() + '\n')

        self.offdesign_TESPyNetworkError(design_path=tmp_path)

    def test_Network_get_comp_without_connections_added(self):
        """Test if components are found prior to initialization."""
        pi = Pipe('pipe')
        a = Connection(self.source, 'out1', pi, 'in1')
        Connection(pi, 'out1', self.sink, 'in1')
        self.nw.add_conns(a)
        msg = (
            "A component with the label 'sink' has been created but must not "
            "be part of the network as the respective connection has not "
            "been added."
        )
        assert self.nw.get_comp("sink") is None, msg

    def test_Network_get_comp_before_initialization(self):
        """Test if components are found prior to initialization."""
        pi = Pipe('pipe')
        a = Connection(self.source, 'out1', pi, 'in1')
        b = Connection(pi, 'out1', self.sink, 'in1')
        self.nw.add_conns(a, b)
        msg = (
            "A component with the label 'pipe' is part of the network "
            "and therefore must be found in the DataFrame."
        )
        assert self.nw.get_comp("pipe") == pi, msg


class TestNetworkIndividualOffdesign:

    def setup_Network_individual_offdesign(self):
        """Set up network for individual offdesign tests."""
        self.nw = Network(T_unit='C', p_unit='bar', v_unit='m3 / s')

        so = Source('source')
        sp = Splitter('splitter', num_out=2)
        self.pump1 = Pump('pump 1')
        self.sc1 = SolarCollector('collector field 1')
        v1 = Valve('valve1')
        self.pump2 = Pump('pump 2')
        self.sc2 = SolarCollector('collector field 2')
        v2 = Valve('valve2')
        me = Merge('merge', num_in=2)
        si = Sink('sink')

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
        inlet = Connection(so, 'out1', sp, 'in1', T=50, p=3, fluid=fl)
        outlet = Connection(me, 'out1', si, 'in1', p=3)

        self.sp_p1 = Connection(sp, 'out1', self.pump1, 'in1')
        self.p1_sc1 = Connection(self.pump1, 'out1', self.sc1, 'in1')
        self.sc1_v1 = Connection(self.sc1, 'out1', v1, 'in1', p=3.1, T=90)
        v1_me = Connection(v1, 'out1', me, 'in1')

        self.sp_p2 = Connection(sp, 'out2', self.pump2, 'in1')
        self.p2_sc2 = Connection(self.pump2, 'out1', self.sc2, 'in1')
        self.sc2_v2 = Connection(self.sc2, 'out1', v2, 'in1', p=3.1, m=0.1)
        v2_me = Connection(v2, 'out1', me, 'in2')

        self.nw.add_conns(inlet, outlet, self.sp_p1, self.p1_sc1, self.sc1_v1,
                          v1_me, self.sp_p2, self.p2_sc2, self.sc2_v2, v2_me)

    def test_individual_design_path_on_connections_and_components(self, tmp_path_factory):
        """Test individual design path specification."""
        tmp_path1 = tmp_path_factory.mktemp("tmp1")
        tmp_path2 = tmp_path_factory.mktemp("tmp2")
        self.setup_Network_individual_offdesign()
        self.nw.solve('design')
        self.nw._convergence_check()
        self.sc2_v2.set_attr(m=0)
        self.nw.solve('design')
        self.nw._convergence_check()
        self.nw.save(tmp_path1)
        v1_design = self.sc1_v1.v.val_SI
        zeta_sc1_design = self.sc1.zeta.val

        self.sc2_v2.set_attr(T=95, state='l', m=None)
        self.sc1_v1.set_attr(m=0.001, T=None)
        self.nw.solve('design')
        self.nw._convergence_check()
        self.nw.save(tmp_path2)
        v2_design = self.sc2_v2.v.val_SI
        zeta_sc2_design = self.sc2.zeta.val

        self.sc1_v1.set_attr(m=None)
        self.sc1_v1.set_attr(design=['T'], offdesign=['v'], state='l')
        self.sc2_v2.set_attr(design=['T'], offdesign=['v'], state='l')

        self.sc2.set_attr(design_path=tmp_path2)
        self.pump2.set_attr(design_path=tmp_path2)
        self.sp_p2.set_attr(design_path=tmp_path2)
        self.p2_sc2.set_attr(design_path=tmp_path2)
        self.sc2_v2.set_attr(design_path=tmp_path2)
        self.nw.solve('offdesign', design_path=tmp_path1)
        self.nw._convergence_check()

        self.sc1.set_attr(E=500)
        self.sc2.set_attr(E=950)

        self.nw.solve('offdesign', design_path=tmp_path1)
        self.nw._convergence_check()
        self.sc2_v2.set_attr(design_path=None)

        # volumetric flow comparison
        msg = f"Design path was set to None, is {self.sc2_v2.design_path}."
        assert self.sc2_v2.design_path is None, msg

        # volumetric flow comparison
        msg = (
            f"Value of volumetric flow must be {v1_design}, is "
            f"{self.sc1_v1.v.val_SI}."
        )
        assert round(v1_design, 5) == round(self.sc1_v1.v.val_SI, 5), msg

        msg = (
            f"Value of volumetric flow must be {v2_design}, is "
            f"{self.sc2_v2.v.val_SI}."
        )
        assert round(v2_design, 5) == round(self.sc2_v2.v.val_SI, 5), msg

        # zeta value of solar collector comparison
        msg = (
            f"Value of zeta must be {zeta_sc1_design}, is {self.sc1.zeta.val}."
        )
        assert round(zeta_sc1_design, 0) == round(self.sc1.zeta.val, 0), msg

        msg = (
            f"Value of zeta must be {zeta_sc2_design}, is {self.sc2.zeta.val}."
        )
        assert round(zeta_sc2_design, 0) == round(self.sc2.zeta.val, 0), msg

    def test_local_offdesign_on_connections_and_components(self, tmp_path_factory):
        """Test local offdesign feature."""
        tmp_path1 = tmp_path_factory.mktemp("tmp1")
        tmp_path2 = tmp_path_factory.mktemp("tmp2")
        self.setup_Network_individual_offdesign()
        self.nw.solve('design')
        self.nw._convergence_check()
        self.sc2_v2.set_attr(m=0)
        self.nw.solve('design')
        self.nw._convergence_check()
        self.nw.save(tmp_path1)

        self.sc1_v1.set_attr(design=['T'], offdesign=['v'], state='l')
        self.sc2_v2.set_attr(design=['T'], offdesign=['v'], state='l')

        self.sc1.set_attr(local_offdesign=True, design_path=tmp_path1)
        self.pump1.set_attr(local_offdesign=True, design_path=tmp_path1)
        self.sp_p1.set_attr(local_offdesign=True, design_path=tmp_path1)
        self.p1_sc1.set_attr(local_offdesign=True, design_path=tmp_path1)
        self.sc1_v1.set_attr(local_offdesign=True, design_path=tmp_path1)
        self.sc1.set_attr(E=500)

        self.sc2_v2.set_attr(T=95, m=None)
        self.nw.solve('design')
        self.nw._convergence_check()
        self.nw.save(tmp_path2)

        # connections and components on side 1 must have switched to offdesign

        msg = ('Solar collector outlet temperature must be different from ' +
               'design value ' + str(round(self.sc1_v1.T.design - 273.15, 1)) +
               ', is ' + str(round(self.sc1_v1.T.val, 1)) + '.')
        assert self.sc1_v1.T.design > self.sc1_v1.T.val, msg

        msg = "Parameter eta_s_char must be set for pump one."
        assert self.pump1.eta_s_char.is_set, msg

        msg = (
            "Parameter v must be set for connection from solar collector1 to "
            "pump1."
        )
        assert self.sc1_v1.v.is_set, msg

    def test_missing_design_path_local_offdesign_on_connections(self, tmp_path):
        """Test missing design path on connections in local offdesign mode."""
        self.setup_Network_individual_offdesign()
        self.nw.solve('design')
        self.nw._convergence_check()
        self.sc2_v2.set_attr(m=0)
        self.nw.solve('design')
        self.nw._convergence_check()
        self.nw.save(tmp_path)

        self.sc1_v1.set_attr(design=['T'], offdesign=['v'], state='l')
        self.sc2_v2.set_attr(design=['T'], offdesign=['v'], state='l')

        self.sc1.set_attr(local_offdesign=True, design_path=tmp_path)
        self.pump1.set_attr(local_offdesign=True, design_path=tmp_path)
        self.sp_p1.set_attr(local_offdesign=True, design_path=tmp_path)
        self.p1_sc1.set_attr(local_offdesign=True, design_path=tmp_path)
        self.sc1_v1.set_attr(local_offdesign=True)
        self.sc1.set_attr(E=500)

        self.sc2_v2.set_attr(T=95, m=None)
        try:
            self.nw.solve('design', init_only=True)
        except TESPyNetworkError:
            pass

class TestNetworkPreprocessing:

    def setup_method(self):
        self.nwk = Network(T_unit="C", p_unit="bar", h_unit='kJ / kg')

    def _create_linear_branch(self):
        a = Source("source")
        b = Pipe("pipe")
        c = Sink("sink")

        c1 = Connection(a, "out1", b, "in1", label="1")
        c2 = Connection(b, "out1", c, "in1", label="2")

        self.nwk.add_conns(c1, c2)

    def _create_recirculation(self):
        self.nwk = Network()

        self.nwk.set_attr(T_unit='C', p_unit='bar', h_unit='kJ / kg')

        source = Source('source')
        merge = Merge('merge')
        component1 = SimpleHeatExchanger('comp1', pr=1)
        splitter = Splitter('splitter')
        component2 = SimpleHeatExchanger('comp2')
        sink = Sink('sink')

        c1 = Connection(source, 'out1', merge, 'in1', label="1")
        c2 = Connection(merge, 'out1', component1, 'in1', label="2")
        c3 = Connection(component1, 'out1', splitter, 'in1', label="3")
        c4 = Connection(splitter, 'out1', component2, 'in1', label="4")
        c5 = Connection(component2, 'out1', merge, 'in2', label="5")
        c6 = Connection(splitter, 'out2', sink, 'in1', label="6")

        self.nwk.add_conns(c1, c2, c3, c4, c5, c6)

        c1.set_attr(p=1, h=200, m=10)
        c3.set_attr(h=180)
        c4.set_attr(m=1)
        c5.set_attr(h=170)

    @mark.skip("Not implemented")
    def test_fluid_linear_branch_distribution(self):
        raise NotImplementedError()

    @mark.skip("Not implemented")
    def test_fluid_connected_branches_distribution(self):
        raise NotImplementedError()

    @mark.skip("Not implemented")
    def test_fluid_independent_branches_distribution(self):
        raise NotImplementedError()

    def test_linear_branch_massflow_presolve(self):
        self._create_linear_branch()
        c1, c2 = self.nwk.get_conn(["1", "2"])
        b = self.nwk.get_comp("pipe")
        c1.set_attr(fluid={"O2": 1}, p=10, T=-20, m=5)
        c2.set_attr(T=-15)
        b.set_attr(pr=1)
        self.nwk.solve("design")
        self.nwk._convergence_check()
        variables = [data["obj"].get_attr(data["variable"]) for data in self.nwk.variables_dict.values()]
        # no mass flow is variable
        assert c1.m not in variables
        assert c2.m not in variables
        # first connection pressure and enthalpy not variable
        assert c1.p not in variables
        assert c1.h not in variables
        # second connection pressure and enthalpy are variable
        assert c2.p in variables
        assert c2.h in variables

    @mark.skip("Not implemented")
    def test_splitting_branch_massflow_presolve(self):
        raise NotImplementedError()

    def test_double_massflow_specification_linear_branch(self):
        self._create_linear_branch()
        c1, c2 = self.nwk.get_conn(["1", "2"])
        c1.set_attr(fluid={"O2": 1}, m=5)
        c2.set_attr(m=5)
        with raises(TESPyNetworkError):
            self.nwk.solve("design")

    def test_missing_fluid_information(self):
        self._create_linear_branch()
        with raises(TESPyNetworkError):
            self.nwk.solve("design")

    def test_referencing_massflow_specification_linear_branch(self):
        self._create_linear_branch()
        c1, c2 = self.nwk.get_conn(["1", "2"])
        c1.set_attr(fluid={"O2": 1})
        c2.set_attr(m=Ref(c1, 1, 0))
        with raises(TESPyNetworkError):
            self.nwk.solve("design")

    @mark.skip("Not implemented")
    def test_linear_branch_fluid_presolve(self):
        raise NotImplementedError()

    @mark.skip("Not implemented")
    def test_splitting_branch_fluid_presolve(self):
        raise NotImplementedError()

    @mark.skip("Not implemented")
    def test_independent_branch_fluid_presolve(self):
        raise NotImplementedError()

    def test_recirculation_structure_two_fluids_without_starting(self):
        self._create_recirculation()
        c1, c2 = self.nwk.get_conn(["1", "2"])
        c1.set_attr(fluid={"water": 0, "R134a": 1})
        self.nwk.solve("design", init_only=True)
        assert c2.fluid.val["R134a"] == 0.5
        assert c2.fluid.val["water"] == 0.5

    def test_recirculation_structure_two_fluids_with_starting(self):
        self._create_recirculation()
        c1, c2, c6 = self.nwk.get_conn(["1", "2", "6"])
        c1.set_attr(fluid={"water": 0, "R134a": 1})
        c2.set_attr(fluid0={"water": 0, "R134a": 1})
        self.nwk.solve("design", init_only=True)
        assert c6.fluid.val["R134a"] == 1

    def test_recirculation_structure_single_fluids(self):
        self._create_recirculation()
        c1, c2 = self.nwk.get_conn(["1", "2"])
        c1.set_attr(fluid={"R134a": 1})
        self.nwk.solve("design", init_only=True)
        assert c2.fluid.val["R134a"] == 1

def test_use_cuda_without_it_being_installed():
    nw = Network()

    a = Source("turbine")
    b = Sink("Compressor")

    c1 = Connection(a, "out1", b, "in1")

    nw.add_conns(c1)

    c1.set_attr(m=1, p=1e5, T=300, fluid={"INCOMP::Water": 1})
    nw.solve("design", use_cuda=True)
    nw._convergence_check()
    assert not nw.use_cuda

def test_component_not_found():
    nw = Network()

    a = Turbine("turbine")
    b = Compressor("Compressor")

    c1 = Connection(a, "out1", b, "in1")
    c2 = Connection(b, "out1", a, "in1")

    nw.add_conns(c1, c2)
    assert nw.get_comp("Turbine") is None

def test_connection_not_found():
    nw = Network()

    a = Turbine("turbine")
    b = Compressor("Compressor")

    c1 = Connection(a, "out1", b, "in1")
    c2 = Connection(b, "out1", a, "in1")

    nw.add_conns(c1, c2)
    assert nw.get_conn("1") is None

def test_missing_source_sink_cycle_closer():
    nw = Network()

    a = Turbine("turbine")
    b = Compressor("Compressor")

    c1 = Connection(a, "out1", b, "in1")
    c2 = Connection(b, "out1", a, "in1")

    nw.add_conns(c1, c2)
    with raises(TESPyNetworkError):
        nw.check_network()
