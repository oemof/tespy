# -*- coding: utf-8

"""Module for testing Network class.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_network.py

SPDX-License-Identifier: MIT
"""
import json
import os
import warnings
from copy import deepcopy

import numpy as np
from pytest import approx
from pytest import mark
from pytest import raises

from tespy.components import Compressor
from tespy.components import Merge
from tespy.components import MovingBoundaryHeatExchanger
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
from tespy.components import WaterElectrolyzer
from tespy.connections import Connection
from tespy.connections import Ref
from tespy.networks import Network
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.fluid_properties import conductivity_mix_ph
from tespy.tools.fluid_properties.wrappers import IncompressibleFluidWrapper
from tespy.tools.helpers import TESPyNetworkError
from tespy.tools.helpers import UserDefinedEquation
from tespy.tools.helpers import _numeric_deriv


class TestNetworks:
    def setup_method(self):
        self.nw = Network()
        self.nw.units.set_defaults(**{
            "pressure": "bar", "pressure_difference": "bar",
            "temperature": "degC"
        })
        self.source = Source('source')
        self.sink = Sink('sink')

    def offdesign_TESPyNetworkError(self, **kwargs):
        with raises(TESPyNetworkError):
            self.nw.solve('offdesign', **kwargs)

    def test_Network_overdetermined_fluid_state(self):
        """Test network linear dependency."""
        a = Connection(
            self.source, 'out1', self.sink, 'in1', p=5, x=1, T=7, fluid={"H2": 1}
        )
        self.nw.add_conns(a)
        with raises(TESPyNetworkError):
            self.nw.solve('design')

    def test_Network_linear_dependency(self):
        pump = Pump("pump")
        a = Connection(self.source, "out1", pump, "in1")
        b = Connection(pump, "out1", self.sink, "in1")
        self.nw.add_conns(a, b)
        a.set_attr(fluid={"water": 1}, p=1, T=25)
        b.set_attr(p=2, T=27)
        pump.set_attr(eta_s=0.9)

        self.nw.solve("design")
        assert self.nw.status == 3

    def test_Network_no_progress(self):
        """Test no convergence progress."""
        pi = Pipe('pipe', pr=1, Q=-100e3)
        a = Connection(
            self.source, 'out1', pi, 'in1', m=1, p=1, T=7, fluid={'water': 1}
        )
        b = Connection(pi, 'out1', self.sink, 'in1')
        self.nw.add_conns(a, b)
        self.nw.solve('design')
        msg = (
            'Test must result in a calculation making no progress, as the '
            'pipe\'s outlet enthalpy is below fluid property range.'
        )
        assert self.nw.status == 2, msg

    def test_Network_max_iter(self):
        """Test reaching maximum iteration count."""
        pi = Pipe('pipe', pr=1, Q=100e3)
        a = Connection(
            self.source, 'out1', pi, 'in1', m=1, p=1, T=7, fluid={'water': 1}
        )
        b = Connection(pi, 'out1', self.sink, 'in1')
        self.nw.add_conns(a, b)
        self.nw.solve('design', max_iter=2)
        assert self.nw.status == 2
        msg = 'Test must result in the itercount being equal to max_iter.'
        assert self.nw.max_iter == self.nw.iter + 1, msg

    def test_Network_delete_conns(self):
        """Test deleting a network's connection."""
        a = Connection(self.source, 'out1', self.sink, 'in1')
        self.nw.add_conns(a)
        self.nw.check_topology()
        msg = ('After the network check, the .checked-property must be True.')
        assert self.nw.checked, msg

        self.nw.del_conns(a)
        msg = (
            'A connection has been deleted, the network consistency check '
            'must be repeated (.checked-property must be False).'
        )
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
        self.nw.assert_convergence()

    def test_Network_missing_connection_in_init_path(self):
        """Test debug message for missing connection in init_path."""
        IF = SubsystemInterface('IF')
        a = Connection(self.source, 'out1', self.sink, 'in1')
        a.set_attr(fluid={"Air": 1})
        self.nw.add_conns(a)
        self.nw.solve('design', init_only=True)
        design_state = self.nw.save(as_dict=True)
        msg = ('After the network check, the .checked-property must be True.')
        assert self.nw.checked, msg

        self.nw.del_conns(a)
        a = Connection(self.source, 'out1', IF, 'in1')
        b = Connection(IF, 'out1', self.sink, 'in1')
        a.set_attr(fluid={"Air": 1})
        self.nw.add_conns(a, b)
        self.nw.solve('design', init_path=design_state, init_only=True)
        msg = ('After the network check, the .checked-property must be True.')
        assert self.nw.checked, msg

    def test_Network_from_json_checked(self, tmp_path):
        """Test state of network if loaded successfully from export."""
        tmp_path = f"{tmp_path}.json"
        a = Connection(self.source, 'out1', self.sink, 'in1')
        self.nw.add_conns(a)
        a.set_attr(fluid={"H2O": 1})
        self.nw.solve('design', init_only=True)
        self.nw.export(tmp_path)
        imported_nwk = Network.from_json(tmp_path)
        imported_nwk.solve('design', init_only=True)
        msg = (
            'If the network import was successful the network check '
            'should have been successful, too, but it is not.'
        )
        assert imported_nwk.checked, msg

    def test_Network_from_dict(self):
        """Test state of network if loaded successfully from export."""
        a = Connection(self.source, 'out1', self.sink, 'in1')
        self.nw.add_conns(a)
        a.set_attr(fluid={"H2O": 1})
        self.nw.solve('design', init_only=True)
        serialization = self.nw.export()
        imported_nwk = Network.from_dict(serialization)
        imported_nwk.solve('design', init_only=True)
        msg = (
            'If the network import was successful the network check '
            'should have been successful, too, but it is not.'
        )
        assert imported_nwk.checked, msg

    def test_Network_import_with_component_parameter_as_variable(self):
        """Test if component variables are retained after import."""
        pipe = Pipe("pipe")
        c1 = Connection(self.source, 'out1', pipe, 'in1', label="c1")
        c2 = Connection(pipe, 'out1', self.sink, 'in1', label="c2")
        self.nw.add_conns(c1, c2)
        c1.set_attr(fluid={"H2O": 1}, m=1, T=25, p=2)
        c2.set_attr(p=1.9)
        pipe.set_attr(Q=0, D="var", ks=0.00005, L=100)
        self.nw.solve("design")
        self.nw.assert_convergence()
        serialization = self.nw.export()
        imported_nwk = Network.from_dict(serialization)
        imported_nwk.solve("design")
        imported_nwk.assert_convergence()
        assert approx(pipe.D.val) == imported_nwk.get_comp("pipe").D.val

    def test_Network_deserialze_component_with_default_charmap(self):
        """Test if component variables are retained after import."""
        pump = Pump("pump")
        c1 = Connection(self.source, 'out1', pump, 'in1', label="c1")
        c2 = Connection(pump, 'out1', self.sink, 'in1', label="c2")
        self.nw.add_conns(c1, c2)
        c1.set_attr(fluid={"H2O": 1}, m=1, T=25, p=2)
        c2.set_attr(p=3)
        pump.set_attr(eta=0.7)
        self.nw.solve("design")
        self.nw.assert_convergence()
        serialization = self.nw.export()
        # this deserialization fails if the CharMap is not correctly loaded
        imported_nwk = Network.from_dict(serialization)
        len(imported_nwk.get_comp("pump").head_flow_map.char_func.x.flatten()) == 2
        len(imported_nwk.get_comp("pump").head_flow_map.char_func.y.flatten()) == 4
        len(imported_nwk.get_comp("pump").head_flow_map.char_func.z.flatten()) == 4


    def test_Network_reader_unknown_component_class(self, tmp_path):
        """Test notsupported component."""
        tmp_path = f"{tmp_path}.json"
        a = Connection(self.source, 'out1', self.sink, 'in1')
        self.nw.add_conns(a)
        a.set_attr(fluid={"H2O": 1})
        self.nw.solve('design', init_only=True)
        self.nw.export(tmp_path)

        with open(tmp_path, "r") as f:
            data = json.load(f)

        data["Component"]["TestComponent"] = {}

        with open(tmp_path, "w") as f:
            json.dump(data, f)

        with raises(TESPyNetworkError):
            Network.from_json(tmp_path)


    def test_Network_missing_data_in_individual_design_case_file(self):
        """Test for missing data in individual design case files."""
        pi = Pipe('pipe', Q=0, pr=0.95, design=['pr'], offdesign=['zeta'])
        a = Connection(self.source, 'out1', pi, 'in1')
        a.set_attr(m=1, p=1, T=293.15, fluid={'water': 1})
        b = Connection(pi, 'out1', self.sink, 'in1')
        self.nw.add_conns(a, b)
        self.nw.solve('design')
        design_state1 = self.nw.save(as_dict=True)
        design_state2 = deepcopy(design_state1)
        design_state2["Connection"] = {}

        b.set_attr(design_path=design_state2)
        self.offdesign_TESPyNetworkError(design_path=design_state1, init_only=True)

    def test_Network_missing_connection_in_design_path(self):
        """Test for missing connection data in design case files."""
        pi = Pipe('pipe', Q=0, pr=0.95, design=['pr'], offdesign=['zeta'])
        a = Connection(
            self.source, 'out1', pi, 'in1', m=1, p=1, T=293.15,
            fluid={'water': 1}
        )
        b = Connection(pi, 'out1', self.sink, 'in1')
        self.nw.add_conns(a, b)
        self.nw.solve('design')
        data = self.nw.save(as_dict=True)

        data["Connection"] = {}

        self.offdesign_TESPyNetworkError(design_path=data)

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

    def test_Network_access_converged_before_solve(self):
        with raises(AttributeError):
            self.nw.converged


class TestNetworkPreprocessing:

    def setup_method(self):
        self.nwk = Network()
        self.nwk.units.set_defaults(**{
            "pressure": "bar", "pressure_difference": "bar",
            "temperature": "degC", "enthalpy": "kJ/kg"
        })

    def _create_linear_branch(self):
        a = Source("source")
        b = Pipe("pipe")
        c = Sink("sink")

        c1 = Connection(a, "out1", b, "in1", label="1")
        c2 = Connection(b, "out1", c, "in1", label="2")

        self.nwk.add_conns(c1, c2)

    def _create_recirculation(self):
        self.nwk = Network()

        self.nwk = Network()
        self.nwk.units.set_defaults(**{
            "pressure": "bar", "pressure_difference": "bar",
            "temperature": "degC", "enthalpy": "kJ/kg"
        })

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
        self.nwk.assert_convergence()
        variables = [
            data["obj"].get_attr(data["variable"])
            for data in self.nwk.variables_dict.values()
        ]
        # no variable at all, everything must have been presolved
        assert c1.m not in variables
        assert c2.m not in variables
        # first connection pressure and enthalpy not variable
        assert c1.p not in variables
        assert c1.h not in variables
        # second connection pressure and enthalpy not variable
        assert c2.p not in variables
        assert c2.h not in variables

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


def test_temperature_unit_C_raises():
    nw = Network()
    with raises(ValueError):
        nw.units.set_defaults(temperature="C")


def test_use_cuda_without_it_being_installed():
    nw = Network()

    a = Source("turbine")
    b = Sink("Compressor")

    c1 = Connection(a, "out1", b, "in1")

    nw.add_conns(c1)

    c1.set_attr(m=1, p=1e5, T=300, fluid={"INCOMP::Water": 1})
    nw.solve("design", use_cuda=True)
    nw.assert_convergence()
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


def _make_simple_network():
    nw = Network()
    so = Source("source")
    si = Sink("sink")
    c = Connection(so, "out1", si, "in1", label="c1")
    nw.add_conns(c)
    return nw, so, si, c


def test_get_comp_found():
    nw, so, si, c = _make_simple_network()
    assert nw.get_comp("source") is so


def test_get_comp_not_found_warns():
    nw, *_ = _make_simple_network()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        result = nw.get_comp("nonexistent")
    assert result is None
    assert len(w) == 1
    assert issubclass(w[0].category, FutureWarning)


def test_get_conn_found():
    nw, so, si, c = _make_simple_network()
    assert nw.get_conn("c1") is c


def test_get_conn_not_found_warns():
    nw, *_ = _make_simple_network()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        result = nw.get_conn("nonexistent")
    assert result is None
    assert len(w) == 1
    assert issubclass(w[0].category, FutureWarning)


def test_get_ude_found():
    nw, so, si, c = _make_simple_network()
    ude = UserDefinedEquation("my_ude", lambda u: 0, lambda u: [], conns=[c])
    nw.add_ude(ude)
    assert nw.get_ude("my_ude") is ude


def test_get_ude_not_found_warns():
    nw, *_ = _make_simple_network()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        result = nw.get_ude("nonexistent")
    assert result is None
    assert len(w) == 1
    assert issubclass(w[0].category, FutureWarning)


def test_missing_source_sink_cycle_closer():
    nw = Network()

    a = Turbine("turbine")
    b = Compressor("Compressor")

    c1 = Connection(a, "out1", b, "in1")
    c2 = Connection(b, "out1", a, "in1")

    nw.add_conns(c1, c2)
    with raises(TESPyNetworkError):
        nw.solve("design")


def test_dublicated_linear_dependent_variables():
    nw = Network()
    nw.units.set_defaults(**{
        "pressure": "bar", "pressure_difference": "bar", "temperature": "degC"
    })

    so = Source("source")
    heater = SimpleHeatExchanger("heater")
    compressor = Compressor("compressor")
    turbine = Turbine("turbine")
    si = Sink("sink")

    c1 = Connection(so, "out1", compressor, "in1", label="c1")
    c2 = Connection(compressor, "out1", heater, "in1", label="c2")
    c3 = Connection(heater, "out1", turbine, "in1", label="c3")
    c4 = Connection(turbine, "out1", si, "in1", label="c4")

    nw.add_conns(c1, c2, c3, c4)

    # fluid has to be specified, otherwise crash due to other issue
    c1.set_attr(fluid={"air": 1}, p=Ref(c4, 1, 0))
    c4.set_attr(p=Ref(c1, 1, 0))

    with raises(TESPyNetworkError):
        nw.solve("design", init_only=True)


def test_cyclic_linear_dependent_variables():
    nw = Network()
    nw.units.set_defaults(**{
        "pressure": "bar", "pressure_difference": "bar", "temperature": "degC"
    })

    so = Source("source")
    heater = SimpleHeatExchanger("heater")
    compressor = Compressor("compressor")
    turbine = Turbine("turbine")
    si = Sink("sink")

    c1 = Connection(so, "out1", compressor, "in1", label="c1")
    c2 = Connection(compressor, "out1", heater, "in1", label="c2")
    c3 = Connection(heater, "out1", turbine, "in1", label="c3")
    c4 = Connection(turbine, "out1", si, "in1", label="c4")

    nw.add_conns(c1, c2, c3, c4)

    # fluid has to be specified, otherwise crash due to other issue
    c1.set_attr(fluid={"air": 1}, p=Ref(c2, 1, 0))
    c2.set_attr(p=Ref(c4, 1, 0))
    c4.set_attr(p=Ref(c1, 1, 0))

    with raises(TESPyNetworkError):
        nw.solve("design", init_only=True)

    adjacency_list, _, _, _ = (
        nw._build_graph(nw._structure_matrix, nw._rhs)
    )
    # Detect cycles (to check for circular dependencies)
    cycle = nw._find_cycles_in_graph(
        {k: [x[0] for x in v] for k, v in adjacency_list.items()}
    )
    # checksum for the variable numbers
    assert sum(cycle) == 19


def test_cyclic_linear_dependent_with_merge_and_split():
    nw = Network()
    nw.units.set_defaults(**{
        "pressure": "bar", "pressure_difference": "bar", "temperature": "degC"
    })

    so = Source("source")
    splitter = Splitter("splitter")
    heater1 = SimpleHeatExchanger("heater 1")
    heater2 = SimpleHeatExchanger("heater 2")
    merge = Merge("merge")
    si = Sink("sink")

    c1 = Connection(so, "out1", splitter, "in1", label="c1")
    c2 = Connection(splitter, "out1", heater1, "in1", label="c2")
    c3 = Connection(heater1, "out1", merge, "in1", label="c3")
    c4 = Connection(splitter, "out2", heater2, "in1", label="c4")
    c5 = Connection(heater2, "out1", merge, "in2", label="c5")
    c6 = Connection(merge, "out1", si, "in1", label="c6")

    nw.add_conns(c1, c2, c3, c4, c5, c6)

    # fluid has to be specified, otherwise crash due to other issue
    c1.set_attr(fluid={"air": 1}, p=1)
    heater1.set_attr(pr=0.98)
    heater2.set_attr(pr=0.98)

    with raises(TESPyNetworkError):
        nw.solve("design", init_only=True)

    adjacency_list, _, _, _ = (
        nw._build_graph(nw._structure_matrix, nw._rhs)
    )
    # Detect cycles (to check for circular dependencies)
    cycle = nw._find_cycles_in_graph(
        {k: [x[0] for x in v] for k, v in adjacency_list.items()}
    )
    # checksum for the variable numbers
    assert sum(cycle) == 45


def test_missing_cyclecloser_but_no_missing_source():
    nw = Network()
    nw.units.set_defaults(**{
        "pressure": "bar", "pressure_difference": "bar", "temperature": "degC",
        "enthalpy": "kJ/kg"
    })

    # Components
    source = Source("source")
    sink = Sink("sink")
    comp = Compressor("compressor")
    cond = MovingBoundaryHeatExchanger("condenser")
    valve = Valve("valve")
    evap = SimpleHeatExchanger("evaporator")

    # Connections (closed loop but no CycleCloser)
    c1 = Connection(evap, "out1", comp, "in1")
    c2 = Connection(comp, "out1", cond, "in1")
    c3 = Connection(cond, "out1", valve, "in1")
    c4 = Connection(valve, "out1", evap, "in1")
    c5 = Connection(source, "out1", cond, "in2")
    c6 = Connection(cond, "out2", sink, "in1")

    nw.add_conns(c1, c2, c3, c4,c5,c6)

    # Set fluid and boundary conditions
    c2.set_attr(fluid={"PROPANE": 1})
    c5.set_attr(fluid={'WATER':1}, p=1, T=50)

    # Component parameters
    comp.set_attr(eta_s=0.7)
    cond.set_attr(td_pinch=3, Q=-15e3)
    evap.set_attr(pr=1, Tamb = 5)

    # This will fail with a fluid key error (instead of warning for the
    # absence of cycle closer)
    with raises(TESPyNetworkError):
        nw.solve(mode="design")


def test_two_phase_in_supercritical_starting_pressure_convergence():
    nw = Network()
    nw.units.set_defaults(**{
        "pressure": "bar", "pressure_difference": "bar", "temperature": "degC"
    })

    so = Source("source")
    heater = SimpleHeatExchanger("heater")
    si = Sink("sink")

    c1 = Connection(so, "out1", heater, "in1", label="c1")
    c2 = Connection(heater, "out1", si, "in1", label="c2")

    nw.add_conns(c1, c2)

    c1.set_attr(fluid={"water": 1}, m=1, p=250, T=400)
    c2.set_attr(x=1, p0=250)

    heater.set_attr(Q=0)

    nw.solve("design")
    nw.assert_convergence()
    assert approx(c2.h.val_SI) == c1.h.val_SI
    assert approx(c2.p.val) == 160.67964


def test_two_phase_in_supercritical_pressure_non_convergence():
    nw = Network()
    nw.units.set_defaults(**{
        "pressure": "bar", "pressure_difference": "bar", "temperature": "degC"
    })

    so = Source("source")
    heater = SimpleHeatExchanger("heater")
    si = Sink("sink")

    c1 = Connection(so, "out1", heater, "in1", label="c1")
    c2 = Connection(heater, "out1", si, "in1", label="c2")

    nw.add_conns(c1, c2)

    c1.set_attr(fluid={"water": 1}, m=1, p=500, T=400)
    c2.set_attr(x=1, p0=250)

    heater.set_attr(Q=0)

    nw.solve("design")
    assert nw.status == 99


def test_postprocessing_supercritical():
    nw = Network()
    nw.units.set_defaults(**{
        "pressure": "bar", "pressure_difference": "bar", "temperature": "degC"
    })

    so = Source("source")
    si = Sink("sink")

    c1 = Connection(so, "out1", si, "in1", label="c1")

    nw.add_conns(c1)

    c1.set_attr(fluid={"water": 1}, m=1, p=500, T=400)
    nw.solve("design")
    assert np.isnan(c1.td_dew.val)
    assert np.isnan(c1.x.val)


def test_nonconverged_simulation_does_not_overwrite_component_specification_1():
    """This creates a result, that shows as converged but actually it did not
    because of internal convergence helpers. It tests, that the user specified
    input is not overwritten by the erroneous result
    """
    nw = Network()
    nw.units.set_defaults(
        pressure="bar", pressure_difference="bar", temperature="degC"
    )

    inflow = Source("inflow")
    outflow = Sink("outflow")
    instance = SimpleHeatExchanger("heat exchanger")

    c1 = Connection(inflow, "out1", instance, "in1")
    c2 = Connection(instance, "out1", outflow, "in1")

    nw.add_conns(c1, c2)
    c1.set_attr(m=0.1, fluid={"N2": 0.7, "O2": 0.15, "Water": 0.15})
    c2.set_attr(p=1, T=20)
    instance.set_attr(Q=1e4, zeta_d4=1e6)

    nw.solve("design")
    assert nw.status == 2
    assert np.isnan(c1.T.val_SI)
    assert instance.zeta_d4.val == 1e6
    assert np.isnan(instance.zeta_d4.val_SI)


def test_nonconverged_simulation_does_not_overwrite_component_specification_2():
    """This creates a result, that shows as converged but actually it did not
    because of internal convergence helpers. It tests, that the user specified
    input is not overwritten by the erroneous result
    """
    nw = Network()
    nw.units.set_defaults(
        pressure="bar", pressure_difference="bar", temperature="degC"
    )

    inflow = Source("source")
    outflow = Sink("sink")
    instance = SimpleHeatExchanger("heatexchanger")

    c1 = Connection(inflow, "out1", instance, "in1", label="c1")
    c2 = Connection(instance, "out1", outflow, "in1", label="c2")

    nw.add_conns(c1, c2)


    c1.set_attr(fluid={"H2O": 1}, T=30, p=1)
    c2.set_attr(T=19.5)
    instance.set_attr(Tamb=20, kA=500, pr=1)
    nw.solve("design")

    assert nw.residual < 1e-3  # residual shows convergence
    assert nw.status == 2  # status shows non-convergence

    assert np.isnan(instance.kA.val_SI)  # calculated SI value is not equal to inputted value
    assert instance.kA.val == 500  # inputted value stays the same

    # recalculation works, because old kA input is correctly retained
    c2.set_attr(T=20.2)
    nw.solve("design")
    assert nw.status == 0


def test_offdesign_of_component_parameter_group():

    nw = Network()
    nw.units.set_defaults(
        pressure="bar", pressure_difference="bar", temperature="degC"
    )

    inflow = Source("source")
    outflow = Sink("sink")
    instance = Pipe("pipe")

    c1 = Connection(inflow, "out1", instance, "in1", label="c1")
    c2 = Connection(instance, "out1", outflow, "in1", label="c2")

    nw.add_conns(c1, c2)

    c1.set_attr(fluid={"H2O": 1}, T=30, p=1, m=1)
    c2.set_attr(T=19.5, p=0.9, design=["p"])
    instance.set_attr(D=0.1, ks=0.0001, L=50, offdesign=["darcy_group"])
    nw.solve("design")
    nw.assert_convergence()
    assert not instance.darcy_group.is_set

    design_state = nw.save(as_dict=True)
    nw.solve("offdesign", design_path=design_state)
    nw.assert_convergence()
    assert instance.darcy_group.is_set


def test_design_of_component_parameter_group():

    nw = Network()
    nw.units.set_defaults(
        pressure="bar", pressure_difference="bar", temperature="degC"
    )

    inflow = Source("source")
    outflow = Sink("sink")
    instance = Pipe("pipe")

    c1 = Connection(inflow, "out1", instance, "in1", label="c1")
    c2 = Connection(instance, "out1", outflow, "in1", label="c2")

    nw.add_conns(c1, c2)

    c1.set_attr(fluid={"H2O": 1}, T=30, p=1, m=1)
    c2.set_attr(T=19.5, offdesign=["p"])
    instance.set_attr(D=0.1, ks=0.0001, L=50, design=["darcy_group"])
    nw.solve("design")
    nw.assert_convergence()
    assert instance.darcy_group.is_set

    design_state = nw.save(as_dict=True)
    nw.solve("offdesign", design_path=design_state)
    nw.assert_convergence()
    assert not instance.darcy_group.is_set


class WaterElectrolyzer(WaterElectrolyzer):

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        constraints['mass_flow_constraints'] = dc_cmc(**{
            'func': self.reactor_mass_flow_func,
            'deriv': self.reactor_mass_flow_deriv,
            'dependents': self.reactor_mass_flow_dependents,
            'num_eq_sets': 2
        })
        return constraints

    def reactor_mass_flow_func(self):
        r"""
        Equations for mass conservation.

        Returns
        -------
        residual : list
            Residual values of equation.

            .. math::

                O_2 = \frac{M_{O_2}}{M_{O_2} + 2 \cdot M_{H_2}}\\
                0 =\dot{m}_\mathrm{in,1}-\dot{m}_\mathrm{out,1}\\
                0=O_2\cdot\dot{m}_\mathrm{H_{2}O,in,2}-
                \dot{m}_\mathrm{O_2,out,2}\\
                0 = \left(1 - O_2\right) \cdot \dot{m}_\mathrm{H_{2}O,in,2} -
                \dot{m}_\mathrm{H_2,out,3}
        """
        # calculate the ratio of o2 in water
        M_o2 = self.outl[1].fluid.wrapper[self.o2]._molar_mass
        M_h2 = self.outl[2].fluid.wrapper[self.h2]._molar_mass

        o2 = M_o2 / (M_o2 + 2 * M_h2)
        # equations for mass flow balance electrolyzer
        residual = []
        residual += [o2 * self.inl[1].m.val_SI - self.outl[1].m.val_SI]
        residual += [(1 - o2) * self.inl[1].m.val_SI - self.outl[2].m.val_SI]
        return np.array(residual)

    def reactor_mass_flow_deriv(self, increment_filter, k, dependents=None):
        r"""
        Calculate the partial derivatives for all mass flow balance equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the mass flow equations.
        """
        deriv_to_m_in = _numeric_deriv(
            self.inl[1].m._reference_container, self.reactor_mass_flow_func
        )
        deriv_to_m_out1 = _numeric_deriv(
            self.outl[1].m._reference_container, self.reactor_mass_flow_func
        )
        deriv_to_m_out2 = _numeric_deriv(
            self.outl[2].m._reference_container, self.reactor_mass_flow_func
        )

        self._partial_derivative(self.inl[1].m, k, deriv_to_m_in[0])
        self._partial_derivative(self.outl[1].m, k, deriv_to_m_out1[0])

        k += 1
        self._partial_derivative(self.inl[1].m, k, deriv_to_m_in[1])
        self._partial_derivative(self.outl[2].m, k, deriv_to_m_out2[1])


def test_component_with_numpy_array_in_residual():
    nw = Network()
    nw.units.set_defaults(**{
        "pressure": "bar", "pressure_difference": "bar",
        "temperature": "degC", "power": "MW"
    })
    instance = WaterElectrolyzer('electrolyzer')

    fw = Source('feed water')
    cw_in = Source('cooling water')
    o2 = Sink('oxygen sink')
    h2 = Sink('hydrogen sink')
    cw_out = Sink('cooling water sink')

    instance.set_attr(pr=0.99, eta=1)

    cw_el = Connection(
        cw_in, 'out1', instance, 'in1', fluid={'H2O': 1}, T=20, p=1
    )
    el_cw = Connection(instance, 'out1', cw_out, 'in1', T=45)

    nw.add_conns(cw_el, el_cw)

    fw_el = Connection(fw, 'out1', instance, 'in2', label='h2o')
    el_o2 = Connection(instance, 'out2', o2, 'in1')
    el_h2 = Connection(instance, 'out3', h2, 'in1', label='h2')

    nw.add_conns(fw_el, el_o2, el_h2)

    fw_el.set_attr(T=25, p=1)
    el_h2.set_attr(T=25)
    instance.set_attr(P=2.5)

    nw.solve('design')
    nw.assert_convergence()


def test_generates_fluid_wrapper_branches_with_inherited_component():

    class FakeSource(Source):
        pass

    nw = Network()

    so = FakeSource('feed water')
    si = Sink('oxygen sink')

    c1 = Connection(so, 'out1', si, 'in1', fluid={'H2O': 1}, T=20 + 273.15, p=1e5)
    nw.add_conns(c1)

    nw.solve("design", init_only=True)


def test_fluid_kwargs_propagation():
    nw = Network()
    nw.units.set_defaults(
        temperature="°C", pressure="bar", pressure_difference="bar"
    )

    pipe = SimpleHeatExchanger("pipe")

    so = Source("source")
    si = Sink("sink")

    c1 = Connection(so, "out1", pipe, "in1", label="c1")
    c2 = Connection(pipe, "out1", si, "in1", label="c2")

    nw.add_conns(c1, c2)

    fluid_kwargs = {
        "temperature_data": np.array([273.15, 373.15]),
        "density_data": np.array([1000, 1100]),
        "heat_capacity_data": np.array([4000, 4100]),
        "viscosity_data": np.array([0.05, 0.00025]),
        "conductivity_data": np.array([0.1425, 0.135])
    }

    c1.set_attr(
        fluid={"f": 1},
        fluid_engines={"f": IncompressibleFluidWrapper},
        fluid_wrapper_kwargs={"f": fluid_kwargs},
        p=1, T=30
    )
    c2.set_attr(p=0.9, T=50)
    pipe.set_attr(Q=1500)

    nw.solve("design")

    # 50 °C is exactly half of range
    # heat capacity is implicitly tested, as it is required to find the
    # temperature from the enthalpy passed into the function
    assert approx(1 / c2.calc_vol()) == 1050
    assert approx(
        conductivity_mix_ph(c2.p.val_SI, c2.h.val_SI, c2.fluid_data)
    ) == 0.13875


def test_skip_postprocessing():
    nw = Network()
    nw.units.set_defaults(temperature="°C", pressure="bar", pressure_difference="bar")

    pipe = SimpleHeatExchanger("pipe")

    so = Source("source")
    si = Sink("sink")

    c1 = Connection(so, "out1", pipe, "in1", label="c1")
    c2 = Connection(pipe, "out1", si, "in1", label="c2")

    nw.add_conns(c1, c2)

    fluid_kwargs = {
        "temperature_data": np.array([273.15, 373.15]),
        "density_data": np.array([1000, 1100]),
        "heat_capacity_data": np.array([4000, 4100]),
        "viscosity_data": np.array([0.05, 0.00025]),
        "conductivity_data": np.array([0.1425, 0.135])
    }

    c1.set_attr(
        fluid={"f": 1},
        fluid_engines={"f": IncompressibleFluidWrapper},
        fluid_wrapper_kwargs={"f": fluid_kwargs},
        p=1, T=30
    )
    c2.set_attr(p=0.9, T=50)
    pipe.set_attr(Q=1500)

    nw.solve("design", skip_postprocess=True)
    nw.assert_convergence()

    assert np.isnan(pipe.pr.val)
    assert np.isnan(pipe.dp.val)
    assert np.isnan(c2.v.val)
    assert np.isnan(c1.s.val)


def test_setting_ref_on_hex_leads_to_linear_dependency():
    nw = Network()
    nw.units.set_defaults(
        temperature="°C",
        pressure="bar", pressure_difference="bar"
    )

    so1 = Source("source1")
    si1 = Sink("sink1")
    so2 = Source("source2")
    si2 = Sink("sink2")

    hex = MovingBoundaryHeatExchanger("hex")

    c1 = Connection(so1, "out1", hex, "in1", label="c1")
    c2 = Connection(hex, "out1", si1, "in1", label="c2")
    c3 = Connection(so2, "out1", hex, "in2", label="c3")
    c4 = Connection(hex, "out2", si2, "in1", label="c4")

    nw.add_conns(c1, c2, c3, c4)

    c1.set_attr(fluid={"water": 1}, td_dew=10, T=70, m=1)
    c2.set_attr(x=0.5)

    c3.set_attr(fluid={"air": 1}, T=40, p=1)
    c4.set_attr(T=Ref(c3, 1, 5))

    hex.set_attr(dp1=0, dp2=0)
    nw.solve("design")
    c3.set_attr(T=c4.T.val)
    nw.solve("design")
    assert nw.status == 0


def test_export_creates_nonexistent_directory(tmp_path):
    nw = Network()
    nw.units.set_defaults(
        temperature="°C", pressure="bar", pressure_difference="bar"
    )
    nw.iterinfo = False
    so = Source("source")
    si = Sink("sink")
    c = Connection(so, "out1", si, "in1")
    nw.add_conns(c)
    c.set_attr(fluid={"water": 1}, T=25, p=1, m=1)
    nw.solve("design")
    nw.assert_convergence()

    export_path = tmp_path / "new_subdir" / "network.json"
    nw.export(str(export_path))
    assert export_path.exists()
    design_path = tmp_path / "new_subdir" / "design.json"
    nw.save(design_path)
    assert design_path.exists()


class TestBackwardsCompatibility:
    """Verify that save/export files written by v0.9.x are still readable."""

    _HERE = os.path.dirname(os.path.abspath(__file__))

    def test_v09_export_design_and_offdesign(self):
        """v0.9 export + flat design state: import, design, offdesign all work."""
        nw = Network.from_json(os.path.join(self._HERE, "_exported_nwk.json"))
        nw.iterinfo = False
        nw.solve("design")
        nw.assert_convergence()
        nw.get_comp('compressor').set_attr(igva='var')
        nw.solve("offdesign", design_path=os.path.join(self._HERE, "_design_state.json"))
        nw.assert_convergence()

    # integrated this test here because it has quite a few different components
    def test_export_folder_csv_files(self, tmp_path):
        nw = Network.from_json(os.path.join(self._HERE, "_exported_nwk.json"))
        nw.solve("design")
        nw.save_csv(tmp_path)
        assert (tmp_path / "Component" / "TurboCompressor.csv").exists()
        assert (tmp_path / "Connection" / "Connection.csv").exists()
