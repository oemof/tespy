# -*- coding: utf-8

"""Module for testing network properties.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_exergy_and_entropy_analysis.py

SPDX-License-Identifier: MIT
"""

from exerpy import ExergyAnalysis
from pytest import raises

from tespy.components import Compressor
from tespy.components import CycleCloser
from tespy.components import Generator
from tespy.components import HeatSource
from tespy.components import Merge
from tespy.components import Motor
from tespy.components import PowerBus
from tespy.components import PowerSink
from tespy.components import PowerSource
from tespy.components import Pump
from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Splitter
from tespy.components import Turbine
from tespy.components import Valve
from tespy.connections import Connection
from tespy.connections import HeatConnection
from tespy.connections import PowerConnection
from tespy.networks import Network
from tespy.tools.global_vars import ERR


class TestClausiusRankine:

    def setup_method(self):
        """Set up clausius rankine cycle with turbine driven feed water pump."""
        self.Tamb = 20
        self.pamb = 1
        self.nw = Network()
        self.nw.units.set_defaults(**{
            "pressure": "bar", "pressure_difference": "bar",
            "temperature": "degC", "enthalpy": "kJ/kg"
        })

        # create components
        splitter1 = Splitter("splitter 1")
        merge1 = Merge("merge 1")
        turb = Turbine("turbine")
        fwp_turb = Turbine("feed water pump turbine")
        condenser = SimpleHeatExchanger("condenser", dissipative=True)
        fwp = Pump("pump")
        steam_generator = SimpleHeatExchanger("steam generator", dissipative=False)
        cycle_close = CycleCloser("cycle closer")

        # power components
        gen = Generator("generator")
        grid = PowerSink("grid")
        fwp_shaft = Generator("fwp shaft")
        heat_src = HeatSource("heat source")

        gen.set_attr(eta=1)

        # create connections
        fs_in = Connection(cycle_close, "out1", splitter1, "in1", label="fs")
        fs_fwpt = Connection(splitter1, "out1", fwp_turb, "in1")
        fs_t = Connection(splitter1, "out2", turb, "in1")
        fwpt_ws = Connection(fwp_turb, "out1", merge1, "in1")
        t_ws = Connection(turb, "out1", merge1, "in2")
        ws = Connection(merge1, "out1", condenser, "in1")
        cond = Connection(condenser, "out1", fwp, "in1", label="cond")
        fw = Connection(fwp, "out1", steam_generator, "in1", label="fw")
        fs_out = Connection(steam_generator, "out1", cycle_close, "in1")
        self.nw.add_conns(
            fs_in, fs_fwpt, fs_t, fwpt_ws, t_ws, ws, cond, fw, fs_out
        )

        fwp_shaft.set_attr(eta=1)

        # power connections: main turbine -> generator -> grid
        e1 = PowerConnection(turb, "power", gen, "power_in", label="e1")
        e2 = PowerConnection(gen, "power_out", grid, "power", label="e2")
        # power connections: fwp turbine -> shaft -> fwp
        e3 = PowerConnection(fwp_turb, "power", fwp_shaft, "power_in", label="e3")
        e4 = PowerConnection(fwp_shaft, "power_out", fwp, "power", label="e4")
        self.nw.add_conns(e1, e2, e3, e4)

        # heat connection: heat source -> steam generator
        h1 = HeatConnection(heat_src, "heat", steam_generator, "heat", label="h1")
        self.nw.add_conns(h1)

        # component parameters
        turb.set_attr(eta_s=1)
        fwp_turb.set_attr(eta_s=1)
        condenser.set_attr(pr=1, dissipative=True)
        fwp.set_attr(eta_s=1)
        steam_generator.set_attr(pr=1, dissipative=False)

        # connection parameters
        fs_in.set_attr(m=10, p=120, T=600, fluid={"water": 1})
        cond.set_attr(T=self.Tamb, x=0.0)

        # solve network
        self.nw.solve("design")
        self.nw.print_results()
        self.nw.assert_convergence()

    def test_exergy_analysis_perfect_cycle(self):
        """Test exergy analysis in the perfect clausius rankine cycle."""
        E_F = {"inputs": ["h1"]}
        E_P = {"inputs": ["e2"]}
        ean = ExergyAnalysis.from_tespy(self.nw, self.Tamb + 273.15, self.pamb * 1e5)
        ean.analyse(E_F=E_F, E_P=E_P)

        msg = (
            f"Exergy destruction of this network must be 0 (smaller than "
            f"{ERR ** 0.5}) for this test but is {round(abs(ean.E_D), 4)} .")
        assert abs(ean.E_D) <= 1e-2, msg

        msg = f"Exergy efficiency of this network must be 1 for this test but is {round(ean.epsilon, 4)} ."
        assert round(ean.epsilon, 4) == 1, msg

        exergy_balance = ean.E_F - ean.E_P - ean.E_L - ean.E_D
        msg = (
            f"Exergy balance must be closed (residual value smaller than "
            f"{ERR ** 0.5}) for this test but is {round(abs(exergy_balance), 4)} .")
        assert abs(exergy_balance) <= ERR ** 0.5, msg

        msg = (
            f"Fuel exergy and product exergy must be identical for this test. "
            f"Fuel exergy value: {round(ean.E_F, 2)}. Product exergy value: {round(ean.E_P, 2)}.")
        delta = round(abs(ean.E_F - ean.E_P), 2)
        assert delta < 1e-2, msg

    def test_exergy_with_non_ideal_components(self):
        """Test exergy analysis with non-ideal components."""
        self.nw.get_comp("steam generator").set_attr(pr=0.9)
        self.nw.get_comp("turbine").set_attr(eta_s=0.9)
        self.nw.get_comp("feed water pump turbine").set_attr(eta_s=0.85)
        self.nw.get_comp("pump").set_attr(eta_s=0.75)
        self.nw.get_comp("generator").set_attr(eta=0.98)
        self.nw.get_comp("fwp shaft").set_attr(eta=0.98)
        self.nw.get_conn("cond").set_attr(T=self.Tamb + 5)

        self.nw.solve("design")
        self.nw.assert_convergence()

        E_F = {"inputs": ["h1"]}
        E_P = {"inputs": ["e2"]}
        ean = ExergyAnalysis.from_tespy(self.nw, self.Tamb + 273.15, self.pamb * 1e5)
        ean.analyse(E_F=E_F, E_P=E_P)

        exergy_balance = ean.E_F - ean.E_P - ean.E_L - ean.E_D
        msg = (
            f"Exergy balance must be closed (residual value smaller than "
            f"{ERR ** 0.5}) for this test but is {round(abs(exergy_balance), 4)} .")
        assert abs(exergy_balance) <= ERR ** 0.5, msg

        df, _, _ = ean.exergy_results(print_results=False)
        assert len(df) > 0, "Exergy results must not be empty."

    def test_exergy_analysis_missing_E_F_E_P_information(self):
        """Test exergy analysis errors with invalid connection labels."""
        ean = ExergyAnalysis.from_tespy(self.nw, self.Tamb + 273.15, self.pamb * 1e5)
        with raises(ValueError):
            ean.analyse(
                E_F={"inputs": ["nonexistent_connection"]},
                E_P={"inputs": ["e2"]},
            )
        with raises(ValueError):
            ean.analyse(
                E_F={"inputs": ["h1"]},
                E_P={"inputs": ["nonexistent_connection"]},
            )


class TestRefrigerator:

    def setup_method(self):
        """Set up simple refrigerator."""
        self.Tamb = 20
        self.pamb = 1
        self.nw = Network()
        self.nw.units.set_defaults(**{
            "pressure": "bar", "pressure_difference": "bar",
            "temperature": "degC", "enthalpy": "kJ/kg"
        })

        # create components
        va = Valve("expansion valve")
        cp = Compressor("compressor")
        cond = SimpleHeatExchanger("condenser", dissipative=True)
        eva = SimpleHeatExchanger("evaporator", dissipative=False)
        cc = CycleCloser("cycle closer")

        # power components
        grid = PowerSource("grid")
        motor = Motor("motor")
        cold_space = HeatSource("cold space")

        motor.set_attr(eta=1)

        # create connections
        cc_cp = Connection(cc, "out1", cp, "in1", label="from eva")
        cp_cond = Connection(cp, "out1", cond, "in1", label="to cond")
        cond_va = Connection(cond, "out1", va, "in1", label="from cond")
        va_eva = Connection(va, "out1", eva, "in1", label="to eva")
        eva_cc = Connection(eva, "out1", cc, "in1")
        self.nw.add_conns(cc_cp, cp_cond, cond_va, va_eva, eva_cc)

        # power connections: grid -> motor -> compressor
        e1 = PowerConnection(grid, "power", motor, "power_in", label="e1")
        e2 = PowerConnection(motor, "power_out", cp, "power", label="e2")
        self.nw.add_conns(e1, e2)

        # heat connection: cold space -> evaporator
        h1 = HeatConnection(cold_space, "heat", eva, "heat", label="h1")
        self.nw.add_conns(h1)

        # component parameters
        cp.set_attr(eta_s=0.9)
        cond.set_attr(pr=0.97)
        eva.set_attr(pr=0.96)

        # connection parameters
        cc_cp.set_attr(m=1, x=1, T=-25, fluid={"R134a": 1})
        cond_va.set_attr(x=0, T=self.Tamb + 1)

        # solve network
        self.nw.solve("design")
        self.nw.assert_convergence()

    def test_exergy_analysis(self):
        """Test exergy analysis at product exergy with T < Tamb."""
        # E_P sign: h1 carries heat from cold_space into evaporator at T < T0,
        # which has negative exergy as an input, so it appears in outputs.
        E_F = {"inputs": ["e1"]}
        E_P = {"outputs": ["h1"]}
        ean = ExergyAnalysis.from_tespy(self.nw, self.Tamb + 273.15, self.pamb * 1e5)
        ean.analyse(E_F=E_F, E_P=E_P)

        exergy_balance = ean.E_F - ean.E_P - ean.E_L - ean.E_D
        msg = (
            f"Exergy balance must be closed (residual value smaller than "
            f"{ERR ** 0.5}) for this test but is {round(abs(exergy_balance), 4)} .")
        assert abs(exergy_balance) <= ERR ** 0.5, msg


class TestCompressedAirIn:

    def setup_method(self):
        """Set up air compressor."""
        self.Tamb = 20
        self.pamb = 1

        self.nw = Network()
        self.nw.units.set_defaults(**{
            "pressure": "bar", "pressure_difference": "bar",
            "temperature": "degC", "enthalpy": "kJ/kg"
        })

        # components
        amb = Source("air intake")
        cp = Compressor("compressor")
        cooler = SimpleHeatExchanger("cooling", dissipative=True)
        cas = Sink("compressed air storage")

        # power components
        grid = PowerSource("grid")
        motor = Motor("motor")
        motor.set_attr(eta=1)

        # create connections
        amb_cp = Connection(amb, "out1", cp, "in1", label="amb_cp")
        cp_cool = Connection(cp, "out1", cooler, "in1")
        cool_cas = Connection(cooler, "out1", cas, "in1", label="cool_cas")
        self.nw.add_conns(amb_cp, cp_cool, cool_cas)

        # power connections: grid -> motor -> compressor
        e1 = PowerConnection(grid, "power", motor, "power_in", label="e1")
        e2 = PowerConnection(motor, "power_out", cp, "power", label="e2")
        self.nw.add_conns(e1, e2)

        # component parameters
        cp.set_attr(eta_s=1)
        cooler.set_attr(pr=1)

        # connection parameters
        amb_cp.set_attr(m=2, T=self.Tamb, p=self.pamb, fluid={"Air": 1})
        cool_cas.set_attr(T=self.Tamb, p=10)

        # solve network
        self.nw.solve("design")
        self.nw.assert_convergence()

    def test_exergy_analysis(self):
        """Test exergy analysis at product exergy with T < Tamb."""
        E_F = {"inputs": ["e1"]}
        E_P = {"inputs": ["cool_cas"], "outputs": ["amb_cp"]}
        ean = ExergyAnalysis.from_tespy(self.nw, self.Tamb + 273.15, self.pamb * 1e5)
        ean.analyse(E_F=E_F, E_P=E_P)

        exergy_balance = ean.E_F - ean.E_P - ean.E_L - ean.E_D
        msg = (
            f"Exergy balance must be closed (residual value smaller than "
            f"{ERR ** 0.5}) for this test but is {round(abs(exergy_balance), 4)} .")
        assert abs(exergy_balance) <= ERR ** 0.5, msg


class TestCompressedAirOut:

    def setup_method(self):
        """Set up air compressed air turbine."""
        self.Tamb = 20
        self.pamb = 1

        self.nw = Network()
        self.nw.units.set_defaults(**{
            "pressure": "bar", "pressure_difference": "bar",
            "temperature": "degC", "enthalpy": "kJ/kg"
        })

        # components
        cas = Source("compressed air storage")
        reheater = SimpleHeatExchanger("reheating", dissipative=False)
        turb = Turbine("turbine")
        amb = Sink("air outlet")

        # power components
        gen = Generator("generator")
        grid = PowerSink("grid")
        heat_src = HeatSource("heat source")
        gen.set_attr(eta=1)

        # create connections
        cas_reheater = Connection(cas, "out1", reheater, "in1", label="cas_in")
        reheater_turb = Connection(reheater, "out1", turb, "in1")
        turb_amb = Connection(turb, "out1", amb, "in1", label="outlet")
        self.nw.add_conns(cas_reheater, reheater_turb, turb_amb)

        # power connections: turbine -> generator -> grid
        e1 = PowerConnection(turb, "power", gen, "power_in", label="e1")
        e2 = PowerConnection(gen, "power_out", grid, "power", label="e2")
        self.nw.add_conns(e1, e2)

        # heat connection: heat source -> reheater
        h1 = HeatConnection(heat_src, "heat", reheater, "heat", label="h1")
        self.nw.add_conns(h1)

        # component parameters
        turb.set_attr(eta_s=1)
        reheater.set_attr(pr=1)

        # connection parameters
        cas_reheater.set_attr(m=2, T=self.Tamb, p=10, fluid={"Air": 1})
        turb_amb.set_attr(p=self.pamb, T=self.Tamb)

        # solve network
        self.nw.solve("design")
        self.nw.assert_convergence()

    def test_exergy_analysis_no_destruction_just_loss(self):
        """Test exergy analysis with all """
        E_F = {"inputs": ["cas_in", "h1"]}
        E_P = {"inputs": ["e2"]}
        E_L = {"inputs": ["outlet"]}

        ean = ExergyAnalysis.from_tespy(self.nw, self.Tamb + 273.15, self.pamb * 1e5)
        ean.analyse(E_F=E_F, E_P=E_P, E_L=E_L)

        msg = (
            "Exergy efficiency must be equal to 1.0 for this test but is "
            f"{round(ean.epsilon, 4)}."
        )
        assert round(ean.epsilon, 4) == 1, msg

        msg = (
            "Exergy destruction must be equal to 0.0 for this test but is "
            f"{round(ean.E_D, 4)}."
        )
        assert round(ean.E_D, 4) == 0, msg

        c = self.nw.get_conn("outlet")
        c.set_attr(T=self.Tamb - 20)
        self.nw.solve("design")
        self.nw.assert_convergence()

        ean = ExergyAnalysis.from_tespy(self.nw, self.Tamb + 273.15, self.pamb * 1e5)
        ean.analyse(E_F=E_F, E_P=E_P, E_L=E_L)

        msg = (
            "Exergy destruction must be equal to 0.0 for this test but is "
            f"{round(ean.E_D, 4)}."
        )
        assert round(ean.E_D, 4) == 0, msg

        msg = (
            f"Exergy loss must be equal to {round(c.Ex_physical, 4)} "
            f"for this test but is {round(ean.E_L, 4)}."
        )
        assert round(ean.E_L, 4) == round(c.Ex_physical, 4), msg


class TestCompression:

    def setup_method(self):
        self.Tamb = 20
        self.pamb = 1

        self.nw = Network()
        self.nw.units.set_defaults(**{
            "pressure": "bar", "pressure_difference": "bar",
            "temperature": "degC", "enthalpy": "kJ/kg"
        })

        # components
        so = Source("inlet")
        cp = Compressor("compressor")
        si = Sink("outlet")

        # power components: grid -> motor (eta=0.9) -> compressor
        grid = PowerSource("grid")
        motor = Motor("motor")
        motor.set_attr(eta=0.9)

        # create connections
        c1 = Connection(so, "out1", cp, "in1", "1")
        c2 = Connection(cp, "out1", si, "in1", "2")
        self.nw.add_conns(c1, c2)

        e1 = PowerConnection(grid, "power", motor, "power_in", label="e1")
        e2 = PowerConnection(motor, "power_out", cp, "power", label="e2")
        self.nw.add_conns(e1, e2)

        # component parameters
        cp.set_attr(eta_s=0.85, pr=5)

        # connection parameters
        c1.set_attr(m=2, T=self.Tamb, p=self.pamb, fluid={"Air": 1})

        # solve network
        self.nw.solve("design")

    def test_larger_T0(self):
        self.nw.get_conn("1").set_attr(T=self.Tamb + 10)
        self.nw.solve("design")
        self.run_analysis()

    def test_T0_cross(self):
        self.nw.get_conn("1").set_attr(T=self.Tamb - 30)
        self.nw.solve("design")
        self.run_analysis()

    def test_smaller_T0(self):
        self.nw.get_conn("1").set_attr(T=None)
        self.nw.get_conn("2").set_attr(T=self.Tamb - 10)
        self.nw.solve("design")
        self.run_analysis()

    def run_analysis(self):
        E_F = {"inputs": ["e1"]}
        E_P = {"inputs": ["2"], "outputs": ["1"]}
        ean = ExergyAnalysis.from_tespy(self.nw, self.Tamb + 273.15, self.pamb * 1e5)
        ean.analyse(E_F=E_F, E_P=E_P)

        exergy_balance = ean.E_F - ean.E_P - ean.E_L - ean.E_D
        msg = (
            "Exergy balance must be closed (residual value smaller than "
            f"{ERR ** 0.5}) for this test but is "
            f"{round(abs(exergy_balance), 4)}."
        )
        assert abs(exergy_balance) <= ERR ** 0.5, msg

        df, _, _ = ean.exergy_results(print_results=False)
        E_D_components = df.loc[df["Component"] != "TOT", "E_D [kW]"].sum() * 1e3
        E_D_nw = ean.E_D
        msg = (
            "The exergy destruction of the aggregated components "
            f"({round(E_D_components)}) must be equal to the exergy "
            f"destruction of the network ({round(E_D_nw)})."
        )
        assert abs(E_D_components - E_D_nw) < 1, msg


class TestExpansion:

    def setup_method(self):
        self.Tamb = 20
        self.pamb = 1

        self.nw = Network()
        self.nw.units.set_defaults(**{
            "pressure": "bar", "pressure_difference": "bar",
            "temperature": "degC", "enthalpy": "kJ/kg"
        })

        # components
        so = Source("inlet")
        tu = Turbine("compressor")
        si = Sink("outlet")

        # power components: turbine -> generator (eta=0.9) -> grid
        gen = Generator("generator")
        grid = PowerSink("grid")
        gen.set_attr(eta=0.9)

        # create connections
        c1 = Connection(so, "out1", tu, "in1", "1")
        c2 = Connection(tu, "out1", si, "in1", "2")
        self.nw.add_conns(c1, c2)

        e1 = PowerConnection(tu, "power", gen, "power_in", label="e1")
        e2 = PowerConnection(gen, "power_out", grid, "power", label="e2")
        self.nw.add_conns(e1, e2)

        # component parameters
        tu.set_attr(eta_s=0.85, pr=1 / 5)

        # connection parameters
        c1.set_attr(m=2, p=10, fluid={"Air": 1})
        c2.set_attr(T=self.Tamb)

        # solve network
        self.nw.solve("design")

    def test_larger_T0(self):
        self.nw.get_conn("2").set_attr(T=self.Tamb + 10)
        self.nw.solve("design")
        self.run_analysis()

    def test_T0_cross(self):
        self.nw.get_conn("2").set_attr(T=self.Tamb - 30)
        self.nw.solve("design")
        self.run_analysis()

    def test_smaller_T0(self):
        self.nw.get_conn("1").set_attr(T=self.Tamb - 10)
        self.nw.get_conn("2").set_attr(T=None)
        self.nw.solve("design")
        self.run_analysis()

    def run_analysis(self):
        E_F = {"inputs": ["1"], "outputs": ["2"]}
        E_P = {"inputs": ["e2"]}
        ean = ExergyAnalysis.from_tespy(self.nw, self.Tamb + 273.15, self.pamb * 1e5)
        ean.analyse(E_F=E_F, E_P=E_P)

        exergy_balance = ean.E_F - ean.E_P - ean.E_L - ean.E_D
        msg = (
            "Exergy balance must be closed (residual value smaller than "
            f"{ERR ** 0.5}) for this test but is "
            f"{round(abs(exergy_balance), 4)}."
        )
        assert abs(exergy_balance) <= ERR ** 0.5, msg

        df, _, _ = ean.exergy_results(print_results=False)
        E_D_components = df.loc[df["Component"] != "TOT", "E_D [kW]"].sum() * 1e3
        E_D_nw = ean.E_D
        msg = (
            "The exergy destruction of the aggregated components "
            f"({round(E_D_components)}) must be equal to the exergy "
            f"destruction of the network ({round(E_D_nw)})."
        )
        assert abs(E_D_components - E_D_nw) < 1, msg
