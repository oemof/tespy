# -*- coding: utf-8 -*-

"""Module for testing two tespy simulation against each other.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_models/test_CGAM_model.py

SPDX-License-Identifier: MIT
"""
import os

import pandas as pd
from CoolProp.CoolProp import PropsSI as CPSI
from exerpy import ExergyAnalysis

from tespy.components import Compressor
from tespy.components import DiabaticCombustionChamber
from tespy.components import Drum
from tespy.components import Generator
from tespy.components import HeatExchanger
from tespy.components import PowerBus
from tespy.components import PowerSink
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Turbine
from tespy.connections import Connection
from tespy.connections import PowerConnection
from tespy.networks import Network


class TestCGAM:

    def setup_method(self):
        self.nwk = Network()
        self.nwk.units.set_defaults(**{
            "pressure": "bar", "pressure_difference": "bar",
            "temperature": "degC"
        })

        air_molar = {
            "O2": 0.2059, "N2": 0.7748, "CO2": 0.0003, "H2O": 0.019, "CH4": 0
        }
        molar_masses = {key: CPSI("M", key) * 1000 for key in air_molar}
        M_air = sum([air_molar[key] * molar_masses[key] for key in air_molar])

        air = {
            key: value / M_air * molar_masses[key]
            for key, value in air_molar.items()
        }

        water = {f: (0 if f != "H2O" else 1) for f in air}
        fuel = {f: (0 if f != "CH4" else 1) for f in air}

        amb = Source("ambient air")
        ch4 = Source("methane")
        fw = Source("feed water")

        ch = Sink("chimney")
        ls = Sink("live steam")

        cmp = Compressor("compressor")
        aph = HeatExchanger("air preheater")
        cb = DiabaticCombustionChamber("combustion chamber")
        tur = Turbine("gas turbine")

        eva = HeatExchanger("evaporator")
        eco = HeatExchanger("economizer")
        dr = Drum("drum")

        c1 = Connection(amb, "out1", cmp, "in1", label="1")
        c2 = Connection(cmp, "out1", aph, "in2", label="2")
        c3 = Connection(aph, "out2", cb, "in1", label="3")
        c10 = Connection(ch4, "out1", cb, "in2", label="10")

        self.nwk.add_conns(c1, c2, c3, c10)

        c4 = Connection(cb, "out1", tur, "in1", label="4")
        c5 = Connection(tur, "out1", aph, "in1", label="5")
        c6 = Connection(aph, "out1", eva, "in1", label="6")
        c6p = Connection(eva, "out1", eco, "in1", label="6p")
        c7 = Connection(eco, "out1", ch, "in1", label="7")

        self.nwk.add_conns(c4, c5, c6, c6p, c7)

        c8 = Connection(fw, "out1", eco, "in2", label="8")
        c8p = Connection(eco, "out2", dr, "in1", label="8p")
        c11 = Connection(dr, "out1", eva, "in2", label="11")
        c11p = Connection(eva, "out2", dr, "in2", label="11p")
        c9 = Connection(dr, "out2", ls, "in1", label="9")

        self.nwk.add_conns(c8, c8p, c11, c11p, c9)

        c8.set_attr(p=20, T=25, m=14, fluid=water)
        c1.set_attr(p=1.013, T=25, fluid=air)
        c10.set_attr(T=25, fluid=fuel, p=12)
        c7.set_attr(p=1.013)
        c3.set_attr(T=850 - 273.15)
        c8p.set_attr(td_bubble=15)
        c11p.set_attr(x=0.5)

        cmp.set_attr(pr=10, eta_s=0.86)
        cb.set_attr(eta=0.98, pr=0.95)
        tur.set_attr(eta_s=0.86)
        aph.set_attr(pr1=0.97, pr2=0.95)
        eva.set_attr(pr1=0.95 ** 0.5)
        eco.set_attr(pr1=0.95 ** 0.5, pr2=1)

        shaft = PowerBus("shaft", num_in=1, num_out=2)
        grid = PowerSink("grid")

        e1 = PowerConnection(tur, "power", shaft, "power_in1", label="e1")
        e2 = PowerConnection(shaft, "power_out1", cmp, "power", label="e2")
        e3 = PowerConnection(shaft, "power_out2", grid, "power", label="e3")

        self.nwk.add_conns(e1, e2, e3)

        e3.set_attr(E=30e6)
        c4.set_attr(T=1520 - 273.15)
        self.nwk.solve("design")
        self.nwk.assert_convergence()

        self.result = self.nwk.results["Connection"].copy()

        for idx in self.result.index:
            c = self.nwk.get_conn(idx)

            molarflow = {
                key: value / molar_masses[key] * c.m.val_SI * 1000
                for key, value in c.fluid.val.items()
            }
            molarflow_sum = sum(molarflow.values())
            molar = {
                key: value / molarflow_sum
                for key, value in molarflow.items()
            }

            self.result.loc[idx, molar.keys()] = molar

        self.result.loc["AC", "P"] = cmp.P.val
        self.result.loc["EXP", "P"] = tur.P.val

    def test_ebsilon(self):
        """Test the deviation with to an Ebsilon model"""
        path = os.path.dirname(__file__)
        ebsilon = pd.read_csv(path + "/cgam-ebsilon-results.csv", index_col=0)
        tespy = self.result.loc[ebsilon.index, ebsilon.columns]

        tespy.loc[:, ["h", "s"]] /= 1e3
        tespy.loc[:, ["T"]] += 273.15

        # set reference enthalpies and entropies

        # air
        air = ["1", "2", "3"]
        tespy.loc[air, ["h", "s"]] -= tespy.loc[air[0], ["h", "s"]]
        ebsilon.loc[air, ["h", "s"]] -= ebsilon.loc[air[0], ["h", "s"]]

        # CH4
        ch4 = ["10"]
        tespy.loc[ch4, ["h", "s"]] -= tespy.loc[ch4[0], ["h", "s"]]
        ebsilon.loc[ch4, ["h", "s"]] -= ebsilon.loc[ch4[0], ["h", "s"]]

        # flue gas
        fg = ["4", "5", "6", "6p", "7"]
        tespy.loc[fg, ["h", "s"]] -= tespy.loc[fg[0], ["h", "s"]]
        ebsilon.loc[fg, ["h", "s"]] -= ebsilon.loc[fg[0], ["h", "s"]]

        delta_abs = (tespy - ebsilon)
        delta_rel = (delta_abs / ebsilon).fillna(0)

        msg = "The deviation for all values must be lower than 0.005."
        assert (delta_rel.abs() < 5e-3).values.all(), msg

    def test_exergy_analysis(self):
        ean = ExergyAnalysis.from_tespy(self.nwk, 298.15, 1.013e5, chemExLib="Ahrendts")
        E_F = {"inputs": ["1", "10"], "outputs": []}
        E_P = {"inputs": ["e3", "9"], "outputs": ["8"]}
        E_L = {"inputs": ["7"], "outputs": []}
        ean.analyse(E_F=E_F, E_P=E_P, E_L=E_L)

        exergy_balance = ean.E_F - ean.E_P - ean.E_L - ean.E_D
        msg = (
            f"Exergy balance must be satisfied. Deviation is "
            f"{round(abs(exergy_balance), 4)}, but must be below {1e-6 ** 0.5}."
        )
        assert abs(exergy_balance) < 1e-6 ** 0.5, msg
