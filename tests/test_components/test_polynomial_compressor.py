# -*- coding: utf-8

"""Module for testing components of type turbomachinery.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_turbomachinery.py

SPDX-License-Identifier: MIT
"""
import numpy as np
import pandas as pd
import pytest
from CoolProp.CoolProp import PropsSI as PSI

from tespy.components import PolynomialCompressor
from tespy.components import PowerSource
from tespy.components import Sink
from tespy.components import Source
from tespy.components.displacementmachinery.polynomial_compressor import calc_EN12900
from tespy.components.displacementmachinery.polynomial_compressor import fit_EN12900
from tespy.components.displacementmachinery.polynomial_compressor import generate_eta_polys_from_data
from tespy.components.displacementmachinery.polynomial_compressor import generate_eta_polys_from_power_and_cooling_polys
from tespy.connections import Connection
from tespy.connections import PowerConnection
from tespy.networks import Network


@pytest.fixture
def power_data():
    power = pd.DataFrame(
        columns=[10,7.5,5,0,-5,-10], index=[30, 40, 50], dtype=float
    )

    power.loc[30] = [62.0,61.8,61.8,61.8,61.7,61.3]
    power.loc[40] = [78.0,78.0,78.0,78.0,77.7,76.8]
    power.loc[50] = [99.2,99.2,99.2,98.9,98.1,96.5]

    return power * 1000


@pytest.fixture
def cooling_data():
    cooling = pd.DataFrame(
        columns=[10,7.5,5,0,-5,-10], index=[30, 40, 50], dtype=float
    )

    cooling.loc[30] = [465600,424100,385500,316700,257900,208000]
    cooling.loc[40] = [418900,380400,344800,281400,227400,181600]
    cooling.loc[50] = [365900,331300,299200,242100,193700,152900]

    return cooling


@pytest.fixture
def reference_state():
    return {
        "T_sh": 20,
        "T_sc": 0,
        "rpm_poly": 50 * 60,
        "rpm_displacement": 20 * 60,
        "displacement": 214
    }


@pytest.fixture
def power_poly(power_data):
    return fit_EN12900(power_data.columns, power_data.index, power_data.values)


@pytest.fixture
def cooling_poly(cooling_data):
    return fit_EN12900(cooling_data.columns, cooling_data.index, cooling_data.values)


def test_fit_en12900(power_data, power_poly):
    t_evap, t_cond = np.meshgrid(power_data.columns, power_data.index)
    np.testing.assert_allclose(
        power_data.values.flatten(),
        calc_EN12900(power_poly, t_evap.flatten(), t_cond.flatten()),
        rtol=1e-3
    )


def test_integration_eta_polys(power_data, cooling_data, power_poly, cooling_poly, reference_state):
    fluid = "R134a"
    eta_s_poly, eta_vol_poly = generate_eta_polys_from_data(power_data, cooling_data, fluid, reference_state)
    t_evap = power_data.columns
    t_cond = power_data.index
    eta_s_poly2, eta_vol_poly2 = generate_eta_polys_from_power_and_cooling_polys(power_poly, cooling_poly, t_evap, t_cond, fluid, reference_state)

    np.testing.assert_allclose(eta_s_poly, eta_s_poly2)
    np.testing.assert_allclose(eta_vol_poly, eta_vol_poly2)


class TestPolynomialCompressor:

    def setup_network(self, instance):
        self.nw = Network()
        self.nw.units.set_defaults(**{
            "pressure": "bar", "temperature": "degC",
            "volumetric_flow": "m3/s"
        })
        self.source = Source('source')
        self.sink = Sink('sink')
        self.c1 = Connection(self.source, 'out1', instance, 'in1')
        self.c2 = Connection(instance, 'out1', self.sink, 'in1')
        self.nw.add_conns(self.c1, self.c2)

    def test_eta_s(self):
        """Test component properties of compressors."""
        instance = PolynomialCompressor('compressor')
        self.setup_network(instance)

        # compress NH3, other fluids in network are for turbine, pump, ...
        fl = {'R134a': 1}
        self.c1.set_attr(fluid=fl, m=1, x=1, T=5)
        self.c2.set_attr(p=6)
        instance.set_attr(eta_s=0.8, dissipation_ratio=0.1)
        self.nw.solve('design')
        self.nw.assert_convergence()

        assert round(instance.eta_s.val, 2) == 0.8

    def test_eta_vol(self, reference_state):
        """Test component properties of compressors."""
        instance = PolynomialCompressor('compressor')
        self.setup_network(instance)

        # compress NH3, other fluids in network are for turbine, pump, ...
        fl = {'R134a': 1}
        self.c1.set_attr(fluid=fl, x=1, T=5)
        self.c2.set_attr(p=6, T=150)
        instance.set_attr(
            eta_vol=0.8, rpm=1500, reference_state=reference_state,
            dissipation_ratio=0
        )
        self.nw.solve('design')
        self.nw.assert_convergence()

        assert round(self.c1.v.val, 2) == round(
            reference_state["displacement"] / 3600
            * instance.rpm.val / reference_state["rpm_displacement"]
            * instance.eta_vol.val,
            2
        )

    def test_eta_vol_with_Q_diss_rel(self, reference_state):
        """Test component properties of compressors."""
        instance = PolynomialCompressor('compressor')
        self.setup_network(instance)

        # compress NH3, other fluids in network are for turbine, pump, ...
        fl = {'R134a': 1}
        self.c1.set_attr(fluid=fl, x=1, T=5)
        self.c2.set_attr(p=6, T=150)
        instance.set_attr(
            eta_vol=0.8, rpm=1500, reference_state=reference_state,
            Q_diss_rel=0
        )
        self.nw.solve('design')
        self.nw.assert_convergence()

        assert round(self.c1.v.val, 2) == round(
            reference_state["displacement"] / 3600
            * instance.rpm.val / reference_state["rpm_displacement"]
            * instance.eta_vol.val,
            2
        )

    def test_eta_vol_var_rpm(self, reference_state):
        """Test component properties of compressors."""
        instance = PolynomialCompressor('compressor')
        self.setup_network(instance)

        # compress NH3, other fluids in network are for turbine, pump, ...
        fl = {'R134a': 1}
        self.c1.set_attr(fluid=fl, x=1, m=1, T=5)
        self.c2.set_attr(p=6, T=150)
        instance.set_attr(
            eta_vol=0.8, rpm="var", reference_state=reference_state,
            dissipation_ratio=0
        )
        self.nw.solve('design')
        self.nw.assert_convergence()

        assert round(self.c1.v.val, 2) == round(
            reference_state["displacement"] / 3600
            * instance.rpm.val / reference_state["rpm_displacement"]
            * instance.eta_vol.val,
            2
        )

    def test_power(self):
        instance = PolynomialCompressor('compressor')
        self.setup_network(instance)

        # compress NH3, other fluids in network are for turbine, pump, ...
        fl = {'R134a': 1}
        self.c1.set_attr(fluid=fl, x=0, T=0)
        p = PSI("P", "Q", 1, "T", 50 + 273.15, "R134a") / 1e5
        self.c2.set_attr(p=p)
        instance.set_attr(P=1e6, dissipation_ratio=0.1, eta_s=0.8)
        self.nw.solve("design")
        self.nw.assert_convergence()

        assert (
            round(instance.P.val, 3)
            ==
            round(
                (self.c2.h.val_SI - self.c1.h.val_SI)
                / (1 - instance.dissipation_ratio.val) * self.c1.m.val_SI,
                3
            )
        )

    def test_with_powerconnection(self):
        instance = PolynomialCompressor('compressor')
        self.setup_network(instance)

        # compress NH3, other fluids in network are for turbine, pump, ...
        fl = {'R134a': 1}
        self.c1.set_attr(fluid=fl, x=0, T=0)
        p = PSI("P", "Q", 1, "T", 50 + 273.15, "R134a") / 1e5
        self.c2.set_attr(p=p)
        instance.set_attr(dissipation_ratio=0.1, eta_s=0.8)

        grid = PowerSource("grid")
        e = PowerConnection(grid, "power", instance, "power", label="e")
        e.set_attr(E=1e6)

        self.nw.add_conns(e)

        self.nw.solve("design")
        self.nw.assert_convergence()

        assert round(e.E.val, 3) == round(instance.P.val, 3)

    def test_eta_poly(self, power_data, cooling_data, reference_state):
        """Test component properties of compressors."""
        instance = PolynomialCompressor('compressor')
        self.setup_network(instance)

        # compress NH3, other fluids in network are for turbine, pump, ...
        fl = {'R134a': 1}
        self.c1.set_attr(fluid=fl, td_dew=reference_state["T_sh"], T=20)
        p = PSI("P", "Q", 1, "T", 50 + 273.15, "R134a") / 1e5
        self.c2.set_attr(p=p)
        eta_s_poly, eta_vol_poly = generate_eta_polys_from_data(
            power_data, cooling_data, "R134a", reference_state
        )
        instance.set_attr(
            eta_vol_poly=eta_vol_poly, rpm=50 * 60,
            eta_s_poly=eta_s_poly, dissipation_ratio=0,
            reference_state=reference_state
        )
        self.nw.solve('design')
        self.nw.assert_convergence()

        assert round(self.c1.v.val, 2) == round(
            reference_state["displacement"] / 3600
            * instance.rpm.val / reference_state["rpm_displacement"]
            * instance.eta_vol.val,
            2
        )

        assert (
            abs((instance.P.val- power_data.loc[50, 0.0]))
            / power_data.loc[50, 0.0] <= 1e-2
        )
