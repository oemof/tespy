# -*- coding: utf-8 -*-

"""Cross-check TESPy MovingBoundaryHeatExchanger against the Bell (2015)
oracle (HX.py) results stored in bell2015_five_cases.json.

For each case the test:
  1. Builds a MovingBoundaryHeatExchanger network with the oracle's fluid and
     flow inputs, without area_hot (so area_zones does not activate yet).
  2. Pins Q to the oracle value and solves for a warm initial state.
  3. Releases Q, sets area_hot, and resolves with oscillation_damping so
     area_zones drives the solution.
  4. Asserts that Q, T_h_out and T_c_out match the oracle within tolerance.

Q sign convention: TESPy stores Q = mdot_h * (h_out - h_in) < 0 for the hot
side.  The oracle reports positive heat transferred, so the warm-start pins
:code:`hx.Q = -Q_oracle` and assertions use :code:`abs(hx.Q.val_SI)`.
"""
import json
from pathlib import Path

import pytest

from tespy.components import MovingBoundaryHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.networks import Network

_ORACLE_PATH = Path(__file__).parent / "bell2015_five_cases.json"
_ORACLE = json.loads(_ORACLE_PATH.read_text())

Q_RTOL   = 1e-3   # 0.1 % relative tolerance on total heat transfer
T_ATOL_K = 0.5    # K  absolute tolerance on outlet temperatures


def _build_and_solve(
    fluid_h, T_h_in, p_h_in, mdot_h,
    fluid_c, T_c_in, p_c_in, mdot_c,
    alpha_sp, alpha_tp, area, R_cond,
    Q_oracle,
):
    nw = Network()
    nw.units.set_defaults(temperature="K", pressure="Pa", pressure_difference="Pa")
    nw.iterinfo = False

    src_h = Source("hot source")
    snk_h = Sink("hot sink")
    src_c = Source("cold source")
    snk_c = Sink("cold sink")
    hx = MovingBoundaryHeatExchanger("hx")

    c1 = Connection(src_h, "out1", hx, "in1")
    c2 = Connection(hx,    "out1", snk_h, "in1")
    c3 = Connection(src_c, "out1", hx, "in2")
    c4 = Connection(hx,    "out2", snk_c, "in1")
    nw.add_conns(c1, c2, c3, c4)

    c1.set_attr(fluid={fluid_h: 1}, T=T_h_in, p=p_h_in, m=mdot_h)
    c3.set_attr(fluid={fluid_c: 1}, T=T_c_in, p=p_c_in, m=mdot_c)

    hx.set_attr(
        pr1=1, pr2=1,
        alpha1_l=alpha_sp, alpha1_tp=alpha_tp, alpha1_g=alpha_sp, alpha1_sc=alpha_sp,
        alpha2_l=alpha_sp, alpha2_tp=alpha_tp, alpha2_g=alpha_sp, alpha2_sc=alpha_sp,
        area_ratio=1.0, R_cond=R_cond,
    )

    # Step 1: warm start — pin Q to the oracle value (area_zones inactive).
    hx.set_attr(Q=-Q_oracle)
    nw.solve("design")

    # Step 2: release Q, activate area_zones via area_hot, resolve.
    hx.set_attr(Q=None, area_hot=area)
    nw.solve("design", oscillation_damping=True)

    return nw, hx, c2, c4


# ── Section A: Water (hot) / n-Propane (cold), varying area ──────────────
_inp_A  = _ORACLE["section_A"]["inputs"]
_cases_A = _ORACLE["section_A"]["cases"]


@pytest.mark.parametrize(
    "case", _cases_A, ids=[f"A={c['A_m2']}_m2" for c in _cases_A]
)
def test_section_A_water_propane(case):
    _, hx, c2, c4 = _build_and_solve(
        _inp_A["Fluid_h"], _inp_A["T_h_in_K"], _inp_A["p_h_in_Pa"], _inp_A["mdot_h_kg_s"],
        _inp_A["Fluid_c"], _inp_A["T_c_in_K"], _inp_A["p_c_in_Pa"], _inp_A["mdot_c_kg_s"],
        _inp_A["alpha_liquid_vapor_W_m2_K"], _inp_A["alpha_two_phase_W_m2_K"],
        case["A_m2"], _inp_A["R_cond"],
        case["Q_W"],
    )
    assert abs(hx.Q.val_SI) == pytest.approx(case["Q_W"], rel=Q_RTOL)
    assert c2.T.val_SI == pytest.approx(case["T_h_out_K"], abs=T_ATOL_K)
    assert c4.T.val_SI == pytest.approx(case["T_c_out_K"], abs=T_ATOL_K)


# ── Section B: n-Propane (hot) / n-Propane (cold), 3 cases ───────────────
_inp_B  = _ORACLE["section_B"]["inputs"]
_cases_B = _ORACLE["section_B"]["cases"]


@pytest.mark.parametrize(
    "case", _cases_B, ids=[c["label"] for c in _cases_B]
)
def test_section_B_propane_propane(case):
    _, hx, c2, c4 = _build_and_solve(
        _inp_B["Fluid_h"], _inp_B["T_h_in_K"], _inp_B["p_h_in_Pa"], case["mdot_h_kg_s"],
        _inp_B["Fluid_c"], _inp_B["T_c_in_K"], _inp_B["p_c_in_Pa"], case["mdot_c_kg_s"],
        _inp_B["alpha_liquid_vapor_W_m2_K"], _inp_B["alpha_two_phase_W_m2_K"],
        case["A_m2"], _inp_B["R_cond"],
        case["Q_W"],
    )
    assert abs(hx.Q.val_SI) == pytest.approx(case["Q_W"], rel=Q_RTOL)
    assert c2.T.val_SI == pytest.approx(case["T_h_out_K"], abs=T_ATOL_K)
    assert c4.T.val_SI == pytest.approx(case["T_c_out_K"], abs=T_ATOL_K)
