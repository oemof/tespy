"""
Basic vapor-compression refrigeration cycle simulation using TeSPy.

Refrigerant: R134a
Cycle: Evaporator → Compressor → Condenser → Expansion Valve → Evaporator

State points:
  c0: CycleCloser → Evaporator  (two-phase mixture, low pressure)
  c1: Evaporator  → Compressor  (saturated vapor, -10°C)
  c2: Compressor  → Condenser   (superheated vapor, high pressure)
  c3: Condenser   → Valve       (saturated liquid, 40°C)
  c4: Valve       → CycleCloser (two-phase mixture, low pressure)
"""

# ── 1. Network ────────────────────────────────────────────────────────────────
from tespy.networks import Network

nw = Network()
nw.units.set_defaults(
    temperature="degC",
    pressure="bar",
    enthalpy="kJ/kg",
    heat="kW",
    power="kW",
)

# ── 2. Components ─────────────────────────────────────────────────────────────
from tespy.components import CycleCloser, Compressor, Valve, SimpleHeatExchanger

cc = CycleCloser("cycle closer")
ev = SimpleHeatExchanger("evaporator")
cp = Compressor("compressor")
co = SimpleHeatExchanger("condenser")
va = Valve("expansion valve")

# ── 3. Connections ────────────────────────────────────────────────────────────
from tespy.connections import Connection

c0 = Connection(cc, "out1", ev, "in1",  label="0")
c1 = Connection(ev, "out1", cp, "in1",  label="1")
c2 = Connection(cp, "out1", co, "in1",  label="2")
c3 = Connection(co, "out1", va, "in1",  label="3")
c4 = Connection(va, "out1", cc, "in1",  label="4")

nw.add_conns(c0, c1, c2, c3, c4)

# ── 4. Boundary conditions ────────────────────────────────────────────────────
# Component parameters
ev.set_attr(pr=0.98, Q=10)        # 10 kW cooling duty; slight pressure drop
co.set_attr(pr=0.98)              # slight pressure drop; heat determined by energy balance
cp.set_attr(eta_s=0.85)           # 85 % isentropic efficiency

# State-point specifications
c1.set_attr(
    fluid={"R134a": 1},           # pure R134a throughout the cycle
    T=-10,                        # evaporation temperature
    x=1,                          # saturated vapor at compressor inlet
)
c3.set_attr(
    T=40,                         # condensation temperature
    x=0,                          # saturated liquid at valve inlet
)

# ── 5. Solve design point ─────────────────────────────────────────────────────
nw.solve(mode="design")
nw.print_results()

# ── 6. Performance summary ────────────────────────────────────────────────────
Q_evap   = ev.Q.val          # refrigeration capacity  (positive = heat in)
W_comp   = cp.P.val          # compressor shaft power   (positive = work in)
Q_cond   = co.Q.val          # heat rejected at condenser (negative = heat out)

COP      = Q_evap / W_comp
T_evap_K = c1.T.val + 273.15
T_cond_K = c3.T.val + 273.15
COP_carnot = T_evap_K / (T_cond_K - T_evap_K)

print("\n" + "=" * 50)
print("  REFRIGERATION CYCLE — PERFORMANCE SUMMARY")
print("=" * 50)
print(f"  Refrigerant           : R134a")
print(f"  Evaporation temp      : {c1.T.val:.1f} °C  ({c1.p.val:.3f} bar)")
print(f"  Condensation temp     : {c3.T.val:.1f} °C  ({c3.p.val:.3f} bar)")
print(f"  Pressure ratio        : {c2.p.val / c1.p.val:.2f}")
print(f"  Mass flow rate        : {c1.m.val:.4f} kg/s")
print(f"  Refrigeration capacity: {Q_evap:.2f} kW")
print(f"  Compressor power      : {W_comp:.2f} kW")
print(f"  Heat rejected         : {abs(Q_cond):.2f} kW")
print(f"  COP (actual)          : {COP:.3f}")
print(f"  COP (Carnot)          : {COP_carnot:.3f}")
print(f"  Carnot efficiency     : {COP / COP_carnot * 100:.1f} %")
print("=" * 50)

# ── 7. Parametric study: COP vs evaporation temperature ──────────────────────
import numpy as np
import matplotlib.pyplot as plt

T_evap_range = np.linspace(-20, 5, 11)
COP_values   = []

for T_evap in T_evap_range:
    c1.set_attr(T=T_evap)
    nw.solve(mode="design")
    COP_values.append(ev.Q.val / cp.P.val)

# Reset to design point
c1.set_attr(T=-10)
nw.solve(mode="design")

# Plot
plt.rc("font", size=13)
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(T_evap_range, COP_values, marker="o", color="#1f567d", linewidth=2)
ax.set_xlabel("Evaporation temperature (°C)")
ax.set_ylabel("COP  (—)")
ax.set_title("R134a Refrigeration Cycle: COP vs Evaporation Temperature\n"
             "(Condensation temp = 40 °C, η_s = 0.85, Q_ref = 10 kW)")
ax.grid(True, alpha=0.4)
plt.tight_layout()
fig.savefig("refrigeration_COP_parametric.svg")
print("\nParametric plot saved → refrigeration_COP_parametric.svg")
