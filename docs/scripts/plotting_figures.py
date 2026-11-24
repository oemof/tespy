"""Create the figures of the plotting page of the documentation.

Run from the repository root:

    python docs/scripts/plotting_figures.py
"""
import os

import matplotlib

matplotlib.use("Agg")

from matplotlib import pyplot as plt  # noqa: E402

from tespy.components import Compressor  # noqa: E402
from tespy.components import CycleCloser  # noqa: E402
from tespy.components import MovingBoundaryHeatExchanger  # noqa: E402
from tespy.components import SimpleHeatExchanger  # noqa: E402
from tespy.components import Sink  # noqa: E402
from tespy.components import Source  # noqa: E402
from tespy.components import Valve  # noqa: E402
from tespy.connections import Connection  # noqa: E402
from tespy.networks import Network  # noqa: E402
from tespy.plotting import plot_characteristic  # noqa: E402
from tespy.plotting import plot_property_diagram  # noqa: E402
from tespy.plotting import plot_QT_diagram  # noqa: E402
from tespy.plotting import plot_sweep  # noqa: E402
from tespy.tools.characteristics import CharMap  # noqa: E402

PATH = os.path.join("docs", "_static", "images", "modules", "plotting")


def heat_pump():
    nw = Network(iterinfo=False)
    nw.units.set_defaults(
        temperature="°C", pressure="bar", pressure_difference="bar",
        heat="kW", power="kW", enthalpy="kJ/kg", entropy="kJ/kg/K"
    )
    cp = Compressor("compressor")
    cc = CycleCloser("cycle_closer")
    cd = MovingBoundaryHeatExchanger("condenser")
    va = Valve("expansion valve")
    ev = SimpleHeatExchanger("evaporator")
    so = Source("water source")
    si = Sink("water sink")
    c1 = Connection(cc, "out1", cd, "in1", label="c1")
    c2 = Connection(cd, "out1", va, "in1", label="c2")
    c3 = Connection(va, "out1", ev, "in1", label="c3")
    c4 = Connection(ev, "out1", cp, "in1", label="c4")
    c5 = Connection(cp, "out1", cc, "in1", label="c5")
    a1 = Connection(so, "out1", cd, "in2", label="a1")
    a2 = Connection(cd, "out2", si, "in1", label="a2")
    nw.add_conns(c1, c2, c3, c4, c5, a1, a2)
    cd.set_attr(dp1=0, dp2=0, Q=-1e3, td_pinch=5)
    ev.set_attr(dp=0)
    cp.set_attr(eta_s=0.8)
    c1.set_attr(fluid={"R290": 1})
    c2.set_attr(td_bubble=5)
    c4.set_attr(td_dew=5, T=15)
    a1.set_attr(fluid={"water": 1}, p=1, T=50)
    a2.set_attr(T=65)
    nw.solve("design")
    nw.assert_convergence()
    return nw


def save(ax, name, variant):
    ax.figure.savefig(os.path.join(PATH, f"{name}_{variant}.png"), dpi=110)
    plt.close(ax.figure)


def main():
    os.makedirs(PATH, exist_ok=True)
    nw = heat_pump()
    char_map = CharMap(
        x=[0.9, 1.0, 1.1],
        y=[[0.5, 0.8, 1.0], [0.6, 0.9, 1.1], [0.7, 1.0, 1.2]],
        z=[[0.8, 0.84, 0.85], [0.82, 0.87, 0.88], [0.84, 0.89, 0.9]],
    )
    sweep = {
        "T_cond": [35, 45, 55, 65] * 2,
        "COP": [6.0, 4.6, 3.7, 3.0, 6.1, 4.8, 3.9, 3.2],
        "fluid": ["R290"] * 4 + ["R600a"] * 4,
    }
    for variant, style in (("light", "tespy"), ("dark", "tespy-dark")):
        save(plot_property_diagram(nw, "c1", style=style), "logph", variant)
        save(plot_QT_diagram(nw.get_comp("condenser"), style=style), "QT", variant)
        save(plot_characteristic(
            char_map, operating_points=[(0.95, 0.85)], style=style,
            title="Compressor map"
        ), "char_map", variant)
        save(plot_sweep(
            sweep, "T_cond", "COP", group="fluid", units={"T_cond": "°C"},
            style=style
        ), "sweep", variant)
    save(plot_property_diagram(nw, "c1", style="tespy-paper"), "logph", "paper")


if __name__ == "__main__":
    main()
