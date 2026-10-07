# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt

from heat_pump_model import HeatPumpModel

# %% ambient conditions and design parameters

pamb = 1.013  # bar
Tamb = 2.8    # °C
Tgeo = 9.5    # °C  (mean geothermal temperature)

# %% set up model and run design calculation

hp = HeatPumpModel("NH3", Tamb=Tamb, pamb=pamb, Tgeo=Tgeo)
hp.solve_model_design()
hp.nw.assert_convergence()
hp.nw.print_results()
hp.save_design()

# %%[logph] log(p)-h diagram
fig, ax = hp.plot_logph_diagram_matplotlib("c1", save_dir=".")
plt.tight_layout()
fig.savefig("NH3_logph.svg")
plt.close()

# %%[exergy] exergy analysis
E_F = {"inputs": ["e1", "c11"], "outputs": ["c13"]}
E_P = {"inputs": ["c23"], "outputs": ["c21"]}

ean = hp.run_exergy_analysis(Tamb, pamb, E_F, E_P)
print("\n##### EXERGY ANALYSIS #####\n")
ean.exergy_results()
