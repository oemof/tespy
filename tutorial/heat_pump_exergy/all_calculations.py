# -*- coding: utf-8 -*-
from itertools import product
import pandas as pd
import matplotlib.pyplot as plt

from heat_pump_model import HeatPumpModel

# %% ambient conditions and design parameters

pamb = 1.013  # bar
Tamb = 1.0    # °C
Tgeo = 10    # °C  (mean geothermal temperature)

# %% set up model and run design calculation
for fluid in ["NH3", "R290"]:

    hp = HeatPumpModel(fluid, Tamb=Tamb, pamb=pamb, Tgeo=Tgeo)
    hp.nw.print_results()
    hp.save_design()

    # %% log(p)-h diagram

    fig, ax = hp.plot_logph_diagram_matplotlib("c1")
    plt.tight_layout()
    fig.savefig(f"{fluid}_logph.svg")
    plt.close()

    # %% exergy analysis at design point

    E_F = {"inputs": ["e1", "c11"], "outputs": ["c13"]}
    E_P = {"inputs": ["c23"], "outputs": ["c21"]}
    exergy_kwargs = {"Tamb": Tamb, "pamb": pamb, "E_F": E_F, "E_P": E_P}
    run_exergy = lambda model: model.run_exergy_analysis(**exergy_kwargs)

    ean = hp.run_exergy_analysis(**exergy_kwargs)
    print("\n##### EXERGY ANALYSIS #####\n")
    ean.exergy_results()

    ean.plot_exergy_waterfall(title=f"{fluid} Heat Pump Exergy Analysis", show_plot=False)
    plt.tight_layout()
    plt.savefig(f"{fluid}_waterfall.svg")
    plt.close()

    # %%[ed_export] exergy destruction waterfall data export
    df_comp, _, _ = ean.exergy_results(print_results=False)
    df_comp = df_comp[(df_comp["Component"] != "TOT") & df_comp["E_F [kW]"].notna()].copy()

    comps = ["E_F"]
    E_D_list = [0]
    running = ean.E_F          # W
    E_P_list = [running]

    for _, row in df_comp.iterrows():
        e_d = row["E_D [kW]"] * 1e3   # W
        if e_d > 1:
            comps.append(row["Component"])
            E_D_list.append(e_d)
            running -= e_d
            E_P_list.append(running)

    comps.append("E_P")
    E_D_list.append(0)
    E_P_list.append(running)

    pd.DataFrame([E_D_list, E_P_list], columns=comps, index=["E_D", "E_P"]).to_csv(f"{fluid}_E_D.csv")

    # %%[parametric] parametric studies

    hp.nw.iterinfo = False

    Tamb_design = Tamb
    Tgeo_design = Tgeo

    # --- varying ambient temperature (exergy only, no re-solve) ---

    Tamb_range = [4, 8, 12, 16, 20]
    eps_Tamb = [
        hp.run_exergy_analysis(T, pamb, E_F, E_P).epsilon
        for T in Tamb_range
    ]
    pd.DataFrame([eps_Tamb], columns=Tamb_range, index=[Tgeo_design]).to_csv(f"{fluid}_eps_Tamb.csv")

    # --- varying mean geothermal temperature (offdesign) ---

    Tgeo_range = [14, 13, 12, 11, 10, 9, 8]
    results_Tgeo = hp.sensitivity_analysis(
        param_dict={
            "T_geo":  [T for T in Tgeo_range]
        },
        result_param_list=["epsilon"],
        mode="offdesign",
        postproc_func=run_exergy,
    )
    pd.DataFrame(
        [results_Tgeo["epsilon"].values],
        columns=Tgeo_range, index=[Tamb_design],
    ).to_csv(f"{fluid}_eps_Tgeo.csv")

    # reset geo temperatures to design values
    hp.set_parameters(T_geo=Tgeo_design)

    # %%[tgeo_ths] varying geothermal and heating system temperatures

    Tgeo_range = [14, 12, 10]
    Ths_range = [45, 40, 35]

    T_geo_list = []
    T_hs_list = []

    for tgeo, ths in product(Tgeo_range, Ths_range):
        T_geo_list.append(tgeo)
        T_hs_list.append(ths)

    result = {
        "T_geo": T_geo_list,
        "T_hs": T_hs_list
    }

    results = hp.sensitivity_analysis(
        param_dict=result,
        result_param_list=["COP", "epsilon"],
        mode="offdesign",
        postproc_func=run_exergy,
    )

    results.pivot(index="T_geo", columns="T_hs", values="COP").to_csv(f"{fluid}_cop_Tgeo_Ths.csv")
    results.pivot(index="T_geo", columns="T_hs", values="epsilon").to_csv(f"{fluid}_eps_Tgeo_Ths.csv")

    # reset heating system temperatures to design values
    hp.set_parameters(T_hs=40)

    # %%[tgeo_q] varying geothermal temperature and heating load

    Q_range = [4.3e3, 4e3, 3.7e3, 3.4e3, 3.1e3, 2.8e3]

    T_geo_list = []
    Q_list = []

    for tgeo, q in product(Tgeo_range, Q_range):
        T_geo_list.append(tgeo)
        Q_list.append(-q)

    result = {
        "T_geo": T_geo_list,
        "Q": Q_list
    }

    results = hp.sensitivity_analysis(
        param_dict=result,
        result_param_list=["COP", "epsilon"],
        mode="offdesign",
        postproc_func=run_exergy,
    )

    results.pivot(index="T_geo", columns="Q", values="COP").to_csv(f"{fluid}_cop_Tgeo_Q.csv")
    results.pivot(index="T_geo", columns="Q", values="epsilon").to_csv(f"{fluid}_eps_Tgeo_Q.csv")
