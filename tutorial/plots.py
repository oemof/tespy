# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# color range
colors = ['#00395b', '#74adc1', '#b54036', '#ec6707',
          '#bfbfbf', '#999999', '#010101']


# %% figure 1: plot component exergy destruction

# color range for E_F and E_P bars
E_F_colors = ['#3a9dce', '#f08e2b', '#f08e2b', '#f08e2b', '#f08e2b', '#db5252']

# read data
df_ED_NH3 = pd.read_csv('NH3_E_D.csv', index_col=0)
df_ED_R410A = pd.read_csv('R410A_E_D.csv', index_col=0)

# get data from dataframes
y_NH3 = df_ED_NH3.columns.values.tolist()
y_NH3_pos = np.arange(len(y_NH3))
E_P_NH3 = df_ED_NH3.loc["E_P"]
E_D_NH3 = df_ED_NH3.loc["E_D"]
y_R410A = df_ED_R410A.columns.values.tolist()
y_R410A_pos = np.arange(len(y_R410A))
E_P_R410A = df_ED_R410A.loc["E_P"]
E_D_R410A = df_ED_R410A.loc["E_D"]

# create bar diagram of absolute exergy destruction
fig, axs = plt.subplots(2, 2, constrained_layout=True,
                        sharex='row', sharey='row')
axs[0, 0].barh(y_NH3_pos, E_P_NH3, align='center', color=E_F_colors)
axs[0, 0].barh(y_NH3_pos, E_D_NH3, align='center', left=E_P_NH3, label='E_D',
               color='#6ed880')
axs[0, 0].set_xlabel('Exergy in W')
axs[0, 0].set_yticks(y_NH3_pos)
axs[0, 0].set_yticklabels(y_NH3)
axs[0, 0].invert_yaxis()
axs[0, 0].set_title('NH3')

axs[0, 1].barh(y_R410A_pos, E_P_R410A, align='center', color=E_F_colors)
axs[0, 1].barh(y_R410A_pos, E_D_R410A, align='center', left=E_P_R410A,
               color='#6ed880')
axs[0, 1].set_xlabel('Exergy in W')
axs[0, 1].set_yticks(y_R410A_pos)
axs[0, 1].set_yticklabels(y_R410A)
axs[0, 1].set_title('R410A')
axs[0, 1].set_xlim(right=1000)

# create bar diagram of percentage exergy destruction
E_F_NH3 = df_ED_NH3.loc["E_P"][0]
E_F_R410A = df_ED_R410A.loc["E_P"][0]


axs[1, 0].barh(y_NH3_pos, E_P_NH3/E_F_NH3, align='center', color=E_F_colors)
axs[1, 0].barh(y_NH3_pos, E_D_NH3/E_F_NH3, align='center',
               left=E_P_NH3 / E_F_NH3, color='#6ed880')
axs[1, 0].set_xlabel('$\epsilon$')
axs[1, 0].set_yticks(y_NH3_pos)
axs[1, 0].set_yticklabels(y_NH3)
axs[1, 0].invert_yaxis()

axs[1, 1].barh(y_R410A_pos, E_P_R410A/E_F_R410A, align='center',
               color=E_F_colors)
axs[1, 1].barh(y_R410A_pos, E_D_R410A/E_F_R410A, align='center',
               left=E_P_R410A / E_F_R410A, color='#6ed880')
axs[1, 1].set_xlabel('$\epsilon$')
axs[1, 1].set_yticks(y_R410A_pos)
axs[1, 1].set_yticklabels(y_R410A)
axs[1, 1].set_xlim(right=1)

fig.suptitle('Component Exergy Destruction', fontsize=14)
fig.legend(bbox_to_anchor=(0.08, 0, 0.9, .0), loc='center',
           ncol=3, borderaxespad=0.)
plt.show()
fig.savefig('diagram_E_D.svg', bbox_inches='tight')


# %% figure 2: epsilon depending on ambient Temperature Tamb
#              and mean geothermal temperature Tgeo

Tamb_design = 2.8
Tgeo_design = 9.5

Tamb_range = [1, 4, 8, 12, 16, 20]
Tgeo_range = [11.5, 10.5, 9.5, 8.5, 7.5, 6.5]

# read data
df_eps_Tamb_NH3 = pd.read_csv('NH3_eps_Tamb.csv', index_col=0)
df_eps_Tgeo_NH3 = pd.read_csv('NH3_eps_Tgeo.csv', index_col=0)
df_eps_Tamb_R410A = pd.read_csv('R410A_eps_Tamb.csv', index_col=0)
df_eps_Tgeo_R410A = pd.read_csv('R410A_eps_Tgeo.csv', index_col=0)

# create plot
fig, axs = plt.subplots(1, 2, constrained_layout=True,
                        sharex='col', sharey='all')
axs[0].plot(Tamb_range, df_eps_Tamb_NH3.loc[Tgeo_design], 'x',
            label='NH3', color=colors[0], markersize=7)
axs[0].plot(Tamb_range, df_eps_Tamb_R410A.loc[Tgeo_design], 'x',
            label='R410A', color=colors[3], markersize=7)
axs[0].set_title('ambient Temperature')
axs[0].set_ylabel('$\epsilon$')
axs[0].set_ylim([0.25, 0.6])
axs[0].set_xlabel('$T_{amb}$ in °C ($T_{geo}$ = ' + str(Tgeo_design) + '°C)')
axs[0].set_ylabel('$\epsilon$')
axs[0].legend(loc='lower left')

axs[1].plot(Tgeo_range, df_eps_Tgeo_NH3.loc[Tamb_design], 'x',
            label='NH3', color=colors[0], markersize=7)
axs[1].plot(Tgeo_range, df_eps_Tgeo_R410A.loc[Tamb_design], 'x',
            label='R410A', color=colors[3], markersize=7)
axs[1].set_title('geothermal Temperature')
axs[1].set_xlabel('$T_{geo}$ in °C ($T_{amb}$ = 2.8°C)')
axs[1].legend(loc='lower left')

fig.suptitle('Exergetic Efficency depending on Temperature', fontsize=14)
fig.savefig('diagram_eps_Tamb_Tgeo.svg', bbox_inches='tight')


# %% figure 3: epsilon and COP depending on mean geothermal temperature Tgeo
#              and mean heating system temperature Ths

Ths_range = [42.5, 37.5, 32.5]
Tgeo_range = [10.5, 8.5, 6.5]

# read data
df_cop_Tgeo_Ths_NH3 = pd.read_csv('NH3_cop_Tgeo_Ths.csv', index_col=0)
df_eps_Tgeo_Ths_NH3 = pd.read_csv('NH3_eps_Tgeo_Ths.csv', index_col=0)
df_cop_Tgeo_Ths_R410A = pd.read_csv('R410A_cop_Tgeo_Ths.csv', index_col=0)
df_eps_Tgeo_Ths_R410A = pd.read_csv('R410A_eps_Tgeo_Ths.csv', index_col=0)

# create plot
fig, axs = plt.subplots(2, 2, constrained_layout=True,
                        sharex='all', sharey='row')
i = 0
for Tgeo in Tgeo_range:
    axs[0, 0].plot(Ths_range, df_eps_Tgeo_Ths_NH3.loc[Tgeo], 'x',
                   color=colors[i], label='$T_{geo}$ = ' + str(Tgeo) +
                   ' °C', markersize=7, linewidth=2)
    axs[1, 0].plot(Ths_range, df_cop_Tgeo_Ths_NH3.loc[Tgeo], 'x',
                   color=colors[i], markersize=7, linewidth=2)
    axs[0, 1].plot(Ths_range, df_eps_Tgeo_Ths_R410A.loc[Tgeo], 'x',
                   color=colors[i], markersize=7, linewidth=2)
    axs[1, 1].plot(Ths_range, df_cop_Tgeo_Ths_R410A.loc[Tgeo], 'x',
                   color=colors[i], markersize=7, linewidth=2)
    i += 1

axs[0, 0].set_ylabel('$\epsilon$')
axs[0, 0].set_ylim([0.42, 0.55])
axs[0, 0].set_title('NH3')

axs[1, 0].set_ylabel('COP')
axs[1, 0].set_xlabel('$T_{heating system}$ in °C')
axs[1, 0].set_ylim([3.5, 6.5])

axs[0, 1].set_title('R410A')

axs[1, 1].set_xlabel('$T_{heating system}$ in °C')

fig.legend(bbox_to_anchor=(0.08, -0.1, 0.9, .0), loc='lower left',
           ncol=3, mode="expand", borderaxespad=0.)
fig.suptitle(
    'COP and $\epsilon$ depending heating and geothermal mean temperature',
    fontsize=12)
fig.savefig('diagram_cop_eps_Tgeo_Ths.svg', bbox_inches='tight')


# %% figure 4: epsilon and COP depending on mean geothermal temperature Tgeo
#              and heating load Q_cond

Q_range = np.array([4.3e3, 4e3, 3.7e3, 3.4e3, 3.1e3, 2.8e3])
Tgeo_range = [10.5, 8.5, 6.5]

# read data
df_cop_Tgeo_Q_NH3 = pd.read_csv('NH3_cop_Tgeo_Q.csv', index_col=0)
df_eps_Tgeo_Q_NH3 = pd.read_csv('NH3_eps_Tgeo_Q.csv', index_col=0)
df_cop_Tgeo_Q_R410A = pd.read_csv('R410A_cop_Tgeo_Q.csv', index_col=0)
df_eps_Tgeo_Q_R410A = pd.read_csv('R410A_eps_Tgeo_Q.csv', index_col=0)

# create plot
fig, axs = plt.subplots(2, 2, constrained_layout=True,
                        sharex='all', sharey='row')

i = 0
for Tgeo in Tgeo_range:
    axs[0, 0].plot(Q_range, df_eps_Tgeo_Q_NH3.loc[Tgeo], 'x',
                   color=colors[i], label='$T_{geo}$ = ' + str(Tgeo) +
                   ' °C', markersize=7, linewidth=2)
    axs[1, 0].plot(Q_range, df_cop_Tgeo_Q_NH3.loc[Tgeo], 'x',
                   color=colors[i], markersize=7, linewidth=2)
    axs[0, 1].plot(Q_range, df_eps_Tgeo_Q_R410A.loc[Tgeo], 'x',
                   color=colors[i], markersize=7, linewidth=2)
    axs[1, 1].plot(Q_range, df_cop_Tgeo_Q_R410A.loc[Tgeo], 'x',
                   color=colors[i], markersize=7, linewidth=2)
    i += 1

axs[0, 0].set_ylabel('$\epsilon$')
axs[0, 0].set_ylim([0.4, 0.55])
axs[0, 0].set_title('NH3')

axs[1, 0].set_ylabel('COP')
axs[1, 0].set_xlabel(r'$Q_{cond}$ in kW')
axs[1, 0].set_ylim([3, 6])
axs[1, 0].set_xlim([2500, 4500])

axs[0, 1].set_title('R410A')

axs[1, 1].set_xlabel(r'$Q_{cond}$ in kW')

plt.xticks(np.arange(2500, 5000, step=1000), np.arange(2.5, 5, step=1))
fig.legend(bbox_to_anchor=(0.08, -0.1, 0.9, .0), loc='lower left',
           ncol=3, mode="expand", borderaxespad=0.)
fig.suptitle(
    'COP and $\epsilon$ depending on load and geothermal mean temperature',
    fontsize=12)
fig.savefig('diagram_cop_eps_Tgeo_Q.svg', bbox_inches='tight')
