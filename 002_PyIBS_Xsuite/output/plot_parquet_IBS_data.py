#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to plot parquet files generated from Xsuite_IBS_allmodels.py
"""
import pandas as pd
import matplotlib.pylab as plt
import yaml
import numpy as np

# Load simulation parameters
with open("../config.yaml", "r") as f:
    config = yaml.safe_load(f)

energy=config['energy']
nemitt_x = config["emit_x"]*1e-6
nemitt_y = config["emit_y"]*1e-6
RF_Voltage = config["V0max"]
Harmonic_Num = config["h"]
bunch_intensity= float(config["bunch_intensity"])
sigma_z = config['sigma_z']
n_part = int(config['n_part'])
n_turns = int(config['n_turns'])
save_to = config['save_to']
nr_kinetic_runs = config['nr_kinetic_runs']  # how many kinetic runs to average over

destination_folder = 'Plots'

# Load dataframes - kinetic
dfs_kinetic = []
if nr_kinetic_runs == 1:
    plot_average = False
    j = 0
    df_kinetic = pd.read_parquet("xsuite_run{}_kinetic_{}.parquet".format(j, n_turns))
    dfs_kinetic.append(pd.read_parquet("xsuite_run{}_kinetic_{}.parquet".format(j, n_turns)))
elif nr_kinetic_runs ==3: 
    plot_average = True
    for j in range(nr_kinetic_runs):
        dfs_kinetic.append(pd.read_parquet("xsuite_run{}_kinetic_{}.parquet".format(j, n_turns)))
        
# Load dataframe - analytical
df_analytical = pd.read_parquet("xsuite_analytical_{}.parquet".format(n_turns))

# ------------------------ DEFINE PLOT PARAMETERS --------------------------------------
SMALL_SIZE = 20
MEDIUM_SIZE = 20
BIGGER_SIZE = 20
plt.rcParams["font.family"] = "serif"
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)    # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)   # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)   # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)   # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


# ----------------------- PLOT THE FIGURE ---------------------------------------------
fl_kinetic=True
fl_analytical=True

if fl_analytical:
    analytical = df_analytical

f, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize = (16,5))
f.suptitle('SPS PB ion tracking: $N_{{b}}$ = {:.2e}, $\sigma_{{z}}$ = {:.3f}'.format(bunch_intensity, sigma_z), fontsize=20)

# Find the average kinetic kick
if nr_kinetic_runs == 1:
    eps_x_mean = np.array(df_kinetic['eps_x'])
    eps_y_mean = np.array(df_kinetic['eps_y'])
    sig_delta_mean = np.array(df_kinetic['sig_delta'])
elif nr_kinetic_runs ==3: 
    All_epsx = np.array([dfs_kinetic[0]['eps_x'], dfs_kinetic[1]['eps_x'], dfs_kinetic[2]['eps_x']])
    All_epsy = np.array([dfs_kinetic[0]['eps_y'], dfs_kinetic[1]['eps_y'], dfs_kinetic[2]['eps_y']])
    All_sig_delta = np.array([dfs_kinetic[0]['sig_delta'], dfs_kinetic[1]['sig_delta'], dfs_kinetic[2]['sig_delta']])
    All_BL = np.array([dfs_kinetic[0]['bl'], dfs_kinetic[1]['bl'], dfs_kinetic[2]['bl']])
    eps_x_mean = np.mean(All_epsx, axis=0)
    eps_y_mean = np.mean(All_epsy, axis=0)
    sig_delta_mean = np.mean(All_sig_delta, axis=0)
    bl_mean = np.mean(All_BL, axis=0)
else:
    print("Set nr of kinetic runs to 3")

# ax1.plot(nag[0], 'r')
plt.sca(ax1)
if fl_kinetic:
    if plot_average:
        plt.plot(eps_x_mean, alpha=0.7, label='Mean Kinetic {} runs'.format(nr_kinetic_runs), c='b')
    else:
        for k, kinetic in enumerate(dfs_kinetic):
            plt.plot(kinetic['eps_x'].values, alpha=0.7, label=f'Kinetic run {k+1}')
if fl_analytical:
    plt.plot(analytical['eps_x'].values, c='k', label='Analytical')
plt.legend(fontsize=12)

plt.sca(ax2)
if fl_kinetic:
    if plot_average:
        plt.plot(eps_y_mean, alpha=0.7, c='b')
    else:
        for k, kinetic in enumerate(dfs_kinetic):
            plt.plot(kinetic['eps_y'].values, alpha=0.7)
if fl_analytical:
    plt.plot(analytical['eps_y'].values, c='k')

plt.sca(ax3)
if fl_kinetic:
    if plot_average:
        plt.plot(sig_delta_mean*1e3, alpha=0.7, c='b')
    else:
        for k, kinetic in enumerate(dfs_kinetic):
            plt.plot(kinetic['sig_delta'].values*1e3, alpha=0.7)
if fl_analytical:
    plt.plot(analytical['sig_delta'].values*1e3, c='k')

ax1.set_ylabel(r'$\varepsilon_x$ [m]')
ax1.set_xlabel('Turns')

ax2.set_ylabel(r'$\varepsilon_y$ [m]')
ax2.set_xlabel('Turns')

ax3.set_ylabel(r'$\sigma_{\delta}$ [$10^{-3}$]')
ax3.set_xlabel('Turns')

plt.tight_layout()
plt.show()
f.savefig("{}/Emittances_SPS_PB_IBS_tracking_{}_turns.png".format(destination_folder, n_turns), dpi=250)
    
# Plot the integral evolution
fig, (ax11, ax22, ax33) = plt.subplots(1, 3, figsize = (16,5))
fig.suptitle('SPS PB ion tracking. IBS: Growth Rates')

ax11.plot(df_kinetic["kinTx"][1:], marker='o', linestyle=None, markerfacecolor='none', label='Kinetic: Tx (last run)')
ax11.plot(df_analytical['Ixx'][1:], c='b', linewidth=4, label='Analytical: Ixx')
ax11.set_xlabel('Turns')
ax11.legend(fontsize=12)

ax22.plot(df_kinetic["kinTy"][1:], marker='o', markersize=8, linestyle=None, markerfacecolor='none', label='Kinetic: Ty (last run)')
ax22.plot(df_analytical['Iyy'][1:], c='g', linewidth=4, label='Analytical: Iyy')
ax22.set_xlabel('Turns')
ax22.legend(fontsize=12)

ax33.plot(df_kinetic["kinTz"][1:], marker='o', markersize=8,  linestyle=None, markerfacecolor='none', label='Kinetic: Tz (last run)')
ax33.plot(df_analytical['Ipp'][1:], c='k', linewidth=4, label='Ipp')
ax33.set_xlabel('Turns')
ax33.legend(fontsize=12)

fig.savefig("{}/IBS_Growth_SPS_Pb_IBS_tracking_{}_turns.png".format(destination_folder, n_turns), dpi=250)

# Plot the bunch length evolution 
fig2, ax = plt.subplots(1, 1, figsize = (10,5))
fig2.suptitle('SPS PB ion tracking. IBS: Bunch lengths')
ax.plot(df_analytical['bl_proton_version'][1:], c='g', marker='o', markersize=8,  linestyle=None, markerfacecolor='none', label='Analytical (proton version)')
ax.plot(df_analytical['bl'][1:], c='r', label='Analytical')
if plot_average:
    ax.plot(bl_mean[1:], c='b', label='Mean Kinetic')
else:
    ax.plot(df_kinetic['bl'][1:], c='b', label='Kinetic (from particles object)')
ax.set_xlabel('Turns')
ax.set_ylabel('Bunch Length')
ax.legend()
fig2.savefig("{}/BunchLength_SPS_Pb_IBS_tracking_{}_turns.png".format(destination_folder, n_turns), dpi=250)