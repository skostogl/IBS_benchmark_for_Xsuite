import pandas as pd
import matplotlib.pylab as plt
import yaml
import numpy as np
import tree_maker

# Load simulation parameters
with open("config.yaml", "r") as f:
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

# Activate the different modes
modes = []
if config['mode_kinetic']: 
    modes.append('kinetic')
if config['mode_analytical']: 
    modes.append('analytical')
if config['mode_simple']: 
    modes.append('simple')

# Iterate over the different modes 
for mode in modes:
    if mode == 'analytical':
        pd_analytical = pd.read_parquet(f"{save_to}/xsuite_{mode}.parquet")
    elif mode == 'kinetic':
        pd_kinetic = pd.read_parquet(f"{save_to}/xsuite_{mode}.parquet")
    elif mode == 'simple':
        pd_simple = pd.read_parquet(f"{save_to}/xsuite_{mode}.parquet")

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

f, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize = (16,5))
f.suptitle('SPS PB ion tracking: $N_{{b}}$ = {:.2e}, $\sigma_{{z}}$ = {:.3f} m'.format(bunch_intensity, sigma_z), fontsize=20)

# ax1.plot(nag[0], 'r')
plt.sca(ax1)
if 'kinetic' in modes:
    plt.plot(pd_kinetic['epsn_x'].values, alpha=0.7, label='Kinetic run')
if 'analytical' in modes:
    plt.plot(pd_analytical['epsn_x'].values, c='k', label='Analytical')
if 'simple' in modes:
    plt.plot(pd_simple['epsn_x'].values, c='r--', label='Simple')
plt.legend(fontsize=12)

plt.sca(ax2)
if 'kinetic' in modes:
    plt.plot(pd_kinetic['epsn_y'].values, alpha=0.7, label='Kinetic run')
if 'analytical' in modes:
    plt.plot(pd_analytical['epsn_y'].values, c='k', label='Analytical')
if 'simple' in modes:
    plt.plot(pd_simple['epsn_y'].values, c='r--', label='Simple')

plt.sca(ax3)
if 'kinetic' in modes:
    plt.plot(1e3*pd_kinetic['sig_delta'].values, alpha=0.7, label='Kinetic run')
if 'analytical' in modes:
    plt.plot(1e3*pd_analytical['sig_delta'].values, c='k', label='Analytical')
if 'simple' in modes:
    plt.plot(1e3*pd_simple['sig_delta'].values, c='r--', label='Simple')

ax1.set_ylabel(r'$\varepsilon_{x,n}$ [m]')
ax1.set_xlabel('Turns')

ax2.set_ylabel(r'$\varepsilon_{y,n}$ [m]')
ax2.set_xlabel('Turns')

ax3.set_ylabel(r'$\sigma_{\delta}$ [$10^{-3}$]')
ax3.set_xlabel('Turns')

plt.tight_layout()
f.savefig("{}/Emittances_SPS_PB_IBS_tracking_{}_turns.png".format(save_to, n_turns), dpi=250)
    
# Plot the integral evolution
fig, (ax11, ax22, ax33) = plt.subplots(1, 3, figsize = (16,5))
fig.suptitle('SPS PB ion tracking. IBS: Growth Rates')
if 'kinetic' in modes:
    ax11.plot(pd_kinetic["kinTx"][1:], marker='o', linestyle=None, markerfacecolor='none', label='Kinetic: Tx (last run)')
if 'analytical' in modes:
    ax11.plot(pd_analytical['Ixx'][1:], c='b', linewidth=4, label='Analytical: Ixx')
if 'simple' in modes:
    ax11.plot(pd_simple['Ixx'][1:], c='b', linewidth=4, label='Simple: Ixx')
ax11.set_xlabel('Turns')
ax11.legend(fontsize=12)

if 'kinetic' in modes:
    ax22.plot(pd_kinetic["kinTy"][1:], marker='o', markersize=8, linestyle=None, markerfacecolor='none', label='Kinetic: Ty (last run)')
if 'analytical' in modes:
    ax22.plot(pd_analytical['Iyy'][1:], c='g', linewidth=4, label='Analytical: Iyy')
if 'simple' in modes:
    ax22.plot(pd_simple['Iyy'][1:], c='b', linewidth=4, label='Simple: Iyy')
ax22.set_xlabel('Turns')
ax22.legend(fontsize=12)

if 'kinetic' in modes:
    ax33.plot(pd_kinetic["kinTz"][1:], marker='o', markersize=8,  linestyle=None, markerfacecolor='none', label='Kinetic: Tz (last run)')
if 'analytical' in modes:
    ax33.plot(pd_analytical['Ipp'][1:], c='k', linewidth=4, label='Ipp')
if 'simple' in modes:
    ax33.plot(pd_simple['Ipp'][1:], c='b', linewidth=4, label='Simple: Ipp')
ax33.set_xlabel('Turns')
ax33.legend(fontsize=12)

fig.savefig("{}/IBS_Growth_SPS_Pb_IBS_tracking_{}_turns.png".format(save_to, n_turns), dpi=250)

# Plot the bunch length evolution 
fig2, ax = plt.subplots(1, 1, figsize = (10,5))
fig2.suptitle('SPS PB ion tracking. IBS: Bunch lengths')
if 'analytical' in modes:
    ax.plot(pd_analytical['bl'][1:], c='r', label='Analytical Ion Bunch length')
if 'kinetic' in modes:
    ax.plot(pd_kinetic['bl'][1:], c='b', label='Kinetic (from particles object)')
if 'simple' in modes:
    ax.plot(pd_simple['bl'][1:], c='b', label='Simple (from particles object)')
ax.set_xlabel('Turns')
ax.set_ylabel('Bunch Length')
ax.legend()
fig2.savefig("{}/BunchLength_SPS_Pb_IBS_tracking_{}_turns.png".format(save_to, n_turns), dpi=250)

# Log as completed
tree_maker.tag_json.tag_it(config['log_file'], 'completed')  
