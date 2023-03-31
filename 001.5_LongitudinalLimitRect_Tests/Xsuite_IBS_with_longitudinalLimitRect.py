"""
Tester to check effect with or without LongitudinalLimitRect, to kill particles outside bucket 
"""

# IBS module from M. Zampetakis in https://github.com/MichZampetakis/IBS_for_Xsuite
import sys
import pathlib
import json
import numpy as np
import xobjects as xo
import xtrack as xt
import xpart as xp
import matplotlib.pyplot as plt
import pandas as pd
from cpymad.madx import Madx
from scipy.constants import e as qe
from scipy.constants import m_p, m_e
import yaml

# Load simulation parameters
with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

import sys
sys.path.append(config['ibs_lib_path'])
from lib.IBSfunctions import *

# Check whether we should use longitudinalLimitRect or not
longitudinalLimit = config['longitudinalLimit']

xsuite_line = config['xsuite_line']
energy=config['energy']
nemitt_x = config["emit_x"]*1e-6
nemitt_y = config["emit_y"]*1e-6
RF_Voltage = config["V0max"]
Harmonic_Num = config["h"]
bunch_intensity= float(config["bunch_intensity"])
sigma_z = config['sigma_z']
n_part = int(config['n_part'])
n_turns = int(config['n_turns'])
IBS_step = int(config['IBS_step'])
Energy_loss = 0
###########################################################

## Load xsuite line
with open(xsuite_line) as fid:
    dd=json.load(fid)
line = xt.Line.from_dict(dd)

p0 = line.particle_ref

# Add longitudinal limit rectangle 
bucket_length = line.get_length()/Harmonic_Num
if longitudinalLimit:
    mode_limit = 'LongitudinalLimitRect'
    line.unfreeze() # if you had already build the tracker
    line.append_element(element=xt.LongitudinalLimitRect(min_zeta=-bucket_length/2, max_zeta=bucket_length/2), name='long_limit')
    line.build_tracker()
else:
    mode_limit = 'no_LongitudinalLimitRect'

## Choose a context
context = xo.ContextCpu()         # For CPU

## Transfer lattice on context and compile tracking code
tracker = xt.Tracker(_context=context, line=line,
                    extra_headers=['#define XTRACK_MULTIPOLE_NO_SYNRAD'])

tracker.optimize_for_tracking()

particles0 = xp.generate_matched_gaussian_bunch(
         num_particles=n_part, total_intensity_particles=bunch_intensity,
         nemitt_x=nemitt_x, nemitt_y=nemitt_y, sigma_z=sigma_z,
         particle_ref=p0, tracker=tracker)

tw = tracker.twiss(particle_ref = p0)

# ----- Initialize IBS -----
IBS = NagaitsevIBS()
IBS.set_beam_parameters(particles0)
IBS.set_optic_functions(tw)

particles = particles0.copy()

## Initialize dictionaries for tbt checks 
columns = {'sig_x', 'sig_y', 'sig_delta', 'mean_delta', 'sig_z', 'BL_theoretical', 'NrParticles'}
tbt_checks = {nn: np.zeros((n_turns), dtype = float) for nn in columns}


# Track and investigate stability of bunch length
for i in range(n_turns):
    print('Turn = ', i)
    print('N_part = ',len(particles.x[particles.state > 0]))

    # Calculate RF bucket length here, and "kill" particles that are outside if longitudinal rect 
    tbt_checks['sig_x'][i] = np.std(particles.x[particles.state > 0])
    tbt_checks['sig_y'][i] = np.std(particles.y[particles.state > 0])
    tbt_checks['sig_delta'][i] = np.std(particles.delta[particles.state > 0])
    tbt_checks['mean_delta'][i] = np.mean(particles.delta[particles.state > 0])
    tbt_checks['sig_z'][i]  = np.std(particles.zeta[particles.state > 0])
    tbt_checks['NrParticles'][i]  = len(particles.zeta[particles.state > 0])
    
    # Check theoretical bunch length 
    Sigma_E = tbt_checks['sig_delta'][i]*IBS.betar**2
    BunchL = ion_BunchLength(IBS.Circu, Harmonic_Num, IBS.EnTot, IBS.slip, 
               Sigma_E, IBS.betar, RF_Voltage*1e-3, Energy_loss, IBS.Ncharg)
    tbt_checks['BL_theoretical'][i] = BunchL

    # Apply the kinetic kick from the calculated coefficient 
    if (i % IBS_step == 0) or (i==1):
       IBS.calculate_kinetic_coefficients(particles)
    IBS.apply_kinetic_kick(particles)
    print(f"Dx: {IBS.Dx}, Dy: {IBS.Dy}, Dz: {IBS.Dz}")
    print(f"Fx: {IBS.Fx}, Fy: {IBS.Fy}, Fz: {IBS.Fz}")

    tracker.track(particles)

# Then track one last turn without any IBS kick to "clean" remaining particles
tracker.track(particles)

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
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# Plot the figures
fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (12,7))
fig.suptitle('SPS Pb ions: Delta mean and spread', fontsize=18)
ax1.plot(tbt_checks['mean_delta'], 'b', label='Mean($\delta$)')
ax1.legend() 
ax2.plot(tbt_checks['sig_delta'], 'r', label='$\\sigma_{\delta}$')
ax2.legend() 
fig.savefig(f'Plots/{mode_limit}_SPS_Pb_ions_Nominal_parameter_check_DELTA_{n_turns}turns.png', dpi=250)

fig2, (ax3, ax4) = plt.subplots(2, 1, figsize = (10,7))
fig2.suptitle('SPS Pb ions - nominal parameters: Bunch Length', fontsize=18)
ax3.plot(tbt_checks['sig_z'], 'g', label='$\\sigma_{z}$ from tracking')
ax3.legend() 
ax4.plot(tbt_checks['BL_theoretical'], 'k', label='IBS theoretical Bunch Length')
ax4.legend() 
fig2.savefig(f'Plots/{mode_limit}_SPS_Pb_ions_BL_{n_turns}turns.png', dpi=250)

fig3, ax5 = plt.subplots(1, 1, figsize = (10,5))
fig3.suptitle('SPS PB ion tracking. Last turn Longitudinal Phase Space')
ax5.plot(particles.zeta[particles.state > 0], particles.delta[particles.state > 0]*1000, '.', markersize=3)
ax5.axvline(x=bucket_length/2, color='r', linestyle='dashed')
ax5.axvline(x=-bucket_length/2, color='r', linestyle='dashed')
ax5.set_xlabel(r'z [-]')
ax5.set_ylabel(r'$\delta$ [1e-3]')
fig3.savefig(f'Plots/{mode_limit}_SPS_Pb_ions_BUCKET_{n_turns}turns.png', dpi=250)

fig4, ax6 = plt.subplots(1, 1, figsize = (10,5))
fig4.suptitle('SPS PB ion tracking. Particles survival')
ax6.plot(tbt_checks['NrParticles'], 'g')
fig4.savefig(f'Plots/{mode_limit}_SPS_Pb_ions_SURVIVAL_{n_turns}turns.png', dpi=250)





