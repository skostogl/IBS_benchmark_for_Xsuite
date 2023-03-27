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

# Initiate values from config file 
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
save_to = config['save_to']
nr_kinetic_runs = config['nr_kinetic_runs']  # how many kinetic runs to average over

# Check that saving location exists 
pathlib.Path(config['save_to']).mkdir(parents=True, exist_ok=True)
###########################################################

## Load xsuite line
with open(xsuite_line) as fid:
    dd=json.load(fid)
line = xt.Line.from_dict(dd)

p0 = line.particle_ref

# Add longitudinal limit rectangle 
bucket_length = line.get_length()/Harmonic_Num
line.unfreeze() # if you had already build the tracker
line.append_element(element=xt.LongitudinalLimitRect(min_zeta=-bucket_length/2, max_zeta=bucket_length/2), name='long_limit')
line.build_tracker()

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

# Initialize benchmark dictionary for the kinetic mode, for the last kinetic run  
check_keys = {'epsn_x', 'epsn_y', 'particles', 'kinTx', 'kinTy', 'kinTz'}
tbt_checks = {nn: np.zeros((n_turns), dtype = float) for nn in check_keys}

# Also initiate for dataframes
dfs_kinetic = []
dfs_analytical = []

#%% First investigate KINETIC mode 
mode = 'kinetic'
for j in range(nr_kinetic_runs):
    print(f"\nRUN: {j+1}\n")
    
    particles = particles0.copy()
    
    # ----- Initialize IBS -----
    IBS = NagaitsevIBS()
    IBS.set_beam_parameters(particles)
    IBS.set_optic_functions(tw)
    
    # Write simulation value parameters for kinetic 
    if j == 0:
        sim_params = ''
        print()
        print("IBS class parameters:\n")
        attributes = {attr: getattr(IBS, attr) for attr in dir(IBS) if not callable(getattr(IBS, attr)) and not attr.startswith('__')}
        for attr, value in attributes.items():
            print(f'{attr}: {value}')
            sim_params+=f'{attr}: {value}\n'
        print()
    
        with open(f"{save_to}/sim_params.txt", "w") as file:
          file.write(sim_params)

    
    dt = 1./IBS.frev # consecutive turns/frev, used only for analytical, 
    
    ## Initialize dictionaries
    record_emit = {'eps_x', 'eps_y', 'sig_delta', 'bl', 'kinTx', 'kinTy', 'kinTz'}
    turn_by_turn = {nn: np.zeros((n_turns), dtype = float) for nn in record_emit}
    
    # --- Initialize 
    sig_x = np.std(particles.x[particles.state > 0])
    sig_y = np.std(particles.y[particles.state > 0])
    sig_delta = np.std(particles.delta[particles.state > 0])
    turn_by_turn['bl'][0]        = np.std(particles.zeta[particles.state > 0])
    turn_by_turn['sig_delta'][0] = sig_delta
    
    # Decide whether to use manual calculation of emittance or the StatisticalEmittance 
    turn_by_turn['eps_x'][0] = (sig_x**2 - (tw['dx'][0] * sig_delta)**2) / tw['betx'][0]
    turn_by_turn['eps_y'][0] = sig_y**2 / tw['bety'][0] 
    
    tbt_checks['epsn_x'][0] = (sig_x**2 - (tw['dx'][0] * sig_delta)**2) / tw['betx'][0]
    tbt_checks['epsn_y'][0] = sig_y**2 / tw['bety'][0] 
    tbt_checks['particles'][0] = sum(particles.state > 0)

    # Track over the given number of turns 
    for i in range(1, n_turns):
    
        print(f'Run {j+1}: Turn = {i}')
        print('N_part = ',len(particles.x[particles.state > 0]))
    
        # Calculate normalized emittances 
        sig_x = np.std(particles.x[particles.state > 0])
        sig_y = np.std(particles.y[particles.state > 0])
        sig_delta = np.std(particles.delta[particles.state > 0])
        turn_by_turn['bl'][i]        = np.std(particles.zeta[particles.state > 0])
        turn_by_turn['sig_delta'][i] = sig_delta
        turn_by_turn['eps_x'][i]     = (sig_x**2 - (tw['dx'][0] * sig_delta)**2) / tw['betx'][0]
        turn_by_turn['eps_y'][i]     = sig_y**2 / tw['bety'][0] 
        
        tbt_checks['epsn_x'][i] = turn_by_turn['eps_x'][i]  
        tbt_checks['epsn_y'][i] = turn_by_turn['eps_y'][i]
        tbt_checks['particles'][i] = sum(particles.state > 0)
        
        # Apply the kinetic kick from the calculated coefficient 
        if (i % IBS_step == 0) or (i==1):
           IBS.calculate_kinetic_coefficients(particles)
        IBS.apply_kinetic_kick(particles)
        print(f"Dx: {IBS.Dx}, Dy: {IBS.Dy}, Dz: {IBS.Dz}")
        print(f"Fx: {IBS.Fx}, Fy: {IBS.Fy}, Fz: {IBS.Fz}")

        # Check the growth rates 
        tbt_checks['kinTx'][i] = IBS.Dx - IBS.Fx
        tbt_checks['kinTy'][i] = IBS.Dy - IBS.Fy
        tbt_checks['kinTz'][i] = IBS.Dz - IBS.Fz
    
        tracker.track(particles)

    Emitt = []
    Emitt.append(turn_by_turn['eps_x'])
    Emitt.append(turn_by_turn['eps_y'])
    Emitt.append(turn_by_turn['sig_delta'])
    Emitt.append(turn_by_turn['bl'])
    Emitt.append(tbt_checks['kinTx'])
    Emitt.append(tbt_checks['kinTy'])
    Emitt.append(tbt_checks['kinTz'])
    
    df_kinetic = pd.DataFrame(np.array(Emitt).T, columns=["eps_x", "eps_y", "sig_delta", "bl", "kinTx", "kinTy", "kinTz"])
    df_kinetic.index.name = 'Turn'
    df_kinetic.to_parquet("{}/xsuite_run{}_{}_{}.parquet".format(save_to, j, mode, n_turns))
    dfs_kinetic.append(df_kinetic)
  
    
#%% Then investigate ANALYTICAL kick
mode = 'analytical'
print(f"Model: {mode}")

particles = particles0.copy()

# ----- Initialize IBS -----
IBS = NagaitsevIBS()
IBS.set_beam_parameters(particles)
IBS.set_optic_functions(tw)

dt = 1./IBS.frev # consecutive turns/frev, used only for analytical, 

## Initialize dictionaries
record_emit = {'eps_x', 'eps_y', 'sig_delta', 'bl', 'bl_proton_version', 'Ixx', 'Iyy', 'Ipp'}
turn_by_turn_analytical = {nn: np.zeros((n_turns), dtype = float) for nn in record_emit}

# --- Initialize 
sig_x = np.std(particles.x[particles.state > 0])
sig_y = np.std(particles.y[particles.state > 0])
sig_delta = np.std(particles.delta[particles.state > 0])
turn_by_turn_analytical['bl'][0] = np.std(particles.zeta[particles.state > 0])
turn_by_turn_analytical['bl_proton_version'][0]  = np.std(particles.zeta[particles.state > 0])
turn_by_turn_analytical['sig_delta'][0] = sig_delta

# Decide whether to use manual calculation of emittance or the StatisticalEmittance 
turn_by_turn_analytical['eps_x'][0] = (sig_x**2 - (tw['dx'][0] * sig_delta)**2) / tw['betx'][0]
turn_by_turn_analytical['eps_y'][0] = sig_y**2 / tw['bety'][0] 

for i in range(1, n_turns):

    print('Turn = ', i)
    print('N_part = ',len(particles.x[particles.state > 0]))

    # Calculate the IBS integrals and add the kick 
    if (i % IBS_step == 0) or (i==1):
         IBS.calculate_integrals(
             turn_by_turn_analytical['eps_x'][i-1],
             turn_by_turn_analytical['eps_y'][i-1],
             turn_by_turn_analytical['sig_delta'][i-1],
             turn_by_turn_analytical['bl'][i-1]
             )
    Emit_x, Emit_y, Sig_M = IBS.emit_evol(turn_by_turn_analytical['eps_x'][i-1],
                                          turn_by_turn_analytical['eps_y'][i-1],
                                          turn_by_turn_analytical['sig_delta'][i-1],
                                          turn_by_turn_analytical['bl'][i-1], 
                                          dt
                                          )
    
    Sigma_E = Sig_M*IBS.betar**2
    BunchL = ion_BunchLength(IBS.Circu, Harmonic_Num, IBS.EnTot, IBS.slip, 
               Sigma_E, IBS.betar, RF_Voltage*1e-3, Energy_loss, IBS.Ncharg)
    BunchL_proton_version = BunchLength(IBS.Circu, Harmonic_Num, IBS.EnTot, IBS.slip, 
                   Sigma_E, IBS.betar, RF_Voltage*1e-3, Energy_loss, IBS.Ncharg)
    
    turn_by_turn_analytical['bl'][i]        = BunchL
    turn_by_turn_analytical['bl_proton_version'][i] = BunchL_proton_version
    turn_by_turn_analytical['sig_delta'][i] = Sig_M
    turn_by_turn_analytical['eps_x'][i]     = Emit_x
    turn_by_turn_analytical['eps_y'][i]     = Emit_y
  
    turn_by_turn_analytical['Ixx'][i] = IBS.Ixx
    turn_by_turn_analytical['Iyy'][i] = IBS.Iyy
    turn_by_turn_analytical['Ipp'][i] = IBS.Ipp  

    tracker.track(particles)

Emitt = []
Emitt.append(turn_by_turn_analytical['eps_x'])
Emitt.append(turn_by_turn_analytical['eps_y'])
Emitt.append(turn_by_turn_analytical['sig_delta'])
Emitt.append(turn_by_turn_analytical['bl'])
Emitt.append(turn_by_turn_analytical['bl_proton_version'])
Emitt.append(turn_by_turn_analytical['Ixx'])
Emitt.append(turn_by_turn_analytical['Iyy'])
Emitt.append(turn_by_turn_analytical['Ipp'])

df_analytical = pd.DataFrame(np.array(Emitt).T, columns=["eps_x", "eps_y", "sig_delta", "bl", "bl_proton_version", 'Ixx', 'Iyy', 'Ipp'])
df_analytical.index.name = 'Turn'
df_analytical.to_parquet("{}/xsuite_{}_{}.parquet".format(save_to, mode, n_turns))
        
# Also save the TBT check for the kinetic growth rate integrals 
df_tbt = pd.DataFrame(tbt_checks)
df_tbt.to_parquet("{}/tbt_checks_{}_turns.parquet".format(save_to, n_turns))

