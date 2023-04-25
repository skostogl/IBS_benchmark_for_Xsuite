"""
Main script to compare different IBS modes: analytical, kinetic, simple - combining with Xsuite tracking 
IBS module from M. Zampetakis in https://github.com/MichZampetakis/IBS_for_Xsuite
"""

import sys
import pathlib
import json
import numpy as np
import xobjects as xo
import xtrack as xt
import xpart as xp
import pandas as pd
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

# Activate the different modes
modes = []
if config['mode_kinetic']: 
    modes.append('kinetic')
if config['mode_analytical']: 
    modes.append('analytical')
if config['mode_simple']: 
    modes.append('simple')

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

# ----- Initialize IBS object -----
IBS = NagaitsevIBS()
IBS.set_beam_parameters(particles0)
IBS.set_optic_functions(tw)

# Write used simulation parameters 
sim_params = ''
print("\nIBS class parameters:\n")
attributes = {attr: getattr(IBS, attr) for attr in dir(IBS) if not callable(getattr(IBS, attr)) and not attr.startswith('__')}
for attr, value in attributes.items():
    print(f'{attr}: {value}')
    sim_params+=f'{attr}: {value}\n'
print()

with open("sim_params.txt", "w") as file:
  file.write(sim_params)

###### ITERATE OVER THE CHOSEN MODES #######
for mode in modes:
    
    print(f"\nModel: {mode}\n")
    particles = particles0.copy()
    dt = 1./IBS.frev # consecutive turns/frev, used only for analytical, 
    
    ## Initialize dictionaries
    turn_by_turn = {}
    record_emit = ['eps_x', 'eps_y', 'epsn_x', 'epsn_y', 'sig_delta', 'bl']
    for nn in record_emit:
        turn_by_turn[nn] = np.zeros((n_turns), dtype = float)
    
    # --- Initialize first turn-by-turn data for all modes 
    sig_x = np.std(particles.x[particles.state > 0])
    sig_y = np.std(particles.y[particles.state > 0])
    sig_delta = np.std(particles.delta[particles.state > 0])
    turn_by_turn['bl'][0]        = np.std(particles.zeta[particles.state > 0])
    turn_by_turn['sig_delta'][0] = sig_delta
    turn_by_turn['eps_x'][0]     = (sig_x**2 - (tw['dx'][0] * sig_delta)**2) / tw['betx'][0]
    turn_by_turn['eps_y'][0]     = sig_y**2 / tw['bety'][0] 

    # Initialize dictionaries for the integrals and growth rates
    if mode == 'analytical':
        turn_by_turn_analytical_integrals = {}
        integral_keys = ['Ixx', 'Iyy', 'Ipp']
        for nn in integral_keys:
            turn_by_turn_analytical_integrals[nn] = np.zeros((n_turns), dtype = float)
    
    elif mode == 'kinetic':
        turn_by_turn_kinetic_integrals = {}
        integral_keys = ['kinTx', 'kinTy', 'kinTz']
        for nn in integral_keys:
            turn_by_turn_kinetic_integrals[nn] = np.zeros((n_turns), dtype = float)
            
    elif mode == 'simple':
        turn_by_turn_simple_integrals = {}
        integral_keys = ['Ixx', 'Iyy', 'Ipp']
        for nn in integral_keys:
            turn_by_turn_simple_integrals[nn] = np.zeros((n_turns), dtype = float)


    # TRACK the particles over the given number of turns 
    for i in range(1, n_turns):
    
        print('Turn = ', i)
        print('N_part active = ',len(particles.x[particles.state > 0]))
    
        # For the non-analytical cases, extract beam properties from particle object 
        if mode == 'kinetic' or mode == 'simple':
            sig_x = np.std(particles.x[particles.state > 0])
            sig_y = np.std(particles.y[particles.state > 0])
            sig_delta = np.std(particles.delta[particles.state > 0])
            turn_by_turn['bl'][i]        = np.std(particles.zeta[particles.state > 0])
            turn_by_turn['sig_delta'][i] = sig_delta
            turn_by_turn['eps_x'][i]     = (sig_x**2 - (tw['dx'][0] * sig_delta)**2) / tw['betx'][0]
            turn_by_turn['eps_y'][i]     = sig_y**2 / tw['bety'][0] 
        
        # For the analytical case, calculate the integrals for each X turns and calculate theoretical emittance evolution
        if mode == 'analytical':
            
            # Calculate the IBS integrals and add the kick 
            if (i % IBS_step == 0) or (i==1):
                 IBS.calculate_integrals(
                     turn_by_turn['eps_x'][i-1],
                     turn_by_turn['eps_y'][i-1],
                     turn_by_turn['sig_delta'][i-1],
                     turn_by_turn['bl'][i-1]
                     )
            Emit_x, Emit_y, Sig_M = IBS.emit_evol(
                                                turn_by_turn['eps_x'][i-1],
                                                turn_by_turn['eps_y'][i-1],
                                                turn_by_turn['sig_delta'][i-1],
                                                turn_by_turn['bl'][i-1], 
                                                dt
                                                  )
            # Calculate energy spread and bunch length 
            Sigma_E = Sig_M*IBS.betar**2
            BunchL = ion_BunchLength(IBS.Circu, Harmonic_Num, IBS.EnTot, IBS.slip, 
                       Sigma_E, IBS.betar, RF_Voltage*1e-3, Energy_loss, IBS.Ncharg)
            
            # Save the data 
            turn_by_turn['bl'][i]        = BunchL
            turn_by_turn['sig_delta'][i] = Sig_M
            turn_by_turn['eps_x'][i]     = Emit_x
            turn_by_turn['eps_y'][i]     = Emit_y
          
            turn_by_turn_analytical_integrals['Ixx'][i] = IBS.Ixx
            turn_by_turn_analytical_integrals['Iyy'][i] = IBS.Iyy
            turn_by_turn_analytical_integrals['Ipp'][i] = IBS.Ipp  

        # Apply the kicks and save the integrals 
        elif mode == "kinetic":
            if (i % IBS_step == 0) or (i==1):
                IBS.calculate_kinetic_coefficients(particles)
            IBS.apply_kinetic_kick(particles)
            print(f"Dx: {IBS.Dx}, Dy: {IBS.Dy}, Dz: {IBS.Dz}")
            print(f"Fx: {IBS.Fx}, Fy: {IBS.Fy}, Fz: {IBS.Fz}")

            # Check the growth rates 
            turn_by_turn_kinetic_integrals['kinTx'][i] = IBS.Dx - IBS.Fx
            turn_by_turn_kinetic_integrals['kinTy'][i] = IBS.Dy - IBS.Fy
            turn_by_turn_kinetic_integrals['kinTz'][i] = IBS.Dz - IBS.Fz
            
        elif mode == "simple":
            if (i % IBS_step == 0) or (i==1):
                IBS.calculate_simple_kick(particles)
            IBS.apply_simple_kick(particles)
            
            turn_by_turn_simple_integrals['Ixx'][i] = IBS.Ixx
            turn_by_turn_simple_integrals['Iyy'][i] = IBS.Iyy
            turn_by_turn_simple_integrals['Ipp'][i] = IBS.Ipp  
            
        # Track the particles
        tracker.track(particles)

    # Add the normalized emittances in one go
    turn_by_turn['epsn_x'] = IBS.betar*IBS.gammar*turn_by_turn['eps_x']
    turn_by_turn['epsn_y'] = IBS.betar*IBS.gammar*turn_by_turn['eps_y']

    # Save the turn-by-turn_data with the integrals 
    if mode == 'analytical':
        # Merge dictionaries 
        z1 = turn_by_turn | turn_by_turn_analytical_integrals
        pd_analytical = pd.DataFrame(z1)
        pd_analytical.index.name = 'Turn'
        pathlib.Path(config['save_to']).mkdir(parents=True, exist_ok=True)
        pd_analytical.to_parquet(f"{save_to}/xsuite_{mode}.parquet")
        
    elif mode == 'kinetic':
        z2 = turn_by_turn | turn_by_turn_kinetic_integrals
        pd_kinetic = pd.DataFrame(z2)
        pd_kinetic.index.name = 'Turn'
        pathlib.Path(config['save_to']).mkdir(parents=True, exist_ok=True)
        pd_kinetic.to_parquet(f"{save_to}/xsuite_{mode}.parquet")
        
    elif mode == 'simple':
        z3 = turn_by_turn | turn_by_turn_simple_integrals
        pd_simple = pd.DataFrame(z3)
        pd_simple.index.name = 'Turn'
        pathlib.Path(config['save_to']).mkdir(parents=True, exist_ok=True)
        pd_simple.to_parquet(f"{save_to}/xsuite_{mode}.parquet")
        
