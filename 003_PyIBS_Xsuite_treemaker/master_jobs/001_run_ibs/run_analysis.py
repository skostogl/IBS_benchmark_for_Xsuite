# IBS module from M. Zampetakis in https://github.com/MichZampetakis/IBS_for_Xsuite
import tree_maker
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

tree_maker.tag_json.tag_it(config['log_file'], 'started')  

import sys
sys.path.append(config['ibs_lib_path'])
from lib.IBSfunctions import *

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

if config['modes'] == ['all']:
  modes = ['analytical', 'kinetic', 'simple']
else:
  modes = config['modes']
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

# ----- Initialize IBS -----
IBS = NagaitsevIBS()
IBS.set_beam_parameters(particles0)
#IBS.Npart = bunch_intensity
IBS.set_optic_functions(tw)


sim_params = ''
print()
print("IBS class parameters:\n")
attributes = {attr: getattr(IBS, attr) for attr in dir(IBS) if not callable(getattr(IBS, attr)) and not attr.startswith('__')}
for attr, value in attributes.items():
    print(f'{attr}: {value}')
    sim_params+=f'{attr}: {value}\n'
sim_params+=f"modes: {modes}\n"
print()

with open("sim_params.txt", "w") as file:
  file.write(sim_params)


#for mode in ['analytical', 'kinetic', 'simple']:
for mode in modes:
  print(f"Model: {mode}")

  particles = particles0.copy()

  dt = 1./IBS.frev # consecutive turns/frev, used only for analytical, 
  
    ## Initialize dictionaries
  turn_by_turn = {}
  
  record_emit = ['eps_x', 'eps_y', 'sig_delta', 'bl']
  for nn in record_emit:
      turn_by_turn[nn] = np.zeros((n_turns), dtype = float)

  # --- Initialize 
  sig_x = np.std(particles.x[particles.state > 0])
  sig_y = np.std(particles.y[particles.state > 0])
  sig_delta = np.std(particles.delta[particles.state > 0])
  turn_by_turn['bl'][0]        = np.std(particles.zeta[particles.state > 0])
  turn_by_turn['sig_delta'][0] = sig_delta
  turn_by_turn['eps_x'][0]     = (sig_x**2 - (tw['dx'][0] * sig_delta)**2) / tw['betx'][0]
  turn_by_turn['eps_y'][0]     = sig_y**2 / tw['bety'][0] 

  for i in range(1, n_turns):
      #import time
      #start = time.time()
  
      print('Turn = ', i)
      print('N_part = ',len(particles.x[particles.state > 0]))
  
      
      if mode != 'analytical':
        sig_x = np.std(particles.x[particles.state > 0])
        sig_y = np.std(particles.y[particles.state > 0])
        sig_delta = np.std(particles.delta[particles.state > 0])
        turn_by_turn['bl'][i]        = np.std(particles.zeta[particles.state > 0])
        turn_by_turn['sig_delta'][i] = sig_delta
        turn_by_turn['eps_x'][i]     = (sig_x**2 - (tw['dx'][0] * sig_delta)**2) / tw['betx'][0]
        turn_by_turn['eps_y'][i]     = sig_y**2 / tw['bety'][0] 
      
      if mode == 'analytical':
        if (i % IBS_step == 0) or (i==1):
             IBS.calculate_integrals(turn_by_turn['eps_x'][i-1],turn_by_turn['eps_y'][i-1],turn_by_turn['sig_delta'][i-1],turn_by_turn['bl'][i-1])
        Emit_x, Emit_y, Sig_M = IBS.emit_evol(turn_by_turn['eps_x'][i-1],turn_by_turn['eps_y'][i-1],turn_by_turn['sig_delta'][i-1],turn_by_turn['bl'][i-1], dt)
        
        Sigma_E = Sig_M*IBS.betar**2
        #BunchL = BunchLength(IBS.Circu, Harmonic_Num, IBS.EnTot, IBS.slip, 
        #           Sigma_E, IBS.betar, RF_Voltage*1e-3, Energy_loss, IBS.Ncharg)
        BunchL = ion_BunchLength(IBS.Circu, Harmonic_Num, IBS.EnTot, IBS.slip, 
                   Sigma_E, IBS.betar, RF_Voltage*1e-3, Energy_loss, IBS.Ncharg)
        
        turn_by_turn['bl'][i]        = BunchL
        turn_by_turn['sig_delta'][i] = Sig_M
        turn_by_turn['eps_x'][i]     = Emit_x
        turn_by_turn['eps_y'][i]     = Emit_y
  
      elif mode == "kinetic":
        if (i % IBS_step == 0) or (i==1):
             IBS.calculate_kinetic_coefficients(particles)
        IBS.apply_kinetic_kick(particles)
      elif mode == "simple":
        if (i % IBS_step == 0) or (i==1):
            IBS.calculate_simple_kick(particles)
        IBS.apply_simple_kick(particles)
  
      tracker.track(particles)
      #end = time.time()
      #print(end - start)
  
  Emitt = []
  Emitt.append(turn_by_turn['eps_x'])
  Emitt.append(turn_by_turn['eps_x']*IBS.betar*IBS.gammar)
  Emitt.append(turn_by_turn['eps_y'])
  Emitt.append(turn_by_turn['eps_y']*IBS.betar*IBS.gammar)
  Emitt.append(turn_by_turn['sig_delta'])
  Emitt.append(turn_by_turn['bl'])
  
  df = pd.DataFrame(np.array(Emitt).T, columns=["eps_x", "epsn_x","eps_y", "epsn_y","sig_delta", "bl"])
  df.index.name = 'Turn'

  pathlib.Path(config['save_to']).mkdir(parents=True, exist_ok=True)

  df.to_parquet(f"{save_to}/xsuite_{mode}.parquet")

tree_maker.tag_json.tag_it(config['log_file'], 'completed')  
