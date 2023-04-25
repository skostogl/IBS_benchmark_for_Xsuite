# This is a modified example of M. Zampetakis in https://github.com/MichZampetakis/IBS_for_Xsuite/blob/main/CLIC_DRs_IBS_allmodels.py

import matplotlib.pyplot as plt
# pip install xsuite
import xtrack as xt
# pip install cpymad
from cpymad.madx import Madx
import xtrack as xt
import xpart as xp
import sys
import pandas as pd
import yaml

with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

sys.path.append(config['ibs_lib_path'])
from lib.IBSfunctions import *


madx = Madx()

# Simulation parameters
sequence = config['sequence']
energy=config['energy']
exn = config["emit_x"]*1e-6
eyn = config["emit_y"]*1e-6
RF_voltage = config["V0max"]
harm_number = config["h"]
cc_name_knobs = config["cc_name_knobs"]
bunch_intensity= float(config["bunch_intensity"])
particle = config["particle"]
sequence_name = config["sequence_name"]
mass = config["mass"] 
radius = config["radius"]
flag_IBS = config["flag_IBS"]
flag_SR = config["flag_SR"]
dt = config["ibs_step"]
duration = config["duration"]

madx.input(f'''
beam, particle={particle}, energy={energy}, mass={mass}, charge={config['charge']}, npart={bunch_intensity}, exn={exn}, eyn={eyn};
''')

madx.call(file=sequence)

madx.use(sequence=sequence_name)
madx.input(f'''
twiss,centre;
''')

twiss = madx.table.twiss.dframe()
summ = madx.table.summ.dframe()

energy_GEV = madx.table.twiss.summary.energy
pc_GEV = madx.table.twiss.summary.pc

p0 = xp.Particles(mass0=mass*1e9, q0=config['charge'], p0c=pc_GEV*1e9)
gt = summ["gammatr"].values[0]
gamma = madx.table.twiss.summary.gamma
my_eta = abs(1./gt**2-1./gamma**2)
betar = np.sqrt(1.0-1.0/gamma**2)
bunch_length = config["blns"]*1e-9/4.0*(c*betar)
bunch_length_level = config["bl_lev"]*1e-9/4.0*(c*betar)
madx.beam.sigt = bunch_length

# Initialize IBS class
IBS = NagaitsevIBS()
IBS.set_beam_parameters(p0)
IBS.Npart = bunch_intensity

if not np.isclose(energy_GEV, IBS.EnTot):
    IBS.EnTot=energy_GEV
    IBS.gammar = gamma
    IBS.betar = betar

twiss = twiss.to_dict(orient='list')
twiss={key: np.array(values) for key, values in twiss.items()}
twiss["slip_factor"] = my_eta
IBS.set_optic_functions(twiss)

emit_x_all   = []
emit_y_all   = []
sigma_z_all  = []
sigma_ns_all = []
txh          = []
tyh          = []
tlh          = []
time         = []

emit_x_all.append(exn/(IBS.betar*IBS.gammar))
emit_y_all.append(eyn/(IBS.betar*IBS.gammar))

sigma_z_all.append(bunch_length) 
sigma_ns_all.append(config["blns"])
txh.append(0.0)
tyh.append(0.0)
tlh.append(0.0)
time.append(0.0)

#madx.input("show, beam;")
# Check IBS attributes
print()
print("IBS class parameters:\n")
attributes = {attr: getattr(IBS, attr) for attr in dir(IBS) if not callable(getattr(IBS, attr)) and not attr.startswith('__')}
for attr, value in attributes.items():
    print(f'{attr}: {value}')
print()

for time_step in range(int(duration/dt)):
  emit_x         = emit_x_all[-1]
  emit_y         = emit_y_all[-1]
  
  #Sigma_E = EnergySpread(IBS.Circu, harm_number, IBS.EnTot*1e9, IBS.slip, sigma_z_all[-1], IBS.betar, RF_voltage*1e6, 0, IBS.Ncharg)
  Sigma_E = ion_EnergySpread(IBS.Circu, harm_number, IBS.EnTot*1e9, IBS.slip, sigma_z_all[-1], IBS.betar, RF_voltage*1e6, 0, IBS.Ncharg)

  Sig_M = Sigma_E/IBS.betar**2
  IBS.calculate_integrals(emit_x,emit_y,Sig_M,sigma_z_all[-1])

  Emit_x, Emit_y, Sig_M = IBS.emit_evol(emit_x,emit_y,Sig_M,sigma_z_all[-1], dt)
  Sigma_E = Sig_M*IBS.betar**2
  #BunchL = BunchLength(IBS.Circu, harm_number, IBS.EnTot, IBS.slip,
  #                  Sigma_E, IBS.betar, RF_voltage*1e-3, 0.0, IBS.Ncharg)
  BunchL = ion_BunchLength(IBS.Circu, harm_number, IBS.EnTot*1e9, IBS.slip, Sigma_E, IBS.betar, RF_voltage*1e6,0.0, IBS.Ncharg)
  print(BunchL)
  
  print(f"Step {time_step}: exn={Emit_x*IBS.betar*IBS.gammar}, eyn={Emit_y*IBS.betar*IBS.gammar}, bl={BunchL}")
  emit_x_all.append(Emit_x)
  emit_y_all.append(Emit_y)
  sigma_z_all.append(BunchL)
  sigma_ns_all.append(BunchL/(c*betar)*4.*1e9)
  time.append( (time_step+1)*dt)
  txh.append(1.0/float(IBS.Ixx)/3600.0)
  tyh.append(1.0/float(IBS.Iyy)/3600.0)
  tlh.append(1.0/float(IBS.Ipp)/3600.0)

  #print("..........")
  #print(IBS.Ixx, IBS.Iyy, IBS.Ipp)
  #print(txh[-1], tyh[-1], tlh[-1])
  #print("..........")

df = pd.DataFrame({'exin': np.array(emit_x_all)*IBS.betar*IBS.gammar, 'eyin':np.array(emit_y_all)*IBS.betar*IBS.gammar, 'sigma_z': sigma_z_all, 'Txh': txh, 'Tyh': tyh, 'Tlh':tlh, 'tt':time, 'npbb': bunch_intensity, 'bl_ns':sigma_ns_all, 'emit_x0':emit_x,'emit_y0':emit_y})
print(df)
df.to_parquet("IBS_output_python.parquet")
