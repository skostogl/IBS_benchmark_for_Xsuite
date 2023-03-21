import cpymad
import numpy as np
import pandas as pd
import math
from cpymad.madx import Madx
madx = Madx()
import cpymad
from cpymad.madx import Madx
from scipy.constants import c
import yaml

with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

def dpp_to_bl(RC, En, nc, RF_voltage, U0, betar, harm_number, etap, dpp):
    return c*RC*math.acos((En*nc*(RF_voltage-U0)*betar**2-dpp**2*En**2*harm_number*np.pi*betar**4*abs(etap))/(En*nc*(RF_voltage-U0)*betar**2))/(2*c*harm_number*np.pi*betar)

def bl_to_dpp(RC, En, nc, RF_voltage, U0, betar, harm_number, etap, bl_i):
    return np.sqrt(2/np.pi)*np.sqrt(nc*(RF_voltage-U0)*(np.sin(bl_i*harm_number*np.pi*betar/RC))**2)/np.sqrt(En*harm_number*abs(etap))/betar

madx = Madx()

# Define beam parameters
sequence = config['sequence']
energy = config["energy"]
exn = config["emit_x"]*1e-6
eyn = config["emit_y"]*1e-6
RF_voltage = config["V0max"]*1e-3
harm_number = config["h"]
cc_name_knobs = config["cc_name_knobs"]
bunch_intensity= config["bunch_intensity"]
particle = config["particle"]
sequence_name = config["sequence_name"]
mass = config["mass"] #938.2723*1e-3 # [GeV/c^2] me 0.51099906*1e-3;
radius = config["radius"] #1.534698e-18 #; ![m] re = 2.8179409e-18;
nc = config['charge']

flag_IBS = config["flag_IBS"]
flag_SR = config["flag_SR"]

ibs_step = config["ibs_step"]
duration = config["duration"]
###############################################################
###############################################################


madx.input(f'''
beam, particle={particle}, energy={energy}, mass={mass}, charge={config['charge']};
''')

madx.call(file=sequence)

madx.use(sequence=sequence_name)
madx.twiss()

twiss = madx.table.twiss.dframe()
summ = madx.table.summ.dframe()

RC = summ["length"].values[0]
gt = summ["gammatr"].values[0]
gamma = madx.table.twiss.summary.gamma
etap = abs(1./gt**2-1./gamma**2)
betar = np.sqrt(1.0-1.0/gamma**2)
En = madx.table.twiss.summary.energy

E0 = En/gamma;


madx.input(f"{cc_name_knobs}={RF_voltage*1e3}*{nc};")
U0 = 0
bunch_length = config["blns"]*1e-9/4.0*(c*betar) 
bunch_length_level = config["bl_lev"]*1e-9/4.0*(c*betar) 
bl_i = bunch_length
exi = exn/(gamma*betar);
eyi = eyn/(gamma*betar);

frev = betar*c/RC
fRF = harm_number*frev
TRF=1./fRF

dpp = bl_to_dpp(RC, En, nc, RF_voltage, U0, betar, harm_number, etap, bl_i)
dee = dpp*betar**2

madx.use(sequence=sequence_name)
madx.twiss()
twiss = madx.table.twiss.dframe()

madx.beam.ex = exi
madx.beam.ey = eyi
madx.beam.sigt = bunch_length
madx.beam.sige = dee
madx.beam.npart = bunch_intensity

if flag_SR:
    madx.input(f'''
    ! synchrotron radiation considered in all dipole magnets
    beam,radiate;
    ! when beam, radiate, cmp equilibrium emittances
    emit,deltap=0;
    value, beam->ex, beam->ey, beam->sigt, beam->sige, beam->et;
    ex0=beam->ex;
    ey0=beam->ey;
    sp0=beam->sige;
    ss0=beam->sigt;
    En=beam->energy;
    ! synchrotron radiation integrals
    twiss,chrom;
    ''')

    ex0 = madx.globals.ex0
    ey0 = madx.globals.ey0
    sp0 = madx.globals.sp0
    ss0 = madx.globals.ss0

    summ = madx.table.summ.dframe()

    I1 = summ["synch_1"].values[0]
    I2 = summ["synch_2"].values[0]
    I3 = summ["synch_3"].values[0]
    I4 = summ["synch_4"].values[0]
    I5 = summ["synch_5"].values[0]
 
    # Energy loss U0 per turn https://cas.web.cern.ch/sites/default/files/lectures/zakopane-2006/rivkin-dynamics.pdf slide 4
    Cgg = 4*np.pi/3*radius/(mass)**3
    U0 = Cgg/2/np.pi*En**4*I2
    
    # https://arxiv.org/pdf/1507.02213.pdf
    taux=2*En*(RC/c)/U0; # Beam size damping time
    tauy=2*En*(RC/c)/U0; # Beam size damping time
    taul=2*En*(RC/c)/U0/2;   # Bunch length damping time
    
    madx.input("beam,radiate=false;")
    #madx.globals["cc_name_knobs"] = 0
    madx.input(f"{cc_name_knobs}=0;")
    
    madx.beam.ex = exi
    madx.beam.ey = eyi
    madx.beam.sigt = bunch_length
    dee = dpp*betar**2
    madx.beam.sige = dee
    madx.beam.npart = bunch_intensity

else:
    U0, taux, tauy, taul, ex0, ey0, sp0, ss0= (0,)*8

madx.input(f'''
twiss,centre;
''')

results = []
current_time = 0
counter = 0
while current_time < duration + ibs_step:
    
    print(f"Currently at (mins): {current_time/60.0}/{duration/60.0}")

    madx.input(f'''
    init_ex = beam->Ex;
    init_ey = beam->Ey;
    bl_ns= beam->sigt/(clight*{betar})*4*1e9;
    ''')

    init_ex = madx.globals.init_ex 
    init_ey = madx.globals.init_ey
    init_bl_ns = madx.globals.bl_ns
    
    if flag_IBS:
        madx.input(f'''
        !show, beam;
        ibs;

        ! Beam size growth time
        Tx=1/ibs.tx/2;
        Ty=1/ibs.ty/2;
        Tl=1/ibs.tl/2;
        ''')
        Tx = madx.globals.Tx
        Ty = madx.globals.Ty
        Tl = madx.globals.Tl
        Txh = (1/madx.globals.Tx)/3600/2;
        Tyh = (1/madx.globals.Ty)/3600/2;
        Tlh = (1/madx.globals.Tl)/3600/2;
        #print(2.0*Tx, 2.0*Ty, 2.0*Tl) # to compare with IBS.Ixx from python IBS
        #print(Txh, Tyh, Tlh)
    else:
        Tx, Ty, Tl, Txh, Tyh, Tlh= (0, )*6
    
    if flag_SR:
        exi = (-ex0+np.exp(2*ibs_step*(Tx-1/taux))*(ex0+exi*(-1+Tx*taux)))/(-1+Tx*taux)
        eyi = (-ey0+np.exp(2*ibs_step*(Ty-1/tauy))*(ey0+eyi*(-1+Ty*tauy)))/(-1+Ty*tauy)
        dpp = (-sp0+np.exp(ibs_step*(Tl-1/taul))*(sp0+dpp*(-1+Tl*taul)))/(-1+Tl*taul)
    else:
        exi = exi*np.exp(2*ibs_step*Tx)
        eyi = eyi*np.exp(2*ibs_step*Ty)
        dpp = dpp*np.exp(ibs_step*Tl)
        U0 = 0

    dee = dpp*betar**2
    bl_i = dpp_to_bl(RC, En, nc, RF_voltage, U0, betar, harm_number, etap, dpp)
    
    if bl_i<bunch_length_level:
        bl_i = bunch_length_level
        dpp = bl_to_dpp(RC, En, nc, RF_voltage, U0, betar, harm_number, etap, bl_i)
        dee = dpp*betar**2

    madx.input(f'''
    beam,EX={exi};
    beam,EY={eyi};
    beam,sigt={bl_i};
    beam,sige={dee};
    beam,pc={np.sqrt(En**2-E0**2)};
    beam,energy={En};
    beam,gamma={gamma};
    beam,npart={bunch_intensity}; 
    show, beam;
    ''')
    counter+=1
    print(f"Step {counter}: exn={exi*betar*gamma}, eyn={eyi*betar*gamma}, bl={bl_i}")
    results.append([current_time, bunch_intensity, ex0, ey0, sp0, ss0, init_ex*(betar*gamma), init_ey*(betar*gamma) , dpp,init_bl_ns, En, Txh, Tyh, Tlh, RF_voltage, taux, tauy, taul])
    
    current_time += ibs_step
    
df = pd.DataFrame(results, columns = ["tt", "npbb", "ex0", "ey0", "sp0", "ss0", "exin", "eyin", "dpp", "bl_ns", "En", "Txh", "Tyh", "Tlh", "V0", "taux", "tauy", "taul"])
print(df)
df.to_parquet("IBS_output_madx.parquet")
#df.to_parquet(config["save_to"])
