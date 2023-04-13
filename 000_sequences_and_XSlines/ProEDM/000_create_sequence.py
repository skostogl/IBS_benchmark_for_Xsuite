# An example from S. Siddique
from cpymad.madx import Madx
import numpy as np
from scipy.constants import c

madx = Madx()

mp      = 938.272046 #MeV/c^2
def gammarel(EGeV,m0=mp):
  return (EGeV*1.e3)/m0

def betarel(EGeV,m0=mp):
  g=gammarel(EGeV,m0=m0)
  return np.sqrt(1-1/g**2)

def cmp_frev(EGeV, C, m0=mp):
  b = betarel(EGeV)
  return (b*c/C)


CC_MV           = 4.0 
CC_LAG          = 0 # below transition
#CC_LAG          = 0.5 # above transition
CC_H            = 6
bunch_intensity = 1e9
energy_GEV      = 0.968
#energy_GEV      = 9.68

emit_x          = 1e-6*betarel(energy_GEV)*gammarel(energy_GEV)  # normalised emit
emit_y          = 1e-6*betarel(energy_GEV)*gammarel(energy_GEV) # normalised emit
sigma_z         = 692e-2
sigma_ns        = sigma_z/(c*betarel(energy_GEV))*4.
#print(sigma_z, sigma_ns*1e9)

C               = 123.778
frev            = cmp_frev(energy_GEV, C) # in Hz
dt              = 1e3/frev # time step (s) for IBS computation, every 5 minutes
duration        = 1e5/frev # total duration (s)

print("Creating sequence with the following parameters: \n")
print(f"RF voltage: {CC_MV}, RF LAG: {CC_LAG}, RF Harmonics number: {CC_H}, Bunch intensity: {bunch_intensity}, Energy in GeV: {energy_GEV}, Normalized emit (um): {emit_x*1e6}, {emit_y*1e6}, Bunch length (m): {sigma_z}, Bunch length in ns: {sigma_ns*1e9}, Circumference: {C}, Revolution frequency: {frev}, Revolution period: {1.0/frev}, beta: {betarel(energy_GEV)}, gamma: {gammarel(energy_GEV)}")

madx.input(f'''
option,-echo,threader;

	QDk=-0.3;
	Qfk=0.05;


	N.Ebend=16;
	l.Ebend=4.809875;
	l.Quad=0.200;
	
	angle=2*pi/(N.Ebend);
	
	k.Qf:=Qfk;
	k.Qd:=QDk;
	k.Qss:=Qss;
	
	Qf: QUADRUPOLE, L= l.Quad, K1=k.Qf;
	dq:  drift, l = 0.2;
	ds:  drift, l = 1.6;
	BPM: monitor, l=0.08;
	dEbend: drift, l=0.42;
	d: Marker;
	SSQ:  QUADRUPOLE, L= l.Quad, K1=k.Qss;
	Qd: QUADRUPOLE, L= l.Quad, K1=k.Qd;
 
    my_vrf = {CC_MV};
    CAVITY:   RFCAVITY, L=0.1, VOLT:=my_vrf,LAG={CC_LAG},HARMON={CC_H};

	index=0.199;
	
	
	EB: MATRIX, L =l.Ebend, 
	RM11 = 0.85418, RM12 = 3.30871, RM13 = 0, RM14 =0, RM15 =0, RM16 = 1.29205,

	RM21 = -0.0817166, RM22 = 0.85418, RM23 = 0, RM24 =0, RM25 =0, RM26 = 0.724056,

	RM31 = 0 , RM32 = 0 , RM33 = 1., RM34  = 3.47954 ,RM35 =0, RM36 = 0,

	RM41 = 0 , RM42 = 0 , RM43 = 0., RM44  = 1. ,RM45 =0, RM46 = 0,

	RM51 = -0.724056, RM52 = -1.29205, RM53 = 0, RM54 =0, RM55 =1, RM56 = 2.94856,

	RM61 = 0, RM62 = 0 , RM63 = 0, RM64 =0, RM65 =0, RM66 = 1;
	
	onecell: Line=(d,SSQ,dq,BPM,ds,ds,dq,Qf,d,Qf,dq,BPM,dEbend,EB,EB,dEbend,BPM,dq,Qd,d,Qd,dq,BPM,dEbend,EB,EB,dEbend,BPM,dq,Qf,d,Qf,dq,ds,ds,dq,SSQ,d);
	proedm: Line=(4*onecell + CAVITY);
	
	beam, PARTICLE=PROTON,MASS=pmass,npart = {bunch_intensity}, CHARGE=1,ENERGY={energy_GEV},EXN={emit_x},EYN={emit_y},bunched=True;
	use, sequence=proedm;
	
	select, flag=twiss, clear;
	select, flag=twiss, column= NAME,S,BETX,ALFX,MUX,X,PX,DX,DPX,BETY,ALFY,MUY,Y,PY,DY,DPY ;
	twiss,save,file="PTR-33M.tfs";
	
	
	!plot,table=twiss,haxis=s,hmin=0,hmax=110,vaxis1=BETX,BETY,DX,colour=100,Range=#S/#E,Noversion=True,noline=True,title="PTR Lattice", file="lattice";
	
   save,sequence=proedm,file="proedm.seq";

''')
