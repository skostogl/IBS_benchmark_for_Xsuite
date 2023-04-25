import os
import xtrack as xt
import xobjects as xo
from cpymad.madx import Madx
import json
import sys
import numpy as np
from cpymad.madx import Madx

seq_name = 'sps'

optics = "/afs/cern.ch/eng/acc-models/sps/2021"

if not os.path.isdir("sps"):
    os.symlink(optics, "sps")

mad = Madx()

# Modified example from /afs/cern.ch/eng/acc-models/sps/2021/examples/job_lhc_ion.madx
mad.input('''
call,file="sps/sps.seq";
call,file="sps/strengths/lhc_ion.str";

Beam, particle=ion, mass=193.7, charge=82, energy = 1415.72;

use,sequence=sps;
''')

twthick = mad.twiss().dframe()

# Slice and match tunes
n_slice_per_element = 4
mad.input(f'''
select, flag=MAKETHIN, SLICE={n_slice_per_element}, thick=false;
MAKETHIN, SEQUENCE={seq_name}, MAKEDIPEDGE=true;
use, sequence={seq_name};

use,sequence=sps;
twiss;

qx0=26.30;
qy0=26.25;

call,file="sps/toolkit/macro.madx";
call,file="sps/toolkit/match_tune.madx";
''')


# Chromaticity matching - use values generated under Test4_Chromaticity matching here 
# as values are not displayed on acc-models 
# https://gitlab.cern.ch/elwaagaa/xsuite-sps-ps-sequence-benchmarker
dq1 = -3.460734474533172e-09  # basically zero 
dq2 = -3.14426538905229e-09
mad.use('sps')
mad.exec("sps_define_sext_knobs();")
mad.exec("sps_set_chroma_weights_q26();")
mad.input(f"""match;
global, dq1={dq1};
global, dq2={dq2};
vary, name=qph_setvalue;
vary, name=qpv_setvalue;
jacobian, calls=10, tolerance=1e-25;
endmatch;""")

####### SET CAVITY VOLTAGE - with info from Hannes #######
# 6x200 MHz cavities: actcse, actcsf, actcsh, actcsi (3 modules), actcsg, actcsj (4 modules)
# acl 800 MHz cavities
# acfca crab cavities
# Ions: all 200 MHz cavities: 1.7 MV, h=4653
harmonic_nb = 4653
nn = 'actcse.31632'
mad.sequence.sps.elements[nn].lag = 0
mad.sequence.sps.elements[nn].volt = 1.7*82 # different convention between madx and xsuite
mad.sequence.sps.elements[nn].freq = mad.sequence[seq_name].beam.freq0*harmonic_nb

twthin = mad.twiss().dframe()

line = xt.Line.from_madx_sequence(mad.sequence[seq_name])
mad_beam = mad.sequence['sps'].beam
import xpart as xp
line.particle_ref = xp.Particles(
        p0c = mad_beam.pc*1e9,
        q0 = mad_beam.charge,
        mass0 = mad_beam.mass*1e9)

nn = 'actcse.31632'
harmonic_nb = 4653
line[nn] 
line[nn].lag = 0.0  # below transition
line[nn].voltage = 1.7e6
line[nn].frequency = mad.sequence[seq_name].beam.freq0*1e6*harmonic_nb

tracker = xt.Tracker(line=line)
tw_xtrack = tracker.twiss()

# Save line for tracking
folder_name = 'xsuite_line'
os.makedirs(folder_name, exist_ok=True)
with open(folder_name +'/sps_line_ions_for_tracking.json', 'w') as fid:
  json.dump(line.to_dict(), fid, cls=xo.JEncoder)


mad.input(f'''
save,sequence=sps,file=SPS_Q26_injection_ions.seq;
twiss, file="filename.tfs";
''')

