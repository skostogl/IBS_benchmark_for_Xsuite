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
call,file="sps/strengths/lhc_q26.str";

Beam, particle=proton, pc = 451;

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

qx0=26.13;
qy0=26.18;

call,file="sps/toolkit/macro.madx";
call,file="sps/toolkit/match_tune.madx";
''')

# Match the chromaticities
# - values from https://acc-models.web.cern.ch/acc-models/sps/2021/scenarios/lhc_proton/
dq1 = 0.0083385932924699; 
dq2 = -0.0019057430908735074;
mad.exec("sps_define_sext_knobs();")
mad.exec("sps_set_chroma_weights_q26();")
mad.input(f"""match;
global, dq1={dq1};
global, dq2={dq2};
vary, name=qph_setvalue;
vary, name=qpv_setvalue;
jacobian, calls=10, tolerance=1e-25;
endmatch;""")


harmonic_nb = 4653
nn = 'actcse.31637'
mad.sequence.sps.elements[nn].lag = 0.5
mad.sequence.sps.elements[nn].volt = 3 # different convention between madx and xsuite
mad.sequence.sps.elements[nn].freq = mad.sequence[seq_name].beam.freq0*harmonic_nb

twthin = mad.twiss().dframe()

line = xt.Line.from_madx_sequence(mad.sequence[seq_name])
mad_beam = mad.sequence['sps'].beam
import xpart as xp
line.particle_ref = xp.Particles(
        p0c = mad_beam.pc*1e9,
        q0 = mad_beam.charge,
        mass0 = mad_beam.mass*1e9)


nn = 'actcse.31637'
harmonic_nb = 4653
line[nn] 
line[nn].lag = 180.0
line[nn].voltage = 3e6
line[nn].frequency = mad.sequence[seq_name].beam.freq0*1e6*harmonic_nb

tracker = xt.Tracker(line=line)
tw_xtrack = tracker.twiss()

# Save line for tracking
folder_name = 'xsuite_line'
os.makedirs(folder_name, exist_ok=True)
with open(folder_name +'/sps_line_protons_for_tracking.json', 'w') as fid:
  json.dump(line.to_dict(), fid, cls=xo.JEncoder)


mad.input(f'''
save,sequence=sps,file=SPS_Q26_top_protons.seq;
''')
