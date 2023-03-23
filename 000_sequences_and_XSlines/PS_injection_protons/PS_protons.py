#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generator of sequence files for PS protons at FLAT BOTTOM 
- also matching tunes and chromaticities 
"""
import xobjects as xo
import xtrack as xt
import xpart as xp

from cpymad.madx import Madx
import json
import os

# Specify optics repository 
optics = '/afs/cern.ch/eng/acc-models/ps/2022'
if not os.path.isdir("ps"):
    os.symlink(optics, "ps")

# Initiate MADX sequence, and call the relevant madx files and macros  
madx = Madx()
madx.call("{}/_scripts/macros.madx".format(optics))
madx.call("{}/scenarios/lhc/1_flat_bottom/ps_fb_lhc.beam".format(optics))   
madx.input('''BRHO      := BEAM->PC * 3.3356;''')
madx.call("{}/ps_mu.seq".format(optics))
madx.call("{}/ps_ss.seq".format(optics))
madx.call("{}/scenarios/lhc/1_flat_bottom/ps_fb_lhc.str".format(optics))
madx.call('PS_match_tunes_and_chroma.madx')
madx.use(sequence='ps')


############ MATCHING TUNES AND CHROMATICITY ############
# Put in values for PS protons to LHC: https://acc-models.web.cern.ch/acc-models/ps/2022/scenarios/lhc/
# When we perform matching, recall the MADX convention that all chromatic functions are multiplied by relativistic beta factor
# Thus, when we match for a given chromaticity, need to account for this factor to get correct value in Xsuite and PTC
beta0 = madx.sequence['ps'].beam.beta 
madx.input("qx = 6.21")
madx.input("qy = 6.245")
madx.input(f"qpx = 0.7273839232/{beta0}")
madx.input(f"qpy = -2.871405047/{beta0}")

madx.use("ps")
madx.input("seqedit, sequence=PS;")
madx.input("flatten;")
madx.input("endedit;")
madx.use("ps")
madx.input("select, flag=makethin, slice=5, thick=false;")
madx.input("makethin, sequence=ps, style=teapot, makedipedge=True;")
madx.use('ps')
madx.input("exec, match_tunes_and_chroma(qx, qy, qpx, qpy);")

madx.use(sequence = 'ps')
twiss_thin = madx.twiss()


############ XSUITE TRACKER ############   
context = xo.ContextCpu()
buf = context.new_buffer()   # using default initial capacity

# Perform Twiss command with MADX
madx.use(sequence='ps')
line = xt.Line.from_madx_sequence(madx.sequence['ps'])
madx_beam = madx.sequence['ps'].beam

particle_sample = xp.Particles(
        p0c = madx_beam.pc*1e9,
        q0 = madx_beam.charge,
        mass0 = madx_beam.mass*1e9)

line.particle_ref = particle_sample
tracker = xt.Tracker(_context=context, _buffer=buf, line=line) 


############ SET CAVITY VOLTAGE - with info from Alexandre Lasheen ############
# Ions: 10 MHz cavities: 200 kV for LHCPILOT, h=16
# In the ps_ss.seq, the 10 MHz cavities are PR_ACC10 - there seems to be 12 of them in the straight sections
harmonic_nb = 16
nn = 'pa.c10.11'  # for now test the first of the RF cavities 
V_RF = 200  # kV - in reality too high for one cavity but should suffice for our tracking purposes

# MADX sequence 
madx.sequence.ps.elements[nn].lag = 0
madx.sequence.ps.elements[nn].volt = V_RF*1e-3*particle_sample.q0 # different convention between madx and xsuite
madx.sequence.ps.elements[nn].freq = madx.sequence['ps'].beam.freq0*harmonic_nb

# Xsuite sequence 
line[nn].lag = 0  # 0 if below transition
line[nn].voltage =  V_RF*1e3 # In Xsuite for ions, do not multiply by charge as in MADX
line[nn].frequency = madx.sequence['ps'].beam.freq0*1e6*harmonic_nb


# Check Twiss command
madx.use(sequence = 'ps')
twiss_thin_RF = madx.twiss()
twiss_xtrack_RF = tracker.twiss()  

# Save MADX sequence
madx.command.save(sequence='ps', file='PS_2022_Protons_matched_with_RF.seq', beam=True)

# Save Xsuite sequence
with open('PS_2022_Protons_matched_with_RF.json', 'w') as fid:
    json.dump(line.to_dict(), fid, cls=xo.JEncoder)
    
# Print sanity checks 
print("\nPROTONS WITH RF: XTRACK vs MADX sequence:")
print("MAD-X thin:   " f"Qx  = {twiss_thin_RF.summary['q1']:.8f}"            f"   Qy  = {twiss_thin_RF.summary['q2']:.8f}")
print("Xsuite:       " f"Qx  = {twiss_xtrack_RF['qx']:.8f}"                  f"   Qy  = {twiss_xtrack_RF['qy']:.8f}\n")
print("MAD-X thin:   " f"Q'x = {twiss_thin_RF.summary['dq1']*beta0:.7f}"     f"   Q'y = {twiss_thin_RF.summary['dq2']*beta0:.7f}")
print("Xsuite:       " f"Q'x = {twiss_xtrack_RF['dqx']:.7f}"                 f"   Q'y = {twiss_xtrack_RF['dqy']:.7f}\n")
print("MAD-X thin:   " f"alpha_p = {twiss_thin_RF.summary.alfa:.7f}")
print("Xsuite:       " f"alpha_p = {twiss_xtrack_RF['momentum_compaction_factor']:.7f}")