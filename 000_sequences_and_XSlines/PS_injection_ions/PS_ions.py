#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generator of sequence files for PS Pb ions at FLAT BOTTOM 
- also matching tunes and chromaticities 
"""
import xobjects as xo
import xtrack as xt
import xpart as xp

from cpymad.madx import Madx
import json
import os

# Specify optics repository 
optics = '/afs/cern.ch/eng/acc-models/ps/2022/'
if not os.path.isdir("ps"):
    os.symlink(optics, "ps")

# Pb ion beam parameters for PS 
m_u = 931.49410242e6  # 1 Dalton in eV/c^2 -- atomic mass unit 
atomic_mass_pb208 = 207.9766519 # Pb208 atomic mass from The AME2016 atomic mass evaluation (II).
m_ion = atomic_mass_pb208*m_u  
Z = 54
E_kin = 0.0722*1e9*atomic_mass_pb208 # total kinetic energy in eV per particle at injection, 
#   from LIU parameter table https://edms.cern.ch/ui/file/1420286/2/LIU-Ions_beam_parameter_table.pdf
E_tot = m_ion + E_kin

# Initiate MADX sequence, and call the relevant madx files and macros  
madx = Madx()
madx.call("{}/_scripts/macros.madx".format(optics))

# Correct beam command tp obtain same Brho and beta0 as Hannes' and Isabelle Table 5 (https://cds.cern.ch/record/2749453)
madx.input(" \
           Beam, particle=ion, mass={}, charge={}, energy = {}; \
           DPP:=BEAM->SIGE*(BEAM->ENERGY/BEAM->PC)^2;  \
           ".format(m_ion/1e9, Z, E_tot/1e9))

madx.call("{}/ps_mu.seq".format(optics))
madx.call("{}/ps_ss.seq".format(optics))
madx.call("{}/scenarios/lhc_ion/1_flat_bottom/ps_fb_ion.str".format(optics))
madx.call('PS_match_tunes_and_chroma.madx')
madx.use(sequence='ps')


############ MATCHING OF TUNES AND CHROMATICITY ############  
# Put in Twiss summary values for PS ions: https://acc-models.web.cern.ch/acc-models/ps/2022/scenarios/lhc_ion/

# When we perform matching, recall the MADX convention that all chromatic functions are multiplied by relativistic beta factor
# Thus, when we match for a given chromaticity, need to account for this factor to get correct value in Xsuite and PTC
beta0 = madx.sequence['ps'].beam.beta 
madx.input("qx = 6.21")
madx.input("qy = 6.245")
madx.input(f"qpx = -5.26716824/{beta0}")
madx.input(f"qpy = -7.199251093/{beta0}")

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
# Ions: 10 MHz cavities: 1.7 MV, h=16
# In the ps_ss.seq, the 10 MHz cavities are PR_ACC10 - there seems to be 12 of them in the straight sections
harmonic_nb = 16
nn = 'pa.c10.11'  # for now, activate the first of the RF cavities 
V_RF = 38.0958  # kV

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
madx.command.save(sequence='ps', file='PS_2022_Pb_ions_matched_with_RF.seq', beam=True)

# Save Xsuite sequence
with open('PS_2022_Pb_ions_matched_with_RF.json', 'w') as fid:
    json.dump(line.to_dict(), fid, cls=xo.JEncoder)
    
# Print sanity checks 
print("\nPB IONS WITH RF: XTRACK vs MADX sequence:")
print("MAD-X thin:   " f"Qx  = {twiss_thin_RF.summary['q1']:.8f}"            f"   Qy  = {twiss_thin_RF.summary['q2']:.8f}")
print("Xsuite:       " f"Qx  = {twiss_xtrack_RF['qx']:.8f}"                  f"   Qy  = {twiss_xtrack_RF['qy']:.8f}\n")
print("MAD-X thin:   " f"Q'x = {twiss_thin_RF.summary['dq1']*beta0:.7f}"     f"   Q'y = {twiss_thin_RF.summary['dq2']*beta0:.7f}")
print("Xsuite:       " f"Q'x = {twiss_xtrack_RF['dqx']:.7f}"                 f"   Q'y = {twiss_xtrack_RF['dqy']:.7f}\n")
print("MAD-X thin:   " f"alpha_p = {twiss_thin_RF.summary.alfa:.7f}")
print("Xsuite:       " f"alpha_p = {twiss_xtrack_RF['momentum_compaction_factor']:.7f}")