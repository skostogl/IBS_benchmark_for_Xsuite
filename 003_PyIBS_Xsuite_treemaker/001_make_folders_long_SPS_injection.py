import pandas as pd
import tree_maker
from tree_maker import NodeJob
from tree_maker import initialize
import time
import os
from pathlib import Path
import itertools
import numpy as np
import yaml
from user_defined_functions import generate_run_sh
from user_defined_functions import generate_run_sh_htc

# Import the configuration
config=yaml.safe_load(open('config.yaml'))
study_name        = f"SPS_injection_ions_long_run"

# Loop over the following parameters
n_part_s = [5000]
n_turns_s = [800000]  # corresponds to slightly less than 20 seconds in the SPS
IBS_step_s = [50]
xsuite_line_s = ["/afs/cern.ch/work/e/elwaagaa/public/IBS/IBS_benchmark_for_Xsuite/000_sequences_and_XSlines/SPS_injection_ions/xsuite_line/sps_line_ions_for_tracking.json"]
bunch_intensity_s = [3.5e8]  # nominal bunch intensity 
emit_x_s = [1.2612]
emit_y_s = [0.9081]
sigma_z_s = [0.039, 0.23]  # short bunch, and nominal length
V0max_s = [3.0]  # voltage as read from the CCC in early 2023

# Fixed parameters
ibs_lib_path = "/afs/cern.ch/work/e/elwaagaa/public/IBS/IBS_for_Xsuite"

# Define what modes to run
mode_kinetic = 1
mode_analytical = 1
mode_simple = 0

#modes = ['all']
energy = 1415.72
h = 4653
mypython = "/afs/cern.ch/work/e/elwaagaa/public/IBS/IBS_benchmark_for_Xsuite/003_PyIBS_Xsuite_treemaker/miniconda/bin/activate"

save_to = f"/eos/user/e/elwaagaa/PhD/Projects/IBS/IBSresults_{study_name}"  # sofia uses "/eos/user/s/skostogl/IBSresults_{study_name}"
from pathlib import Path
Path(save_to).mkdir(parents=True, exist_ok=True)

children={}
for child, (n_part, n_turns, IBS_step, xsuite_line, bunch_intensity, emit_x, emit_y, sigma_z, V0max) in enumerate(itertools.product(n_part_s, n_turns_s, IBS_step_s, xsuite_line_s, bunch_intensity_s, emit_x_s, emit_y_s, sigma_z_s, V0max_s)):
    print("Study: ", n_part, n_turns, IBS_step, xsuite_line, bunch_intensity, emit_x, emit_y, sigma_z, V0max)
    children[f"{study_name}/{child:03}"] = {
                                    'ibs_lib_path':ibs_lib_path,
                                    'n_part':n_part,
                                    'n_turns': n_turns,
                                    'IBS_step': IBS_step,
                                    'xsuite_line': xsuite_line,
                                    'energy': energy,
                                    'bunch_intensity': bunch_intensity,
                                    'emit_x': emit_x,
                                    'emit_y': emit_y,
                                    'sigma_z': sigma_z,
                                    'V0max': V0max,
                                    'h': h,
                                    'mode_kinetic': mode_kinetic,
				    'mode_analytical': mode_analytical,
				    'mode_simple': mode_simple,
                                    #'modes': modes,
                                    'save_to':f"{save_to}/{child:03}",
                                    'log_file': f"{os.getcwd()}/{study_name}/{child:03}/tree_maker.log"
                                    }
config['root']['children'] = children
config['root']['setup_env_script'] = mypython

# Create tree object
start_time = time.time()
root = initialize(config)
print('Done with the tree creation.')
print("--- %s seconds ---" % (time.time() - start_time))

# From python objects we move the nodes to the file-system.
start_time = time.time()
#root.make_folders(generate_run_sh)
root.make_folders(generate_run_sh_htc)
print('The tree folders are ready.')
print("--- %s seconds ---" % (time.time() - start_time))

import shutil
shutil.move("tree_maker.json", f"tree_maker_{study_name}.json")
shutil.move("tree_maker.log", f"tree_maker_{study_name}.log")
