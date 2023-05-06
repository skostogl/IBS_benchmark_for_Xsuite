# Treemaker for parametric studies

- tree_maker from G.Sterbini: https://gitlab.cern.ch/abpcomputing/sandbox/tree_maker
- example based on https://github.com/xsuite/example_DA_study

Including::
- setup.sh: to install python and tree_maker
- 001_make_folders.py: to create the tree
- 002_chronjob.py: to submit in htcondor

In order to facilitate parametric scans for the [IBS benchmarking](https://github.com/skostogl/IBS_benchmark_for_Xsuite), we use the [Tree maker](https://gitlab.cern.ch/abpcomputing/sandbox/tree_maker)1 in the `003_PyIBS_Xsuite_treemaker` folder 

## Main set-up instructions
1) First run [source setup.sh](https://github.com/skostogl/IBS_benchmark_for_Xsuite/blob/main/003_PyIBS_Xsuite_treemaker/setup.sh) on LXPLUS (comment and uncomment appropriate lines depending on operating system)
	1) If some lines do not work, run them manually
	2) For example, I had to do this with 
	```
	python -m pip install jupyterlab matplotlib pandas scipy numpy'
	
	git clone https://github.com/xsuite/tree_maker.git
	cd tree_maker
	git checkout release/v0.1.0
	python -m pip install .
	```
2) Then check `which python` --> this should give something like `/afs/cern.ch/work/e/elwaagaa/public/IBS/IBS_benchmark_for_Xsuite/003_PyIBS_Xsuite_treemaker/miniconda/bin/python`,  NOT `usr/bin/python`
	1) Running `ipython` and then `import tree_maker` should work fine
	2) Can add `source /afs/cern.ch/work/.../bin/activate` in the .bashrc file when opening afs
3) Add correct paths on AFS and EOS for saving data in the [001_make_folders.py](https://github.com/skostogl/IBS_benchmark_for_Xsuite/blob/main/003_PyIBS_Xsuite_treemaker/001_make_folders.py) 
	1) Define machine and parameters desired for the scan - include these correctly in the arrays
	2) If using Xsuite lines and sequences, make sure these paths exist first
	3) Then run `python 001_make_folders.py` in the correct Python environment - this should create a folder structure each containing `config.yaml` files and scripts `run_analysis.py` and `plot_analysis.py` adapted for each case. 
	4) Check that job folders e.g. 'SPS_injection' has been created. In these, try to run `./run.sh` after having commented away the output such that it shows in the terminal: `python {python_command} # > output_ht.txt 2> error_ht.txt Edited`  
	5) If the `run_analysis.py` starts, then it should work well 
4) Run `python 002_chronjob.py` to submit jobs to HTCondor

Remember 
- To install all necessary Python packages in the correct python path, e.g `python -m pip install cpymad`
- Make sure to run all scripts to generate the correct sequence files 

More [AFS info](https://twiki.cern.ch/twiki/bin/view/Main/NikhefDiskResources)
