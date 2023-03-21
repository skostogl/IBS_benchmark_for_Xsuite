A series of examples including IBS with the aim of integrating IBS effects in Xsuite:

- IBS MADX code: modified version on March 2023 wrapped in cpymad. Original code from: S. Papadopoulou and F. Antoniou in https://github.com/stpapado/calculateIBS_MADX.git
- IBS Python code: IBS library by M. Zampetakis in https://github.com/MichZampetakis/IBS_for_Xsuite. No local copies of these libraries in the current project.
- Xsuite by G. Iadarola et al in https://github.com/xsuite

----------------------

1. 000_sequences_and_XSlines: save sequences and create Xsuite lines. Examples for LHC protons and ions, SPS protons and ions, injection and top energy. Parameters might not exactly correspond to operational scenario. For the LHC, 000_pymask.py and optics_specific_tools.py can be found in https://github.com/lhcopt/lhcmask/tree/master/python_examples 

2. 001_PythonIBS_vs_MadxIBS: benchmarking of analytical approaches with MADX IBS and PyIBS. Formulas for bunch length and energy deviation are valid only for small amplitude oscillations.
    - master: contains ibs_madx.py and ibs_python.py which are the template scripts to use. ibs_python.py is based on the example of https://github.com/MichZampetakis/IBS_for_Xsuite/blob/main/CLIC_DRs_IBS_allmodels.py
    - examples: A series of examples that benchmarks the analytical approach of MADX-IBS and PyIBS. Input parameters are passed through config.yaml and sequences are read from 000_sequences_and_XSlines.

3. 002_PyIBS_Xsuite: PyIBS with Xsuite example for SPS ions. Based on CLIC example of https://github.com/MichZampetakis/IBS_for_Xsuite/blob/main/CLIC_DRs_IBS_allmodels.py
