- 000_pymask.py and optics_specific_tools.py can be found in https://github.com/lhcopt/lhcmask/tree/master/python_examples
- to save the sequence apart from the xsuite_line do:
    mad_track.input(f'''
    save,sequence=lhcb1,file=LHC.seq;
    ''')
