# -*- coding: utf-8 -*-
"""
Command line convertion tool for simba3d outputs

This will convert npz outputs from simba3d to some other outputs

Created on Tue Jan 16 10:04:17 2018

@author: Michael Rosenthal
"""

import sys
import numpy as np
import os
import scipy.io
from simba3d.mp_manager import convert

def printhelp():
    """
    print the help text for the convertion tool
    """
    print('Convert simba3d report output from npz to .mat, or .json')
    print('simba3d-convertion --ext_out [.mat, .json, .txt, or .pdb] <list of files>')

def main(args=None):
    """
    main function executing the convertion tool
    """
    if args is None:
        args=sys.argv[:]
    if ('--help' in args) | ('-h'in args): # check if the help option was called
        printhelp()
        return 0; # do not continue if the help was called
    ext_out='.json' # set default output extention
    if '--ext_out' in args: # check for the output extention
        ext_out=args[args.index('--ext_out')+1]
        del args[args.index('--ext_out')+1]
        del args[args.index('--ext_out')]
    del args[0] # the first arg is just the function call
    for filename in args: # loop through all the arguments and convert them
        print("Converting "+filename+" to "+ext_out)
        convert(filename,ext_out)

if __name__ == "__main__":
   main()
