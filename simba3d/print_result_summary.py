# -*- coding: utf-8 -*-
"""
command line tool to print parts of a simba3d output report file

Created on Tue Jan 16 13:11:09 2018

@author: Michael Rosenthal
"""


import sys
import numpy as np
import os
import scipy.io
from simba3d.mp_manager import load_result

def printhelp():
    """
    print the document options for the simba3d print utility
    """
    print('Print specific keys from .mat or .npz simba3d report file')
    print('simba3d-print -i path/to/file [keys]')

def main(args=None):
    """
    main function call for the simba3d print utility
    """
    if args is None:
        args=sys.argv[:]
    if ('--help' in args) | ('-h'in args):
        printhelp()
        return 0;
    filename=None;
    if '-i' in args:
        filename=args[args.index('-i')+1]
        del args[args.index('-i')+1]
        del args[args.index('-i')]
    del args[0]
    summary=[]
    if filename is None:
        printhelp()
    else:
        try:
            summary=load_result(filename)
        except:
            print("\tError loading "+filename)
        if args:
            for key in args:
                if key in summary:
                    print(key+"="+str(summary[key]))
        else:
            print('keys='+str(summary.keys()))
if __name__ == "__main__":
   main()
