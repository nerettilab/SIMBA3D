# -*- coding: utf-8 -*-
"""
A simple tool for creating initalized curves for SIMBA3D

Currently it only uses gaussain random noise, but in the future more
initialization methods may be implemented.

@author: Michael Rosenthal
"""

import os
import sys
import numpy as np
from scipy.io import loadmat, savemat
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from difflib import SequenceMatcher
try: import simplejson as json # try to import simplejson
except ImportError: import json #otherwise import json

import simba3d.plotting_tools as pt
import simba3d.srvf_open_curve_Rn as srvf
import simba3d.mp_manager as mp
from simba3d.matrixlabish import significant_figures,keyboard
from simba3d.latex_reports import make_latex_report_header,make_latex_table,make_latex_report_footer

#%%
def printhelp():
    """
    print the document options for the simba3d command line utility
    """
    print('Create a random initialization')
    print('[options] -o <result files>')
    print('\t-o or --output-files  <result_file_name.csv> csv filename')
    print('[Options]')
    print('-n or --number-of-nodes <number of nodes>')
    print('-m or --method <gaussian> specify the method of initialization (currently only has gaussian random noise')
    print('-s or --seed <a non-negative integer> you can specify the seed if you want to')
    #print '-f or --filter <string parameter name> <minimum value> <maximum value>'

def initialize_gaussian(n,d,rng=None):
    if rng is None:
        return np.random.rand(n*d).reshape([d,n])
    else:
        return rng.normal(0,1,n*d).reshape([d,n])

# add your method here
def initialize_template(n,d,rng=None):
    '''
    This is just a placeholder for adding future methods
    '''
    if rng is None:
        '''
        If the random number generator is none, then the user did not set the
        seed.
        '''
        return np.random.rand(n*d).reshape([d,n])
    else:
        '''
        If the random number generator is not none, then the user did not set
        the seed.
        '''
        return rng.normal(0,1,n*d).reshape([d,n])


def main(args=None):
    """
    main function for the simba3d display command line utility
    """
    if args is None:
       args=sys.argv[:]
    ii=1;
    n=5000
    d=3
    rng=None
    method='gaussian'
    suported_outputs=['.csv','.json']
    outputfiles=[]
    while ii < len(args):
        if (args[ii]== '-h')|(args[ii]== '--help'):
           printhelp()
           sys.exit()
        if (args[ii]== '-s')|(args[ii]== '--seed'):
            ii+=1
            seed=int(args[ii])
            rng=np.random.RandomState(seed)
        if (args[ii]== '-n')|(args[ii]== '--number-of-nodes'):
            ii+=1
            n=int(args[ii])
        if (args[ii]== '-o')|(args[ii]== '--output-files'):
           outputfiles=[]
           ii+=1
           while ii < len(args):
               outputfiles.append(args[ii])
               ii+=1
        ii+=1
    #print inputfiles
    if not outputfiles:
        print('No output provided')
        printhelp()
    for outputfile in outputfiles:
        print(outputfile)
        basename=os.path.basename(outputfile)
        filepath=os.path.dirname(outputfile)
        if filepath =="":
            filepath='.'
        if not os.path.exists(filepath):
            print("file path does not exist")
            os.mkdir(filepath)
            print("Directory ", filepath, " Created ")
        ext=os.path.splitext(basename)[-1]
        if ext not in suported_outputs:
            print(ext+' is not a supported output format')
        else:
            if method is 'gaussian':
                curve=initialize_gaussian(n,d,rng)
            # add your method here
            if method is 'template':
                curve=initialize_template(n,d,rng)
            #

            # write output to file
            if ext == '.csv':
                np.savetxt(outputfile,curve.T,delimiter=',')
            if ext == '.json':
                with open(outputfile, 'w') as jsonfile:
                     json.dump(mp.jsonify({'initialized_curve':curve}),jsonfile)
                     #

if __name__ == "__main__":
   main()
