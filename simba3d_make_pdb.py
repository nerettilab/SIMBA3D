'''
stand alone script for converting simba3d results into a pdb format which can be read by chimera.
'''

import sys
import numpy as np
import os
from scipy.io import loadmat
try: import simplejson as json # try to import simplejson
except ImportError: import json #otherwise import json

def make_pdb_file_00(outputfilepath,X):
    '''
    taks a 3xn array and converts it into a pdb format which can be read by chimera.
    '''
    smallest=abs(X).min()
    largest=abs(X).max()
    X*=100.0/largest
    with open(outputfilepath, 'w') as result:
        result.write("Structure Inference from Multiscale BAyes in 3 Dimensions (SIMBA3D)\n")
        for ii in range(len(X[0])):
            result.write('ATOM  ')
            result.write('{: 5d}'.format(ii+1))
            result.write('   CA MET A'+str(ii+1).ljust(8))
            result.write('{: 8.3f}'.format(X[0,ii]))
            result.write('{: 8.3f}'.format(X[1,ii]))
            result.write('{: 8.3f}'.format(X[2,ii]))
            result.write('  0.20 10.00\n')
        for ii in range(len(X[0])-1):
            result.write('CONECT')
            result.write('{: 5d}'.format(ii+1))
            result.write('{: 5d}'.format(ii+2)+'\n')

def printhelp_00():
    """
    print the help text for the convertion tool
    """
    print('Convert simba3d report output from npz to .mat, or .json')
    print('simba3d-make-pdb <list of files>')


def load_result_00(outputfile):
    """
    This will load a simba3d result
    """
    results_dir='.'
    extension = os.path.splitext(outputfile)[1]
    summary={}
    try:

        if extension == '.npy':
            data=np.load( os.path.join(results_dir,str(outputfile)))
            if 'summary' in data:
              summary=data['summary'].item()
            else:
              summary={'X_evol':[data]}
        if extension == '.npz':
            data=np.load( os.path.join(results_dir,str(outputfile)))
            summary=data['summary'].item()
        if extension == '.mat':
            summary=loadmat(os.path.join(results_dir,str(outputfile)))
        if extension == '.json':
            with open(os.path.join(results_dir,str(outputfile)),'r') as result:
                summary=dict(json.load(result))
            summary['X_evol']=[np.array(X) for X in summary['X_evol']]
            summary['initialized_curve']=np.array(summary['initialized_curve'])
    except:
        summary={}
    return summary

def save_data_00(outputfilepath,summary):
    """
    This will save a simba3d result

    Summary is the result of the experiment
    """
    ext=os.path.splitext(outputfilepath)[-1].lower()
    if ext=='.pdb':
        make_pdb_file_00(outputfilepath,np.array(summary['X_evol'][-1]))

def convert_00(filename,ext_out=".json"):
    """
    convert .npz simba3d result into some other supported format
    """
    summary=[]
    try:
        ext_in = os.path.splitext(filename)[1]
        if (ext_in == ext_out):
            print("\tError input extention matches output extention")
        else:
            summary=load_result_00(filename)
    except:
            print("\tError loading "+filename)
    if summary:
        try:
            print("\t saving :"+os.path.splitext(filename)[0]+ext_out)
            save_data_00(os.path.splitext(filename)[0]+ext_out,summary)
        except:
            print("\tError converting "+filename+" to "+ext_out)
    else:
        print("\tCould not load "+filename)

def main(args=None):
    """
    main function executing the convertion tool
    """
    if args is None:
        args=sys.argv[:]
    if ('--help' in args) | ('-h'in args): # check if the help option was called
        printhelp_00()
        return 0; # do not continue if the help was called
    ext_out='.pdb' #
    del args[0] # the first arg is just the function call
    for filename in args: # loop through all the arguments and convert them
        print("Converting "+filename+" to "+ext_out)
        convert_00(filename,ext_out)

if __name__ == "__main__":
   main()
