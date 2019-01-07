# -*- coding: utf-8 -*-
"""
This contains the command line interface for simba3d.

Created on Thu Sep 07 06:48:42 2017

@author: Michael Rosenthal
"""

# import simplejson or json
#
# This is useful because the syntax is pythonic, it's human readable, and it
# can be read by other programs in other languages.
#
# json is potentially and older version of simplejson which is updated more
# frequently than Python
#
# A good practice is to use one or the other as a fallback.
try: import simplejson as json # try to import simplejson
except ImportError: import json #otherwise import json
import sys
import os
import time
import simba3d.mp_manager as mp
import uuid
import numpy as np

from simba3d.matrixlabish import significant_figures,keyboard

#print('Number of arguments:', len(sys.argv), 'arguments.'
#print('Argument List:', str(sys.argv)

#%%
def printhelp():
    '''
    print(document options for the simba3d taskrun command
    '''
    print('To apply multiprocessing to a list of tasks specify # of cores and a json tasklist')
    print('\t-r <json_tasklist_file>')
    print('')
    print('Alternatively, to specify a single run, specify file and parameter options:')
    filename_option_descriptions={
            '--inputdir':
                'local dircetory for data input files (for simplier relative paths)',
            '--outputdir'  :
                'local dircetory for data output files (for simplier relative paths)',
            '--cm_file or --pairwise_contact_matrix':
                'contact matrix file (currently supports mat and npy files)',
            '--population_contact_matrix':
                'popullation data contact matrix file (currently supports mat and npy files)',
            '--init_structure_file --initialized_curve':
                'initialized 3D structure file (currently supports mat and npy files)',
            '--o or --output_filename':
                'The output filename',
            '--override-output_filename':
                'overide the output file if the file already exists (otherwise it will skip it'
            }
    parameter_option_descriptions={
            '--lambda_1_value or --term_weight_uniform_spacing <value>':
                'set value weight for the uniform spacing penalty',
            '--lambda_2_value or --term_weight_smoothing <value>':
                'set value weight for the uniform smoothing penalty',
            '--lambda_3_value or --term_weight_population_prior <value>':
                'set value weight for the population prior penalty'
                }
    print('File and Directory options:')
    for option in filename_option_descriptions.keys():
        print('\t'+option)
        print('\t\t'+filename_option_descriptions[option])
    print('Parameter options:')
    for option in parameter_option_descriptions.keys():
        print('\t'+option)
        print('\t\t'+parameter_option_descriptions[option]   )
#%%
def main(args=None):
   """
   main simba3d task run call
   """
   if args is None:
       args=sys.argv[:]
       #print(args
   inputfile = None
   cores=1;
   ii=1;
   override_output_filename=False
   UUID=uuid.uuid4()
   task={
            # task parameters
            'taskname'  : str(UUID),
            'uuid'  : str(UUID),#  assign a unique identifier to resume tasklists
            # initialization settings
            'randomize_initialization'  : False, # Randomize the initialized curve
            #'seed'                      : int(np.random.rand()*4294967295),   # fixed seed for randomizing inititialization

            # debug options
            'check_jacobian'    : False, # check analytical gradient with numerical gradient 9

             # data file names to load
            'file_names'  :
                {
                'inputdir'  : './', # location of data input files
                'outputdir' : './', # directory where to output results
                # the following are relative to the inputdir directory
                #"initialized_curve"         : '',
                #"pairwise_contact_matrix"   : '',
                #"population_contact_matrix" : '',
                #"prior_shape_model"         : 'prior_shape_model.npy',
                "output_filename"           :  str(UUID),
                },

            # data
            # if all the tasks use the same data, it may be better to load it once
            'data'  :
                {
                #'initialized_curve'         : initialized_curve,
                #'t'                         : t,   # optional
                #'pairwise_contact_matrix'   : [],
                #'index'                     : index,
                #'population_contact_matrix' : pairwise_contact_pop,
                #'prior_shape_model'         : prior_shape_model,
                },


            # parameter settings
            'parameters'    :
                {
                'a'              : -3.0,
                'b'              : 1.0, # not identifiable with unconstrained scale
                'term_weights'   :
                    {
                    'data'             : 1.0e0,      # weight for the data term
                    'uniform spacing'  : 0.0, # scaled first order penalty
                    'smoothing'        : 0.0, # scaled second order penalty
                    'population prior' : 0.0,    # weight for the population matrix prior
                    'shape prior'      : 0.0,    # weight for the shape prior

                    # below are unsupported penalties
                    #'firstroughness'   :0.0e3,   # weight for the fist order roughness
                    #'secondroughness'  :0.0e-16,   # weight for the second order roughness
                    #'scaledfirstroughness'  :0.0,   # weight for the scaled second order roughness
                    #'scaledsecondroughness' :0.0,   # weight for the scaled second order roughness
                    #'parameterization' :0, # not implemented
                    },
                },
            # options
            'options'   :
                {
                'maxitr'    : 100000, # set maximum number of iterations
                'display'   : True, # display function values at each iteration
                'store'     : False, # store iterative curves
                'method'    :'BFGS',# Broyden, Fletcher, Goldfarb, and Shanno (quasi-Newton, only first derivatives are used)
                #'method'    :'Nelder-Mead',# robust but slow numerical approach
                #'method'    :'CG',# nonlinear conjugate gradient algorithm by Polak and Ribiere (only first derivatives are used)
                #'method'    :'Newton-CG',# (truncated Newton method)
                },
            }
   if len(args)==1:
       printhelp()
       sys.exit()
   while ii < len(args):
       if args[ii]== '-h':
           printhelp()
           sys.exit()
       elif args[ii]== '-c':
           ii+=1
           #cores=int(args[ii] )
       elif args[ii]== '-r':
           ii+=1
           inputfiles=[args[ii] for ii in range(ii,len(args))]
       #FILE SETTINGS
       elif args[ii]== '--inputdir':
           ii+=1
           task['file_names']['inputdir']=args[ii]
       elif args[ii]== '--outputdir':
           ii+=1
           task['file_names']['outputdir']=args[ii]
       elif args[ii]== '--cm_file':
           ii+=1
           task['file_names']['pairwise_contact_matrix']=args[ii]
       elif args[ii]== '--pairwise_contact_matrix':
           ii+=1
           task['file_names']['pairwise_contact_matrix']=args[ii]
       elif args[ii]== '--population_contact_matrix':
           ii+=1
           task['file_names']['population_contact_matrix']=args[ii]
       elif args[ii]== '--init_structure_file':
           ii+=1
           task['file_names']['initialized_curve']=args[ii]
       elif args[ii]== '--initialized_curve':
           ii+=1
           task['file_names']['initialized_curve']=args[ii]
       elif args[ii]== '--o':
           ii+=1
           task['file_names']['output_filename']=args[ii]
       elif args[ii]== '--output_filename':
           ii+=1
           task['file_names']['output_filename']=args[ii]
       elif args[ii]=='--override-output_filename':
           override_output_filename=True
       #PARAMETER SETTINGS
       elif args[ii]== '--param_a':
           ii+=1
           task['parameters']['a']=np.float(args[ii]  )
       elif args[ii]== '--param_b':
           ii+=1
           task['parameters']['b']=np.float(args[ii]  )
       elif args[ii]== '--term_weight_data':
           ii+=1
           task['parameters']['term_weights']['data']=np.float(args[ii] )
       elif args[ii]== '--lambda_1_value':
           ii+=1
           task['parameters']['term_weights']['uniform spacing']=np.float(args[ii] )
       elif args[ii]== '--term_weight_uniform_spacing':
           ii+=1
           task['parameters']['term_weights']['uniform spacing']=np.float(args[ii])
       elif args[ii]== '--lambda_2_value':
           ii+=1
           task['parameters']['term_weights']['smoothing']=np.float(args[ii] )
       elif args[ii]== '--term_weight_smoothing':
           ii+=1
           task['parameters']['term_weights']['smoothing']=np.float(args[ii] )
       elif args[ii]== '--lambda_3_value':
           ii+=1
           task['parameters']['term_weights']['population prior']=np.float(args[ii])
       elif args[ii]== '--term_weight_population_prior':
           ii+=1
           task['parameters']['term_weights']['population prior']=np.float(args[ii]   )
       ii+=1
   print('Input file is ', inputfiles)
   print('number of cores', str(cores))
   if not inputfiles:
       # check that the file already exist

       if override_output_filename:
           mp.run_task(task)
       else:
           outputfilepath=mp.get_outputfilepath(task)
           ext=os.path.splitext(outputfilepath)[-1].lower()
           if (ext != '.npz')|(ext=='.mat')|(ext=='.json'):
               outputfilepath+='.npz'
           if not mp.load_result(outputfilepath):
               mp.run_task(task)
           else:
               print("Skipping task with same output file")
   else:
       for inputfile in inputfiles:
        with open(inputfile,'r') as tasklist:
         tasks=json.load(tasklist)
         t= time.time()
         mp.mp_worker(tasks)
         elapsed=time.time()-t
         print('Total '+str(elapsed)+' seconds elapsed')

if __name__ == "__main__":
   main()
