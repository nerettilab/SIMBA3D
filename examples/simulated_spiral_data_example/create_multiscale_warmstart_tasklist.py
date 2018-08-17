# -*- coding: utf-8 -*-
"""
Created on Thu Sep 07 05:53:57 2017

@author: Michael Rosenthal
"""
import os
import numpy as np
import uuid
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

import simba3d.mp_manager as mp

# the user needs to specify the location of the data and the results here
inputdir='data/'   # root of all the data inputs
outputdir='results/' # directory where the files will be outputed

    
#% write the down sampled matrices to file

# create the downsampled population matrices
#res,population_contact_matricies=mp.save_matrices(inputdir,'population_mESC_M19_a.npy')

# create the downsampled single cell matrices
res,pairwise_contact_matricies=mp.save_matrices(inputdir,'simulate_data_400.npy')

# specify the matrix index range
number_of_resampled_datasets=len(pairwise_contact_matricies)
# use the six highest resolution resampled matrices (YOU! make sure you does get out of index range)
downsampled_matrix_index=range(number_of_resampled_datasets-4,number_of_resampled_datasets)
# manually enter which downsample matrices you want to use
# downsampled_matrix_index=[4,5,6,7]
#%% create the task list
taskss=[]
for jj in range(1):
  tasks=[]
  UUID=str(uuid.uuid4()) # generated for keeping track of individual tasks
  # ii indexes which downsampled matrices to include in the multiscale warmstart
  for ii in downsampled_matrix_index: 
      tasks.append(
          {
              # task parameters        
              'taskname'  : 'run on '+str(res[ii])+'x'+str(res[ii]),
              'uuid'  : UUID+"_"+str(ii),#  assign a unique identifier to resume tasklists
              # initialization settings        
              'randomize_initialization'  : False, # Randomize the initialized curve
              #'seed'                      : int(np.random.rand()*4294967295),   # fixed seed for randomizing inititialization
              'usewarmstarts': True, # this tell pleases not to run in parallel and the use previous iteration to start the next iteration
              # output settings
              'plotfigs'  : False,    # to plot the figures after each run
              'animate'   : False,    # to animate the iterative energy 
              
              # debug options
              'check_jacobian'    : False, # check analytical gradient with numerical gradient 9
              
               # data file names to load
              'file_names'  :
                  {
                  'inputdir'  : inputdir, # location of data input files
                  'outputdir' : outputdir, # directory where to output results
                  # the following are relative to the inputdir directory
                  "pairwise_contact_matrix"   : pairwise_contact_matricies[ii],
                  #"initialized_curve"         : 'initialization.npy',
                  #"population_contact_matrix" : population_contact_matricies[ii],   
                  #"prior_shape_model"         : 'prior_shape_model.npy',       
                  'output_filename':"simulated_data_"+UUID+'_'+str(res[ii])+'.npz' # manually assign name to output file
                  },
              'index_parameters' : # specify which matrix indices are included in optimization
                  { 
                  'offdiagonal' : 1 # which off diagonal entry to start at?
                  #'missing_rowsum_threshold' : 0 # automatically set rows that sum up to less than the threshold in the log-likelihood term at missing indexes (for stability)
                  # below are advanced index parameters
                  #'index_mississing': np.array([],dtype=int) # manually specify index of missing entries (usuall because certain regions are unobservable)         
                  #'pairwise_combinations': np.array([[0,1]],dtype=int) # manually specify which pairwise indexes will be used in the log-likelihood term 
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
                      'data'             : 1.0,      # weight for the data term
                      'uniform spacing'  : 0.1, # scaled first order penalty
                      'smoothing'        : 0.1, # scaled second order penalty
                       #  'population prior' : 0.0,    # weight for the population matrix prior
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
                  #'gradient tolerance':1e-3,
                  'maxitr'    : 100000, # set maximum number of iterations
                  'display'   : True, # display function values at each iteration
                  'store'     : False, # store iterative curves
                  #'method'    : 'default'
                  'method'    :'BFGS',# Broyden, Fletcher, Goldfarb, and Shanno (quasi-Newton, only first derivatives are used)
                  
                  # below are alternative blackbox standard optimization algorithms 
                  # (Make sure you comment out BFGS if you want to use one of these)
                  #'method'    :'Nelder-Mead',# robust but slow numerical approach
                  #'method'    :'CG',# nonlinear conjugate gradient algorithm by Polak and Ribiere (only first derivatives are used)
                  #'method'    :'Newton-CG',# (truncated Newton method)
                  },
              }    
      )
  taskss.append(tasks)        
#%%
# create index



#%% Use JSON to write the tasklist to file
taskfilepath='simulated_data_multiscale_warmstart_tasklist.json'
with open(taskfilepath, 'w') as tasklist:
    json.dump(taskss, tasklist,indent=4,sort_keys=True) # pretty print the JSON to make it human readible
    #json.dump(tasks, tasklist) # if you don't care about pretty printing the JSON just dump it this way

