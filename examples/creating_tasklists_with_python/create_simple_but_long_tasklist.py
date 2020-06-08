# -*- coding: utf-8 -*-
"""
An example of creating tasklists

@author: Michael Rosenthal
"""
import os
import numpy as np
import uuid
try: import simplejson as json # try to import simplejson
except ImportError: import json #otherwise import json

from simba3d.mp_manager import significant_figures

taskfilepath="simple_but_long_tasklist.json"

tasks=[]

lambda_2=0.5;
lambda_3=0.5;
UUID0=str(uuid.uuid4()) # generate a unique id
iii=0
for lambda_1 in np.linspace(0,1,11):
	UUID=UUID0+"_"+str(iii).zfill(5) # unique id for each subtask
	tasks.append(
		{
		    # task parameters        
		    "usewarmstarts":False,
		    "taskname"  : str(UUID),
		    "uuid"  : str(UUID),#  assign a unique identifier to resume tasklists
		    # initialization settings        
		    "randomize_initialization"  : False, # Randomize the initialized curve
		    #"seed"                      : int(np.random.rand()*4294967295),   # fixed seed for randomizing inititialization                    
		    # debug options
		    "check_jacobian"    : False, # check analytical gradient with numerical gradient 9           
		     # data file names to load
		    "file_names_inputdir"  : "data/", # location of data input files
		    "file_names_outputdir" : "results/", # directory where to output results
		     # the following are relative to the inputdir directory
		    "file_names_initialized_curve"         : "initialized_curve.npy",
		    "file_names_pairwise_contact_matrix"   : "simulated_contact_matrix.npy",         
		    "file_names_output_filename"           :  str(UUID)+".json",
		    # parameter settings
		    "parameters_a"              : -3.0,
		    "parameters_b"              : 1.0, # not identifiable with unconstrained scale
		    "parameters_term_weights_data"             : 1.0,      # weight for the data term
		    "parameters_term_weights_uniform_spacing"  : lambda_1, # scaled first order penalty
		    "parameters_term_weights_smoothing"        : lambda_2, # scaled second order penalty
		    "parameters_term_weights_population_prior" : lambda_3,    # weight for the population matrix prior
		    "parameters_term_weights_shape_prior"      : 0.0,    # weight for the shape_prior
		    # options
		    "options_maxitr"    : 100000, # set maximum number of iterations
		    "options_display"   : True, # display function values at each iteration
		    "options_store"     : False, # store iterative curves
		    "options_method"    : "BFGS",# Broyden, Fletcher, Goldfarb, and Shanno (quasi-Newton, only first derivatives are used)
		    #"options_method"    :"Nelder-Mead",# robust but slow numerical approach
		    #"options_method"    :"CG",# nonlinear conjugate gradient algorithm by Polak and Ribiere (only first derivatives are used)
		    #"options_method"    :"Newton-CG",# (truncated Newton method)
		    "options_gradient tolerance"    : 1e-5,
	    }    
	)
	iii+=1
#Use JSON to write the tasklist to file

with open(taskfilepath, "w") as tasklist:
    json.dump(tasks, tasklist,indent=4,sort_keys=True) # pretty print the JSON to make it human readible
    #json.dump(tasks, tasklist) # if you don"t care about pretty printing the JSON just dump it this way
#%%
