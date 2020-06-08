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

# the user needs to specify the location of the data and the results here
inputdir='../data'
matrixdir='haploid/matrices/sparse/'   # root of all the data inputs
outputdir='results/' # directory where the files will be outputed
'''
cells=['NXT-1014_sup']
'''
cells=['NXT-1014_sup',
    'NXT-1124_sup',
    'NXT-2187_sup',
    'NXT-827_sup',
    'NXT-102_sup',
    'NXT-1125_sup',
    'NXT-2194_sup',
    'NXT-828_sup',
    'NXT-1034_sup',
    'NXT-1133_sup',
    'NXT-222_sup',
    'NXT-832_sup',
    'NXT-1035_sup',
    'NXT-1136_sup',
    'NXT-232_sup',
    'NXT-833_sup',
    'NXT-1037_sup',
    'NXT-1142_sup',
    'NXT-2374_sup',
    'NXT-835_sup',
    'NXT-1038_sup',
    'NXT-1144_sup',
    'NXT-2392_sup',
    'NXT-836_sup',
    'NXT-1042_sup',
    'NXT-1147_sup',
    'NXT-2441_sup',
    'NXT-841_sup',
    'NXT-1045_sup',
    'NXT-1154_sup',
    'NXT-248_sup',
    'NXT-846_sup',
    'NXT-1046_sup',
    'NXT-1162_sup',
    'NXT-442_sup',
    'NXT-847_sup',
    'NXT-1047_sup',
    'NXT-1163_sup',
    'NXT-466_sup',
    'NXT-853_sup',
    'NXT-1048_sup',
    'NXT-1171_sup',
    'NXT-471_sup',
    'NXT-85_sup',
    'NXT-1052_sup',
    'NXT-1173_sup',
    'NXT-476_sup',
    'NXT-864_sup',
    'NXT-1057_sup',
    'NXT-1175_sup',
    'NXT-477_sup',
    'NXT-902_sup',
    'NXT-1061_sup',
    'NXT-1184_sup',
    'NXT-482_sup',
    'NXT-906_sup',
    'NXT-1063_sup',
    'NXT-1186_sup',
    'NXT-483_sup',
    'NXT-912_sup',
    'NXT-1064_sup',
    'NXT-1187_sup',
    'NXT-762_sup',
    'NXT-915_sup',
    'NXT-1066_sup',
    'NXT-1192_sup',
    'NXT-766_sup',
    'NXT-924_sup',
    'NXT-1071_sup',
    'NXT-1193_sup',
    'NXT-772_sup',
    'NXT-926_sup',
    'NXT-1072_sup',
    'NXT-1194_sup',
    'NXT-774_sup',
    'NXT-932_sup',
    'NXT-1075_sup',
    'NXT-1202_sup',
    'NXT-782_sup',
    'NXT-935_sup',
    'NXT-1076_sup',
    'NXT-1205_sup',
    'NXT-785_sup',
    'NXT-937_sup',
    'NXT-107_sup',
    'NXT-1216_sup',
    'NXT-787_sup',
    'NXT-941_sup',
    'NXT-1086_sup',
    'NXT-157_sup',
    'NXT-792_sup',
    'NXT-944_sup',
    'NXT-1091_sup',
    'NXT-2073_sup',
    'NXT-793_sup',
    'NXT-945_sup',
    'NXT-1092_sup',
    'NXT-2085_sup',
    'NXT-796_sup',
    'NXT-955_sup',
    'NXT-1102_sup',
    'NXT-2087_sup',
    'NXT-803_sup',
    'NXT-962_sup',
    'NXT-1103_sup',
    'NXT-2092_sup',
    'NXT-812_sup',
    'NXT-964_sup',
    'NXT-1106_sup',
    'NXT-2095_sup',
    'NXT-813_sup',
    'NXT-973_sup',
    'NXT-1108_sup',
    'NXT-2102_sup',
    'NXT-814_sup',
    'NXT-982_sup',
    'NXT-1113_sup',
    'NXT-2114_sup',
    'NXT-817_sup',
    'NXT-994_sup',
    'NXT-1114_sup',
    'NXT-2136_sup',
    'NXT-822_sup',
    'NXT-1117_sup',
    'NXT-2167_sup',
    'NXT-824_sup']

    
'''
chromosomes=[]
for ii in range(1,20):
    chromosomes.append(str(ii))
chromosomes.append('X')
'''
chromosomes=['5','10','19']
resolutions=['1000000','500000','100000']
#%% create the task list

l1_penalties=[10.0];
l2_penalties=[0.0];
l3_penalties=[0.1];
ind=500;
for l1 in l1_penalties:
 for l2 in l2_penalties:
  for l3 in l3_penalties:  
   for cell in cells:
    for chromosome in chromosomes:
     for treate_zeros_missing in [True]:
        first_run= True
        tasks=[]
        ind+=1
        UUID=str(ind).zfill(6)+'_'+str(uuid.uuid4())[:8] # generated for keeping track of individual tasks
        ind2=0
        for resolution in resolutions:
            if (not os.path.exists('results/'+chromosome)):
                os.mkdir('results/'+chromosome)
            inputdir0=os.path.join(inputdir,resolution+'_reindexed_contact_matrices')
            pairwise_contact_matrix=cell+'_'+resolution+'_'+chromosome+'_'+chromosome+'_contact_matrix.json'
            population_contact_matrix='population_'+chromosome+'_'+chromosome+'_contact_matrix.json'
            print(os.path.join(inputdir,pairwise_contact_matrix)+'\n')
            if (1):
                ind2+=1
                SUBUUID=UUID+"_"+str(ind2).zfill(6)
                tasks.append({
                  # task parameters
                  'taskname'  : 'run on '+resolution,
                  'uuid'  : SUBUUID,#  assign a unique identifier to resume tasklists
                  # initialization settings
                  'randomize_initialization'  : False, # Randomize the initialized curve
                  #'seed'                      : int(np.random.rand()*4294967295),   # fixed seed for randomizing inititialization
                  'usewarmstarts': True, # this tell pleases not to run in parallel and the use previous iteration to start the next iteration
                  'nonspecified_zeros_as_missing':treate_zeros_missing,
                  # debug options
                  'check_jacobian'    : False, # check analytical gradient with numerical gradient 9
                   # data 15file names to load
                  'file_names_inputdir'  : inputdir, # location of data input files
                  'file_names_outputdir' : os.path.join(outputdir,chromosome), # directory where to output results
                  # the following are relative to the inputdir directory
                  "file_names_sparse_pairwise_contact_matrix"   : os.path.join(matrixdir,resolution+'_reindexed_contact_matrices',pairwise_contact_matrix),
                  #"file_names_initialized_curve"         : 'initialized_curves/init000.csv',
                  "file_names_sparse_population_contact_matrix" : os.path.join(matrixdir,resolution+'_reindexed_contact_matrices',population_contact_matrix),
                  #"prior_shape_model"         : 'prior_shape_model.npy',
                  'file_names_output_filename':SUBUUID+"_"+resolution+"_"+cell+'_'+chromosome+'.json', # manually assign name to output file
                  # specify which matrix indices are included in optimization
                  'index_parameters_offdiagonal' : 1 ,# which off diagonal entry to start at?
                  #'index_parameters_missing_rowsum_threshold' : 0, # automatically set rows that sum up to less than the threshold in the log-likelihood term at missing indexes (for stability)
                  # below are advanced index parameters
                  #'index_parameters_index_mississing': np.array([],dtype=int), # manually specify index of missing entries (usuall because certain regions are unobservable)
                  #'index_parameters_pairwise_combinations': np.array([[0,1]],dtype=int), # manually specify which pairwise indexes will be used in the log-likelihood term
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
                  'parameters_a'              : -3.0,
                  'parameters_b'              : 1.0, # not identifiable with unconstrained scale
                  'parameters_term_weights_data'             : 1.0,      # weight for the data term
                  'parameters_term_weights_uniform spacing'  : l1, # scaled first order penalty
                  'parameters_term_weights_smoothing'        : l2, # scaled second order penalty
                  'parameters_term_weights_population prior' : l3,    # weight for the population matrix prior
                  'parameters_term_weights_shape prior'      : 0.0,    # weight for the shape prior
                  # options
                  #'options_gradient tolerance':1e-3,
                  'options_maxitr'    : 1000000, # set maximum number of iterations
                  'options_display'   : True, # display function values at each iteration
                  'options_store'     : False, # store iterative curves
                  #'options_method'    : 'default'
                  'options_method'    :'L-BFGS-B',# Limited Memory Broyden, Fletcher, Goldfarb, and Shanno (quasi-Newton, only first derivatives are used)
                  # below are alternative blackbox standard optimization algorithms
                  # (Make sure you comment out BFGS if you want to use one of these)
                  #'options_method'    :'Nelder-Mead',# robust but slow numerical approach
                  #'options_method'    :'CG',# nonlinear conjugate gradient algorithm by Polak and Ribiere (only first derivatives are used)
                  #'options_method'    :'Newton-CG',# (truncated Newton method)
                  })
                if first_run:
                      tasks[-1]['file_names_initialized_curve']='initialized_curves/init'+str(ind).zfill(4)+'.csv'
                      first_run= False
                taskfilepath='tasks/'+UUID+'_tasklist.json'
                with open(taskfilepath, 'w') as tasklist:
                    json.dump([tasks], tasklist,indent=4,sort_keys=True) # pretty print the JSON to make it human readible
                    #json.dump(tasks, tasklist) # if you don't care about pretty printing the JSON just dump it this way
