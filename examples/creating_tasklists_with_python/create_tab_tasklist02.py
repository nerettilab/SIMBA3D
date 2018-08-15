# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 10:49:52 2018

@author: star
"""

# -*- coding: utf-8 -*-
"""
An example of creating tasklists
Created on Thu Sep 07 05:53:57 2017

@author: Michael Rosenthal
"""
import os
import numpy as np
import uuid
import matplotlib.pyplot as plt
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

from simba3d.mp_manager import significant_figures

from simba3d.plotting_tools import matshow
import simba3d.matrixlabish as ml
taskfilepath='tab_tasklist.txt'


tasks=[]
lambda_1=0.5; 
lambda_2=0.5;
lambda_3=0.0;

inputdir='data/'
data_matrix_files=['simulated_contact_matrix.npy']
tab_indexes=[
    [
    [0,35],
    [34,45],
    [44,80]
    ]
]
plt.close('all')
subdata_matrix_files=[]

for data_id in range(len(data_matrix_files)):
    data=np.load(os.path.join(inputdir,datamatrix[data_id]))
    n,n=shape(data)
    for ii in range(len(tab_indexes[data_id])):
        b1=tab_indexes[data_id][ii][0]
        b2=tab_indexes[data_id][ii][1]
        
        submatrix=data[range(b1,b2),:]
        submatrix=submatrix[:,range(b1,b2)]
        fig2=plt.figure();fig2.clf()
        matshow(submatrix,fig2)
        
        filename=data_matrix_files[data_id].split(".")[0]+"_"+str(b1)+"_"+str(b2)+'.npy'
        subdata_matrix_files.append(filename)        
        filename=os.path.join(inputdir,filename)
        np.save(filename,submatrix)
        print filename
#%% graphical aid in selecting tabs
if 0:
    data_id=0
    data=np.load(os.path.join(inputdir,datamatrix[data_id]))
    n,n=shape(data)
    fig=plt.figure(1)
    fig.clf()
    matshow(data,fig)
    clicked=ml.getpts(fig)
    #%%
    values=clicked.x
    values.extend(clicked.y)
    b1=int(floor(min(values)))
    b2=int(ceil(max(values)))

    fig2=plt.figure(2);fig2.clf()
    submatrix=data[range(b1,b2),:]
    submatrix=submatrix[:,range(b1,b2)]
    matshow(submatrix,fig2)

    print "["+str(b1)+","+str(b2)+"]"
#%%

for ii in range(len(subdata_matrix_files)):
    UUID=uuid.uuid4()
    tasks.append(
        {
            # task parameters        
            'usewarmstarts':True,
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
                'inputdir'  : inputdir, # location of data input files
                'outputdir' : 'results/', # directory where to output results
                # the following are relative to the inputdir directory
                "initialized_curve"         : 'initialized_curve.npy',
                "pairwise_contact_matrix"   : subdata_matrix_files[ii],
                #"population_contact_matrix" : 'population_mESC_M19_a.npy',   
                #"prior_shape_model"         : 'prior_shape_model.npy',           
                "output_filename"           :  str(UUID)+'.npz',
                },
            
            # data                
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
                    'uniform spacing'  : lambda_1, # scaled first order penalty
                    'smoothing'        : lambda_2, # scaled second order penalty
                    'population prior' : lambda_3,    # weight for the population matrix prior
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
                'method'    : 'BFGS',# Broyden, Fletcher, Goldfarb, and Shanno (quasi-Newton, only first derivatives are used)
                #'method'    :'Nelder-Mead',# robust but slow numerical approach
                #'method'    :'CG',# nonlinear conjugate gradient algorithm by Polak and Ribiere (only first derivatives are used)
                #'method'    :'Newton-CG',# (truncated Newton method)
                'gradient tolerance'    : 1e-5,
                },
            }    
    )
#Use JSON to write the tasklist to file

with open(taskfilepath, 'w') as tasklist:
    json.dump(tasks, tasklist,indent=4,sort_keys=True) # pretty print the JSON to make it human readible
    #json.dump(tasks, tasklist) # if you don't care about pretty printing the JSON just dump it this way
#%%
