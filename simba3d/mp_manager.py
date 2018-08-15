# -*- coding: utf-8 -*-
"""

Simplified multi processing wrapper for simba3d

The mp_handler can be called by inputing a list of tasks and specifying the
desired number of cores. Selecting the optimal number of cores can be a bit
subjective. Generally, one uses the number of virtual cores available on the
system.

Created on Wed Sep 06 18:59:03 2017

@author: Michael Rosenthal
"""
from simba3d.matrixlabish import significant_figures,keyboard
from pprint import pprint
import os
import multiprocessing
import time
import uuid
import numpy as np
from scipy.io import loadmat, savemat
# A good practice is to use one or the other as a fallback.
try: import simplejson as json # try to import simplejson
except ImportError: import json #otherwise import json

import simba3d.optimizer as opt
#import simba3d.plotting_tools as pt
import simba3d.srvf_open_curve_Rn as srvf
import simba3d.uuid_check as uc

def create_downsampled_matrices(data):
    """
    This takes a matrix and returns a list of downsampled matrixes that recursively
    reduces the dimension of the matrix in roughly half each time

    Each time, the matrix is zero-padded to make the number of entries even, then 
    entries in the added together in groups of 4.

    It assumes the matrix is square, and that the diagonal entries are ignored.
    """
    n,m=np.shape(data)
    data[np.eye(n)==1]=0
    C=[]
    C.append(data)
    k=0
    while n>=8: # keep going until the dimension of the data is less than 8
        if np.mod(n,2)==1:
            C[k]=np.concatenate((C[k],np.zeros((1,n))),axis=0)
            C[k]=np.concatenate((C[k],np.zeros((n+1,1))),axis=1)
            n+=1
        
        n=int(n/2.0);
        C.append(np.zeros((n,n)))
        for i in range(n-1):
            for j in range(i+1,n):
                C[k+1][i,j]=C[k][2*i,2*j]
                C[k+1][i,j]+=C[k][2*i,2*j+1]
                C[k+1][i,j]+=C[k][2*i+1,2*j]
                C[k+1][i,j]+=C[k][2*i+1,2*j+1]
        A=np.triu(C[k+1],1)
        C[k+1]= C[k+1]+A.T
        k=k+1
    return C;
    
def save_matrices(inputdir,data_matrix_file,outputdir=None):
    """
    This saves a downsampled versions of a matrix npy file.

    You tell it the directory and the name of the file, and this will create the 
    down sampled matrix files for the multiresolution warm-starts.
    """
    if (outputdir==None):
      """
      if the user does not specify the outputdir, assume the ouput and input 
      directory are the same.
      """
      outputdir=inputdir
    # load the data matrix to down sample
    data=np.load(os.path.join(inputdir,data_matrix_file))
    # create the down sampled matrices
    Crev=create_downsampled_matrices(data)
    multires_filenames=[]
    res=[]
    reversed_index=list(reversed(range(len(Crev))))
    C=[]
    # loop through the downsampled matrices in reverse order
    for ii in range(len(Crev)):
        C.append(Crev[reversed_index[ii]])
        # store the resolutions for automating the tasklist generation
        res.append(len(C[ii]))
        file_name,file_extension=os.path.splitext(data_matrix_file)
        # store the file names for automating the tasklist generation
        multires_filenames.append(file_name+'_'+str(res[ii])+'x'+str(res[ii])+file_extension)
        # save the multiresolution matrices to the desired directory
        print('Saving '+multires_filenames[ii])
        np.save(os.path.join(outputdir,multires_filenames[ii]),C[ii])
    return res,multires_filenames


def jsonify(data):
    """
    This is a small serializer for saving output in human readable json format
    """
    if isinstance(data,dict):
        serialized_summary= dict() 
        for key,value in data.iteritems():
          if isinstance(value, list):
            value = [jsonify(item) for item in value]
          elif isinstance(value, list):
            value = jsonify(value)
          elif type(value).__module__=='numpy':
            value=value.tolist()
          serialized_summary[key]=value    
    elif type(data).__module__=='numpy':
        serialized_summary=data.tolist()
    else: 
        serialized_summary=data
    return serialized_summary

def load_result(outputfile):
    """
    This will load a simba3d result
    """
    results_dir='.'
    extension = os.path.splitext(outputfile)[1]
    summary={}
    try:
        
        if extension == '.npy':
            data=np.load( os.path.join(results_dir,str(outputfile)))
            summary=data['summary'].item()
        if extension == '.npz':
            data=np.load( os.path.join(results_dir,str(outputfile)))
            summary=data['summary'].item()
        if extension == '.mat':
            summary=loadmat(os.path.join(results_dir,str(outputfile)))
        if extension == '.json':
            with open(os.path.join(results_dir,str(outputfile)),'r') as result:
                summary=dict(json.load(result))
    except:
        summary={}
    return summary
def save_data(outputfilepath,summary):
    """
    This will save a simba3d result
    
    Summary is the result of the experiment
    """
    ext=os.path.splitext(outputfilepath)[-1].lower()
    if ext=='.mat':
        savemat(outputfilepath,summary)
        return outputfilepath        
    elif ext=='.json':
        with open(outputfilepath, 'w') as result:
            # serialize
            json.dump(jsonify(summary), result,indent=4,sort_keys=True)  
    elif ext=='.txt':
        with open(outputfilepath, 'w') as result:
            pprint(summary, result)      
    else:
        if ext != '.npz':
            outputfilepath+='.npz'
        np.savez(outputfilepath,summary=summary) 
    return outputfilepath         
def convert(filename,ext_out=".npz"):
    """
    convert .npz simba3d result into some other supported format
    """
    summary=[]
    try:
        ext_in = os.path.splitext(filename)[1]
        if (ext_in == ext_out):
            print "\tError input extention matches output extention"
        elif  (ext_in!='.npz'):
            print "\t Only convertion from .npz file is supported"
        else:
            summary=load_result(filename)
    except:
            print "\tError loading "+filename
    if summary:
        save_data(os.path.splitext(filename)[0]+ext_out,summary)
        try:
            print "\t saving :"+os.path.splitext(filename)[0]+ext_out
            save_data(os.path.splitext(filename)[0]+ext_out,summary)
        except:
            print "\tError converting "+filename+" to "+ext_out
    else:
        print "\tCould not load "+filename
def get_tag_name(task):
    '''
    Gets a tag name (this may go away in future versions)
    '''
    tag=''
    for penalty in task['parameters']['term_weights'].keys():
        
        f=task['parameters']['term_weights'][penalty]
        if f !=0.0:
            if penalty!='data':
                tag+='_'+penalty+''+significant_figures(f)
            elif f!=1.0:
                tag+='_'+penalty+''+significant_figures(f)
    if 'seed' in task:
        tag+='_seed'+str(task['seed'])
    return tag 
def run_tasks(tasks,resume=True):  
    '''
    plural task runner.
    
    This runs a list of independent tasks and divides the work across multiple
    processes if desired.
    
    Ther resume option tells simba3d to skip over tasks with uuid's that have
    already been ran.
    '''
    if type(tasks) is dict: # if a dicitonary is passed in, but it in a list so we know the layer
        tasks=[tasks]
    if resume:
        index_remaining=uc.check_tasks_index(tasks) # find which tasks still need to run
        print index_remaining
    else:
        index_remaining=range(len(tasks))
    usewarmstarts=False # default setting for warmstarts

    if 'uuid' in tasks[0]:
      UUID=tasks[0]["uuid"] # if the task has a uuid use that
    else:
      UUID=str(uuid.uuid4()) # generate a uuid
    for ii in index_remaining: # loop through remaining tasks
        if 'uuid' not in tasks[ii]: # make sure the task has a uuid
          tasks[ii]["uuid"]=UUID+'_'+str(ii) 
        if ii==0:
            outputfilepath=run_task(tasks[0])   # run the first task
        else:        
            if 'usewarmstarts' in tasks[ii]:
                usewarmstarts=tasks[ii]['usewarmstarts'] # if usewarmstarts parameter is set then use that value
            if usewarmstarts:
                print "Using warmstarts sequentially from previous task"
                if 'uuid' not in tasks[ii]:
                    tasks[ii]["uuid"]=UUID+'_'+str(ii) # make sure the task has a uuid
                outputfilepath=get_outputfilepath(tasks[ii-1]) # get the outputfile path from the previous run
                print outputfilepath                
                '''
                ext=os.path.splitext(outputfilepath)[-1].lower() 
                if ext=='.mat':
                    data=loadmat(outputfilepath)
                else:
                    data=np.load(outputfilepath) 
                summary=data['summary'].item()
                '''
                summary=load_result(outputfilepath)
                X0=np.array(summary['X_evol'][-1]) # get the last curve from previous run
                # get dimension of previous run
                if 'd' in summary:
                    d=summary['d']
                else:
                    d,n0=np.shape(X0)
                if 'n' in summary:
                    n0=summary['n']
                else:
                    d,n0=np.shape(X0)
                X0=X0.reshape((d,n0))
                # get dimension of next run
                data=load_data(tasks[ii])
                (m,n)=np.shape( data['pairwise_contact_matrix'])     
                # interpolate previous curve to match dimension
                t0=np.linspace(0,1,n0)
                t=np.linspace(0,1,n)
                initialized_curve=srvf.interp(t0,X0,t)
                initialized_curve=srvf.interp(t0,X0,t) # interpolate to initialize the next curve
                if 'data' in tasks[ii]:
                    tasks[ii]['data']['initialized_curve']=initialized_curve
                else:
                    tasks[ii]['data']={'initialized_curve' : initialized_curve}          
            outputfilepath=run_task(tasks[ii])        
def run_task(task):  
    """
    run a singular task
    """
    UUID=uuid.uuid4()# generate a unique task identifier
    # set defaults
    if 'uuid' not in task:
        task['uuid']= str(UUID) # a unique identifier to the specific task
    if 'taskname' not in task:
        task['taskname']= "unnamed task" # give the task a more intuitive name       
    if 'randomize_initialization' not in task:
        task['randomize_initialization']=False
    if  'check_jacobian' in task: # to animate the iterative energy
         check_jacobian= task['check_jacobian']
    else:         
        check_jacobian=False # check analytical gradient with numerical gradient
        #s= int(np.random.rand()*4294967295)
    if  'sigma' in task: # to animate the iterative energy
         sig= task['sigma']
    else:
        sig=.005;
    if 'parameters' not in task:
        task['parameters']= {
                'a'              : -3.0,
                'b'              : 1.0, # not identifiable with unconstrained scale
                    }
    if 'term_weights' not in task['parameters']:
        task['parameters']['term_weights']= {
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
                    }
    if 'store_task' not in task: # set default for store_task
        task['store_task']=True # store the task in the output file

    parameters=task['parameters']
    term_types=parameters['term_weights'].keys()        
    term_weights=[]             
    for term in term_types:        
        term_weights.append(np.array(parameters['term_weights'][term],dtype=np.double))
    if 'inputdir' not in task['file_names']:
        task['file_names']['inputdir']='./'
    inputdir=task['file_names']['inputdir']
    if 'outputdir' not in task['file_names']:
         task['file_names']['outputdir']='./'
    outputdir=task['file_names']['outputdir']   
    # check that output directory exists and create it is needed
    directory = os.path.dirname(outputdir)
    try:
        os.stat(directory)
    except:
        print "Please make sure the output directory exists"
        return 0
    outputfilepath=get_outputfilepath(task)
    

    # load data ###############################################################   
    data=load_data(task);
        
                  
    if 'initialized_curve' in data: 
        if task['randomize_initialization']:  
            # randomize initialized curve with gaussian noise
            if 'seed' in task:
                rng=np.random.RandomState(task['seed'])
                data['initialized_curve']=rng.normal(data['initialized_curve'],sig)     
            else:
                data['initialized_curve']=np.random.normal(data['initialized_curve'],sig) 
    # check that the initialized curve matches dimension of data
    
    if 'index_parameters' in task:      
        if 'missing_rowsum_threshold' in task['index_parameters']:
            # specifies a threshold matrix row sum to treat an entry as missing data (missing nodes are ignored in the optimization)
            data['missing_rowsum_threshold']=task['index_parameters']['missing_rowsum_threshold']
        if 'index_mississing' in task['index_parameters']:
            # specifies which entries are treated as missing (missing nodes are ignored in the optimization)
            data['index_mississing']=task['index_parameters']['index_mississing']    
        if 'off_diagonal' in task['index_parameters']:
            # offset off of the diagonal entries which is treated as missing (and ignored)
            data['off_diagonal']=task['index_parameters']['off_diagonal']               
        if 'pairwise_combinations' in task['index_parameters']:
            # optionally specify specific pairwise combination to use
            data['pairwise_combinations']=task['index_parameters']['pairwise_combinations']                       
    #%
    
    if 'options' in task:
        options=task['options']
    else:
        options = {
                "method": "BFGS", 
                "display": True, 
                "maxitr": 100000
                }
    
    start=time.time()
                
    result=opt.opt_E(data,parameters,options)     
    print "Running ..."
    print "\tTaskname: "+task['taskname']
    if 'uuid' in task:
        print "\tUUID: "+task['uuid']
    if 'title' in task:
        print "\ttitle: "+task['title']        
    if 'description' in task:
        print "\tdescription: "+task['description']             
    result.run()
    end=time.time()   
    comptime=  end-start  
    print "Time ellapse: "+str((comptime)/(60))+" minutes"
    if check_jacobian:
        # Check analytical jacobian is correct ########################################
        jac_error=result.check_jacobian() # this should be a small number
        # Check the correctness of a gradient function by comparing it against a 
        # (forward) finite-difference approximation of the gradient.
        print 'Error between analytical and numerical gradient at initialized curve: '+str(jac_error)
        # norm of the difference between analytical jacobian and the finite difference approximation of jacobian at initialized curve.
    if result.XK is None:
        result.XK=[result.estimated_curve]
    
    summary=summarize_results(result,task)
    
    print "saving results to file:"+outputfilepath
    return save_data(outputfilepath,summary)
def get_outputfilepath(task):
    """
    returns the output file path of a task
    """
    outputfilepath=None
    if 'outputdir' in task['file_names']:
        outputdir= task['file_names']['outputdir']
    else:
        outputdir='.'
    if 'output_filename' in task['file_names']:
        outputfilepath=outputdir+task['file_names']['output_filename']
    else:
        tag=get_tag_name(task)
        output_string   =  task['taskname'].replace(' ','_')+"_"+task['uuid']+'.npz'
        outputfilepath  = os.path.join(outputdir,output_string )   
    return outputfilepath
def summarize_results(result,task):
    """
    creates a summary of the result as dict
    """
    summary={
            #'Csc':data['pairwise_contact_matrix'],
            #'Cpop':data['population_contact_matrix'],
            #'Xmean':data['prior_shape_model'],
            'n':result.n,
            'd':result.d,
            'a':result.a,
            'b':result.b,   
            'E_evol':result.e,
            'X_evol':result.XK,
            'computation_time':result.computation_time,
            'method':result.method,
            'initialized_curve':result.initialized_curve,
            }
    if task['store_task']:
        summary['json_task']=task
    if 'uuid'in task:
        summary['uuid']=task['uuid']         
    for term in  result.term_weights.keys():
        summary['weight_'+term.replace(' ','_')]= result.term_weights[term]
        print term+' = '+str(summary['weight_'+term.replace(' ','_')])
    return summary
    
def mp_worker(task):
    if 'taskname' in task:
        taskname=task['taskname']
    else:
        taskname='task'
    #output_string=taskname+get_tag_name(task)
    print '\nProcess %s \tStarted\t %s\n' % (taskname,time.strftime('%d%b%Y %H:%M:%S',time.gmtime()))
    # place task here and pass in parameters
    t= time.time()
    run_tasks(task)
    elapsed=time.time()-t    
    print '\nProcess %s \tDone \t %s seconds elapsed' % (taskname, elapsed)
    
def load_data(task):
    '''
    
    potential problem here, this only loads numpy arrays npy files
    '''
    if 'inputdir' not in task['file_names']:
        task['file_names']['inputdir']='./'
    inputdir=task['file_names']['inputdir'] 
    datatypes=list(set(task['file_names'].keys())- set(['inputdir','outputdir','output_filename']))
    if 'data' in task:
        data=task['data'];
    else:
        data={};
    for datatype in datatypes:
        if datatype not in data: # don't load it if the user specified the data at startup
             #print "\nLoading "+ datatype+'\n'
             data[datatype]=np.load(inputdir+task['file_names'][datatype])
    return data 

def mp_handler(tasks,cores):
    """
    This handles the multi processing
    """
    p = multiprocessing.Pool(cores)
    p.map(mp_worker, tasks,1)
 #%%
