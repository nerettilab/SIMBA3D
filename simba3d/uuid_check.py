# -*- coding: utf-8 -*-
"""
Tools for checking if a file with a specific uuid already exists

Created on Tue Oct  3 08:51:56 2017

@author: Michael Rosenthal
"""

from simba3d.matrixlabish import significant_figures,keyboard
import numpy as np
import os
import scipy.io
from scipy.io import loadmat
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
        outputfilepath=os.path.join(outputdir,task['file_names']['output_filename'])
    else:
        output_string   =  task['taskname'].replace(' ','_')+"_"+task['uuid']+'.json'
        outputfilepath  = os.path.join(outputdir,output_string )
    return outputfilepath

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

def get_uuid_from_directory(results_dir):
    '''
    load all the files in the specified results directory and extract the uuid.

    This is useful for determining which files have already been run.
    '''
    files=os.listdir(results_dir)
    uuid=[]
    for filename in files:
        extension = os.path.splitext(filename)[1]
        try:
            summary=load_result(os.path.join(results_dir,str(filename)))
            if 'uuid' in summary:
                uuid.append(summary['uuid'])
        except:
            '''
            '''
    return uuid;
def get_uuid_from_tasklist(taskfilepath):
    '''
    Load a tasklist and get all the uuid

    This is useful for determining which files have already been run.
    '''
    with open(taskfilepath,'r') as tasklist:
        tasks=json.load(tasklist)
    uuid=[]
    for task in tasks:
        if 'uuid' in task:
            uuid.append(task['uuid'])
    return uuid;
def get_uuid_remaining_tasks(taskfilepath,results_dir):
    '''
    Find uuid of tasks in the taskfilepath that are not in the
    results_dir directory.
    '''
    uuid_finished_tasks=get_uuid_from_directory(results_dir)
    uuid_total_tasks=get_uuid_from_tasklist(taskfilepath)
    uuid_remaining_tasks=[]
    for uuid in uuid_total_tasks:
        if uuid not in uuid_finished_tasks:
            print( uuid)
            uuid_remaining_tasks.append(uuid)
    return uuid_remaining_tasks
def get_uuid_from_tasks(tasks):
    '''
    Find uuid of tasks in the already load list of tasks
    '''
    uuid=[]
    for task in tasks:
        if 'uuid' in task:
            uuid.append(task['uuid'])
    return uuid
def get_remaining_tasks(taskfilepath,results_dir):
    '''
    Find tasks in the taskfilepath that are not in the results_dir directory.
    '''
    uuid_finished_tasks=get_uuid_from_directory(results_dir)
    with open(taskfilepath,'r') as tasklist:
        tasks=json.load(tasklist)
    tasks_remaining=[]
    for task in tasks:
        if 'uuid' in task:
            if task['uuid'] not in uuid_finished_tasks:
                print( task['uuid'])
                tasks_remaining.append(task)
    return tasks_remaining
def check_tasks(tasks):
    '''
    Checks each task in a list of task has not already been completed
    '''
    tasks_remaining=[]
    uuid_finished_tasks=[]
    results_dirs=[]
    for task in tasks:
        
        get_outputfilepath(task);
        results_dir='.'
        if 'outputdir' in task['file_names']:
            results_dir= os.path.join(results_dir,task['file_names']['outputdir'])
        if results_dir not in results_dirs:
            results_dirs.append(results_dir)
            local_finished_tasks=get_uuid_from_directory(results_dir)
            uuid_finished_tasks.extend(local_finished_tasks)
            print( 'Searching for UUID in '+results_dir)
            print( '\tFound '+str(len(local_finished_tasks))+' completed tasks in '+results_dir)
        if 'uuid' in task:
            if task['uuid'] not in uuid_finished_tasks:
                print( task['uuid']+' imcomplete')
                tasks_remaining.append(task)
    print( 'Progress: '+str(len(tasks)-len(tasks_remaining))+' of '+str(len(tasks))+' tasks are completed')
    return tasks_remaining

def check_tasks_index0(tasks):
    '''
    Checks each task in a list of task has not already been completed and returns integer index of remaining tasks
    '''
    uuid_finished_tasks=[]
    results_dirs=[]
    ii=0;
    index=[];
    for task in tasks:
        results_dir='.'
        if 'outputdir' in task['file_names']:
            results_dir= os.path.join(results_dir,task['file_names']['outputdir'])
        if results_dir not in results_dirs:
            results_dirs.append(results_dir)
            local_finished_tasks=get_uuid_from_directory(results_dir)
            uuid_finished_tasks.extend(local_finished_tasks)
            print( 'Searching for UUID in '+results_dir)
            print( '\tFound '+str(len(local_finished_tasks))+' completed tasks in '+results_dir)
        if 'uuid' in task:
            if task['uuid'] not in uuid_finished_tasks:
                print( task['uuid']+' Not found')
                index.append(ii)
        else:
            index.append(ii)
        ii+=1;

    print( 'Progress: '+str(len(tasks)-len(index))+' of '+str(len(tasks))+' tasks are completed')
    return index

def check_tasks_index(tasks):
    '''
    Checks each task in a list of task has not already been completed and returns integer index of remaining tasks
    '''
    uuid_finished_tasks=[]
    results_dirs=[]
    ii=0;
    index=[];
    for task in tasks:
        results_dir='.'
        
        if type(task) != type(list()):
            task=[task]
        number_of_subtasks=len(task)
        number_of_subtasks_complete=0
        # check to see if they give a unique file name or uuid
        for subtask in task: 
            output=get_outputfilepath(subtask)
            if os.path.isfile(output):
                summary=load_result(output)
                '''
                If the file exists and you can load it, then say the task is
                complete. I don't want to force people to use the uuid anymore.
                It is a difficult sell and it makes it more difficult to 
                explain. I want to more gently suggest it because I think it
                will make it easier for people to search for a specific result
                later. Also, if we want to coordinate tasks detributed over a 
                network it will be really important to have that in there so
                that tasks do not get duplicated.
                '''
                number_of_subtasks_complete+=1
                
                # check for the uuid too
                if 'uuid' in subtask:
                    if 'uuid' in summary:
                        if subtask['uuid']==summary['uuid']:
                           print(subtask['uuid']+': subtask completed')
                        else:
                           print(subtask['uuid']+': uuid did not match located file')
                    else:
                        print(subtask['uuid']+': uuid did not match located file')
                else:
                    print('No uuid provided')
        print(str(number_of_subtasks_complete)+' of '+str(number_of_subtasks)+' completed')
        if (number_of_subtasks==number_of_subtasks_complete):
            '''
            the task is complete
            '''
        else:
            '''
            the task is not complete
            '''
            index.append(ii)
        ii+=1
    print( 'Progress: '+str(len(tasks)-len(index))+' of '+str(len(tasks))+' tasks are completed')
    return index

def check_tasks_indexs(tasks):
    '''
    Checks each task in a list of task has not already been completed and returns integer index of remaining tasks
    '''
    uuid_finished_tasks=[]
    results_dirs=[]
    ii=0;
    index=[];
    if type(tasks)==dict:
        tasks=[tasks]
    for task in tasks:
        results_dir='.'
        if 'outputdir' in task['file_names']:
            results_dir= os.path.join(results_dir,task['file_names']['outputdir'])
        if results_dir not in results_dirs:
            results_dirs.append(results_dir)
            local_finished_tasks=get_uuid_from_directory(results_dir)
            uuid_finished_tasks.extend(local_finished_tasks)
            #print 'Searching for UUID in '+results_dir
            #print '\tFound '+str(len(local_finished_tasks))+' completed tasks in '+results_dir
        if 'uuid' in task:
            if task['uuid'] not in uuid_finished_tasks:
                #print task['uuid']+' Not found'
                index.append(ii)
        ii+=1;
    #print 'Progress: '+str(len(tasks)-len(index))+' of '+str(len(tasks))+' tasks are completed'
    return index
def check_status(tasks):
    '''
    Checks each task in a list of task has not already been completed.

    returns the number of task finish and the total number of task for each
    subtask
    '''
    uuid_finished_tasks=[]
    results_dirs=[]
    if type(tasks)==dict:
        tasks=[tasks]
    lengths=np.zeros(len(tasks))
    individual_summary=np.zeros((2,len(tasks)))
    for ii in range(len(tasks)):
        if type(tasks[ii])==dict:
            n=1
        else:
            n=len(tasks[ii])
        index=check_tasks_indexs(tasks[ii])
        k=len(index)
        individual_summary[0,ii]=n-k
        individual_summary[1,ii]=n
    individual_percentages=individual_summary[0,:]/individual_summary[1,:]
    total_summary=np.sum(individual_summary,axis=1)
    total_percentage=total_summary[0]/total_summary[1]
    return individual_summary,individual_percentages,total_summary,total_percentage
#check_status(tasks)
#%%
#get_uuid_remaining_tasks('tasklist.txt','results')
#get_remaining_tasks('tasklist.txt','results')
