# -*- coding: utf-8 -*-
"""
Display graphical results from simba3d.opt_E object. 

These used to be in the opt_E class itself, but I separated it so that 
matplotlib is not required to run the optimizer.

Created on Thu Sep 07 11:40:43 2017

@author: Michael Rosenthal
"""

import matplotlib.pyplot as plt
import matplotlib.animation as ani
import mpl_toolkits.mplot3d.axes3d as p3

import simba3d.optimizer as opt
import simba3d.plotting_tools as pt
import simba3d.mp_manager as mp

def plot_energy_evolution(result,fig=None):
    '''
    plot the evolution
    ''' 
    if fig is None:
        fig=plt.figure();
    plt.plot(result.e)
    return fig

def plot_estimated_curve(result,ax=None,fig=None):
    '''
    plot the estimated curve stored in the instance
    '''
    if fig is None:
        fig=plt.figure();
    if ax is None:
            ax=fig.add_subplot(111,projection='3d')
    pt.plot_curve(result.estimated_curve,result.t,ax,fig)       
    return ax,fig

def plot_initialized_curve(result,ax=None,fig=None):
    '''
    plot the curve used to initiailize the gradient ascent algorithm
    '''
    if fig is None:
        fig=plt.figure();     
    if ax is None:
            ax=fig.add_subplot(111,projection='3d')         
    pt.plot_curve(result.initialized_curve,None,ax,fig)  
    return ax,fig
def create_images_from_task(task):
    '''
    Will load the individual task (if possible) and plot the relavent images
    '''
    outputfile=get_outputfilepath(task)
    if outputfile is not None:
        summary=mp.load_result(outputfile)
