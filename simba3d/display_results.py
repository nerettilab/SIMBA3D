# -*- coding: utf-8 -*-
"""
A quick graphical display of a simba3d .npz result

Created on Thu Sep 14 14:40:13 2017

@author: Michael Rosenthal
"""

import os
import sys
import numpy as np
from scipy.io import loadmat, savemat
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import simba3d.plotting_tools as pt
import simba3d.srvf_open_curve_Rn as srvf
import simba3d.mp_manager as mp
from simba3d.matrixlabish import significant_figures,keyboard


#%%
def printhelp():
    """
    print the document options for the simba3d display command line utility
    """
    print 'Graphically display results collected by simba3d'
    print 'simba3d-disp [options] -i <result files>'
    print '\t-i or --input-files  <result files> '
    print '[Options]'
    print '-c or --center-to <result filename to center plots to>'
    #print '-f or --filter <string parameter name> <minimum value> <maximum value>'
def result_passes_filter(summary,param_name,param_min,param_max):
    """
    This will check if the parameter is between the two specified values
    """
    passes_filter=True
    for ii in range(len(param_name)):
        if param_name[ii] in summary:
            passes_filter*=summary[param_name[ii]]>=param_min[ii]
            passes_filter*=summary[param_name[ii]]<=param_max[ii]
    return passes_filter
def main(args=None):
    """
    main function for the simba3d display command line utility
    """
    if args is None:
       args=sys.argv[:]
    ii=1;
    param_name=[]
    param_min=[]
    param_max=[]
    center_to_filename=None;
    while ii < len(args):
        if (args[ii]== '-h')|(args[ii]== '--help'):
           printhelp()
           sys.exit()
        if (args[ii]== '-f')|(args[ii]== '--filter'):
            ii+=1
            param_name.append(str(args[ii]))
            ii+=1
            param_min.append(np.float(args[ii]))
            ii+=1
            param_max.append(np.float(args[ii])   )
        if (args[ii]== '-c')|(args[ii]== '--center-to'):
            ii+=1
            center_to_filename=str(args[ii])
        if (args[ii]== '-i')|(args[ii]== '--input-files'):
           inputfiles=[]
           ii+=1
           while ii < len(args):
               inputfiles.append(args[ii])
               ii+=1

        ii+=1
    #print inputfiles
    curves=[]
    scurves=[]
    scurves_init=[]
    length=[]
    energy=[]

    # specify the file all the other figures will rotate to
    if center_to_filename is None:
        center_to_filename=inputfiles[0]
    summary=mp.load_result(center_to_filename)
    X0=np.array(summary['X_evol'][-1]) # get the last curve
    # get dimension
    if 'd' in summary:
        d=summary['d']
    else:
        d,n0=np.shape(X0)
    if 'n' in summary:
        n=summary['n']
    else:
        d,n=np.shape(X0)
    X0=X0.reshape((d,n)) # correct dimensions
    center_curve,mu=srvf.center_curve(X0)
    scenter_curve,scale=srvf.scale_curve(center_curve)

    for inputfile in inputfiles:
        summary=mp.load_result(inputfile)
        if result_passes_filter(summary,param_name,param_min,param_max):
            #print summary.keys()
            if 'uuid' in summary:
              print summary['uuid']
            if 'E_evol'in summary:
                energy.append(summary['E_evol'])
                #print summary.keys()
            if 'X_evol'in summary:
                X0=np.array(summary['X_evol'][-1]) # get the last curve
                # get dimension
                if 'd' in summary:
                    d=summary['d']
                else:
                    d,n=np.shape(X0)
                if 'n' in summary:
                    n=summary['n']
                else:
                    d,n=np.shape(X0)
                X0=X0.reshape((d,n)) # correct dimensions
                curves.append(X0)
                length.append(n)

                curve,mu=srvf.center_curve(X0)
                scurve,scale=srvf.scale_curve(curve)
                scurves.append(scurve)
                scurve,rot=srvf.find_best_rotation(scenter_curve,scurve)
                scurves[-1]=scurve
            if 'initialized_curve'in summary:
                X0=np.array(summary['initialized_curve']) # get the last curve
                #keyboard()
                # get dimension
                if 'd' in summary:
                    d=summary['d']
                else:
                    d,n=np.shape(X0)
                if 'n' in summary:
                    n=summary['n']
                else:
                    d,n=np.shape(X0)
                X0=X0.reshape((d,n)) # correct dimensions
                curves.append(X0)
                length.append(n)

                curve,mu=srvf.center_curve(X0)
                scurve,scale=srvf.scale_curve(curve)
                scurves_init.append(scurve)
                scurve,rot=srvf.find_best_rotation(scenter_curve,scurve)
                scurves_init[-1]=scurve


    #plt.subplots_adjust(left=0.0,right=1.0,bottom=0.0,top=1.0,wspace=0.0,hspace=0.0)
    if energy:
        fig1=plt.figure()
        for ii in range(len(energy)):
            plt.figure(fig1.number)
            plt.plot(energy[ii])
        plt.title("Energy Evolution")

        fig3=plt.figure()
        plt.subplots_adjust(left=0.0,right=1.0,bottom=0.0,top=1.0,wspace=0.0,hspace=0.0)
        #fig2.tight_layout()
        ax3=fig3.add_subplot(111,projection='3d')
        ax3.set_axis_off()
        (m,n)= np.shape(scurves_init[-1])
        t=np.linspace(0,1,n)
        pt.plot_curves(scurves_init,t,ax=ax3,fig=fig3)
        #pt.plot_pts(summary['X_evol'][-1],t,ax=ax2,fig=fig2)
        #print summary.keys()
        for ii in range(len(scurves_init)):
            pt.plot_pts(scurves_init[ii],t,ax=ax3,fig=fig3)
        plt.title("Initialized Curve")
        plt.figure(fig3.number)

        fig2=plt.figure()
        plt.subplots_adjust(left=0.0,right=1.0,bottom=0.0,top=1.0,wspace=0.0,hspace=0.0)
        #fig2.tight_layout()
        ax2=fig2.add_subplot(111,projection='3d')
        ax2.set_axis_off()
        (m,n)= np.shape(scurves[-1])
        t=np.linspace(0,1,n)
        pt.plot_curves(scurves,t,ax=ax2,fig=fig2)
        #pt.plot_pts(summary['X_evol'][-1],t,ax=ax2,fig=fig2)
        #print summary.keys()
        for ii in range(len(scurves)):
            pt.plot_pts(scurves[ii],t,ax=ax2,fig=fig2)
        plt.figure(fig2.number)
        plt.title("Estimated Curve")



        plt.show(block=True)

    else:
        print "No file specified matches filter requirements"
if __name__ == "__main__":
   main()
