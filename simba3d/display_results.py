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
from difflib import SequenceMatcher

import simba3d.plotting_tools as pt
import simba3d.srvf_open_curve_Rn as srvf
import simba3d.mp_manager as mp
from simba3d.matrixlabish import significant_figures,keyboard
from simba3d.latex_reports import make_latex_report_header,make_latex_table,make_latex_report_footer

#%%
def printhelp():
    """
    print the document options for the simba3d display command line utility
    """
    print('Graphically display results collected by simba3d')
    print('simba3d-disp [options] -i <result files>')
    print('\t-i or --input-files  <result files> ')
    print('[Options]')
    print('-c or --center-to <result filename to center plots to>')
    print('-o of --output-directory <output/directory> specify where to create the summary report ( ./ is the default)' )
    print('-f or --format <format as png,eps, or jpg> format for outputed images (png is the default)')
    print('-p or --print-each <report name> print each result individually')

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
    report_directory='.'
    image_ext='png'
    summary_name='results_summary'
    print_each=False
    latex_print=False
    inputfiles=[]
    while ii < len(args):
        print(args[ii])
        if (args[ii]== '-h')|(args[ii]== '--help'):
           printhelp()
           sys.exit()
        if (args[ii]== '-p')|(args[ii]== '--print-each'):
            ii+=1
            summary_name = str(args[ii])
            print('\t'+args[ii])
            print_each=True
            latex_print=True
        if (args[ii]== '-o')|(args[ii]== '--output-directory'):
            ii+=1
            report_directory = str(args[ii])
            print('\t'+args[ii])
            print_each=True
            latex_print=True
        if (args[ii]== '-f')|(args[ii]== '--format'):
            ii+=1
            image_ext = str(args[ii])
            print('\t'+args[ii])
            print_each=True
            latex_print=True
        if (args[ii]== '-c')|(args[ii]== '--center-to'):
            ii+=1
            center_to_filename=str(args[ii])
            print('\t'+args[ii])
        if (args[ii]== '-i')|(args[ii]== '--input-files'):
           inputfiles=[]
           ii+=1
           while ii < len(args):
               inputfiles.append(args[ii])
               print('\t'+args[ii])
               ii+=1
        ii+=1
    image_directory = os.path.join(report_directory,'figures')
    # create the report directory if it does not exist
    if print_each:
        try:
            if not os.path.exists(report_directory):
                os.mkdir(report_directory)
                print("Directory ", report_directory, " Created ")
            if not os.path.exists(image_directory):
                os.mkdir(image_directory)
                print("Directory ", image_directory, " Created ")
        except:
            print("Potentially failed to create report directories")
    #print inputfiles
    curves=[]
    scurves=[]
    scurves_init=[]
    length=[]
    energy=[]
    if not inputfiles:
        print('No inputs provided')
        printhelp()
    else:
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

    latex=make_latex_report_header('figures')

    for inputfile in inputfiles:
        summary=mp.load_result(inputfile)
        if result_passes_filter(summary,param_name,param_min,param_max):
            if 'uuid' in summary:
              print(summary['uuid'])
            if 'E_evol'in summary:
                energy.append(summary['E_evol'])
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
            weight_uniform_spacing=None
            if "weight_uniform_spacing" in summary:
                weight_uniform_spacing=summary['weight_uniform_spacing']
            weight_smoothing=None
            if "weight_smoothing" in summary:
                weight_smoothing=summary['weight_smoothing']
            weight_population_prior=None
            if "weight_population_prior" in summary:
                weight_population_prior=summary['weight_population_prior']
            computation_time=None
            if "computation_time" in summary:
                computation_time=summary['computation_time']
            nonspecified_zeros_as_missing=None
            if "nonspecified_zeros_as_missing" in summary:
                nonspecified_zeros_as_missing=summary['nonspecified_zeros_as_missing']
            if 'initialized_curve'in summary:
                X0=np.array(summary['initialized_curve']) # get the last curve
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
                plt.close('all')
                fig1=plt.figure()
                plt.figure(fig1.number)
                plt.plot(energy[-1])
                plt.title("Energy Evolution")

                fig3=plt.figure()
                plt.subplots_adjust(left=0.0,right=1.0,bottom=0.0,top=1.0,wspace=0.0,hspace=0.0)
                #fig2.tight_layout()
                ax3=fig3.add_subplot(111,projection='3d')
                ax3.set_axis_off()
                (m,n)= np.shape(scurves_init[-1])
                t=np.linspace(0,1,n)
                pt.plot_curve(scurves_init[-1],t,ax=ax3,fig=fig3)
                pt.plot_pts(scurves_init[-1],t,ax=ax3,fig=fig3)
                plt.title("Initialized Curve")
                plt.figure(fig3.number)

                fig2=plt.figure()
                plt.subplots_adjust(left=0.0,right=1.0,bottom=0.0,top=1.0,wspace=0.0,hspace=0.0)
                #fig2.tight_layout()
                ax2=fig2.add_subplot(111,projection='3d')
                ax2.set_axis_off()
                (m,n)= np.shape(scurves[-1])
                t=np.linspace(0,1,n)
                pt.plot_curve(scurves[-1],t,ax=ax2,fig=fig2)
                pt.plot_pts(scurves[-1],t,ax=ax2,fig=fig2)
                plt.figure(fig2.number)
                plt.title("Estimated Curve")
            if not print_each:
                plt.show()
            if print_each:
                if image_ext not in ["png",'pdf','svg','ps','eps']:
                    print('invalid image format:'+image_ext)
                else:
                        base=os.path.basename(inputfile);
                        base=base.split('.')[0]
                        image_name_tmp=os.path.join(image_directory,base)
                        print(image_name_tmp)
                        try:
                            fig1.savefig(image_name_tmp+'_energies.'+ image_ext)
                            fig3.savefig(image_name_tmp+'_initial_curves.'+ image_ext)
                            fig2.savefig(image_name_tmp+'_estimated_curves.'+ image_ext)
                            if latex_print:
                                params=dict()
                                params['inputfile']=inputfile
                                params['table width']=3
                                params['images']=[
                                        base+'_energies.'+ image_ext,
                                        base+'_initial_curves.'+ image_ext,
                                        base+'_estimated_curves.'+ image_ext
                                      ]
                                params['statistics']={
                                        'Final Energy' : energy[-1][-1],
                                        'Total Iterations':len(energy[-1]),
                                        'Total Computation Time': computation_time,
                                        'Uniform Spacing Penalty':weight_uniform_spacing,
                                        'Smoothing Penalty':weight_smoothing,
                                        'Population Penalty':weight_population_prior,
                                        'Nonspecified Zeros As Missing':nonspecified_zeros_as_missing
                                        }

                                latex_table=make_latex_table(params)
                                latex+=latex_table
                        except:
                            print("unable to create image for:"+inputfile)
    if latex_print:
        latex+=make_latex_report_footer()
        print(latex)
        with open(os.path.join(report_directory, summary_name+ '.tex'), 'w') as result:
            result.write(latex)


if __name__ == "__main__":
   main()
