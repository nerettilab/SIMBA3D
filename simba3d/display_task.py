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
from scipy.sparse import coo_matrix
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from difflib import SequenceMatcher
import uuid

from simba3d.mp_manager import load_data
from simba3d.uuid_check import check_tasks_index, load_result,get_outputfilepath
import simba3d.plotting_tools as pt
import simba3d.srvf_open_curve_Rn as srvf
import simba3d.mp_manager as mp
from simba3d.matrixlabish import significant_figures,keyboard
import simba3d.latex_reports as lr

try: import simplejson as json # try to import simplejson
except ImportError: import json #otherwise import json

#%%
def printhelp():
    """
    print the document options for the simba3d display command line utility
    """
    print('Graphically display input matrices read by simba3d')
    print('simba3d-disp [options] -i <result files>')
    print('\t-i or --input-files  <task_file_name.json> ')
    print('[Options]')
    print('-c or --colormap <colormap name> pink, jet, bone, winter, summer, ...')
    print('-d or --dense flag which indicates a dense matrix input is provided')
    print('-o of --output-directory <output/directory> specify where to create the summary report ( ./ is the default)' )
    print('-f or --format <format as png,eps, or jpg> format for outputed images (png is the default)')
    print('-p or --print-each <report name> print each result individually')
    print('-q or --quantile-color-scale <value> value is an integer between 0 and 100 specifying the upper threshold quantile for the color limit')

    #print '-f or --filter <string parameter name> <minimum value> <maximum value>'

def main(args=None):
    """
    main function for the simba3d display command line utility
    """
    if args is None:
       args=sys.argv[:]
    ii=1;
    report_directory='.'
    image_ext='png'
    summary_name='tasks_summary'
    print_each=False
    latex_print=False
    colormap_name='pink'
    p=100
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
        if (args[ii]== '-c')|(args[ii]== '--colormap'):
            ii+=1
            colormap_name=str(args[ii])
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

    latex=lr.make_latex_report_header('figures')
    tasks=[]
    scenter_curve=None
    for inputfile in inputfiles:
        curves=[]
        scurves=[]
        scurves_init=[]
        length=[]
        energy=[]
        print('Reading:'+inputfile)
        basename=os.path.basename(inputfile)
        filepath=os.path.dirname(inputfile)
        latex+=r'\pagebreak' +u'\n'
        latex+=lr.make_latex_section('Task Summary',inputfile)

        with open(inputfile,'r') as tasklist:
            tasks=json.load(tasklist)
            number_of_tasks=str(len(tasks))
            number_of_tasks_remaining=str(len(check_tasks_index(tasks))) # find which tasks still need to run
            for task in tasks:
                index_remaining=check_tasks_index(task)
                number_of_subtasks=str(len(task))
                number_of_subtasks_remaining=str(len(check_tasks_index(task))) # find which tasks still need to run
                for subtask in task:
                    UUID=None
                    if 'uuid' in subtask:
                        UUID=subtask['uuid']

                    latex+=lr.make_latex_subsection('Subtask Summary',UUID)
                    if UUID is None:
                        UUID=str(uuid.uuid4())
                    matrix_files=[]
                    is_sparse=[]
                    initialized_curve_files=[]
                    data=load_data(subtask)
                    #
                    # get the file input names
                    latex+=r'\subsubsection{Inputs}' +u'\n'
                    if 'file_names' in subtask:
                        inputdir='.'
                        if 'inputdir' in subtask['file_names']:
                            inputdir=subtask['file_names']['inputdir']
                        input_parameters= '';
                        for key in list(subtask):
                            skip='data' in key
                            skip+='file_names' in key
                            if skip ==0:
                                if type(subtask[key])==type(dict()):
                                    for subkey in list(subtask[key]):
                                        if type(subtask[key][subkey])==type(dict()):
                                            for subsubkey in list(subtask[key][subkey]):
                                                name=key+' '+subkey+' '+subsubkey
                                                input_parameters+=  name+':'+str(subtask[key][subkey][subsubkey])+u'\n'
                                        else:
                                            name=key+' '+subkey
                                            input_parameters+=  name+':'+str(subtask[key][subkey])+u'\n'
                                else:
                                   name=key
                                   input_parameters+= name+':'+str(subtask[key])+u'\n'

                        latex+=lr.make_latex_table({'inputfile':input_parameters})
                        for key in list(subtask['file_names']):
                            if 'output' not in key:
                                params=dict()
                                params['inputfile']=key+':'+subtask['file_names'][key]
                                params['table width']=2
                                if key in data:
                                    filename=UUID+key
                                    if 'initialized_curve' in key:

                                        curve,mu=srvf.center_curve(data[key])
                                        scurve,scale=srvf.scale_curve(curve)
                                        if scenter_curve is None:
                                            scenter_curve=scurve
                                        else:
                                            scurve,rot=srvf.find_best_rotation(scenter_curve,scurve)
                                        fig3=plt.figure()
                                        plt.subplots_adjust(left=0.0,right=1.0,bottom=0.0,top=1.0,wspace=0.0,hspace=0.0)
                                        #fig2.tight_layout()
                                        ax3=fig3.add_subplot(111,projection='3d')
                                        ax3.set_axis_off()
                                        (m,n)= np.shape(scurve)
                                        t=np.linspace(0,1,n)
                                        pt.plot_curve(scurve,t,ax=ax3,fig=fig3)
                                        pt.plot_pts(scurve,t,ax=ax3,fig=fig3)
                                        plt.title("Initialized Curve")
                                        if not print_each:
                                            plt.show()
                                        else:
                                            imagename=os.path.join(image_directory,filename+'.'+image_ext)
                                            fig3.savefig(imagename)

                                            params['images']=[
                                                    filename+'.'+image_ext
                                                  ]
                                            params['statistics']={
                                                    'm' : str(m),
                                                    'n' : str(n),
                                                    }
                                    if '_matrix' in key:
                                        # you have a matrix
                                        if 'sparse_' in key:
                                            # it is sparse
                                            matrix=data[key].todense()
                                        else:
                                            # it is not sparse
                                            matrix=data[key]
                                        q=np.percentile(matrix,p)
                                        plt.close('all')
                                        fig1=plt.figure()
                                        plt.figure(fig1.number)
                                        plt.imshow(matrix,colormap_name)
                                        plt.clim([0,q])
                                        plt.colorbar()
                                        if not print_each:
                                            plt.show()
                                        else:
                                            imagename=os.path.join(image_directory,filename+'.'+image_ext)
                                            fig1.savefig(imagename)


                                            (m,n)=np.shape(matrix)

                                            params['images']=[
                                                    filename+'.'+image_ext
                                                  ]
                                            params['statistics']={
                                                    'Row Length' : str(m),
                                                    'Column Length':str(n),
                                                    'Total Entries':str(n*m),
                                                    'Non-zero Entries':str(np.sum(matrix>0)),
                                                    'Sum': str(np.sum(matrix)),
                                                    'Max':str(np.max(matrix)),
                                                    'Non-trivial Rows':str(np.sum(sum(matrix>0)>0))
                                                    }
                                latex+=lr.make_latex_table(params)
                        # summarize input files\
                    latex+=r'\subsubsection{Results}' +u'\n'
                    summary=mp.load_result(get_outputfilepath(subtask))
                    if not summary:
                         latex+=u'No output result file found.\n'
                    else:
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
                            else:
                                fig1.savefig(os.path.join(image_directory,UUID+'_energies.'+ image_ext))
                                fig3.savefig(os.path.join(image_directory,UUID+'_initial_curves.'+ image_ext))
                                fig2.savefig(os.path.join(image_directory,UUID+'_estimated_curves.'+ image_ext))
                            if latex_print:
                                params=dict()
                                params['inputfile']=get_outputfilepath(subtask)
                                params['table width']=3
                                params['images']=[
                                        UUID+'_energies.'+ image_ext,
                                        UUID+'_initial_curves.'+ image_ext,
                                        UUID+'_estimated_curves.'+ image_ext
                                      ]
                                params['statistics']={
                                        'Final Energy' : energy[-1][-1],
                                        'Total Iterations':len(energy[-1]),
                                        'Total Computation Time': computation_time,
                                        }
                                latex_table=lr.make_latex_table(params)
                                latex+=latex_table

    if latex_print:
        latex+=r'\end{document}'+u'\n'
        print(latex)
        with open(os.path.join(report_directory, summary_name+ '.tex'), 'w') as result:
            result.write(latex)
if __name__ == "__main__":
   main()
