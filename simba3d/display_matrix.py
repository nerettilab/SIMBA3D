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
    print('\t-i or --input-files  <result files> ')
    print('[Options]')
    print('-c or --colormap <colormap name> pink, jet, bone, winter, summer, ...')
    print('-d or --dense flag which indicates a dense matrix input is provided')
    print('-q or --quantile-color-scale <value> value is an integer between 0 and 100 specifying the upper threshold quantile for the color limit')
    print('-p or --print-each <outputdir> <summaryname> <format as png,eps, or jpg> print each image individually')'

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
    summaryname='summary'
    image_name=''
    image_ext=''
    print_each=False
    latex_print=False
    colormap_name='pink'
    issparse=True
    p=100
    while ii < len(args):
        if (args[ii]== '-h')|(args[ii]== '--help'):
           printhelp()
           sys.exit()
        if (args[ii]== '-p')|(args[ii]== '--print-each'):
            ii+=1
            report_directory = str(args[ii])
            image_directory = os.path.join(report_directory,'figures')
            ii+=1
            summary_name = str(args[ii])
            ii+=1
            image_ext = str(args[ii])
            print_each=True
            latex_print=True
            # create the report directory if it does not exist
            try:
                if not os.path.exists(report_directory):
                    os.mkdir(report_directory)
                    print("Directory ", report_directory, " Created ")
                if not os.path.exists(image_directory):
                    os.mkdir(image_directory)
                    print("Directory ", image_directory, " Created ")
            except:
                print("Potentially failed to create report directories")
        if (args[ii]== '-c')|(args[ii]== '--colormap'):
            ii+=1
            colormap_name = str(args[ii])
        if (args[ii]== '-d')|(args[ii]== '--dense'):
            issparse=False
        if (args[ii]== '-q')|(args[ii]== '--quantile-color-scale'):
            ii+=1
            p = float(args[ii])
        if (args[ii]== '-i')|(args[ii]== '--input-files'):
           inputfiles=[]
           ii+=1
           while ii < len(args):
               inputfiles.append(args[ii])
               ii+=1
        ii+=1


    latex=lr.make_latex_report_header('figures')
    for inputfile in inputfiles:
        ext=os.path.splitext(inputfile)[-1]
        filename=os.path.splitext(os.path.basename(inputfile))[0]+'_'+ext.split('.')[1]
        print(filename)
        if issparse:
            data=mp.load_sparse_matrix(inputfile)
            q=np.percentile(data.data,p)
            matrix=data.todense()
        else:
             #print "\nLoading "+ datatype+'\n'



             if (ext=='.npy')| (ext=='.npz'):
                 matrix=np.load(inputfile)
             elif (ext=='.csv'):
                 matrix=np.loadtxt(inputfile,delimiter=',',dtype=np.float)
                 #keyboard()
             elif (ext=='.mat'):
                 data=loadmat(inputfile)
                 for key in list(data):
                     if '__' not in key:
                         matrix=data[key]
             elif (ext=='.json'):
                 with open(inputfile, 'r') as jsonfile:
                     data=json.load(jsonfile)
                 for key in list(data):
                     if '__' not in key:
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
         if image_ext not in ["png",'pdf','svg','ps','eps']:
             print('invalid image format:'+image_ext)
         else:
            imagename=os.path.join(image_directory,filename+'.'+image_ext)
            fig1.savefig(imagename)
            if latex_print:
                (m,n)=np.shape(matrix)
                params=dict()
                params['inputfile']=inputfile
                params['table width']=1
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
                latex_table=lr.make_latex_table(params)
                latex+=latex_table
                '''
                latex_table=u'\n'E
                latex_table+=r'\begin{lstlisting}' +u'\n'
                latex_table+=inputfile +u'\n'
                latex_table+=r'\end{lstlisting}' +u'\n'
                latex_table+=r'\includegraphics[width=.5\textwidth]{'+imagename+r'}\\ '+u'\n'
                latex_table+=r'\begin{tabular}{cc}'+u'\n'
                latex_table+="Row dimension& "+str(m)+r"\\"+u'\n'
                latex_table+="Column dimension& "+str(n)+r"\\"+u'\n'
                latex_table+="Sum& "+str(np.sum(matrix))+r"\\"+u'\n'
                latex_table+="Max& "+str(np.max(matrix))+r"\\"+u'\n'
                latex_table+="Non-zero entries& "+str(np.sum(matrix>0))+r"\\"+u'\n'
                latex_table+="Total entries& "+str(n*m)+r"\\"+u'\n'
                latex_table+="Non-trivial rows & "+str(np.sum(sum(matrix>0)>0))+r"\\"+u'\n'
                latex_table+=r'\end{tabular}'+u'\n'

                latex_table+=u'\n'
                '''
                print(latex_table)
                latex+=latex_table
    if latex_print:
        latex+=r'\end{document}'+u'\n'
        print(latex)
        with open(os.path.join(report_directory, summary_name+ '.tex'), 'w') as result:
            result.write(latex)

if __name__ == "__main__":
   main()
