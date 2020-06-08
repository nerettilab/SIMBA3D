#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 09:36:32 2019

@author: star
"""
import numpy as np
from simba3d.matrixlabish import significant_figures,keyboard


def make_latex_report_header(image_path):
    latex=""
    latex+=r'\documentclass{report}'+u'\n'
    latex+=r'\usepackage{verbatim}'+u'\n'
    latex+=r'\usepackage{listings}'+u'\n'
    latex+=r'\lstset{breaklines=true,basicstyle=\ttfamily}'+u'\n'
    latex+=r'\usepackage[margin=2cm]{geometry}'+u'\n'
    latex+=r'\usepackage{graphicx}'+u'\n'
    latex+=r'\graphicspath{{'+image_path+r'/}}'+u'\n'
    latex+=r'\begin{document}'+u'\n'
    latex+=r'\small'+u'\n'
    return latex
def make_latex_table(params):
    latex_table=u'\n'
    latex_table+=r'\begin{lstlisting}' +u'\n'
    if 'inputfile' in params:
        latex_table+=params['inputfile'] +u'\n'
    latex_table+=r'\end{lstlisting}' +u'\n'
    if 'table width' in params:
        ncols=params['table width']
    else:
        ncols=1.0;

    if 'images' in params:
        latex_table+=r'\begin{tabular}{'
        for ii in range(ncols):
            latex_table+='c'
        latex_table+='}'+u'\n'
        nrows=np.ceil(float(len(params['images']))/float( ncols))
        width_percentage=str(round(1.0/float(ncols),2)-0.01)
        itr=0;
        for ii in range(int(nrows)):
            for jj in range(int(ncols)):
                if itr< len(params['images']):
                    latex_table+=r'\includegraphics[width='
                    latex_table+=width_percentage
                    latex_table+=r'\textwidth]{'
                    latex_table+=params['images'][itr]
                    latex_table+=r'}'
                    latex_table+=u'\n'
                if jj < (ncols-1):
                    latex_table+=r'&'
                else:
                    latex_table+=r'\\'
                latex_table+=u'\n'
                itr+=1;
        latex_table+=r'\end{tabular}'+u'\n'

    if 'statistics' in params:
        latex_table+=r'\begin{tabular}{cc}'+u'\n'
        for key in list(params['statistics']):
            latex_table+=key
            latex_table+=r'&'
            latex_table+=str(params['statistics'][key])
            latex_table+=r'\\'
            latex_table+=u'\n'
        latex_table+=r'\end{tabular}'+u'\n'
        latex_table+=u'\n'
    #print(latex_table)
    return latex_table
def make_latex_section(section_name,exact_text=None):
    latex=r'\section{'
    latex+=section_name
    latex+=r'}'+u'\n'
    if exact_text is not None:
        latex+=r'\begin{lstlisting}' +u'\n'
        latex+=exact_text+u'\n'
        latex+=r'\end{lstlisting}' +u'\n'
    return latex
def make_latex_subsection(section_name,exact_text=None):
    latex=r'\subsection{'
    latex+=section_name
    latex+=r'}'+u'\n'
    if exact_text is not None:
        latex+=r'\begin{lstlisting}' +u'\n'
        latex+=exact_text+u'\n'
        latex+=r'\end{lstlisting}' +u'\n'
    return latex
def make_latex_report_footer():
    return r'\end{document}'+u'\n'
