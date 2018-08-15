# -*- coding: utf-8 -*-
"""

I created some wrapers using matplotlib to make it feel more like a matrix 
laboratory.

Created on Fri Apr 21 12:36:46 2017

@author: Michael Rosenthal
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import SymLogNorm,Normalize

def set_axis_equal(ax,xlimit=None,ylimit=None,zlimit=None):
    '''
    set axis equal
    '''
    w=np.array([np.diff(ax.get_xlim()),
    np.diff(ax.get_ylim()),
    np.diff(ax.get_zlim())])
    s=max(w)/w
    if xlimit is None:
        xlimit=ax.get_xlim()
    if ylimit is None:
        ylimit=ax.get_ylim()      
    if zlimit is None:
        zlimit=ax.get_zlim()            
    limits=np.array([xlimit,
                ylimit,
                zlimit])
    m=np.array([np.mean(ax.get_xlim()),
    np.mean(ax.get_ylim()),
    np.mean(ax.get_zlim())])
    newlimits=np.array([s[:,0]*(limits[:,0]-m)+m,
                    s[:,0]*(limits[:,1]-m)+m])
    ax.set_xlim(newlimits[:,0])
    ax.set_ylim(newlimits[:,1])
    ax.set_zlim(newlimits[:,2])
    return ax

def plot_curve(curve,t=None,ax=None,fig=None):
    '''
    plot a 2D or 3d parametrized curve
    '''

    [d,n] = np.shape(curve)
    if d>n:
        curve=curve.T
        [d,n] = np.shape(curve)
    if fig is None:
        fig=plt.figure()        
    if ax is None:
        if d==2:# if the curve is 2 dimensional
            ax=fig.add_subplot(111)
        if d==3:# if the curve is 3 dimensional
            ax=fig.add_subplot(111,projection='3d')
    line=None;
    if t is None: # if temperal parametrization is not given
        if d==2:# if the curve is 2 dimensional
            line,=ax.plot(curve[0,:],curve[1,:])
        if d==3:# if the curve is 3 dimensional
            line,=ax.plot(curve[0,:],curve[1,:],curve[2,:])            
    else:  # if temperal parametrization is given
        cn=Normalize(min(t),max(t))
        if d==2:# if the curve is 2 dimensional
            for i in xrange(n-1):
                line,=ax.plot(curve[0,i:i+2],curve[1,i:i+2],color =plt.cm.jet(cn(t[i])));
        if d==3:# if the curve is 3 dimensional
            for i in xrange(n-1):
                line,=ax.plot(curve[0,i:i+2],curve[1,i:i+2],curve[2,i:i+2],color =plt.cm.jet(cn(t[i])));                
            set_axis_equal(ax)
    
    return fig,ax
def plot_pts(curve,t=None,ax=None,fig=None):
    '''
    plot a 2D or 3d parametrized curve
    '''

    [d,n] = np.shape(curve)
    if d>n:
        curve=curve.T
        [d,n] = np.shape(curve)
    if fig is None:
        fig=plt.figure()        
    if ax is None:
        if d==2:# if the curve is 2 dimensional
            ax=fig.add_subplot(111)
        if d==3:# if the curve is 3 dimensional
            ax=fig.add_subplot(111,projection='3d')
    line=None;
    if t is None: # if temperal parametrization is not given
        if d==2:# if the curve is 2 dimensional
            line=ax.scatter(curve[0,:],curve[1,:],marker='.')
        if d==3:# if the curve is 3 dimensional
            line=ax.scatter(curve[0,:],curve[1,:],curve[2,:],marker='.')            
    else:  # if temperal parametrization is given
        cn=Normalize(min(t),max(t))
        if d==2:# if the curve is 2 dimensional
            for i in xrange(n-1):
                line=ax.scatter(curve[0,i:i+2],curve[1,i:i+2],markeredgecolor='w',color =plt.cm.jet(cn(t[i])),marker='.');
        if d==3:# if the curve is 3 dimensional
            line=ax.scatter(curve[0,:],curve[1,:],curve[2,:],c=t);                
            set_axis_equal(ax)
    
    return fig,ax    
def plot_curves(curves,t=None,ax=None,fig=None):
    '''
    plot a list of 2D or 3d parametrized curves
    '''
    [d,n] = np.shape(curves[0])
    if fig is None:
        fig=plt.figure()        
    if ax is None:
        if d==2:# if the curve is 2 dimensional
            ax=fig.add_subplot(111)
        if d==3:# if the curve is 3 dimensional
            ax=fig.add_subplot(111,projection='3d')
    n=len(curves)
    if t is None:
        t=np.linspace(0,1,len(curves[ii]))
    for ii in range(n):    
        plot_curve(curves[ii],t,ax,fig);
    set_axis_equal(ax)
    return fig,ax    
    
def matshow(matrix,fig=None,cmap="pink_r"):
    '''
    plot a matrix as an image
    '''
    if fig is None:
        fig=plt.figure()
    [n1,n2] = np.shape(matrix)
    fig.clf()
    ax = fig.add_subplot(111)
    ax.set_xlim((0), n1)
    ax.set_ylim(n2,0)
    m = ax.matshow(matrix,
               cmap=cmap, norm=SymLogNorm(1), origin="bottom")
    
    return fig,ax,m
    
