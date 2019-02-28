"""
An assortment of low level functions that I use. Some of them resemble 
functions similiar to a popular matrix laboratory.
"""
import code
import sys
import os
import numpy as np
from numpy import array,shape,matrix,ones,transpose,sqrt,mod,arctan2,pi
import matplotlib.pyplot as plt
from itertools import count

def significant_figures(x,sig=4):
    """
    create string without decimal points for file names
    """
    if x==0:
        return "0"
    else:
        dec=int(np.floor(np.log10(abs(x))))-sig
        xround=np.round(x,-dec)
        return str(int(xround*pow(10,-dec)))+'e'+str(dec)
         
def float_range(start, step, stop):
    ''' Function that mimics the start:step:stop '''
    # like a:b:c 
    arg=[]
    for i in count():
        yld = start + step*i        
        if yld > stop:
            break
        else:
            arg.append(yld)
    return array(arg)

def keyboard(banner=None):
    ''' Function that mimics the keyboard command '''
    # use exception trick to pick up the current frame
    try:
        raise None
    except:
        frame = sys.exc_info()[2].tb_frame.f_back
    print( "# Use quit() to exit :) Happy debugging!")
    # evaluate commands in current namespace
    namespace = frame.f_globals.copy()
    namespace.update(frame.f_locals)
    try:
        code.interact(banner=banner, local=namespace)
    except SystemExit:
        return 
def distancematrix(p0,p1=None,getangles=False):
    """given two sets of points p0 (mX2 array) and p1 (nX2 array), compute the pairwise distance matrix D (mXn arra)"""
    if p1 is None:
        p1=p0
    m=shape(p0)[0]
    n=shape(p1)[0]
    p0=matrix(p0)
    p1=matrix(p1)
    xdiff=p0[:,0]*ones((1,n))-ones((m,1))*transpose(p1[:,0])
    ydiff=p0[:,1]*ones((1,n))-ones((m,1))*transpose(p1[:,1])
    xsqrdist=pow(array(xdiff),2)
    ysqrdist=pow(array(ydiff),2)
    distance=sqrt(xsqrdist+ysqrdist)
    if getangles:
        return matrix(distance),matrix(mod(arctan2(xdiff,ydiff),2*pi))
    else:
        return matrix(distance)
def differencematrix(p0,p1=None):
    """given two sets of points p0 (mX2 array) and p1 (nX2 array), compute the pairwise differance matrixes xdiff,ydiff (mXn arrays)"""
    if p1 is None:
        p1=p0
    m=shape(p0)[0]
    n=shape(p1)[0]
    p0=matrix(p0)
    p1=matrix(p1)
    xdiff=p0[:,0]*ones((1,n))-ones((m,1))*transpose(p1[:,0])
    ydiff=p0[:,1]*ones((1,n))-ones((m,1))*transpose(p1[:,1])
    return xdiff,ydiff
cls=lambda: os.system('cls')
whos=lambda: [v for v in globals().keys() if not v.startswith('_')]

class getlines:    
    """ 
    like getpts, but with line in between
    """
    def __init__(self,fig):
        self.line = plt.plot([0],[0])[0]
        self.x = []
        self.y = []
        self.fig=fig
        self.cid1 = self.line.figure.canvas.mpl_connect('button_press_event', self.on_click)        
        self.cid2 = self.line.figure.canvas.mpl_connect('key_press_event', self.on_key)
        self.go=True
        print('Press Enter when you have finished getting points.')
    def on_click(self,event):
        self.x.append(event.xdata)
        self.y.append(event.ydata)
        if len(self.x)>0:
            plt.scatter(self.x,self.y ,s=100,c='black',edgecolor='white',linewidth=2)
        if len(self.x)>1:
            plt.plot(self.x,self.y,c='black',linewidth=2)
            plt.draw()
    def on_key(self,event):
        if event.key=='enter':
            print('Exited getpts')
            self.fig.canvas.mpl_disconnect(self.cid1)
            self.fig.canvas.mpl_disconnect(self.cid2)
            self.go=False
    def __call__(self, event):
        if event.inaxes!=self.line.axes: return
class getpts:    
    """ 
    like getpts
    """
    def __init__(self,fig):
        self.line = plt.plot([0],[0])[0]
        self.x = []
        self.y = []
        self.fig=fig
        self.cid1 = self.line.figure.canvas.mpl_connect('button_press_event', self.on_click)        
        self.cid2 = self.line.figure.canvas.mpl_connect('key_press_event', self.on_key)
        self.go=True
        print('Press Enter when you have finished getting points.')
    def on_click(self,event):
        self.x.append(event.xdata)
        self.y.append(event.ydata)
        if len(self.x)>0:
            plt.scatter(self.x,self.y ,s=100,c='black',edgecolor='white',linewidth=2)
    def on_key(self,event):
        if event.key=='enter':
            print('Exited getpts')
            self.fig.canvas.mpl_disconnect(self.cid1)
            self.fig.canvas.mpl_disconnect(self.cid2)
            self.go=False
    def __call__(self, event):
        if event.inaxes!=self.line.axes: return
class getptstk:    
    """
    like getpts
    """
    def __init__(self,fig,tk):
        self.line = plt.plot([0],[0])[0]
        self.x = []
        self.y = []
        self.tk=tk
        self.fig=fig
        self.cid1 = self.line.figure.canvas.mpl_connect('button_press_event', self.on_click)        
        self.cid2 = self.line.figure.canvas.mpl_connect('key_press_event', self.on_key)
        self.go=True
        print('Press Enter when you have finished getting points.')
    def on_click(self,event):
        self.x.append(event.xdata)
        self.tk.longitude=self.x
        self.y.append(event.ydata)
        self.tk.latitude=self.y
        if len(self.x)>0:
            plt.scatter(self.x,self.y ,s=100,c='black',edgecolor='white',linewidth=2)
    def on_key(self,event):
        if event.key=='enter':
            print('Exited getpts')
            self.fig.canvas.mpl_disconnect(self.cid1)
            self.fig.canvas.mpl_disconnect(self.cid2)
            self.go=False
    def __call__(self, event):
        if event.inaxes!=self.line.axes: 
            return
def file_len(fname):
    """
    This counts the length of a file. Just got threwn in here because it is low level
    """
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
