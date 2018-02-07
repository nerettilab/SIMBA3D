# -*- coding: utf-8 -*-
"""
Square Root Velocity Function tools for open curve in R^n

Incomplete module 

@author: Michael Rosenthal
"""

import numpy as np
import numpy.linalg as npla
from scipy.integrate import cumtrapz



def extrap(x,xp,yp):
    """np.interp with linear extrapolation"""
    y= np.interp(x,xp,yp)
    y=np.where(x<xp[0],yp[0]+(x-xp[0])*(yp[0]-yp[1])/(xp[0]-xp[1]),y)
    y=np.where(x>xp[-1],yp[-1]+(x-xp[-1])*(yp[-1]-yp[-2])/(xp[-1]-xp[-2]),y)
    return y    

def interp(t,X,tnew=None):
    """
    interpolation with linear extrapoltion
    """
    n=len(t)
    [d,nX]=np.shape(X)
    transpose=False
    if nX != n:
        if d==n:
            X=X.T
            [d,nX]=np.shape(X)
            transpose=True
        else:
           return -1;
    if tnew  is None:
        tnew= t
    n2=len(tnew)
    Xnew=np.zeros([d,n2])
    for ii in range(d):
        Xnew[ii,:]=extrap(tnew,t,X[ii,:])
    if transpose:
        return Xnew.T
    else:
        return Xnew        
        
def gradient(curve):
    '''
    row-wise gradient for unifromly sampled curve
    '''
    [n,T] = np.shape(curve)
    v=np.NaN*np.ones([n,T]);
    for ii in range(n):
        v[ii,:]=np.gradient(curve[ii,:],1.0/(T-1))
    return v;


def length(curve,t=None):
    """
    Return the length of a parametrized curve
    """
    [n,T]=np.shape(curve)
    if n>T:
        curve=curve.T
        T=n
    if t is None:
        t=np.linspace(0,1,T)        
    curvedot=gradient(curve)
    normcurvedot=np.sqrt(np.sum(pow(curvedot,2.0),0))
    return np.trapz(normcurvedot,np.linspace(0,1,T))
    
def scale_curve(curve,t=None):
    '''
    Scale a curve to unit length
    '''
    [M,T] =np.shape(curve)
    if M>T:
        curve=curve.T
        T=M
    if t is None:
        t=np.linspace(0,1,T)         
    scale=length(curve,t)
    return curve/scale,scale

    
def curve_to_psi(curve):
    """
    Curve to square root speed function
    """
    curvedot=gradient(curve)
    normcurvedot=np.sqrt(np.sum(pow(curvedot,2.0),0))
    return np.sqrt(normcurvedot)
    
def roughness2(curve,t=None):
    """
    Return the first (same as length) and second order roughness measure of a 
    parametrized curve
    
    The first and second order gradient are computed numerically, and the 
    roughness measure is simply the Frobenious norm of first and second order
    gradients.
    
    If the temperal parametrization t is not specified then the curve is 
    assumed to be parametrized uniformly over the interval [0,1]
    """    
    [n,T]=np.shape(curve)
    if n>T:
        curve=curve.T
        T=n
    if t is None:
        t=np.linspace(0,1,T)
    curvedot=gradient(curve)
    normcurvedot=np.sqrt(np.sum(pow(curvedot,2.0),0))
    R1=np.trapz(normcurvedot,t)
    
    curvedotdot=gradient(curvedot)
    normcurvedotdot=np.sqrt(np.sum(pow(curvedotdot,2.0),0))
    R2=np.trapz(normcurvedotdot,t)
    
    return R1,R2
    
    
def InnerProd_Q(q1,q2):
    '''
    Inner product in the Q space
    '''
    [n,T]=np.shape(q1)
    return np.trapz(np.sum(q1*q2,0),np.linspace(0,1,T))
    
def curve_to_q(curve):
    ''' 
    Convert function from square root velocity space back to the original function space
    '''
    [n,T]=np.shape(curve)
    transpose=False;
    if T<n:
        curve=curve.T
        [n,T] = np.shape(curve)  
        transpose=True;
    v = gradient(curve)
    length = np.sum(np.sqrt(np.sum(v*v)))/n
    v=v/length
    q=np.zeros([n,T]);
    eps=0.0001;
    5
    
    for i in range(T):
        L = np.sqrt(np.linalg.linalg.norm(v[:,i]));
        if L > eps:
            q[:,i] = v[:,i]/L;
        else:
            q[:,i] = v[:,i]*eps;
    if transpose:
        return q.T
    else:
        return q

def q_to_curve(q):
    ''' 
    Convert function from square root velocity space back to the original function space
    '''
    [n,T]=np.shape(q);
    transpose=False;
    if T<n:
        q=q.T
        [n,T] = np.shape(q)  
        transpose=True;    
    qnorm =np.ones(T)
    for ii in range(T):
        qnorm[ii] =np.linalg.linalg.norm(q[:,ii])
    integrand=np.zeros([n,T])
    for ii in range(n):
        integrand[ii,:]= q[ii,:]*qnorm
    p=cumtrapz(integrand,axis=1,initial=0)/T
    if transpose:
        return p.T
    else:
        return p
        
def resample_curve(curve,N):
    '''
    resample a parametrized curve uniformly at N time points
    '''
    [n,T]=np.shape(curve)
    transpose=False;
    if T<n:
        curve=curve.T
        [n,T] = np.shape(curve)  
        transpose=True;  
    delt= np.zeros(T)
    for r in range(1,T):
        delt[r]=np.linalg.norm(curve[:,r]-curve[:,r-1])
    cumdelt=np.cumsum(delt)/np.sum(delt)
    newdelt=np.linspace(0,1,N)
    Xn=np.zeros([n,N]);
    Xn= interp(cumdelt,curve,newdelt)
    '''
    # for closed curves
    for i in range(n):
        Xn[i,:]= spline(cumdelt,curve[i,:],newdelt)
    '''
    if transpose:
        return Xn.T
    else:
        return Xn   

def center_curve(curve):
    '''
    center the curve to the origin
    '''
    [n,T]=np.shape(curve)
    transpose=False;
    if T<n:
        curve=curve.T
        [n,T] = np.shape(curve)  
        transpose=True;  
    centered=curve;
    mu=np.zeros([n,1])
    for d in range(n):
        mu[d]=np.mean(curve[d,:])
        centered[d,:]-=mu[d]
    if transpose:
        return (centered.T,mu.T)
    else:
        return (centered,mu)


def standarize_curve(curve):
    '''
    center the curve to the origin
    '''
    (centered,mu)=center_curve(curve);
    (standardized,scale)=scale_curve(centered);
    return (standardized,scale,mu)
    
def find_best_rotation(curve1,curve2):
    '''
    Optimally rotate the second curve to the first curve
    '''
    [n,T]=np.shape(curve1)

    transpose=False;
    if T<n:
        curve1=curve1.T
        curve2=curve2.T
        [n,T] = np.shape(curve1)  
        transpose=True
    [T,n1] = np.shape(curve1)  
    [T,n2] = np.shape(curve2) 
    if n1!=n2: # if dimensions do not align then resampleW
        n= max([n1,n2])
        t1=np.linspace(0,1,n1)
        t2=np.linspace(0,1,n2)
        t=np.linspace(0,1,n)
        recurve1=interp(t1,curve1,t)   
        recurve2=interp(t2,curve2,t)
        A=np.matrix(recurve1)*np.matrix(recurve2.T)
    else:        
        A=np.matrix(curve1)*np.matrix(curve2.T)
    U,S,V=npla.svd(A)
    S=np.eye(n)
    if (npla.det(A)<0):
        S[:,-1]=-S[:,-1]
    R=U*S*V
    Rcurve2=np.array(R*np.matrix(curve2))        
    if transpose:
        return (Rcurve2.T,R)
    else:
        return (Rcurve2,R)
        
def find_best_ortho(curve1,curve2):
    '''
    Optimally rotate or reflect the second curve to the first curve
    '''

    [n,T]=np.shape(curve1)
    transpose=False;
    if T<n:
        curve1=curve1.T
        curve2=curve2.T
        transpose=True
    [T,n1] = np.shape(curve1)  
    [T,n2] = np.shape(curve2) 
    if n1!=n2: # if dimensions do not align then resampleW
        n= max([n1,n2])
        t1=np.linspace(0,1,n1)
        t2=np.linspace(0,1,n2)
        t=np.linspace(0,1,n)
        recurve1=interp(t1,curve1,t)   
        recurve2=interp(t2,curve2,t)
        A=np.matrix(recurve1)*np.matrix(recurve2.T)
    else:
        A=np.matrix(curve1)*np.matrix(curve2.T)
    U,S,V=npla.svd(A)
    R=U*V
    Rcurve2=np.array(R*np.matrix(curve2))        
    if transpose:
        return (Rcurve2.T,R)
    else:
        return (Rcurve2,R) 
        
def standardize_curves(curves):
    '''
    First removes translation and scale, and then rotates the remaining curves to first curve.
    '''
    n=len(curves)
    mu_curve,mu_scale,mu_center=standarize_curve(curves[0]);
    standarized_curves=[mu_curve]
    for ii in range(1,n):
        ii_curve,ii_scale,ii_center=standarize_curve(curves[ii]);
        #Rest_curve,R=find_best_rotation(mu_curve,ii_curve)
        Rest_curve,R=find_best_ortho(mu_curve,ii_curve)
        standarized_curves.append(Rest_curve)
    return standarized_curves;