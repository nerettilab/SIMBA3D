# -*- coding: utf-8 -*-
"""
Optimization tools are are defined here.

A ptyhon class called opt_E manages the scipy optimize options and creates a
friendlier user interface.

Created on Thu Apr 20 10:26:05 2017

@author: Michael Rosenthal
"""

import time

import numpy as np
import scipy.spatial.distance as ssd
from scipy.spatial.distance import pdist
from scipy.sparse import csr_matrix
from scipy.optimize import minimize,check_grad
from scipy.misc import comb

from simba3d.matrixlabish import keyboard

#import plotting_tools as pt
import simba3d.srvf_open_curve_Rn as srvf # most of this is not used right now

# set the data precision used within the algorithm
# !incomplete task! (only done for h1 and data term)
# This was started as a debug process, but I do not see a reason to finish it
# or undo it.
DTYPE=np.float # there is no point in changing this now,

class stop_optimizing_exception(Exception):
    pass

def listIndices(n,offdiag):
    '''
    List the pairs of index values (i,j) that indicate the upper triangular
    entries of an n x n matrix.

    Inputs:
        n - The size of the square matrix.
        offdiag - An integer value indicating the first off-diagonal to contain
        usable data in the matrix C (usually offdiag = 1).

    Outputs:
        idx - An nchoosek(n,2) x 2 matrix consisting of the (i,j) index
        combinations such that i<j.
'''
    kk=int(0)
    numberofpairs=int(comb(n-offdiag+1,2))
    idx=np.zeros((numberofpairs,2),dtype='int')
    for ii in range(n-offdiag):
        for jj in range(ii+offdiag,n):
            idx[kk,0]=ii
            idx[kk,1]=jj
            kk=kk+1
    return idx


def scale_initialized_curve(curve,pairwise_contact_matrix,a,b,index,idxmiss):
    '''
    optimally scale the curve to the data for fewer iterations
    '''
    X,mu=srvf.center_curve(curve)
    C=pairwise_contact_matrix
    if idxmiss.size:
        X=np.delete(X,idxmiss,1)
        C=np.delete(C,idxmiss,0)
        C=np.delete(C,idxmiss,1)
    n=len(C)
    index=np.tril(np.ones((n,n)),-1)==1
    a=np.double(a);
    b=np.double(b);
    D=np.array(compute_pairwise_distance(X),dtype=np.double)
    c=np.array(C[index],dtype=np.double)
    d=D[index]
    scale= pow(sum(c)/(b*sum(pow(d,a))),1/a)

    return scale*curve

def compute_pairwise_distance(curve):
    """
    Compute the pairwise distance matrix of a parametrized curve

    This function assumes that the longer dimension is the dimension of the
    parameterization. The reason for this is to avoid errors caused by
    transposing the matrix.
    """

    [n,T]=np.shape(curve)
    if n<T: # make sure the larger dimension is used
        curve=curve.T
    return np.array(ssd.squareform(ssd.pdist(curve,'euclidean')))

def loglikelihood_Varoquaux(pairwise_distance_matrix,a,b,pairwise_contact_matrix):
    """
    the loglikelihood function using the poisson model proposed in the
    Varoquaux paper.
    """
    L=np.sum(-b*pow(pairwise_distance_matrix,a)+pairwise_contact_matrix*(np.log(b)*a*np.log(pairwise_distance_matrix)))
    return L

def h1penalty(curve,S1):
    '''
    Compute the h1 penalty and its gradient.

    The penalty can be interpreted as the sample variance of the normalized
    curve
    '''
    #increase precision here

    # makes sure the curve is oriented correctly
    n,p=np.shape(curve)
    if p>n:
        curve=curve.T
        n,p=np.shape(curve)


    # initialized computation parameters

    term=np.zeros([n,p],dtype=DTYPE);
    d= np.array(curve[range(1,n),:]-curve[range(n-1),:],dtype=DTYPE)
    w=np.sqrt(np.sum(pow(d,2.0),1,dtype=DTYPE))
    term[0,:]=-d[0,:]/w[0]
    for ii in range(1,n-1):
        term[ii,:]=d[ii-1,:]/w[ii-1]-d[ii,:]/w[ii]
    term[n-1,:]=d[n-2,:]/w[n-2]
    s=np.sum(pow(w,2.0),dtype=DTYPE)
    L=np.sum(w,dtype=DTYPE)

    # penalty value
    num= (n-1.0)*s
    den=pow(L,2.0)
    h1= num/den-1.0
    # gradient
    mat=np.array(S1*np.matrix(curve),dtype=DTYPE)
    gradh1=(n-1.0)*(mat/pow(L,2.0)-2.0*term*s/pow(L,3.0))
    #print h1'
    return h1,gradh1.T

def h2penalty(curve):
    '''
    Compute the h2 penalty and its gradient.

    the average \cos(\theta_i), where \theta_i is the angle between the
    triplet (x_{i-1},x_i,x_{i+1}).
    '''

    # makes sure the curve is oriented correctly
    n,p=np.shape(curve)
    transpose=False
    if p>n:
        curve=curve.T
        n,p=np.shape(curve)
        transpose=True

    # initialized computation parameters
    h2=0
    gradh2=np.zeros([n,p]);

   # compute differences between ajacent nodes
    d= curve[range(1,n),:]-curve[range(n-1),:]

    # compute weights according to gaussian
    w=np.sqrt(np.sum(pow(d,2.0),1))

    # loop through the terms of the gradient
    for ii in range(n-1):
        if ii>0:

            # update the second order penalty
            h2 -= sum(d[ii-1,:]*d[ii,:])/(w[ii-1]*w[ii])

            # comput the gradient
            if ii==1:
                # compute the first gradient term
                gradh2[0,:]=d[0,:]*sum(-d[0,:]*d[1,:])/(pow(w[0],3.0)*w[1])
                gradh2[0,:]+=d[1,:]/(w[0]*w[1])
            elif ii==2:
                    # compute the second gradient term
                    gradh2[1,:]=d[0,:]*(1/(w[0]*w[1])+np.sum(d[0,:]*d[1,:])/(pow(w[0],3.0)*w[1]))
                    gradh2[1,:]+=d[1,:]*(-1/(w[0]*w[1])-np.sum(d[0,:]*d[1,:])/(w[0]*pow(w[1],3))-np.sum(d[1,:]*d[2,:])/(pow(w[1],3.0)*w[2]))
                    gradh2[1,:]+=d[2,:]*(1/(w[1]*w[2]))
            else:
                    # compute the i^th gradient term (excluding the first and last two terms)
                    gradh2[ii-1,:]=d[ii-3,:]*(-1/(w[ii-3]*w[ii-2]))
                    gradh2[ii-1,:]+=d[ii-2,:]*(np.sum(d[ii-3,:]*d[ii-2,:])/(w[ii-3]*pow(w[ii-2],3.0))+1/(w[ii-2]*w[ii-1])+np.sum(d[ii-2,:]*d[ii-1,:])/(pow(w[ii-2],3.0)*w[ii-1]))
                    gradh2[ii-1,:]+=d[ii-1,:]*(-1/(w[ii-2]*w[ii-1])-np.sum(d[ii-2,:]*d[ii-1,:])/(w[ii-2]*pow(w[ii-1],3.0))-np.sum(d[ii-1,:]*d[ii,:])/(pow(w[ii-1],3.0)*w[ii]))
                    gradh2[ii-1,:]+=d[ii,:]*(1/(w[ii-1]*w[ii]));

    # compute the next to last gradient term
    gradh2[n-2,:]=d[n-4,:]*(-1/(w[n-4]*w[n-3]))
    gradh2[n-2,:]+=d[n-3,:]*(np.sum(d[n-4,:]*d[n-3,:])/(w[n-4]*pow(w[n-3],3.0))+1/(w[n-3]*w[n-2])+np.sum(d[n-3,:]*d[n-2,:])/(pow(w[n-3],3.0)*w[n-2]))
    gradh2[n-2,:]+=d[n-2,:]*(-1/(w[n-3]*w[n-2])-np.sum(d[n-3,:]*d[n-2,:])/(w[n-3]*pow(w[n-2],3.0)));

    # compute the last gradient term
    gradh2[n-1,:]=d[n-3,:]*(-1/(w[n-3]*w[n-2]))
    gradh2[n-1,:]+=d[n-2,:]*(np.sum(d[n-3,:]*d[n-2,:])/(w[n-3]*pow(w[n-2],3.0)));

    # final scaling adjustment to terms
    h2=h2/(n-2)
    gradh2=gradh2/(n-2)

    if transpose: # if the curve was transposed, to fit the function, transpose it again the match the original input
        gradh2=gradh2.T
    #print h2
    return h2,gradh2

def h4penalty(curve,radius=1.0):
    '''
    Compute the h4 lamina penalty and its gradient.

    '''
    # makes sure the curve is oriented correctly
    n,p=np.shape(curve)
    transpose=False
    if p>n:
        curve=curve.T
        n,p=np.shape(curve)
        transpose=True
    normed=np.sqrt(np.sum(pow(curve,2.0),1))
    beyondind=normed>radius
    h4=0.0
    gradh4=np.zeros([n,p])
    return h4,gradh4

def parametrization_error(curve,true_square_root_speed=None,t=None):
    """
    Error of the curve parameterization given a square roor speed.

    If the true square root speed is not specified, then uniform speed is used
    by default.
    """
    M,T=np.shape(curve)
    if M>T:
        curve=curve.T
        T=M

    if t is None:
        t=np.linspace(0,1,T)
    if true_square_root_speed is None:
        true_square_root_speed=1.0/T
    curve,scale=srvf.scale_curve(curve,t)
    psi_curve=srvf.curve_to_psi(curve)
    return np.trapz(pow(psi_curve-true_square_root_speed,2.0),t)

def penalized_log_likelihood(curve,t,pairwise_contact_matrix,a,b,term_weights,square_root_speed=None,pairwise_distance_matrix=None):
    """
    penalized log likelihood
    """
    if pairwise_distance_matrix is None: #if the do not already have the pairwise distance matrix computed, then compute it
        pairwise_distance_matrix=compute_pairwise_distance(curve)

    L=0     # initialize log likelihood term
    R1=0    # initialize first order term
    R2=0    # initialize second order term
    Q=0     # initialize parametrization penalty term
    S=0     # initialize shape prior term
    if term_weights[0]!=0:
        L=term_weights[0]*loglikelihood_Varoquaux(pairwise_distance_matrix,a,b,pairwise_contact_matrix)
    if (term_weights[1]!=0)&(term_weights[2]==0):
        R1=term_weights[1]*srvf.length(curve,t)
    if (term_weights[2]!=0):
        R1,R2=srvf.roughness2(curve,t)
        R1=term_weights[1]*R1
        R2=term_weights[2]*R2
    if (term_weights[3]!=0):
        Q=term_weights[3]*parametrization_error(curve,square_root_speed,t)
    if (term_weights[4]!=0):
        S=term_weights[4]*0 # not implemented yet
    return L-R1-R2-Q-S



def dif(x1,x2):
    return x1-x2

def analytical_data_gradient0(curve,pairwise_contact_matrix,a,b,index):
    '''
    compute log-likelihood gradient analytically
    '''

    m,n=np.shape(curve);
    x=np.zeros([n,1],dtype=np.float64)
    x[:,0]=curve[0,:]
    y=np.zeros([n,1],dtype=np.float64)
    y[:,0]=curve[1,:]
    z=np.zeros([n,1],dtype=np.float64)
    z[:,0]=curve[2,:]
    diffx=np.array(ssd.squareform(pdist(x,dif)),dtype=np.float64)
    diffx[index]=-diffx[index]
    diffy=np.array(ssd.squareform(pdist(y,dif)),dtype=np.float64)
    diffy[index]=-diffy[index]
    diffz=np.array(ssd.squareform(pdist(z,dif)),dtype=np.float64)
    diffz[index]=-diffz[index]
    dsqr=pow(diffx,2.0)+pow(diffy,2.0)+pow(diffz,2.0)
    d=np.sqrt(dsqr[index])
    c=pairwise_contact_matrix[index]
    #c=np.array(expected_contact_matrix[index],dtype=np.float64)
    coef=a*(c-b*pow(dsqr[index],a/2.0))/(dsqr[index])
    gradx=np.zeros([n,n])
    grady=np.zeros([n,n])
    gradz=np.zeros([n,n])
    gradx[index]=coef*diffx[index]
    grady[index]=coef*diffy[index]
    gradz[index]=coef*diffz[index]
    nabla_data=np.array([   np.sum(gradx,0)-np.sum(gradx,1),
                            np.sum(grady,0)-np.sum(grady,1),
                            np.sum(gradz,0)-np.sum(gradz,1)])
    e_data=np.sum(a*c*np.log(d)-b*pow(d,a))

    #R_1,R_2=srvf.roughness2(curve,t)
    #e=lambda[0]*e_data+lambda[1]*R_1+lambda[2]*R_2
    return e_data,nabla_data
def analytical_data_gradient(curve,pairwise_contact_matrix,a,b,offdiag=1,idxmiss=None,idxlist=None):
    '''
    compute log-likelihood gradient analytically
    '''
    if idxlist is not None:
        idxmiss = None
    if idxmiss is None:
        idxmiss=np.array([],dtype='int')
    X=np.delete(curve,idxmiss,1)
    C=np.delete(pairwise_contact_matrix,idxmiss,0)
    C=np.delete(C,idxmiss,1)
    [p,n]=np.shape(X)
    
    if idxlist is None:
        idxlist=listIndices(n,offdiag)
    
    differences=np.array(X[:,idxlist[:,0]]-X[:,idxlist[:,1]],dtype=DTYPE)
    #np.nansum(differences)
    squared_distances=DTYPE(0.0);
    for ii in range(p):
        squared_distances+=differences[ii,:]*differences[ii,:]
    distances=np.sqrt(squared_distances)
    #np.nansum(distances)
    a=DTYPE(a)
    b=DTYPE(b)
    bda=b*pow(distances,a)
    #ind=np.tril(np.ones((n,n))==1,-offdiag)
    coef=(C[idxlist[:,0],idxlist[:,1]]-bda)/squared_distances
    gradE_data=np.zeros((p,n),dtype=DTYPE)
    for ii in range(p):
        gradmatrix=np.zeros((n,n),dtype=DTYPE)
        gradmatrix[idxlist[:,0],idxlist[:,1]]=coef*differences[ii,:]
        gradE_data[ii,:]=-a*(np.sum(gradmatrix,0).T-np.sum(gradmatrix,1))

    if len(idxmiss)>0:
        n+=len(idxmiss)
        temp=np.zeros((p,n),dtype=DTYPE)
        list0=np.array(range(n),dtype='int')
        list0=np.delete(list0,idxmiss)
        temp[:,list0]=gradE_data
        gradE_data=temp
    E_data=np.sum(a*C[idxlist[:,0],idxlist[:,1]]*np.log(distances)-bda,dtype=DTYPE)

    return E_data,gradE_data
def create_S1(n):
    '''
    create sparce matrix to quickly compute first order roughness gradient
    '''
    S1=csr_matrix(([2,2],([0,n-1],[0,n-1])),[n,n])
    S1+=csr_matrix((4*np.ones(n-1-1),(range(1,n-1),range(1,n-1))),[n,n])
    S1+=csr_matrix((-2*np.ones(n-1),(range(0,n-1),range(1,n))),[n,n])
    S1+=csr_matrix((-2*np.ones(n-1),(range(1,n),range(0,n-1))),[n,n])
    return S1
def create_S2(n):
    '''
    create sparce matrix to quickly compute second order roughness gradient
    '''
    S2=csr_matrix(([2,-4,-4,10],([0,0,1,1],[0,1,0,1])),[n,n])
    S2+=csr_matrix(([2,-4,-4,10],([n-1,n-2,n-1,n-2],[n-1,n-1,n-2,n-2])),[n,n])
    S2+=csr_matrix((2*np.ones(n-2),(range(2,n),range(n-2))),[n,n])
    S2+=csr_matrix((2*np.ones(n-2),(range(n-2),range(2,n))),[n,n])
    S2+=csr_matrix((12*np.ones(n-4),(range(2,n-2),range(2,n-2))),[n,n])
    S2+=csr_matrix((-8*np.ones(n-3),(range(1,n-2),range(2,n-1))),[n,n])
    S2+=csr_matrix((-8*np.ones(n-3),(range(2,n-1),range(1,n-2))),[n,n])
    return S2
    #S=term_weights[1]*S1+term_weights[2]*S2

def energy_roughness(curve):
    '''
    compute roughness 1rst and 2nd order measure
    '''
    E1=0;
    [d,n]=np.shape(curve)
    for ii in range(1,n):
        v=curve[:,ii]-curve[:,ii-1]
        E1+=sum(v*v)
    E2=0;
    for ii in range(1,n-1):
        v=curve[:,ii-1]-2*curve[:,ii]+curve[:,ii+1]
        E2+=sum(v*v)
    return E1,E2

def energy_shape_prior(curve,shape_model):
    '''
    compute energy and gradient of the shape prior
    '''
    curve,m=srvf.center_curve(curve)
    shape_model,m=srvf.center_curve(shape_model)
    shape_model=shape_model*np.sqrt(np.sum(pow(curve,2)))/np.sqrt(np.sum(pow(shape_model,2)))
    shape_model,R=srvf.find_best_ortho(curve,shape_model)
    diff=curve-shape_model
    e_shape=np.sum(pow(diff,2))
    gradient_shape=2*diff
    return e_shape,gradient_shape
def analytical_roughness_gradient(curve,S,term_weights):
    '''
    analytically compute the gradient of the roughness using the sparce S
    matrix.
    '''
    return np.array(S*np.matrix(curve.T).T)

def energy_scale_invariant_first_roughness_gradient(curve,S=None):
    '''
    scale invariant roughness penalty gradient
    '''
    [d,n]=np.shape(curve)
    if S is None:
        S=create_S1(n)
    w=np.zeros(n-1)
    diff=np.zeros([n-1,d])
    term=np.zeros([n,d])
    s=0
    for i in range(1,n):
        diff[i-1,:]=curve[:,i]-curve[:,i-1]
        w[i-1]=np.linalg.norm(diff[i-1,:])
        s+=w[i-1]*w[i-1]
        if i==1:
            term[i-1,:]=-diff[i-1,:]/w[i-1]
        else:
            term[i-1,:]=diff[i-2,:]/w[i-2] -diff[i-1]/w[i-1]
            #s=s+pow(np.linalg.norm(diff[i-1,:]-diff[i-2,:]),2.0)
    term[n-1,:]=diff[n-2,:]/w[n-2]
    L=sum(w)
    LL=L*L
    E=s/LL
    gradE=np.array(S*np.matrix(curve.T)).T/LL-2*term.T*s/(LL*L)
    return E,gradE;

def energy_scale_invariant_second_roughness_gradient(curve,S=None):
    '''
    scale invariant roughness penalty gradient
    '''
    [d,n]=np.shape(curve)
    if S is None:
        S=create_S2(n)
    w=np.zeros(n-1)
    diff=np.zeros([n-1,d])
    term=np.zeros([n,d])
    s=0
    for i in range(1,n):
        diff[i-1,:]=curve[:,i]-curve[:,i-1]
        w[i-1]=np.linalg.norm(diff[i-1,:])
        if i==1:
            term[i-1,:]=-diff[i-1,:]/w[i-1]
        else:
            term[i-1,:]=diff[i-2,:]/w[i-2] -diff[i-1]/w[i-1]
            s=s+pow(np.linalg.norm(diff[i-1,:]-diff[i-2,:]),2.0)
    term[n-1,:]=diff[n-2,:]/w[n-2]
    L=sum(w)
    LL=L*L
    E=s/LL
    gradE=np.array(S*np.matrix(curve.T)).T/LL-2*term.T*s/(LL*L)
    return E,gradE;
class opt_E():
    '''
    Primary optimization manager for simba3d

    This class manages defaults and data storage for the results. Basically,
    it is a wraper for the scipy optimizer, but it stores results and manages
    input parameters.

    '''
    def __init__(self,data,parameters=None,options=None):
        '''
        initialize the energy optimization object
        '''

        self.n=10;  # number of node points this should be overided by user data
        self.d=3;   # this might be overided by user data, but right now I think it only works for dimension d=3

        # set parameters and defaults #########################################
        # parameters
        if parameters is None: # if parameters were not passed, then use default
            parameters={}
            print("Default parameter settings will be used:")
        # a parameter default check
        if 'a' not in parameters: # check for a
            parameters['a']=DTYPE(-3); # set default
            print("Default a ="+str(parameters['a']))
        self.a=DTYPE(parameters['a'])
        # b parameter default check
        if 'b' not in parameters: # check for b
            parameters['b']=DTYPE(1);  # set default (b is not identifiable with scale of points)
            print( "Default b ="+str(parameters['b']))
        self.b=DTYPE(parameters['b'])
        # lamina_radius parameter default check
        if 'lamina_radius' not in parameters: # check for b
            parameters['lamina_radius']=DTYPE(1);  # set default (b is not identifiable with scale of points)
            print( "Default lamina_radius ="+str(parameters['lamina_radius']))
        self.lamina_radius=DTYPE(parameters['lamina_radius'])
        # term_weights parameter default check
        if 'term_weights' not in parameters: # check for term weights
            parameters['term_weights']= {
                                 }
        # set the supported penalties weights
        self.term_weights=parameters['term_weights']
        if "data" not in self.term_weights:
            self.term_weights['data']=DTYPE(1.0) # weight for the data term
            #print "Default data weight  ="
            #print self.term_weights['data']
        if "h1" not in self.term_weights: # scale and intensity invariant "first order"-like penalty
            if "uniform spacing" not in self.term_weights:
                self.term_weights['h1']=DTYPE(0.0 )   # weight for the h1 penalty
                print("Default population h1 weight  =")
                print(self.term_weights['h1'])
            else:
                self.term_weights['h1']=DTYPE(self.term_weights['uniform spacing'])
        if "h2" not in self.term_weights: # scale and intensity invariant "second order"-like penalty
            if "smoothing" not in self.term_weights:
                self.term_weights['h2']=DTYPE(0.0)    # weight for the h2 penalty
                print("Default population h2 weight  =")
                print(self.term_weights['h2'])
            else:
                self.term_weights['h2']=DTYPE(self.term_weights['smoothing'])
        if "h4" not in self.term_weights: # scale and intensity invariant "second order"-like penalty
            if "lamina" not in self.term_weights:
                self.term_weights['h4']=DTYPE(0.0)    # weight for the h2 penalty
                print( "Default lamina h4 weight  =")
                print( self.term_weights['h4'])
            else:
                self.term_weights['h4']=DTYPE(self.term_weights['lamina'])
        if "population prior" not in self.term_weights:
            self.term_weights['population prior']=DTYPE(0.0)    # weight for the population matrix prior
            print( "Default population prior weight  =")
            print( self.term_weights['population prior'])
        if "shape prior" not in self.term_weights:
            self.term_weights['shape prior']=DTYPE(0.0)    # weight for the shape prior
            print( "Default shape prior weight  =")
            print( self.term_weights['shape prior'])

        # set defaults for unsuported penalties
        # scale and intensity dependent penalties are no longer supported
        if "firstroughness" not in self.term_weights: # depends on scale of curve and total number of interaction in data
            self.term_weights['firstroughness']=0.0    # weight for the second order roughness
            #print "Default firstroughness' weight  ="
            #print self.term_weights['firstroughness']
        if "secondroughness" not in self.term_weights:
            # depends on scale of curve and total number of interaction in data
            self.term_weights['secondroughness']=0.0    # weight for the second order roughness
            #print "Default secondroughness weight  ="
            #print self.term_weights['secondroughness']
        if "scaledfirstroughness" not in self.term_weights:
            self.term_weights['scaledfirstroughness']=0.0    # weight for the second order roughness invaraint scale of curve but dependent on intesity of interaction
            #print "Default scaledfirstroughness weight  ="
            #print self.term_weights['scaledfirstroughness']
        if "scaledsecondroughness" not in self.term_weights:
            self.term_weights['scaledsecondroughness']=0.0    # weight for the second order roughness invaraint scale of curve but dependent on intesity of interaction
            #print "Default scaledsecondroughness weight  ="
            #print self.term_weights['scaledsecondroughness']


        # set data and defaults ###############################################
        # pairwise contact matrix
        if 'sparse_pairwise_contact_matrix' in data:
          print( data['sparse_pairwise_contact_matrix'])
          if 'pairwise_contact_matrix' not in data:
              # this is done the dumb way for now
              data['pairwise_contact_matrix'] = data['sparse_pairwise_contact_matrix'].toarray()
              data['pairwise_contact_matrix'] = data['pairwise_contact_matrix']+data['pairwise_contact_matrix'].transpose()

        if 'pairwise_contact_matrix' not in data: # if no pairwise contact matrix is given
            data['pairwise_contact_matrix'] =0 # set to zero
            if self.term_weights['data']!=0:
                self.term_weights['data']=0
                print("\nNo pairwise interaction matrix proivded for data term!\n")
        else:
            self.pairwise_contact_matrix=self.term_weights['data']*data['pairwise_contact_matrix']
            (m,n)=np.shape(self.pairwise_contact_matrix)
            self.n=min([m,n])
        # population contact matrix
        if 'sparse_population_contact_matrix' in data:
          print( data['sparse_population_contact_matrix'])
          if 'population_contact_matrix' not in data:
              # this is done the dumb way for now
              data['population_contact_matrix'] =data['sparse_population_contact_matrix'].toarray()
              data['population_contact_matrix'] = data['population_contact_matrix']+data['population_contact_matrix'].transpose()
        if 'population_contact_matrix' not in data: # if no population contact matrix is given
            self.population_contact_matrix=0
            if self.term_weights['population prior']!=0:
                self.term_weights['population prior']=0
                print("No population matrix provided for data prior!\n")
        else:
            self.population_contact_matrix=data['population_contact_matrix']
            (m,n)=np.shape(self.pairwise_contact_matrix)
            self.n=min([m,n])
        # shape prior model data
        if 'prior_shape_model' not in data:
            self.prior_shape_model=None
            if self.term_weights['shape prior']!=0:
                self.term_weights['shape prior']=0
                print( "No prior shape model provided for shape prior!\n")
        else:
            self.prior_shape_model=data['prior_shape_model']
        if 'initialized_curve' not in data:
            if 'seed' in options:
                rng=np.random.RandomState(options['seed'])
                data['initialized_curve']=rng.normal(0,1,self.n*self.d).reshape([self.d,self.n])
            else:
                data['initialized_curve']=np.random.rand(self.n*self.d).reshape([self.d,self.n])
            print("No initialized_curve provided, using randomized curve!\n")
            print('Number of points: '+str(self.n))
            print('Curve dimension: '+str(self.d))
        self.initialized_curve=data['initialized_curve']

        # check that initialized curve has correct dimension
        (mm,nn)=np.shape(self.initialized_curve)
        if (nn!=n): # if dimension does not match then interpolate
            print("\nInitialized curve has "+str(nn)+" nodes, "+
                "but data has "+str(n)+" nodes!")
            print("Interpolating initialized curve to match data dimension!\n")
            t0=np.linspace(0,1,nn)
            t=np.linspace(0,1,n)
            self.initialized_curve=srvf.interp(t0,self.initialized_curve,t)

        # time parameterization
        if 't' not in data:
            data['t']=np.linspace(0,1,self.n)
        self.t=data['t']


        '''index is no longer used, it has been replaced with pairwise_combinations
        # index of contact matrix to use in calculation of log-likelihood term
        if 'index' not in data:
            data['index']=np.triu(np.ones((self.n,self.n)),1)==1
            # print "No index provided for log-likelihood term!\n"
        '''

        # temporarily retroactively assign equvialent index until a fix can be made
        self.index=np.triu(np.ones((self.n,self.n)),1)==1

        # population matrix prior  adjustment #################################
        Ndata=sum(data['pairwise_contact_matrix'][self.index] )+10e-16
        Nprior=1
        if self.term_weights['population prior']:
            Nprior=sum(self.population_contact_matrix[self.index] )+10e-16
            self.pairwise_contact_matrix+=self.term_weights['population prior']*(Ndata/Nprior)*self.population_contact_matrix

        if 1:
            self.b=(self.term_weights['data']+self.term_weights['population prior']*Ndata/Nprior)*self.b
            self.Ndata=Ndata
        else:
            self.b=(self.term_weights['population prior']*Ndata/Nprior)*self.b
            # adjust weights
            self.term_weights['scaledfirstroughness']*=Ndata
            self.term_weights['scaledsecondroughness']*=Ndata
            self.term_weights['data']*=Ndata # scale by amount of data observed to keep tuning parameters comparible
            self.Ndata=1
        self.pairwise_contact_matrix[np.identity(n)==1]=0;

        # set missing indexes #################################################
        # threshold value for the row sum to consider the index as a missing value
        # this is used only if the user does not manually specify the missing indexes
        if 'missing_rowsum_threshold' not in data:
            data['missing_rowsum_threshold']=0
        self.missing_rowsum_threshold=data['missing_rowsum_threshold']
        # Index the nodes that are missing
        # specifies which nodes do not have enough observed data
        if 'index_mississing' not in data:
            ans=np.where((np.sum(self.pairwise_contact_matrix,0)<=self.missing_rowsum_threshold))
            data['index_mississing']=ans[0]
        self.index_mississing=data['index_mississing']
        # How many off from the diagonal can we start reading
        if 'off_diagonal' not in data:
            data['off_diagonal']=1
        self.off_diagonal=data['off_diagonal']
        # A list of all pairwise combinations added in the data term
        # This speeds up the calculation considerably to store it.
        if 'pairwise_combinations' not in data:
            data['pairwise_combinations']=listIndices(self.n-len(data['index_mississing']),self.off_diagonal)
        self.pairwise_combinations=data['pairwise_combinations']


        # set options #########################################################
        if 'maxitr' not in options:
            options['maxitr']=100000
            print( 'Maximum number of iterations: '+str(options['maxitr']))
        self.maxitr=options['maxitr']
        if 'display' not in options:
            options['display']=False
            print( 'display output: '+str(options['display']))
        self.display=options['display']
        if 'store' not in options:
            options['store']=False
            print( 'store iterative values: '+str(options['store']))
        self.store=options['store']
        if 'method' not in options:
            #options['method']='Nelder-Mead'
            options['method']='BFGS'
            #options['method']='CG'
            #options['method']='Newton-CG'
            print( 'Optimization Method: '+options['method'])
        self.method= options['method']
        if 'gradient tolerance' not in options:
            options['gradient tolerance']=1e-5   # default for MatLab and Scipy Optimize
            print( 'store iterative values: '+str(options['store']))
        self.gtol=options['gradient tolerance']
        if 'minimum energy' not in options:
            options['minimum energy']=-np.inf  # stop if minimum energy is obtained
        self.minimum_energy= options['minimum energy']
        # store useful variables ##############################################
        self.shp=np.shape(self.initialized_curve)
        print( self.shp)
        self.XK=None;
        if self.store:
            self.XK=np.zeros([self.maxitr,np.prod(self.shp)])
        self.e_current=None;
        self.itr=0;
        self.e=np.zeros(self.maxitr)
        self.n=max(self.shp)
        self.S1=0
        if (self.term_weights["firstroughness"]!=0)|(self.term_weights["scaledfirstroughness"]!=0)|(self.term_weights["h1"]!=0):
            self.S1=create_S1(self.n)
        self.S2=0
        if (self.term_weights["secondroughness"]!=0)|(self.term_weights["scaledsecondroughness"]!=0):
            self.S2=create_S2(self.n)
        self.S=self.term_weights["firstroughness"]*self.S1+self.term_weights["secondroughness"]*self.S2

        # initalize results variables
        self.result=None;
        self.estimated_curve=None
        self.ax=None;
        self.line=None;

        self.aniitr=0

        # rescale initialized curve
        if self.term_weights['data']!=0:
            self.initialized_curve=scale_initialized_curve(self.initialized_curve,self.pairwise_contact_matrix,self.a,self.b,self.index,self.index_mississing)

        self.start_time=None
        self.end_time=None
        self.computation_time=None
        # run optimizer when initilaized
        if 0:
            self.run();
    def energy_jacobian(self,curve,*args):
        '''
        analytically compute energy and jacobian energy
        '''
        shp=self.shp
        # compute data term
        e_data=0
        gradient_data=0
        if (self.term_weights['data']!=0)|(self.term_weights['population prior']!=0):
            e_data,gradient_data=analytical_data_gradient(curve.reshape(shp),self.pairwise_contact_matrix,self.a,self.b,offdiag=self.off_diagonal,idxmiss=self.index_mississing,idxlist=self.pairwise_combinations)
            gradient_data=gradient_data.flatten()

        # compute shape terms
        e_shape=0 # initialize the energy
        gradient_shape=0 # initialize the gradient
        if (self.term_weights["shape prior"]!=0):
            e_shape,gradient_shape=energy_shape_prior(curve.reshape(shp),self.prior_shape_model)
            gradient_shape=gradient_shape.flatten()

        # compute improved roughness terms (h1 and h2 penallties)
        e_h1=0   # initialize the energy
        gradient_h1=0 # initialize the gradient
        if (self.term_weights["h1"]!=0):
            [e_h1,gradient_h1]=h1penalty(curve.reshape(shp),self.S1)
            gradient_h1=gradient_h1.flatten()
        e_h2=0   # initialize the energy
        gradient_h2=0 # initialize the gradient
        if (self.term_weights["h2"]!=0):
            [e_h2,gradient_h2]=h2penalty(curve.reshape(shp))
            gradient_h2=gradient_h2.flatten()
        # compute the lamina penalty
        e_h4=0   # initialize the energy
        gradient_h4=0 # initialize the gradient
        if (self.term_weights["h4"]!=0):
            [e_h4,gradient_h4]=h4penalty(curve.reshape(shp),self.lamina_radius)
            gradient_h4=gradient_h4.flatten()
        # compute unsuported roughness terms (these have been replaced with better ones)
        gradient_roughness=0 # initialize the gradient
        e_r1=0 # initialize the energy
        e_r2=0 # initialize the energy
        if (self.term_weights["firstroughness"]!=0)|(self.term_weights["secondroughness"]!=0):
            e_r1,e_r2=energy_roughness(curve.reshape(shp))
            gradient_roughness= np.array(self.S*np.matrix(curve.reshape(shp)).T).T
            gradient_roughness= gradient_roughness.flatten()
        # compute unsuported scaled roughness (these have been replaced with better ones)
        e_sr1=0 # initialize the energy
        gradient_scaled_first_roughness=0 # initialize the gradient
        if (self.term_weights["scaledfirstroughness"]!=0):
            e_sr1,gradient_scaled_first_roughness=energy_scale_invariant_first_roughness_gradient(curve.reshape(shp),self.S1)
            gradient_scaled_first_roughness= gradient_scaled_first_roughness.flatten()
        e_sr2=0 # initialize the energy
        gradient_scaled_second_roughness=0 # initialize the gradient
        if (self.term_weights["scaledsecondroughness"]!=0):
            e_sr2,gradient_scaled_second_roughness=energy_scale_invariant_second_roughness_gradient(curve.reshape(shp),self.S2)
            gradient_scaled_second_roughness= gradient_scaled_second_roughness.flatten()

        # compute total energy
        e=DTYPE(0.0);
        e-=e_data/self.Ndata # data term
        e+=self.term_weights["h1"]*e_h1# scaled first order roughness terms
        e+=self.term_weights["h2"]*e_h2# scaled second order roughness terms
        e+=self.term_weights["shape prior"]*e_shape# first order roughness terms
        e+=self.term_weights["h4"]*e_h4# lamina penalty
        # add unsupported penalties
        e+=self.term_weights["firstroughness"]*e_r1# first order roughness terms
        e+=self.term_weights["secondroughness"]*e_r2# second order roughness terms
        e+=self.term_weights["scaledfirstroughness"]*e_sr1# scaled first order roughness terms
        e+=self.term_weights["scaledsecondroughness"]*e_sr2# scaled second order roughness terms

        self.e_current=e; # store the current energy

        # compute total gradient
        gradient=DTYPE(0.0);
        gradient-=gradient_data/self.Ndata# data term
        gradient+=self.term_weights["h1"]*gradient_h1 # the h1 term
        gradient+=self.term_weights["h2"]*gradient_h2 # the h2 term
        gradient+=self.term_weights["h4"]*gradient_h4 # the h4 term
        gradient+=self.term_weights["shape prior"]*gradient_shape # shape prior (needs scaled)
        # add unsupported penalties
        gradient+=gradient_roughness # roughness terms
        gradient+=self.term_weights["scaledfirstroughness"]*gradient_scaled_first_roughness # scaled second order roughness term
        gradient+=self.term_weights["scaledsecondroughness"]*gradient_scaled_second_roughness # scaled second order roughness term
        return e,gradient

    def accelerate_gradient_ascent_algorithm(self,cbfun=None,options=None,delta0=None):
        if options is None:
            options={   'maxiter':self.maxitr,
                'gtol':1e-5,
                }
        backtracking=True
        curve=self.initialized_curve
        tol=options['gtol']
        maxitr=options['maxiter']
        if delta0 is None:
            delta0=1

        ii=0

        deltamin=1e-15
        tau=.05
        c=.1

        E=np.zeros((maxitr,1))
        normgrad=np.zeros((maxitr,1))

        self.e_current,grad=self.energy_jacobian(curve)
        E[ii]=self.e_current

        normgrad[ii]=np.linalg.norm(grad)
        #print 'Energy: '+str(E[ii])+'\tIteratiuon:'+str(ii)

        current_curve=curve
        previous_curve=current_curve.flatten()
        if cbfun is not None:
            cbfun(current_curve.flatten())
        delta=delta0

        # first to a regular gradient descent step
        candidate_curve=previous_curve-delta*grad
        candidate_e,candidate_grad=self.energy_jacobian(candidate_curve)
        if backtracking: # Armijo Backtracking
            ag_condition=E[ii]-candidate_e
            while (ag_condition<delta*c*pow(normgrad[ii],2.0)) & (delta>deltamin):
                delta=tau*delta
                candidate_curve=previous_curve-delta*grad
                candidate_e,candidate_grad=self.energy_jacobian(candidate_curve)
                ag_condition=E[ii]-candidate_e
        current_curve=candidate_curve
        grad=candidate_grad
        E[ii+1]=candidate_e
        self.e_current=E[ii+1]
        normgrad[ii+1]=np.linalg.norm(grad)
        ii+=1
        delta=delta0

        X=current_curve.reshape(self.shp)
        if cbfun is not None:
            cbfun(current_curve.flatten())
        while (normgrad[ii]>tol) & ((ii<(maxitr-1))|(self.n<200)):
            # Accelerated gradient descent step
            Y=current_curve+(ii-1.0)/(ii+2.0)*(current_curve-previous_curve)
            e_current,grad=self.energy_jacobian(Y)
            candidate_curve=Y-delta*grad

            #if backtracking & (ii<10): # just do it for the first few iterations
            if backtracking & ((ii<10)): # just do it for the first few iterations
            #if backtracking:
                delta=delta0
                candidate_e,candidate_grad=self.energy_jacobian(candidate_curve)
                while (candidate_e> (e_current+np.trace(np.matrix(grad).T*np.matrix(candidate_curve-Y))+1.0/(2.0*delta)*pow(np.linalg.norm(candidate_curve-Y),2.0) ))& (delta>deltamin):
                    delta=tau*delta
                    candidate_curve= Y-delta*grad
                    candidate_e,candidate_grad=self.energy_jacobian(candidate_curve)
                E[ii+1] = candidate_e
                #print 'step size: '+str(delta)
            #update
            previous_curve=current_curve
            current_curve=candidate_curve
            normgrad[ii+1]=np.sqrt(np.sum(pow(grad,2.0)))
            self.e_current=e_current
            if cbfun is not None:
                cbfun(current_curve.flatten())
            ii+=1
            #print 'Energy: '+str(E[ii])+'\tIteratiuon:'+str(ii)+'\tnorm of gradient:'+str(normgrad[ii]) +'\tdelta:'+str(delta)
            X=current_curve.reshape(self.shp)
        return X,E
    def callbackdispstore(self,xk):
        '''
        display iteration during callback and store iterative result
        '''
        print( "Iteration "+str(self.itr)+" e="+str(self.e_current))

        if xk is not None:
            self.estimated_curve=xk.reshape(self.shp) ;
        self.e[self.itr]=self.e_current
        self.XK[self.itr]=xk
        self.itr+=1;

        if self.e_current<self.minimum_energy:
            raise stop_optimizing_exception()
    def callbackdisp(self,xk):
        '''
        display iteration during callback
        '''
        print( "Iteration "+str(self.itr)+" e="+str(self.e_current))
        self.e[self.itr]=self.e_current
        if xk is not None:
            self.estimated_curve=xk.reshape(self.shp) ;
        self.itr+=1;
        if self.e_current<self.minimum_energy:
            raise stop_optimizing_exception()
    def callback(self,xk):
        '''
        lightweight callback without storing iterative result or printing iterative value
        '''
        #print "Iteration "+str(self.itr)+" e="+str(self.e_current)
        if xk is not None:
            self.estimated_curve=xk.reshape(self.shp) ;
        self.e[self.itr]=self.e_current
        self.itr+=1;
        if self.e_current<self.minimum_energy:
            raise stop_optimizing_exception()
    def callbackstore(self,xk):
        '''
        store iterative result without displaying during iterations
        '''
        #print "Iteration "+str(self.itr)+" e="+str(self.e_current)
        if xk is not None:
            self.estimated_curve=xk.reshape(self.shp) ;
        self.e[self.itr]=self.e_current
        self.XK[self.itr]=xk
        self.itr+=1;
        if self.e_current<self.minimum_energy:
            raise stop_optimizing_exception()
    def run(self):

        options={   'maxiter':self.maxitr,
                    'disp':self.display,
                    #'eps':1.5e-7,
                    #'norm':np.inf,
                    'gtol':self.gtol,
                    }
        if self.display:
            if self.store:
                cbfun=self.callbackdispstore
            else:
                cbfun=self.callbackdisp
        else:
            if self.store:
                cbfun=self.callbackstore
            else:
                cbfun=self.callback
        self.start_time=time.time()
        if (self.method== 'default'):
            self.method='BFGS'
        try:
            if (self.method== 'AGM'):
                self.accelerate_gradient_ascent_algorithm(cbfun,options)
            else:
                self.result=minimize(self.energy_jacobian,self.initialized_curve.flatten(),args=(self),method=self.method,tol=None,jac=True,options=options,callback=cbfun)
                #self.estimated_curve=self.result.x.reshape(self.shp)
        except stop_optimizing_exception:
            self.result=None;

        self.end_time=time.time()
        self.computation_time=self.end_time-self.start_time
        self.e=self.e[range(self.itr)]
        if self.store:
            self.XK=self.XK[range(self.itr)]
    def just_energy(self,curve,*args):
         e,gradient=self.energy_jacobian(curve,*args)
         return e;
    def just_gradient(self,curve,*args):
         e,gradient=self.energy_jacobian(curve,*args)
         return gradient;
    def check_jacobian(self):
        return check_grad(self.just_energy,self.just_gradient,self.initialized_curve.flatten())
    def get_data(self,i):
        '''
        get_data function call for creating animation
        '''
        i=int(np.floor((np.double(self.itr)/np.double(self.aniitr))*i))
        return self.XK[i].reshape(self.shp)
