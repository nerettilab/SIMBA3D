# -*- coding: utf-8 -*-
"""

Simulate 

@author: Michael Rosenthal
"""
display_the_result=False
""" to display the results... -->
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
display_the_result=True
<-- ... move this line to the other arrow above """ 


import numpy as np
import os
import scipy.io

from simba3d.optimizer import compute_pairwise_distance,listIndices
import simba3d.srvf_open_curve_Rn as srvf
import simba3d.plotting_tools as pt

#%% load the ground truth curve
ground_truth_curve=np.load("data/ground_truth_curves/double_spiral_curve.npy")
t=np.linspace(0,1,len(ground_truth_curve)) # parametrized the original curve

#%% resample the ground truth curve at the desired number of points
a=-3 # the power parameter relating the # of pairwise interactions and distance
b=.01 # scaling parameter relating the # of pairwise interactions and distance
#n=1280  # number of resampled points
n=400  # number of resampled points
tn=np.linspace(0,1,n) # curve parametrization
# resample the gorund ture curve with repsect the the tn parametrization
resampled_ground_truth_curve=srvf.interp(t,ground_truth_curve,tn).T
# save the resampled ground truth curve
np.save("data/ground_truth_curves/double_spiral_"+str(n)+".npy",resampled_ground_truth_curve)

#%% simulate a data matrix according to the Varoquaux model
# compute pairwise distance
D=compute_pairwise_distance(resampled_ground_truth_curve) 
ind=listIndices(n,1) # get the upper triangular index of the matrix
# compute expect number of interactions
MU=b*pow(D,a)
# simulate from independent poisson distribution
c=np.random.poisson(MU[ind[:,0],ind[:,1]])
# construct the pairwise interaction matrix
C=np.zeros([n,n])
C[ind[:,0],ind[:,1]]=c # only the upper triangular part
C=C+C.T # make the matrix symetric
# save the simulated matrix
np.save("data/simulate_data_"+str(n)+".npy",C)


        
fig1,ax1=pt.plot_curve(resampled_ground_truth_curve,tn)
pt.plot_pts(resampled_ground_truth_curve,t,ax=ax1,fig=fig1)   
fig1.suptitle('Ground Truth Curve',fontsize=18)
fig1.savefig('images/ground_truth_curve.png');

fig2,ax2,m=pt.matshow(1.0*(C>0))
fig2.suptitle('Non-Zero Simulated Data',fontsize=18)
fig2.savefig('images/simulated_data_nonzeros.png');

fig3,ax3,m=pt.matshow(C)
fig3.suptitle('Simulated Data',fontsize=18)   
fig3.savefig('images/simulated_data.png');
        
if display_the_result:
  plt.show(block=True) # If you want to display the results uncomment the thing above
