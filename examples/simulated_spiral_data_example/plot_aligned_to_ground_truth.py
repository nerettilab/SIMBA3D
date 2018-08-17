# -*- coding: utf-8 -*-
"""
Compute the total computation time from a list of simba3d output report files.

@author: Michael Rosenthal
"""
import os


import numpy as np

import simba3d.mp_manager as mp
import simba3d.srvf_open_curve_Rn as srvf
import simba3d.plotting_tools as pt

#%% load ground truth curve
ground_truth_curve=np.load("data/ground_truth_curves/double_spiral_400.npy")
#remove translation and scale of the ground truth curve
center_curve,mu=srvf.center_curve(ground_truth_curve) # remove translation
center_curve,scale=srvf.scale_curve(center_curve) # remove scale

#%% load the output report file curve
report_file ='results/simulated_data_b7b1fe84-8b80-4bb4-9af8-8f0f9e7aa3b7_400.npz'
summary=mp.load_result(report_file)
curves=summary['X_evol']
#remove translation and scale of the ground truth curve
curve,mu=srvf.center_curve(curves[-1]) # remove translation
curve,scale=srvf.scale_curve(curve) # remove scale
#align estimated curve to the ground truth curve
curve,rot=srvf.find_best_rotation(center_curve,curve)# align to the first curve

# plot the estimated curve
d,n=np.shape(curve);
t=np.linspace(0,1,n)
fig1,ax1=pt.plot_curve(curve,t)
pt.plot_pts(curve,t,ax=ax1,fig=fig1)   
fig1.suptitle('Estimated Curve',fontsize=18)
fig1.savefig('images/estimated_curve.png');

# add the ground truth curve on top
d,n=np.shape(center_curve);
t=np.linspace(0,1,n)
fig1,ax1=pt.plot_curve(center_curve,ax=ax1,fig=fig1)
pt.plot_pts(center_curve,ax=ax1,fig=fig1)   
fig1.suptitle('Estimated Curve with Ground Truth',fontsize=18)
fig1.savefig('images/estimated_and_groundtruth_curve.png');

# print the weight settings
print("uniform spacing penalty:"+str(summary['weight_uniform_spacing']))
print("smoothing penalty:"+str(summary['weight_smoothing']))
print("population prior penalty:"+str(summary['weight_population_prior']))

