"""



@author: Michael Rosenthal
"""

# import the gradients


import numpy as np

def gradient_h1(
				penalty_weight,
				length,
				sum_adjacent_distance,
				sum_of_squared_adjacent_distance,
				x_differences,
				y_differences,
				z_differences,
				distances,
				x_gradient,
				y_gradient,
				z_gradient
				):
	'''
	uniform spacing penalty
	'''
	L_squared=sum_adjacent_distance*sum_adjacent_distance;
	L_cubed=sum_adjacent_distance*sum_adjacent_distance*sum_adjacent_distance;
	S=sum_of_squared_adjacent_distance;
	term0=2.0/L_squared;
	coef=(length-1)*(term0-2.0*S/(L_cubed*distances[0]));
	x_gradient[0]+=penalty_weight*coef*x_differences[0];
	y_gradient[0]+=penalty_weight*coef*y_differences[0];
	z_gradient[0]+=penalty_weight*coef*z_differences[0];

	for ii in range(1,length-1):
		adj_ind=int((ii+2)*(ii+1)/2 -1);#get the adjacent index
		adj_ind0=int((ii+1)*(ii)/2 -1);#get the previous adjacent index
		coef0=coef;
		#coef=(length-1)*(term0+2.0*S/(L_cubed*distances[adj_ind0]));
		x_gradient[ii]-=penalty_weight*x_differences[adj_ind0]*coef0;
		y_gradient[ii]-=penalty_weight*y_differences[adj_ind0]*coef0;
		z_gradient[ii]-=penalty_weight*z_differences[adj_ind0]*coef0;
		coef=(length-1)*(term0-2.0*S/(L_cubed*distances[adj_ind]));
		x_gradient[ii]+=penalty_weight*x_differences[adj_ind]*coef;
		y_gradient[ii]+=penalty_weight*y_differences[adj_ind]*coef;
		z_gradient[ii]+=penalty_weight*z_differences[adj_ind]*coef;
	ii=length-2;
	adj_ind=int((ii+2)*(ii+1)/2 -1);#get the adjacent index
	coef=(length-1)*(term0-2.0*S/(L_cubed*distances[adj_ind]));
	x_gradient[length-1]-=penalty_weight*x_differences[adj_ind]*coef;
	y_gradient[length-1]-=penalty_weight*y_differences[adj_ind]*coef;
	z_gradient[length-1]-=penalty_weight*z_differences[adj_ind]*coef;