"""



@author: Michael Rosenthal
"""

# import the gradients


import numpy as np

from simba3d.pairwise_computations import triu_ij_to_ind

def gradient_poisson(   penalty_weight,
						contact_row_index,
						contact_col_index,
						contact_value,
						poisson_a_parameter,
						poisson_b_parameter,
						contact_sum,
						pairwise_difference_x,
						pairwise_difference_y,
						pairwise_difference_z,
						pairwise_distance,
						x_gradient,
						y_gradient,
						z_gradient,
						series
						):
	'''
	'''
	series[0]=0.0;
	d=0.0
	differences=[]
	for ind in range(len(contact_value)):
		ind2=triu_ij_to_ind(contact_row_index[ind],contact_col_index[ind]);
		differences.append(pairwise_distance[ind2]);
	for ind in range(len(contact_value)):
		ind2=triu_ij_to_ind(contact_row_index[ind],contact_col_index[ind]);
		series[0]-=poisson_a_parameter*contact_value[ind]*np.log(pairwise_distance[ind2]);
		series[0]+=poisson_b_parameter*pow(pairwise_distance[ind2],poisson_a_parameter);
		scale=-penalty_weight*poisson_a_parameter*(contact_value[ind]-poisson_b_parameter*pow(pairwise_distance[ind2],poisson_a_parameter));
		scale/=(pairwise_distance[ind2]*pairwise_distance[ind2]);
		scale/=contact_sum;
		x_gradient[contact_row_index[ind]]+=scale*pairwise_difference_x[ind2];
		y_gradient[contact_row_index[ind]]+=scale*pairwise_difference_y[ind2];
		z_gradient[contact_row_index[ind]]+=scale*pairwise_difference_z[ind2];
		x_gradient[contact_col_index[ind]]-=scale*pairwise_difference_x[ind2];
		y_gradient[contact_col_index[ind]]-=scale*pairwise_difference_y[ind2];
		z_gradient[contact_col_index[ind]]-=scale*pairwise_difference_z[ind2];
	series[0]/=contact_sum;

