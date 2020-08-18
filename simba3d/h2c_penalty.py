"""



@author: Michael Rosenthal
"""

# import the gradients


import numpy as np

from simba3d.pairwise_computations import triu_ind_to_ij,triu_ij_to_ind
def run_h2c_computations(
							length,
							number_of_pairs,
							h2c_radius,
							average_adjacent_distance,
							pairwise_distance,
							h2b_F,
							series
							):
	'''
	'''
	series[0]=0.0
	for ind in range(number_of_pairs):
		ii,jj=triu_ind_to_ij(ind)
		u_jk= pairwise_distance[ind]/average_adjacent_distance
		h2b_F[ind]=1.0-u_jk
		if (pairwise_distance[ind]/average_adjacent_distance<h2c_radius):
			series[0]+=pow(h2c_radius-pairwise_distance[ind]/average_adjacent_distance,4.0);
	series[0]*=2.0/np.float((length*(length-1.0)));

def gradient_h2c(
				penalty_weight,
				length,
				number_of_pairs,
				h2b_alpha,
				h2c_F,
				h2c_radius,
				average_adjacent_distance,
				pairwise_difference_x,
				pairwise_difference_y,
				pairwise_difference_z,
				pairwise_distance,
				x_gradient,
				y_gradient,
				z_gradient
				):
	'''
	pairwise repulsion penalty
	'''
	for ii in range(length):
		# partial derivative of mu
		if (ii==0):
			ind2=triu_ij_to_ind(0,1);
			partial_mu_xi_x=pairwise_difference_x[ind2]/((length-1)*pairwise_distance[ind2]);
			partial_mu_xi_y=pairwise_difference_y[ind2]/((length-1)*pairwise_distance[ind2]);
			partial_mu_xi_z=pairwise_difference_z[ind2]/((length-1)*pairwise_distance[ind2]);
		elif (ii==length-1):
			ind2=triu_ij_to_ind(ii-1,ii);
			partial_mu_xi_x=-pairwise_difference_x[ind2]/((length-1)*pairwise_distance[ind2]);
			partial_mu_xi_y=-pairwise_difference_y[ind2]/((length-1)*pairwise_distance[ind2]);
			partial_mu_xi_z=-pairwise_difference_z[ind2]/((length-1)*pairwise_distance[ind2]);
		else:
			ind2=triu_ij_to_ind(ii,ii+1);
			partial_mu_xi_x=pairwise_difference_x[ind2]/((length-1)*pairwise_distance[ind2]);
			partial_mu_xi_y=pairwise_difference_y[ind2]/((length-1)*pairwise_distance[ind2]);
			partial_mu_xi_z=pairwise_difference_z[ind2]/((length-1)*pairwise_distance[ind2]);
			ind2=triu_ij_to_ind(ii-1,ii);
			partial_mu_xi_x-=pairwise_difference_x[ind2]/((length-1)*pairwise_distance[ind2]);
			partial_mu_xi_y-=pairwise_difference_y[ind2]/((length-1)*pairwise_distance[ind2]);
			partial_mu_xi_z-=pairwise_difference_z[ind2]/((length-1)*pairwise_distance[ind2]);	
		for ind in range(number_of_pairs):
			jj,kk=triu_ind_to_ij(ind);
			if ((1.0-h2c_F[ind])<h2c_radius):
				# partial derivative of pairwise_distance
				partial_djk_xi_x=0.0;
				partial_djk_xi_y=0.0;
				partial_djk_xi_z=0.0;
				if (jj==ii):
					partial_djk_xi_x=pairwise_difference_x[ind]/pairwise_distance[ind];
					partial_djk_xi_y=pairwise_difference_y[ind]/pairwise_distance[ind];
					partial_djk_xi_z=pairwise_difference_z[ind]/pairwise_distance[ind];
				elif (kk==ii):
					partial_djk_xi_x=-pairwise_difference_x[ind]/pairwise_distance[ind];
					partial_djk_xi_y=-pairwise_difference_y[ind]/pairwise_distance[ind];
					partial_djk_xi_z=-pairwise_difference_z[ind]/pairwise_distance[ind];
				# partial derivative of ujk
				partial_ujk_xi_x=(partial_djk_xi_x+(h2c_F[ind]-1.0)*partial_mu_xi_x)/average_adjacent_distance;
				partial_ujk_xi_y=(partial_djk_xi_y+(h2c_F[ind]-1.0)*partial_mu_xi_y)/average_adjacent_distance;
				partial_ujk_xi_z=(partial_djk_xi_z+(h2c_F[ind]-1.0)*partial_mu_xi_z)/average_adjacent_distance;
				# compute gradient
				scale=-penalty_weight*2/(length*(length-1));
				F=(h2c_radius-pairwise_distance[ind]/average_adjacent_distance);
				term=scale*4*F*F*F;
				#
				x_gradient[ii]+=term*partial_ujk_xi_x;
				y_gradient[ii]+=term*partial_ujk_xi_y;
				z_gradient[ii]+=term*partial_ujk_xi_z;